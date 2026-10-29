#include "Cats3DModel.h"

#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <sstream>
#include <utility>

#include "TAxis.h"
#include "TFile.h"
#include "TH1.h"
#include "TH3.h"
#include "TTree.h"
#include "exp_femto_3d/Logging.h"
#ifdef EXP_FEMTO_3D_HAS_CATS3D
#include "CATS3D.h"
#include "CATS3DFastKernelCache.h"
#endif

namespace exp_femto_3d {
  namespace {
    [[maybe_unused]] constexpr double kPionMassMeV = 139.57;  // Match the supplied CATS3D reference.
    struct Axis {
      int count = 0;
      double low = 0, high = 0;
      double Center(int i) const {
        return low + (i + .5) * (high - low) / count;
      }
      int Bin(double q) const {
        int i = static_cast<int>(std::floor((q - low) * count / (high - low)));
        if (i < 0 || i >= count) {
          throw Cats3DDomainError("q coordinate outside computed CATS3D grid.");
        }
        return i;
      }
    };
    Axis SelectAxis(const TAxis &axis, double qmax) {
      const double width = axis.GetBinWidth(1);
      if (!std::isfinite(width) || width <= 0) {
        throw std::runtime_error("Invalid CATS3D q axis.");
      }
      int first = 0, last = 0;
      for (int bin = 1; bin <= axis.GetNbins(); ++bin) {
        if (std::abs(axis.GetBinWidth(bin) - width) > 1e-10 * width) {
          throw std::runtime_error("CATS3D requires uniform q axes; no implicit rebinning is performed.");
        }
        if (std::abs(axis.GetBinCenter(bin)) <= qmax) {
          if (!first) {
            first = bin;
          }
          last = bin;
        }
      }
      if (!first) {
        throw std::runtime_error("CATS3D fit box contains no q bins.");
      }
      return {last - first + 1, axis.GetBinLowEdge(first), axis.GetBinUpEdge(last)};
    }
    std::string AxisKey(const std::array<Axis, 3> &axes) {
      std::ostringstream out;
      out << std::setprecision(17);
      for (const auto &a : axes) {
        out << a.count << ':' << a.low << ':' << a.high << '|';
      }
      return out.str();
    }
    std::array<Axis, 3> Axes(const TH3 &h, double qmax) {
      return {SelectAxis(*h.GetXaxis(), qmax), SelectAxis(*h.GetYaxis(), qmax), SelectAxis(*h.GetZaxis(), qmax)};
    }

    struct KernelBudget {
      std::size_t bytes;
      double peak_gib;
    };

    /** The native ByQ loader owns a staging vector until CATS3D copies it.
     * Allow both full float arrays, source weights, and 1 GiB of ROOT/input
     * overhead. This estimate is a pre-allocation budget, not a hard RSS cap.
     */
    KernelBudget CheckBudget(const std::array<Axis, 3> &axes, const Cats3DConfig &cfg, const Logger *logger = nullptr) {
      if (cfg.grid_bins < 2 || cfg.subdivisions < 1 || !std::isfinite(cfg.max_kernel_gib) || cfg.max_kernel_gib <= 0
          || !std::isfinite(cfg.max_peak_gib) || cfg.max_peak_gib <= 0) {
        throw std::runtime_error("Invalid CATS3D memory/grid configuration.");
      }
      constexpr long double gib = 1ULL << 30;
      const long double cells = std::pow(static_cast<long double>(cfg.grid_bins) * cfg.subdivisions, 3);
      const long double nq = static_cast<long double>(axes[0].count) * axes[1].count * axes[2].count;
      const long double bytes = cells * nq * sizeof(float);
      const long double peak = 2 * bytes + 16 * cells + gib;
      const std::string estimate = "kernel GiB=" + std::to_string(double(bytes / gib))
                                   + ", estimated process peak GiB=" + std::to_string(double(peak / gib));
      if (logger)
        logger->Info("CATS3D " + estimate);
      if (!std::isfinite(bytes) || bytes >= std::ldexp(1.0L, 8 * sizeof(std::size_t))
          || cells > std::numeric_limits<unsigned>::max() || nq > std::numeric_limits<unsigned>::max()
          || bytes > cfg.max_kernel_gib * gib || peak > cfg.max_peak_gib * gib) {
        throw std::runtime_error("CATS3D memory budget exceeded before allocation: " + estimate
                                 + ", kernel limit GiB=" + std::to_string(cfg.max_kernel_gib)
                                 + ", process peak limit GiB=" + std::to_string(cfg.max_peak_gib));
      }
      return {static_cast<std::size_t>(bytes), double(peak / gib)};
    }

#ifdef EXP_FEMTO_3D_HAS_CATS3D
    // Readable automatic names are conveniences only. Native file metadata,
    // rather than a library hash or an ExpFemto-specific marker, governs loading.
    std::string NativeCacheName(const std::array<Axis, 3> &axes, const Cats3DConfig &cfg) {
      std::ostringstream out;
      out << std::setprecision(12) << "cats3d";
      for (const auto &axis : axes) out << "_q" << axis.count << '_' << 1000 * axis.low << '_' << 1000 * axis.high;
      out << "_r" << cfg.grid_bins << '_' << cfg.grid_min_fm << '_' << cfg.grid_max_fm << "_sub" << cfg.subdivisions
          << (cfg.with_coulomb ? "_coulomb" : "_plane") << "_pw" << cfg.max_pw;
      return out.str();
    }

    void CheckCacheOrder(const std::filesystem::path &path, const std::array<Axis, 3> &axes) {
      TFile file(path.c_str(), "READ");
      auto *meta = dynamic_cast<TTree *>(file.Get("CATS3DFastKernel6DByQ_meta"));
      auto *tree = dynamic_cast<TTree *>(file.Get("CATS3DFastKernel6DByQ"));
      const auto nq = std::size_t(axes[0].count) * axes[1].count * axes[2].count;
      if (file.IsZombie() || !meta || meta->GetEntries() != 1 || !tree
          || tree->GetEntries() != static_cast<Long64_t>(nq)) {
        throw std::runtime_error("CATS3D native cache is incomplete or has incompatible q dimensions: "
                                 + path.string());
      }
      unsigned index = 0, qo = 0, qs = 0, ql = 0;
      tree->SetBranchStatus("*", false);
      for (const auto &item :
           {std::pair<const char *, unsigned *>{"qIndex", &index}, {"iqOut", &qo}, {"iqSide", &qs}, {"iqLong", &ql}}) {
        tree->SetBranchStatus(item.first, true);
        if (tree->SetBranchAddress(item.first, item.second) < 0)
          throw std::runtime_error("CATS3D native cache missing q index branch.");
      }
      for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
        if (tree->GetEntry(i) <= 0 || index != static_cast<unsigned>(i) || qo >= unsigned(axes[0].count)
            || qs >= unsigned(axes[1].count) || ql >= unsigned(axes[2].count)
            || index != (qo * axes[1].count + qs) * axes[2].count + ql) {
          throw std::runtime_error("CATS3D native cache has missing/duplicate/misordered q entries.");
        }
      }
    }
#endif
  }  // namespace

  Cats3DLambda::Cats3DLambda(const Cats3DLambdaInput &input) : input_(input) {
    if (input.meaning == "pair_purity") {
      if (!input.weak_decay_fraction || !std::isfinite(*input.weak_decay_fraction) || *input.weak_decay_fraction < 0
          || *input.weak_decay_fraction > 1) {
        throw std::runtime_error("pair_purity requires an explicit valid weak_decay_fraction.");
      }
      scale_ = std::pow(1.0 - *input.weak_decay_fraction, 2);
    } else if (input.meaning != "final_lambda" || input.weak_decay_fraction) {
      throw std::runtime_error("Invalid CATS3D lambda meaning/correction.");
    }
    if (input.kind == "constant") {
      if (!input.value || !std::isfinite(*input.value) || *input.value < 0 || *input.value > 1) {
        throw std::runtime_error("CATS3D constant lambda/pair purity must be in [0,1].");
      }
    } else if (input.kind == "histogram") {
      if (input.axis_unit != "MeV" && input.axis_unit != "GeV") {
        throw std::runtime_error("External lambda histogram requires explicit MeV or GeV axis_unit.");
      }
      TFile file(input.file.c_str(), "READ");
      if (file.IsZombie()) {
        throw std::runtime_error("Cannot open external lambda file: " + input.file);
      }
      const auto *h = dynamic_cast<TH1 *>(file.Get(input.histogram.c_str()));
      if (!h || h->GetDimension() != 1 || h->GetNbinsX() < 1) {
        throw std::runtime_error("Missing/invalid external lambda TH1: " + input.histogram);
      }
      histogram_.reset(static_cast<TH1 *>(h->Clone("LambdaInput")));
      histogram_->SetDirectory(nullptr);
      for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
        const double value = h->GetBinContent(bin), error = h->GetBinError(bin);
        if (!std::isfinite(value) || value < 0 || value > 1 || !std::isfinite(error) || error < 0) {
          throw std::runtime_error("Invalid external lambda/purity histogram value or error.");
        }
      }
    } else {
      throw std::runtime_error("Unknown external lambda kind.");
    }
  }
  Cats3DLambda::~Cats3DLambda() = default;
  double Cats3DLambda::Evaluate(double k_mev) const {
    if (!std::isfinite(k_mev) || k_mev < 0) {
      throw std::runtime_error("Invalid lambda-table momentum.");
    }
    if (!histogram_) {
      return *input_.value * scale_;
    }
    const double x = input_.axis_unit == "GeV" ? k_mev / 1000.0 : k_mev;
    const auto &axis = *histogram_->GetXaxis();
    if (x < axis.GetXmin() || x > axis.GetXmax()) {
      throw std::runtime_error("External lambda histogram does not cover required k_table=|q|/2.");
    }
    const int n = histogram_->GetNbinsX();
    if (x <= axis.GetBinCenter(1)) {
      return histogram_->GetBinContent(1) * scale_;
    }
    if (x >= axis.GetBinCenter(n)) {
      return histogram_->GetBinContent(n) * scale_;
    }
    const int bin = axis.FindFixBin(x);
    const int lo = x < axis.GetBinCenter(bin) ? bin - 1 : bin;
    const double fraction = (x - axis.GetBinCenter(lo)) / (axis.GetBinCenter(lo + 1) - axis.GetBinCenter(lo));
    return scale_ * ((1 - fraction) * histogram_->GetBinContent(lo) + fraction * histogram_->GetBinContent(lo + 1));
  }

  struct Cats3DModel::Impl {
    Cats3DConfig config;
    std::array<Axis, 3> axes;
    std::string geometry, contract, key, cache_path;
    std::size_t bytes = 0, fold_count = 0;
    double peak_gib = 0;
    std::array<double, 6> last{}, inverse{};
    double source_norm = 0;
    bool prepared = false;
    std::vector<double> core, nominal_lambda, lambda;
    std::shared_ptr<const Cats3DLambda> input;
    std::optional<double> constant_scan;
    double lambda_scale = 1;
#ifdef EXP_FEMTO_3D_HAS_CATS3D
    std::unique_ptr<CATS3D> kitty;
    static double Source(void *context, double *p) {
      const auto &self = *static_cast<Impl *>(context);
      const double x = p[3], y = p[4], z = p[5];
      const auto &a = self.inverse;
      return self.source_norm
             * std::exp(-.25
                        * (a[0] * x * x + a[1] * y * y + a[2] * z * z + 2 * a[3] * x * y + 2 * a[4] * x * z
                           + 2 * a[5] * y * z));
    }
#endif
    std::size_t Index(int x, int y, int z) const {
      return (std::size_t(x) * axes[1].count + y) * axes[2].count + z;
    }
    std::unique_ptr<TH3D> Histogram(const char *name, const std::vector<double> &values) const {
      auto h = std::make_unique<TH3D>(name,
                                      ";q_{out} (GeV/c);q_{side} (GeV/c);q_{long} (GeV/c)",
                                      axes[0].count,
                                      axes[0].low,
                                      axes[0].high,
                                      axes[1].count,
                                      axes[1].low,
                                      axes[1].high,
                                      axes[2].count,
                                      axes[2].low,
                                      axes[2].high);
      h->SetDirectory(nullptr);
      h->Sumw2();  // Evaluated model/map bins have zero conditional numerical errors.
      for (int x = 0; x < axes[0].count; ++x) {
        for (int y = 0; y < axes[1].count; ++y) {
          for (int z = 0; z < axes[2].count; ++z) {
            h->SetBinContent(x + 1, y + 1, z + 1, values.at(Index(x, y, z)));
          }
        }
      }
      return h;
    }
  };

  bool Cats3DModel::Available() {
#ifdef EXP_FEMTO_3D_HAS_CATS3D
    return true;
#else
    return false;
#endif
  }
  void Cats3DModel::ValidateInputGrid(const TH3 &reference,
                                      double qmax,
                                      const Cats3DConfig &config,
                                      const Cats3DLambda &input) {
    const auto axes = Axes(reference, qmax);
    (void)CheckBudget(axes, config);
    // The table interval is contiguous and interpolation is convex. Endpoints of
    // the actual grid's norm range therefore prove coverage for every q bin.
    double min2 = 0, max2 = 0;
    for (const auto &axis : axes) {
      double minimum = std::numeric_limits<double>::infinity(), maximum = 0;
      for (int i = 0; i < axis.count; ++i) {
        minimum = std::min(minimum, std::abs(axis.Center(i)));
        maximum = std::max(maximum, std::abs(axis.Center(i)));
      }
      min2 += minimum * minimum;
      max2 += maximum * maximum;
    }
    (void)input.Evaluate(500 * std::sqrt(min2));
    (void)input.Evaluate(500 * std::sqrt(max2));
  }
  Cats3DModel::Cats3DModel(const TH3 &reference, double qmax, const Cats3DConfig &config, const Logger *logger)
      : impl_(std::make_unique<Impl>()) {
    auto &s = *impl_;
    s.config = config;
    s.axes = Axes(reference, qmax);
    s.geometry = AxisKey(s.axes);
    if (config.grid_bins < 2 || config.subdivisions < 1 || config.max_pw < 1 || config.max_pw > 65535
        || !std::isfinite(config.grid_min_fm) || !std::isfinite(config.grid_max_fm) || config.grid_min_fm >= 0
        || config.grid_max_fm <= 0 || !std::isfinite(config.max_kernel_gib) || config.max_kernel_gib <= 0) {
      throw std::runtime_error("Invalid CATS3D numerical configuration.");
    }
    if (!config.kernel_cache_file.empty() && !config.use_disk_cache)
      throw std::runtime_error("CATS3D kernel_cache_file requires use_disk_cache=true.");
    const auto budget = CheckBudget(s.axes, config, logger);
    s.bytes = budget.bytes;
    s.peak_gib = budget.peak_gib;
#ifdef EXP_FEMTO_3D_HAS_CATS3D
    const std::size_t nq = std::size_t(s.axes[0].count) * s.axes[1].count * s.axes[2].count;
    std::ostringstream contract;
    contract << std::setprecision(17) << "native_cats3d_byq_v1|q_GeV_to_MeV_no_half|k_table_qnorm_half|" << s.geometry
             << config.grid_bins << '|' << config.grid_min_fm << '|' << config.grid_max_fm << '|' << config.subdivisions
             << '|' << config.max_pw << '|' << config.with_coulomb << "|pdg=211,211|qs=1|mass=" << kPionMassMeV;
    s.contract = contract.str();
    s.key = NativeCacheName(s.axes, config);
    s.kitty = std::make_unique<CATS3D>();
    auto &cat = *s.kitty;
    cat.SetMomBinsOut(s.axes[0].count, 1000 * s.axes[0].low, 1000 * s.axes[0].high);
    cat.SetMomBinsSide(s.axes[1].count, 1000 * s.axes[1].low, 1000 * s.axes[1].high);
    cat.SetMomBinsLong(s.axes[2].count, 1000 * s.axes[2].low, 1000 * s.axes[2].high);
    cat.SetGridX(config.grid_bins, config.grid_min_fm, config.grid_max_fm);
    cat.SetGridY(config.grid_bins, config.grid_min_fm, config.grid_max_fm);
    cat.SetGridZ(config.grid_bins, config.grid_min_fm, config.grid_max_fm);
    cat.SetPdgId(211, 211);
    cat.SetQuantumStatistics(1);
    cat.SetQ1Q2(config.with_coulomb ? 1 : 0);
    cat.SetRedMass(.5 * kPionMassMeV);
    cat.SetMaxPw(config.max_pw);
    cat.SetWfMode(config.with_coulomb ? CATS3D::wfCoulomb : CATS3D::wfPlaneWave);
    cat.SetFastKernelMode(true);
    cat.SetFastKernelSubdivisions(config.subdivisions);
    // Install once. Updating only analytic parameters must not invalidate the WF kernel.
    cat.SetAnaSource(Impl::Source, &s, 6);
    if (config.use_disk_cache) {
      const bool external = !config.kernel_cache_file.empty();
      const std::filesystem::path dir =
          config.kernel_cache_directory.empty() ? EXP_FEMTO_3D_CATS3D_CACHE_DIR : config.kernel_cache_directory;
      const std::filesystem::path path =
          external ? std::filesystem::path(config.kernel_cache_file) : dir / (s.key + ".root");
      s.cache_path = path.string();
      if (external && !std::filesystem::is_regular_file(path))
        throw std::runtime_error("CATS3D explicit kernel_cache_file does not exist: " + path.string());
      if (std::filesystem::exists(path)) {
        CheckCacheOrder(path, s.axes);
        // Native CATS3D validates axes, r centers, mass, charge, QS and partial
        // waves before allocating a kernel. Arbitrary filenames are supported.
        if (!LoadCATS3DFastKernelByQ(cat, path.string().c_str()))
          throw std::runtime_error("CATS3D native cache grid/physics mismatch or corruption: " + path.string());
        const auto &kernel = cat.GetFastWaveFunctionKernel();
        if (kernel.size() * sizeof(float) != s.bytes || !std::all_of(kernel.begin(), kernel.end(), [](float v) {
              return std::isfinite(v) && v >= 0;
            }))
          throw std::runtime_error("CATS3D native cache contains invalid kernel values.");
      } else {
        cat.BuildFastKernelOnly();
        const auto &kernel = cat.GetFastWaveFunctionKernel();
        if (!std::all_of(kernel.begin(), kernel.end(), [](float v) {
              return std::isfinite(v) && v >= 0;
            }))
          throw std::runtime_error("CATS3D generated invalid wave-function kernel.");
        std::filesystem::create_directories(dir);
        const auto temp = std::filesystem::path(path.string() + ".tmp." + std::to_string(getpid()));
        try {
          // The saved file is the unmodified native ByQ format. No custom
          // fingerprint or library binary signature is required to read it.
          if (!SaveCATS3DFastKernelByQ(cat, temp.string().c_str(), "CATS3DFastKernel6DByQ", false))
            throw std::runtime_error("Cannot save CATS3D kernel cache.");
          std::filesystem::rename(temp, path);
        } catch (...) {
          std::filesystem::remove(temp);
          throw;
        }
      }
    }
    s.core.resize(nq);
#else
    (void)logger;
    throw std::runtime_error("correlation_backend=cats3d requires a build with CATS3D support.");
#endif
  }
  Cats3DModel::~Cats3DModel() = default;
  bool Cats3DModel::Matches(const TH3 &h, double qmax) const {
    return impl_->geometry == AxisKey(Axes(h, qmax));
  }
  void Cats3DModel::SetLambda(std::shared_ptr<const Cats3DLambda> input) {
    auto &s = *impl_;
    if (!input) {
      throw std::runtime_error("Missing CATS3D lambda input.");
    }
    s.input = std::move(input);
    s.nominal_lambda.clear();
    for (int x = 0; x < s.axes[0].count; ++x) {
      for (int y = 0; y < s.axes[1].count; ++y) {
        for (int z = 0; z < s.axes[2].count; ++z) {
          const double qo = s.axes[0].Center(x), qs = s.axes[1].Center(y), ql = s.axes[2].Center(z);
          s.nominal_lambda.push_back(s.input->Evaluate(500 * std::sqrt(qo * qo + qs * qs + ql * ql)));
        }
      }
    }
    SetScan(std::nullopt);
  }
  void Cats3DModel::SetScan(std::optional<double> constant, double scale) {
    auto &s = *impl_;
    if (!s.input) {
      throw std::runtime_error("Set lambda input before a scan.");
    }
    if (!std::isfinite(scale) || scale < 0
        || (constant && (!std::isfinite(*constant) || *constant < 0 || *constant > 1))) {
      throw Cats3DDomainError("Invalid fixed lambda scan coordinate.");
    }
    if (!constant && scale != 1 && s.input->Input().kind != "histogram") {
      throw std::runtime_error("Lambda scale scan requires a histogram input.");
    }
    auto values = s.nominal_lambda;
    for (auto &value : values) {
      value = constant ? *constant : scale * value;
      if (!std::isfinite(value) || value < 0 || value > 1) {
        throw Cats3DDomainError("Lambda scan produces values outside [0,1]; no clipping is applied.");
      }
    }
    s.constant_scan = constant;
    s.lambda_scale = scale;
    s.lambda = std::move(values);
  }
  void Cats3DModel::Prepare(const double *p, bool full) {
    auto &s = *impl_;
    const std::array<double, 6> a = {p[2], p[3], p[4], full ? p[5] : 0, full ? p[6] : 0, full ? p[7] : 0};
    if (s.prepared && a == s.last) {
      return;
    }
    s.prepared = false;
    // Sylvester/Cholesky pivots enforce a strictly positive definite source matrix.
    const double minor = a[0] * a[1] - a[3] * a[3];
    const double det =
        a[0] * (a[1] * a[2] - a[5] * a[5]) - a[3] * (a[3] * a[2] - a[4] * a[5]) + a[4] * (a[3] * a[5] - a[4] * a[1]);
    if (!std::all_of(a.begin(),
                     a.end(),
                     [](double v) {
                       return std::isfinite(v);
                     })
        || !(a[0] > 0) || !(minor > 0) || !(det > 0) || !std::isfinite(det)) {
      throw Cats3DDomainError("Gaussian R2 matrix is not positive definite.");
    }
    s.inverse = {(a[1] * a[2] - a[5] * a[5]) / det,
                 (a[0] * a[2] - a[4] * a[4]) / det,
                 minor / det,
                 (a[4] * a[5] - a[3] * a[2]) / det,
                 (a[3] * a[5] - a[4] * a[1]) / det,
                 (a[3] * a[4] - a[0] * a[5]) / det};
    s.source_norm = 1.0 / (std::pow(4 * kPi, 1.5) * std::sqrt(det));
    if (!std::isfinite(s.source_norm) || !std::all_of(s.inverse.begin(), s.inverse.end(), [](double v) {
          return std::isfinite(v);
        })) {
      throw Cats3DDomainError("Gaussian source matrix inversion is not finite.");
    }
#ifdef EXP_FEMTO_3D_HAS_CATS3D
    for (unsigned i = 0; i < 6; ++i) {
      s.kitty->SetAnaSource(i, a[i], true);
    }
    s.kitty->KillTheCat();
    for (int x = 0; x < s.axes[0].count; ++x) {
      for (int y = 0; y < s.axes[1].count; ++y) {
        for (int z = 0; z < s.axes[2].count; ++z) {
          const double value = s.kitty->GetCorrFun(x, y, z);
          if (!std::isfinite(value) || value < 0 || s.kitty->GetCorrFunErr(x, y, z) >= 1e6) {
            throw Cats3DDomainError("CATS3D produced an empty/invalid source convolution.");
          }
          s.core[s.Index(x, y, z)] = value;
        }
      }
    }
#endif
    s.last = a;
    s.prepared = true;
    ++s.fold_count;
  }
  double Cats3DModel::Evaluate(double qo, double qs, double ql, const double *p, bool full, bool baseline) {
    Prepare(p, full);
    auto &s = *impl_;
    if (!s.input) {
      throw std::runtime_error("CATS3D lambda input is not set.");
    }
    const auto i = s.Index(s.axes[0].Bin(qo), s.axes[1].Bin(qs), s.axes[2].Bin(ql));
    const double b = baseline ? 1 + p[full ? 9 : 6] * (qo * qo + qs * qs + ql * ql) : 1;
    const double value = p[0] * b * (1 + s.lambda.at(i) * (s.core.at(i) - 1));
    if (!std::isfinite(value) || value <= 0) {
      throw Cats3DDomainError("CATS3D model ratio is non-positive or non-finite.");
    }
    return value;
  }
  double Cats3DModel::ScalarLambda() const {
    const auto &s = *impl_;
    return s.constant_scan
               ? *s.constant_scan
               : (s.input && s.input->Input().kind == "constant" ? s.lambda_scale * s.input->Evaluate(0)
                                                                 : std::numeric_limits<double>::quiet_NaN());
  }
  std::unique_ptr<TH3D> Cats3DModel::CoreHistogram() const {
    if (!impl_->prepared) {
      throw std::runtime_error("CATS3D core has not been evaluated.");
    }
    return impl_->Histogram("CATS3DCoreCF3D", impl_->core);
  }
  std::unique_ptr<TH3D> Cats3DModel::LambdaHistogram() const {
    return impl_->Histogram("LambdaMap3D", impl_->lambda);
  }
  std::unique_ptr<TH3D> Cats3DModel::FinalHistogram(const double *p, bool full, bool baseline) {
    auto &s = *impl_;
    std::vector<double> values;
    for (int x = 0; x < s.axes[0].count; ++x) {
      for (int y = 0; y < s.axes[1].count; ++y) {
        for (int z = 0; z < s.axes[2].count; ++z) {
          values.push_back(Evaluate(s.axes[0].Center(x), s.axes[1].Center(y), s.axes[2].Center(z), p, full, baseline));
        }
      }
    }
    return s.Histogram("ModelCF3D", values);
  }
  const Cats3DLambda &Cats3DModel::LambdaInput() const {
    if (!impl_->input) {
      throw std::runtime_error("Missing lambda input.");
    }
    return *impl_->input;
  }
  const std::string &Cats3DModel::Contract() const {
    return impl_->contract;
  }
  const std::string &Cats3DModel::CacheKey() const {
    return impl_->key;
  }
  std::size_t Cats3DModel::KernelBytes() const {
    return impl_->bytes;
  }
  const std::string &Cats3DModel::KernelCachePath() const {
    return impl_->cache_path;
  }
  double Cats3DModel::EstimatedPeakGiB() const {
    return impl_->peak_gib;
  }
  bool Cats3DModel::WithCoulomb() const {
    return impl_->config.with_coulomb;
  }
  std::size_t Cats3DModel::FoldCount() const {
    return impl_->fold_count;
  }
  int Cats3DModel::PmlMaxCalls() const {
    return impl_->config.pml_max_calls;
  }
}  // namespace exp_femto_3d
