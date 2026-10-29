#include "exp_femto_3d/CfNormalization.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "TAxis.h"
#include "TBranch.h"
#include "TH3.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

namespace exp_femto_3d {
  namespace {
    constexpr const char *kAxes[] = {"q_out", "q_side", "q_long"};
    std::array<const TAxis *, 3> Axes(const TH3 &h) {
      return {h.GetXaxis(), h.GetYaxis(), h.GetZaxis()};
    }
    bool Close(double a, double b) {
      return std::abs(a - b) <= 1.e-12 * std::max({std::abs(a), std::abs(b), 1.e-300});
    }
    void CheckIntegral(double value) {
      if (!std::isfinite(value) || value < 0.)
        throw std::runtime_error("CF normalization integral is negative or non-finite.");
    }
  }  // namespace

  std::string ToString(CfNormalizationMode mode) {
    if (mode == CfNormalizationMode::kFullIntegral)
      return "full_integral";
    if (mode == CfNormalizationMode::kBoxShell)
      return "box_shell";
    throw std::runtime_error("Unknown CF normalization mode.");
  }

  void ValidateCfNormalizationAxes(const std::array<const TAxis *, 3> &se, const std::array<const TAxis *, 3> &me) {
    for (int a = 0; a < 3; ++a) {
      if (!se[a] || !me[a] || se[a]->GetNbins() != me[a]->GetNbins())
        throw std::runtime_error("SE/ME normalization axes have different bin counts.");
      for (int b = 1; b <= se[a]->GetNbins() + 1; ++b) {
        const double se_edge = se[a]->GetBinLowEdge(b), me_edge = me[a]->GetBinLowEdge(b);
        if (!std::isfinite(se_edge) || !std::isfinite(me_edge)
            || (b > 1 && (se_edge <= se[a]->GetBinLowEdge(b - 1) || me_edge <= me[a]->GetBinLowEdge(b - 1))))
          throw std::runtime_error("CF normalization axes must have finite increasing edges.");
        if (std::abs(se_edge - me_edge) > kCfNormalizationTolerance)
          throw std::runtime_error(std::string("SE/ME normalization axis edges disagree: ") + kAxes[a]);
      }
    }
  }

  void ValidateCfNormalizationConfig(const CfNormalizationConfig &config) {
    ToString(config.mode);
    for (const auto &range : config.ranges) {
      if (config.mode == CfNormalizationMode::kFullIntegral) {
        if (range[0] != 0. || range[1] != 0.)
          throw ConfigError("full_integral does not accept normalization axis ranges.");
      } else if (!std::isfinite(range[0]) || !std::isfinite(range[1]) || range[0] < 0. || range[1] <= range[0]) {
        throw ConfigError("box_shell ranges require finite 0 <= min < max (GeV/c).");
      }
    }
  }

  CfNormalizationResult ResolveCfNormalizationRegion(const CfNormalizationConfig &config,
                                                     const std::array<const TAxis *, 3> &axes) {
    ValidateCfNormalizationConfig(config);
    CfNormalizationResult result;
    result.config = config;
    long double outer_count = 1, inner_count = 1;
    for (int a = 0; a < 3; ++a) {
      if (!axes[a])
        throw std::runtime_error("Missing CF normalization axis.");
      const auto &axis = *axes[a];
      int of = 0, ol = 0, inf = 0, inl = 0;
      const double low = config.ranges[a][0], high = config.ranges[a][1];
      if (config.mode == CfNormalizationMode::kBoxShell) {
        if (!std::isfinite(low) || !std::isfinite(high) || low < 0. || high <= low)
          throw std::runtime_error("CF box_shell range requires finite 0 <= min < max.");
        if (axis.GetBinLowEdge(1) > -high + kCfNormalizationTolerance
            || axis.GetBinUpEdge(axis.GetNbins()) < high - kCfNormalizationTolerance)
          throw std::runtime_error(std::string("CF normalization outer box is outside axis ") + kAxes[a]);
      } else if (config.mode != CfNormalizationMode::kFullIntegral) {
        throw std::runtime_error("Unknown CF normalization mode.");
      }
      for (int b = 1; b <= axis.GetNbins(); ++b) {
        const double q = std::abs(axis.GetBinCenter(b));
        if (config.mode == CfNormalizationMode::kFullIntegral || q <= high + kCfNormalizationTolerance) {
          if (!of)
            of = b;
          ol = b;
          // A center at an inner boundary belongs to the shell, not the excluded box.
          if (config.mode == CfNormalizationMode::kBoxShell && q < low - kCfNormalizationTolerance) {
            if (!inf)
              inf = b;
            inl = b;
          }
        }
      }
      result.outer_first[a] = of;
      result.outer_last[a] = ol;
      result.inner_first[a] = inf;
      result.inner_last[a] = inl;
      outer_count *= of ? ol - of + 1 : 0;
      inner_count *= inf ? inl - inf + 1 : 0;
    }
    if (outer_count <= inner_count)
      throw std::runtime_error("CF normalization region contains no normal bins.");
    return result;
  }

  double CfNormalizationIntegral(const TH3 &hist, const CfNormalizationResult &region) {
    // Use ROOT's original accumulation for full_integral to preserve old contents and errors exactly.
    if (region.config.mode == CfNormalizationMode::kFullIntegral) {
      const double sum = hist.Integral(1, hist.GetNbinsX(), 1, hist.GetNbinsY(), 1, hist.GetNbinsZ(), "width");
      CheckIntegral(sum);
      return sum;
    }
    double sum = 0.;
    for (int x = region.outer_first[0]; x <= region.outer_last[0]; ++x)
      for (int y = region.outer_first[1]; y <= region.outer_last[1]; ++y)
        for (int z = region.outer_first[2]; z <= region.outer_last[2]; ++z) {
          if (x >= region.inner_first[0] && x <= region.inner_last[0] && y >= region.inner_first[1]
              && y <= region.inner_last[1] && z >= region.inner_first[2] && z <= region.inner_last[2])
            continue;
          const double content = hist.GetBinContent(x, y, z);
          if (!std::isfinite(content))
            throw std::runtime_error("CF normalization contains a non-finite bin.");
          sum += content * hist.GetXaxis()->GetBinWidth(x) * hist.GetYaxis()->GetBinWidth(y)
                 * hist.GetZaxis()->GetBinWidth(z);
        }
    CheckIntegral(sum);
    return sum;
  }

  CfNormalizationResult ComputeCfNormalization(const TH3 &se, const TH3 &me, const CfNormalizationConfig &config) {
    ValidateCfNormalizationAxes(Axes(se), Axes(me));
    auto result = ResolveCfNormalizationRegion(config, Axes(se));
    result.se_integral = CfNormalizationIntegral(se, result);
    result.me_integral = CfNormalizationIntegral(me, result);
    if (result.se_integral <= 0. || result.me_integral <= 0.)
      throw std::runtime_error("Cannot normalize CF with zero SE/ME integral.");
    result.cf_scale = result.me_integral / result.se_integral;
    result.present = true;
    ValidateCfNormalizationResult(result);
    return result;
  }

  void ValidateCfNormalizationResult(const CfNormalizationResult &r) {
    if (!r.present)
      return;
    ToString(r.config.mode);
    if (!(r.se_integral > 0.) || !(r.me_integral > 0.) || !(r.cf_scale > 0.) || !std::isfinite(r.se_integral)
        || !std::isfinite(r.me_integral) || !std::isfinite(r.cf_scale)
        || !Close(r.cf_scale, r.me_integral / r.se_integral))
      throw std::runtime_error("Corrupt CF normalization integrals or cf_scale.");
    long double outer = 1, inner = 1;
    for (int a = 0; a < 3; ++a) {
      const auto &bounds = r.config.ranges[a];
      if (r.config.mode == CfNormalizationMode::kBoxShell) {
        if (!std::isfinite(bounds[0]) || !std::isfinite(bounds[1]) || bounds[0] < 0. || bounds[1] <= bounds[0])
          throw std::runtime_error("Corrupt CF normalization ranges.");
      } else if (bounds[0] != 0. || bounds[1] != 0.) {
        throw std::runtime_error("full_integral metadata must not contain axis ranges.");
      }
      if (r.outer_first[a] < 1 || r.outer_last[a] < r.outer_first[a]
          || ((r.inner_first[a] == 0) != (r.inner_last[a] == 0))
          || (r.inner_first[a] != 0
              && (r.inner_first[a] < r.outer_first[a] || r.inner_last[a] > r.outer_last[a]
                  || r.inner_last[a] < r.inner_first[a]))
          || (r.config.mode == CfNormalizationMode::kFullIntegral && r.inner_first[a] != 0))
        throw std::runtime_error("Corrupt CF normalization bin bounds.");
      outer *= r.outer_last[a] - r.outer_first[a] + 1;
      inner *= r.inner_first[a] ? r.inner_last[a] - r.inner_first[a] + 1 : 0;
    }
    if (outer <= inner)
      throw std::runtime_error("Corrupt empty CF normalization region.");
  }

  double ResolvePMLRawScale(const TH3 &se, const TH3 &me, const CfNormalizationResult &stored) {
    if (!stored.present) {
      ValidateCfNormalizationAxes(Axes(se), Axes(me));
      const auto region = ResolveCfNormalizationRegion({}, Axes(se));
      const double same = CfNormalizationIntegral(se, region), mixed = CfNormalizationIntegral(me, region);
      // Preserve the legacy minimizer setup-error path for zero-count inputs.
      return same > 0. && mixed > 0. ? same / mixed : 0.;
    }
    ValidateCfNormalizationResult(stored);
    ValidateCfNormalizationAxes(Axes(se), Axes(me));
    const auto resolved = ResolveCfNormalizationRegion(stored.config, Axes(se));
    if (resolved.outer_first != stored.outer_first || resolved.outer_last != stored.outer_last
        || resolved.inner_first != stored.inner_first || resolved.inner_last != stored.inner_last)
      throw std::runtime_error("CF normalization metadata disagrees with raw histogram axes.");
    // C_model is normalized; the count likelihood needs raw S/B = C_model / cf_scale.
    return stored.se_integral / stored.me_integral;
  }

  CfNormalizationTreeWriter::CfNormalizationTreeWriter(TTree &tree) {
    tree.Branch("normalization_version", &version_);
    tree.Branch("normalization_mode", &mode_);
    for (int a = 0; a < 3; ++a) {
      const std::string name = std::string("normalization_") + kAxes[a];
      tree.Branch(name.c_str(), value_.config.ranges[a].data(), (name + "[2]/D").c_str());
    }
    const auto branch_bins = [&](const char *name, std::array<int, 3> &bins) {
      tree.Branch(name, bins.data(), (std::string(name) + "[3]/I").c_str());
    };
    branch_bins("normalization_outer_first", value_.outer_first);
    branch_bins("normalization_outer_last", value_.outer_last);
    branch_bins("normalization_inner_first", value_.inner_first);
    branch_bins("normalization_inner_last", value_.inner_last);
    tree.Branch("normalization_se_integral", &value_.se_integral);
    tree.Branch("normalization_me_integral", &value_.me_integral);
    tree.Branch("cf_scale", &value_.cf_scale);
  }
  void CfNormalizationTreeWriter::Set(const CfNormalizationResult &result) {
    if (!result.present)
      throw std::runtime_error("Cannot mix legacy and new CF normalization metadata.");
    ValidateCfNormalizationResult(result);
    value_ = result;
    mode_ = ToString(result.config.mode);
  }

  struct CfNormalizationTreeReader::Impl {
    TTreeReaderValue<int> version;
    TTreeReaderValue<std::string> mode;
    std::vector<std::unique_ptr<TTreeReaderArray<double>>> ranges;
    std::vector<std::unique_ptr<TTreeReaderArray<int>>> bins;
    TTreeReaderValue<double> se, me, scale;
    explicit Impl(TTreeReader &r)
        : version(r, "normalization_version"),
          mode(r, "normalization_mode"),
          se(r, "normalization_se_integral"),
          me(r, "normalization_me_integral"),
          scale(r, "cf_scale") {
      for (const auto *axis : kAxes)
        ranges.emplace_back(
            std::make_unique<TTreeReaderArray<double>>(r, (std::string("normalization_") + axis).c_str()));
      for (const auto *name : {"normalization_outer_first",
                               "normalization_outer_last",
                               "normalization_inner_first",
                               "normalization_inner_last"})
        bins.emplace_back(std::make_unique<TTreeReaderArray<int>>(r, name));
    }
  };

  CfNormalizationTreeReader::CfNormalizationTreeReader(TTree &tree, TTreeReader &reader) {
    const std::vector<std::string> names = {"normalization_version",
                                            "normalization_mode",
                                            "normalization_q_out",
                                            "normalization_q_side",
                                            "normalization_q_long",
                                            "normalization_outer_first",
                                            "normalization_outer_last",
                                            "normalization_inner_first",
                                            "normalization_inner_last",
                                            "normalization_se_integral",
                                            "normalization_me_integral",
                                            "cf_scale"};
    const auto count = std::count_if(names.begin(), names.end(), [&](const auto &n) {
      return tree.GetBranch(n.c_str()) != nullptr;
    });
    if (!count)
      return;  // Only complete absence is legacy, not partial metadata.
    if (count != static_cast<int>(names.size()))
      throw std::runtime_error("Incomplete CF normalization metadata.");
    for (const auto &name : names) {
      if (name == "normalization_mode") {
        const std::string type = tree.GetBranch(name.c_str())->GetClassName();
        if (type != "string" && type != "std::string")
          throw std::runtime_error("Corrupt CF normalization mode branch type.");
        continue;
      }
      const auto *leaf = tree.GetLeaf(name.c_str());
      const bool bin = name.find("_first") != std::string::npos || name.find("_last") != std::string::npos;
      const bool range =
          name == "normalization_q_out" || name == "normalization_q_side" || name == "normalization_q_long";
      const char *type = bin || name == "normalization_version" ? "Int_t" : "Double_t";
      if (!leaf || std::string(leaf->GetTypeName()) != type || leaf->GetLenStatic() != (bin ? 3 : range ? 2 : 1))
        throw std::runtime_error("Corrupt CF normalization branch: " + name);
    }
    impl_ = std::make_unique<Impl>(reader);
  }
  CfNormalizationTreeReader::~CfNormalizationTreeReader() = default;
  CfNormalizationResult CfNormalizationTreeReader::Read() const {
    if (!impl_)
      return {};
    auto &v = *impl_;
    if (*v.version != 1)
      throw std::runtime_error("Unsupported CF normalization metadata version.");
    CfNormalizationResult result;
    result.present = true;
    if (*v.mode == "box_shell")
      result.config.mode = CfNormalizationMode::kBoxShell;
    else if (*v.mode != "full_integral")
      throw std::runtime_error("Corrupt CF normalization mode.");
    for (int a = 0; a < 3; ++a) {
      if (v.ranges[a]->GetSize() != 2)
        throw std::runtime_error("Corrupt normalization range length.");
      result.config.ranges[a] = {(*v.ranges[a])[0], (*v.ranges[a])[1]};
    }
    std::array<int, 3> *bins[] = {&result.outer_first, &result.outer_last, &result.inner_first, &result.inner_last};
    for (int b = 0; b < 4; ++b) {
      if (v.bins[b]->GetSize() != 3)
        throw std::runtime_error("Corrupt normalization bin array length.");
      for (int a = 0; a < 3; ++a) (*bins[b])[a] = (*v.bins[b])[a];
    }
    result.se_integral = *v.se;
    result.me_integral = *v.me;
    result.cf_scale = *v.scale;
    ValidateCfNormalizationResult(result);
    return result;
  }
}  // namespace exp_femto_3d
