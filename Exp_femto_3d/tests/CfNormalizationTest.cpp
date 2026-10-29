#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include "TFile.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TROOT.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "exp_femto_3d/CfNormalization.h"
#include "exp_femto_3d/Config.h"
#include "exp_femto_3d/Workflow.h"

using namespace exp_femto_3d;
namespace {
  void Require(bool value, const std::string &message) {
    if (!value)
      throw std::runtime_error(message);
  }
  void Near(double a, double b, const std::string &message, double tolerance = 1.e-10) {
    Require(
        std::isfinite(a) && std::isfinite(b) && std::abs(a - b) <= tolerance * std::max({1., std::abs(a), std::abs(b)}),
        message);
  }
  template <class F>
  void Reject(F action, const std::string &message) {
    bool rejected = false;
    try {
      action();
    } catch (const std::exception &) {
      rejected = true;
    }
    Require(rejected, message);
  }
  CfNormalizationConfig Shell(double low = .15, double high = .20) {
    CfNormalizationConfig c;
    c.mode = CfNormalizationMode::kBoxShell;
    for (auto &range : c.ranges) range = {low, high};
    return c;
  }
  std::array<const TAxis *, 3> Axes(const TH3 &h) {
    return {h.GetXaxis(), h.GetYaxis(), h.GetZaxis()};
  }
  void Fill(TH3 &h, double value) {
    for (int x = 1; x <= h.GetNbinsX(); ++x)
      for (int y = 1; y <= h.GetNbinsY(); ++y)
        for (int z = 1; z <= h.GetNbinsZ(); ++z) {
          h.SetBinContent(x, y, z, value);
          h.SetBinError(h.GetBin(x, y, z), std::sqrt(value));
        }
  }

  void CheckGeometryAndAlgebra() {
    TH3D se("se", "", 5, -.25, .25, 5, -.25, .25, 5, -.25, .25);
    TH3D me("me", "", 5, -.25, .25, 5, -.25, .25, 5, -.25, .25);
    se.Sumw2();
    me.Sumw2();
    Fill(se, 100);
    Fill(me, 100);
    se.SetBinContent(3, 3, 3, 1000);
    se.SetBinContent(0, 0, 0, 1.e12);
    me.SetBinContent(6, 6, 6, 1.e12);
    const auto full = ComputeCfNormalization(se, me, {}), shell = ComputeCfNormalization(se, me, Shell());
    Near(shell.cf_scale, 1., "central enhancement contaminated shell normalization");
    Near(shell.se_integral, 9.8, "shell must contain 125-27=98 bins exactly once");
    Near(full.cf_scale, 12500. / 13400., "legacy integral or under/overflow handling changed");
    Require(shell.outer_first == std::array<int, 3>{1, 1, 1} && shell.inner_first == std::array<int, 3>{2, 2, 2},
            "resolved bin bounds wrong");
    const auto raw_full = ResolvePMLRawScale(se, me, full), raw_shell = ResolvePMLRawScale(se, me, shell);
    Near(ResolvePMLRawScale(se, me, {}), raw_full, "legacy fallback changed");
    // Norm changes with the CF convention; the raw-count prediction and likelihood must not change.
    for (double model : {.8, 1., 1.7})
      Near(model * raw_full,
           model * (shell.cf_scale / full.cf_scale) * raw_shell,
           "PML normalization change altered raw prediction");
    for (int axis = 0; axis < 3; ++axis) {
      TH3D point("point", "", 5, -.25, .25, 5, -.25, .25, 5, -.25, .25);
      int bin[] = {3, 3, 3};
      bin[axis] = 1;
      point.SetBinContent(bin[0], bin[1], bin[2], 1.);
      Near(CfNormalizationIntegral(point, shell), .001, "negative face omitted");
      point.Reset();
      bin[axis] = 5;
      point.SetBinContent(bin[0], bin[1], bin[2], 1.);
      Near(CfNormalizationIntegral(point, shell), .001, "positive face omitted");
    }
    auto unequal = Shell();
    unequal.ranges = {{{.15, .20}, {.05, .10}, {.05, .20}}};
    TH3D unit("unit", "", 5, -.25, .25, 5, -.25, .25, 5, -.25, .25);
    Fill(unit, 1.);
    Near(CfNormalizationIntegral(unit, ResolveCfNormalizationRegion(unequal, Axes(unit))),
         .072,
         "independent axis bounds wrong");
    Near(CfNormalizationIntegral(unit, ResolveCfNormalizationRegion(Shell(.10, .20), Axes(unit))),
         .124,
         "inner endpoint must be included");
    Near(CfNormalizationIntegral(unit, ResolveCfNormalizationRegion(Shell(.05, .10), Axes(unit))),
         .026,
         "outer endpoint must be included");
    Near(CfNormalizationIntegral(unit, ResolveCfNormalizationRegion(Shell(.155, .205), Axes(unit))),
         .098,
         "non-aligned bounds must use centers");
    Near(CfNormalizationIntegral(unit, ResolveCfNormalizationRegion(Shell(0., .20), Axes(unit))),
         .125,
         "empty inner box is valid");
    Reject(
        [&] {
          ResolveCfNormalizationRegion(Shell(.15, .26), Axes(unit));
        },
        "out-of-axis range accepted");
    Reject(
        [&] {
          ResolveCfNormalizationRegion(Shell(.151, .199), Axes(unit));
        },
        "empty bin-center shell accepted");
    auto bad = shell;
    bad.cf_scale *= 2;
    Reject(
        [&] {
          ValidateCfNormalizationResult(bad);
        },
        "inconsistent stored scale accepted");
    bad = shell;
    bad.outer_first[0] = 0;
    Reject(
        [&] {
          ValidateCfNormalizationResult(bad);
        },
        "underflow in metadata accepted");
    TH3D wrong("wrong", "", 4, -.25, .25, 5, -.25, .25, 5, -.25, .25);
    Reject(
        [&] {
          ComputeCfNormalization(unit, wrong, Shell());
        },
        "SE/ME axis mismatch accepted");
    Fill(me, 0.);
    Near(ResolvePMLRawScale(se, me, {}), 0., "legacy zero-count setup-error path changed");
    Reject(
        [&] {
          ComputeCfNormalization(se, me, Shell());
        },
        "zero integral silently normalized");
    me.SetBinContent(1, 1, 1, -1.);
    Reject(
        [&] {
          ComputeCfNormalization(se, me, Shell());
        },
        "negative shell count accepted");
  }

  void CheckNonuniformBins() {
    const double edges[] = {-.30, -.21, -.12, -.03, .02, .07, .14, .22, .30};
    TH3D h("nonuniform", "", 8, edges, 8, edges, 8, edges);
    double expected = 0., counts = 0.;
    for (int x = 1; x <= 8; ++x)
      for (int y = 1; y <= 8; ++y)
        for (int z = 1; z <= 8; ++z) {
          const double content = 10 * x + 2 * y + z;
          h.SetBinContent(x, y, z, content);
          const double qo = std::abs(h.GetXaxis()->GetBinCenter(x));
          const double qs = std::abs(h.GetYaxis()->GetBinCenter(y));
          const double ql = std::abs(h.GetZaxis()->GetBinCenter(z));
          if (std::max({qo, qs, ql}) <= .25 && std::max({qo, qs, ql}) >= .10) {
            counts += content;
            expected += content * (edges[x] - edges[x - 1]) * (edges[y] - edges[y - 1]) * (edges[z] - edges[z - 1]);
          }
        }
    const auto region = ResolveCfNormalizationRegion(Shell(.10, .25), Axes(h));
    Near(CfNormalizationIntegral(h, region), expected, "shell omitted bin volume weighting");
    Require(std::abs(expected - counts) > 1., "nonuniform fixture cannot distinguish count/width conventions");
    Near(CfNormalizationIntegral(h, ResolveCfNormalizationRegion({}, Axes(h))),
         h.Integral(1, 8, 1, 8, 1, 8, "width"),
         "legacy width sum changed");
  }

  void CheckMetadata(const std::filesystem::path &dir) {
    TH3D h("h", "", 5, -.25, .25, 5, -.25, .25, 5, -.25, .25);
    Fill(h, 2.);
    const auto original = ComputeCfNormalization(h, h, Shell());
    const auto path = dir / "metadata.root";
    {
      TFile file(path.c_str(), "RECREATE");
      TTree tree("catalog", "");
      tree.SetDirectory(nullptr);
      CfNormalizationTreeWriter writer(tree);
      writer.Set(original);
      tree.Fill();
      tree.Write();
      TTree legacy("legacy", "");
      legacy.SetDirectory(nullptr);
      int id = 1;
      legacy.Branch("id", &id);
      legacy.Fill();
      legacy.Write();
      TTree partial("partial", "");
      partial.SetDirectory(nullptr);
      partial.Branch("normalization_version", &id);
      partial.Fill();
      partial.Write();
      TTree corrupt("corrupt", "");
      corrupt.SetDirectory(nullptr);
      CfNormalizationTreeWriter bad_writer(corrupt);
      bad_writer.Set(original);
      double bad_scale = -1.;
      corrupt.GetBranch("cf_scale")->SetAddress(&bad_scale);
      corrupt.Fill();
      corrupt.Write();
    }
    TFile file(path.c_str(), "READ");
    auto *tree = dynamic_cast<TTree *>(file.Get("catalog"));
    Require(tree, "metadata tree missing");
    TTreeReader reader(tree);
    CfNormalizationTreeReader normal(*tree, reader);
    Require(reader.Next(), "normalization metadata row missing");
    const auto read = normal.Read();
    Require(read.present && read.config.ranges == original.config.ranges && read.inner_last == original.inner_last,
            "metadata round trip lost geometry");
    Near(read.cf_scale, original.cf_scale, "metadata round trip changed scale");
    auto *legacy = dynamic_cast<TTree *>(file.Get("legacy"));
    TTreeReader lr(legacy);
    CfNormalizationTreeReader old(*legacy, lr);
    Require(lr.Next() && !old.Read().present, "old metadata must retain legacy semantics");
    auto *partial = dynamic_cast<TTree *>(file.Get("partial"));
    TTreeReader pr(partial);
    Reject(
        [&] {
          CfNormalizationTreeReader incomplete(*partial, pr);
        },
        "partial metadata treated as legacy");
    auto *corrupt = dynamic_cast<TTree *>(file.Get("corrupt"));
    TTreeReader cr(corrupt);
    CfNormalizationTreeReader bad(*corrupt, cr);
    Require(cr.Next(), "corrupt fixture missing");
    Reject(
        [&] {
          bad.Read();
        },
        "corrupt scale accepted");
  }

  void WriteToy(const std::filesystem::path &path) {
    TFile file(path.c_str(), "RECREATE");
    auto *task = file.mkdir("task");
    auto *sd = task->mkdir("Same");
    auto *md = task->mkdir("Mixed");
    const int bins[] = {5, 5, 5, 1, 1, 1, 2};
    const double lo[] = {-.25, -.25, -.25, .2, 0., 0., 0.}, hi[] = {.25, .25, .25, .4, 10., 1., kPi};
    THnSparseF se("sparse", "", 7, bins, lo, hi), me("sparse", "", 7, bins, lo, hi);
    se.Sumw2();
    me.Sumw2();
    for (int phi = 0; phi < 2; ++phi)
      for (int x = 0; x < 5; ++x)
        for (int y = 0; y < 5; ++y)
          for (int z = 0; z < 5; ++z) {
            double q[] = {-.20 + .1 * x, -.20 + .1 * y, -.20 + .1 * z, .3, 5., .5, (phi + .5) * kPi / 2};
            const double base = 100000. + 10000 * phi;
            const double exponent = (4 * q[0] * q[0] + 6 * q[1] * q[1] + 9 * q[2] * q[2]) / std::pow(.1973269804, 2);
            const double same = (1.2 + .2 * phi) * base * (1. + .6 * std::exp(-exponent));
            se.Fill(q, same);
            me.Fill(q, base);
            se.SetBinError(se.GetBin(q), std::sqrt(same));
            me.SetBinError(me.GetBin(q), std::sqrt(base));
          }
    sd->cd();
    se.Write();
    md->cd();
    me.Write();
  }

  // Independent count-likelihood check: a wrong normalization can be hidden by a free Norm,
  // so verify the stored objective against raw counts using the FILE's scale, not the build config.
  void CheckFitObjective(const std::filesystem::path &fit_path,
                         const std::filesystem::path &cf_path,
                         bool expect_metadata = true) {
    TFile fit(fit_path.c_str(), "READ"), cf(cf_path.c_str(), "READ");
    auto *tree = dynamic_cast<TTree *>(fit.Get("meta/FitCatalog"));
    Require(tree, "fit catalog missing");
    TTreeReader reader(tree);
    CfNormalizationTreeReader nr(*tree, reader);
    TTreeReaderValue<std::string> slice(reader, "slice_id");
    TTreeReaderValue<double> norm(reader, "norm"), lambda(reader, "lambda"), out(reader, "rout2"),
        side(reader, "rside2"), lng(reader, "rlong2"), objective(reader, "fit_statistic");
    int rows = 0;
    while (reader.Next()) {
      ++rows;
      auto info = nr.Read();
      Require(info.present == expect_metadata, "fit normalization provenance changed");
      const auto prefix = "slices/" + *slice + "/";
      auto *se = dynamic_cast<TH3D *>(cf.Get((prefix + "SE_raw3d").c_str()));
      auto *me = dynamic_cast<TH3D *>(cf.Get((prefix + "ME_raw3d").c_str()));
      Require(se && me, "fit raw fixture missing");
      if (!info.present)
        info = ComputeCfNormalization(*se, *me, {});
      double expected = 0.;
      for (int x = 1; x <= 5; ++x)
        for (int y = 1; y <= 5; ++y)
          for (int z = 1; z <= 5; ++z) {
            const double qo = se->GetXaxis()->GetBinCenter(x), qs = se->GetYaxis()->GetBinCenter(y),
                         ql = se->GetZaxis()->GetBinCenter(z);
            const double model =
                *norm
                * (1.
                   + *lambda
                         * std::exp(-(*out * qo * qo + *side * qs * qs + *lng * ql * ql) / std::pow(.1973269804, 2)));
            const double ratio = model / info.cf_scale, s = se->GetBinContent(x, y, z), m = me->GetBinContent(x, y, z);
            expected +=
                -2. * (s * std::log(ratio * (s + m) / (s * (ratio + 1.))) + m * std::log((s + m) / (m * (ratio + 1.))));
          }
      Near(*objective, expected, "PML objective used a normalization inconsistent with CF metadata", 1.e-7);
    }
    Require(rows == 3, "expected integrated plus two phi fit rows");
  }

  void CheckWorkflow(const std::filesystem::path &dir) {
    WriteToy(dir / "input.root");
    ApplicationConfig config;
    config.input.input_root = (dir / "input.root").string();
    config.input.task_name = "task";
    config.input.same_event_subtask = "Same";
    config.input.mixed_event_subtask = "Mixed";
    config.input.sparse_object_name = "sparse";
    config.output.output_directory = dir.string();
    config.output.fit_report_directory = dir.string();
    config.output.cf_root_name = "full.root";
    config.output.fit_root_name = "fit.root";
    config.output.fit_summary_name = "fit.tsv";
    config.output.fit_report_root_name = "report.root";
    config.centrality_bins = {{0., 10., ""}};
    config.mt_bins = {{.2, .4, ""}};
    config.build.split_mixed_event_by_phi = true;
    config.build.progress = ProgressMode::kDisabled;
    config.fit.progress = ProgressMode::kDisabled;
    config.fit.model = FitModel::kDiag;
    config.fit.options.use_pml = true;
    config.fit.options.fit_q_max = .21;
    config.fit.options.parameters.alpha.fixed_value = 2.;
    config.fit.options.parameters.lambda.fixed_value = .6;
    config.fit.options.parameters.rout2.initial = 4.;
    config.fit.options.parameters.rside2.initial = 6.;
    config.fit.options.parameters.rlong2.initial = 9.;
    ValidateApplicationConfig(config);
    const Logger logger(LogLevel::kError);
    Require(RunBuildCf(config, logger).stored_slices == 3, "full CF build failed");
    const auto full = LoadSliceCatalog((dir / "full.root").string());
    auto shell = config;
    shell.output.cf_root_name = "shell.root";
    shell.build.normalization = Shell();
    Require(RunBuildCf(shell, logger).stored_slices == 3, "shell CF build failed");
    const auto selected = LoadSliceCatalog((dir / "shell.root").string());
    Require(full.size() == selected.size(), "normalization changed slice selection");
    TFile f((dir / "full.root").c_str(), "READ");
    TFile shell_file((dir / "shell.root").c_str(), "READ");
    for (std::size_t n = 0; n < full.size(); ++n) {
      const auto &a = full[n];
      const auto &b = selected[n];
      Require(a.normalization.present && b.normalization.config.mode == CfNormalizationMode::kBoxShell,
              "build metadata missing");
      auto *old = dynamic_cast<TH3D *>(f.Get(a.cf_object_path.c_str()));
      auto *now = dynamic_cast<TH3D *>(shell_file.Get(b.cf_object_path.c_str()));
      auto *raw = dynamic_cast<TH3D *>(f.Get(a.se_object_path.c_str()));
      auto *mixed = dynamic_cast<TH3D *>(f.Get(a.me_object_path.c_str()));
      Require(old && now && raw && mixed, "CF toy objects missing");
      std::unique_ptr<TH3D> expected(static_cast<TH3D *>(raw->Clone("expected")));
      std::unique_ptr<TH3D> denominator(static_cast<TH3D *>(mixed->Clone("denominator")));
      expected->SetDirectory(nullptr);
      denominator->SetDirectory(nullptr);
      expected->Scale(1. / raw->Integral(1, 5, 1, 5, 1, 5, "width"));
      denominator->Scale(1. / mixed->Integral(1, 5, 1, 5, 1, 5, "width"));
      expected->Divide(denominator.get());
      const double scale = b.normalization.cf_scale / a.normalization.cf_scale;
      for (int x = 1; x <= 5; ++x)
        for (int y = 1; y <= 5; ++y)
          for (int z = 1; z <= 5; ++z) {
            const auto bin = old->GetBin(x, y, z);
            Near(old->GetBinContent(bin), expected->GetBinContent(bin), "full CF contents changed");
            Near(old->GetBinError(bin), expected->GetBinError(bin), "full CF errors changed");
            Near(now->GetBinContent(bin), scale * old->GetBinContent(bin), "shell changed CF shape");
            Near(now->GetBinError(bin), scale * old->GetBinError(bin), "shell errors did not scale consistently");
          }
      auto weighted = std::unique_ptr<TH3D>(static_cast<TH3D *>(now->Clone("weighted")));
      weighted->SetDirectory(nullptr);
      auto *shell_me = dynamic_cast<TH3D *>(shell_file.Get(b.me_object_path.c_str()));
      weighted->Multiply(shell_me);
      Near(CfNormalizationIntegral(*weighted, b.normalization) / b.normalization.me_integral,
           1.,
           "ME-weighted shell mean is not unity");
    }
    // Build config remains full_integral: fit-only must nevertheless honor shell.root metadata.
    Require(RunFit(config, logger, std::nullopt, (dir / "shell.root").string()).fitted_slices == 3,
            "fit-only shell input failed");
    CheckFitObjective(dir / "fit.root", dir / "shell.root");
    config.output.fit_root_name = "full_fit.root";
    config.output.fit_summary_name = "full_fit.tsv";
    config.output.fit_report_root_name = "full_report.root";
    Require(RunFit(config, logger).fitted_slices == 3, "legacy-mode fit failed");
    CheckFitObjective(dir / "full_fit.root", dir / "full.root");

    // Recreate a genuine old CF catalog without ANY normalization branches; raw histograms stay intact.
    std::filesystem::copy_file(dir / "full.root", dir / "legacy.root");
    {
      TFile legacy((dir / "legacy.root").c_str(), "UPDATE");
      auto *catalog = dynamic_cast<TTree *>(legacy.Get("meta/SliceCatalog"));
      Require(catalog, "legacy source catalog missing");
      catalog->SetBranchStatus("normalization_*", 0);
      catalog->SetBranchStatus("cf_scale", 0);
      legacy.cd("meta");
      std::unique_ptr<TTree> copy(catalog->CloneTree(-1, "fast"));
      Require(copy != nullptr, "cannot create legacy catalog fixture");
      copy->SetDirectory(nullptr);
      copy->Write("SliceCatalog", TObject::kOverwrite);
    }
    const auto legacy_entries = LoadSliceCatalog((dir / "legacy.root").string());
    Require(legacy_entries.size() == 3 && !legacy_entries.front().normalization.present,
            "legacy file must not fabricate normalization metadata");
    config.output.fit_root_name = "legacy_fit.root";
    config.output.fit_summary_name = "legacy_fit.tsv";
    config.output.fit_report_root_name = "legacy_report.root";
    Require(RunFit(config, logger, std::nullopt, (dir / "legacy.root").string()).fitted_slices == 3,
            "legacy CF PML fit failed");
    CheckFitObjective(dir / "legacy_fit.root", dir / "legacy.root", false);

    // A populated central region must not rescue a zero shell integral.
    WriteToy(dir / "zero_shell.root");
    {
      TFile file((dir / "zero_shell.root").c_str(), "UPDATE");
      auto *mixed = dynamic_cast<THnSparseF *>(file.Get("task/Mixed/sparse"));
      Require(mixed, "zero-shell sparse missing");
      int coordinates[7];
      for (Long64_t bin = 0; bin < mixed->GetNbins(); ++bin) {
        mixed->GetBinContent(bin, coordinates);
        double max_q = 0.;
        for (int axis = 0; axis < 3; ++axis)
          max_q = std::max(max_q, std::abs(mixed->GetAxis(axis)->GetBinCenter(coordinates[axis])));
        if (max_q >= .15)
          mixed->SetBinContent(bin, 0.);
      }
      file.cd("task/Mixed");
      mixed->Write("sparse", TObject::kOverwrite);
    }
    auto zero = shell;
    zero.input.input_root = (dir / "zero_shell.root").string();
    zero.output.cf_root_name = "zero_output.root";
    zero.build.split_mixed_event_by_phi = false;
    const auto zero_stats = RunBuildCf(zero, logger);
    Require(zero_stats.stored_slices == 0 && zero_stats.skipped_zero_mixed_event_groups == 1,
            "zero shell ME must skip its group without full-integral fallback");

    const auto sentinel = dir / "preflight.root";
    {
      std::ofstream out(sentinel);
      out << "preserve this output";
    }
    shell.output.cf_root_name = "preflight.root";
    shell.build.normalization = Shell(.1, .3);
    Reject(
        [&] {
          (void)RunBuildCf(shell, logger);
        },
        "out-of-axis preflight accepted");
    std::ifstream saved(sentinel);
    std::string line;
    std::getline(saved, line);
    Require(line == "preserve this output", "bad normalization truncated existing output");
  }
}  // namespace
int main() {
  try {
    gROOT->SetBatch(true);
    TH1::AddDirectory(false);
    const auto dir =
        std::filesystem::temp_directory_path() / ("exp_femto_3d_normalization_" + std::to_string(getpid()));
    std::filesystem::create_directories(dir);
    CheckGeometryAndAlgebra();
    CheckNonuniformBins();
    CheckMetadata(dir);
    CheckWorkflow(dir);
    std::filesystem::remove_all(dir);
    std::cout << "CF normalization geometry, metadata, build, errors and PML regression passed\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
