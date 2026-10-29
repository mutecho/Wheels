#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

#include "Cats3DModel.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TNamed.h"
#include "TROOT.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "exp_femto_3d/Config.h"
#include "exp_femto_3d/Logging.h"
#include "exp_femto_3d/Workflow.h"
using namespace exp_femto_3d;
namespace {
  void Require(bool value, const std::string &message) {
    if (!value) {
      throw std::runtime_error(message);
    }
  }
  template <class F>
  void Reject(F fn, const char *message) {
    bool failed = false;
    try {
      fn();
    } catch (const std::exception &) {
      failed = true;
    }
    Require(failed, message);
  }
  std::string Bytes(const std::string &path) {
    std::ifstream f(path, std::ios::binary);
    std::ostringstream s;
    s << f.rdbuf();
    return s.str();
  }
  void WriteToy(const std::string &path) {
    TFile file(path.c_str(), "RECREATE");
    auto *task = file.mkdir("task");
    auto *same_dir = task->mkdir("Same");
    auto *mixed_dir = task->mkdir("Mixed");
    const int bins[] = {4, 4, 4, 2, 1, 1, 3};
    const double lo[] = {-.12, -.12, -.12, .2, 0, -.5, 0}, hi[] = {.12, .12, .12, .6, 10, .5, kPi};
    THnSparseF same("sparse", "", 7, bins, lo, hi), mixed("sparse", "", 7, bins, lo, hi);
    same.Sumw2();
    mixed.Sumw2();
    for (int mt = 0; mt < 2; ++mt) {
      for (int phi = 0; phi < 3; ++phi) {
        for (int x = 0; x < 4; ++x) {
          for (int y = 0; y < 4; ++y) {
            for (int z = 0; z < 4; ++z) {
              double v[] = {-.09 + .06 * x, -.09 + .06 * y, -.09 + .06 * z, .3 + .2 * mt, 5, 0, (phi + .5) * kPi / 3};
              const double lambda = mt == 0 ? .7 : .55;
              const double argument = ((16 + mt) * v[0] * v[0] + (20.25 + mt) * v[1] * v[1] + (25 + mt) * v[2] * v[2])
                                      / std::pow(.1973269804, 2);
              const double m = 200000 + 10000 * phi, s = std::round(m * (1 + lambda * std::exp(-argument)));
              same.Fill(v, s);
              mixed.Fill(v, m);
              same.SetBinError(same.GetBin(v), std::sqrt(s));
              mixed.SetBinError(mixed.GetBin(v), std::sqrt(m));
            }
          }
        }
      }
    }
    same_dir->cd();
    same.Write();
    mixed_dir->cd();
    mixed.Write();
  }
  void CheckFit(const std::string &path, int expected, bool table, bool full) {
    TFile file(path.c_str(), "READ");
    auto *tree = dynamic_cast<TTree *>(file.Get("meta/FitCatalog"));
    Require(tree && tree->GetEntries() == expected, "missing fit catalog rows");
    TTreeReader reader(tree);
    TTreeReaderValue<std::string> slice(reader, "slice_id"), backend(reader, "correlation_backend");
    TTreeReaderValue<int> fixed(reader, "lambda_is_external_fixed"), valid(reader, "point_estimate_valid"),
        errors(reader, "parameter_errors_valid"), mt(reader, "mt_index"), ndf(reader, "ndf");
    TTreeReaderValue<double> rout(reader, "rout2"), rside(reader, "rside2"), rlong(reader, "rlong2"),
        lambda(reader, "lambda");
    while (reader.Next()) {
      Require(*backend == "cats3d" && *fixed == 1, "backend or external-fixed provenance missing");
      Require(*valid == 1 && *errors == 1, "toy minimum/covariance invalid for " + *slice);
      Require(std::abs(*rout - (16 + *mt)) < .35 && std::abs(*rside - (20.25 + *mt)) < .35
                  && std::abs(*rlong - (25 + *mt)) < .35,
              "integer-count Gaussian toy radius recovery failed");
      Require(*ndf == 64 - (full ? 7 : 4), "external lambda/alpha wrongly counted as free parameters");
      Require(table ? std::isnan(*lambda) : std::isfinite(*lambda), "table input fabricated a scalar lambda");
      const std::string dir = "fits/" + *slice + "/";
      auto *model = dynamic_cast<TH3D *>(file.Get((dir + "ModelCF3D").c_str()));
      auto *core = dynamic_cast<TH3D *>(file.Get((dir + "CATS3DCoreCF3D").c_str()));
      Require(model && core && file.Get((dir + "LambdaMap3D").c_str()), "evaluated CATS3D artifacts missing");
      Require(!file.Get((dir + "LevyFit3D").c_str()), "CATS3D persisted a misleading Bowler callback");
      auto *graph = dynamic_cast<TGraph *>(file.Get((dir + "Fit_ProjX").c_str()));
      Require(graph && graph->GetN() == 4, "projection does not use actual computed bins");
      const double qwindow = .04;
      for (int i = 0; i < graph->GetN(); ++i) {
        double qo = 0, value = 0;
        graph->GetPoint(i, qo, value);
        double sum = 0;
        int count = 0;
        for (int y = 1; y <= 4; ++y) {
          if (std::abs(model->GetYaxis()->GetBinCenter(y)) <= qwindow) {
            for (int z = 1; z <= 4; ++z) {
              if (std::abs(model->GetZaxis()->GetBinCenter(z)) <= qwindow) {
                sum += model->GetBinContent(i + 1, y, z);
                ++count;
              }
            }
          }
        }
        Require(std::abs(value - sum / count) < 1e-10, "projection curve differs from persisted CATS3D model");
      }
    }
  }
}  // namespace
int main() {
  try {
    gROOT->SetBatch(true);
    Logger logger(LogLevel::kWarn);
    const auto root = std::filesystem::path(EXP_FEMTO_3D_SOURCE_DIR) / "build/tests/cats3d_workflow";
    std::filesystem::create_directories(root);
    const auto source = std::filesystem::path(EXP_FEMTO_3D_SOURCE_DIR) / "config/examples/cats3d_gaussian.toml";
    auto cfg = LoadApplicationConfig(source.string());
    Require(
        cfg.fit.correlation_backend == CorrelationBackend::kCats3D && cfg.fit.options.parameters.alpha.fixed_value == 2,
        "CATS3D example config did not parse");
    const auto original = Bytes(source.string());
    const auto invalid_path = (root / "invalid.toml").string();
    const auto invalid = [&](std::string text) {
      std::ofstream f(invalid_path);
      f << text;
      f.close();
      Reject(
          [&] {
            (void)LoadApplicationConfig(invalid_path);
          },
          "invalid CATS3D configuration accepted");
    };
    auto bad = original;
    bad.replace(bad.find("meaning = \"final_lambda\""),
                std::string("meaning = \"final_lambda\"").size(),
                "meaning = \"pair_purity\"");
    invalid(bad);
    bad = original;
    bad.insert(bad.find("correlation_backend"), "lambda_mode = \"shared_phi\"\n");
    invalid(bad);
    bad = original;
    bad.insert(bad.find("correlation_backend"), "coulomb_mode = \"none\"\n");
    invalid(bad);
    bad = original;
    bad += "\n[fit.parameters.alpha]\nfixed_value=1.4\n";
    invalid(bad);
    bad = original;
    bad.replace(bad.find("max_peak_gib = 8.0"), std::string("max_peak_gib = 8.0").size(), "max_peak_gib = 0");
    invalid(bad);
    bad = original;
    bad.insert(bad.find("with_coulomb"), "kernel_cache_file=123\n");
    invalid(bad);
    bad = original;
    bad.replace(bad.find("use_disk_cache = true"),
                std::string("use_disk_cache = true").size(),
                "use_disk_cache = false\nkernel_cache_file=\"native.root\"");
    invalid(bad);
    auto portable_config = original;
    portable_config.insert(portable_config.find("with_coulomb"), "kernel_cache_file=\"../native.root\"\n");
    {
      std::ofstream f(invalid_path);
      f << portable_config;
    }
    const auto parsed_portable = LoadApplicationConfig(invalid_path);
    Require(parsed_portable.fit.cats3d.kernel_cache_file == (root / "../native.root").lexically_normal().string(),
            "native kernel input path did not resolve relative to its configuration");
    auto production_example = LoadApplicationConfig(
        (std::filesystem::path(EXP_FEMTO_3D_SOURCE_DIR) / "config/oo_build_and_fit_6bins_cats3d.toml").string());
    (void)production_example;
    cfg.input.input_root = (root / "raw_counts.root").string();
    cfg.input.task_name = "task";
    cfg.input.same_event_subtask = "Same";
    cfg.input.mixed_event_subtask = "Mixed";
    cfg.input.sparse_object_name = "sparse";
    cfg.output.output_directory = root.string();
    cfg.output.fit_report_directory = root.string();
    cfg.output.cf_root_name = "cf.root";
    cfg.output.fit_root_name = "nominal.root";
    cfg.output.fit_summary_name = "nominal.tsv";
    cfg.output.fit_report_root_name = "report.root";
    cfg.centrality_bins = {{0, 10, ""}};
    cfg.fit_centrality_bins = cfg.centrality_bins;
    cfg.mt_bins = {{.2, .4, ""}, {.4, .6, ""}};
    cfg.fit_mt_bins = cfg.mt_bins;
    cfg.phi_bins = {{0, kPi / 3, ""}, {kPi / 3, 2 * kPi / 3, ""}, {2 * kPi / 3, kPi, ""}};
    cfg.build.mt_rebin = {};
    cfg.build.phi_rebin = {};
    cfg.build.phi_rebin.configured = true;
    cfg.build.phi_rebin.enabled = true;
    cfg.build.phi_rebin.mode = RebinMode::kRanges;
    cfg.build.split_mixed_event_by_phi = true;
    cfg.build.progress = cfg.fit.progress = ProgressMode::kDisabled;
    cfg.fit.model = FitModel::kDiag;
    cfg.fit.options.use_q2_baseline = false;
    cfg.fit.options.fit_q_max = .12;
    auto &cats = cfg.fit.cats3d;
    cats.with_coulomb = false;
    cats.grid_bins = 31;
    cats.grid_min_fm = -30;
    cats.grid_max_fm = 30;
    cats.use_disk_cache = false;
    cats.max_kernel_gib = .1;
    if (!Cats3DModel::Available()) {
      Reject(
          [&] {
            (void)RunFit(cfg, logger);
          },
          "no-CATS3D runtime silently fell back to Bowler");
      std::cout << "CATS3D configuration validation and disabled backend rejection passed\n";
      return 0;
    }
    // Exercise nominal, table-lambda and scan paths with the file's shell convention.
    cfg.build.normalization.mode = CfNormalizationMode::kBoxShell;
    for (auto &range : cfg.build.normalization.ranges) range = {.06, .12};
    WriteToy(cfg.input.input_root);
    ValidateApplicationConfig(cfg);
    (void)RunBuildCf(cfg, logger);
    const auto lambda_file = (root / "lambda.root").string();
    {
      TFile f(lambda_file.c_str(), "RECREATE");
      TH1D a("lambda_0", "", 4, 0, .2), b("lambda_1", "", 4, 0, .2);
      for (int i = 1; i <= 4; ++i) {
        a.SetBinContent(i, .7);
        a.SetBinError(i, .01);
        b.SetBinContent(i, .55);
        b.SetBinError(i, .02);
      }
      a.Write();
      b.Write();
    }
    Cats3DLambdaInput input;
    input.kind = "histogram";
    input.file = lambda_file;
    input.histogram = "lambda_0";
    input.axis_unit = "GeV";
    cats.lambda = input;
    {
      TFile f((root / "cf.root").c_str(), "READ");
      auto *tree = dynamic_cast<TTree *>(f.Get("meta/SliceCatalog"));
      TTreeReader r(tree);
      TTreeReaderValue<std::string> group(r, "group_id");
      TTreeReaderValue<int> mt(r, "mt_index");
      while (r.Next()) {
        if (*mt == 1) {
          auto override = input;
          override.group_id = *group;
          override.histogram = "lambda_1";
          cats.lambda_groups = {override};
          break;
        }
      }
    }
    cats.lambda_scan.enabled = true;
    cats.lambda_scan.constant_values = {.65, .7, 1.2};
    cats.lambda_scan.scale_values = {.95, 1, 2};
    cats.lambda_scan.root_name = "scan.root";
    ValidateApplicationConfig(cfg);
    const auto stats = RunFit(cfg, logger);
    Require(stats.fitted_slices == 8, "not all toy slices fitted");
    CheckFit((root / "nominal.root").string(), 8, true, false);
    {
      TFile f((root / "scan.root").c_str(), "READ");
      TFile cf((root / "cf.root").c_str(), "READ");
      auto *catalog = dynamic_cast<TTree *>(cf.Get("meta/SliceCatalog"));
      TTreeReader r(catalog);
      TTreeReaderValue<std::string> group(r, "group_id"), slice(r, "slice_id");
      TTreeReaderValue<int> integrated(r, "is_phi_integrated");
      std::map<std::string, std::vector<std::string>> expected;
      while (r.Next()) {
        if (!*integrated) {
          expected[*group].push_back(*slice);
        }
      }
      for (const auto &[id, ids] : expected) {
        for (const char *kind : {"constant", "scale"}) {
          const std::string base = "lambda_scan/" + std::string(kind) + "/" + id + "/";
          auto *tree = dynamic_cast<TTree *>(f.Get((base + "GroupPoints").c_str()));
          Require(tree && tree->GetEntries() == 3, "group scan points missing");
          TTreeReader g(tree);
          TTreeReaderValue<double> coord(g, "coordinate"), objective(g, "objective");
          TTreeReaderValue<bool> valid(g, "valid");
          TTreeReaderValue<std::vector<std::string>> members(g, "member_slice_ids");
          while (g.Next()) {
            Require(members->size() == 3, "integrated phi was included in group sum");
            if (*coord > 1.1) {
              Require(!*valid && std::isnan(*objective), "invalid lambda coordinate was clipped or omitted");
              continue;
            }
            double sum = 0;
            for (const auto &sid : ids) {
              auto *t = dynamic_cast<TTree *>(f.Get((base + sid + "/ScanPoints").c_str()));
              Require(t && t->GetEntries() == 3, "slice scan points missing");
              auto *attempts = dynamic_cast<TTree *>(f.Get((base + sid + "/AttemptPoints").c_str()));
              Require(attempts && attempts->GetEntries() > 3, "scan retries missing");
              TTreeReader q(t);
              TTreeReaderValue<double> x(q, "coordinate"), o(q, "objective");
              TTreeReaderValue<bool> v(q, "valid");
              while (q.Next()) {
                if (*x == *coord) {
                  Require(*v, "valid toy scan failed");
                  sum += *o;
                }
              }
            }
            Require(*valid && std::abs(*objective - sum) < 1e-8, "group curve is not the complete phi objective sum");
          }
        }
      }
    }
    const auto formal = Bytes((root / "nominal.root").string()), tsv = Bytes((root / "nominal.tsv").string()),
               report = Bytes((root / "report.root").string());
    cats.lambda_scan.execution_mode = "scan_only";
    cats.lambda_scan.constant_values = {.7};
    cats.lambda_scan.scale_values.clear();
    (void)RunFit(cfg, logger);
    Require(formal == Bytes((root / "nominal.root").string()) && tsv == Bytes((root / "nominal.tsv").string())
                && report == Bytes((root / "report.root").string()),
            "scan_only changed formal outputs");
    cats.lambda_scan.root_name = cfg.output.fit_root_name;
    Reject(
        [&] {
          (void)RunFit(cfg, logger);
        },
        "scan-only output/input collision accepted");
    Require(formal == Bytes((root / "nominal.root").string()), "failed scan preflight truncated production output");
    cats.lambda_scan.root_name = "scan.root";
    cats.kernel_cache_file = (root / "nominal.root").string();
    Reject(
        [&] {
          (void)RunFit(cfg, logger);
        },
        "formal output could overwrite native kernel input");
    Require(formal == Bytes((root / "nominal.root").string()), "kernel path collision truncated input");
    cats.kernel_cache_file.clear();
    cats.lambda_scan = {};
    cats.lambda_groups.clear();
    cats.lambda = input;
    cfg.fit_mt_bins = {{.2, .4, ""}};
    cfg.fit.model = FitModel::kFull;
    cfg.output.fit_root_name = "full.root";
    cfg.output.fit_summary_name = "full.tsv";
    cfg.output.fit_report_root_name = "full_report.root";
    (void)RunFit(cfg, logger);
    CheckFit((root / "full.root").string(), 4, true, true);
    cfg.fit.model = FitModel::kDiag;
    cfg.fit.options.use_pml = false;
    cfg.output.fit_root_name = "chi2.root";
    cfg.output.fit_summary_name = "chi2.tsv";
    cfg.output.fit_report_root_name = "chi2_report.root";
    (void)RunFit(cfg, logger);
    CheckFit((root / "chi2.root").string(), 4, true, false);
    cfg.fit.options.use_pml = true;
    cats.pml_max_calls = 1;
    cats.lambda_scan.enabled = true;
    cats.lambda_scan.constant_values = {.7};
    cats.lambda_scan.root_name = "nonconverged_scan.root";
    cfg.output.fit_root_name = "nonconverged.root";
    cfg.output.fit_summary_name = "nonconverged.tsv";
    cfg.output.fit_report_root_name = "nonconverged_report.root";
    const auto failed = RunFit(cfg, logger);
    Require(failed.fitted_slices == 0, "explicit call limit did not retain nonconverged fits");
    {
      TFile f((root / "nonconverged.root").c_str(), "READ");
      auto *t = dynamic_cast<TTree *>(f.Get("meta/FitCatalog"));
      Require(t && t->GetEntries() == 4, "nonconverged fit rows disappeared");
      TTreeReader r(t);
      TTreeReaderValue<int> valid(r, "point_estimate_valid"), errors(r, "parameter_errors_valid");
      TTreeReaderValue<double> objective(r, "fit_statistic"), error(r, "rout2_err");
      while (r.Next()) {
        Require(!*valid && !*errors && std::isnan(*objective) && std::isnan(*error),
                "failed fit acquired a fake objective/error");
      }
    }
    std::cout << "CATS3D PML/full/chi2 integer-count recovery, group lambda inputs, scans, "
                 "artifacts and scan-only isolation passed\n";
    std::cout << "ARTIFACTS: " << root << '\n';
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
