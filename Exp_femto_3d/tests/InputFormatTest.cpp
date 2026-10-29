#include <algorithm>
#include <cmath>
#include <filesystem>
#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TNamed.h"
#include "TPad.h"
#include "TTree.h"
#include "exp_femto_3d/Config.h"
#include "exp_femto_3d/Workflow.h"

namespace {
  using namespace exp_femto_3d;
  namespace fs = std::filesystem;
  const Logger logger(LogLevel::kError);
  constexpr const char *object = "hQoutVsQsideVsQlongVsMtVsCentVsEventPlaneAngleVsQvector";

  void Expect(bool ok, const std::string &message) {
    if (!ok) throw std::runtime_error(message);
  }

  void Near(double a, double b, const std::string &message) {
    if (std::isnan(a) && std::isnan(b)) return;
    Expect(std::isfinite(a) && std::isfinite(b) &&
               std::abs(a - b) <= 1.e-8 * std::max({1., std::abs(a), std::abs(b)}),
           message + ": " + std::to_string(a) + " vs " + std::to_string(b));
  }

  // Identical weighted fills in two producer layouts. Splitting each weight into
  // two fills makes stored sum(w^2) observably different from a Poisson fallback.
  void WriteInput(const fs::path &path, bool femto, int fault = 0) {
    TFile out(path.c_str(), "RECREATE");
    auto *task = out.mkdir(femto ? "femto-pair-track-track" : "task");
    auto *base = femto ? task->mkdir("TrackTrack") : task;
    const int phi_axis = femto ? 5 : 6, qn_axis = femto ? 6 : 5;
    for (bool mixed : {false, true}) {
      auto *dir = base->mkdir(femto ? (mixed ? "ME" : "SE") : (mixed ? "Mixed" : "Same"));
      if (femto) dir = dir->mkdir("Analysis");
      int bins[7] = {4, 4, 4, 4, 2, 2, 2};
      double low[7] = {-.2, -.2, -.2, .2, 0, 0, 0};
      double high[7] = {.2, .2, .2, 1.4, 20, 2, 2};
      bins[phi_axis] = fault == 3 ? 7 : 8;
      high[phi_axis] = fault == 5 ? kPi : 2 * kPi;
      if (fault == 1 && mixed) high[qn_axis] = 3;
      const int nd = fault == 2 ? 6 : 7;
      THnSparseF h(object, "producer counts", nd, bins, low, high);
      const double mt_edges[] = {.2, .3, .5, .8, 1.4};
      h.SetBinEdges(3, mt_edges);
      if (fault == 4) {
        const double edges[] = {0,           kPi / 4,     kPi / 2, 3 * kPi / 4, kPi, 5 * kPi / 4 + .1,
                                3 * kPi / 2, 7 * kPi / 4, 2 * kPi};
        h.SetBinEdges(phi_axis, edges);
      }
      h.Sumw2();
      auto fill = [&](double qo, double qs, double ql, double mt, double cent, double phi, double qn,
                      double w) {
        double v[7] = {qo, qs, ql, mt, cent, 0, 0};
        v[phi_axis] = phi;
        v[qn_axis] = qn;
        h.Fill(v, .4 * w);
        h.Fill(v, .6 * w);
      };
      if (nd == 7) {
        for (int im = 0; im < 4; ++im)
          for (int ic = 0; ic < 2; ++ic)
            for (int ip = 0; ip < 8; ++ip)
              for (int iq = 0; iq < 2; ++iq)
                for (int x = 0; x < 4; ++x)
                  for (int y = 0; y < 4; ++y)
                    for (int z = 0; z < 4; ++z) {
                      const double qo = -.15 + .1 * x, qs = -.15 + .1 * y, ql = -.15 + .1 * z;
                      const double correlation =
                          1 + (.3 + .05 * ip) *
                                  std::exp(-20 * (qo * qo + qs * qs + ql * ql) / (.1973269804 * .1973269804));
                      const double w = (100 + 10 * im + 3 * ic + 2 * ip + iq) * (mixed ? 1 : correlation);
                      fill(qo, qs, ql, .5 * (mt_edges[im] + mt_edges[im + 1]), 5 + 10 * ic,
                           (ip + .5) * kPi / 4, iq + .5, w);
                    }
        // Auxiliary flows must not leak into phi_all or qn_all normal-bin selections.
        fill(.05, .05, .05, .25, 5, -.1, .5, 1.e6);
        fill(.05, .05, .05, .25, 5, 2 * kPi + .1, .5, 1.e6);
        fill(.05, .05, .05, .25, 5, .1, 3.5, 1.e6);
        fill(.05, .05, .05, 2, 5, .1, .5, 1.e6);
      }
      dir->cd();
      h.Write();
      TNamed wrong("wrong_type", "not a histogram");
      wrong.Write();
    }
  }

  ApplicationConfig Config(const fs::path &base, const fs::path &input, bool femto, const std::string &name) {
    ApplicationConfig c;
    c.input.input_root = input.string();
    c.input.format = femto ? InputFormat::kFemto : InputFormat::kFemtoDream;
    c.input.task_name = femto ? "femto-pair-track-track" : "task";
    c.input.same_event_subtask = femto ? "TrackTrack/SE/Analysis" : "Same";
    c.input.mixed_event_subtask = femto ? "TrackTrack/ME/Analysis" : "Mixed";
    c.input.sparse_object_name = object;
    c.output.output_directory = base.string();
    c.output.cf_root_name = name + ".root";
    c.output.fit_root_name = name + "_fit.root";
    c.output.fit_summary_name = name + "_fit.tsv";
    c.output.fit_report_root_name = name + "_report.root";
    c.build.reopen_output_file_per_slice = false;
    c.build.progress = ProgressMode::kDisabled;
    c.fit.progress = ProgressMode::kDisabled;
    c.fit.model = FitModel::kDiag;
    c.build.mt_rebin.configured = true;
    c.build.phi_rebin.configured = true;
    c.centrality_bins = {{0, 10, ""}};
    ValidateApplicationConfig(c);
    return c;
  }

  void CompareEntry(const SliceCatalogEntry &a, const SliceCatalogEntry &b) {
    Expect(std::tie(a.slice_id, a.group_id, a.slice_directory, a.se_object_path, a.me_object_path,
                    a.cf_object_path, a.projection_x_path, a.projection_y_path, a.projection_z_path,
                    a.centrality_index, a.mt_index, a.qn_index, a.phi_index, a.qn_label,
                    a.build_uses_symmetric_phi_range, a.split_mixed_event_by_phi, a.split_mixed_event_by_qn,
                    a.mt_rebin_enabled, a.mt_rebin_mode, a.phi_rebin_enabled, a.phi_rebin_mode,
                    a.is_qn_integrated, a.is_phi_integrated, a.phi_mode) ==
               std::tie(b.slice_id, b.group_id, b.slice_directory, b.se_object_path, b.me_object_path,
                        b.cf_object_path, b.projection_x_path, b.projection_y_path, b.projection_z_path,
                        b.centrality_index, b.mt_index, b.qn_index, b.phi_index, b.qn_label,
                        b.build_uses_symmetric_phi_range, b.split_mixed_event_by_phi,
                        b.split_mixed_event_by_qn, b.mt_rebin_enabled, b.mt_rebin_mode, b.phi_rebin_enabled,
                        b.phi_rebin_mode, b.is_qn_integrated, b.is_phi_integrated, b.phi_mode),
           "catalog identity/policy mismatch");
    const double av[] = {a.cent_low,       a.cent_high,       a.mt_low,           a.mt_high,
                         a.qn_low,         a.qn_high,         a.raw_phi_low,      a.raw_phi_high,
                         a.raw_phi_center, a.display_phi_low, a.display_phi_high, a.display_phi_center};
    const double bv[] = {b.cent_low,       b.cent_high,       b.mt_low,           b.mt_high,
                         b.qn_low,         b.qn_high,         b.raw_phi_low,      b.raw_phi_high,
                         b.raw_phi_center, b.display_phi_low, b.display_phi_high, b.display_phi_center};
    for (int i = 0; i < 12; ++i) Near(av[i], bv[i], "catalog physical coordinates mismatch");
    const auto &n = a.normalization, &m = b.normalization;
    Expect(n.present == m.present && n.config.mode == m.config.mode && n.config.ranges == m.config.ranges &&
               n.outer_first == m.outer_first && n.outer_last == m.outer_last &&
               n.inner_first == m.inner_first && n.inner_last == m.inner_last,
           "normalization metadata mismatch");
    Near(n.se_integral, m.se_integral, "SE normalization mismatch");
    Near(n.me_integral, m.me_integral, "ME normalization mismatch");
    Near(n.cf_scale, m.cf_scale, "CF scale mismatch");
  }

  void CompareOutputs(const fs::path &old_path, const fs::path &new_path) {
    const auto a = LoadSliceCatalog(old_path.string()), b = LoadSliceCatalog(new_path.string());
    Expect(!a.empty() && a.size() == b.size(), "catalog size mismatch");
    TFile old_file(old_path.c_str(), "READ"), new_file(new_path.c_str(), "READ");
    bool second_half = false;
    for (std::size_t i = 0; i < a.size(); ++i) {
      CompareEntry(a[i], b[i]);
      second_half = second_half || (!b[i].is_phi_integrated && b[i].raw_phi_center > kPi);
      for (const auto &path : {a[i].se_object_path, a[i].me_object_path, a[i].cf_object_path,
                               a[i].projection_x_path, a[i].projection_y_path, a[i].projection_z_path}) {
        auto *s = dynamic_cast<TH1 *>(old_file.Get(path.c_str()));
        auto *m = dynamic_cast<TH1 *>(new_file.Get(path.c_str()));
        Expect(s && m && s->GetNcells() == m->GetNcells(), "missing output or histogram size mismatch");
        Expect(s->GetSumw2N() > 0 && m->GetSumw2N() > 0, "stored errors lost");
        for (int bin = 0; bin < s->GetNcells(); ++bin) {
          Near(s->GetBinContent(bin), m->GetBinContent(bin), "raw/CF/projection content mismatch: " + path);
          Near(s->GetBinError(bin), m->GetBinError(bin), "raw/CF/projection error mismatch: " + path);
        }
      }
      if (b[i].is_phi_integrated) Near(b[i].raw_phi_high, 2 * kPi, "phi_all does not cover full input");
    }
    Expect(second_half, "phi > pi was lost or relabeled");
  }

  // Compare against independently built, unmerged half-period histograms. The
  // oracle sums raw counts/variances, not already normalized CFs or fitted radii.
  void TestFoldedCounts(const fs::path &base, const fs::path &input) {
    ApplicationConfig native_fold;
    bool differs_from_cf_average = false;
    for (int scenario = 0; scenario < 3; ++scenario) {
      auto separate = Config(base, input, true, "fold_reference_" + std::to_string(scenario));
      separate.build.split_mixed_event_by_phi = true;
      separate.build.split_same_event_by_qn = scenario > 0;
      separate.build.split_mixed_event_by_qn = scenario > 0;
      if (scenario > 0) separate.qn_bins = {{0, 1, "low"}, {1, 2, "high"}};
      if (scenario == 1) {
        for (auto *r : {&separate.build.mt_rebin, &separate.build.phi_rebin}) {
          r->enabled = true;
          r->mode = RebinMode::kFactor;
          r->factor = 2;
        }
        separate.build.normalization.mode = CfNormalizationMode::kBoxShell;
        for (auto &range : separate.build.normalization.ranges) range = {.1, .2};
      }
      if (scenario == 2) {
        separate.build.mt_rebin.enabled = separate.build.phi_rebin.enabled = true;
        separate.build.mt_rebin.mode = separate.build.phi_rebin.mode = RebinMode::kRanges;
        separate.mt_bins = {{.2, .5, ""}, {.5, 1.4, ""}};
        separate.phi_bins = {
            {0, kPi / 2, ""}, {kPi / 2, kPi, ""}, {kPi, 3 * kPi / 2, ""}, {3 * kPi / 2, 2 * kPi, ""}};
      }
      auto folded = separate;
      folded.build.phi_mode = PhiMode::kFoldPi;
      folded.output.cf_root_name = "folded_" + std::to_string(scenario) + ".root";
      if (scenario == 2) folded.phi_bins = {{0, kPi / 2, ""}, {kPi / 2, kPi, ""}};
      ValidateApplicationConfig(separate);
      ValidateApplicationConfig(folded);
      (void)RunBuildCf(separate, logger);
      const auto stats = RunBuildCf(folded, logger);
      Expect(stats.phi_input_bins == 8 && stats.phi_output_bins == (scenario == 0 ? 4U : 2U),
             "folded native/rebinned angle count is wrong");
      const auto ref = LoadSliceCatalog((base / separate.output.cf_root_name).string());
      const auto out = LoadSliceCatalog((base / folded.output.cf_root_name).string());
      TFile original((base / separate.output.cf_root_name).c_str(), "READ");
      TFile merged((base / folded.output.cf_root_name).c_str(), "READ");
      for (const auto &entry : out) {
        Expect(entry.phi_mode == "fold_pi" && entry.raw_phi_high <= kPi + 1.e-9, "folded metadata lost");
        const auto find = [&](double shift) -> const SliceCatalogEntry & {
          auto it = std::find_if(ref.begin(), ref.end(), [&](const auto &r) {
            return r.group_id == entry.group_id && r.is_phi_integrated == entry.is_phi_integrated &&
                   (entry.is_phi_integrated ||
                    (std::abs(r.raw_phi_low - entry.raw_phi_low - shift) < 1.e-9 &&
                     std::abs(r.raw_phi_high - entry.raw_phi_high - shift) < 1.e-9));
          });
          Expect(it != ref.end(), "missing unmerged reference slice");
          return *it;
        };
        const auto &first = find(0), &second = find(kPi);
        const double se_norm = first.normalization.se_integral +
                               (entry.is_phi_integrated ? 0 : second.normalization.se_integral);
        const double me_norm = first.normalization.me_integral +
                               (entry.is_phi_integrated ? 0 : second.normalization.me_integral);
        Near(entry.normalization.cf_scale, me_norm / se_norm,
             "folded normalization was averaged instead of recomputed");
        auto *s = dynamic_cast<TH3D *>(merged.Get(entry.se_object_path.c_str()));
        auto *b = dynamic_cast<TH3D *>(merged.Get(entry.me_object_path.c_str()));
        auto *cf = dynamic_cast<TH3D *>(merged.Get(entry.cf_object_path.c_str()));
        auto *s1 = dynamic_cast<TH3D *>(original.Get(first.se_object_path.c_str()));
        auto *b1 = dynamic_cast<TH3D *>(original.Get(first.me_object_path.c_str()));
        auto *c1 = dynamic_cast<TH3D *>(original.Get(first.cf_object_path.c_str()));
        auto *s2 = dynamic_cast<TH3D *>(original.Get(second.se_object_path.c_str()));
        auto *b2 = dynamic_cast<TH3D *>(original.Get(second.me_object_path.c_str()));
        auto *c2 = dynamic_cast<TH3D *>(original.Get(second.cf_object_path.c_str()));
        Expect(s && b && cf && s1 && b1 && c1 && s2 && b2 && c2, "missing folding oracle histogram");
        for (int x = 1; x <= 4; ++x)
          for (int y = 1; y <= 4; ++y)
            for (int z = 1; z <= 4; ++z) {
              const int bin = s->GetBin(x, y, z);
              const double sw =
                  s1->GetBinContent(bin) + (entry.is_phi_integrated ? 0 : s2->GetBinContent(bin));
              const double bw =
                  b1->GetBinContent(bin) + (entry.is_phi_integrated ? 0 : b2->GetBinContent(bin));
              const double sv = std::pow(s1->GetBinError(bin), 2) +
                                (entry.is_phi_integrated ? 0 : std::pow(s2->GetBinError(bin), 2));
              const double bv = std::pow(b1->GetBinError(bin), 2) +
                                (entry.is_phi_integrated ? 0 : std::pow(b2->GetBinError(bin), 2));
              Near(s->GetBinContent(bin), sw, "folded SE is not the raw sum");
              Near(b->GetBinContent(bin), bw, "folded ME is not the raw sum");
              Near(std::pow(s->GetBinError(bin), 2), sv, "SE weighted variances did not add");
              Near(std::pow(b->GetBinError(bin), 2), bv, "ME weighted variances did not add");
              const double expected = entry.normalization.cf_scale * sw / bw;
              Near(cf->GetBinContent(bin), expected, "folded CF was not normalized after summing counts");
              Near(cf->GetBinError(bin),
                   entry.normalization.cf_scale * std::sqrt(sv / (bw * bw) + sw * sw * bv / std::pow(bw, 4)),
                   "folded CF error propagation changed");
              if (!entry.is_phi_integrated &&
                  std::abs(expected - .5 * (c1->GetBinContent(bin) + c2->GetBinContent(bin))) > 1.e-6)
                differs_from_cf_average = true;
            }
      }
      if (scenario == 0) native_fold = folded;
    }
    Expect(differs_from_cf_average, "folding oracle cannot distinguish wrong averaging of CFs");
    native_fold.fit_mt_bins = {{.2, .3, ""}};
    native_fold.output.fit_root_name = "folded_fit.root";
    native_fold.output.fit_report_root_name = "folded_report.root";
    native_fold.output.fit_summary_name = "folded_fit.tsv";
    const auto fit = RunFit(native_fold, logger);
    Expect(fit.selected_slices == 5, "folded fit should read four phi bins plus phi_all");
    TFile result((base / native_fold.output.fit_root_name).c_str(), "READ");
    const auto catalog = LoadSliceCatalog((base / native_fold.output.cf_root_name).string());
    auto *harmonic = dynamic_cast<TF1 *>(
        result.Get(("summary/R2_vs_phi/" + catalog.front().group_id + "/Rout2_phi_fit").c_str()));
    Expect(harmonic, "folded R2 summary did not retain the 0..pi interface");
    Near(harmonic->GetXmin(), 0, "folded summary low edge");
    Near(harmonic->GetXmax(), kPi, "folded summary high edge");
  }

  // Saved CATS summaries must apply completeness separately to each half, so
  // an invalid first-half error cannot suppress a valid second-half harmonic.
  void TestCatsPeriodReports(const fs::path &base, const std::vector<SliceCatalogEntry> &catalog) {
    std::vector<LevyFitResult> results;
    for (const auto &e : catalog) {
      if (e.mt_low != .2 || e.is_phi_integrated) continue;
      LevyFitResult r;
      r.slice_id = e.slice_id;
      r.group_id = e.group_id;
      r.centrality_index = e.centrality_index;
      r.mt_index = e.mt_index;
      r.mt_low = e.mt_low;
      r.mt_high = e.mt_high;
      r.cent_low = e.cent_low;
      r.cent_high = e.cent_high;
      r.qn_index = e.qn_index;
      r.qn_low = e.qn_low;
      r.qn_high = e.qn_high;
      r.qn_label = e.qn_label;
      r.phi_index = e.phi_index;
      r.phi = e.raw_phi_center;
      r.raw_phi_low = e.raw_phi_low;
      r.raw_phi_high = e.raw_phi_high;
      r.display_phi_low = e.display_phi_low;
      r.display_phi_high = e.display_phi_high;
      r.point_estimate_valid = true;
      r.parameter_errors_valid = true;
      r.correlation_backend = "cats3d";
      r.lambda_is_external_fixed = true;
      r.lambda_input_kind = "constant";
      r.lambda_input_meaning = "final_lambda";
      r.lambda_err = 0;
      r.alpha_err = 0;
      r.rout2 = 20 + 2 * std::cos(2 * r.phi);
      r.rside2 = 18 + std::cos(2 * r.phi);
      r.rlong2 = 22;
      r.rout2_err = r.rside2_err = r.rlong2_err = .1;
      r.alpha = 2;
      r.lambda = .5;
      results.push_back(r);
    }
    Expect(results.size() == 8, "synthetic CATS half-period report needs eight input members");
    for (bool incomplete : {false, true}) {
      const auto path = base / (incomplete ? "cats_halves_incomplete.root" : "cats_halves.root");
      if (fs::exists(path)) fs::remove(path);
      if (incomplete) results.front().parameter_errors_valid = false;
      WriteCats3DFitReportRoot(path.string(), results);
      TFile file(path.c_str(), "READ");
      const std::string tail = "/summary/R2_vs_phi/" + results.front().group_id + "/";
      auto *first = dynamic_cast<TF1 *>(file.Get(("phi_0_pi" + tail + "Rout2_phi_fit").c_str()));
      auto *second = dynamic_cast<TF1 *>(file.Get(("phi_pi_2pi" + tail + "Rout2_phi_fit").c_str()));
      Expect(second, "valid second-half CATS summary was lost");
      Near(second->GetXmin(), kPi, "second-half CATS fit low limit");
      Near(second->GetXmax(), 2 * kPi, "second-half CATS fit high limit");
      Expect(incomplete ? first == nullptr : first != nullptr, "CATS completeness was not checked per half");
      if (incomplete)
        Expect(file.Get(("phi_0_pi" + tail + "SummaryStatus").c_str()), "missing incomplete-half status");
    }
  }

  void Sentinel(const fs::path &path) {
    TFile file(path.c_str(), "RECREATE");
    TNamed mark("sentinel", "preserve");
    mark.Write();
  }

  void RejectPreserving(const fs::path &path, const std::function<void()> &run,
                        const std::string &diagnostic) {
    Sentinel(path);
    bool rejected = false;
    try {
      run();
    } catch (const std::exception &e) {
      rejected = true;
      Expect(std::string(e.what()).find(diagnostic) != std::string::npos,
             "unexpected diagnostic: " + std::string(e.what()));
    }
    Expect(rejected, "invalid input was accepted");
    TFile file(path.c_str(), "READ");
    Expect(file.Get("sentinel"), "validation truncated output");
  }
}  // namespace

int main() {
  try {
    const fs::path base = fs::temp_directory_path() / "exp_femto_3d_input_format_test";
    fs::create_directories(base);
    const auto old_input = base / "old.root", new_input = base / "new.root";
    WriteInput(old_input, false);
    WriteInput(new_input, true);
    for (int scenario = 0; scenario < 4; ++scenario) {
      auto a = Config(base, old_input, false, "old_" + std::to_string(scenario));
      auto b = Config(base, new_input, true, "new_" + std::to_string(scenario));
      for (auto *c : {&a, &b}) {
        c->build.split_same_event_by_qn = scenario > 0;
        c->build.split_mixed_event_by_phi = scenario == 1 || scenario == 3;
        c->build.split_mixed_event_by_qn = scenario == 1 || scenario == 2;
        if (scenario > 0) c->qn_bins = {{0, 1, "low"}, {1, 2, "high"}};
        if (scenario == 1) {
          for (auto *r : {&c->build.mt_rebin, &c->build.phi_rebin}) {
            r->enabled = true;
            r->mode = RebinMode::kFactor;
            r->factor = 2;
          }
          c->build.normalization.mode = CfNormalizationMode::kBoxShell;
          for (auto &range : c->build.normalization.ranges) range = {.1, .2};
        }
        if (scenario == 2) {
          c->build.mt_rebin.enabled = true;
          c->build.phi_rebin.enabled = true;
          c->build.mt_rebin.mode = RebinMode::kRanges;
          c->build.phi_rebin.mode = RebinMode::kRanges;
          c->mt_bins = {{.2, .5, ""}, {.5, 1.4, ""}};
          c->phi_bins = {{0, kPi, ""}, {kPi, 2 * kPi, ""}};
        }
        ValidateApplicationConfig(*c);
        const auto stats = RunBuildCf(*c, logger);
        Expect(stats.stored_slices > 0 && stats.phi_input_bins == 8, "layout resolved wrong phi axis");
        const std::size_t mt_groups = (scenario == 1 || scenario == 2) ? 2 : 4;
        const std::size_t phi_groups = scenario == 1 ? 4 : (scenario == 2 ? 2 : 8);
        Expect(stats.stored_slices == mt_groups * (scenario > 0 ? 3 : 1) * (phi_groups + 1),
               "missing layout slices");
      }
      CompareOutputs(base / a.output.cf_root_name, base / b.output.cf_root_name);
    }
    // Independent weighted-error oracle for one known native slice (qn integrated).
    auto c = Config(base, new_input, true, "new_0");
    auto entries = LoadSliceCatalog((base / c.output.cf_root_name).string());
    TFile file((base / c.output.cf_root_name).c_str(), "READ");
    bool checked_weighted_oracle = false;
    for (const auto &e : entries)
      if (e.mt_low == .2 && !e.is_phi_integrated && e.phi_index == 0) {
        auto *me = dynamic_cast<TH3D *>(file.Get(e.me_object_path.c_str()));
        Expect(me, "missing weighted oracle histogram");
        // Scenario 0 deliberately shares a phi-integrated ME denominator: 100..115.
        Near(me->GetBinContent(1, 1, 1), 1720, "qn/phi integration included flows or wrong axis");
        Near(std::pow(me->GetBinError(1, 1, 1), 2), .52 * 185240, "weighted variance was replaced");
        checked_weighted_oracle = true;
      }
    Expect(checked_weighted_oracle, "weighted oracle slice not found");
    file.Close();
    auto bad = c;
    bad.output.cf_root_name = "preserve.root";
    for (const auto &name : {"missing", "wrong_type"}) {
      bad.input.sparse_object_name = name;
      RejectPreserving(
          base / bad.output.cf_root_name, [&] { (void)RunBuildCf(bad, logger); }, "input.format=femto");
    }
    bad = c;
    bad.output.cf_root_name = "preserve.root";
    for (int fault : {1, 2}) {
      const auto input = base / ("bad_" + std::to_string(fault) + ".root");
      WriteInput(input, true, fault);
      bad.input.input_root = input.string();
      RejectPreserving(
          base / bad.output.cf_root_name, [&] { (void)RunBuildCf(bad, logger); },
          fault == 1 ? "q-vector" : "seven-dimensional");
    }
    bad = c;
    bad.output.cf_root_name = "preserve.root";
    bad.build.mt_rebin.enabled = true;
    bad.build.mt_rebin.mode = RebinMode::kRanges;
    bad.mt_bins = {{.21, .5, ""}};
    RejectPreserving(base / bad.output.cf_root_name, [&] { (void)RunBuildCf(bad, logger); }, "align");
    for (int fault : {3, 4, 5}) {
      const auto input = base / ("unfoldable_" + std::to_string(fault) + ".root");
      WriteInput(input, true, fault);
      bad = c;
      bad.input.input_root = input.string();
      bad.output.cf_root_name = "preserve.root";
      bad.build.phi_mode = PhiMode::kFoldPi;
      RejectPreserving(base / bad.output.cf_root_name, [&] { (void)RunBuildCf(bad, logger); }, "fold_pi");
    }
    bad = c;
    bad.output.cf_root_name = "preserve.root";
    bad.build.phi_rebin.enabled = true;
    bad.build.phi_rebin.mode = RebinMode::kFactor;
    bad.build.phi_rebin.factor = 8;
    RejectPreserving(base / bad.output.cf_root_name, [&] { (void)RunBuildCf(bad, logger); }, "straddles");
    bad = c;
    bad.output.cf_root_name = "preserve.root";
    bad.build.map_pair_phi_to_symmetric_range = true;
    RejectPreserving(base / bad.output.cf_root_name, [&] { (void)RunBuildCf(bad, logger); }, "map_pair_phi");
    bad = c;
    bad.fit.map_pair_phi_to_symmetric_range = true;
    RejectPreserving(base / bad.output.fit_root_name, [&] { (void)RunFit(bad, logger); }, "map_pair_phi");
#ifdef EXP_FEMTO_3D_HAS_CATS3D
    bad.fit.correlation_backend = CorrelationBackend::kCats3D;
    RejectPreserving(base / bad.output.fit_root_name, [&] { (void)RunFit(bad, logger); }, "map_pair_phi");
#endif
    // A small synthetic fit exercises the unchanged catalog interface over 0..2pi.
    c.fit_mt_bins = {{.2, .3, ""}};
    const auto fit = RunFit(c, logger);
    Expect(fit.selected_slices == 9, "fit did not read all eight native phi slices plus phi_all");
    TFile fitted((base / c.output.fit_root_name).c_str(), "READ");
    auto *tree = dynamic_cast<TTree *>(fitted.Get("meta/FitCatalog"));
    Expect(tree && tree->GetEntries() == 9, "fit catalog smoke failed");
    const auto group =
        std::find_if(entries.begin(), entries.end(), [](const auto &e) { return e.mt_low == .2; });
    Expect(group != entries.end(), "missing smoke-fit group");
    TFile report((base / c.output.fit_report_root_name).c_str(), "READ");
    for (const auto &period : {std::string("phi_0_pi"), std::string("phi_pi_2pi")}) {
      const double low = period == "phi_0_pi" ? 0 : kPi;
      auto *harmonic = dynamic_cast<TF1 *>(
          fitted.Get((period + "/summary/R2_vs_phi/" + group->group_id + "/Rout2_phi_fit").c_str()));
      auto *points = dynamic_cast<TGraph *>(
          fitted.Get((period + "/summary/R2_vs_phi/" + group->group_id + "/Rout2_vs_phi").c_str()));
      Expect(harmonic && points && points->GetN() == 4, "missing independent half-period harmonic");
      Near(harmonic->GetXmin(), low, "half-period summary low edge changed");
      Near(harmonic->GetXmax(), low + kPi, "half-period summary high edge changed");
      Near(points->GetXaxis()->GetXmin(), low, "half-period graph low limit");
      Near(points->GetXaxis()->GetXmax(), low + kPi, "half-period graph high limit");
      auto *overview = dynamic_cast<TCanvas *>(report.Get(
          (period + "/source_parameters/cent_0.00-10.00/mt_0.20-0.30/source_parameters_overview_canvas")
              .c_str()));
      Expect(overview && overview->GetPad(3), "missing half-period overview canvas");
      auto *trend = dynamic_cast<TF1 *>(
          overview->GetPad(3)->GetListOfPrimitives()->FindObject("Rout2_vs_phi_overview_fit"));
      Expect(trend, "missing independent half-period overview trend");
      Near(trend->GetXmin(), low, "overview half-period low edge changed");
      Near(trend->GetXmax(), low + kPi, "overview half-period high edge changed");
    }
    TestFoldedCounts(base, new_input);
    TestCatsPeriodReports(base, entries);
    std::cout << "INPUT_FORMAT_TEST_OK\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
