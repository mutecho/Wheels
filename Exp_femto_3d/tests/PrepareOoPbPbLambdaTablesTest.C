#include "../scripts/prepare_oo_pbpb_lambda_tables.C"

#include <functional>
#include "TH2D.h"
#include "TMemFile.h"

namespace {
  void Expect(bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error(message);
  }

  void Near(double actual, double expected, const std::string &message) {
    if (!std::isfinite(actual) || !std::isfinite(expected)
        || std::abs(actual - expected) > 2.e-12 * std::max(1.0, std::abs(expected))) {
      throw std::runtime_error(message + ": actual=" + std::to_string(actual)
                               + " expected=" + std::to_string(expected));
    }
  }

  void Rejects(const std::function<void()> &operation, const std::string &reason) {
    try {
      operation();
    } catch (const std::exception &error) {
      Expect(std::string(error.what()).find(reason) != std::string::npos,
             "Wrong failure reason: " + std::string(error.what()));
      return;
    }
    throw std::runtime_error("Expected rejection: " + reason);
  }

  void CheckSyntheticCases() {
    TDirectory::TContext context;
    TH1D p1("p1", "", 3, 0., .3), p2("p2", "", 3, 0., .3);
    TH1D w1("w1", "", 3, 0., .3), w2("w2", "", 3, 0., .3);
    for (TH1 *h : {&p1, &p2, &w1, &w2}) { h->SetDirectory(nullptr); h->Sumw2(); }
    for (int bin = 1; bin <= 2; ++bin) {
      p1.SetBinContent(bin, .6); p1.SetBinError(bin, .1);
      p2.SetBinContent(bin, .9); p2.SetBinError(bin, .2);
      w1.SetBinContent(bin, 1.); w1.SetBinError(bin, 1.);
      w2.SetBinContent(bin, 3.); w2.SetBinError(bin, std::sqrt(3.));
    }
    auto merge = [&]() { return oo_pbpb_lambda::Merge({&p1, &p2}, {&w1, &w2}, "synthetic", .15); };
    auto table = merge();
    Near(table.purity->GetBinContent(1), .825, "Known weighted purity");
    Near(table.purity->GetBinError(1), std::sqrt(.37) / 4., "Fixed-weight purity error");
    Near(table.weight->GetBinError(1), 2., "Preserved weight uncertainty");
    Expect(table.purity->GetNbinsX() == 2, "Only unsupported tail should be trimmed");
    Near(table.purity->GetXaxis()->GetXmax(), .2, "Trimmed high edge");
    // Flow bins cannot contribute to the merge or extend its support.
    w1.SetBinContent(0, 1.e6); w1.SetBinContent(4, 1.e6);
    Near(merge().purity->GetBinContent(1), .825, "Ignore flow weights");
    w1.SetBinContent(1, 0); w2.SetBinContent(1, 0);
    Rejects([&]() { merge(); }, "Zero total weight");
    w1.SetBinContent(1, 1); w2.SetBinContent(1, 3);
    w1.SetBinContent(2, 0); w2.SetBinContent(2, 0);
    w1.SetBinContent(3, 1);
    Rejects([&]() { merge(); }, "Zero total weight");
    w1.SetBinContent(3, 0);
    Rejects([&]() { merge(); }, "required fit k range");
    w1.SetBinContent(2, 1); w2.SetBinContent(2, 3);
    TH1D wrong_axis("wrong_axis", "", 3, 0., .6);
    wrong_axis.SetDirectory(nullptr);
    Rejects([&]() { oo_pbpb_lambda::Merge({&p1, &p2}, {&w1, &wrong_axis}, "bad", .15); }, "Incompatible k axes");
    p1.SetBinContent(1, 1.1);
    Rejects([&]() { merge(); }, "Invalid purity/weight");
    p1.SetBinContent(1, .6); w1.SetBinContent(1, -1.);
    Rejects([&]() { merge(); }, "Invalid purity/weight");
    TMemFile missing("missing.root", "RECREATE");
    Rejects([&]() { oo_pbpb_lambda::ReadHistogram(missing, "absent"); }, "Missing/invalid");
    TH2D two_d("two_d", "", 2, 0., 1., 2, 0., 1.);
    two_d.Write();
    Rejects([&]() { oo_pbpb_lambda::ReadHistogram(missing, "two_d"); }, "Missing/invalid");
    std::cout << "PASS synthetic weighting, error formula, tail trimming, flow exclusion, missing objects, axis mismatch, invalid values, zero weights, and coverage\n";
  }

  // Independent readback uses the stored source objects, not Merge(), and checks
  // every output bin against its defining weighted sums and unchanged k edges.
  void CheckPreparedFile(const char *source_directory, const char *output_path) {
    TFile prepared(output_path, "READ");
    Expect(!prepared.IsZombie(), "Prepared file is unreadable");
    auto *status = dynamic_cast<TObjString *>(prepared.Get("meta/PreparationStatus"));
    Expect(status && status->GetString() == "complete", "Missing completion metadata");
    Expect(prepared.GetDirectory("tables")->GetNkeys() == 18, "Expected nine purity and nine weight tables");
    const std::vector<std::vector<std::string>> members = {
        {"0.20to0.30", "0.30to0.40"}, {"0.40to0.50", "0.50to0.60"},
        {"0.60to0.70", "0.70to0.90", "0.90to1.10"}};
    const std::array<std::string, 3> targets = {"0.20to0.40", "0.40to0.60", "0.60to1.10"};
    int checked = 0;
    for (const std::string cent : {"0to10", "10to30", "30to50"}) {
      const std::string original_path = std::string(source_directory) + "/kStarDepLambda_cent" + cent + "_pion_PbPb.root";
      TFile original(original_path.c_str(), "READ");
      Expect(!original.IsZombie(), "Cannot read original file");
      const std::string saved_prefix = "sources/cent" + cent + "/";
      auto md5 = std::unique_ptr<TMD5>(TMD5::FileChecksum(original_path.c_str()));
      const auto *saved_md5 = dynamic_cast<TObjString *>(prepared.Get((saved_prefix + "SourceMD5").c_str()));
      const auto *saved_path = dynamic_cast<TObjString *>(prepared.Get((saved_prefix + "SourcePath").c_str()));
      Expect(md5 && saved_md5 && saved_md5->GetString() == md5->AsString(), "Source MD5 mismatch");
      Expect(saved_path && saved_path->GetString() == original_path.c_str(), "Source path mismatch");
      TIter next(original.GetListOfKeys());
      while (auto *key = static_cast<TKey *>(next())) {
        auto *a = dynamic_cast<TH1 *>(original.Get(key->GetName()));
        auto *b = dynamic_cast<TH1 *>(prepared.Get((saved_prefix + key->GetName()).c_str()));
        Expect(a && b && a->GetNcells() == b->GetNcells(), "Missing source copy");
        for (int bin = 0; bin < a->GetNcells(); ++bin) {
          Near(b->GetBinContent(bin), a->GetBinContent(bin), "Source content preservation");
          Near(b->GetBinError(bin), a->GetBinError(bin), "Source error preservation");
        }
      }
      for (std::size_t g = 0; g < targets.size(); ++g) {
        const std::string name = "tables/hPairPurity_mT" + targets[g] + "_cent" + cent;
        const auto *p = dynamic_cast<TH1 *>(prepared.Get(name.c_str()));
        const auto *w = dynamic_cast<TH1 *>(prepared.Get((name + "_weights").c_str()));
        Expect(p && w, "Missing merged tables");
        std::vector<const TH1 *> source_p, source_w;
        for (const auto &member : members[g]) {
          source_p.push_back(dynamic_cast<TH1 *>(original.Get(("hPurityKstar_mT" + member + "_cent" + cent).c_str())));
          source_w.push_back(dynamic_cast<TH1 *>(original.Get(("hWeightKstar_mT" + member + "_cent" + cent).c_str())));
          Expect(source_p.back() && source_w.back(), "Missing defining source table");
        }
        for (int bin = 1; bin <= source_p.front()->GetNbinsX(); ++bin) {
          double den = 0., numerator = 0., variance = 0., weight_variance = 0.;
          for (std::size_t s = 0; s < source_p.size(); ++s) {
            const double weight = source_w[s]->GetBinContent(bin);
            den += weight;
            numerator += weight * source_p[s]->GetBinContent(bin);
            variance += std::pow(weight * source_p[s]->GetBinError(bin), 2);
            weight_variance += std::pow(source_w[s]->GetBinError(bin), 2);
          }
          if (bin > p->GetNbinsX()) { Near(den, 0., "Only zero-weight tails may be removed"); continue; }
          Expect(den > 0, "Unsupported output bin");
          Near(p->GetBinContent(bin), numerator / den, "Weighted purity readback");
          Near(p->GetBinError(bin), std::sqrt(variance) / den, "Purity error readback");
          Near(w->GetBinContent(bin), den, "Weight readback");
          Near(w->GetBinError(bin), std::sqrt(weight_variance), "Weight error readback");
          Near(p->GetXaxis()->GetBinLowEdge(bin), source_p.front()->GetXaxis()->GetBinLowEdge(bin), "Unchanged k edge");
        }
        Near(p->GetXaxis()->GetXmax(), source_p.front()->GetXaxis()->GetBinUpEdge(p->GetNbinsX()), "Unchanged final edge");
        Expect(p->GetXaxis()->GetXmin() == 0 && p->GetXaxis()->GetXmax() >= .5 * std::sqrt(3.) * .15,
               "Conservative full fit-grid coverage");
        ++checked;
      }
    }
    Rejects([&]() { oo_pbpb_lambda::Prepare(source_directory, output_path); }, "already exists");
    std::cout << "PASS independent readback of " << checked << " tables, all source copies, provenance, coverage, and overwrite refusal\n";
  }
}  // namespace

void PrepareOoPbPbLambdaTablesTest(const char *source_directory, const char *output_path) {
  try {
    CheckSyntheticCases();
    CheckPreparedFile(source_directory, output_path);
    std::cout << "ALL_CHECKS_PASSED\n";
  } catch (const std::exception &error) {
    std::cerr << "[test failure] " << error.what() << '\n';
    gSystem->Exit(1);
  }
}
