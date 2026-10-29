#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TKey.h"
#include "TMD5.h"
#include "TObjString.h"
#include "TSystem.h"

namespace oo_pbpb_lambda {
  struct MtGroup {
    std::string label;
    std::vector<std::string> members;
  };

  const std::array<std::string, 3> centralities = {"0to10", "10to30", "30to50"};
  const std::array<MtGroup, 3> mt_groups = {{
      {"0.20to0.40", {"0.20to0.30", "0.30to0.40"}},
      {"0.40to0.60", {"0.40to0.50", "0.50to0.60"}},
      {"0.60to1.10", {"0.60to0.70", "0.70to0.90", "0.90to1.10"}},
  }};

  std::string HistogramName(const std::string &prefix, const std::string &mt, const std::string &cent) {
    return prefix + "_mT" + mt + "_cent" + cent;
  }

  const TH1 &ReadHistogram(TDirectory &directory, const std::string &name) {
    const auto *histogram = dynamic_cast<TH1 *>(directory.Get(name.c_str()));
    if (!histogram || histogram->GetDimension() != 1 || histogram->GetNbinsX() < 1) {
      throw std::runtime_error("Missing/invalid one-dimensional histogram: " + name);
    }
    return *histogram;
  }

  void RequireSameAxis(const TH1 &reference, const TH1 &histogram) {
    if (histogram.GetDimension() != 1 || reference.GetNbinsX() != histogram.GetNbinsX()) {
      throw std::runtime_error("Incompatible k axes: " + std::string(histogram.GetName()));
    }
    for (int bin = 1; bin <= reference.GetNbinsX() + 1; ++bin) {
      const double a = reference.GetXaxis()->GetBinLowEdge(bin);
      const double b = histogram.GetXaxis()->GetBinLowEdge(bin);
      if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a - b) > 1.e-12
          || (bin > 1 && a <= reference.GetXaxis()->GetBinLowEdge(bin - 1))) {
        throw std::runtime_error("Incompatible k axes: " + std::string(histogram.GetName()));
      }
    }
  }

  struct MergedTables {
    std::unique_ptr<TH1D> purity;
    std::unique_ptr<TH1D> weight;
  };

  /**
   * Merge complete mT intervals at unchanged k* bin edges (GeV/c). W is the
   * supplied PbPb pair weight, not an mT-bin width or an OO pair weight.
   * Source errors are treated as independent and W as fixed for purity errors.
   * Only normal bins contribute; unsupported trailing bins are removed.
   */
  MergedTables Merge(const std::vector<const TH1 *> &purities,
                     const std::vector<const TH1 *> &weights,
                     const std::string &name,
                     double required_k_max_gev) {
    if (purities.empty() || purities.size() != weights.size() || !purities.front()
        || purities.front()->GetDimension() != 1 || purities.front()->GetNbinsX() < 1
        || !std::isfinite(required_k_max_gev) || required_k_max_gev <= 0) {
      throw std::runtime_error("Invalid merge inputs or required k range.");
    }
    const TH1 &reference = *purities.front();
    const int bins = reference.GetNbinsX();
    std::vector<double> sum_w(bins), sum_wp(bins), sum_w2_error2(bins), sum_weight_error2(bins);
    for (std::size_t member = 0; member < purities.size(); ++member) {
      if (!purities[member] || !weights[member]) throw std::runtime_error("Null source histogram.");
      RequireSameAxis(reference, *purities[member]);
      RequireSameAxis(reference, *weights[member]);
      for (int bin = 1; bin <= bins; ++bin) {
        const double p = purities[member]->GetBinContent(bin);
        const double ep = purities[member]->GetBinError(bin);
        const double w = weights[member]->GetBinContent(bin);
        const double ew = weights[member]->GetBinError(bin);
        if (!std::isfinite(p) || p < 0 || p > 1 || !std::isfinite(ep) || ep < 0
            || !std::isfinite(w) || w < 0 || !std::isfinite(ew) || ew < 0) {
          throw std::runtime_error("Invalid purity/weight content or error: "
                                   + std::string(purities[member]->GetName()));
        }
        const int i = bin - 1;
        sum_w[i] += w;
        sum_wp[i] += w * p;
        sum_w2_error2[i] += w * w * ep * ep;
        sum_weight_error2[i] += ew * ew;
        if (!std::isfinite(sum_w[i]) || !std::isfinite(sum_wp[i])
            || !std::isfinite(sum_w2_error2[i]) || !std::isfinite(sum_weight_error2[i])) {
          throw std::runtime_error("Nonfinite weighted sum.");
        }
      }
    }
    int last = bins;
    while (last > 0 && sum_w[last - 1] == 0) --last;
    if (last == 0) throw std::runtime_error("No positive total weight.");
    for (int bin = 1; bin <= last; ++bin) {
      if (sum_w[bin - 1] <= 0) throw std::runtime_error("Zero total weight inside supported k range.");
    }
    if (reference.GetXaxis()->GetBinLowEdge(1) > 0
        || reference.GetXaxis()->GetBinUpEdge(last) < required_k_max_gev) {
      throw std::runtime_error("Weighted table does not cover the required fit k range.");
    }
    std::vector<double> edges(last + 1);
    for (int bin = 1; bin <= last + 1; ++bin) edges[bin - 1] = reference.GetXaxis()->GetBinLowEdge(bin);
    MergedTables result;
    result.purity = std::make_unique<TH1D>(name.c_str(), "PbPb pair purity; k* (GeV/c); Pair purity", last, edges.data());
    result.weight = std::make_unique<TH1D>((name + "_weights").c_str(),
                                         "Sum of PbPb input weights; k* (GeV/c); Sum of weights", last, edges.data());
    result.purity->SetDirectory(nullptr);
    result.weight->SetDirectory(nullptr);
    result.purity->Sumw2();
    result.weight->Sumw2();
    // The weak-decay survival factor is deliberately NOT applied here. The
    // consuming TOML declares pair_purity and applies (1 - 0.12)^2 exactly once.
    for (int bin = 1; bin <= last; ++bin) {
      const int i = bin - 1;
      const double p = sum_wp[i] / sum_w[i];
      const double ep = std::sqrt(sum_w2_error2[i]) / sum_w[i];
      if (!std::isfinite(p) || p < 0 || p > 1 || !std::isfinite(ep)) {
        throw std::runtime_error("Invalid merged purity/error.");
      }
      result.purity->SetBinContent(bin, p);
      result.purity->SetBinError(bin, ep);
      result.weight->SetBinContent(bin, sum_w[i]);
      result.weight->SetBinError(bin, std::sqrt(sum_weight_error2[i]));
    }
    return result;
  }

  void WriteText(TDirectory &directory, const std::string &key, const std::string &text) {
    TDirectory::TContext context(&directory);
    TObjString object(text.c_str());
    if (object.Write(key.c_str()) <= 0) throw std::runtime_error("Failed to write metadata: " + key);
  }

  /** Prepare all nine tables before creating a new, non-overwriting output. */
  void Prepare(const std::string &input_directory, const std::string &output_root) {
    namespace fs = std::filesystem;
    if (output_root.empty() || fs::exists(output_root)) {
      throw std::runtime_error("Output is empty or already exists: " + output_root);
    }
    TDirectory::TContext context;
    // Conservatively covers every q_i in [-0.15,0.15] GeV/c, including edges;
    // the actual OO bin-center queries span about 0.00433--0.12557 GeV/c.
    const double required_k_max = 0.5 * std::sqrt(3.0) * 0.15;
    std::vector<std::unique_ptr<TFile>> inputs;
    std::vector<std::string> source_paths, checksums;
    std::vector<MergedTables> merged;
    for (const auto &cent : centralities) {
      const auto source = fs::absolute(fs::path(input_directory)
                                      / ("kStarDepLambda_cent" + cent + "_pion_PbPb.root")).lexically_normal();
      auto file = std::unique_ptr<TFile>(TFile::Open(source.c_str(), "READ"));
      if (!file || file->IsZombie()) throw std::runtime_error("Cannot open source: " + source.string());
      if (file->TestBit(TFile::kRecovered)) throw std::runtime_error("Recovered source file: " + source.string());
      auto checksum = std::unique_ptr<TMD5>(TMD5::FileChecksum(source.c_str()));
      if (!checksum) throw std::runtime_error("Cannot checksum source: " + source.string());
      source_paths.push_back(source.string());
      checksums.emplace_back(checksum->AsString());
      for (const auto &group : mt_groups) {
        std::vector<const TH1 *> p, w;
        for (const auto &member : group.members) {
          p.push_back(&ReadHistogram(*file, HistogramName("hPurityKstar", member, cent)));
          w.push_back(&ReadHistogram(*file, HistogramName("hWeightKstar", member, cent)));
        }
        merged.push_back(Merge(p, w, HistogramName("hPairPurity", group.label, cent), required_k_max));
      }
      inputs.push_back(std::move(file));
    }
    const fs::path destination(output_root);
    if (!destination.parent_path().empty()) fs::create_directories(destination.parent_path());
    auto output = std::unique_ptr<TFile>(TFile::Open(output_root.c_str(), "CREATE"));
    if (!output || output->IsZombie()) throw std::runtime_error("Cannot create output: " + output_root);
    try {
      auto *tables = output->mkdir("tables");
      auto *sources = output->mkdir("sources");
      auto *metadata = output->mkdir("meta");
      if (!tables || !sources || !metadata) throw std::runtime_error("Cannot create output directories.");
      WriteText(*metadata, "Contract",
                "pair_purity, axis_unit=GeV; normal bins only; P=sum(W_i*P_i)/sum(W_i). "
                "Weights are the supplied PbPb weights; no mT-width factor, no OO reweighting. "
                "Purity errors: sqrt(sum(W_i^2*sigma_i^2))/sum(W_i), fixed weights and independent input errors. "
                "No input covariance is available. Weight errors are stored separately and not propagated into purity. "
                "Only zero-weight trailing bins are trimmed; internal/leading zero weights are errors. "
                "No weak-decay correction is applied to stored tables. Consumer uses (1-0.12)^2 exactly once. "
                "Conditional fit errors do not propagate these external table errors.");
      WriteText(*metadata, "FitCoverage", "fit_q_max=0.15 GeV/c; k_table=|q|/2; required_k_max_GeV=0.1299038105676658");
      for (std::size_t c = 0; c < centralities.size(); ++c) {
        auto *source_directory = sources->mkdir(("cent" + centralities[c]).c_str());
        if (!source_directory) throw std::runtime_error("Cannot create source directory.");
        WriteText(*source_directory, "SourcePath", source_paths[c]);
        WriteText(*source_directory, "SourceMD5", checksums[c]);
        TIter next(inputs[c]->GetListOfKeys());
        while (auto *key = static_cast<TKey *>(next())) {
          std::unique_ptr<TObject> object(key->ReadObj());
          TDirectory::TContext source_context(source_directory);
          if (!object || object->Write(key->GetName()) <= 0) {
            throw std::runtime_error("Failed to preserve source object: " + std::string(key->GetName()));
          }
        }
        for (std::size_t m = 0; m < mt_groups.size(); ++m) {
          auto &table = merged[c * mt_groups.size() + m];
          TDirectory::TContext table_context(tables);
          if (table.purity->Write() <= 0 || table.weight->Write() <= 0) throw std::runtime_error("Failed to write merged table.");
          std::string members;
          for (const auto &member : mt_groups[m].members) members += (members.empty() ? "" : ",") + member;
          WriteText(*metadata, table.purity->GetName(), "cent=" + centralities[c] + "; mT inputs=" + members);
          std::cout << "TABLE tables/" << table.purity->GetName() << " bins=" << table.purity->GetNbinsX()
                    << " k_max_GeV=" << table.purity->GetXaxis()->GetXmax() << '\n';
        }
      }
      WriteText(*metadata, "PreparationStatus", "complete");
      if (output->Write() < 0 || output->TestBit(TFile::kWriteError)) throw std::runtime_error("ROOT output write failed.");
      output->Close();
      if (output->TestBit(TFile::kWriteError)) throw std::runtime_error("ROOT output close failed.");
    } catch (...) {
      output.reset();
      fs::remove(destination);
      throw;
    }
    std::cout << "PREPARED_TABLES=9 OUTPUT=" << fs::absolute(destination).string() << '\n';
  }
}  // namespace oo_pbpb_lambda

/** ROOT entry point: input directory and a new output ROOT path are required. */
void prepare_oo_pbpb_lambda_tables(const char *input_directory, const char *output_root) {
  try {
    oo_pbpb_lambda::Prepare(input_directory, output_root);
  } catch (const std::exception &error) {
    std::cerr << "[error] " << error.what() << '\n';
    gSystem->Exit(1);
  }
}
