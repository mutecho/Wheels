#include "Comparison.h"

#include <TAxis.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH2.h>
#include <THnBase.h>
#include <THnSparse.h>
#include <TKey.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace femtocompare {
namespace {

Axis ReadAxis(const TAxis& a, const std::string& name) {
  Axis result{name, a.GetTitle(), {}};
  for (int i = 1; i <= a.GetNbins() + 1; ++i) result.edges.push_back(a.GetBinLowEdge(i));
  return result;
}

std::string Interval(const Axis& a, int bin) {
  std::ostringstream s;
  s << std::setprecision(4) << a.edges[bin - 1] << "-" << a.edges[bin];
  return s.str();
}

std::unique_ptr<TFile> Open(const std::string& path) {
  auto result = std::unique_ptr<TFile>(TFile::Open(path.c_str(), "READ"));
  if (!result || result->IsZombie()) throw std::runtime_error("Cannot open " + path);
  if (result->TestBit(TFile::kRecovered)) throw std::runtime_error("Recovered input requires inspection: " + path);
  return result;
}

std::unique_ptr<THnBase> ReadSparse(TFile& f, const std::string& task, bool mixed) {
  const std::string key = task + (mixed ? "/MixedEvent_3Dqn/" : "/SameEvent_3Dqn/") + sparseName;
  auto* obj = f.Get(key.c_str());
  auto* h = dynamic_cast<THnBase*>(obj);
  if (!h) throw std::runtime_error("Missing THnBase: " + std::string(f.GetName()) + ":" + key);
  return std::unique_ptr<THnBase>(h);
}

Map2D ReadMap(const TH2& h) {
  Map2D result{ReadAxis(*h.GetXaxis(), "delta_eta"), ReadAxis(*h.GetYaxis(), "delta_phi"), {}, {}};
  result.count.resize(result.x.bins() * result.y.bins());
  result.variance.resize(result.count.size());
  for (int y = 1; y <= h.GetNbinsY(); ++y) {
    for (int x = 1; x <= h.GetNbinsX(); ++x) {
      const auto i = (y - 1) * h.GetNbinsX() + x - 1;
      result.count[i] = h.GetBinContent(x, y);
      result.variance[i] = std::pow(h.GetBinError(x, y), 2);
    }
  }
  return result;
}

Histogram ProjectMap(const Map2D& m, bool yAxis) {
  Histogram result(yAxis ? m.y : m.x);
  for (int y = 0; y < m.y.bins(); ++y) {
    for (int x = 0; x < m.x.bins(); ++x) {
      const int i = y * m.x.bins() + x;
      const int bin = yAxis ? y : x;
      result.count[bin] += m.count[i];
      result.variance[bin] += m.variance[i];
    }
  }
  return result;
}

void Require(bool ok, const char* text) {
  if (!ok) throw std::runtime_error(std::string("Self-test failed: ") + text);
}

}  // namespace

Histogram::Histogram(Axis a) : axis(std::move(a)) {
  if (!axis.edges.empty()) {
    count.resize(axis.bins());
    variance.resize(axis.bins());
  }
}

double Histogram::total() const { return std::accumulate(count.begin(), count.end(), 0.); }

void RequireSameAxis(const Axis& a, const Axis& b) {
  if (a.edges.size() != b.edges.size()) throw std::runtime_error("Different bin counts: " + a.name);
  for (std::size_t i = 0; i < a.edges.size(); ++i) {
    if (std::abs(a.edges[i] - b.edges[i]) > 1.e-11 * std::max({1., std::abs(a.edges[i]), std::abs(b.edges[i])})) {
      throw std::runtime_error("Different bin edges: " + a.name);
    }
  }
}

Plan MakePlan(const THnBase& source) {
  if (source.GetNdimensions() != dimensions) throw std::runtime_error("Expected seven sparse axes");
  Plan p;
  const std::array<std::string, dimensions> names{"dk_out", "dk_side", "dk_long", "mt", "centrality", "qn", "pair_phi"};
  for (int d = 0; d < dimensions; ++d) p.axes[d] = ReadAxis(*source.GetAxis(d), names[d]);
  p.selections.push_back({"inclusive", "All normal bins"});
  // The +-0.05 GeV/c cube is only legal when these are exact source edges.
  bool aligned = true;
  for (int d = 0; d < 3; ++d) for (double edge : {-0.05, 0.05}) {
    aligned &= std::any_of(p.axes[d].edges.begin(), p.axes[d].edges.end(), [edge](double v) { return std::abs(v - edge) < 1.e-10; });
  }
  if (aligned) p.selections.push_back({"low_q_cube", "|DK_{out,side,long}| < 0.05 GeV/c", 0, 0, true});
  for (int c = 1; c <= p.axes[4].bins(); ++c)
    p.selections.push_back({"cent_" + std::to_string(c), "Centrality " + Interval(p.axes[4], c) + "%", c});
  for (int m = 1; m <= p.axes[3].bins(); ++m)
    p.selections.push_back({"mt_" + std::to_string(m), "m_{T} " + Interval(p.axes[3], m) + " (native axis units)", 0, m});
  for (int c = 1; c <= p.axes[4].bins(); ++c) for (int m = 1; m <= p.axes[3].bins(); ++m)
    p.selections.push_back({"cent_" + std::to_string(c) + "_mt_" + std::to_string(m),
                            "Centrality " + Interval(p.axes[4], c) + "%, m_{T} " + Interval(p.axes[3], m), c, m});
  return p;
}

Sample Accumulate(const THnBase& source, const Plan& p) {
  if (source.GetNdimensions() != dimensions) throw std::runtime_error("Sparse dimension mismatch");
  for (int d = 0; d < dimensions; ++d) {
    auto axis = ReadAxis(*source.GetAxis(d), p.axes[d].name);
    RequireSameAxis(p.axes[d], axis);
    if (axis.title != p.axes[d].title) throw std::runtime_error("Axis title mismatch: " + axis.name);
  }
  Sample out;
  out.entries = source.GetEntries();
  out.filledBins = source.GetNbins();
  out.storedSumw2 = source.GetCalculateErrors();
  out.spectra.resize(p.selections.size());
  for (auto& selection : out.spectra) for (int d = 0; d < dimensions; ++d) selection[d] = Histogram(p.axes[d]);
  for (auto pair : p.mapAxes) {
    Map2D m{p.axes[pair[0]], p.axes[pair[1]], {}, {}};
    m.count.resize(m.x.bins() * m.y.bins());
    m.variance.resize(m.count.size());
    out.maps.push_back(std::move(m));
  }
  // Native mT/centrality bins determine membership without splitting any bin.
  std::vector<std::vector<int>> membership((p.axes[4].bins() + 1) * (p.axes[3].bins() + 1));
  for (int c = 1; c <= p.axes[4].bins(); ++c) for (int m = 1; m <= p.axes[3].bins(); ++m) {
    auto& list = membership[c * (p.axes[3].bins() + 1) + m];
    for (std::size_t s = 0; s < p.selections.size(); ++s) {
      auto sel = p.selections[s];
      if ((!sel.centralityBin || sel.centralityBin == c) && (!sel.mtBin || sel.mtBin == m)) list.push_back(s);
    }
  }
  std::array<int, dimensions> coord{};
  for (Long64_t i = 0; i < source.GetNbins(); ++i) {
    const double weight = source.GetBinContent(i, coord.data());
    if (!std::isfinite(weight) || weight < 0) throw std::runtime_error("Nonfinite/negative input count");
    out.allWeight += weight;
    bool normal = true;
    for (int d = 0; d < dimensions; ++d) {
      if (coord[d] == 0) { out.underflow[d] += weight; normal = false; }
      if (coord[d] > p.axes[d].bins()) { out.overflow[d] += weight; normal = false; }
    }
    if (!normal || weight == 0) continue;
    out.normalWeight += weight;
    const double variance = out.storedSumw2 ? source.GetBinError2(i) : weight;
    if (!std::isfinite(variance) || variance < 0) throw std::runtime_error("Invalid source variance");
    bool low = true;
    for (int d = 0; d < 3; ++d)
      low &= p.axes[d].edges[coord[d] - 1] >= -0.05 - 1.e-10 && p.axes[d].edges[coord[d]] <= 0.05 + 1.e-10;
    for (int s : membership[coord[4] * (p.axes[3].bins() + 1) + coord[3]]) {
      if (p.selections[s].lowMomentum && !low) continue;
      for (int d = 0; d < dimensions; ++d) {
        auto& h = out.spectra[s][d];
        h.count[coord[d] - 1] += weight;
        h.variance[coord[d] - 1] += variance;
      }
    }
    for (std::size_t k = 0; k < p.mapAxes.size(); ++k) {
      auto& h = out.maps[k];
      const auto pair = p.mapAxes[k];
      const int bin = (coord[pair[1]] - 1) * h.x.bins() + coord[pair[0]] - 1;
      h.count[bin] += weight;
      h.variance[bin] += variance;
    }
  }
  // These are independent conservation checks on the one-pass projections.
  for (int d = 0; d < dimensions; ++d)
    if (std::abs(out.spectra[0][d].total() - out.normalWeight) > 1.e-8 * std::max(1., out.normalWeight))
      throw std::runtime_error("Projection failed count conservation");
  if (!(out.normalWeight > 0)) throw std::runtime_error("No pairs inside all normal bins");
  return out;
}

Point NormalizedBin(const Histogram& h, int bin, bool density) {
  const double n = h.total();
  if (n <= 0) return {};
  const double p = h.count[bin] / n;
  const double v = std::accumulate(h.variance.begin(), h.variance.end(), 0.);
  // Includes Cov(n_i,N)=Var(n_i). Does not treat the measured total as fixed.
  const double variance = std::max(0., ((1 - 2 * p) * h.variance[bin] + p * p * v) / (n * n));
  const double width = density ? h.axis.edges[bin + 1] - h.axis.edges[bin] : 1.;
  return {p / width, variance / (width * width), true};
}

Point Divide(Point a, Point b) {
  if (!a.valid || !b.valid || b.value <= 0) return {};
  const double r = a.value / b.value;
  const double v = a.variance / (b.value * b.value) + a.value * a.value * b.variance / std::pow(b.value, 4);
  return {r, v, std::isfinite(r) && std::isfinite(v)};
}

double TotalVariation(const Histogram& a, const Histogram& b) {
  RequireSameAxis(a.axis, b.axis);
  if (a.total() <= 0 || b.total() <= 0) return -1;
  double result = 0;
  for (int i = 0; i < a.axis.bins(); ++i) result += std::abs(a.count[i] / a.total() - b.count[i] / b.total());
  return 0.5 * result;
}

void RunSelfTests() {
  Histogram h({"toy", "toy", {0., 1., 3.}});
  h.count = {25, 75}; h.variance = h.count;
  auto p = NormalizedBin(h, 0);
  Require(std::abs(p.value - 0.25) < 1.e-12 && std::abs(p.variance - 0.001875) < 1.e-12, "normalized covariance");
  Require(std::abs(NormalizedBin(h, 1, true).value - 0.375) < 1.e-12, "variable-width density");
  Require(!Divide({1, 1, true}, {0, 0, true}).valid, "zero denominator masked");
  bool rejected = false;
  try { RequireSameAxis(h.axis, {"toy", "toy", {0, 1.1, 3}}); } catch (...) { rejected = true; }
  Require(rejected, "complete axis edge check");
  const int bins[dimensions] = {4, 4, 4, 2, 2, 1, 2};
  const double lo[dimensions] = {-0.1, -0.1, -0.1, 0.2, 0, 0, 0};
  const double hi[dimensions] = {0.1, 0.1, 0.1, 1.0, 80, 10, 3.141592653589793};
  THnSparseD sparse("toy_sparse", "toy_sparse", dimensions, bins, lo, hi);
  sparse.Sumw2();
  double pos[dimensions] = {0.025, 0.025, 0.025, 0.3, 5, 1, 0.2};
  sparse.Fill(pos, 2); sparse.Fill(pos, 3);
  pos[0] = 0.075; sparse.Fill(pos, 4);
  pos[0] = 0.2; sparse.Fill(pos, 7);
  const auto plan = MakePlan(sparse);
  const auto sample = Accumulate(sparse, plan);
  Require(sample.allWeight == 16 && sample.normalWeight == 9 && sample.overflow[0] == 7, "all-axis flow exclusion");
  Require(sample.spectra[1][0].total() == 5, "aligned low-momentum selection");
  Require(std::accumulate(sample.spectra[0][0].variance.begin(), sample.spectra[0][0].variance.end(), 0.) == 29,
          "preserve weighted Sumw2");
  std::cout << "Self-tests passed: axes, normalization covariance, variable widths, ratio masks, sparse selections, Sumw2.\n";
}

void RunComparison(const Options& o) {
  std::filesystem::create_directories(o.output);
  std::ofstream inventory(o.output / "inventory.csv"), flow(o.output / "flows.csv"), axes(o.output / "axes.csv");
  inventory << "sample,entries,filled_sparse_bins,stored_sumw2,sum_all_bins,sum_all_normal_bins,normal_fraction\n" << std::setprecision(16);
  flow << "sample,axis,underflow_weight,overflow_weight\n" << std::setprecision(16);
  axes << "axis,title,edge_index,edge\n" << std::setprecision(16);
  Plan plan;
  std::array<Sample, 4> samples;
  std::vector<FourMaps> qaMaps;
  std::vector<FourSpectra> qa;
  const std::vector<std::string> qaNames{"detadphidetadphi0Before_0", "detadphidetadphi0After_0", "at_PV_before", "at_PV_after"};
  for (const auto& name : qaNames) {
    qaMaps.push_back({"qa_" + name, "Close-pair QA: " + name, {}});
    for (const auto& axis : {"x", "y"}) qa.push_back({"qa_" + name + "_" + axis, "Close-pair QA: " + name + " (" + axis + ")", {}});
  }
  std::ofstream inputs(o.output / "inputs.txt");
  inputs << "A=" << o.fileA << "\nB=" << o.fileB << "\ntaskA=" << o.taskA << "\ntaskB=" << o.taskB << '\n';
  for (int file = 0; file < 2; ++file) {
    auto input = Open(file == 0 ? o.fileA : o.fileB);
    const auto task = file == 0 ? o.taskA : o.taskB;
    inputs << "\nFile " << (file == 0 ? 'A' : 'B') << " pair tasks present:\n";
    TIter next(input->GetListOfKeys());
    while (auto* key = dynamic_cast<TKey*>(next())) {
      if (std::string(key->GetName()).find("femto-dream-pair-task") == 0) inputs << key->GetName() << '\n';
    }
    for (int event = 0; event < 2; ++event) {
      const int index = file * 2 + event;
      const std::string label = std::string(file == 0 ? "A" : "B") + (event ? "_ME" : "_SE");
      std::cout << "Reading " << label << ": " << input->GetName() << std::endl;
      {
        auto sparse = ReadSparse(*input, task, event);
        if (index == 0) {
          plan = MakePlan(*sparse);
          for (const auto& axis : plan.axes) for (std::size_t b = 0; b < axis.edges.size(); ++b)
            axes << axis.name << ",\"" << axis.title << "\"," << b << ',' << axis.edges[b] << '\n';
        }
        samples[index] = Accumulate(*sparse, plan);
      }
      const auto& sample = samples[index];
      inventory << label << ',' << sample.entries << ',' << sample.filledBins << ',' << sample.storedSumw2 << ','
                << sample.allWeight << ',' << sample.normalWeight << ',' << sample.normalWeight / sample.allWeight << '\n';
      for (int d = 0; d < dimensions; ++d)
        flow << label << ',' << plan.axes[d].name << ',' << sample.underflow[d] << ',' << sample.overflow[d] << '\n';
      for (std::size_t q = 0; q < qaNames.size(); ++q) {
        const auto path = task + "/kTrack_kTrack/" + qaNames[q] + (event ? "_ME" : "_SE");
        auto* h = dynamic_cast<TH2*>(input->Get(path.c_str()));
        if (!h) throw std::runtime_error("Missing QA histogram: " + path);
        qaMaps[q].histogram[index] = ReadMap(*h);
        for (int d = 0; d < 2; ++d) qa[q * 2 + d].histogram[index] = ProjectMap(qaMaps[q].histogram[index], d == 1);
      }
      std::cout << label << " all-bin sum=" << sample.allWeight << ", all-normal sum=" << sample.normalWeight << std::endl;
    }
  }
  // Validate auxiliary maps before writing any derived ratios.
  for (const auto& q : qaMaps) for (int i = 1; i < 4; ++i) {
    RequireSameAxis(q.histogram[0].x, q.histogram[i].x);
    RequireSameAxis(q.histogram[0].y, q.histogram[i].y);
  }
  WritePlots(o, plan, samples, qa, qaMaps);
}

}  // namespace femtocompare
