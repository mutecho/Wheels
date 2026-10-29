#pragma once

#include <array>
#include <filesystem>
#include <string>
#include <vector>

class THnBase;

namespace femtocompare {

constexpr int dimensions = 7;
inline constexpr const char* sparseName = "relPair3dRmTMultPercentileQnPairphi";

struct Axis {
  std::string name, title;
  std::vector<double> edges;
  int bins() const { return static_cast<int>(edges.size()) - 1; }
};

/// Independent source-bin variances are summed before normalization. Input
/// Sumw2 is retained; absent Sumw2 uses ROOT's unweighted-count convention.
struct Histogram {
  Axis axis;
  std::vector<double> count, variance;
  Histogram() = default;
  explicit Histogram(Axis a);
  double total() const;
};

struct Point {
  double value = 0, variance = 0;
  bool valid = false;
};

struct Selection {
  std::string name, title;
  int centralityBin = 0, mtBin = 0;
  bool lowMomentum = false;
};

struct Map2D {
  Axis x, y;
  std::vector<double> count, variance;
};

struct Plan {
  std::array<Axis, dimensions> axes;
  std::vector<Selection> selections;
  // Ordered (x,y) pairs; all maps use the identical all-normal-bin sample.
  std::vector<std::array<int, 2>> mapAxes{{0, 1}, {0, 2}, {1, 2}, {3, 4}, {6, 4}, {6, 3}};
};

/// A single SE or ME input. No source object survives file close.
struct Sample {
  std::vector<std::array<Histogram, dimensions>> spectra;
  std::vector<Map2D> maps;
  double entries = 0, allWeight = 0, normalWeight = 0;
  long long filledBins = 0;
  bool storedSumw2 = false;
  std::array<double, dimensions> underflow{}, overflow{};
};

struct FourSpectra {
  std::string name, title;
  // A same, A mixed, B same, B mixed throughout the program.
  std::array<Histogram, 4> histogram;
};

struct FourMaps {
  std::string name, title;
  std::array<Map2D, 4> histogram;
};

struct Options {
  std::string fileA = "/Users/allenzhou/ALICE/alidata/hyperloop_res/femtoep/OO/AnalysisResults_3d_25ae_pass2_loosecuttest.root";
  std::string fileB = "/Users/allenzhou/ALICE/alidata/hyperloop_setup/femtoep_oo/run/res_20261020/AnalysisResults.root";
  std::string taskA = "femto-dream-pair-task-track-track", taskB = taskA;
  std::string labelA = "A: 25ae loosecuttest", labelB = "B: local 20261020";
  std::filesystem::path output = "comparison_25ae_vs_local_20261020";
};

/// Reject different complete edge arrays instead of aligning by bin centers.
void RequireSameAxis(const Axis& a, const Axis& b);
Plan MakePlan(const THnBase& source);
/// One sparse pass fills all inclusive and native centrality/mT slices. Every
/// accepted pair must be in normal bins on ALL seven axes. No rebinning occurs.
Sample Accumulate(const THnBase& source, const Plan& plan);
Point NormalizedBin(const Histogram& h, int bin, bool density = false);
/// Ratios of positive denominators only; validity is separate from the value.
Point Divide(Point numerator, Point denominator);
double TotalVariation(const Histogram& a, const Histogram& b);
void RunSelfTests();
void RunComparison(const Options& options);

/// Writes individual PNGs, a multipage PDF, raw/derived ROOT objects and CSVs.
void WritePlots(const Options& options, const Plan& plan, const std::array<Sample, 4>& samples,
                const std::vector<FourSpectra>& qa, const std::vector<FourMaps>& qaMaps);

}  // namespace femtocompare
