#pragma once

#include <array>
#include <memory>
#include <string>

#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH3D.h"
#include "femto3d/AnalysisConfig.h"

namespace femto3d {

  struct SliceHistograms {
    std::unique_ptr<TH3D> source_histogram;
    std::array<std::unique_ptr<TH1D>, 6> projection_histograms;
  };

  [[nodiscard]] std::unique_ptr<TH3D> CreateTH3D(const std::string &name,
                                                 const AxisSpec &x_axis,
                                                 const AxisSpec &y_axis,
                                                 const AxisSpec &z_axis);

  void FillToTH3D(TH3D &histogram, const PairSeparationLCMS &separation, bool warn_on_overflow = true);

  [[nodiscard]] std::unique_ptr<TH1D> CreateProjectionHistogram(const std::string &name,
                                                                const AxisSpec &axis,
                                                                const ProjectionDefinition &direction);

  void FillProjectionHistograms(const PairSeparationLCMS &separation,
                                const std::array<ProjectionDefinition, 6> &directions,
                                std::array<std::unique_ptr<TH1D>, 6> &projection_histograms,
                                bool warn_on_overflow = true);

  [[nodiscard]] SliceHistograms CreateSliceHistograms(const std::string &slice_name,
                                                      const AnalysisConfig &config,
                                                      const std::array<ProjectionDefinition, 6> &directions);

  void WriteSliceHistograms(TDirectory &directory, SliceHistograms &histograms);

}  // namespace femto3d
