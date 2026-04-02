#include "femto3d/Histogramming.h"

#include <iostream>
#include <stdexcept>

namespace femto3d {

namespace {

void WarnIfOutsideHistogram(const TH3D& histogram,
                            const PairSeparationLCMS& separation) {
  const bool outside_x = separation.rho_out < histogram.GetXaxis()->GetXmin() ||
                         separation.rho_out >= histogram.GetXaxis()->GetXmax();
  const bool outside_y = separation.rho_side < histogram.GetYaxis()->GetXmin() ||
                         separation.rho_side >= histogram.GetYaxis()->GetXmax();
  const bool outside_z = separation.rho_long < histogram.GetZaxis()->GetXmin() ||
                         separation.rho_long >= histogram.GetZaxis()->GetXmax();

  if (outside_x || outside_y || outside_z) {
    std::cerr << "Warning: value (" << separation.rho_out << ", "
              << separation.rho_side << ", " << separation.rho_long
              << ") exceeds TH3D range for histogram " << histogram.GetName()
              << ".\n";
  }
}

void WarnIfOutsideHistogram(const TH1D& histogram, const double value) {
  if (value < histogram.GetXaxis()->GetXmin() ||
      value >= histogram.GetXaxis()->GetXmax()) {
    std::cerr << "Warning: value " << value
              << " exceeds TH1D range for histogram " << histogram.GetName()
              << ".\n";
  }
}

}  // namespace

std::unique_ptr<TH3D> CreateTH3D(const std::string& name,
                                 const AxisSpec& x_axis,
                                 const AxisSpec& y_axis,
                                 const AxisSpec& z_axis) {
  auto histogram = std::make_unique<TH3D>(name.c_str(),
                                          name.c_str(),
                                          x_axis.NumBins(),
                                          x_axis.min,
                                          x_axis.max,
                                          y_axis.NumBins(),
                                          y_axis.min,
                                          y_axis.max,
                                          z_axis.NumBins(),
                                          z_axis.min,
                                          z_axis.max);
  histogram->SetDirectory(nullptr);
  histogram->Sumw2();
  histogram->GetXaxis()->SetTitle(x_axis.title.c_str());
  histogram->GetYaxis()->SetTitle(y_axis.title.c_str());
  histogram->GetZaxis()->SetTitle(z_axis.title.c_str());
  return histogram;
}

void FillToTH3D(TH3D& histogram,
                const PairSeparationLCMS& separation,
                const bool warn_on_overflow) {
  if (warn_on_overflow) {
    WarnIfOutsideHistogram(histogram, separation);
  }

  histogram.Fill(separation.rho_out, separation.rho_side, separation.rho_long);
}

std::unique_ptr<TH1D> CreateProjectionHistogram(
    const std::string& name,
    const AxisSpec& axis,
    const ProjectionDefinition& direction) {
  auto histogram = std::make_unique<TH1D>(name.c_str(),
                                          direction.title.c_str(),
                                          axis.NumBins(),
                                          axis.min,
                                          axis.max);
  histogram->SetDirectory(nullptr);
  histogram->Sumw2();
  histogram->GetXaxis()->SetTitle(direction.title.c_str());
  histogram->GetYaxis()->SetTitle("Counts");
  return histogram;
}

void FillProjectionHistograms(
    const PairSeparationLCMS& separation,
    const std::array<ProjectionDefinition, 6>& directions,
    std::array<std::unique_ptr<TH1D>, 6>& projection_histograms,
    const bool warn_on_overflow) {
  for (std::size_t index = 0; index < directions.size(); ++index) {
    TH1D* histogram = projection_histograms[index].get();
    if (histogram == nullptr) {
      throw std::runtime_error("Projection histogram container is not initialized.");
    }

    const double projected_value = ProjectOntoDirection(separation, directions[index]);
    if (warn_on_overflow) {
      WarnIfOutsideHistogram(*histogram, projected_value);
    }

    histogram->Fill(projected_value);
  }
}

SliceHistograms CreateSliceHistograms(
    const std::string& slice_name,
    const AnalysisConfig& config,
    const std::array<ProjectionDefinition, 6>& directions) {
  SliceHistograms histograms;
  histograms.source_histogram =
      CreateTH3D("source_3d",
                 config.histograms.rho_out_axis,
                 config.histograms.rho_side_axis,
                 config.histograms.rho_long_axis);
  histograms.source_histogram->SetTitle(slice_name.c_str());

  for (std::size_t index = 0; index < directions.size(); ++index) {
    const std::string histogram_name =
        "proj_" + directions[index].short_name;
    histograms.projection_histograms[index] =
        CreateProjectionHistogram(histogram_name,
                                  config.histograms.projection_axis,
                                  directions[index]);
  }

  return histograms;
}

void WriteSliceHistograms(TDirectory& directory, SliceHistograms& histograms) {
  directory.cd();

  if (histograms.source_histogram != nullptr) {
    histograms.source_histogram->Write();
  }

  for (const auto& histogram : histograms.projection_histograms) {
    if (histogram != nullptr) {
      histogram->Write();
    }
  }
}

}  // namespace femto3d
