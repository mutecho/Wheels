#include "femto3d/Workflow.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TTree.h"
#include "femto3d/EventPlane.h"
#include "femto3d/Histogramming.h"
#include "femto3d/InputReader.h"
#include "femto3d/ProjectionFit.h"
#include "femto3d/SourceExtraction.h"

namespace femto3d {

  namespace {

    struct RadiusSummaryPoint {
      double phi_center = 0.0;
      double phi_error = 0.0;
      double value = 0.0;
      double error = 0.0;
      bool valid = false;
    };

    using RadiusSummaryMap =
        std::map<std::pair<std::size_t, std::size_t>, std::array<std::vector<RadiusSummaryPoint>, 6>>;

    struct SummaryGraphStats {
      int point_count = 0;
      double mean_y = 0.0;
      double rms_y = 0.0;
    };

    int FindBinIndex(const std::vector<RangeBin> &bins, const double value) {
      for (std::size_t index = 0; index < bins.size(); ++index) {
        if (bins[index].Contains(value)) {
          return static_cast<int>(index);
        }
      }

      if (!bins.empty() && std::abs(value - bins.back().max) < 1.0e-12) {
        return static_cast<int>(bins.size() - 1U);
      }

      return -1;
    }

    TDirectory *GetOrCreateDirectory(TDirectory &parent, const std::string &name) {
      if (TDirectory *existing = parent.GetDirectory(name.c_str())) {
        return existing;
      }

      TDirectory *created = parent.mkdir(name.c_str());
      if (created == nullptr) {
        throw std::runtime_error("Failed to create directory: " + name);
      }

      return created;
    }

    std::string BuildSliceName(const RangeBin &centrality_bin, const RangeBin &mt_bin, const RangeBin &phi_bin) {
      return centrality_bin.label + "__" + mt_bin.label + "__" + phi_bin.label;
    }

    std::string BuildRelativePhiAxisTitle(const int harmonic_order) {
      return "#phi_{pair} - #Psi_{" + std::to_string(harmonic_order) + "}";
    }

    std::string FormatSummaryStat(const double value) {
      std::ostringstream stream;
      stream << std::fixed << std::setprecision(3) << value;
      return stream.str();
    }

    SummaryGraphStats ComputeSummaryGraphStats(const std::vector<RadiusSummaryPoint> &points) {
      SummaryGraphStats stats;
      double sum_y = 0.0;
      double sum_y2 = 0.0;

      for (const RadiusSummaryPoint &point : points) {
        if (!point.valid || !std::isfinite(point.value)) {
          continue;
        }

        ++stats.point_count;
        sum_y += point.value;
        sum_y2 += point.value * point.value;
      }

      if (stats.point_count <= 0) {
        return stats;
      }

      const double normalization = static_cast<double>(stats.point_count);
      stats.mean_y = sum_y / normalization;
      const double mean_y2 = sum_y2 / normalization;
      stats.rms_y = std::sqrt(std::max(0.0, mean_y2 - stats.mean_y * stats.mean_y));
      return stats;
    }

    void ApplySummaryGraphStyle(TGraphErrors &graph) {
      constexpr int kSummaryGraphColor = 602;
      constexpr int kSummaryMarkerStyle = 20;

      graph.SetMarkerStyle(kSummaryMarkerStyle);
      graph.SetMarkerSize(1.2);
      graph.SetMarkerColor(kSummaryGraphColor);
      graph.SetLineColor(kSummaryGraphColor);
      graph.SetLineWidth(2);
      graph.SetFillStyle(0);
      graph.SetDrawOption("ALPE1");
    }

    TPaveText *CreateSummaryInfoBox(const std::string &graph_name, const SummaryGraphStats &stats) {
      auto *info_box = new TPaveText(0.63, 0.72, 0.89, 0.89, "NDC");
      info_box->SetName((graph_name + "_summary_box").c_str());
      info_box->SetBorderSize(1);
      info_box->SetFillColor(0);
      info_box->SetFillStyle(1001);
      info_box->SetShadowColor(0);
      info_box->SetTextAlign(12);
      info_box->SetTextFont(42);
      info_box->AddText(("Points: " + std::to_string(stats.point_count)).c_str());
      info_box->AddText(("Mean(Y): " + FormatSummaryStat(stats.mean_y)).c_str());
      info_box->AddText(("RMS(Y): " + FormatSummaryStat(stats.rms_y)).c_str());
      return info_box;
    }

    bool PassesFemtoSelection(const ParticleState &particle, const ParticleSelectionConfig &selection) {
      if (particle.pdg != selection.target_pdg) {
        return false;
      }

      const double eta = ComputeParticleEta(particle);
      if (!std::isfinite(eta) || eta < selection.femto_eta_min || eta >= selection.femto_eta_max) {
        return false;
      }

      const double pt = ComputeParticlePt(particle);
      return pt >= selection.femto_pt_min && pt < selection.femto_pt_max;
    }

    EventPlaneResult ResolveEventPlane(const AnalysisConfig &config,
                                       const bool use_internal_event_plane,
                                       const bool use_input_event_plane,
                                       const std::vector<ParticleState> &event_plane_candidates,
                                       const EventData &event_data) {
      EventPlaneResult result;
      result.failure_reason = EventPlaneFailureReason::kDisabled;

      if (use_internal_event_plane) {
        result = ReconstructEventPlane(event_plane_candidates, config.event_plane);
        if (result.success) {
          return result;
        }
      }

      if (use_input_event_plane) {
        if (event_data.has_input_event_plane) {
          return MakeInputEventPlaneResult(event_data.event_plane_psi, config.event_plane.harmonic_order);
        }

        if (!result.success && result.failure_reason == EventPlaneFailureReason::kDisabled) {
          result.failure_reason = EventPlaneFailureReason::kMissingInputBranch;
        }
      }

      return result;
    }

    void RecordEventPlaneFailure(const EventPlaneResult &event_plane_result, AnalysisStatistics &statistics) {
      switch (event_plane_result.failure_reason) {
        case EventPlaneFailureReason::kNoCandidates:
        case EventPlaneFailureReason::kInsufficientCandidates:
          ++statistics.events_rejected_no_event_plane_candidates;
          break;
        case EventPlaneFailureReason::kQTooSmall:
          ++statistics.events_rejected_small_q_vector;
          break;
        case EventPlaneFailureReason::kMissingInputBranch:
          ++statistics.events_rejected_missing_input_event_plane;
          break;
        case EventPlaneFailureReason::kNone:
        case EventPlaneFailureReason::kDisabled:
          break;
      }
    }

    void FillResultArrays(const std::array<ProjectionFitResult, 6> &results,
                          double values[6],
                          double errors[6],
                          int success_flags[6]) {
      for (std::size_t index = 0; index < results.size(); ++index) {
        values[index] = results[index].r2;
        errors[index] = results[index].r2_error;
        success_flags[index] = results[index].success ? 1 : 0;
      }
    }

    void FillBoolFlags(const std::array<bool, 6> &flags, int output_flags[6]) {
      for (std::size_t index = 0; index < flags.size(); ++index) {
        output_flags[index] = flags[index] ? 1 : 0;
      }
    }

    void WriteAnalysisStatistics(TDirectory &directory, const AnalysisStatistics &statistics) {
      long long events_read = static_cast<long long>(statistics.events_read);
      long long selected_particles = static_cast<long long>(statistics.selected_particles);
      long long events_with_valid_event_plane = static_cast<long long>(statistics.events_with_valid_event_plane);
      long long events_rejected_no_event_plane_candidates =
          static_cast<long long>(statistics.events_rejected_no_event_plane_candidates);
      long long events_rejected_small_q_vector = static_cast<long long>(statistics.events_rejected_small_q_vector);
      long long events_rejected_missing_input_event_plane =
          static_cast<long long>(statistics.events_rejected_missing_input_event_plane);
      long long events_rejected_insufficient_femto_particles =
          static_cast<long long>(statistics.events_rejected_insufficient_femto_particles);
      long long candidate_pairs = static_cast<long long>(statistics.candidate_pairs);
      long long accepted_pairs = static_cast<long long>(statistics.accepted_pairs);
      long long rejected_close_pairs = static_cast<long long>(statistics.rejected_close_pairs);
      long long rejected_invalid_pair_kinematics = static_cast<long long>(statistics.rejected_invalid_pair_kinematics);
      long long r2_summary_points_skipped_invalid_hbt_error =
          static_cast<long long>(statistics.r2_summary_points_skipped_invalid_hbt_error);

      TTree statistics_tree("analysis_statistics", "analysis_statistics");
      statistics_tree.Branch("events_read", &events_read, "events_read/L");
      statistics_tree.Branch("selected_particles", &selected_particles, "selected_particles/L");
      statistics_tree.Branch(
          "events_with_valid_event_plane", &events_with_valid_event_plane, "events_with_valid_event_plane/L");
      statistics_tree.Branch("events_rejected_no_event_plane_candidates",
                             &events_rejected_no_event_plane_candidates,
                             "events_rejected_no_event_plane_candidates/L");
      statistics_tree.Branch(
          "events_rejected_small_q_vector", &events_rejected_small_q_vector, "events_rejected_small_q_vector/L");
      statistics_tree.Branch("events_rejected_missing_input_event_plane",
                             &events_rejected_missing_input_event_plane,
                             "events_rejected_missing_input_event_plane/L");
      statistics_tree.Branch("events_rejected_insufficient_femto_particles",
                             &events_rejected_insufficient_femto_particles,
                             "events_rejected_insufficient_femto_particles/L");
      statistics_tree.Branch("candidate_pairs", &candidate_pairs, "candidate_pairs/L");
      statistics_tree.Branch("accepted_pairs", &accepted_pairs, "accepted_pairs/L");
      statistics_tree.Branch("rejected_close_pairs", &rejected_close_pairs, "rejected_close_pairs/L");
      statistics_tree.Branch(
          "rejected_invalid_pair_kinematics", &rejected_invalid_pair_kinematics, "rejected_invalid_pair_kinematics/L");
      statistics_tree.Branch("r2_summary_points_skipped_invalid_hbt_error",
                             &r2_summary_points_skipped_invalid_hbt_error,
                             "r2_summary_points_skipped_invalid_hbt_error/L");
      statistics_tree.Fill();
      directory.cd();
      statistics_tree.Write();
    }

    void WriteEventPlaneMetadata(TDirectory &directory, const AnalysisConfig &config) {
      int harmonic_order = config.event_plane.harmonic_order;
      int enabled = config.event_plane.enabled ? 1 : 0;
      int use_internal_reconstruction = config.event_plane.use_internal_reconstruction ? 1 : 0;
      int fallback_to_input_branch = config.event_plane.fallback_to_input_branch ? 1 : 0;
      int weight_mode = static_cast<int>(config.event_plane.weight_mode);
      double eta_min = config.event_plane.eta_min;
      double eta_max = config.event_plane.eta_max;
      long long min_candidates = static_cast<long long>(config.event_plane.min_candidates);
      double min_q_magnitude = config.event_plane.min_q_magnitude;

      TTree metadata_tree("event_plane_metadata", "event_plane_metadata");
      metadata_tree.Branch("harmonic_order", &harmonic_order, "harmonic_order/I");
      metadata_tree.Branch("enabled", &enabled, "enabled/I");
      metadata_tree.Branch(
          "use_internal_reconstruction", &use_internal_reconstruction, "use_internal_reconstruction/I");
      metadata_tree.Branch("fallback_to_input_branch", &fallback_to_input_branch, "fallback_to_input_branch/I");
      metadata_tree.Branch("eta_min", &eta_min, "eta_min/D");
      metadata_tree.Branch("eta_max", &eta_max, "eta_max/D");
      metadata_tree.Branch("weight_mode", &weight_mode, "weight_mode/I");
      metadata_tree.Branch("min_candidates", &min_candidates, "min_candidates/L");
      metadata_tree.Branch("min_q_magnitude", &min_q_magnitude, "min_q_magnitude/D");
      metadata_tree.Fill();
      directory.cd();
      metadata_tree.Write();
    }

    void WriteMtSlicingMetadata(TDirectory &directory, const AnalysisConfig &config) {
      int uses_standard_femto_mt = 1;
      int pair_system_mt_used_for_lcms = 1;
      int uses_explicit_reference_mass = config.selection.femto_mt_reference_mass >= 0.0 ? 1 : 0;
      int fit_summary_fields_are_femto_mt = 1;
      int directory_labels_are_femto_mt = 1;
      int target_pdg = config.selection.target_pdg;
      double femto_mt_reference_mass = config.selection.femto_mt_reference_mass;
      double femto_mt_mass_tolerance = config.selection.femto_mt_mass_tolerance;

      TTree metadata_tree("mt_slicing_metadata", "mt_slicing_metadata");
      metadata_tree.Branch("uses_standard_femto_mt", &uses_standard_femto_mt, "uses_standard_femto_mt/I");
      metadata_tree.Branch(
          "pair_system_mt_used_for_lcms", &pair_system_mt_used_for_lcms, "pair_system_mt_used_for_lcms/I");
      metadata_tree.Branch(
          "uses_explicit_reference_mass", &uses_explicit_reference_mass, "uses_explicit_reference_mass/I");
      metadata_tree.Branch(
          "fit_summary_fields_are_femto_mt", &fit_summary_fields_are_femto_mt, "fit_summary_fields_are_femto_mt/I");
      metadata_tree.Branch(
          "directory_labels_are_femto_mt", &directory_labels_are_femto_mt, "directory_labels_are_femto_mt/I");
      metadata_tree.Branch("target_pdg", &target_pdg, "target_pdg/I");
      metadata_tree.Branch("femto_mt_reference_mass", &femto_mt_reference_mass, "femto_mt_reference_mass/D");
      metadata_tree.Branch("femto_mt_mass_tolerance", &femto_mt_mass_tolerance, "femto_mt_mass_tolerance/D");
      metadata_tree.Fill();
      directory.cd();
      metadata_tree.Write();
    }

    void WriteProjectionFitMetadata(TDirectory &directory, const AnalysisConfig &config) {
      int use_adaptive_integration = config.projection_fit.use_adaptive_integration ? 1 : 0;
      int accept_forced_posdef_covariance_as_valid =
          config.projection_fit.accept_forced_posdef_covariance_as_valid ? 1 : 0;
      int fail_alpha_result_when_error_invalid = config.projection_fit.fail_alpha_result_when_error_invalid ? 1 : 0;
      int fail_hbt_results_when_error_invalid = config.projection_fit.fail_hbt_results_when_error_invalid ? 1 : 0;
      int fail_directional_results_when_error_invalid =
          config.projection_fit.fail_directional_results_when_error_invalid ? 1 : 0;
      int accept_hbt_central_value_only_for_summary =
          config.projection_fit.accept_hbt_central_value_only_for_summary ? 1 : 0;

      TTree metadata_tree("projection_fit_metadata", "projection_fit_metadata");
      metadata_tree.Branch("use_adaptive_integration", &use_adaptive_integration, "use_adaptive_integration/I");
      metadata_tree.Branch("accept_forced_posdef_covariance_as_valid",
                           &accept_forced_posdef_covariance_as_valid,
                           "accept_forced_posdef_covariance_as_valid/I");
      metadata_tree.Branch("fail_alpha_result_when_error_invalid",
                           &fail_alpha_result_when_error_invalid,
                           "fail_alpha_result_when_error_invalid/I");
      metadata_tree.Branch("fail_hbt_results_when_error_invalid",
                           &fail_hbt_results_when_error_invalid,
                           "fail_hbt_results_when_error_invalid/I");
      metadata_tree.Branch("fail_directional_results_when_error_invalid",
                           &fail_directional_results_when_error_invalid,
                           "fail_directional_results_when_error_invalid/I");
      metadata_tree.Branch("accept_hbt_central_value_only_for_summary",
                           &accept_hbt_central_value_only_for_summary,
                           "accept_hbt_central_value_only_for_summary/I");
      metadata_tree.Fill();
      directory.cd();
      metadata_tree.Write();
    }

    void WriteR2SummaryGraphs(TDirectory &directory, const AnalysisConfig &config, RadiusSummaryMap &summary_points) {
      const std::array<std::string, 6> graph_names = {
          "Rout2_vs_phi", "Rside2_vs_phi", "Rlong2_vs_phi", "Ros2_vs_phi", "Rol2_vs_phi", "Rsl2_vs_phi"};
      const std::array<std::string, 6> graph_titles = {
          "R_{out}^{2}", "R_{side}^{2}", "R_{long}^{2}", "R_{os}^{2}", "R_{ol}^{2}", "R_{sl}^{2}"};
      const std::string relative_phi_title = BuildRelativePhiAxisTitle(config.event_plane.harmonic_order);

      for (auto &entry : summary_points) {
        const std::size_t cent_index = entry.first.first;
        const std::size_t mt_index = entry.first.second;
        TDirectory *cent_directory = GetOrCreateDirectory(directory, config.centrality_bins[cent_index].label);
        TDirectory *mt_directory = GetOrCreateDirectory(*cent_directory, config.mt_bins[mt_index].label);

        for (std::size_t radius_index = 0; radius_index < entry.second.size(); ++radius_index) {
          auto &points = entry.second[radius_index];
          std::sort(points.begin(), points.end(), [](const RadiusSummaryPoint &left, const RadiusSummaryPoint &right) {
            return left.phi_center < right.phi_center;
          });

          const auto valid_count =
              static_cast<int>(std::count_if(points.begin(), points.end(), [](const RadiusSummaryPoint &point) {
                return point.valid;
              }));
          if (valid_count == 0) {
            continue;
          }

          auto graph = std::make_unique<TGraphErrors>(valid_count);
          graph->SetName(graph_names[radius_index].c_str());
          graph->SetTitle((graph_titles[radius_index] + " vs " + relative_phi_title).c_str());
          graph->GetXaxis()->SetTitle(relative_phi_title.c_str());
          graph->GetYaxis()->SetTitle(graph_titles[radius_index].c_str());
          ApplySummaryGraphStyle(*graph);

          const SummaryGraphStats stats = ComputeSummaryGraphStats(points);

          int point_index = 0;
          for (const RadiusSummaryPoint &point : points) {
            if (!point.valid) {
              continue;
            }

            graph->SetPoint(point_index, point.phi_center, point.value);
            graph->SetPointError(point_index, point.phi_error, point.error);
            ++point_index;
          }

          mt_directory->cd();
          auto canvas =
              std::make_unique<TCanvas>((graph_names[radius_index] + "_canvas").c_str(), graph->GetTitle(), 800, 600);
          canvas->SetTicks(1, 1);
          canvas->cd();
          graph->Draw("ALPE1");
          CreateSummaryInfoBox(graph_names[radius_index], stats)->Draw();
          canvas->Modified();
          canvas->Update();

          graph->Write();
          canvas->Write();
        }
      }
    }

  }  // namespace

  AnalysisStatistics RunAnalysis(const ApplicationConfig &application_config) {
    const AnalysisConfig &config = application_config.analysis;
    if (config.centrality_bins.empty() || config.mt_bins.empty() || config.phi_bins.empty()) {
      throw std::invalid_argument("Analysis bins must not be empty.");
    }
    if (config.selection.femto_mt_mass_tolerance < 0.0) {
      throw std::invalid_argument("Femto mT mass tolerance must be non-negative.");
    }

    const bool use_internal_event_plane = config.event_plane.enabled && config.event_plane.use_internal_reconstruction;
    const bool use_input_event_plane = config.event_plane.fallback_to_input_branch;
    if (!use_internal_event_plane && !use_input_event_plane) {
      throw std::invalid_argument(
          "Analysis requires an event-plane source: enable internal "
          "reconstruction or input-branch fallback.");
    }

    const std::array<ProjectionDefinition, 6> directions = MakeProjectionDefinitions();

    const std::vector<EventData> events = LoadEventData(application_config);

    AnalysisStatistics statistics;
    std::map<SliceKey, SliceHistograms> slice_histograms;

    for (const EventData &event_data : events) {
      ++statistics.events_read;

      const int cent_index = FindBinIndex(config.centrality_bins, event_data.centrality);
      if (cent_index < 0) {
        continue;
      }

      const std::size_t particle_count = event_data.particles.size();

      std::vector<ParticleState> event_plane_candidates;
      std::vector<ParticleState> selected_particles;
      selected_particles.reserve(particle_count);
      if (use_internal_event_plane) {
        event_plane_candidates.reserve(particle_count);
      }

      for (std::size_t particle_index = 0; particle_index < particle_count; ++particle_index) {
        const ParticleState &particle = event_data.particles[particle_index];

        if (use_internal_event_plane && IsEventPlaneCandidate(particle, config.event_plane)) {
          event_plane_candidates.push_back(particle);
        }

        if (PassesFemtoSelection(particle, config.selection)) {
          selected_particles.push_back(particle);
        }
      }

      statistics.selected_particles += selected_particles.size();

      const EventPlaneResult event_plane_result = ResolveEventPlane(
          config, use_internal_event_plane, use_input_event_plane, event_plane_candidates, event_data);
      if (!event_plane_result.success) {
        RecordEventPlaneFailure(event_plane_result, statistics);
        continue;
      }

      ++statistics.events_with_valid_event_plane;

      if (selected_particles.size() < 2U) {
        ++statistics.events_rejected_insufficient_femto_particles;
        continue;
      }

      for (std::size_t first = 0; first < selected_particles.size(); ++first) {
        for (std::size_t second = first + 1U; second < selected_particles.size(); ++second) {
          ++statistics.candidate_pairs;

          const ParticleState &particle1 = selected_particles[first];
          const ParticleState &particle2 = selected_particles[second];

          try {
            const PairKinematics pair_kinematics = ComputePairKinematics(particle1,
                                                                         particle2,
                                                                         config.selection.femto_mt_reference_mass,
                                                                         config.selection.femto_mt_mass_tolerance);
            const int mt_index = FindBinIndex(config.mt_bins, pair_kinematics.femto_mt);
            if (mt_index < 0) {
              continue;
            }

            const double phi_minus_psi = ComputePairPhiMinusEventPlane(
                pair_kinematics, event_plane_result.psi, config.event_plane.harmonic_order);
            const int phi_index = FindBinIndex(config.phi_bins, phi_minus_psi);
            if (phi_index < 0) {
              continue;
            }

            const PairSeparationLCMS separation =
                ComputeSinglePairDistanceLCMS(particle1, particle2, config.close_pair_distance_tolerance);

            const SliceKey key{static_cast<std::size_t>(cent_index),
                               static_cast<std::size_t>(mt_index),
                               static_cast<std::size_t>(phi_index)};

            auto iterator = slice_histograms.find(key);
            if (iterator == slice_histograms.end()) {
              const std::string slice_name = BuildSliceName(
                  config.centrality_bins[key.cent_index], config.mt_bins[key.mt_index], config.phi_bins[key.phi_index]);
              iterator = slice_histograms.emplace(key, CreateSliceHistograms(slice_name, config, directions)).first;
            }

            FillToTH3D(*iterator->second.source_histogram, separation, config.histograms.warn_on_overflow);
            FillProjectionHistograms(
                separation, directions, iterator->second.projection_histograms, config.histograms.warn_on_overflow);
            ++statistics.accepted_pairs;
          } catch (const std::invalid_argument &error) {
            const std::string message = error.what();
            if (message.find("Close pair rejection") != std::string::npos) {
              ++statistics.rejected_close_pairs;
            } else {
              ++statistics.rejected_invalid_pair_kinematics;
            }
          }
        }
      }
    }

    std::unique_ptr<TFile> output_file(TFile::Open(application_config.output_root_path.c_str(), "RECREATE"));
    if (output_file == nullptr || output_file->IsZombie()) {
      throw std::runtime_error("Failed to create output ROOT file: " + application_config.output_root_path);
    }

    TDirectory *top_directory = GetOrCreateDirectory(*output_file, config.top_directory_name);
    TDirectory *r2_summary_directory = GetOrCreateDirectory(*output_file, config.r2_summary_directory_name);

    int cent_bin_index = 0;
    int femto_mt_bin_index = 0;
    int phi_bin_index = 0;
    double cent_min = 0.0;
    double cent_max = 0.0;
    double femto_mt_min = 0.0;
    double femto_mt_max = 0.0;
    double phi_min = 0.0;
    double phi_max = 0.0;
    double directional_r2[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double directional_r2_error[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int directional_success[6] = {0, 0, 0, 0, 0, 0};
    double hbt_r2[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double hbt_r2_error[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int hbt_success[6] = {0, 0, 0, 0, 0, 0};
    double alpha = 0.0;
    double alpha_error = 0.0;
    int alpha_success = 0;
    int alpha_error_valid = 0;
    int covariance_status = 0;
    int covariance_valid = 0;
    int hbt_error_valid[6] = {0, 0, 0, 0, 0, 0};
    int directional_error_valid[6] = {0, 0, 0, 0, 0, 0};

    TTree fit_summary_tree("fit_summary", "fit_summary");
    fit_summary_tree.Branch("cent_bin_index", &cent_bin_index, "cent_bin_index/I");
    fit_summary_tree.Branch("femto_mt_bin_index", &femto_mt_bin_index, "femto_mt_bin_index/I");
    fit_summary_tree.Branch("phi_bin_index", &phi_bin_index, "phi_bin_index/I");
    fit_summary_tree.Branch("cent_min", &cent_min, "cent_min/D");
    fit_summary_tree.Branch("cent_max", &cent_max, "cent_max/D");
    fit_summary_tree.Branch("femto_mt_min", &femto_mt_min, "femto_mt_min/D");
    fit_summary_tree.Branch("femto_mt_max", &femto_mt_max, "femto_mt_max/D");
    fit_summary_tree.Branch("phi_min", &phi_min, "phi_min/D");
    fit_summary_tree.Branch("phi_max", &phi_max, "phi_max/D");
    fit_summary_tree.Branch("directional_r2", directional_r2, "directional_r2[6]/D");
    fit_summary_tree.Branch("directional_r2_error", directional_r2_error, "directional_r2_error[6]/D");
    fit_summary_tree.Branch("directional_success", directional_success, "directional_success[6]/I");
    fit_summary_tree.Branch("hbt_r2", hbt_r2, "hbt_r2[6]/D");
    fit_summary_tree.Branch("hbt_r2_error", hbt_r2_error, "hbt_r2_error[6]/D");
    fit_summary_tree.Branch("hbt_success", hbt_success, "hbt_success[6]/I");
    fit_summary_tree.Branch("hbt_error_valid", hbt_error_valid, "hbt_error_valid[6]/I");
    fit_summary_tree.Branch("alpha", &alpha, "alpha/D");
    fit_summary_tree.Branch("alpha_error", &alpha_error, "alpha_error/D");
    fit_summary_tree.Branch("alpha_success", &alpha_success, "alpha_success/I");
    fit_summary_tree.Branch("alpha_error_valid", &alpha_error_valid, "alpha_error_valid/I");
    fit_summary_tree.Branch("covariance_status", &covariance_status, "covariance_status/I");
    fit_summary_tree.Branch("covariance_valid", &covariance_valid, "covariance_valid/I");
    fit_summary_tree.Branch("directional_error_valid", directional_error_valid, "directional_error_valid[6]/I");

    RadiusSummaryMap summary_points;

    for (auto &entry : slice_histograms) {
      const SliceKey &key = entry.first;
      SliceHistograms &histograms = entry.second;
      const RangeBin &centrality_bin = config.centrality_bins[key.cent_index];
      const RangeBin &mt_bin = config.mt_bins[key.mt_index];
      const RangeBin &phi_bin = config.phi_bins[key.phi_index];
      const std::string slice_name = BuildSliceName(centrality_bin, mt_bin, phi_bin);

      TDirectory *cent_directory = GetOrCreateDirectory(*top_directory, centrality_bin.label);
      TDirectory *mt_directory = GetOrCreateDirectory(*cent_directory, mt_bin.label);
      TDirectory *phi_directory = GetOrCreateDirectory(*mt_directory, phi_bin.label);

      SliceFitProducts fit_products =
          FitSliceHistograms(slice_name, histograms.projection_histograms, directions, config.projection_fit);

      phi_directory->cd();
      WriteSliceHistograms(*phi_directory, histograms);
      for (const auto &fit_function : fit_products.fit_functions) {
        if (fit_function != nullptr) {
          fit_function->Write();
        }
      }
      if (fit_products.fit_canvas != nullptr) {
        fit_products.fit_canvas->Write();
      }

      cent_bin_index = static_cast<int>(key.cent_index);
      femto_mt_bin_index = static_cast<int>(key.mt_index);
      phi_bin_index = static_cast<int>(key.phi_index);
      cent_min = centrality_bin.min;
      cent_max = centrality_bin.max;
      femto_mt_min = mt_bin.min;
      femto_mt_max = mt_bin.max;
      phi_min = phi_bin.min;
      phi_max = phi_bin.max;
      FillResultArrays(fit_products.directional_results, directional_r2, directional_r2_error, directional_success);
      ApplyDirectionalErrorValidityMask(fit_products.directional_error_valid, directional_r2_error);
      FillResultArrays(fit_products.hbt_radii_results, hbt_r2, hbt_r2_error, hbt_success);
      ApplyHbtErrorValidityMask(fit_products.hbt_error_valid, hbt_r2_error);
      alpha = fit_products.alpha;
      alpha_error = fit_products.alpha_error;
      ApplyAlphaErrorValidityMask(fit_products.alpha_error_valid, &alpha_error);
      alpha_success = fit_products.alpha_success ? 1 : 0;
      alpha_error_valid = fit_products.alpha_error_valid ? 1 : 0;
      covariance_status = fit_products.covariance_status;
      covariance_valid = fit_products.covariance_valid ? 1 : 0;
      FillBoolFlags(fit_products.hbt_error_valid, hbt_error_valid);
      FillBoolFlags(fit_products.directional_error_valid, directional_error_valid);
      fit_summary_tree.Fill();

      auto &radii_per_phi = summary_points[std::make_pair(key.cent_index, key.mt_index)];
      for (std::size_t radius_index = 0; radius_index < fit_products.hbt_radii_results.size(); ++radius_index) {
        const R2SummaryPointDecision decision = DecideHbtR2SummaryPoint(fit_products.hbt_radii_results[radius_index],
                                                                        fit_products.hbt_error_valid[radius_index],
                                                                        config.projection_fit);
        RadiusSummaryPoint point;
        point.phi_center = phi_bin.Center();
        point.phi_error = 0.5 * (phi_bin.max - phi_bin.min);
        point.value = fit_products.hbt_radii_results[radius_index].r2;
        point.error = decision.error;
        point.valid = decision.write_point;
        if (decision.skipped_invalid_hbt_error) {
          ++statistics.r2_summary_points_skipped_invalid_hbt_error;
        }
        radii_per_phi[radius_index].push_back(point);
      }
    }

    output_file->cd();
    fit_summary_tree.Write();
    WriteAnalysisStatistics(*output_file, statistics);
    WriteEventPlaneMetadata(*output_file, config);
    WriteMtSlicingMetadata(*output_file, config);
    WriteProjectionFitMetadata(*output_file, config);
    WriteR2SummaryGraphs(*r2_summary_directory, config, summary_points);
    output_file->Write();

    return statistics;
  }

}  // namespace femto3d
