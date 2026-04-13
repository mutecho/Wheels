#include "femto3d/Workflow.h"

#include "TFile.h"
#include "TTree.h"

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using femto3d::ApplicationConfig;
using femto3d::InputSchema;
using femto3d::MakeDefaultApplicationConfig;

void Expect(const bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void WriteLegacyRoot(const std::string& path) {
  TFile file(path.c_str(), "RECREATE");
  if (file.IsZombie()) {
    throw std::runtime_error("Failed to create ROOT file: " + path);
  }
  file.cd();

  double centrality = 5.0;
  double event_plane_psi = 0.15;
  std::vector<int> pdg;
  std::vector<double> px;
  std::vector<double> py;
  std::vector<double> pz;
  std::vector<double> mass;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> t;

  TTree events("events", "events");
  events.SetDirectory(&file);
  events.Branch("centrality", &centrality, "centrality/D");
  events.Branch("event_plane_psi", &event_plane_psi, "event_plane_psi/D");
  events.Branch("pdg", &pdg);
  events.Branch("px", &px);
  events.Branch("py", &py);
  events.Branch("pz", &pz);
  events.Branch("mass", &mass);
  events.Branch("x", &x);
  events.Branch("y", &y);
  events.Branch("z", &z);
  events.Branch("t", &t);

  pdg = {211, 211};
  px = {0.25, -0.22};
  py = {0.02, 0.04};
  pz = {0.08, -0.06};
  mass = {femto3d::kChargedPionMass, femto3d::kChargedPionMass};
  x = {1.0, 1.2};
  y = {0.0, 0.1};
  z = {0.2, -0.2};
  t = {0.0, 0.15};
  events.Fill();

  events.Write();
  file.Close();
}

ApplicationConfig MakeLegacyConfig(const std::string& input_root,
                                   const std::string& output_root) {
  ApplicationConfig config = MakeDefaultApplicationConfig();
  config.input_schema = InputSchema::kLegacyVectorTree;
  config.input_root_path = input_root;
  config.output_root_path = output_root;
  config.analysis.event_plane =
      femto3d::MakeDefaultEventPlaneConfigForSchema(config.input_schema);
  config.analysis.event_plane.use_internal_reconstruction = false;
  config.analysis.event_plane.fallback_to_input_branch = true;
  config.analysis.centrality_bins = {{0.0, 100.0, "cent_0_100"}};
  config.analysis.mt_bins = {{0.0, 1.0, "femto_mt_0_1"}};
  config.analysis.phi_bins = {{-femto3d::kPi / 2.0, femto3d::kPi / 2.0, "phi_all"}};
  config.analysis.selection.target_pdg = 211;
  config.analysis.selection.femto_eta_min = -1.0;
  config.analysis.selection.femto_eta_max = 1.0;
  config.analysis.selection.femto_pt_min = 0.0;
  config.analysis.selection.femto_pt_max = 10.0;
  return config;
}

}  // namespace

int main() {
  try {
    const std::filesystem::path temp_dir =
        std::filesystem::temp_directory_path() /
        "eventgen_femto_3d_legacy_workflow_smoke";
    std::filesystem::create_directories(temp_dir);

    const std::string input_root = (temp_dir / "legacy_valid.root").string();
    const std::string output_root = (temp_dir / "legacy_output.root").string();
    WriteLegacyRoot(input_root);

    const ApplicationConfig config = MakeLegacyConfig(input_root, output_root);
    const femto3d::AnalysisStatistics statistics = femto3d::RunAnalysis(config);
    Expect(statistics.events_read == 1U, "Expected one legacy event.");
    Expect(statistics.accepted_pairs == 1U, "Expected one accepted legacy pair.");

    TFile output_file(output_root.c_str(), "READ");
    Expect(!output_file.IsZombie(), "Expected readable legacy output ROOT file.");
    Expect(output_file.GetDirectory("Femto3D") != nullptr,
           "Expected Femto3D directory for legacy run.");
    Expect(output_file.GetDirectory("R2Summary") != nullptr,
           "Expected R2Summary directory for legacy run.");
    Expect(output_file.Get("fit_summary") != nullptr,
           "Expected fit_summary tree for legacy run.");
    Expect(output_file.Get("analysis_statistics") != nullptr,
           "Expected analysis_statistics tree for legacy run.");
  } catch (const std::exception& error) {
    std::cerr << "legacy_workflow_smoke_test failed: " << error.what() << "\n";
    return 1;
  }

  return 0;
}
