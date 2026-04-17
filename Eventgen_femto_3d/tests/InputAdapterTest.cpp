#include "femto3d/InputReader.h"
#include "femto3d/Workflow.h"

#include "TFile.h"
#include "TTree.h"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

using femto3d::ApplicationConfig;
using femto3d::InputSchema;
using femto3d::MakeDefaultApplicationConfig;

void Expect(const bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

template <typename Callable>
void ExpectThrows(Callable&& callable, const std::string& message) {
  bool threw = false;
  try {
    callable();
  } catch (const std::exception&) {
    threw = true;
  }
  Expect(threw, message);
}

ApplicationConfig MakeBlastwaveConfig(const std::string& input_root,
                                      const std::string& output_root) {
  ApplicationConfig config = MakeDefaultApplicationConfig();
  config.input_root_path = input_root;
  config.output_root_path = output_root;
  config.input_schema = InputSchema::kBlastwaveFlatTrees;
  config.analysis.event_plane =
      femto3d::MakeDefaultEventPlaneConfigForSchema(config.input_schema);
  config.analysis.centrality_bins = {{0.0, 100.0, "cent_0.00-100.00"}};
  config.analysis.mt_bins = {{0.0, 1.0, "femto_mt_0.00-1.00"}};
  config.analysis.phi_bins = {{-femto3d::kPi / 2.0, femto3d::kPi / 2.0, "phi_all"}};
  config.analysis.selection.target_pdg = 211;
  config.analysis.selection.femto_eta_min = -1.0;
  config.analysis.selection.femto_eta_max = 1.0;
  config.analysis.selection.femto_pt_min = 0.0;
  config.analysis.selection.femto_pt_max = 10.0;
  return config;
}

void WriteValidBlastwaveRoot(const std::string& path) {
  TFile file(path.c_str(), "RECREATE");
  if (file.IsZombie()) {
    throw std::runtime_error("Failed to create ROOT file: " + path);
  }
  file.cd();

  int event_id = 0;
  double centrality = 0.0;
  double psi2 = 0.0;
  TTree events("events", "events");
  events.SetDirectory(&file);
  events.Branch("event_id", &event_id, "event_id/I");
  events.Branch("centrality", &centrality, "centrality/D");
  events.Branch("psi2", &psi2, "psi2/D");

  event_id = 0;
  centrality = 5.0;
  psi2 = 0.25;
  events.Fill();

  event_id = 1;
  centrality = 15.0;
  psi2 = -0.20;
  events.Fill();

  int particle_event_id = 0;
  int pid = 211;
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  double mass = femto3d::kChargedPionMass;
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  double t = 0.0;
  TTree particles("particles", "particles");
  particles.SetDirectory(&file);
  particles.Branch("event_id", &particle_event_id, "event_id/I");
  particles.Branch("pid", &pid, "pid/I");
  particles.Branch("px", &px, "px/D");
  particles.Branch("py", &py, "py/D");
  particles.Branch("pz", &pz, "pz/D");
  particles.Branch("mass", &mass, "mass/D");
  particles.Branch("x", &x, "x/D");
  particles.Branch("y", &y, "y/D");
  particles.Branch("z", &z, "z/D");
  particles.Branch("t", &t, "t/D");

  particle_event_id = 0;
  px = 0.30;
  py = 0.05;
  pz = 0.10;
  x = 1.0;
  y = 0.0;
  z = 0.2;
  t = 0.0;
  particles.Fill();

  particle_event_id = 0;
  px = -0.28;
  py = 0.04;
  pz = -0.08;
  x = 1.3;
  y = 0.1;
  z = -0.1;
  t = 0.2;
  particles.Fill();

  events.Write();
  particles.Write();
  file.Close();
}

void WriteGapBlastwaveRoot(const std::string& path) {
  TFile file(path.c_str(), "RECREATE");
  if (file.IsZombie()) {
    throw std::runtime_error("Failed to create ROOT file: " + path);
  }
  file.cd();

  int event_id = 0;
  double centrality = 0.0;
  double psi2 = 0.0;
  TTree events("events", "events");
  events.SetDirectory(&file);
  events.Branch("event_id", &event_id, "event_id/I");
  events.Branch("centrality", &centrality, "centrality/D");
  events.Branch("psi2", &psi2, "psi2/D");

  event_id = 0;
  centrality = 5.0;
  psi2 = 0.1;
  events.Fill();

  event_id = 2;
  centrality = 25.0;
  psi2 = 0.3;
  events.Fill();

  events.Write();
  file.Close();
}

void WriteOrphanBlastwaveRoot(const std::string& path) {
  TFile file(path.c_str(), "RECREATE");
  if (file.IsZombie()) {
    throw std::runtime_error("Failed to create ROOT file: " + path);
  }
  file.cd();

  int event_id = 0;
  double centrality = 0.0;
  double psi2 = 0.0;
  TTree events("events", "events");
  events.SetDirectory(&file);
  events.Branch("event_id", &event_id, "event_id/I");
  events.Branch("centrality", &centrality, "centrality/D");
  events.Branch("psi2", &psi2, "psi2/D");

  event_id = 0;
  centrality = 5.0;
  psi2 = 0.1;
  events.Fill();

  int particle_event_id = 1;
  int pid = 211;
  double px = 0.2;
  double py = 0.0;
  double pz = 0.0;
  double mass = femto3d::kChargedPionMass;
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  double t = 0.0;
  TTree particles("particles", "particles");
  particles.SetDirectory(&file);
  particles.Branch("event_id", &particle_event_id, "event_id/I");
  particles.Branch("pid", &pid, "pid/I");
  particles.Branch("px", &px, "px/D");
  particles.Branch("py", &py, "py/D");
  particles.Branch("pz", &pz, "pz/D");
  particles.Branch("mass", &mass, "mass/D");
  particles.Branch("x", &x, "x/D");
  particles.Branch("y", &y, "y/D");
  particles.Branch("z", &z, "z/D");
  particles.Branch("t", &t, "t/D");
  particles.Fill();

  events.Write();
  particles.Write();
  file.Close();
}

}  // namespace

int main() {
  try {
    const std::filesystem::path temp_dir =
        std::filesystem::temp_directory_path() /
        "eventgen_femto_3d_input_adapter_test";
    std::filesystem::create_directories(temp_dir);

    const std::string valid_root = (temp_dir / "blastwave_valid.root").string();
    const std::string output_root = (temp_dir / "analysis_output.root").string();
    WriteValidBlastwaveRoot(valid_root);

    ApplicationConfig config = MakeBlastwaveConfig(valid_root, output_root);
    const std::vector<femto3d::EventData> events = femto3d::LoadEventData(config);
    Expect(events.size() == 2U, "Expected two events in blast-wave fixture.");
    Expect(events[0].particles.size() == 2U,
           "Expected first event to aggregate two particles.");
    Expect(events[1].particles.empty(),
           "Expected second event to allow zero particles.");
    Expect(std::abs(events[0].centrality - 5.0) < 1.0e-12,
           "Expected centrality to be propagated from events tree.");
    Expect(std::abs(events[0].event_plane_psi - 0.25) < 1.0e-12,
           "Expected psi2 to be propagated from events tree.");

    const femto3d::AnalysisStatistics statistics = femto3d::RunAnalysis(config);
    Expect(statistics.events_read == 2U, "Expected analysis to read two events.");
    Expect(statistics.events_with_valid_event_plane == 2U,
           "Expected both events to have valid event planes.");
    Expect(statistics.events_rejected_insufficient_femto_particles == 1U,
           "Expected one zero-particle event to be rejected after event-plane handling.");
    Expect(statistics.accepted_pairs == 1U,
           "Expected one accepted pair from the populated event.");

    TFile output_file(output_root.c_str(), "READ");
    Expect(!output_file.IsZombie(), "Expected analysis output ROOT file to be readable.");
    Expect(output_file.GetDirectory("Femto3D") != nullptr,
           "Expected Femto3D output directory.");
    Expect(output_file.GetDirectory("R2Summary") != nullptr,
           "Expected R2Summary output directory.");
    Expect(output_file.Get("analysis_statistics") != nullptr,
           "Expected analysis_statistics tree in output file.");

    const std::string gap_root = (temp_dir / "blastwave_gap.root").string();
    WriteGapBlastwaveRoot(gap_root);
    config.input_root_path = gap_root;
    ExpectThrows([&config]() { (void)femto3d::LoadEventData(config); },
                 "Expected non-contiguous event ids to fail.");

    const std::string orphan_root = (temp_dir / "blastwave_orphan.root").string();
    WriteOrphanBlastwaveRoot(orphan_root);
    config.input_root_path = orphan_root;
    ExpectThrows([&config]() { (void)femto3d::LoadEventData(config); },
                 "Expected orphan particles to fail.");
  } catch (const std::exception& error) {
    std::cerr << "input_adapter_test failed: " << error.what() << "\n";
    return 1;
  }

  return 0;
}
