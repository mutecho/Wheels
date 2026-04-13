#include "femto3d/Workflow.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TList.h"
#include "TPaveText.h"
#include "TText.h"
#include "TTree.h"

#include <cmath>
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

void WriteRichLegacyRoot(const std::string& path) {
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

  for (int event_index = 0; event_index < 8; ++event_index) {
    const double shift = 0.05 * static_cast<double>(event_index);
    const double time_shift = 0.02 * static_cast<double>(event_index);
    pdg = {211, 211, 211, 211};
    px = {0.30, -0.27, 0.26, -0.24};
    py = {0.02, 0.03, -0.04, -0.03};
    pz = {0.08, -0.07, 0.05, -0.06};
    mass = {femto3d::kChargedPionMass,
            femto3d::kChargedPionMass,
            femto3d::kChargedPionMass,
            femto3d::kChargedPionMass};
    x = {1.0 + shift, 1.3 + shift, 0.8 + shift, 1.1 + shift};
    y = {0.0, 0.1 + 0.01 * event_index, -0.2, 0.15 - 0.01 * event_index};
    z = {0.2, -0.1, 0.15, -0.05};
    t = {0.0 + time_shift, 0.2 + time_shift, 0.1 + time_shift, 0.25 + time_shift};
    events.Fill();
  }

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
  config.analysis.centrality_bins = {{0.0, 100.0, "cent_all"}};
  config.analysis.mt_bins = {{0.0, 1.0, "mt_all"}};
  config.analysis.phi_bins = {{-femto3d::kPi / 2.0, femto3d::kPi / 2.0, "phi_all"}};
  config.analysis.selection.target_pdg = 211;
  config.analysis.selection.femto_eta_min = -1.0;
  config.analysis.selection.femto_eta_max = 1.0;
  config.analysis.selection.femto_pt_min = 0.0;
  config.analysis.selection.femto_pt_max = 10.0;
  config.analysis.projection_fit.accept_forced_posdef_covariance_as_valid = true;
  return config;
}

std::string CollectTextLines(TPaveText& text_box) {
  std::string collected;
  TIter next(text_box.GetListOfLines());
  while (TObject* object = next()) {
    auto* text = dynamic_cast<TText*>(object);
    if (text == nullptr) {
      continue;
    }
    if (!collected.empty()) {
      collected.push_back('\n');
    }
    collected += text->GetTitle();
  }
  return collected;
}

}  // namespace

int main() {
  try {
    const std::filesystem::path temp_dir =
        std::filesystem::temp_directory_path() /
        "eventgen_femto_3d_r2summary_visibility_test";
    std::filesystem::create_directories(temp_dir);

    const std::string input_root = (temp_dir / "legacy_rich.root").string();
    const std::string output_root = (temp_dir / "legacy_r2summary.root").string();
    WriteRichLegacyRoot(input_root);

    const ApplicationConfig config = MakeLegacyConfig(input_root, output_root);
    const femto3d::AnalysisStatistics statistics = femto3d::RunAnalysis(config);
    Expect(statistics.accepted_pairs >= 40U,
           "Expected rich legacy fixture to produce enough accepted pairs.");

    TFile output_file(output_root.c_str(), "READ");
    Expect(!output_file.IsZombie(), "Expected readable R2 summary output ROOT file.");

    const std::string base_path = "R2Summary/cent_all/mt_all/";
    const std::vector<std::string> graph_names = {
        "Rout2_vs_phi", "Rside2_vs_phi", "Rlong2_vs_phi",
        "Ros2_vs_phi",  "Rol2_vs_phi",   "Rsl2_vs_phi"};
    for (const std::string& graph_name : graph_names) {
      auto* graph = dynamic_cast<TGraphErrors*>(
          output_file.Get((base_path + graph_name).c_str()));
      Expect(graph != nullptr, "Expected visible R2 summary graph: " + graph_name);
      Expect(graph->GetN() == 1, "Expected one summary point in graph: " + graph_name);
      Expect(graph->GetMarkerStyle() == 20,
             "Expected persisted marker style on graph: " + graph_name);
      Expect(std::abs(graph->GetMarkerSize() - 1.2) < 1.0e-12,
             "Expected persisted marker size on graph: " + graph_name);
      Expect(graph->GetLineWidth() == 2,
             "Expected persisted line width on graph: " + graph_name);

      auto* canvas = dynamic_cast<TCanvas*>(
          output_file.Get((base_path + graph_name + "_canvas").c_str()));
      Expect(canvas != nullptr, "Expected summary canvas for graph: " + graph_name);
      TList* primitives = canvas->GetListOfPrimitives();
      Expect(primitives != nullptr,
             "Expected canvas primitive list for graph: " + graph_name);
      Expect(primitives->FindObject(graph_name.c_str()) != nullptr,
             "Expected canvas to contain drawn graph: " + graph_name);

      auto* info_box = dynamic_cast<TPaveText*>(
          primitives->FindObject((graph_name + "_summary_box").c_str()));
      Expect(info_box != nullptr,
             "Expected canvas to contain summary info box: " + graph_name);
      const std::string info_box_text = CollectTextLines(*info_box);
      Expect(info_box_text.find("Points") != std::string::npos,
             "Expected summary info box to report point count.");
      Expect(info_box_text.find("Mean(Y)") != std::string::npos,
             "Expected summary info box to report mean value.");
      Expect(info_box_text.find("RMS(Y)") != std::string::npos,
             "Expected summary info box to report RMS value.");
    }
  } catch (const std::exception& error) {
    std::cerr << "r2_summary_visibility_test failed: " << error.what() << "\n";
    return 1;
  }

  return 0;
}
