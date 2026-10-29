#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <THnBase.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

// Independent real-input check: ROOT's own Projection is compared with the
// program's manual sparse accumulation. Inputs and comparison output are READ.
void validate_comparison() {
  try {
    const char* files[] = {
      "/Users/allenzhou/ALICE/alidata/hyperloop_res/femtoep/OO/AnalysisResults_3d_25ae_pass2_loosecuttest.root",
      "/Users/allenzhou/ALICE/alidata/hyperloop_setup/femtoep_oo/run/res_20261020/AnalysisResults.root"};
    std::unique_ptr<TFile> result(TFile::Open("comparison_25ae_vs_local_20261020/comparison.root", "READ"));
    if (!result || result->IsZombie()) throw std::runtime_error("Invalid comparison output");
    int checked = 0;
    for (int file = 0; file < 2; ++file) {
      std::unique_ptr<TFile> input(TFile::Open(files[file], "READ"));
      if (!input || input->IsZombie()) throw std::runtime_error("Invalid source input");
      for (int event = 0; event < 2; ++event) {
        const std::string label = std::string(file ? "B" : "A") + (event ? "_ME" : "_SE");
        const std::string path = std::string("femto-dream-pair-task-track-track/") +
            (event ? "MixedEvent_3Dqn/" : "SameEvent_3Dqn/") + "relPair3dRmTMultPercentileQnPairphi";
        std::unique_ptr<THnBase> source(dynamic_cast<THnBase*>(input->Get(path.c_str())));
        if (!source) throw std::runtime_error("Missing source sparse histogram");
        for (int d = 0; d < source->GetNdimensions(); ++d) source->GetAxis(d)->SetRange(1, source->GetAxis(d)->GetNbins());
        for (int d : {0, 4}) {
          std::unique_ptr<TH1> reference(source->Projection(d, "E"));
          reference->SetDirectory(nullptr);
          const auto key = std::string("spectra/inclusive/") + (d == 0 ? "dk_out/" : "centrality/") + label;
          auto* actual = dynamic_cast<TH1*>(result->Get(key.c_str()));
          if (!actual || actual->GetNbinsX() != reference->GetNbinsX()) throw std::runtime_error("Projection schema mismatch");
          for (int b = 1; b <= reference->GetNbinsX(); ++b) {
            const double expected = reference->GetBinContent(b);
            if (std::abs(actual->GetBinContent(b) - expected) > 1.e-10 * std::max(1., std::abs(expected)))
              throw std::runtime_error("ROOT projection mismatch: " + key);
          }
          ++checked;
        }
        std::unique_ptr<TH2> reference(source->Projection(1, 0, "E"));
        reference->SetDirectory(nullptr);
        auto* actual = dynamic_cast<TH2*>(result->Get(("maps/dk_out_vs_dk_side/" + label).c_str()));
        if (!actual) throw std::runtime_error("Missing 2D projection");
        for (int y = 1; y <= reference->GetNbinsY(); ++y) for (int x = 1; x <= reference->GetNbinsX(); ++x) {
          const double expected = reference->GetBinContent(x, y);
          if (std::abs(actual->GetBinContent(x, y) - expected) > 1.e-10 * std::max(1., std::abs(expected)))
            throw std::runtime_error("2D ROOT projection mismatch for " + label);
        }
        ++checked;
        std::cout << "ROOT Projection cross-check passed: " << label << std::endl;
      }
    }
    std::cout << "Validated " << checked << " independent ROOT 1D/2D projections against real inputs." << std::endl;
    gSystem->Exit(0);
  } catch (const std::exception& e) {
    std::cerr << "VALIDATION FAILED: " << e.what() << std::endl;
    gSystem->Exit(1);
  }
}
