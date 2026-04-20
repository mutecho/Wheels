#include <filesystem>
#include <memory>

#include "TFile.h"
#include "TH1.h"
#include "THnSparse.h"

int main() {
  const std::filesystem::path temp_root =
      std::filesystem::temp_directory_path() / "exp_femto_1d_root_runtime_probe.root";

  {
    TFile output(temp_root.string().c_str(), "RECREATE");
    if (output.IsZombie()) {
      return 1;
    }

    const int bins[4] = {4, 4, 4, 4};
    const double min[4] = {0.0, 0.2, 0.0, 0.0};
    const double max[4] = {0.2, 0.4, 10.0, 3.14159265358979323846};
    auto sparse = std::make_unique<THnSparseF>("sparse", "sparse", 4, bins, min, max);
    double values[4] = {0.01, 0.3, 5.0, 0.3};
    sparse->Fill(values, 1.0);
    sparse->Write("sparse");
  }

  {
    TFile input(temp_root.string().c_str(), "READ");
    if (input.IsZombie()) {
      std::filesystem::remove(temp_root);
      return 2;
    }

    auto *sparse = dynamic_cast<THnSparseF *>(input.Get("sparse"));
    if (sparse == nullptr) {
      std::filesystem::remove(temp_root);
      return 3;
    }

    auto *projection = static_cast<TH1D *>(sparse->Projection(0));
    if (projection == nullptr) {
      std::filesystem::remove(temp_root);
      return 4;
    }
    delete projection;
  }

  std::filesystem::remove(temp_root);
  return 0;
}
