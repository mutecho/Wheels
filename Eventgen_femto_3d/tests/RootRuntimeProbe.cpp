#include "TFile.h"
#include "TTree.h"

#include <filesystem>
#include <string>

int main() {
  const std::filesystem::path temp_root =
      std::filesystem::temp_directory_path() /
      "eventgen_femto_3d_root_runtime_probe.root";

  {
    TFile output(temp_root.string().c_str(), "RECREATE");
    if (output.IsZombie()) {
      return 1;
    }

    double value = 1.23;
    TTree tree("events", "events");
    tree.Branch("value", &value, "value/D");
    tree.Fill();
    tree.Write();
  }

  {
    TFile input(temp_root.string().c_str(), "READ");
    if (input.IsZombie()) {
      std::filesystem::remove(temp_root);
      return 2;
    }

    auto* tree = dynamic_cast<TTree*>(input.Get("events"));
    if (tree == nullptr || tree->GetEntries() != 1) {
      std::filesystem::remove(temp_root);
      return 3;
    }
  }

  std::filesystem::remove(temp_root);
  return 0;
}
