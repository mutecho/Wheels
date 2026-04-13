#include "femto3d/InputReader.h"

#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"

namespace femto3d {

  namespace {

    TTree &GetRequiredTree(TFile &file, const std::string &tree_name) {
      TTree *tree = dynamic_cast<TTree *>(file.Get(tree_name.c_str()));
      if (tree == nullptr) {
        throw std::runtime_error("Failed to find input TTree: " + tree_name);
      }
      return *tree;
    }

    void ValidateBranchExists(TTree &tree, const std::string &branch_name) {
      if (tree.GetBranch(branch_name.c_str()) == nullptr) {
        throw std::runtime_error("Missing ROOT branch: " + branch_name);
      }
    }

    void ValidateScalarLeafType(TTree &tree, const std::string &branch_name, const std::string &expected_type) {
      TLeaf *leaf = tree.GetLeaf(branch_name.c_str());
      if (leaf == nullptr) {
        throw std::runtime_error("Missing ROOT leaf: " + branch_name);
      }
      if (expected_type != leaf->GetTypeName()) {
        throw std::runtime_error("ROOT branch '" + branch_name + "' has type '" + std::string(leaf->GetTypeName())
                                 + "', expected '" + expected_type + "'.");
      }
    }

    template <typename ValueType>
    void BindScalarBranch(TTree &tree, const std::string &branch_name, ValueType *value) {
      ValidateBranchExists(tree, branch_name);
      if (tree.SetBranchAddress(branch_name.c_str(), value) < 0) {
        throw std::runtime_error("Failed to bind ROOT branch: " + branch_name);
      }
    }

    template <typename ValueType>
    void BindVectorBranch(TTree &tree, const std::string &branch_name, std::vector<ValueType> **value) {
      ValidateBranchExists(tree, branch_name);
      if (tree.SetBranchAddress(branch_name.c_str(), value) < 0) {
        throw std::runtime_error("Failed to bind ROOT branch: " + branch_name);
      }
    }

    std::size_t ValidateParticleBranches(const std::vector<int> &pdg,
                                         const std::vector<double> &px,
                                         const std::vector<double> &py,
                                         const std::vector<double> &pz,
                                         const std::vector<double> &mass,
                                         const std::vector<double> &x,
                                         const std::vector<double> &y,
                                         const std::vector<double> &z,
                                         const std::vector<double> &t) {
      const std::size_t size = pdg.size();
      const bool consistent = px.size() == size && py.size() == size && pz.size() == size && mass.size() == size
                              && x.size() == size && y.size() == size && z.size() == size && t.size() == size;
      if (!consistent) {
        throw std::runtime_error("Particle branch vector sizes are inconsistent.");
      }

      return size;
    }

    ParticleState MakeParticleState(const std::vector<int> &pdg,
                                    const std::vector<double> &px,
                                    const std::vector<double> &py,
                                    const std::vector<double> &pz,
                                    const std::vector<double> &mass,
                                    const std::vector<double> &x,
                                    const std::vector<double> &y,
                                    const std::vector<double> &z,
                                    const std::vector<double> &t,
                                    const std::size_t particle_index) {
      ParticleState particle;
      particle.pdg = pdg[particle_index];
      particle.px = px[particle_index];
      particle.py = py[particle_index];
      particle.pz = pz[particle_index];
      particle.mass = mass[particle_index];
      particle.x = x[particle_index];
      particle.y = y[particle_index];
      particle.z = z[particle_index];
      particle.t = t[particle_index];
      return particle;
    }

    std::vector<EventData> LoadLegacyVectorTreeData(TFile &input_file, const LegacyInputTreeConfig &config) {
      TTree &tree = GetRequiredTree(input_file, config.tree_name);
      ValidateScalarLeafType(tree, config.centrality_branch, "Double_t");
      if (tree.GetBranch(config.event_plane_branch.c_str()) != nullptr) {
        ValidateScalarLeafType(tree, config.event_plane_branch, "Double_t");
      }

      double centrality = 0.0;
      double event_plane = 0.0;
      std::vector<int> *pdg = nullptr;
      std::vector<double> *px = nullptr;
      std::vector<double> *py = nullptr;
      std::vector<double> *pz = nullptr;
      std::vector<double> *mass = nullptr;
      std::vector<double> *x = nullptr;
      std::vector<double> *y = nullptr;
      std::vector<double> *z = nullptr;
      std::vector<double> *t = nullptr;

      BindScalarBranch(tree, config.centrality_branch, &centrality);
      const bool has_input_event_plane = tree.GetBranch(config.event_plane_branch.c_str()) != nullptr;
      if (has_input_event_plane) {
        BindScalarBranch(tree, config.event_plane_branch, &event_plane);
      }
      BindVectorBranch(tree, config.pdg_branch, &pdg);
      BindVectorBranch(tree, config.px_branch, &px);
      BindVectorBranch(tree, config.py_branch, &py);
      BindVectorBranch(tree, config.pz_branch, &pz);
      BindVectorBranch(tree, config.mass_branch, &mass);
      BindVectorBranch(tree, config.x_branch, &x);
      BindVectorBranch(tree, config.y_branch, &y);
      BindVectorBranch(tree, config.z_branch, &z);
      BindVectorBranch(tree, config.t_branch, &t);

      std::vector<EventData> events;
      events.reserve(static_cast<std::size_t>(tree.GetEntries()));

      for (Long64_t entry = 0; entry < tree.GetEntries(); ++entry) {
        tree.GetEntry(entry);
        if (pdg == nullptr || px == nullptr || py == nullptr || pz == nullptr || mass == nullptr || x == nullptr
            || y == nullptr || z == nullptr || t == nullptr) {
          throw std::runtime_error("Legacy vector-tree reader encountered a null particle branch buffer.");
        }

        const std::size_t particle_count = ValidateParticleBranches(*pdg, *px, *py, *pz, *mass, *x, *y, *z, *t);

        EventData event;
        event.centrality = centrality;
        if (has_input_event_plane && std::isfinite(event_plane)) {
          event.event_plane_psi = event_plane;
          event.has_input_event_plane = true;
        }
        event.particles.reserve(particle_count);
        for (std::size_t particle_index = 0; particle_index < particle_count; ++particle_index) {
          event.particles.push_back(MakeParticleState(*pdg, *px, *py, *pz, *mass, *x, *y, *z, *t, particle_index));
        }
        events.push_back(std::move(event));
      }

      tree.ResetBranchAddresses();
      return events;
    }

    std::vector<EventData> LoadBlastwaveFlatTreeData(TFile &input_file, const BlastwaveInputTreeConfig &config) {
      TTree &events_tree = GetRequiredTree(input_file, config.events_tree);
      ValidateScalarLeafType(events_tree, config.event_id_branch, "Int_t");
      ValidateScalarLeafType(events_tree, config.centrality_branch, "Double_t");
      ValidateScalarLeafType(events_tree, config.event_plane_branch, "Double_t");

      int event_id = 0;
      double centrality = 0.0;
      double event_plane = 0.0;
      BindScalarBranch(events_tree, config.event_id_branch, &event_id);
      BindScalarBranch(events_tree, config.centrality_branch, &centrality);
      BindScalarBranch(events_tree, config.event_plane_branch, &event_plane);

      std::vector<EventData> events;
      events.reserve(static_cast<std::size_t>(events_tree.GetEntries()));

      std::optional<int> first_event_id;
      std::optional<int> previous_event_id;
      for (Long64_t entry = 0; entry < events_tree.GetEntries(); ++entry) {
        events_tree.GetEntry(entry);
        if (first_event_id.has_value() && event_id != *previous_event_id + 1) {
          throw std::runtime_error("Blast-wave events tree requires continuous event_id values. Expected "
                                   + std::to_string(*previous_event_id + 1) + " but found " + std::to_string(event_id)
                                   + ".");
        }
        if (!first_event_id.has_value()) {
          first_event_id = event_id;
        }
        if (!std::isfinite(event_plane)) {
          throw std::runtime_error(
              "Blast-wave events tree contains a non-finite event-plane value for "
              "event_id="
              + std::to_string(event_id) + ".");
        }

        EventData event;
        event.centrality = centrality;
        event.event_plane_psi = event_plane;
        event.has_input_event_plane = true;
        events.push_back(std::move(event));
        previous_event_id = event_id;
      }

      events_tree.ResetBranchAddresses();

      const int first_valid_event_id = first_event_id.value_or(0);
      const int last_valid_event_id = previous_event_id.value_or(-1);

      TTree &particles_tree = GetRequiredTree(input_file, config.particles_tree);
      ValidateScalarLeafType(particles_tree, config.event_id_branch, "Int_t");
      ValidateScalarLeafType(particles_tree, config.pid_branch, "Int_t");
      ValidateScalarLeafType(particles_tree, config.px_branch, "Double_t");
      ValidateScalarLeafType(particles_tree, config.py_branch, "Double_t");
      ValidateScalarLeafType(particles_tree, config.pz_branch, "Double_t");
      ValidateScalarLeafType(particles_tree, config.mass_branch, "Double_t");
      ValidateScalarLeafType(particles_tree, config.x_branch, "Double_t");
      ValidateScalarLeafType(particles_tree, config.y_branch, "Double_t");
      ValidateScalarLeafType(particles_tree, config.z_branch, "Double_t");
      ValidateScalarLeafType(particles_tree, config.t_branch, "Double_t");

      int particle_event_id = 0;
      int pid = 0;
      double px = 0.0;
      double py = 0.0;
      double pz = 0.0;
      double mass = 0.0;
      double x = 0.0;
      double y = 0.0;
      double z = 0.0;
      double t = 0.0;

      BindScalarBranch(particles_tree, config.event_id_branch, &particle_event_id);
      BindScalarBranch(particles_tree, config.pid_branch, &pid);
      BindScalarBranch(particles_tree, config.px_branch, &px);
      BindScalarBranch(particles_tree, config.py_branch, &py);
      BindScalarBranch(particles_tree, config.pz_branch, &pz);
      BindScalarBranch(particles_tree, config.mass_branch, &mass);
      BindScalarBranch(particles_tree, config.x_branch, &x);
      BindScalarBranch(particles_tree, config.y_branch, &y);
      BindScalarBranch(particles_tree, config.z_branch, &z);
      BindScalarBranch(particles_tree, config.t_branch, &t);

      for (Long64_t entry = 0; entry < particles_tree.GetEntries(); ++entry) {
        particles_tree.GetEntry(entry);
        if (particle_event_id < first_valid_event_id || particle_event_id > last_valid_event_id) {
          throw std::runtime_error(
              "Blast-wave particles tree references unknown event_id=" + std::to_string(particle_event_id) + ".");
        }

        const std::size_t event_index = static_cast<std::size_t>(particle_event_id - first_valid_event_id);
        ParticleState particle;
        particle.pdg = pid;
        particle.px = px;
        particle.py = py;
        particle.pz = pz;
        particle.mass = mass;
        particle.x = x;
        particle.y = y;
        particle.z = z;
        particle.t = t;
        events[event_index].particles.push_back(particle);
      }

      particles_tree.ResetBranchAddresses();
      return events;
    }

  }  // namespace

  std::vector<EventData> LoadEventData(const ApplicationConfig &config) {
    std::unique_ptr<TFile> input_file(TFile::Open(config.input_root_path.c_str(), "READ"));
    if (input_file == nullptr || input_file->IsZombie()) {
      throw std::runtime_error("Failed to open input ROOT file: " + config.input_root_path);
    }

    switch (config.input_schema) {
      case InputSchema::kLegacyVectorTree:
        return LoadLegacyVectorTreeData(*input_file, config.analysis.input.legacy);
      case InputSchema::kBlastwaveFlatTrees:
        return LoadBlastwaveFlatTreeData(*input_file, config.analysis.input.blastwave);
    }

    throw std::runtime_error("Unsupported input schema.");
  }

}  // namespace femto3d
