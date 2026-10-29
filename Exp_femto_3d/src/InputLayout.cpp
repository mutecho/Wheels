#include "InputLayout.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "TAxis.h"
#include "exp_femto_3d/Config.h"

namespace exp_femto_3d {

  InputAxisLayout ResolveInputAxisLayout(const InputFormat format) {
    InputAxisLayout layout;
    switch (format) {
      case InputFormat::kFemtoDream:
        return layout;
      case InputFormat::kFemto:
        layout.phi = 5;
        layout.qn = 6;
        return layout;
    }
    throw ConfigError("Cannot resolve unknown input.format.");
  }

  void ValidateInputSparses(const THnSparseF &se, const THnSparseF &me, const InputConfig &config,
                            const InputAxisLayout &layout, const std::string &se_path,
                            const std::string &me_path) {
    const std::string context =
        "input.format=" + ToString(config.format) + " SE=" + se_path + " ME=" + me_path;
    if (se.GetNdimensions() != 7 || me.GetNdimensions() != 7)
      throw std::runtime_error(context + ": expected seven-dimensional THnSparseF inputs.");
    if (config.format != InputFormat::kFemto) return;
    std::array<std::string, 7> names;
    names[layout.q[0]] = "q_out";
    names[layout.q[1]] = "q_side";
    names[layout.q[2]] = "q_long";
    names[layout.mt] = "mT";
    names[layout.centrality] = "centrality";
    names[layout.phi] = "phi";
    names[layout.qn] = "q-vector";
    for (int a = 0; a < 7; ++a) {
      const auto *s = se.GetAxis(a), *m = me.GetAxis(a);
      const std::string label = context + ": axis " + std::to_string(a) + " (" + names[a] + ")";
      if (s->GetNbins() < 1 || s->GetNbins() != m->GetNbins())
        throw std::runtime_error(label + " has incompatible bin counts.");
      for (int b = 1; b <= s->GetNbins() + 1; ++b) {
        const double x = s->GetBinLowEdge(b), y = m->GetBinLowEdge(b);
        const double tolerance = 1.e-9 * std::max({1., std::abs(x), std::abs(y)});
        if (!std::isfinite(x) || !std::isfinite(y) ||
            (b > 1 && (x <= s->GetBinLowEdge(b - 1) || y <= m->GetBinLowEdge(b - 1))))
          throw std::runtime_error(label + " must have finite increasing edges.");
        if (std::abs(x - y) > tolerance) throw std::runtime_error(label + " has different SE/ME bin edges.");
      }
    }
  }

  void ValidateSymmetricPhiDomain(const double low, const double high, const bool mapped) {
    if (mapped && (!std::isfinite(low) || !std::isfinite(high) || low < -1.e-9 || high > kPi + 1.e-9))
      throw std::runtime_error(
          "map_pair_phi_to_symmetric_range requires raw phi within [0, pi]; "
          "disable it for [0, 2*pi] input. No periodic merging is performed.");
  }

  TAxis ResolvePhiSelectionAxis(const TAxis &source, const PhiMode mode) {
    if (mode == PhiMode::kSeparate) return source;
    const int n = source.GetNbins(), half = n / 2;
    const double low = source.GetBinLowEdge(1), high = source.GetBinUpEdge(n);
    // O2 JSON configurations commonly round pi to six decimal places. Allow
    // that endpoint precision, but require paired bin edges to match tightly.
    if (n % 2 != 0 || half < 1 || std::abs(low) > 1.e-5 || std::abs(high - 2 * kPi) > 1.e-5)
      throw std::runtime_error("build.phi_mode=fold_pi requires an even, paired full 0..2pi phi axis.");
    const double period = .5 * (high - low);
    std::vector<double> edges;
    for (int b = 1; b <= half + 1; ++b) {
      const double edge = source.GetBinLowEdge(b);
      if (std::abs(source.GetBinLowEdge(b + half) - edge - period) > 1.e-9)
        throw std::runtime_error("build.phi_mode=fold_pi requires matching phi and phi+pi bin edges.");
      edges.push_back(edge);
    }
    return TAxis(half, edges.data());
  }

  std::unique_ptr<TH3D> ProjectInput3D(THnSparseF &source, const InputAxisLayout &layout,
                                       const PhiMode mode) {
    auto *phi = source.GetAxis(layout.phi);
    const int first = phi->GetFirst(), last = phi->GetLast();
    auto result = std::unique_ptr<TH3D>(source.Projection(layout.q[0], layout.q[1], layout.q[2]));
    result->SetDirectory(nullptr);
    if (mode == PhiMode::kSeparate) return result;
    const int half = phi->GetNbins() / 2;
    struct RestoreRange {
      TAxis *axis;
      int first, last;
      ~RestoreRange() { axis->SetRange(first, last); }
    } restore{phi, first, last};
    phi->SetRange(first + half, last + half);
    auto partner = std::unique_ptr<TH3D>(source.Projection(layout.q[0], layout.q[1], layout.q[2]));
    partner->SetDirectory(nullptr);
    result->Add(partner.get());
    return result;
  }

}  // namespace exp_femto_3d
