#include "Cats3DModel.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TNamed.h"
#include "TROOT.h"
#include "TTree.h"
#include "exp_femto_3d/Config.h"
#ifdef EXP_FEMTO_3D_HAS_CATS3D
#include "CATS3D.h"
#endif
#include <unistd.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
using namespace exp_femto_3d;
static void Require(bool value, const char *message) {
  if (!value) {
    throw std::runtime_error(message);
  }
}
template <class F>
void Reject(F fn, const char *message) {
  bool failed = false;
  try {
    fn();
  } catch (const std::exception &) {
    failed = true;
  }
  Require(failed, message);
}
#ifdef EXP_FEMTO_3D_HAS_CATS3D
static double DirectGaussian(double *p) {
  return std::exp(-p[3] * p[3] / 64 - p[4] * p[4] / 100 - p[5] * p[5] / 144) / (std::pow(4 * kPi, 1.5) * 4 * 5 * 6);
}
#endif
int main() {
  try {
    gROOT->SetBatch(true);
    const auto scratch = std::filesystem::temp_directory_path() / ("exp_cats3d_model_" + std::to_string(getpid()));
    std::filesystem::create_directories(scratch);
    const auto file = (scratch / "lambda.root").string();
    {
      TFile f(file.c_str(), "RECREATE");
      const double edges[] = {0, .1, .3, .5}, mev[] = {0, 100, 300, 500};
      TH1D h("lambda", "", 3, edges), m("lambda_mev", "", 3, mev);
      for (int b = 1; b <= 3; ++b) {
        h.SetBinContent(b, .2 + .2 * b);
        h.SetBinError(b, .01 * b);
        m.SetBinContent(b, h.GetBinContent(b));
        m.SetBinError(b, h.GetBinError(b));
      }
      h.Write();
      m.Write();
    }
    Cats3DLambdaInput input;
    input.kind = "histogram";
    input.file = file;
    input.histogram = "lambda";
    input.axis_unit = "GeV";
    Cats3DLambda table(input);
    Require(std::abs(table.Evaluate(125) - .5) < 1e-12, "lambda interpolation incorrect");
    Require(std::abs(table.Evaluate(0) - .4) < 1e-12 && std::abs(table.Evaluate(500) - .8) < 1e-12,
            "lambda edge-center treatment incorrect");
    Require(table.Histogram()->GetBinError(2) == .02, "external histogram error was lost");
    Reject(
        [&] {
          table.Evaluate(501);
        },
        "lambda out-of-range was silently clamped");
    auto mev = input;
    mev.axis_unit = "MeV";
    mev.histogram = "lambda_mev";
    Cats3DLambda table_mev(mev);
    Require(std::abs(table.Evaluate(125) - table_mev.Evaluate(125)) < 1e-12, "lambda units inconsistent");
    auto purity = input;
    purity.meaning = "pair_purity";
    purity.weak_decay_fraction = .1;
    Cats3DLambda pure(purity);
    Require(std::abs(pure.Evaluate(500) - .648) < 1e-12, "pair purity weak correction incorrect");
    purity.weak_decay_fraction.reset();
    Reject(
        [&] {
          Cats3DLambda invalid(purity);
        },
        "missing explicit weak fraction accepted");
    auto missing = input;
    missing.histogram = "absent";
    Reject(
        [&] {
          Cats3DLambda invalid(missing);
        },
        "missing histogram accepted");
    Cats3DLambdaInput constant;
    constant.value = .7;
    TH3D reference("reference", "", 4, -.08, .08, 4, -.08, .08, 4, -.08, .08);
    reference.SetDirectory(nullptr);
    Cats3DConfig cfg;
    cfg.grid_bins = 49;
    cfg.grid_min_fm = -36;
    cfg.grid_max_fm = 36;
    cfg.with_coulomb = false;
    cfg.use_disk_cache = false;
    cfg.max_kernel_gib = .2;
    if (!Cats3DModel::Available()) {
      Reject(
          [&] {
            Cats3DModel disabled(reference, .08, cfg);
          },
          "disabled CATS3D silently fell back");
      std::cout << "External lambda and no-CATS3D explicit rejection passed\n";
      return 0;
    }
    double p[] = {1, .7, 16, 25, 36, 3, -2, 1, 2, 0};
    {
      Cats3DModel model(reference, .08, cfg);
      model.SetLambda(std::make_shared<Cats3DLambda>(constant));
      auto predicted = model.FinalHistogram(p, true, false);
      auto core = model.CoreHistogram();
      for (int x = 1; x <= 4; ++x) {
        for (int y = 1; y <= 4; ++y) {
          for (int z = 1; z <= 4; ++z) {
            const double qo = reference.GetXaxis()->GetBinCenter(x), qs = reference.GetYaxis()->GetBinCenter(y),
                         ql = reference.GetZaxis()->GetBinCenter(z);
            const double exponent =
                (16 * qo * qo + 25 * qs * qs + 36 * ql * ql + 6 * qo * qs - 4 * qo * ql + 2 * qs * ql)
                / std::pow(.1973269804, 2);
            const double exact = 1 + std::exp(-exponent);
            Require(std::abs(core->GetBinContent(x, y, z) - exact) < .0005,
                    "full-matrix plane-wave Gaussian closure failed (q/2 or source convention)");
          }
        }
      }
      const auto folds = model.FoldCount();
      p[0] = 1.1;
      p[9] = 2;
      const double q = .02, q2 = 3 * q * q, c = core->GetBinContent(3, 3, 3);
      Require(std::abs(model.Evaluate(q, q, q, p, true, true) - 1.1 * (1 + 2 * q2) * (1 + .7 * (c - 1))) < 1e-12,
              "normalization/baseline lambda mixing wrong");
      Require(model.FoldCount() == folds, "normalization/baseline unnecessarily changed source fold");
      model.SetScan(0);
      Require(std::abs(model.Evaluate(q, q, q, p, true, false) - 1.1) < 1e-12, "lambda=0 closure failed");
      model.SetScan(1);
      Require(std::abs(model.Evaluate(q, q, q, p, true, false) - 1.1 * c) < 1e-12, "lambda=1 closure failed");
      p[2] += 1e-5;
      model.Prepare(p, true);
      Require(model.FoldCount() == folds + 1, "approximate parameter cache reused a distinct source");
      p[5] = 100;
      Reject(
          [&] {
            model.Prepare(p, true);
          },
          "non-positive-definite matrix accepted");
      p[5] = 3;
      model.Prepare(p, true);
      model.SetLambda(std::make_shared<Cats3DLambda>(input));
      Reject(
          [&] {
            model.SetScan(std::nullopt, 3);
          },
          "invalid scale lambda was clipped");
    }
    p[0] = 1;
    p[2] = 16;
    p[5] = p[6] = p[7] = 0;
    p[9] = 0;
    TH3D small("small", "", 2, -.04, .04, 2, -.04, .04, 2, -.04, .04);
    small.SetDirectory(nullptr);
    auto cached = cfg;
    cached.grid_bins = 9;
    cached.grid_min_fm = -12;
    cached.grid_max_fm = 12;
    cached.with_coulomb = true;
    cached.use_disk_cache = true;
    cached.kernel_cache_directory = (scratch / "cache").string();
    std::string key;
    double first = 0;
    {
      Cats3DModel model(small, .04, cached);
      model.SetLambda(std::make_shared<Cats3DLambda>(constant));
      first = model.Evaluate(.02, .02, .02, p, true, false);
      key = model.CacheKey();
      Require(model.Contract().find("|lib=") == std::string::npos,
              "kernel provenance still binds a library binary fingerprint");
#ifdef EXP_FEMTO_3D_HAS_CATS3D
      CATS3D direct;
      direct.SetMomBinsOut(2, -40, 40);
      direct.SetMomBinsSide(2, -40, 40);
      direct.SetMomBinsLong(2, -40, 40);
      direct.SetGridX(9, -12, 12);
      direct.SetGridY(9, -12, 12);
      direct.SetGridZ(9, -12, 12);
      direct.SetPdgId(211, 211);
      direct.SetQuantumStatistics(1);
      direct.SetQ1Q2(1);
      direct.SetRedMass(.5 * 139.57);
      direct.SetMaxPw(16);
      direct.SetWfMode(CATS3D::wfCoulomb);
      direct.SetFastKernelMode(true);
      direct.SetFastKernelSubdivisions(1);
      CATSparameters pars(CATSparameters::tSource3D, 3, true);
      pars.SetParameter(0, 4);
      pars.SetParameter(1, 5);
      pars.SetParameter(2, 6);
      direct.SetAnaSource(DirectGaussian, pars);
      direct.KillTheCat();
      auto core = model.CoreHistogram();
      for (int x = 0; x < 2; ++x) {
        for (int y = 0; y < 2; ++y) {
          for (int z = 0; z < 2; ++z) {
            Require(std::abs(core->GetBinContent(x + 1, y + 1, z + 1) - direct.GetCorrFun(x, y, z)) < 1e-12,
                    "Coulomb adapter differs from direct CATS3D");
          }
        }
      }
#endif
    }
    {
      Cats3DModel model(small, .04, cached);
      model.SetLambda(std::make_shared<Cats3DLambda>(constant));
      Require(first == model.Evaluate(.02, .02, .02, p, true, false), "disk cache altered model");
    }
    auto uncached = cached;
    uncached.use_disk_cache = false;
    {
      Cats3DModel model(small, .04, uncached);
      model.SetLambda(std::make_shared<Cats3DLambda>(constant));
      Require(first == model.Evaluate(.02, .02, .02, p, true, false), "cache on/off altered model");
    }
    const auto cache = (scratch / "cache" / (key + ".root")).string();
    const auto portable_file = (scratch / "native_kernel_with_any_name.root").string();
    std::filesystem::rename(cache, portable_file);
    auto portable = cached;
    portable.kernel_cache_file = portable_file;
    {
      TFile f(portable_file.c_str(), "UPDATE");
      Require(!f.Get("ExpFemto3DContract") && !f.Get("ExpFemto3DKernelDigest"),
              "generated kernel is not the portable native CATS3D format");
      // Legacy private tags are irrelevant; only native grid/physics metadata
      // governs compatibility with another CATS3D producer.
      TNamed meta("ExpFemto3DContract", "a different binary/build fingerprint");
      meta.Write("", TObject::kOverwrite);
      TNamed digest("ExpFemto3DKernelDigest", "legacy optional marker");
      digest.Write("", TObject::kOverwrite);
    }
    {
      Cats3DModel model(small, .04, portable);
      model.SetLambda(std::make_shared<Cats3DLambda>(constant));
      Require(first == model.Evaluate(.02, .02, .02, p, true, false), "portable native file changed the model");
      Require(model.KernelCachePath() == portable_file, "explicit kernel path lost from provenance");
    }
    auto absent = portable;
    absent.kernel_cache_file = (scratch / "absent.root").string();
    Reject(
        [&] {
          Cats3DModel model(small, .04, absent);
        },
        "missing explicit file silently rebuilt");
    Require(!std::filesystem::exists(absent.kernel_cache_file), "load-only input was created");
    auto wrong_physics = portable;
    wrong_physics.with_coulomb = false;
    Reject(
        [&] {
          Cats3DModel model(small, .04, wrong_physics);
        },
        "wrong Coulomb kernel accepted");
    TH3D shifted("shifted", "", 2, -.05, .05, 2, -.05, .05, 2, -.05, .05);
    shifted.SetDirectory(nullptr);
    Reject(
        [&] {
          Cats3DModel model(shifted, .05, portable);
        },
        "wrong q centers accepted");
    {
      TFile f(portable_file.c_str(), "UPDATE");
      f.Delete("CATS3DFastKernel6DByQ;*");
      TTree tree("CATS3DFastKernel6DByQ", "damaged payload");
      unsigned index = 0, qo = 0, qs = 0, ql = 0;
      double x = 0, y = 0, z = 0;
      std::vector<float> values(9 * 9 * 9, 1.f);
      tree.Branch("qIndex", &index);
      tree.Branch("iqOut", &qo);
      tree.Branch("iqSide", &qs);
      tree.Branch("iqLong", &ql);
      tree.Branch("qOut", &x);
      tree.Branch("qSide", &y);
      tree.Branch("qLong", &z);
      tree.Branch("kernelR", &values);
      for (index = 0; index < 8; ++index) {
        qo = index / 4;
        qs = (index / 2) % 2;
        ql = index % 2;
        if (index == 7)
          values.back() = std::numeric_limits<float>::quiet_NaN();
        tree.Fill();
      }
      tree.Write();
    }
    Reject(
        [&] {
          Cats3DModel model(small, .04, portable);
        },
        "non-finite kernel payload accepted");
    auto limited = cfg;
    limited.max_kernel_gib = 1e-9;
    Reject(
        [&] {
          Cats3DModel model(reference, .08, limited);
        },
        "kernel memory limit ignored");
    auto peak_limited = cfg;
    peak_limited.max_peak_gib = .5;
    Reject(
        [&] {
          Cats3DModel model(reference, .08, peak_limited);
        },
        "peak memory limit ignored");
    TH3D production("production", "", 40, -.2, .2, 40, -.2, .2, 40, -.2, .2);
    production.SetDirectory(nullptr);
    auto production_cfg = cfg;
    production_cfg.grid_bins = 31;
    production_cfg.max_kernel_gib = 4;
    production_cfg.max_peak_gib = 6.5;
    Reject(
        [&] {
          Cats3DModel::ValidateInputGrid(production, .15, production_cfg, Cats3DLambda(constant));
        },
        "production-size peak limit did not reject before allocating the kernel");
    production_cfg.max_peak_gib = 8;
    Cats3DModel::ValidateInputGrid(production, .15, production_cfg, Cats3DLambda(constant));
    const double uneven[] = {-.08, -.01, .08}, regular[] = {-.08, 0, .08};
    TH3D nonuniform("nonuniform", "", 2, uneven, 2, regular, 2, regular);
    Reject(
        [&] {
          Cats3DModel model(nonuniform, .08, cfg);
        },
        "nonuniform q axis was silently changed");
    std::filesystem::remove_all(scratch);
    std::cout << "CATS3D full Gaussian analytic closure, direct Coulomb closure, external lambda "
                 "portable native cache and memory budgets passed\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
