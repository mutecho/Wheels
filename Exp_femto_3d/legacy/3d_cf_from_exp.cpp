#include <TAxis.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <THnSparse.h>
#include <TKey.h>
#include <TLegend.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TROOT.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include "TF1.h"
#include "TF3.h"
#include "TH1.h"
#include "TH3.h"
#include "TTree.h"

using namespace std;

// This macro contains two largely independent parts:
// 1. legacy 1D helpers used for older EP-differential studies;
// 2. the current 3D CF workflow:
//      THnSparse(SE/ME) -> normalized 3D CF -> 3D Levy fit -> summary graphs.
//
// The active production entry point is `_3d_cf_from_exp()`.

//==============================================================================
// Basic Data Types And Global Constants
//==============================================================================

enum EventType { kSameEvent, kMixedEvent };

struct EPBin {
  double low1, high1;
  double low2, high2;
  string name;
};

struct Levy3DFitResult {
  string fitModel = "diag";
  string histName;
  string groupKey;
  int centLow = -1;
  int centHigh = -1;
  double mTLow = 0.;
  double mTHigh = 0.;
  double phi = 0.;
  bool isPhiIntegrated = false;
  double norm = 0.;
  double normErr = 0.;
  double lambda = 0.;
  double lambdaErr = 0.;
  double rout2 = 0.;
  double rout2Err = 0.;
  double rside2 = 0.;
  double rside2Err = 0.;
  double rlong2 = 0.;
  double rlong2Err = 0.;
  double routside2 = 0.;
  double routside2Err = 0.;
  double routlong2 = 0.;
  double routlong2Err = 0.;
  double rsidelong2 = 0.;
  double rsidelong2Err = 0.;
  double alpha = 0.;
  double alphaErr = 0.;
  double baselineQ2 = 0.;
  double baselineQ2Err = 0.;
  double chi2 = 0.;
  double edm = -1.;
  int ndf = 0;
  int status = -1;
  int minuitIstat = -1;
  bool hasOffDiagonal = false;
  bool usesCoulomb = false;
  bool usesCoreHaloLambda = true;
  bool usesQ2Baseline = false;
  bool usesPML = false;
};

struct Levy3DFitOptions {
  bool useCoulomb = false;
  bool useCoreHaloLambda = true;
  bool useQ2Baseline = false;
  bool usePML = false;
  double fitQMax = 0.15;
};

// 1D projections of 3D histograms are averaged inside |q_other| < 0.04 GeV/c.
constexpr double kProjection1DWindow = 0.04;

// hbar*c in GeV*fm. Needed because q is stored in GeV/c and R^2 is interpreted
// in fm^2, so the Levy exponent must contain (R^2 q^2)/(hbar*c)^2.
constexpr double kHbarC = 0.1973269804;

// Like-sign pion-pion Bohr radius in fm, used by the simple Gamow-factor
// approximation for Coulomb correction.
constexpr double kPiPiLikeSignBohrRadiusFm = 387.5;

// Numerical guardrails for the full 3D Levy model and the PML objective.
constexpr double kFullR2MatrixTolerance = 1e-10;
constexpr double kLevyArgumentTolerance = 1e-12;
constexpr double kInvalidFullModelCFValue = 1e6;
constexpr double kFitPenaltyValue = 1e30;

//==============================================================================
// Generic ROOT/File Utilities
//==============================================================================

bool FileExists_bool(const string &filename) {
  ifstream file(filename);
  return file.good();
}

void FileExists_warn(const string &filename) {
  if (!FileExists_bool(filename)) {
    cerr << "Error: File: " << filename << " doesn't exist!" << endl;
  }
}

unique_ptr<TFile> GetROOT_unique(const string Readpath, const string ReadFilename, const string Option) {
  string lowerOption = Option;
  transform(lowerOption.begin(), lowerOption.end(), lowerOption.begin(), [](unsigned char c) {
    return tolower(c);
  });  // 全转小写
  if (!(lowerOption == "read" || lowerOption == "write" || lowerOption == "recreate")) {
    cout << "Invalid Option: " << Option << ", should be one of: read, write, recreate" << endl;
  }
  string rfilename = Readpath + "/" + ReadFilename + ".root";  // 加上root后缀
  FileExists_warn(rfilename);

  return make_unique<TFile>(rfilename.c_str(), Option.c_str());
}

TFile *GetROOT(const string Readpath, const string ReadFilename, const string Option) {
  string lowerOption = Option;
  transform(lowerOption.begin(), lowerOption.end(), lowerOption.begin(), [](unsigned char c) {
    return tolower(c);
  });  // 全转小写
  if (!(lowerOption == "read" || lowerOption == "write" || lowerOption == "recreate" || lowerOption == "update")) {
    cout << "Invalid Option: " << Option << ", should be one of: read, write, recreate or update" << endl;
  }
  string rfilename = Readpath + "/" + ReadFilename + ".root";  // 加上root后缀
  FileExists_warn(rfilename);

  return new TFile(rfilename.c_str(), Option.c_str());
}

bool InitializeROOTFile(const string &path, const string &filename) {
  auto wf = GetROOT(path, filename, "recreate");
  if (!wf || wf->IsZombie()) {
    cout << "ERROR: cannot initialize ROOT file " << path << "/" << filename << ".root" << endl;
    delete wf;
    return false;
  }
  wf->Close();
  delete wf;
  return true;
}

std::string doubleToString(double x, int precision = 2) {
  std::stringstream ss;
  ss << std::fixed << std::setprecision(precision) << x;
  return ss.str();
}

TFile *OpenOutputROOTFileForUpdate(const string &path, const string &filename, const string &purpose) {
  auto wf = GetROOT(path, filename, "update");
  if (!wf || wf->IsZombie()) {
    cout << "ERROR: cannot " << purpose << " " << path << "/" << filename << ".root" << endl;
    delete wf;
    return nullptr;
  }
  return wf;
}

TDirectory *GetOrCreateOutputDirectory(TFile *wf, const string &dirName) {
  if (!wf) {
    return nullptr;
  }
  auto dir = wf->GetDirectory(dirName.c_str());
  if (!dir) {
    dir = wf->mkdir(dirName.c_str());
  }
  return dir;
}

//==============================================================================
// Legacy 1D CF Utilities
//==============================================================================

// Construct a 1D correlation function C(k*) from same-event and mixed-event
// histograms using a constant normalization extracted from [normLow, normHigh].
TH1D *GetCFfromSM(TH1D *h_se, TH1D *h_me, double normLow, double normHigh, double kstarMax) {
  TF1 *constant = new TF1("constant", "1", 0, 10);
  cout << "Normalizing SE and ME histograms between k* = " << normLow << " and " << normHigh << " GeV/c." << endl;
  TH1D *h_se_c = (TH1D *)h_se->Clone();
  TH1D *h_me_c = (TH1D *)h_me->Clone();
  cout << "Cloned SE and ME histograms." << endl;

  int binNorm[2];
  binNorm[0] = h_se_c->FindBin(normLow);
  binNorm[1] = h_se_c->FindBin(normHigh);
  double factorN = h_me_c->Integral(binNorm[0], binNorm[1]) / h_se_c->Integral(binNorm[0], binNorm[1]);
  h_me_c->Divide(constant, factorN);
  h_se_c->Divide(h_me_c);
  TH1D *h_cf_re = new TH1D();
  int Upperbin = h_se_c->FindBin(kstarMax);

  for (int nBin = 1; nBin <= Upperbin; nBin++) {
    double content = h_se_c->GetBinContent(nBin);
    double error = h_se_c->GetBinError(nBin);
    h_cf_re->SetBinContent(nBin, content);
    h_cf_re->SetBinError(nBin, error);
  }
  return h_cf_re;
}

TH1D *onefromndHisto(unique_ptr<TFile> rfile, string hpath, bool isReranged, double rerangeMax) {
  auto h4 = (THnSparseF *)rfile->Get(hpath.c_str());

  // 获取轴
  auto ax_kstar = h4->GetAxis(0);

  TH1D *h1 = h4->Projection(0);
  TH1D *h1c = (TH1D *)h1->Clone();

  h1c->SetTitle("CF from ndHisto");
  h1c->GetXaxis()->SetTitle("k* (GeV/c)");
  h1c->GetYaxis()->SetTitle("C(k*)");

  if (!isReranged) {
    return h1c;
  }

  // do cf
  //   TF1 *constant = new TF1("constant", "1", 0, 10);
  //   int binNorm[2];
  //   binNorm[0] = h1sc->FindBin(0.5);
  //   binNorm[1] = h1sc->FindBin(0.8);
  //   double factorN = h1mc->Integral(binNorm[0], binNorm[1]) /
  //                    h1sc->Integral(binNorm[0], binNorm[1]);
  //   h1mc->Divide(constant, factorN);
  //   h1sc->Divide(h1mc);
  //   TH1D *h1sc_re = (TH1D *)h1sc->Clone();

  int binRange[2];
  //   binRange[0] = h1sc->FindBin(0.);
  binRange[0] = h1c->FindBin(0.);
  binRange[1] = h1c->FindBin(rerangeMax);

  TH1D *h1c_re = (TH1D *)h1c->Clone();

  h1c_re->Reset();  // 清空内容

  for (int i = binRange[0]; i <= binRange[1]; i++) {
    double content = h1c->GetBinContent(i);
    double error = h1c->GetBinError(i);
    int newBin = i - binRange[0] + 1;  // 从1开始填
    h1c_re->SetBinContent(newBin, content);
    h1c_re->SetBinError(newBin, error);
  }

  // 重新设置x轴范围
  h1c_re->GetXaxis()->SetRangeUser(0., 0.25);

  h1c_re->SetTitle("reranged CF from ndHisto");
  h1c_re->GetXaxis()->SetTitle("k* (GeV/c)");
  h1c_re->GetYaxis()->SetTitle("C(k*)");
  return h1c_re;
  //   TCanvas *c1 = new TCanvas("c1", "CF from ndHisto", 800, 600);
  //   c1->SetMargin(0.13, 0.05, 0.12, 0.07);
  //   c1->SetGrid();
  //   h1sc->Draw("hist");
  //   c1->Update();

  //   wfFemtoep->cd();
  //   h1sc->Write("CF_from_ndHisto");
  //   h1sc_re->Write("CF_reranged_from_ndHisto");
  //   wfFemtoep->Close();
  //   rfFemtoep->Close();
}

TH1D *ndHistoRead(TFile *rf, string histoname, int centLow, int centHigh, double mTLow, double mTHigh, EPBin epbin) {
  // string fullpath = filepath + "/" +filename;
  THnSparseF *h4 = (THnSparseF *)rf->Get(histoname.c_str());
  //                              "relPairkstarmTMultMultPercentileQn");;
  // if (eventType == kSameEvent) {
  //   cout << "Processing Same Event ndHisto..." << endl;
  //   h4 = (THnSparseF
  //   *)rf->Get("femto-dream-pair-task-track-track/SameEvent_EP/"
  //                              "relPairkstarmTMultMultPercentileQn");
  // } else {
  //   cout << "Processing Mixed Event ndHisto..." << endl;
  //   h4 =
  //       (THnSparseF
  //       *)rf->Get("femto-dream-pair-task-track-track/MixedEvent_EP/"
  //                             "relPairkstarmTMultMultPercentileQn");
  // }

  // 获取轴
  auto ax_kstar = h4->GetAxis(0);
  auto ax_mT = h4->GetAxis(1);
  auto ax_cent = h4->GetAxis(2);
  auto ax_EP = h4->GetAxis(3);

  //   std::vector<EPBin> epBins;
  //   if (eventType == kSameEvent) {
  //     std::cout << "Using Same Event EP bins." << std::endl;
  //     epBins = {
  //         {0, TMath::Pi() / 4, 3 * TMath::Pi() / 4, TMath::Pi(), "In_plane"},
  //         {TMath::Pi() / 4, 3 * TMath::Pi() / 4, -1, -1, "Out_of_plane"}};
  //   } else {
  //     std::cout << "Using Mixed Event EP bins." << std::endl;
  //     epBins = {{0, TMath::Pi(), -1, -1, "Min bias EP"}};
  //   }

  int colorList[] = {kRed + 1, kBlue + 1, kGreen + 2, kOrange + 7};

  // 创建新的画布
  string cname = "cent=" + to_string(static_cast<int>(centLow)) + "-" + to_string(static_cast<int>(centHigh))
                 + "_mT=" + doubleToString(mTLow, 1) + "-" + doubleToString(mTHigh, 1);
  // TCanvas *c = new TCanvas(cname, cname, 800,600);
  //  c->SetMargin(0.13, 0.05, 0.12, 0.07); c->SetGrid();

  TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);

  // std::cout << cname << std::endl;

  // 限制中心度与 mT
  ax_cent->SetRangeUser(centLow, centHigh);
  ax_mT->SetRangeUser(mTLow, mTHigh);

  int iColor = 0;

  // 如果是双区间（in-plane），要合并两部分
  TH1D *h_kstar_total = nullptr;

  if (epbin.low2 > 0) {  // 有镜像部分
    ax_EP->SetRangeUser(epbin.low1, epbin.high1);
    auto h1 = (TH1D *)h4->Projection(0);
    ax_EP->SetRangeUser(epbin.low2, epbin.high2);
    auto h2 = (TH1D *)h4->Projection(0);
    h1->Add(h2);
    h_kstar_total = h1;
  } else {
    ax_EP->SetRangeUser(epbin.low1, epbin.high1);
    h_kstar_total = (TH1D *)h4->Projection(0);
  }
  // 改到这
  //  设置样式
  h_kstar_total->SetLineColor(colorList[iColor++ % 4]);
  h_kstar_total->SetLineWidth(2);
  h_kstar_total->SetTitle(cname.c_str());
  h_kstar_total->GetXaxis()->SetTitle("k* (GeV/c)");
  h_kstar_total->GetYaxis()->SetTitle("Counts");
  string Histname = cname + "_" + epbin.name;
  h_kstar_total->SetName(Histname.c_str());

  return h_kstar_total;

  // 归一化
  // if (h_kstar_total->Integral() > 0)
  //     h_kstar_total->Scale(1.0 / h_kstar_total->Integral());

  //   if (iColor == 1)
  //     h_kstar_total->Draw("hist");
  //   else
  //     h_kstar_total->Draw("hist same");

  // leg->AddEntry(h_kstar_total, ep.name, "l");
}

TH1D *calc_cf_from_sme_rerange(TH1D *h_se,
                               TH1D *h_me,
                               std::pair<double, double> normrange,
                               std::pair<double, double> kstarRange) {
  TF1 *constant = new TF1("constant", "1", 0, 10);

  h_se->Sumw2();
  h_me->Sumw2();

  TH1D *h_se_c = (TH1D *)h_se->Clone();
  TH1D *h_me_c = (TH1D *)h_me->Clone();

  int binNorm[2];
  binNorm[0] = h_se_c->FindBin(normrange.first);
  binNorm[1] = h_se_c->FindBin(normrange.second);
  double factorN = h_me_c->Integral(binNorm[0], binNorm[1]) / h_se_c->Integral(binNorm[0], binNorm[1]);
  h_me_c->Divide(constant, factorN);
  h_se_c->Divide(h_me_c);

  int binRange[2];
  binRange[0] = h_se_c->FindBin(kstarRange.first);
  binRange[1] = h_se_c->FindBin(kstarRange.second);

  int nBins = binRange[1] - binRange[0] + 1;
  TH1D *h_cf_re = new TH1D("", "", nBins, kstarRange.first, kstarRange.second);

  for (int i = binRange[0]; i <= binRange[1]; i++) {
    double content = h_se_c->GetBinContent(i);
    double error = h_se_c->GetBinError(i);
    int newBin = i - binRange[0] + 1;  // 从1开始填
    h_cf_re->SetBinContent(newBin, content);
    h_cf_re->SetBinError(newBin, error);
  }
  return h_cf_re;
}

//==============================================================================
// 3D CF Construction Helpers
//==============================================================================

TH1D *BuildProjectionXWithinWindow(TH3D *hist, const string &name, double qMax);
TH1D *BuildProjectionYWithinWindow(TH3D *hist, const string &name, double qMax);
TH1D *BuildProjectionZWithinWindow(TH3D *hist, const string &name, double qMax);

// Integrate only visible bins, excluding underflow/overflow. With "width" the
// result is a density-like integral in q-space.
double IntegralVisibleRange(TH3D *hist, bool useWidth = false) {
  if (!hist) {
    return 0.0;
  }
  return hist->Integral(1, hist->GetNbinsX(), 1, hist->GetNbinsY(), 1, hist->GetNbinsZ(), useWidth ? "width" : "");
}

// Save three 1D projections of a 3D histogram. When useWindow=true, the
// projection along one axis is the average inside the orthogonal window
// |q| < kProjection1DWindow.
void Write1DProjections(
    TH3D *hCF, TDirectory *dir, const string &baseName, const string &yTitle = "C(q)", bool useWindow = true) {
  TH1D *hProjX = nullptr;
  TH1D *hProjY = nullptr;
  TH1D *hProjZ = nullptr;

  if (useWindow) {
    hProjX = BuildProjectionXWithinWindow(hCF, baseName + "_ProjX", kProjection1DWindow);
    hProjY = BuildProjectionYWithinWindow(hCF, baseName + "_ProjY", kProjection1DWindow);
    hProjZ = BuildProjectionZWithinWindow(hCF, baseName + "_ProjZ", kProjection1DWindow);
  } else {
    hProjX = hCF->ProjectionX((baseName + "_ProjX").c_str(), 1, hCF->GetNbinsY(), 1, hCF->GetNbinsZ());
    hProjY = hCF->ProjectionY((baseName + "_ProjY").c_str(), 1, hCF->GetNbinsX(), 1, hCF->GetNbinsZ());
    hProjZ = hCF->ProjectionZ((baseName + "_ProjZ").c_str(), 1, hCF->GetNbinsX(), 1, hCF->GetNbinsY());
  }

  string hProjXName = baseName + "_ProjX";
  hProjX->SetName(hProjXName.c_str());
  hProjX->SetTitle(hProjXName.c_str());
  hProjX->GetXaxis()->SetTitle("q_{out} (GeV/c)");
  hProjX->GetYaxis()->SetTitle(yTitle.c_str());
  dir->WriteObject(hProjX, hProjXName.c_str());

  string hProjYName = baseName + "_ProjY";
  hProjY->SetName(hProjYName.c_str());
  hProjY->SetTitle(hProjYName.c_str());
  hProjY->GetXaxis()->SetTitle("q_{side} (GeV/c)");
  hProjY->GetYaxis()->SetTitle(yTitle.c_str());
  dir->WriteObject(hProjY, hProjYName.c_str());

  string hProjZName = baseName + "_ProjZ";
  hProjZ->SetName(hProjZName.c_str());
  hProjZ->SetTitle(hProjZName.c_str());
  hProjZ->GetXaxis()->SetTitle("q_{long} (GeV/c)");
  hProjZ->GetYaxis()->SetTitle(yTitle.c_str());
  dir->WriteObject(hProjZ, hProjZName.c_str());

  delete hProjX;
  delete hProjY;
  delete hProjZ;
}

// Build and save 3D correlation functions from same-event / mixed-event
// THnSparse input.
//
// Input sparse-axis convention:
//   axis 0: q_out, axis 1: q_side, axis 2: q_long,
//   axis 3: mT (or kT-like binning used by this task),
//   axis 4: centrality, axis 6: phi_pair - Psi_EP.
//
// For each requested (cent, mT, phi) slice the code computes
//   C(q) = [SE_slice(q) / ∫SE_slice(q)d^3q] / [ME_cent,mT(q) /
//   ∫ME_cent,mT(q)d^3q]
// and stores the raw SE, raw ME, the resulting 3D CF, and the 1D projections of
// the CF into the output ROOT file.
//
// Optionally, the 1D projections of the normalized SE/ME 3D histograms can also
// be written, without storing the normalized 3D histograms themselves.
//
// If mapPairPhiToSymmetricRange=true, phi bins in (pi/2, pi) are labeled as
// phi-pi, so the output naming range becomes [-pi/2, pi/2].
void CFCalc3D(string rpath,
              string rfilename,
              string taskname,
              string subtaskname_se,
              string subtaskname_me,
              string wpath,
              string wfilename,
              std::vector<std::pair<double, double>> centBins,
              std::vector<std::pair<double, double>> mTBins,
              bool mapPairPhiToSymmetricRange = true,
              bool writeNormalizedSEME1DProjections = false,
              bool reopenOutputFilePerSlice = true) {
  const bool oldAddDirectory = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  auto rf = GetROOT(rpath, rfilename, "read");
  if (!InitializeROOTFile(wpath, wfilename)) {
    if (rf) {
      rf->Close();
      delete rf;
    }
    TH1::AddDirectory(oldAddDirectory);
    return;
  }

  string hpath_se = taskname + "/" + subtaskname_se + "/relPair3dRmTMultPercentileQnPairphi";
  string hpath_me = taskname + "/" + subtaskname_me + "/relPair3dRmTMultPercentileQnPairphi";

  auto hSE_sparse_origin = (THnSparseF *)rf->Get(hpath_se.c_str());
  auto hME_sparse_origin = (THnSparseF *)rf->Get(hpath_me.c_str());

  if (!hSE_sparse_origin || !hME_sparse_origin) {
    cout << "ERROR: Cannot find THnSparse!" << endl;
    if (rf) {
      rf->Close();
      delete rf;
    }
    TH1::AddDirectory(oldAddDirectory);
    return;
  }
  auto hSE_sparse = (THnSparseF *)hSE_sparse_origin->Clone();
  auto hSE_sparse_allphi = (THnSparseF *)hSE_sparse_origin->Clone();
  auto hME_sparse = (THnSparseF *)hME_sparse_origin->Clone();
  hSE_sparse->Sumw2();
  hSE_sparse_allphi->Sumw2();
  hME_sparse->Sumw2();

  auto ax_phi = hSE_sparse->GetAxis(6);
  int nPhiBins = ax_phi->GetNbins();
  TFile *sharedOutputFile = nullptr;
  if (!reopenOutputFilePerSlice) {
    sharedOutputFile = OpenOutputROOTFileForUpdate(wpath, wfilename, "reuse output ROOT file during CF building");
    if (!sharedOutputFile) {
      rf->Close();
      delete hSE_sparse;
      delete hSE_sparse_allphi;
      delete hME_sparse;
      delete rf;
      TH1::AddDirectory(oldAddDirectory);
      return;
    }
  }

  cout << "Starting 3D CF calculation..." << endl;

  for (auto [centLow, centHigh] : centBins) {
    for (auto [mTLow, mTHigh] : mTBins) {
      string baseName = "cent=" + to_string((int)centLow) + "-" + to_string((int)centHigh)
                        + "_mT=" + doubleToString(mTLow, 2) + "-" + doubleToString(mTHigh, 2);

      cout << "Processing " << baseName << endl;

      // 设置 cut（只设置一次）
      hSE_sparse->GetAxis(4)->SetRangeUser(centLow, centHigh);
      hSE_sparse->GetAxis(3)->SetRangeUser(mTLow, mTHigh);

      hSE_sparse_allphi->GetAxis(4)->SetRangeUser(centLow, centHigh);
      hSE_sparse_allphi->GetAxis(3)->SetRangeUser(mTLow, mTHigh);

      hME_sparse->GetAxis(4)->SetRangeUser(centLow, centHigh);
      hME_sparse->GetAxis(3)->SetRangeUser(mTLow, mTHigh);

      auto hME_raw = (TH3D *)hME_sparse->Projection(0, 1, 2);
      hME_raw->SetDirectory(nullptr);
      auto hME_norm = (TH3D *)hME_raw->Clone((baseName + "_ME_norm_tmp").c_str());
      hME_norm->SetDirectory(nullptr);
      const double intME = IntegralVisibleRange(hME_norm, true);
      if (intME == 0.0) {
        cout << "WARNING: zero mixed-event integral for " << baseName << endl;
        delete hME_raw;
        delete hME_norm;
        continue;
      }
      hME_norm->Scale(1.0 / intME);

      auto writeStoredSlice = [&](TH3D *hSE_raw, TH3D *hSE_norm, const string &hname) {
        string hSERawName = hname + "_SE_raw3d";
        string hMERawName = hname + "_ME_raw3d";

        hSE_raw->SetName(hSERawName.c_str());
        hSE_raw->SetTitle(hSERawName.c_str());

        auto hMERawStored = (TH3D *)hME_raw->Clone(hMERawName.c_str());
        hMERawStored->SetDirectory(nullptr);
        hMERawStored->SetTitle(hMERawName.c_str());

        auto hCF = (TH3D *)hSE_norm->Clone(hname.c_str());
        hCF->SetDirectory(nullptr);
        hCF->Divide(hME_norm);

        hCF->SetName(hname.c_str());
        hCF->SetTitle(hname.c_str());
        hCF->GetXaxis()->SetTitle("q_{out} (GeV/c)");
        hCF->GetYaxis()->SetTitle("q_{side} (GeV/c)");
        hCF->GetZaxis()->SetTitle("q_{long} (GeV/c)");

        TFile *wf = sharedOutputFile;
        const bool ownsOutputFile = (wf == nullptr);
        if (ownsOutputFile) {
          wf = OpenOutputROOTFileForUpdate(wpath, wfilename, "update output ROOT file");
        }
        if (!wf || wf->IsZombie()) {
          delete hCF;
          delete hMERawStored;
          return;
        }

        auto dir = GetOrCreateOutputDirectory(wf, hname);
        if (!dir) {
          cout << "ERROR: cannot create output directory " << hname << endl;
          if (ownsOutputFile) {
            wf->Close();
            delete wf;
          }
          delete hCF;
          delete hMERawStored;
          return;
        }
        dir->WriteObject(hSE_raw, hSERawName.c_str());
        dir->WriteObject(hMERawStored, hMERawName.c_str());
        if (writeNormalizedSEME1DProjections) {
          Write1DProjections(hSE_norm, dir, hname + "_SE_norm3d", "Normalized density", true);
          Write1DProjections(hME_norm, dir, hname + "_ME_norm3d", "Normalized density", true);
        }
        dir->WriteObject(hCF, hname.c_str());
        Write1DProjections(hCF, dir, hname, "C(q)", true);
        if (ownsOutputFile) {
          wf->Close();
          delete wf;
        } else {
          wf->cd();
        }
        delete hCF;
        delete hMERawStored;
      };

      auto hSE_all_raw = (TH3D *)hSE_sparse_allphi->Projection(0, 1, 2);
      hSE_all_raw->SetDirectory(nullptr);
      auto hSE_all = (TH3D *)hSE_all_raw->Clone((baseName + "_SE_all_norm_tmp").c_str());
      hSE_all->SetDirectory(nullptr);
      const double intSEAll = IntegralVisibleRange(hSE_all, true);
      if (intSEAll == 0.0) {
        cout << "WARNING: zero same-event integral for " << baseName << " phi=all" << endl;
      } else {
        hSE_all->Scale(1.0 / intSEAll);
        writeStoredSlice(hSE_all_raw, hSE_all, baseName + "_phi=all_CF3D");
      }
      delete hSE_all_raw;
      delete hSE_all;

      // 遍历 pair-phi bin；可选地把 (pi/2, pi) 映射到 (-pi/2, 0)
      for (int i = 1; i <= nPhiBins; ++i) {
        double phi = ax_phi->GetBinCenter(i);

        // -------- SE --------
        hSE_sparse->GetAxis(6)->SetRange(i, i);
        auto hSE_a_raw = (TH3D *)hSE_sparse->Projection(0, 1, 2);
        hSE_a_raw->SetDirectory(nullptr);
        auto hSE_a = (TH3D *)hSE_a_raw->Clone((baseName + "_SE_slice_norm_tmp").c_str());
        hSE_a->SetDirectory(nullptr);

        // -------- 归一化（关键：分别归一）--------
        double intSE = IntegralVisibleRange(hSE_a, true);

        if (intSE == 0.0) {
          cout << "WARNING: zero integral, skip bin " << i << endl;
          delete hSE_a_raw;
          delete hSE_a;
          continue;
        }

        hSE_a->Scale(1.0 / intSE);

        // -------- 命名并写出 --------
        double phi_mapped = phi;
        if (mapPairPhiToSymmetricRange && phi > TMath::Pi() / 2.) {
          phi_mapped -= TMath::Pi();
        }
        string hname = baseName + "_phi=" + doubleToString(phi_mapped, 2) + "_CF3D";
        writeStoredSlice(hSE_a_raw, hSE_a, hname);

        delete hSE_a_raw;
        delete hSE_a;
      }

      delete hME_raw;
      delete hME_norm;
    }
  }

  if (sharedOutputFile) {
    sharedOutputFile->Close();
    delete sharedOutputFile;
  }
  rf->Close();
  delete hSE_sparse;
  delete hSE_sparse_allphi;
  delete hME_sparse;
  delete rf;
  cout << "3D CF results have been written to " << wpath << "/" << wfilename << ".root" << endl;
  TH1::AddDirectory(oldAddDirectory);
}

//==============================================================================
// 3D Levy-Fit I/O And Model Helpers
//==============================================================================

bool EndsWith(const string &value, const string &suffix) {
  if (value.size() < suffix.size()) {
    return false;
  }
  return value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// Parse histogram names produced by CFCalc3D:
//   cent=A-B_mT=C-D_phi=E_CF3D
//   cent=A-B_mT=C-D_phi=all_CF3D
// The parsed values are later reused both for bin selection and for the
// phi-dependent summary plots.
bool ParseCF3DHistogramName(const string &histName, Levy3DFitResult &result) {
  int centLow = 0;
  int centHigh = 0;
  double mTLow = 0.;
  double mTHigh = 0.;
  double phi = 0.;

  int matched =
      sscanf(histName.c_str(), "cent=%d-%d_mT=%lf-%lf_phi=%lf_CF3D", &centLow, &centHigh, &mTLow, &mTHigh, &phi);
  if (matched == 5) {
    result.histName = histName;
    result.centLow = centLow;
    result.centHigh = centHigh;
    result.mTLow = mTLow;
    result.mTHigh = mTHigh;
    result.phi = phi;
    result.isPhiIntegrated = false;
    result.groupKey = "cent=" + to_string(centLow) + "-" + to_string(centHigh) + "_mT=" + doubleToString(mTLow, 2) + "-"
                      + doubleToString(mTHigh, 2);
    return true;
  }

  int matchedAll = sscanf(histName.c_str(), "cent=%d-%d_mT=%lf-%lf_phi=all_CF3D", &centLow, &centHigh, &mTLow, &mTHigh);
  if (matchedAll != 4) {
    return false;
  }

  result.histName = histName;
  result.centLow = centLow;
  result.centHigh = centHigh;
  result.mTLow = mTLow;
  result.mTHigh = mTHigh;
  result.phi = -1.0;
  result.isPhiIntegrated = true;
  result.groupKey = "cent=" + to_string(centLow) + "-" + to_string(centHigh) + "_mT=" + doubleToString(mTLow, 2) + "-"
                    + doubleToString(mTHigh, 2);
  return true;
}

// Read one stored 3D CF histogram from the output of CFCalc3D.
// The code supports both the older flat layout and the current directory-based
// layout in which each slice lives inside a folder with the same name.
TH3D *LoadStoredCFHistogram(TFile *rf, const string &objectName) {
  if (!rf) {
    return nullptr;
  }

  if (auto hCF = dynamic_cast<TH3D *>(rf->Get(objectName.c_str()))) {
    auto hClone = (TH3D *)hCF->Clone((objectName + "_data3d").c_str());
    hClone->SetDirectory(nullptr);
    return hClone;
  }

  auto dir = dynamic_cast<TDirectory *>(rf->Get(objectName.c_str()));
  if (!dir) {
    return nullptr;
  }

  auto hCF = dynamic_cast<TH3D *>(dir->Get(objectName.c_str()));
  if (!hCF) {
    return nullptr;
  }

  auto hClone = (TH3D *)hCF->Clone((objectName + "_data3d").c_str());
  hClone->SetDirectory(nullptr);
  return hClone;
}

TH3D *LoadStoredHistogramWithSuffix(TFile *rf, const string &objectName, const string &suffix) {
  if (!rf) {
    return nullptr;
  }

  const string fullName = objectName + suffix;
  if (auto hist = dynamic_cast<TH3D *>(rf->Get(fullName.c_str()))) {
    auto clone = (TH3D *)hist->Clone((fullName + "_clone").c_str());
    clone->SetDirectory(nullptr);
    return clone;
  }

  auto dir = dynamic_cast<TDirectory *>(rf->Get(objectName.c_str()));
  if (!dir) {
    return nullptr;
  }

  auto hist = dynamic_cast<TH3D *>(dir->Get(fullName.c_str()));
  if (!hist) {
    return nullptr;
  }

  auto clone = (TH3D *)hist->Clone((fullName + "_clone").c_str());
  clone->SetDirectory(nullptr);
  return clone;
}

string BuildFitOptionTag(const Levy3DFitOptions &options) {
  return string(options.useCoreHaloLambda ? "lam" : "nolam") + "_" + (options.useCoulomb ? "coul" : "nocoul") + "_"
         + (options.useQ2Baseline ? "baseq2" : "nobaseq2") + "_" + (options.usePML ? "pml" : "chi2");
}

// Point-like like-sign pi-pi Coulomb correction in the Gamow approximation.
// Here |q| is used as an approximate q_inv and k* ≈ |q|/2.
double ComputeLikeSignPiPiGamowFactor(double qOut, double qSide, double qLong) {
  const double qMagnitude = std::sqrt(qOut * qOut + qSide * qSide + qLong * qLong);  // Approximate q_inv.
  const double kStarFm = 0.5 * qMagnitude / kHbarC;
  if (kStarFm <= 1e-12) {
    return 0.0;
  }

  const double eta = 1.0 / (kStarFm * kPiPiLikeSignBohrRadiusFm);
  const double twoPiEta = 2.0 * TMath::Pi() * eta;
  if (twoPiEta > 700.0) {
    return 0.0;
  }

  const double denominator = std::exp(twoPiEta) - 1.0;
  if (denominator <= 0.0) {
    return 0.0;
  }

  const double gamow = twoPiEta / denominator;
  return std::max(0.0, std::min(gamow, 1.0));
}

double ComputeBowlerSinyukovLikeSignPiPiValue(double norm,
                                              double lambda,
                                              double levyExponent,
                                              bool useCoulomb,
                                              bool useCoreHaloLambda,
                                              double qOut,
                                              double qSide,
                                              double qLong) {
  const double lambdaEff = useCoreHaloLambda ? lambda : 1.0;
  const double coulombFactor = useCoulomb ? ComputeLikeSignPiPiGamowFactor(qOut, qSide, qLong) : 1.0;
  const double quantumStatTerm = std::exp(-levyExponent);
  return norm * ((1.0 - lambdaEff) + lambdaEff * coulombFactor * (1.0 + quantumStatTerm));
}

double ComputeQ2BaselineFactor(double qOut, double qSide, double qLong, double baselineQ2, bool useQ2Baseline) {
  if (!useQ2Baseline) {
    return 1.0;
  }
  const double q2 = qOut * qOut + qSide * qSide + qLong * qLong;
  return 1.0 + baselineQ2 * q2;
}

bool IsFullR2MatrixPositiveSemiDefinite(double rout2,
                                        double rside2,
                                        double rlong2,
                                        double routside2,
                                        double routlong2,
                                        double rsidelong2,
                                        double tolerance = kFullR2MatrixTolerance) {
  if (rout2 < -tolerance || rside2 < -tolerance || rlong2 < -tolerance) {
    return false;
  }

  const double detOutSide = rout2 * rside2 - routside2 * routside2;
  const double detOutLong = rout2 * rlong2 - routlong2 * routlong2;
  const double detSideLong = rside2 * rlong2 - rsidelong2 * rsidelong2;
  if (detOutSide < -tolerance || detOutLong < -tolerance || detSideLong < -tolerance) {
    return false;
  }

  const double determinant = rout2 * (rside2 * rlong2 - rsidelong2 * rsidelong2)
                             - routside2 * (routside2 * rlong2 - routlong2 * rsidelong2)
                             + routlong2 * (routside2 * rsidelong2 - routlong2 * rside2);
  return determinant >= -tolerance;
}

bool HasValidFullR2MatrixFromParameterArray(const double *par) {
  if (!par) {
    return false;
  }
  return IsFullR2MatrixPositiveSemiDefinite(par[2], par[3], par[4], par[5], par[6], par[7]);
}

double EvaluateDiagonalLevyCF(double qOut,
                              double qSide,
                              double qLong,
                              double norm,
                              double lambda,
                              double rout2,
                              double rside2,
                              double rlong2,
                              double alpha,
                              double baselineQ2,
                              const Levy3DFitOptions &fitOptions) {
  const double qOut2 = qOut * qOut;
  const double qSide2 = qSide * qSide;
  const double qLong2 = qLong * qLong;
  const double argument = (rout2 * qOut2 + rside2 * qSide2 + rlong2 * qLong2) / (kHbarC * kHbarC);
  const double levyExponent = std::pow(std::max(argument, 0.0), alpha / 2.0);
  const double femtoValue = ComputeBowlerSinyukovLikeSignPiPiValue(
      norm, lambda, levyExponent, fitOptions.useCoulomb, fitOptions.useCoreHaloLambda, qOut, qSide, qLong);
  const double baselineFactor = ComputeQ2BaselineFactor(qOut, qSide, qLong, baselineQ2, fitOptions.useQ2Baseline);
  return femtoValue * baselineFactor;
}

double EvaluateFullLevyCF(double qOut,
                          double qSide,
                          double qLong,
                          double norm,
                          double lambda,
                          double rout2,
                          double rside2,
                          double rlong2,
                          double routside2,
                          double routlong2,
                          double rsidelong2,
                          double alpha,
                          double baselineQ2,
                          const Levy3DFitOptions &fitOptions) {
  if (!IsFullR2MatrixPositiveSemiDefinite(rout2, rside2, rlong2, routside2, routlong2, rsidelong2)) {
    return kInvalidFullModelCFValue;
  }
  const double argument =
      (rout2 * qOut * qOut + rside2 * qSide * qSide + rlong2 * qLong * qLong + 2.0 * routside2 * qOut * qSide
       + 2.0 * routlong2 * qOut * qLong + 2.0 * rsidelong2 * qSide * qLong)
      / (kHbarC * kHbarC);
  if (argument < -kLevyArgumentTolerance) {
    return kInvalidFullModelCFValue;
  }
  const double protectedArgument = argument < 0.0 ? 0.0 : argument;
  const double levyExponent = std::pow(protectedArgument, alpha / 2.0);
  const double femtoValue = ComputeBowlerSinyukovLikeSignPiPiValue(
      norm, lambda, levyExponent, fitOptions.useCoulomb, fitOptions.useCoreHaloLambda, qOut, qSide, qLong);
  const double baselineFactor = ComputeQ2BaselineFactor(qOut, qSide, qLong, baselineQ2, fitOptions.useQ2Baseline);
  return femtoValue * baselineFactor;
}

// Diagonal Bertsch-Pratt Levy fit model:
//   C(q) = N * [(1-lambda_eff) + lambda_eff * K_coul(q)
//                * (1 + exp(-(A(q))^(alpha/2)))]
// with
//   A(q) = (Rout^2 qout^2 + Rside^2 qside^2 + Rlong^2 qlong^2)/(hbar*c)^2.
//
// If core-halo is disabled, lambda_eff is fixed to 1.
// If Coulomb is disabled, K_coul(q) is fixed to 1.
double Levy3DModel(double *x, double *par) {
  Levy3DFitOptions fitOptions;
  fitOptions.useQ2Baseline = par[7] > 0.5;
  fitOptions.useCoulomb = par[8] > 0.5;
  fitOptions.useCoreHaloLambda = par[9] > 0.5;
  const double qOut = x[0];
  const double qSide = x[1];
  const double qLong = x[2];
  return EvaluateDiagonalLevyCF(qOut, qSide, qLong, par[0], par[1], par[2], par[3], par[4], par[5], par[6], fitOptions);
}

// Full symmetric 3D Levy model with six independent R^2 matrix elements:
//   A(q) = [Rout^2 qout^2 + Rside^2 qside^2 + Rlong^2 qlong^2
//           + 2 Routside^2 qout qside
//           + 2 Routlong^2 qout qlong
//           + 2 Rsidelong^2 qside qlong] / (hbar*c)^2 .
double Levy3DFullModel(double *x, double *par) {
  Levy3DFitOptions fitOptions;
  fitOptions.useQ2Baseline = par[10] > 0.5;
  fitOptions.useCoulomb = par[11] > 0.5;
  fitOptions.useCoreHaloLambda = par[12] > 0.5;
  const double qOut = x[0];
  const double qSide = x[1];
  const double qLong = x[2];
  return EvaluateFullLevyCF(
      qOut, qSide, qLong, par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9], fitOptions);
}

// Build the TF3 used for the diagonal Levy fit.
// Parameter convention:
//   [0]=Norm, [1]=lambda, [2]=Rout^2, [3]=Rside^2, [4]=Rlong^2, [5]=alpha,
//   [6]=b_{q^2}.
TF3 *BuildLevyFitFunction(const string &funcName, const Levy3DFitOptions &fitOptions) {
  const double q2Max = 3.0 * fitOptions.fitQMax * fitOptions.fitQMax;
  const double baselineMin = q2Max > 0.0 ? -0.9 / q2Max : -10.0;
  const double baselineMax = q2Max > 0.0 ? 2.0 / q2Max : 10.0;
  auto fitFunc = new TF3(funcName.c_str(),
                         Levy3DModel,
                         -fitOptions.fitQMax,
                         fitOptions.fitQMax,
                         -fitOptions.fitQMax,
                         fitOptions.fitQMax,
                         -fitOptions.fitQMax,
                         fitOptions.fitQMax,
                         10);
  fitFunc->SetParName(0, "Norm");
  fitFunc->SetParName(1, "lambda");
  fitFunc->SetParName(2, "Rout2");
  fitFunc->SetParName(3, "Rside2");
  fitFunc->SetParName(4, "Rlong2");
  fitFunc->SetParName(5, "alpha");
  fitFunc->SetParName(6, "BaselineQ2");
  fitFunc->SetParName(7, "UseQ2Baseline");
  fitFunc->SetParName(8, "UseCoulomb");
  fitFunc->SetParName(9, "UseCoreHaloLambda");

  fitFunc->SetParameters(1.0, 0.5, 25.0, 25.0, 25.0, 1.5, 0.0, 0.0, 0.0, 1.0);
  fitFunc->SetParLimits(0, 0.5, 1.5);
  fitFunc->SetParLimits(1, 0.0, 1.0);
  fitFunc->SetParLimits(2, 0.01, 400.0);
  fitFunc->SetParLimits(3, 0.01, 400.0);
  fitFunc->SetParLimits(4, 0.01, 400.0);
  fitFunc->SetParLimits(5, 0.5, 2.0);
  fitFunc->SetParLimits(6, baselineMin, baselineMax);
  fitFunc->FixParameter(7, fitOptions.useQ2Baseline ? 1.0 : 0.0);
  fitFunc->FixParameter(8, fitOptions.useCoulomb ? 1.0 : 0.0);
  fitFunc->FixParameter(9, fitOptions.useCoreHaloLambda ? 1.0 : 0.0);
  if (!fitOptions.useQ2Baseline) {
    fitFunc->FixParameter(6, 0.0);
  }
  if (!fitOptions.useCoreHaloLambda) {
    fitFunc->FixParameter(1, 1.0);
  }
  fitFunc->SetNpx(60);
  fitFunc->SetNpy(60);
  fitFunc->SetNpz(60);
  return fitFunc;
}

// Build the TF3 used for the full Levy fit.
// Parameter convention:
//   [0]=Norm, [1]=lambda, [2]=Rout^2, [3]=Rside^2, [4]=Rlong^2,
//   [5]=Routside^2, [6]=Routlong^2, [7]=Rsidelong^2, [8]=alpha, [9]=b_{q^2}.
TF3 *BuildFullLevyFitFunction(const string &funcName, const Levy3DFitOptions &fitOptions) {
  const double q2Max = 3.0 * fitOptions.fitQMax * fitOptions.fitQMax;
  const double baselineMin = q2Max > 0.0 ? -0.9 / q2Max : -10.0;
  const double baselineMax = q2Max > 0.0 ? 2.0 / q2Max : 10.0;
  auto fitFunc = new TF3(funcName.c_str(),
                         Levy3DFullModel,
                         -fitOptions.fitQMax,
                         fitOptions.fitQMax,
                         -fitOptions.fitQMax,
                         fitOptions.fitQMax,
                         -fitOptions.fitQMax,
                         fitOptions.fitQMax,
                         13);
  fitFunc->SetParName(0, "Norm");
  fitFunc->SetParName(1, "lambda");
  fitFunc->SetParName(2, "Rout2");
  fitFunc->SetParName(3, "Rside2");
  fitFunc->SetParName(4, "Rlong2");
  fitFunc->SetParName(5, "Routside2");
  fitFunc->SetParName(6, "Routlong2");
  fitFunc->SetParName(7, "Rsidelong2");
  fitFunc->SetParName(8, "alpha");
  fitFunc->SetParName(9, "BaselineQ2");
  fitFunc->SetParName(10, "UseQ2Baseline");
  fitFunc->SetParName(11, "UseCoulomb");
  fitFunc->SetParName(12, "UseCoreHaloLambda");

  const double initialParameters[13] = {1.0, 0.5, 25.0, 25.0, 25.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 1.0};
  fitFunc->SetParameters(initialParameters);
  fitFunc->SetParLimits(0, 0.5, 1.5);
  fitFunc->SetParLimits(1, 0.0, 1.0);
  fitFunc->SetParLimits(2, 0.01, 400.0);
  fitFunc->SetParLimits(3, 0.01, 400.0);
  fitFunc->SetParLimits(4, 0.01, 400.0);
  fitFunc->SetParLimits(8, 0.5, 2.0);
  fitFunc->SetParLimits(9, baselineMin, baselineMax);
  fitFunc->FixParameter(10, fitOptions.useQ2Baseline ? 1.0 : 0.0);
  fitFunc->FixParameter(11, fitOptions.useCoulomb ? 1.0 : 0.0);
  fitFunc->FixParameter(12, fitOptions.useCoreHaloLambda ? 1.0 : 0.0);
  if (!fitOptions.useQ2Baseline) {
    fitFunc->FixParameter(9, 0.0);
  }
  if (!fitOptions.useCoreHaloLambda) {
    fitFunc->FixParameter(1, 1.0);
  }
  fitFunc->SetNpx(60);
  fitFunc->SetNpy(60);
  fitFunc->SetNpz(60);
  return fitFunc;
}

//==============================================================================
// Projection / Drawing Helpers For The 3D Fit
//==============================================================================

std::pair<int, int> GetAxisRangeForWindow(const TAxis *axis, double qMax) {
  const int firstBin = axis->FindBin(-qMax + 1e-9);
  const int lastBin = axis->FindBin(qMax - 1e-9);
  return {std::max(firstBin, 1), std::min(lastBin, axis->GetNbins())};
}

TH1D *BuildProjectionXWithinWindow(TH3D *hist, const string &name, double qMax) {
  auto [yMin, yMax] = GetAxisRangeForWindow(hist->GetYaxis(), qMax);
  auto [zMin, zMax] = GetAxisRangeForWindow(hist->GetZaxis(), qMax);
  auto hProj = hist->ProjectionX(name.c_str(), yMin, yMax, zMin, zMax);
  const int nWindowBins = std::max(0, yMax - yMin + 1) * std::max(0, zMax - zMin + 1);
  if (nWindowBins > 0) {
    hProj->Scale(1.0 / static_cast<double>(nWindowBins));
  }
  return hProj;
}

TH1D *BuildProjectionYWithinWindow(TH3D *hist, const string &name, double qMax) {
  auto [xMin, xMax] = GetAxisRangeForWindow(hist->GetXaxis(), qMax);
  auto [zMin, zMax] = GetAxisRangeForWindow(hist->GetZaxis(), qMax);
  auto hProj = hist->ProjectionY(name.c_str(), xMin, xMax, zMin, zMax);
  const int nWindowBins = std::max(0, xMax - xMin + 1) * std::max(0, zMax - zMin + 1);
  if (nWindowBins > 0) {
    hProj->Scale(1.0 / static_cast<double>(nWindowBins));
  }
  return hProj;
}

TH1D *BuildProjectionZWithinWindow(TH3D *hist, const string &name, double qMax) {
  auto [xMin, xMax] = GetAxisRangeForWindow(hist->GetXaxis(), qMax);
  auto [yMin, yMax] = GetAxisRangeForWindow(hist->GetYaxis(), qMax);
  auto hProj = hist->ProjectionZ(name.c_str(), xMin, xMax, yMin, yMax);
  const int nWindowBins = std::max(0, xMax - xMin + 1) * std::max(0, yMax - yMin + 1);
  if (nWindowBins > 0) {
    hProj->Scale(1.0 / static_cast<double>(nWindowBins));
  }
  return hProj;
}

void StyleProjectionHistogram(TH1D *hist, int color, int markerStyle, const string &xTitle) {
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(markerStyle);
  hist->SetMarkerSize(0.9);
  hist->SetLineWidth(2);
  hist->GetXaxis()->SetTitle(xTitle.c_str());
  hist->GetYaxis()->SetTitle("C(q)");
}

void StyleProjectionCurve(TGraph *graph, int color) {
  graph->SetLineColor(color);
  graph->SetLineWidth(3);
  graph->SetMarkerSize(0.);
}

// The fit curves shown in the 1D projection canvases are not histogram
// projections of a rebinned TF3 object. Instead, they are explicit window
// averages of the continuous TF3 model, sampled on a dense x-grid.
double GetGraphMaximum(const TGraph *graph) {
  if (!graph || graph->GetN() <= 0) {
    return 0.0;
  }
  double x = 0.0;
  double y = 0.0;
  graph->GetPoint(0, x, y);
  double maxValue = y;
  for (int i = 1; i < graph->GetN(); ++i) {
    graph->GetPoint(i, x, y);
    maxValue = std::max(maxValue, y);
  }
  return maxValue;
}

double EvaluateProjectionXWindowAverage(TF3 *fitFunc, const TH3D *referenceHist, double qOut, double qMax) {
  auto [yMin, yMax] = GetAxisRangeForWindow(referenceHist->GetYaxis(), qMax);
  auto [zMin, zMax] = GetAxisRangeForWindow(referenceHist->GetZaxis(), qMax);
  double valueSum = 0.0;
  int nPoints = 0;
  for (int iy = yMin; iy <= yMax; ++iy) {
    const double qSide = referenceHist->GetYaxis()->GetBinCenter(iy);
    for (int iz = zMin; iz <= zMax; ++iz) {
      const double qLong = referenceHist->GetZaxis()->GetBinCenter(iz);
      valueSum += fitFunc->Eval(qOut, qSide, qLong);
      ++nPoints;
    }
  }
  return nPoints > 0 ? valueSum / static_cast<double>(nPoints) : 0.0;
}

double EvaluateProjectionYWindowAverage(TF3 *fitFunc, const TH3D *referenceHist, double qSide, double qMax) {
  auto [xMin, xMax] = GetAxisRangeForWindow(referenceHist->GetXaxis(), qMax);
  auto [zMin, zMax] = GetAxisRangeForWindow(referenceHist->GetZaxis(), qMax);
  double valueSum = 0.0;
  int nPoints = 0;
  for (int ix = xMin; ix <= xMax; ++ix) {
    const double qOut = referenceHist->GetXaxis()->GetBinCenter(ix);
    for (int iz = zMin; iz <= zMax; ++iz) {
      const double qLong = referenceHist->GetZaxis()->GetBinCenter(iz);
      valueSum += fitFunc->Eval(qOut, qSide, qLong);
      ++nPoints;
    }
  }
  return nPoints > 0 ? valueSum / static_cast<double>(nPoints) : 0.0;
}

double EvaluateProjectionZWindowAverage(TF3 *fitFunc, const TH3D *referenceHist, double qLong, double qMax) {
  auto [xMin, xMax] = GetAxisRangeForWindow(referenceHist->GetXaxis(), qMax);
  auto [yMin, yMax] = GetAxisRangeForWindow(referenceHist->GetYaxis(), qMax);
  double valueSum = 0.0;
  int nPoints = 0;
  for (int ix = xMin; ix <= xMax; ++ix) {
    const double qOut = referenceHist->GetXaxis()->GetBinCenter(ix);
    for (int iy = yMin; iy <= yMax; ++iy) {
      const double qSide = referenceHist->GetYaxis()->GetBinCenter(iy);
      valueSum += fitFunc->Eval(qOut, qSide, qLong);
      ++nPoints;
    }
  }
  return nPoints > 0 ? valueSum / static_cast<double>(nPoints) : 0.0;
}

TGraph *BuildProjectionCurveXWithinWindow(
    TF3 *fitFunc, const TH3D *referenceHist, const string &name, double qMax, int nSamples = 400) {
  auto graph = new TGraph(nSamples);
  graph->SetName(name.c_str());
  graph->SetTitle(name.c_str());
  const double xMin = referenceHist->GetXaxis()->GetBinLowEdge(1);
  const double xMax = referenceHist->GetXaxis()->GetBinUpEdge(referenceHist->GetNbinsX());
  for (int i = 0; i < nSamples; ++i) {
    const double x = xMin + (xMax - xMin) * static_cast<double>(i) / (nSamples - 1);
    graph->SetPoint(i, x, EvaluateProjectionXWindowAverage(fitFunc, referenceHist, x, qMax));
  }
  return graph;
}

TGraph *BuildProjectionCurveYWithinWindow(
    TF3 *fitFunc, const TH3D *referenceHist, const string &name, double qMax, int nSamples = 400) {
  auto graph = new TGraph(nSamples);
  graph->SetName(name.c_str());
  graph->SetTitle(name.c_str());
  const double xMin = referenceHist->GetYaxis()->GetBinLowEdge(1);
  const double xMax = referenceHist->GetYaxis()->GetBinUpEdge(referenceHist->GetNbinsY());
  for (int i = 0; i < nSamples; ++i) {
    const double x = xMin + (xMax - xMin) * static_cast<double>(i) / (nSamples - 1);
    graph->SetPoint(i, x, EvaluateProjectionYWindowAverage(fitFunc, referenceHist, x, qMax));
  }
  return graph;
}

TGraph *BuildProjectionCurveZWithinWindow(
    TF3 *fitFunc, const TH3D *referenceHist, const string &name, double qMax, int nSamples = 400) {
  auto graph = new TGraph(nSamples);
  graph->SetName(name.c_str());
  graph->SetTitle(name.c_str());
  const double xMin = referenceHist->GetZaxis()->GetBinLowEdge(1);
  const double xMax = referenceHist->GetZaxis()->GetBinUpEdge(referenceHist->GetNbinsZ());
  for (int i = 0; i < nSamples; ++i) {
    const double x = xMin + (xMax - xMin) * static_cast<double>(i) / (nSamples - 1);
    graph->SetPoint(i, x, EvaluateProjectionZWindowAverage(fitFunc, referenceHist, x, qMax));
  }
  return graph;
}

string FormatParameterLine(const string &label, double value, double error, const string &unit = "") {
  std::stringstream ss;
  ss << std::fixed << std::setprecision(3) << label << " = " << value << " #pm " << error;
  if (!unit.empty()) {
    ss << " " << unit;
  }
  return ss.str();
}

string BuildFitModeTitle(const Levy3DFitResult &fitResult) {
  return fitResult.hasOffDiagonal ? "Full Levy fit" : "Diagonal Levy fit";
}

string BuildFitSwitchLine(const Levy3DFitResult &fitResult) {
  return string("Coulomb(#pi^{#pm}#pi^{#pm}): ") + (fitResult.usesCoulomb ? "on" : "off")
         + ", core-halo: " + (fitResult.usesCoreHaloLambda ? "on" : "off")
         + ", q^{2} baseline: " + (fitResult.usesQ2Baseline ? "on" : "off");
}

string DescribeCovarianceQuality(int istat) {
  switch (istat) {
    case 0:
      return "not available";
    case 1:
      return "approximate";
    case 2:
      return "forced pos-def";
    case 3:
      return "full, accurate";
    default:
      return "not applicable";
  }
}

string CovarianceQualityToken(int istat) {
  switch (istat) {
    case 0:
      return "not_available";
    case 1:
      return "approximate";
    case 2:
      return "forced_pos_def";
    case 3:
      return "full_accurate";
    default:
      return "not_applicable";
  }
}

string BuildFitStatisticLine(const Levy3DFitResult &fitResult) {
  std::stringstream statLine;
  statLine << std::fixed << std::setprecision(2) << (fitResult.usesPML ? "-2 ln L/NDF = " : "#chi^{2}/NDF = ")
           << fitResult.chi2 << "/" << fitResult.ndf;
  return statLine.str();
}

TPaveText *BuildFitParameterBox(const Levy3DFitResult &fitResult, double x1, double y1, double x2, double y2) {
  auto box = new TPaveText(x1, y1, x2, y2, "NDC");
  box->SetFillColor(0);
  box->SetFillStyle(1001);
  box->SetBorderSize(1);
  box->SetTextAlign(12);
  box->SetTextFont(42);
  double textSize = fitResult.hasOffDiagonal ? 0.024 : 0.028;
  if (fitResult.usesPML) {
    textSize = fitResult.hasOffDiagonal ? 0.022 : 0.024;
  }
  box->SetTextSize(textSize);

  box->AddText(BuildFitModeTitle(fitResult).c_str());
  box->AddText(BuildFitSwitchLine(fitResult).c_str());
  box->AddText(FormatParameterLine("N", fitResult.norm, fitResult.normErr).c_str());
  if (fitResult.usesCoreHaloLambda) {
    box->AddText(FormatParameterLine("#lambda", fitResult.lambda, fitResult.lambdaErr).c_str());
  } else {
    box->AddText("#lambda fixed = 1.000");
  }
  box->AddText(FormatParameterLine("#alpha", fitResult.alpha, fitResult.alphaErr).c_str());
  if (fitResult.usesQ2Baseline) {
    box->AddText(
        FormatParameterLine("b_{q^{2}}", fitResult.baselineQ2, fitResult.baselineQ2Err, "(GeV/c)^{-2}").c_str());
  } else {
    box->AddText("b_{q^{2}} fixed = 0.000");
  }
  box->AddText(FormatParameterLine("R_{out}^{2}", fitResult.rout2, fitResult.rout2Err, "fm^{2}").c_str());
  box->AddText(FormatParameterLine("R_{side}^{2}", fitResult.rside2, fitResult.rside2Err, "fm^{2}").c_str());
  box->AddText(FormatParameterLine("R_{long}^{2}", fitResult.rlong2, fitResult.rlong2Err, "fm^{2}").c_str());

  if (fitResult.hasOffDiagonal) {
    box->AddText(FormatParameterLine("R_{outside}^{2}", fitResult.routside2, fitResult.routside2Err, "fm^{2}").c_str());
    box->AddText(FormatParameterLine("R_{outlong}^{2}", fitResult.routlong2, fitResult.routlong2Err, "fm^{2}").c_str());
    box->AddText(
        FormatParameterLine("R_{sidelong}^{2}", fitResult.rsidelong2, fitResult.rsidelong2Err, "fm^{2}").c_str());
  }
  box->AddText((string("Fit method: ") + (fitResult.usesPML ? "PML" : "chi2")).c_str());
  box->AddText(BuildFitStatisticLine(fitResult).c_str());
  if (fitResult.usesPML) {
    std::stringstream edmLine;
    edmLine << std::scientific << std::setprecision(3) << "EDM = " << fitResult.edm;
    box->AddText(edmLine.str().c_str());
    box->AddText((string("istat = ") + std::to_string(fitResult.minuitIstat)).c_str());
    box->AddText((string("Cov quality: ") + DescribeCovarianceQuality(fitResult.minuitIstat)).c_str());
  }

  return box;
}

TCanvas *BuildProjectionCanvas(
    const string &canvasName, TH1D *hData, TGraph *gFit, const string &xTitle, const Levy3DFitResult &fitResult) {
  StyleProjectionHistogram(hData, kBlack, 20, xTitle);
  StyleProjectionCurve(gFit, kRed + 1);

  const double maxValue = std::max(hData->GetMaximum(), GetGraphMaximum(gFit));
  hData->SetMaximum(maxValue * 1.15);

  auto canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600);
  canvas->SetMargin(0.13, 0.05, 0.12, 0.07);
  canvas->SetGrid();
  hData->Draw("E1");
  gFit->Draw("L SAME");

  auto legend = new TLegend(0.62, 0.72, 0.88, 0.88);
  legend->SetBorderSize(0);
  legend->AddEntry(hData, "Data projection", "lep");
  legend->AddEntry(gFit, "Levy fit projection", "l");
  legend->Draw();

  const double y1 = fitResult.hasOffDiagonal ? (fitResult.usesPML ? 0.32 : 0.40) : (fitResult.usesPML ? 0.42 : 0.50);
  auto parameterBox = BuildFitParameterBox(fitResult, 0.16, y1, 0.58, 0.88);
  parameterBox->Draw();

  canvas->Update();
  return canvas;
}

TCanvas *Build3DComparisonCanvas(const string &canvasName,
                                 TH3D *hData,
                                 TF3 *fitFunc,
                                 const Levy3DFitResult &fitResult) {
  auto canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1400, 600);
  canvas->Divide(2, 1);
  canvas->cd(1);
  gPad->SetTheta(24);
  gPad->SetPhi(32);
  hData->Draw("BOX2Z");
  canvas->cd(2);
  gPad->SetTheta(24);
  gPad->SetPhi(32);
  fitFunc->SetTitle((fitResult.histName + "_fit_function").c_str());
  fitFunc->GetXaxis()->SetTitle("q_{out} (GeV/c)");
  fitFunc->GetYaxis()->SetTitle("q_{side} (GeV/c)");
  fitFunc->GetZaxis()->SetTitle("q_{long} (GeV/c)");
  fitFunc->SetLineColor(kRed + 1);
  fitFunc->SetLineWidth(2);
  fitFunc->Draw("ISO");
  const double y1 = fitResult.hasOffDiagonal ? (fitResult.usesPML ? 0.18 : 0.28) : (fitResult.usesPML ? 0.30 : 0.40);
  auto parameterBox = BuildFitParameterBox(fitResult, 0.12, y1, 0.58, 0.88);
  parameterBox->Draw();
  canvas->Update();
  return canvas;
}

void WriteFitResultsSummary(const string &txtPath, const vector<Levy3DFitResult> &results);
void WriteR2Graphs(TFile *wf, const vector<Levy3DFitResult> &results);

//==============================================================================
// Fit Execution Helpers
//==============================================================================

// Lightweight helper that rebuilds a single normalized 3D CF directly from
// sparse input. It is kept as a utility / fallback and is not the main
// production path once CFCalc3D output files already exist.
TH3D *BuildCFHistogramFromSparse(
    THnSparseF *hSE_sparse, THnSparseF *hME_sparse, int phiBin, int phiBinSym, const string &histName) {
  (void)phiBinSym;
  hSE_sparse->GetAxis(6)->SetRange(phiBin, phiBin);
  auto hSE_a = (TH3D *)hSE_sparse->Projection(0, 1, 2);
  hSE_a->SetDirectory(nullptr);

  auto hME_a = (TH3D *)hME_sparse->Projection(0, 1, 2);
  hME_a->SetDirectory(nullptr);

  const double intSE = IntegralVisibleRange(hSE_a, true);
  const double intME = IntegralVisibleRange(hME_a, true);
  if (intSE == 0.0 || intME == 0.0) {
    delete hSE_a;
    delete hME_a;
    return nullptr;
  }

  hSE_a->Scale(1.0 / intSE);
  hME_a->Scale(1.0 / intME);

  auto hCF = (TH3D *)hSE_a->Clone(histName.c_str());
  hCF->SetDirectory(nullptr);
  hCF->SetTitle(histName.c_str());
  hCF->Divide(hME_a);
  hCF->GetXaxis()->SetTitle("q_{out} (GeV/c)");
  hCF->GetYaxis()->SetTitle("q_{side} (GeV/c)");
  hCF->GetZaxis()->SetTitle("q_{long} (GeV/c)");

  delete hSE_a;
  delete hME_a;
  return hCF;
}

struct Levy3DPMLContext {
  TH3D *hSERaw = nullptr;
  TH3D *hMERaw = nullptr;
  bool useFullModel = false;
  Levy3DFitOptions fitOptions;
  double rawSameToMixedIntegralRatio = 1.0;
};

static Levy3DPMLContext gLevy3DPMLContext;

double EvaluateLevyModelFromParameterArray(
    double qOut, double qSide, double qLong, const double *par, bool useFullModel) {
  double x[3] = {qOut, qSide, qLong};
  return useFullModel ? Levy3DFullModel(x, const_cast<double *>(par)) : Levy3DModel(x, const_cast<double *>(par));
}

// The production CF stored by CFCalc3D is
//   C_norm(q) = [SE(q) / I_SE] / [ME(q) / I_ME],
// where I_SE and I_ME are visible-bin Integral("width") values of the slice.
//
// For a Poisson-likelihood fit performed directly on raw same/mixed counts, the
// model entering the likelihood must therefore be
//   R_raw(q) = SE(q) / ME(q) = C_norm(q) * (I_SE / I_ME).
double ComputeRawToNormalizedCFScale(TH3D *hSERaw, TH3D *hMERaw) {
  if (!hSERaw || !hMERaw) {
    return 0.0;
  }

  const double intSE = IntegralVisibleRange(hSERaw, true);
  const double intME = IntegralVisibleRange(hMERaw, true);
  if (intSE <= 0.0 || intME <= 0.0) {
    return 0.0;
  }
  return intSE / intME;
}

double ComputePMLNeg2LogLContribution(double sameCounts, double mixedCounts, double modelRatio) {
  if (sameCounts < 0.0 || mixedCounts < 0.0 || modelRatio <= 0.0 || !std::isfinite(modelRatio)) {
    return kFitPenaltyValue;
  }
  if (sameCounts == 0.0 && mixedCounts == 0.0) {
    return 0.0;
  }
  if (sameCounts == 0.0) {
    return 2.0 * mixedCounts * std::log1p(modelRatio);
  }
  if (mixedCounts == 0.0) {
    const double logTerm =
        modelRatio >= 1.0 ? std::log1p(1.0 / modelRatio) : std::log1p(modelRatio) - std::log(modelRatio);
    return 2.0 * sameCounts * logTerm;
  }

  const double totalCounts = sameCounts + mixedCounts;
  const double arg1 = modelRatio * totalCounts / (sameCounts * (modelRatio + 1.0));
  const double arg2 = totalCounts / (mixedCounts * (modelRatio + 1.0));
  if (arg1 <= 0.0 || arg2 <= 0.0 || !std::isfinite(arg1) || !std::isfinite(arg2)) {
    return kFitPenaltyValue;
  }

  return -2.0 * (sameCounts * std::log(arg1) + mixedCounts * std::log(arg2));
}

void Levy3DPMLFCN(Int_t &npar, Double_t *grad, Double_t &f, Double_t *par, Int_t flag) {
  (void)npar;
  (void)grad;
  (void)flag;
  if (!gLevy3DPMLContext.hSERaw || !gLevy3DPMLContext.hMERaw) {
    f = kFitPenaltyValue;
    return;
  }
  if (gLevy3DPMLContext.useFullModel && !HasValidFullR2MatrixFromParameterArray(par)) {
    f = kFitPenaltyValue;
    return;
  }

  const int nx = gLevy3DPMLContext.hSERaw->GetNbinsX();
  const int ny = gLevy3DPMLContext.hSERaw->GetNbinsY();
  const int nz = gLevy3DPMLContext.hSERaw->GetNbinsZ();

  double neg2LogL = 0.0;
  for (int ix = 1; ix <= nx; ++ix) {
    const double qOut = gLevy3DPMLContext.hSERaw->GetXaxis()->GetBinCenter(ix);
    if (std::abs(qOut) > gLevy3DPMLContext.fitOptions.fitQMax) {
      continue;
    }
    for (int iy = 1; iy <= ny; ++iy) {
      const double qSide = gLevy3DPMLContext.hSERaw->GetYaxis()->GetBinCenter(iy);
      if (std::abs(qSide) > gLevy3DPMLContext.fitOptions.fitQMax) {
        continue;
      }
      for (int iz = 1; iz <= nz; ++iz) {
        const double qLong = gLevy3DPMLContext.hSERaw->GetZaxis()->GetBinCenter(iz);
        if (std::abs(qLong) > gLevy3DPMLContext.fitOptions.fitQMax) {
          continue;
        }

        const double sameCounts = gLevy3DPMLContext.hSERaw->GetBinContent(ix, iy, iz);
        const double mixedCounts = gLevy3DPMLContext.hMERaw->GetBinContent(ix, iy, iz);
        if (sameCounts < 0.0 || mixedCounts < 0.0) {
          f = kFitPenaltyValue;
          return;
        }
        if (sameCounts == 0.0 && mixedCounts == 0.0) {
          continue;
        }

        double modelRatio =
            EvaluateLevyModelFromParameterArray(qOut, qSide, qLong, par, gLevy3DPMLContext.useFullModel);
        modelRatio *= gLevy3DPMLContext.rawSameToMixedIntegralRatio;
        if (modelRatio <= 0.0 || !std::isfinite(modelRatio)) {
          f = kFitPenaltyValue;
          return;
        }

        const double neg2LogLContribution = ComputePMLNeg2LogLContribution(sameCounts, mixedCounts, modelRatio);
        if (!std::isfinite(neg2LogLContribution) || neg2LogLContribution < 0.0
            || neg2LogLContribution >= kFitPenaltyValue) {
          f = kFitPenaltyValue;
          return;
        }

        neg2LogL += neg2LogLContribution;
      }
    }
  }

  f = std::isfinite(neg2LogL) ? neg2LogL : kFitPenaltyValue;
}

int CountPMLUsableBins(TH3D *hSERaw, TH3D *hMERaw, double qMax) {
  if (!hSERaw || !hMERaw) {
    return 0;
  }
  int nPoints = 0;
  for (int ix = 1; ix <= hSERaw->GetNbinsX(); ++ix) {
    const double qOut = hSERaw->GetXaxis()->GetBinCenter(ix);
    if (std::abs(qOut) > qMax) {
      continue;
    }
    for (int iy = 1; iy <= hSERaw->GetNbinsY(); ++iy) {
      const double qSide = hSERaw->GetYaxis()->GetBinCenter(iy);
      if (std::abs(qSide) > qMax) {
        continue;
      }
      for (int iz = 1; iz <= hSERaw->GetNbinsZ(); ++iz) {
        const double qLong = hSERaw->GetZaxis()->GetBinCenter(iz);
        if (std::abs(qLong) > qMax) {
          continue;
        }
        if (hSERaw->GetBinContent(ix, iy, iz) + hMERaw->GetBinContent(ix, iy, iz) > 0.0) {
          ++nPoints;
        }
      }
    }
  }
  return nPoints;
}

double EstimatePMLStepSize(int parameterIndex, bool useFullModel) {
  if ((!useFullModel && (parameterIndex >= 2 && parameterIndex <= 4))
      || (useFullModel && (parameterIndex >= 2 && parameterIndex <= 7))) {
    return 0.1;
  }
  return 0.01;
}

bool IsPMLParameterFixed(int parameterIndex, bool useFullModel, const Levy3DFitOptions &fitOptions) {
  if (!useFullModel) {
    if (parameterIndex >= 7) {
      return true;
    }
    if (parameterIndex == 1 && !fitOptions.useCoreHaloLambda) {
      return true;
    }
    if (parameterIndex == 6 && !fitOptions.useQ2Baseline) {
      return true;
    }
    return false;
  }

  if (parameterIndex >= 10) {
    return true;
  }
  if (parameterIndex == 1 && !fitOptions.useCoreHaloLambda) {
    return true;
  }
  if (parameterIndex == 9 && !fitOptions.useQ2Baseline) {
    return true;
  }
  return false;
}

void ConfigurePMLMinuit(TMinuit &minuit, TF3 *fitFunc, bool useFullModel, const Levy3DFitOptions &fitOptions) {
  const int nPar = fitFunc->GetNpar();
  Int_t ierr = 0;
  for (int i = 0; i < nPar; ++i) {
    double lower = 0.0;
    double upper = 0.0;
    fitFunc->GetParLimits(i, lower, upper);
    const double value = fitFunc->GetParameter(i);
    const double step = EstimatePMLStepSize(i, useFullModel);
    const bool isFixed = IsPMLParameterFixed(i, useFullModel, fitOptions);
    if (isFixed) {
      minuit.mnparm(i, fitFunc->GetParName(i), value, step, 0.0, 0.0, ierr);
      minuit.FixParameter(i);
    } else {
      minuit.mnparm(i, fitFunc->GetParName(i), value, step, lower, upper, ierr);
    }
  }
}

bool RunPMLFit(TF3 *fitFunc,
               TH3D *hSERaw,
               TH3D *hMERaw,
               bool useFullModel,
               const Levy3DFitOptions &fitOptions,
               double &fitStatistic,
               int &ndf,
               int &fitStatus,
               double &edm,
               int &minuitIstat) {
  if (!fitFunc || !hSERaw || !hMERaw) {
    return false;
  }

  const double rawToNormalizedScale = ComputeRawToNormalizedCFScale(hSERaw, hMERaw);
  if (rawToNormalizedScale <= 0.0 || !std::isfinite(rawToNormalizedScale)) {
    return false;
  }

  gLevy3DPMLContext.hSERaw = hSERaw;
  gLevy3DPMLContext.hMERaw = hMERaw;
  gLevy3DPMLContext.useFullModel = useFullModel;
  gLevy3DPMLContext.fitOptions = fitOptions;
  gLevy3DPMLContext.rawSameToMixedIntegralRatio = rawToNormalizedScale;

  const int nPar = fitFunc->GetNpar();
  TMinuit minuit(nPar);
  minuit.SetFCN(Levy3DPMLFCN);
  minuit.SetPrintLevel(-1);
  minuit.SetErrorDef(1.0);

  ConfigurePMLMinuit(minuit, fitFunc, useFullModel, fitOptions);

  Int_t ierr = 0;
  Double_t arglist[2];
  arglist[0] = 100000;
  arglist[1] = 0.1;
  minuit.mnexcm("MIGRAD", arglist, 2, ierr);
  Int_t migradIerr = ierr;
  arglist[0] = 0;
  minuit.mnexcm("HESSE", arglist, 1, ierr);

  Double_t fmin = 0.0;
  Double_t fedm = 0.0;
  Double_t errdef = 0.0;
  Int_t npari = 0;
  Int_t nparx = 0;
  Int_t istat = 0;
  minuit.mnstat(fmin, fedm, errdef, npari, nparx, istat);
  (void)errdef;
  (void)nparx;

  for (int i = 0; i < nPar; ++i) {
    double value = 0.0;
    double error = 0.0;
    minuit.GetParameter(i, value, error);
    fitFunc->SetParameter(i, value);
    fitFunc->SetParError(i, error);
  }

  fitStatistic = fmin;
  edm = fedm;
  ndf = std::max(0, CountPMLUsableBins(hSERaw, hMERaw, fitOptions.fitQMax) - npari);
  minuitIstat = istat;
  fitStatus = migradIerr != 0 ? -migradIerr : istat;
  return true;
}

bool FitAndWriteSingleCFHistogram(TH3D *hCF,
                                  TH3D *hSERaw,
                                  TH3D *hMERaw,
                                  Levy3DFitResult &fitResult,
                                  const string &wpath,
                                  const string &wfilename,
                                  bool useFullModel,
                                  const Levy3DFitOptions &fitOptions,
                                  TFile *sharedOutputFile = nullptr) {
  auto fillResultFromFunction =
      [&](TF3 *fitFunc, double fitStatistic, int fitNdf, int fitStatus, double fitEdm, int fitMinuitIstat) {
        fitResult.fitModel = useFullModel ? "full" : "diag";
        fitResult.hasOffDiagonal = useFullModel;
        fitResult.usesCoulomb = fitOptions.useCoulomb;
        fitResult.usesCoreHaloLambda = fitOptions.useCoreHaloLambda;
        fitResult.usesQ2Baseline = fitOptions.useQ2Baseline;
        fitResult.usesPML = fitOptions.usePML;
        fitResult.norm = fitFunc->GetParameter(0);
        fitResult.normErr = fitFunc->GetParError(0);
        fitResult.lambda = fitOptions.useCoreHaloLambda ? fitFunc->GetParameter(1) : 1.0;
        fitResult.lambdaErr = fitOptions.useCoreHaloLambda ? fitFunc->GetParError(1) : 0.0;
        fitResult.rout2 = fitFunc->GetParameter(2);
        fitResult.rout2Err = fitFunc->GetParError(2);
        fitResult.rside2 = fitFunc->GetParameter(3);
        fitResult.rside2Err = fitFunc->GetParError(3);
        fitResult.rlong2 = fitFunc->GetParameter(4);
        fitResult.rlong2Err = fitFunc->GetParError(4);
        if (useFullModel) {
          fitResult.routside2 = fitFunc->GetParameter(5);
          fitResult.routside2Err = fitFunc->GetParError(5);
          fitResult.routlong2 = fitFunc->GetParameter(6);
          fitResult.routlong2Err = fitFunc->GetParError(6);
          fitResult.rsidelong2 = fitFunc->GetParameter(7);
          fitResult.rsidelong2Err = fitFunc->GetParError(7);
          fitResult.alpha = fitFunc->GetParameter(8);
          fitResult.alphaErr = fitFunc->GetParError(8);
          fitResult.baselineQ2 = fitOptions.useQ2Baseline ? fitFunc->GetParameter(9) : 0.0;
          fitResult.baselineQ2Err = fitOptions.useQ2Baseline ? fitFunc->GetParError(9) : 0.0;
        } else {
          fitResult.alpha = fitFunc->GetParameter(5);
          fitResult.alphaErr = fitFunc->GetParError(5);
          fitResult.baselineQ2 = fitOptions.useQ2Baseline ? fitFunc->GetParameter(6) : 0.0;
          fitResult.baselineQ2Err = fitOptions.useQ2Baseline ? fitFunc->GetParError(6) : 0.0;
        }
        fitResult.chi2 = fitStatistic;
        fitResult.edm = fitOptions.usePML ? fitEdm : -1.0;
        fitResult.ndf = fitNdf;
        fitResult.status = fitStatus;
        fitResult.minuitIstat = fitOptions.usePML ? fitMinuitIstat : -1;
      };

  const string objectName = fitResult.histName;

  TF3 *fitFunc = useFullModel ? BuildFullLevyFitFunction(objectName + "_levy3d_full_fit", fitOptions)
                              : BuildLevyFitFunction(objectName + "_levy3d_fit", fitOptions);
  double fitStatistic = 0.0;
  double fitEdm = -1.0;
  int fitNdf = 0;
  int fitStatus = -1;
  int fitMinuitIstat = -1;
  bool fitSucceeded = false;
  if (fitOptions.usePML) {
    fitSucceeded = RunPMLFit(
        fitFunc, hSERaw, hMERaw, useFullModel, fitOptions, fitStatistic, fitNdf, fitStatus, fitEdm, fitMinuitIstat);
  } else {
    auto fitStatusObject = hCF->Fit(fitFunc, "RSMNQ0");
    fitStatistic = fitFunc->GetChisquare();
    fitNdf = fitFunc->GetNDF();
    fitStatus = static_cast<int>(fitStatusObject);
    fitSucceeded = true;
  }
  fillResultFromFunction(fitFunc, fitStatistic, fitNdf, fitStatus, fitEdm, fitMinuitIstat);
  if (!fitSucceeded) {
    delete fitFunc;
    return false;
  }

  const bool oldBatchMode = gROOT->IsBatch();
  gROOT->SetBatch(kTRUE);

  const string fitHistName = objectName + (useFullModel ? "_full_fit3d" : "_fit3d");
  fitFunc->SetName(fitHistName.c_str());
  fitFunc->SetTitle(fitHistName.c_str());

  auto hProjXData = BuildProjectionXWithinWindow(hCF, objectName + "_data_ProjX", kProjection1DWindow);
  auto hProjYData = BuildProjectionYWithinWindow(hCF, objectName + "_data_ProjY", kProjection1DWindow);
  auto hProjZData = BuildProjectionZWithinWindow(hCF, objectName + "_data_ProjZ", kProjection1DWindow);

  auto hProjXFit = BuildProjectionCurveXWithinWindow(
      fitFunc, hCF, objectName + (useFullModel ? "_full_fit_ProjX" : "_fit_ProjX"), kProjection1DWindow);
  auto hProjYFit = BuildProjectionCurveYWithinWindow(
      fitFunc, hCF, objectName + (useFullModel ? "_full_fit_ProjY" : "_fit_ProjY"), kProjection1DWindow);
  auto hProjZFit = BuildProjectionCurveZWithinWindow(
      fitFunc, hCF, objectName + (useFullModel ? "_full_fit_ProjZ" : "_fit_ProjZ"), kProjection1DWindow);

  auto cProjX = BuildProjectionCanvas(objectName + (useFullModel ? "_canvas_full_ProjX" : "_canvas_ProjX"),
                                      hProjXData,
                                      hProjXFit,
                                      "q_{out} (GeV/c)",
                                      fitResult);
  auto cProjY = BuildProjectionCanvas(objectName + (useFullModel ? "_canvas_full_ProjY" : "_canvas_ProjY"),
                                      hProjYData,
                                      hProjYFit,
                                      "q_{side} (GeV/c)",
                                      fitResult);
  auto cProjZ = BuildProjectionCanvas(objectName + (useFullModel ? "_canvas_full_ProjZ" : "_canvas_ProjZ"),
                                      hProjZData,
                                      hProjZFit,
                                      "q_{long} (GeV/c)",
                                      fitResult);
  auto c3D =
      Build3DComparisonCanvas(objectName + (useFullModel ? "_canvas_full_3D" : "_canvas_3D"), hCF, fitFunc, fitResult);

  TFile *wf = sharedOutputFile;
  const bool ownsOutputFile = (wf == nullptr);
  if (ownsOutputFile) {
    wf = OpenOutputROOTFileForUpdate(wpath, wfilename, "write output ROOT file");
  }
  if (!wf || wf->IsZombie()) {
    gROOT->SetBatch(oldBatchMode);
    delete cProjX;
    delete cProjY;
    delete cProjZ;
    delete c3D;
    delete hProjXData;
    delete hProjYData;
    delete hProjZData;
    delete hProjXFit;
    delete hProjYFit;
    delete hProjZFit;
    delete fitFunc;
    return false;
  }

  auto dir = GetOrCreateOutputDirectory(wf, objectName);
  if (!dir) {
    cout << "ERROR: cannot create output directory " << objectName << endl;
    if (ownsOutputFile) {
      wf->Close();
      delete wf;
    }
    gROOT->SetBatch(oldBatchMode);
    delete cProjX;
    delete cProjY;
    delete cProjZ;
    delete c3D;
    delete hProjXData;
    delete hProjYData;
    delete hProjZData;
    delete hProjXFit;
    delete hProjYFit;
    delete hProjZFit;
    delete fitFunc;
    return false;
  }
  dir->cd();

  hCF->Write();
  fitFunc->Write();
  hProjXData->Write();
  hProjYData->Write();
  hProjZData->Write();
  hProjXFit->Write();
  hProjYFit->Write();
  hProjZFit->Write();
  cProjX->Write();
  cProjY->Write();
  cProjZ->Write();
  c3D->Write();

  if (ownsOutputFile) {
    wf->Close();
    delete wf;
  } else {
    wf->cd();
  }
  delete cProjX;
  delete cProjY;
  delete cProjZ;
  delete c3D;
  delete hProjXData;
  delete hProjYData;
  delete hProjZData;
  delete hProjXFit;
  delete hProjYFit;
  delete hProjZFit;
  delete fitFunc;
  gROOT->SetBatch(oldBatchMode);
  return true;
}

// Match a stored CF histogram against the user-requested fit bin lists.
bool MatchSelectedBin(const Levy3DFitResult &fitResult,
                      const vector<pair<double, double>> &centBins,
                      const vector<pair<double, double>> &mTBins) {
  const bool centMatched = std::any_of(centBins.begin(), centBins.end(), [&](const pair<double, double> &bin) {
    return fitResult.centLow == static_cast<int>(bin.first) && fitResult.centHigh == static_cast<int>(bin.second);
  });
  if (!centMatched) {
    return false;
  }

  return std::any_of(mTBins.begin(), mTBins.end(), [&](const pair<double, double> &bin) {
    return std::abs(fitResult.mTLow - bin.first) < 1e-6 && std::abs(fitResult.mTHigh - bin.second) < 1e-6;
  });
}

// Main fit driver used in production.
//
// Input:
//   rpath/rfilename : ROOT file produced by CFCalc3D.
//   centBins/mTBins : only histograms matching these slices are fitted.
//   useFullModel    : false -> diagonal Levy, true -> full six-component Levy.
//
// Output:
//   1. a ROOT file containing the fitted TF3, 1D projection canvases, 3D
//      comparison canvases and summary TGraphErrors/TF1 objects;
//   2. a txt file containing one row per fitted CF histogram.
void FitCF3DWithSelectedBins(const string &rpath,
                             const string &rfilename,
                             const string &wpath,
                             const string &wfilename,
                             const string &txtFilename,
                             std::vector<std::pair<double, double>> centBins,
                             std::vector<std::pair<double, double>> mTBins,
                             bool useFullModel,
                             const Levy3DFitOptions &fitOptions = Levy3DFitOptions(),
                             bool reopenOutputFilePerSlice = true) {
  const bool oldAddDirectory = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  auto rf = GetROOT(rpath, rfilename, "read");
  if (!rf || rf->IsZombie()) {
    cout << "ERROR: cannot open input ROOT file " << rpath << "/" << rfilename << ".root" << endl;
    delete rf;
    TH1::AddDirectory(oldAddDirectory);
    return;
  }

  if (!InitializeROOTFile(wpath, wfilename)) {
    rf->Close();
    delete rf;
    TH1::AddDirectory(oldAddDirectory);
    return;
  }

  TFile *sharedOutputFile = nullptr;
  if (!reopenOutputFilePerSlice) {
    sharedOutputFile = OpenOutputROOTFileForUpdate(wpath, wfilename, "reuse output ROOT file during fitting");
    if (!sharedOutputFile) {
      rf->Close();
      delete rf;
      TH1::AddDirectory(oldAddDirectory);
      return;
    }
  }

  vector<Levy3DFitResult> fitResults;
  TIter nextKey(rf->GetListOfKeys());
  TKey *key = nullptr;

  while ((key = (TKey *)nextKey())) {
    string objectName = key->GetName();
    if (!EndsWith(objectName, "_CF3D")) {
      continue;
    }

    Levy3DFitResult fitResult;
    if (!ParseCF3DHistogramName(objectName, fitResult)) {
      cout << "WARNING: cannot parse histogram name: " << objectName << endl;
      continue;
    }

    if (!MatchSelectedBin(fitResult, centBins, mTBins)) {
      cout << "Skipping histogram " << objectName << " due to unmatched centrality or mT bin." << endl;
      continue;
    }

    cout << "Fitting selected histogram " << objectName << endl;

    auto hCF = LoadStoredCFHistogram(rf, objectName);
    if (!hCF) {
      cout << "WARNING: cannot read histogram " << objectName << endl;
      continue;
    }

    TH3D *hSERaw = nullptr;
    TH3D *hMERaw = nullptr;
    if (fitOptions.usePML) {
      hSERaw = LoadStoredHistogramWithSuffix(rf, objectName, "_SE_raw3d");
      hMERaw = LoadStoredHistogramWithSuffix(rf, objectName, "_ME_raw3d");
      if (!hSERaw || !hMERaw) {
        cout << "WARNING: PML fit for " << objectName << " requires stored raw SE/ME histograms (_SE_raw3d/_ME_raw3d). "
             << "Please rebuild the CF ROOT file with the current CFCalc3D." << endl;
        delete hSERaw;
        delete hMERaw;
        delete hCF;
        continue;
      }
    }

    if (FitAndWriteSingleCFHistogram(
            hCF, hSERaw, hMERaw, fitResult, wpath, wfilename, useFullModel, fitOptions, sharedOutputFile)) {
      fitResults.push_back(fitResult);
    }

    delete hSERaw;
    delete hMERaw;
    delete hCF;
  }

  TFile *wf = sharedOutputFile;
  const bool ownsSummaryFile = (wf == nullptr);
  if (ownsSummaryFile) {
    wf = OpenOutputROOTFileForUpdate(wpath, wfilename, "update output ROOT file for summary graphs");
  }
  if (!wf || wf->IsZombie()) {
    cout << "ERROR: cannot update output ROOT file " << wpath << "/" << wfilename << ".root for summary graphs" << endl;
  } else {
    WriteR2Graphs(wf, fitResults);
    if (ownsSummaryFile) {
      wf->Close();
      delete wf;
    } else {
      wf->cd();
    }
  }
  if (sharedOutputFile) {
    sharedOutputFile->Close();
    delete sharedOutputFile;
  }

  WriteFitResultsSummary(wpath + "/" + txtFilename + ".txt", fitResults);

  rf->Close();
  delete rf;
  cout << (useFullModel ? "Full 3D Levy fit" : "Diagonal 3D Levy fit") << " results have been written to " << wpath
       << "/" << wfilename << ".root and " << wpath << "/" << txtFilename << ".txt" << endl;
  TH1::AddDirectory(oldAddDirectory);
}

// Write one plain-text summary row per fitted histogram. This is the most
// convenient output for subsequent scripting / plotting outside ROOT.
void WriteFitResultsSummary(const string &txtPath, const vector<Levy3DFitResult> &results) {
  ofstream out(txtPath);
  out << "# fitModel usesCoulomb usesCoreHaloLambda usesQ2Baseline usesPML "
         "histName centLow centHigh "
         "mTLow mTHigh phi isPhiIntegrated norm normErr lambda lambdaErr "
         "Rout2 Rout2Err Rside2 Rside2Err Rlong2 Rlong2Err Routside2 "
         "Routside2Err "
         "Routlong2 Routlong2Err Rsidelong2 Rsidelong2Err alpha alphaErr "
         "baselineQ2 baselineQ2Err "
         "chi2 edm ndf status minuitIstat covarianceQuality\n";
  out << std::fixed << std::setprecision(6);

  for (const auto &result : results) {
    out << result.fitModel << " " << (result.usesCoulomb ? 1 : 0) << " " << (result.usesCoreHaloLambda ? 1 : 0) << " "
        << (result.usesQ2Baseline ? 1 : 0) << " " << (result.usesPML ? 1 : 0) << " " << result.histName << " "
        << result.centLow << " " << result.centHigh << " " << result.mTLow << " " << result.mTHigh << " " << result.phi
        << " " << (result.isPhiIntegrated ? 1 : 0) << " " << result.norm << " " << result.normErr << " "
        << result.lambda << " " << result.lambdaErr << " " << result.rout2 << " " << result.rout2Err << " "
        << result.rside2 << " " << result.rside2Err << " " << result.rlong2 << " " << result.rlong2Err << " "
        << result.routside2 << " " << result.routside2Err << " " << result.routlong2 << " " << result.routlong2Err
        << " " << result.rsidelong2 << " " << result.rsidelong2Err << " " << result.alpha << " " << result.alphaErr
        << " " << result.baselineQ2 << " " << result.baselineQ2Err << " " << result.chi2 << " " << result.edm << " "
        << result.ndf << " " << result.status << " " << result.minuitIstat << " "
        << CovarianceQualityToken(result.minuitIstat) << "\n";
  }
}

// Build R^2(phi) / alpha(phi) / lambda(phi) summary graphs from the per-slice
// fit results. The functional forms follow the usual second-harmonic EP
// expansion:
//   R^2(phi) = R0^2 + 2 R2^2 cos(2phi)     for diagonal terms and Routlong^2,
//   R^2(phi) = R0^2 + 2 R2^2 sin(2phi)     for Routside^2 and Rsidelong^2,
//   alpha(phi) = const, lambda(phi) = const.
void WriteR2Graphs(TFile *wf, const vector<Levy3DFitResult> &results) {
  map<string, vector<Levy3DFitResult>> groupedResults;
  for (const auto &result : results) {
    if (result.isPhiIntegrated) {
      continue;
    }
    groupedResults[result.groupKey].push_back(result);
  }

  auto graphDir = wf->mkdir("R2_vs_phi");
  graphDir->cd();

  for (auto &[groupKey, groupResults] : groupedResults) {
    if (groupResults.empty()) {
      continue;
    }
    sort(groupResults.begin(), groupResults.end(), [](const Levy3DFitResult &lhs, const Levy3DFitResult &rhs) {
      return lhs.phi < rhs.phi;
    });

    const bool usesMappedPhiRange =
        std::any_of(groupResults.begin(), groupResults.end(), [](const Levy3DFitResult &result) {
          return result.phi < -1e-6;
        });
    const double phiFitMin = usesMappedPhiRange ? -TMath::Pi() / 2. : 0.0;
    const double phiFitMax = usesMappedPhiRange ? TMath::Pi() / 2. : TMath::Pi();

    const int nPoints = static_cast<int>(groupResults.size());
    const bool hasOffDiagonal =
        std::any_of(groupResults.begin(), groupResults.end(), [](const Levy3DFitResult &result) {
          return result.hasOffDiagonal;
        });
    const bool usesCoreHaloLambda =
        std::any_of(groupResults.begin(), groupResults.end(), [](const Levy3DFitResult &result) {
          return result.usesCoreHaloLambda;
        });
    const bool usesQ2Baseline =
        std::any_of(groupResults.begin(), groupResults.end(), [](const Levy3DFitResult &result) {
          return result.usesQ2Baseline;
        });
    auto gRout2 = new TGraphErrors(nPoints);
    auto gRside2 = new TGraphErrors(nPoints);
    auto gRlong2 = new TGraphErrors(nPoints);
    TGraphErrors *gRoutside2 = hasOffDiagonal ? new TGraphErrors(nPoints) : nullptr;
    TGraphErrors *gRoutlong2 = hasOffDiagonal ? new TGraphErrors(nPoints) : nullptr;
    TGraphErrors *gRsidelong2 = hasOffDiagonal ? new TGraphErrors(nPoints) : nullptr;
    auto gAlpha = new TGraphErrors(nPoints);
    TGraphErrors *gLambda = usesCoreHaloLambda ? new TGraphErrors(nPoints) : nullptr;
    TGraphErrors *gBaselineQ2 = usesQ2Baseline ? new TGraphErrors(nPoints) : nullptr;

    for (int i = 0; i < nPoints; ++i) {
      const auto &result = groupResults[i];
      gRout2->SetPoint(i, result.phi, result.rout2);
      gRout2->SetPointError(i, 0., result.rout2Err);

      gRside2->SetPoint(i, result.phi, result.rside2);
      gRside2->SetPointError(i, 0., result.rside2Err);

      gRlong2->SetPoint(i, result.phi, result.rlong2);
      gRlong2->SetPointError(i, 0., result.rlong2Err);

      if (hasOffDiagonal) {
        gRoutside2->SetPoint(i, result.phi, result.routside2);
        gRoutside2->SetPointError(i, 0., result.routside2Err);

        gRoutlong2->SetPoint(i, result.phi, result.routlong2);
        gRoutlong2->SetPointError(i, 0., result.routlong2Err);

        gRsidelong2->SetPoint(i, result.phi, result.rsidelong2);
        gRsidelong2->SetPointError(i, 0., result.rsidelong2Err);
      }

      gAlpha->SetPoint(i, result.phi, result.alpha);
      gAlpha->SetPointError(i, 0., result.alphaErr);

      if (usesCoreHaloLambda) {
        gLambda->SetPoint(i, result.phi, result.lambda);
        gLambda->SetPointError(i, 0., result.lambdaErr);
      }
      if (usesQ2Baseline) {
        gBaselineQ2->SetPoint(i, result.phi, result.baselineQ2);
        gBaselineQ2->SetPointError(i, 0., result.baselineQ2Err);
      }
    }

    auto fitCosRout2 = new TF1((groupKey + "_Rout2_phi_fit").c_str(), "[0]+2.0*[1]*cos(2.0*x)", phiFitMin, phiFitMax);
    fitCosRout2->SetParNames("Rout2_0", "Rout2_2");
    fitCosRout2->SetParameters(groupResults.front().rout2, 0.0);

    auto fitCosRside2 = new TF1((groupKey + "_Rside2_phi_fit").c_str(), "[0]+2.0*[1]*cos(2.0*x)", phiFitMin, phiFitMax);
    fitCosRside2->SetParNames("Rside2_0", "Rside2_2");
    fitCosRside2->SetParameters(groupResults.front().rside2, 0.0);

    auto fitCosRlong2 = new TF1((groupKey + "_Rlong2_phi_fit").c_str(), "[0]+2.0*[1]*cos(2.0*x)", phiFitMin, phiFitMax);
    fitCosRlong2->SetParNames("Rlong2_0", "Rlong2_2");
    fitCosRlong2->SetParameters(groupResults.front().rlong2, 0.0);

    TF1 *fitSinRoutside2 = nullptr;
    TF1 *fitCosRoutlong2 = nullptr;
    TF1 *fitSinRsidelong2 = nullptr;
    if (hasOffDiagonal) {
      fitSinRoutside2 =
          new TF1((groupKey + "_Routside2_phi_fit").c_str(), "[0]+2.0*[1]*sin(2.0*x)", phiFitMin, phiFitMax);
      fitSinRoutside2->SetParNames("Routside2_0", "Routside2_2");
      fitSinRoutside2->SetParameters(groupResults.front().routside2, 0.0);

      fitCosRoutlong2 =
          new TF1((groupKey + "_Routlong2_phi_fit").c_str(), "[0]+2.0*[1]*cos(2.0*x)", phiFitMin, phiFitMax);
      fitCosRoutlong2->SetParNames("Routlong2_0", "Routlong2_2");
      fitCosRoutlong2->SetParameters(groupResults.front().routlong2, 0.0);

      fitSinRsidelong2 =
          new TF1((groupKey + "_Rsidelong2_phi_fit").c_str(), "[0]+2.0*[1]*sin(2.0*x)", phiFitMin, phiFitMax);
      fitSinRsidelong2->SetParNames("Rsidelong2_0", "Rsidelong2_2");
      fitSinRsidelong2->SetParameters(groupResults.front().rsidelong2, 0.0);
    }

    auto fitConstAlpha = new TF1((groupKey + "_alpha_phi_fit").c_str(), "[0]", phiFitMin, phiFitMax);
    fitConstAlpha->SetParName(0, "alpha0");
    fitConstAlpha->SetParameter(0, groupResults.front().alpha);

    TF1 *fitConstLambda = nullptr;
    if (usesCoreHaloLambda) {
      fitConstLambda = new TF1((groupKey + "_lambda_phi_fit").c_str(), "[0]", phiFitMin, phiFitMax);
      fitConstLambda->SetParName(0, "lambda0");
      fitConstLambda->SetParameter(0, groupResults.front().lambda);
    }

    TF1 *fitConstBaselineQ2 = nullptr;
    if (usesQ2Baseline) {
      fitConstBaselineQ2 = new TF1((groupKey + "_baselineQ2_phi_fit").c_str(), "[0]", phiFitMin, phiFitMax);
      fitConstBaselineQ2->SetParName(0, "bq2_0");
      fitConstBaselineQ2->SetParameter(0, groupResults.front().baselineQ2);
    }

    string routName = groupKey + "_Rout2_vs_phi";
    gRout2->SetName(routName.c_str());
    gRout2->SetTitle((routName
                      + ";#phi_{pair}-#Psi_{EP} (rad);R_{out}^{2} "
                        "(fm^{2})")
                         .c_str());
    gRout2->SetMarkerStyle(20);
    gRout2->SetMarkerColor(kRed + 1);
    gRout2->SetLineColor(kRed + 1);
    gRout2->Fit(fitCosRout2, "QN");
    gRout2->Write();
    fitCosRout2->Write();

    string rsideName = groupKey + "_Rside2_vs_phi";
    gRside2->SetName(rsideName.c_str());
    gRside2->SetTitle((rsideName
                       + ";#phi_{pair}-#Psi_{EP} (rad);R_{side}^{2} "
                         "(fm^{2})")
                          .c_str());
    gRside2->SetMarkerStyle(21);
    gRside2->SetMarkerColor(kBlue + 1);
    gRside2->SetLineColor(kBlue + 1);
    gRside2->Fit(fitCosRside2, "QN");
    gRside2->Write();
    fitCosRside2->Write();

    string rlongName = groupKey + "_Rlong2_vs_phi";
    gRlong2->SetName(rlongName.c_str());
    gRlong2->SetTitle((rlongName
                       + ";#phi_{pair}-#Psi_{EP} (rad);R_{long}^{2} "
                         "(fm^{2})")
                          .c_str());
    gRlong2->SetMarkerStyle(22);
    gRlong2->SetMarkerColor(kGreen + 2);
    gRlong2->SetLineColor(kGreen + 2);
    gRlong2->Fit(fitCosRlong2, "QN");
    gRlong2->Write();
    fitCosRlong2->Write();

    if (hasOffDiagonal) {
      string routsideName = groupKey + "_Routside2_vs_phi";
      gRoutside2->SetName(routsideName.c_str());
      gRoutside2->SetTitle((routsideName + ";#phi_{pair}-#Psi_{EP} (rad);R_{outside}^{2} (fm^{2})").c_str());
      gRoutside2->SetMarkerStyle(34);
      gRoutside2->SetMarkerColor(kCyan + 1);
      gRoutside2->SetLineColor(kCyan + 1);
      gRoutside2->Fit(fitSinRoutside2, "QN");
      gRoutside2->Write();
      fitSinRoutside2->Write();

      string routlongName = groupKey + "_Routlong2_vs_phi";
      gRoutlong2->SetName(routlongName.c_str());
      gRoutlong2->SetTitle((routlongName + ";#phi_{pair}-#Psi_{EP} (rad);R_{outlong}^{2} (fm^{2})").c_str());
      gRoutlong2->SetMarkerStyle(29);
      gRoutlong2->SetMarkerColor(kViolet + 1);
      gRoutlong2->SetLineColor(kViolet + 1);
      gRoutlong2->Fit(fitCosRoutlong2, "QN");
      gRoutlong2->Write();
      fitCosRoutlong2->Write();

      string rsidelongName = groupKey + "_Rsidelong2_vs_phi";
      gRsidelong2->SetName(rsidelongName.c_str());
      gRsidelong2->SetTitle((rsidelongName + ";#phi_{pair}-#Psi_{EP} (rad);R_{sidelong}^{2} (fm^{2})").c_str());
      gRsidelong2->SetMarkerStyle(47);
      gRsidelong2->SetMarkerColor(kAzure + 2);
      gRsidelong2->SetLineColor(kAzure + 2);
      gRsidelong2->Fit(fitSinRsidelong2, "QN");
      gRsidelong2->Write();
      fitSinRsidelong2->Write();
    }

    string alphaName = groupKey + "_alpha_vs_phi";
    gAlpha->SetName(alphaName.c_str());
    gAlpha->SetTitle((alphaName + ";#phi_{pair}-#Psi_{EP} (rad);#alpha").c_str());
    gAlpha->SetMarkerStyle(23);
    gAlpha->SetMarkerColor(kMagenta + 1);
    gAlpha->SetLineColor(kMagenta + 1);
    gAlpha->Fit(fitConstAlpha, "QN");
    gAlpha->Write();
    fitConstAlpha->Write();

    if (usesCoreHaloLambda) {
      string lambdaName = groupKey + "_lambda_vs_phi";
      gLambda->SetName(lambdaName.c_str());
      gLambda->SetTitle((lambdaName + ";#phi_{pair}-#Psi_{EP} (rad);#lambda").c_str());
      gLambda->SetMarkerStyle(33);
      gLambda->SetMarkerColor(kOrange + 7);
      gLambda->SetLineColor(kOrange + 7);
      gLambda->Fit(fitConstLambda, "QN");
      gLambda->Write();
      fitConstLambda->Write();
    }

    if (usesQ2Baseline) {
      string baselineName = groupKey + "_baselineQ2_vs_phi";
      gBaselineQ2->SetName(baselineName.c_str());
      gBaselineQ2->SetTitle((baselineName
                             + ";#phi_{pair}-#Psi_{EP} (rad);b_{q^{2}} "
                               "((GeV/c)^{-2})")
                                .c_str());
      gBaselineQ2->SetMarkerStyle(27);
      gBaselineQ2->SetMarkerColor(kTeal + 2);
      gBaselineQ2->SetLineColor(kTeal + 2);
      gBaselineQ2->Fit(fitConstBaselineQ2, "QN");
      gBaselineQ2->Write();
      fitConstBaselineQ2->Write();
    }

    delete gRout2;
    delete gRside2;
    delete gRlong2;
    delete gAlpha;
    delete fitCosRout2;
    delete fitCosRside2;
    delete fitCosRlong2;
    delete fitConstAlpha;
    if (usesCoreHaloLambda) {
      delete gLambda;
      delete fitConstLambda;
    }
    if (usesQ2Baseline) {
      delete gBaselineQ2;
      delete fitConstBaselineQ2;
    }

    if (hasOffDiagonal) {
      delete gRoutside2;
      delete gRoutlong2;
      delete gRsidelong2;
      delete fitSinRoutside2;
      delete fitCosRoutlong2;
      delete fitSinRsidelong2;
    }
  }

  wf->cd();
}

//==============================================================================
// Backward-Compatible Fit Wrappers
//==============================================================================

// Legacy wrappers kept for backward compatibility: fit all *_CF3D histograms in
// the input ROOT file without additional cent/mT filtering.
void FitCF3DWithLevy(const string &rpath,
                     const string &rfilename,
                     const string &wpath,
                     const string &wfilename,
                     const string &txtFilename,
                     const Levy3DFitOptions &fitOptions = Levy3DFitOptions(),
                     bool reopenOutputFilePerSlice = true) {
  const bool oldAddDirectory = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  auto rf = GetROOT(rpath, rfilename, "read");
  if (!InitializeROOTFile(wpath, wfilename)) {
    if (rf) {
      rf->Close();
      delete rf;
    }
    TH1::AddDirectory(oldAddDirectory);
    return;
  }

  if (!rf || rf->IsZombie()) {
    cout << "ERROR: cannot open input ROOT file " << rpath << "/" << rfilename << ".root" << endl;
    delete rf;
    TH1::AddDirectory(oldAddDirectory);
    return;
  }

  TFile *sharedOutputFile = nullptr;
  if (!reopenOutputFilePerSlice) {
    sharedOutputFile = OpenOutputROOTFileForUpdate(wpath, wfilename, "reuse output ROOT file during fitting");
    if (!sharedOutputFile) {
      rf->Close();
      delete rf;
      TH1::AddDirectory(oldAddDirectory);
      return;
    }
  }

  vector<Levy3DFitResult> fitResults;
  TIter nextKey(rf->GetListOfKeys());
  TKey *key = nullptr;

  while ((key = (TKey *)nextKey())) {
    string objectName = key->GetName();
    if (!EndsWith(objectName, "_CF3D")) {
      continue;
    }

    auto hCF = LoadStoredCFHistogram(rf, objectName);
    if (!hCF) {
      continue;
    }

    Levy3DFitResult fitResult;
    if (!ParseCF3DHistogramName(objectName, fitResult)) {
      cout << "WARNING: cannot parse histogram name: " << objectName << endl;
      delete hCF;
      continue;
    }

    TH3D *hSERaw = nullptr;
    TH3D *hMERaw = nullptr;
    if (fitOptions.usePML) {
      hSERaw = LoadStoredHistogramWithSuffix(rf, objectName, "_SE_raw3d");
      hMERaw = LoadStoredHistogramWithSuffix(rf, objectName, "_ME_raw3d");
      if (!hSERaw || !hMERaw) {
        cout << "WARNING: PML fit for " << objectName << " requires stored raw SE/ME histograms (_SE_raw3d/_ME_raw3d). "
             << "Please rebuild the CF ROOT file with the current CFCalc3D." << endl;
        delete hSERaw;
        delete hMERaw;
        delete hCF;
        continue;
      }
    }

    cout << "Fitting " << objectName << endl;
    if (FitAndWriteSingleCFHistogram(
            hCF, hSERaw, hMERaw, fitResult, wpath, wfilename, false, fitOptions, sharedOutputFile)) {
      fitResults.push_back(fitResult);
    }
    delete hSERaw;
    delete hMERaw;
    delete hCF;
  }

  TFile *wf = sharedOutputFile;
  const bool ownsSummaryFile = (wf == nullptr);
  if (ownsSummaryFile) {
    wf = OpenOutputROOTFileForUpdate(wpath, wfilename, "update output ROOT file for summary graphs");
  }
  if (!wf || wf->IsZombie()) {
    cout << "ERROR: cannot update output ROOT file " << wpath << "/" << wfilename << ".root for summary graphs" << endl;
  } else {
    WriteR2Graphs(wf, fitResults);
    if (ownsSummaryFile) {
      wf->Close();
      delete wf;
    } else {
      wf->cd();
    }
  }
  if (sharedOutputFile) {
    sharedOutputFile->Close();
    delete sharedOutputFile;
  }
  WriteFitResultsSummary(wpath + "/" + txtFilename + ".txt", fitResults);

  rf->Close();
  delete rf;
  cout << "Diagonal 3D Levy fit results have been written to " << wpath << "/" << wfilename << ".root and " << wpath
       << "/" << txtFilename << ".txt" << endl;
  TH1::AddDirectory(oldAddDirectory);
}

void FitCF3DWithFullLevy(const string &rpath,
                         const string &rfilename,
                         const string &wpath,
                         const string &wfilename,
                         const string &txtFilename,
                         const Levy3DFitOptions &fitOptions = Levy3DFitOptions(),
                         bool reopenOutputFilePerSlice = true) {
  const bool oldAddDirectory = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  auto rf = GetROOT(rpath, rfilename, "read");
  if (!InitializeROOTFile(wpath, wfilename)) {
    if (rf) {
      rf->Close();
      delete rf;
    }
    TH1::AddDirectory(oldAddDirectory);
    return;
  }

  if (!rf || rf->IsZombie()) {
    cout << "ERROR: cannot open input ROOT file " << rpath << "/" << rfilename << ".root" << endl;
    delete rf;
    TH1::AddDirectory(oldAddDirectory);
    return;
  }

  TFile *sharedOutputFile = nullptr;
  if (!reopenOutputFilePerSlice) {
    sharedOutputFile = OpenOutputROOTFileForUpdate(wpath, wfilename, "reuse output ROOT file during fitting");
    if (!sharedOutputFile) {
      rf->Close();
      delete rf;
      TH1::AddDirectory(oldAddDirectory);
      return;
    }
  }

  vector<Levy3DFitResult> fitResults;
  TIter nextKey(rf->GetListOfKeys());
  TKey *key = nullptr;

  while ((key = (TKey *)nextKey())) {
    string objectName = key->GetName();
    if (!EndsWith(objectName, "_CF3D")) {
      continue;
    }

    auto hCF = LoadStoredCFHistogram(rf, objectName);
    if (!hCF) {
      continue;
    }

    Levy3DFitResult fitResult;
    if (!ParseCF3DHistogramName(objectName, fitResult)) {
      cout << "WARNING: cannot parse histogram name: " << objectName << endl;
      delete hCF;
      continue;
    }

    TH3D *hSERaw = nullptr;
    TH3D *hMERaw = nullptr;
    if (fitOptions.usePML) {
      hSERaw = LoadStoredHistogramWithSuffix(rf, objectName, "_SE_raw3d");
      hMERaw = LoadStoredHistogramWithSuffix(rf, objectName, "_ME_raw3d");
      if (!hSERaw || !hMERaw) {
        cout << "WARNING: PML fit for " << objectName << " requires stored raw SE/ME histograms (_SE_raw3d/_ME_raw3d). "
             << "Please rebuild the CF ROOT file with the current CFCalc3D." << endl;
        delete hSERaw;
        delete hMERaw;
        delete hCF;
        continue;
      }
    }

    cout << "Fitting full Levy model for " << objectName << endl;
    if (FitAndWriteSingleCFHistogram(
            hCF, hSERaw, hMERaw, fitResult, wpath, wfilename, true, fitOptions, sharedOutputFile)) {
      fitResults.push_back(fitResult);
    }
    delete hSERaw;
    delete hMERaw;
    delete hCF;
  }

  TFile *wf = sharedOutputFile;
  const bool ownsSummaryFile = (wf == nullptr);
  if (ownsSummaryFile) {
    wf = OpenOutputROOTFileForUpdate(wpath, wfilename, "update output ROOT file for summary graphs");
  }
  if (!wf || wf->IsZombie()) {
    cout << "ERROR: cannot update output ROOT file " << wpath << "/" << wfilename << ".root for summary graphs" << endl;
  } else {
    WriteR2Graphs(wf, fitResults);
    if (ownsSummaryFile) {
      wf->Close();
      delete wf;
    } else {
      wf->cd();
    }
  }
  if (sharedOutputFile) {
    sharedOutputFile->Close();
    delete sharedOutputFile;
  }
  WriteFitResultsSummary(wpath + "/" + txtFilename + ".txt", fitResults);

  rf->Close();
  delete rf;
  cout << "Full 3D Levy fit results have been written to " << wpath << "/" << wfilename << ".root and " << wpath << "/"
       << txtFilename << ".txt" << endl;
  TH1::AddDirectory(oldAddDirectory);
}

//==============================================================================
// Macro Entry Point
//==============================================================================

// Main macro configuration.
//
// Workflow convention:
//   1. doBuildCF3D=true  -> read AnalysisResults and create
//   EP_dependence_CF_*.root
//   2. doFitDiag/full    -> read the CF ROOT file and perform Levy fits
//
// The fit stage intentionally reads the CF ROOT file produced in step 1, so the
// expensive same-event / mixed-event projections do not need to be recomputed
// every time the fit settings are changed.
void _3d_cf_from_exp() {
  // string suffix = "_its";
  // string suffix = "_2400bins";
  string suffix = "_no_epmix";
  // string suffix = "_cent_mix";
  // string data_set = "23zzi_pass5";
  // string data_set = "23zzh_pass5_noepmix";
  // string data_set = "23zzm_pass5_medium";
  // string data_set = "23zzh_pass5_t3_morebin";
  // string data_set = "23zzh_pass5_t4_morebin";
  //   string data_set = "23zzh_pass5_t4_fulltof";
  // string data_set = "23zzh_pass5_t3_morebin";
  // string data_set = "23zzh_pass5_small_t3";
  // string data_set = "23zzh_pass5_small_dbgt4";
  // string data_set = "23zzh_pass5_small_noepmix";
  string data_set = "23zzh_pass5_medium";

  // string data_set = "25ae_pass2";
  // string data_set = "25ae_pass2_test";
  // string data_set = "25ae_pass2_noepmix";
  // string data_set = "25ae_pass2_noepmix_moremult";
  // string data_set = "25ae_pass2_epmix";
  // string data_set = "25ae_pass2_moredepth";

  //--------------------------------------------------------------------------
  // Dataset / path selection
  //--------------------------------------------------------------------------
  string rpath_pbpb = "/Users/allenzhou/ALICE/alidata/hyperloop_res/femtoep/PbPb";
  string rpath_oo = "/Users/allenzhou/ALICE/alidata/hyperloop_res/femtoep/OO";
  string wpath_pbpb = "/Users/allenzhou/ALICE/alidata/femtoep_res/PbPb";
  string wpath_oo = "/Users/allenzhou/ALICE/alidata/femtoep_res/OO";

  string wpath, rpath;

  // enum system = {"PbPb", "OO"}; forbidden
  //   string rpath = rpath_pbpb;
  bool kisoo = 0;
  bool kispbpb = 1;

  if (kisoo && kispbpb) {
    cout << "ERROR: both kisoo and kispbpb cannot be true at the same time." << endl;
    return;
  }
  if (kisoo) {
    rpath = rpath_oo;
    wpath = wpath_oo;
  } else if (kispbpb) {
    rpath = rpath_pbpb;
    wpath = wpath_pbpb;
  }

  string rfilename = "AnalysisResults_3d_" + data_set;
  string rtaskname_main = "femto-dream-pair-task-track-track";
  string rtaskname_sub = "femto-dream-pair-task-track-track" + suffix;
  string rsubtask_se = "SameEvent_3Dqn";
  string rsubtask_me = "MixedEvent_3Dqn";

  string wfilename_main = "EP_dependence_CF_" + data_set;
  string wfilename_sub = "EP_dependence_CF_" + data_set + suffix;

  string wfilename, rtaskname;

  bool kis_subwagon = 0;

  if (kis_subwagon) {
    rtaskname = rtaskname_sub;
    wfilename = wfilename_sub;
  } else {
    rtaskname = rtaskname_main;
    wfilename = wfilename_main;
  }

  //--------------------------------------------------------------------------
  // CF production binning
  //--------------------------------------------------------------------------
  std::vector<std::pair<double, double>> centBins = {{0, 10}, {10, 30}, {30, 50}, {50, 80}, {80, 100}};
  std::vector<std::pair<double, double>> mTBins_kt = {
      {0.244131, 0.331059}, {0.331059, 0.423792}, {0.423792, 0.51923}, {0.51923, 0.713863}, {0.244131, 0.713863}};
  std::vector<std::pair<double, double>> mTBins_mt = {{0.2, 0.3},
                                                      {0.3, 0.4},
                                                      {0.4, 0.5},
                                                      {0.5, 0.6},
                                                      {0.6, 0.7},
                                                      {0.7, 0.8},
                                                      {0.8, 1.0},
                                                      {1.0, 1.2},
                                                      {1.2, 1.6},
                                                      {1.6, 2}};

  //--------------------------------------------------------------------------
  // Stage switches
  //--------------------------------------------------------------------------
  bool doBuildCF3D = false;
  bool doFitDiag = false;
  bool doFitFull = true;
  bool mapPairPhiToSymmetricRange = false;
  bool writeNormalizedSEME1DProjections = false;
  bool reopenOutputFilePerSlice = true;
  bool fitUseCoulomb = true;
  bool fitUseCoreHaloLambda = true;
  bool fitUseQ2Baseline = true;
  bool fitUsePML = true;

  //--------------------------------------------------------------------------
  // Fit binning (can be a subset of the CF production binning)
  //--------------------------------------------------------------------------
  std::vector<std::pair<double, double>> fitCentBins = {{0, 10}, {10, 30}, {30, 50}};
  std::vector<std::pair<double, double>> fitMTBins_mt = {
      {0.20, 0.30},
      {0.30, 0.40},
      {0.40, 0.50},
  };
  std::vector<std::pair<double, double>> fitMTBins_kt = {{0.24, 0.33}, {0.42, 0.52}, {0.52, 0.71}, {0.24, 0.71}};

  std::vector<std::pair<double, double>> fitMTBins = fitMTBins_kt;
  std::vector<std::pair<double, double>> mTBins = mTBins_kt;

  //--------------------------------------------------------------------------
  // Fit model options and output naming
  //--------------------------------------------------------------------------
  Levy3DFitOptions fitOptions;
  fitOptions.useCoulomb = fitUseCoulomb;
  fitOptions.useCoreHaloLambda = fitUseCoreHaloLambda;
  fitOptions.useQ2Baseline = fitUseQ2Baseline;
  fitOptions.usePML = fitUsePML;
  fitOptions.fitQMax = 0.15;
  const string fitOptionTag = BuildFitOptionTag(fitOptions);

  string fitInputPath = wpath;
  string fitInputFilename = wfilename;
  string fitDiagOutput = "EP_dependence_CF_fit_" + fitOptionTag + "_" + data_set + (kis_subwagon ? suffix : "");
  string fitDiagTxt = "EP_dependence_CF_fit_params_" + fitOptionTag + "_" + data_set + (kis_subwagon ? suffix : "");
  string fitFullOutput = "EP_dependence_CF_full_fit_" + fitOptionTag + "_" + data_set + (kis_subwagon ? suffix : "");
  string fitFullTxt =
      "EP_dependence_CF_full_fit_params_" + fitOptionTag + "_" + data_set + (kis_subwagon ? suffix : "");

  if (doBuildCF3D) {
    CFCalc3D(rpath,
             rfilename,
             rtaskname,
             rsubtask_se,
             rsubtask_me,
             wpath,
             wfilename,
             centBins,
             mTBins,
             mapPairPhiToSymmetricRange,
             writeNormalizedSEME1DProjections,
             reopenOutputFilePerSlice);
  }

  if (doFitDiag) {
    FitCF3DWithSelectedBins(fitInputPath,
                            fitInputFilename,
                            wpath,
                            fitDiagOutput,
                            fitDiagTxt,
                            fitCentBins,
                            fitMTBins,
                            false,
                            fitOptions,
                            reopenOutputFilePerSlice);
  }

  if (doFitFull) {
    FitCF3DWithSelectedBins(fitInputPath,
                            fitInputFilename,
                            wpath,
                            fitFullOutput,
                            fitFullTxt,
                            fitCentBins,
                            fitMTBins,
                            true,
                            fitOptions,
                            reopenOutputFilePerSlice);
  }
}
