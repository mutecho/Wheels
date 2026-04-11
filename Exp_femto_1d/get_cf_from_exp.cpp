#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include "TF1.h"
#include "TH1.h"
#include "TTree.h"
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <THnSparse.h>
#include <TLegend.h>
#include <TMath.h>

using namespace std;

enum EventType { kSameEvent, kMixedEvent };

struct EPBin {
  double low1, high1;
  double low2, high2;
  string name;
};

bool FileExists_bool(const string &filename) {
  ifstream file(filename);
  return file.good();
}

void FileExists_warn(const string &filename) {
  if (!FileExists_bool(filename)) {
    cerr << "Error: File: " << filename << " doesn't exist!" << endl;
  }
}

unique_ptr<TFile> GetROOT_unique(const string Readpath,
                                 const string ReadFilename,
                                 const string Option) {
  string lowerOption = Option;
  transform(lowerOption.begin(), lowerOption.end(), lowerOption.begin(),
            [](unsigned char c) { return tolower(c); }); // 全转小写
  if (!(lowerOption == "read" || lowerOption == "write" ||
        lowerOption == "recreate")) {
    cout << "Invalid Option: " << Option
         << ", should be one of: read, write, recreate" << endl;
  }
  string rfilename = Readpath + "/" + ReadFilename + ".root"; // 加上root后缀
  FileExists_warn(rfilename);

  return make_unique<TFile>(rfilename.c_str(), Option.c_str());
}

TFile *GetROOT(const string Readpath, const string ReadFilename,
               const string Option) {
  string lowerOption = Option;
  transform(lowerOption.begin(), lowerOption.end(), lowerOption.begin(),
            [](unsigned char c) { return tolower(c); }); // 全转小写
  if (!(lowerOption == "read" || lowerOption == "write" ||
        lowerOption == "recreate")) {
    cout << "Invalid Option: " << Option
         << ", should be one of: read, write, recreate" << endl;
  }
  string rfilename = Readpath + "/" + ReadFilename + ".root"; // 加上root后缀
  FileExists_warn(rfilename);

  return new TFile(rfilename.c_str(), Option.c_str());
}

std::string doubleToString(double x, int precision = 2) {
  std::stringstream ss;
  ss << std::fixed << std::setprecision(precision) << x;
  return ss.str();
}

TH1D *GetCFfromSM(TH1D *h_se, TH1D *h_me, double normLow, double normHigh,
                  double kstarMax) {
  TF1 *constant = new TF1("constant", "1", 0, 10);
  cout << "Normalizing SE and ME histograms between k* = " << normLow << " and "
       << normHigh << " GeV/c." << endl;
  TH1D *h_se_c = (TH1D *)h_se->Clone();
  TH1D *h_me_c = (TH1D *)h_me->Clone();
  cout << "Cloned SE and ME histograms." << endl;

  int binNorm[2];
  binNorm[0] = h_se_c->FindBin(normLow);
  binNorm[1] = h_se_c->FindBin(normHigh);
  double factorN = h_me_c->Integral(binNorm[0], binNorm[1]) /
                   h_se_c->Integral(binNorm[0], binNorm[1]);
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

TH1D *onefromndHisto(unique_ptr<TFile> rfile, string hpath, bool isReranged,
                     double rerangeMax) {
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

  h1c_re->Reset(); // 清空内容

  for (int i = binRange[0]; i <= binRange[1]; i++) {
    double content = h1c->GetBinContent(i);
    double error = h1c->GetBinError(i);
    int newBin = i - binRange[0] + 1; // 从1开始填
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

TH1D *ndHistoRead(TFile *rf, string histoname, int centLow, int centHigh,
                  double mTLow, double mTHigh, EPBin epbin) {
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
  string cname = "cent=" + to_string(static_cast<int>(centLow)) + "-" +
                 to_string(static_cast<int>(centHigh)) +
                 "_mT=" + doubleToString(mTLow, 1) + "-" +
                 doubleToString(mTHigh, 1);
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

  if (epbin.low2 > 0) { // 有镜像部分
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

TH1D *calc_cf_from_sme_rerange(TH1D *h_se, TH1D *h_me,
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
  double factorN = h_me_c->Integral(binNorm[0], binNorm[1]) /
                   h_se_c->Integral(binNorm[0], binNorm[1]);
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
    int newBin = i - binRange[0] + 1; // 从1开始填
    h_cf_re->SetBinContent(newBin, content);
    h_cf_re->SetBinError(newBin, error);
  }
  return h_cf_re;
}

void CFCalcWith_Cent_Mt_pairphi_full(
    string rpath, string rfilename, string taskname, string subtaskname_se,
    string subtaskname_me, string wpath, string wfilename,
    std::vector<std::pair<double, double>> centBins,
    std::vector<std::pair<double, double>> mTBins, std::vector<EPBin> epBins_se,
    std::vector<EPBin> epBins_me, std::vector<std::string> pairphiBins,
    std::pair<double, double> normrange = {0.5, 0.8},
    std::pair<double, double> kstarRange = {0.0, 0.5}) {
  auto rf = GetROOT(rpath, rfilename, "read");
  auto wfFemtoep = GetROOT(wpath, wfilename, "recreate");

  // 重归一化的区间
  int binNorm[2];
  // 用于计算的fake histo
  double factorN;
  // 最后要显示的区间
  int binRange[2];
  //   binRange[0] = h1sc->FindBin(0.);

  TF1 *constant = new TF1("constant", "1", 0, 10);

  string hpath_se =
      taskname + "/" + subtaskname_se + "/relPairkstarmTMultMultPercentileQn";
  string hpath_me =
      taskname + "/" + subtaskname_me + "/relPairkstarmTMultMultPercentileQn";

  //   TH1D *h1_se = onefromndHisto(rf, hpath_se, false, 0);
  //   TH1D *h1_me = onefromndHisto(rf, hpath_me, false, 0);

  cout << "Starting CF calculation with cent, mT, pairphi bins..." << endl;

  for (auto [centLow, centHigh] : centBins) {
    for (auto [mTLow, mTHigh] : mTBins) {
      string hname = "cent=" + to_string(static_cast<int>(centLow)) + "-" +
                     to_string(static_cast<int>(centHigh)) +
                     "_mT=" + doubleToString(mTLow, 1) + "-" +
                     doubleToString(mTHigh, 1);
      // string hname_me = hname + "/" + "ME_Minbias_EP";
      string hname_me = hname + "_" + "ME_Minbias_EP";
      TH1D *h1_me = ndHistoRead(rf, hpath_me, centLow, centHigh, mTLow, mTHigh,
                                epBins_me[0]);
      if (!h1_me) {
        cout << "ERROR: ME histo not found: " << hname_me << endl;
        continue;
      }
      // h1_me->Sumw2();
      wfFemtoep->WriteObject(h1_me, hname_me.c_str());
      TH1D *h1mc = (TH1D *)(h1_me->Clone());

      string hname_se = hname + "_SE_Minbias_EP";
      // string hname_se = hname + "/SE_Minbias_EP";
      // auto h1_se = (TH1D *)rf_se->Get(hname_se.c_str());
      TH1D *h1_se = ndHistoRead(rf, hpath_se, centLow, centHigh, mTLow, mTHigh,
                                epBins_me[0]);
      if (!h1_se) {
        cout << "ERROR: ME histo not found: " << hname_me << endl;
        continue;
      }
      // h1_se->Sumw2();
      cout << "Processing: " << hname_se << endl;
      TH1D *h1sc = (TH1D *)h1_se->Clone();
      wfFemtoep->WriteObject(h1_se, hname_se.c_str());
      // for (auto phibin_me : epBins_me) {
      //     string hname_me = hname + "/" + "Min bias EP";
      //     TH1D *h1_me = rf_me->Get<TH1D>(hname_me.c_str());
      //     wfFemtoep->WriteObject(h1_me, hname_me.c_str());
      //     if (!h1_me) {
      //         cout << "ERROR: ME histo not found: " << hname_me << endl;
      //         continue;
      //     }
      //     TH1D *h1mc = (TH1D *)(h1_me->Clone());
      // }
      TH1D *h1cf_re_mb =
          calc_cf_from_sme_rerange(h1_se, h1_me, normrange, kstarRange);
      string wHistoname = hname + "_CF_reranged_Minbias_EP";
      // string wHistoname = hname + "/CF_reranged_Minbias_EP";
      h1cf_re_mb->SetTitle("reranged CF from ndHisto");
      h1cf_re_mb->GetXaxis()->SetTitle("k* (GeV/c)");
      h1cf_re_mb->GetYaxis()->SetTitle("C(k*)");
      wfFemtoep->WriteObject(h1cf_re_mb, wHistoname.c_str());

      TH1D *h_cf_cache[2];
      int icache = 0;

      for (auto phibin_se : epBins_se) {
        TH1D *h1mc = (TH1D *)(h1_me->Clone());
        // string hname_se_ep = hname + "/SE_" + phibin_se.name;
        string hname_se_ep = hname + "_SE_" + phibin_se.name;
        // auto h1_se = (TH1D *)rf_se->Get(hname_se.c_str());
        TH1D *h1_se = ndHistoRead(rf, hpath_se, centLow, centHigh, mTLow,
                                  mTHigh, phibin_se);
        if (!h1_se) {
          cout << "ERROR: ME histo not found: " << hname_me << endl;
          continue;
        }
        // h1_se->Sumw2();
        cout << "Processing: " << hname_se_ep << endl;
        TH1D *h1sc = (TH1D *)h1_se->Clone();
        wfFemtoep->WriteObject(h1_se, hname_se_ep.c_str());
        // TH1D *h1cf_re = GetCFfromSM(h1_se, h1_me, 0.5, 0.8, 0.25);

        // // binNorm[0] = h1sc->FindBin(0.5);
        // // binNorm[1] = h1sc->FindBin(0.8);
        // binNorm[0] = h1sc->FindBin(normrange.first);
        // binNorm[1] = h1sc->FindBin(normrange.second);
        // factorN = h1mc->Integral(binNorm[0], binNorm[1]) /
        //           h1sc->Integral(binNorm[0], binNorm[1]);
        // h1mc->Divide(constant, factorN);
        // h1sc->Divide(h1mc);
        // // string wHistoname = hname + "/CF_reranged_" + phibin_se.name;
        // // h1sc->SetTitle("reranged CF from ndHisto");
        // // h1sc->GetXaxis()->SetTitle("k* (GeV/c)");
        // // h1sc->GetXaxis()->SetRange(0, 0.5);
        // // h1sc->GetYaxis()->SetTitle("C(k*)");
        // // wfFemtoep->WriteObject(h1sc, wHistoname.c_str());

        // binRange[0] = h1sc->FindBin(kstarRange.first);
        // binRange[1] = h1sc->FindBin(kstarRange.second);

        // // cout << "binRange: " << binRange[0] << " to " << binRange[1] <<
        // endl;

        // int nBins = binRange[1] - binRange[0] + 1;
        // // cout << "nBins: " << nBins << endl;
        // TH1D *h1cf_re =
        //     new TH1D("", "", nBins, kstarRange.first, kstarRange.second);
        // // TH1D *h1cf_re = new TH1D();

        // for (int i = binRange[0]; i <= binRange[1]; i++) {
        //   double content = h1sc->GetBinContent(i);
        //   double error = h1sc->GetBinError(i);
        //   int newBin = i - binRange[0] + 1; // 从1开始填
        //   // cout << "Bin " << newBin << " content: " << content
        //   //      << ", error: " << error << endl;
        //   h1cf_re->SetBinContent(newBin, content);
        //   h1cf_re->SetBinError(newBin, error);
        // }
        TH1D *h1cf_re =
            calc_cf_from_sme_rerange(h1_se, h1_me, normrange, kstarRange);
        // string wHistoname = hname + "/CF_reranged_" + phibin_se.name;
        string wHistoname = hname + "_CF_reranged_" + phibin_se.name;
        h1cf_re->SetTitle("reranged CF from ndHisto");
        h1cf_re->GetXaxis()->SetTitle("k* (GeV/c)");
        h1cf_re->GetYaxis()->SetTitle("C(k*)");
        wfFemtoep->WriteObject(h1cf_re, wHistoname.c_str());
        h_cf_cache[icache++] = h1cf_re;
        // wfFemtoep->WriteObject(h1sc, wHistoname.c_str());
      }
      h_cf_cache[0]->SetLineColor(kRed);
      h_cf_cache[1]->SetLineColor(kBlue);

      TCanvas *c = new TCanvas("c", "in and out", 800, 600);

      /* 先画第一张 */
      h_cf_cache[0]->Draw("E");

      /* 再叠加第二张 */
      h_cf_cache[1]->Draw("E SAME");

      /* 图例 */
      TLegend *leg = new TLegend(0.65, 0.7, 0.88, 0.88);
      leg->AddEntry(h_cf_cache[0], "in", "l");
      leg->AddEntry(h_cf_cache[1], "out", "l");
      leg->Draw();
      string canvas_name = hname + "_CF_reranged_in_and_out_plane";
      wfFemtoep->WriteObject(c, canvas_name.c_str());
      c->Delete();
    }
  }
}

void get_cf_from_exp() {
  // string suffix = "_its";
  // string suffix = "_2400bins";
  string suffix = "_no_epmix";
  // string suffix = "_cent_mix";
  // string data_set = "23zzi_pass5";
  // string data_set = "23zzh_pass5_noepmix";
  // string data_set = "23zzm_pass5_medium";
  // string data_set = "23zzh_pass5_t3_morebin";
  // string data_set = "23zzh_pass5_t4_morebin";
  string data_set = "23zzh_pass5_t4_fulltof";
  // string data_set = "23zzh_pass5_t3_morebin";
  // string data_set = "23zzh_pass5_small_t3";
  // string data_set = "23zzh_pass5_small_dbgt4";
  // string data_set = "23zzh_pass5_small_noepmix";

  // string data_set = "25ae_pass2";
  // string data_set = "25ae_pass2_noepmix";
  // string data_set = "25ae_pass2_noepmix_moremult";
  // string data_set = "25ae_pass2_epmix";
  // string data_set = "25ae_pass2_moredepth";

  // string rpath1 = "/Users/allenzhou/ALICE/alidata/femtoep_res";
  string rpath_pbpb =
      "/Users/allenzhou/ALICE/alidata/hyperloop_res/femtoep/PbPb";
  string rpath_oo = "/Users/allenzhou/ALICE/alidata/hyperloop_res/femtoep/OO";

  // enum system = {"PbPb", "OO"}; forbidden
  string rpath = rpath_pbpb;
  // string rpath = rpath_oo;

  // string rpath2 = "/Users/allenzhou/ALICE/scripts/femtoep";
  // string rpath2 = rpath1;
  // string rpath2 = "/Users/allenzhou/ALICE/alidata/femtoep_res";
  // string rfilename = "23_zzh_pass5_small_3";
  // string rfilename = "AnalysisResults_23zzf_pass5";
  // string rfilename = "AnalysisResults_" + data_set;
  string rfilename = "AnalysisResults_" + data_set;
  // string rtaskname = "femto-dream-pair-task-track-track";
  string rtaskname_main = "femto-dream-pair-task-track-track";
  string rtaskname_sub = "femto-dream-pair-task-track-track" + suffix;
  string rsubtask_se = "SameEvent_EP";
  string rsubtask_me = "MixedEvent_EP";

  // string wpath = "/Users/allenzhou/ALICE/scripts/femtoep";
  string wpath = "/Users/allenzhou/ALICE/alidata/femtoep_res";
  string wfilename_main = "EP_dependence_CF_" + data_set;
  string wfilename_sub = "EP_dependence_CF_" + data_set + suffix;
  // string wfilename2 = "Check EP mixing 2";

  string wfilename, rtaskname;

  bool kis_subwagon = 1;

  if (kis_subwagon) {
    rtaskname = rtaskname_sub;
    wfilename = wfilename_sub;
  } else {
    rtaskname = rtaskname_main;
    wfilename = wfilename_main;
  }

  // std::vector<std::pair<double, double>> centBins = {
  //     {0, 30}, {30, 70}, {70, 100}};
  // std::vector<std::pair<double, double>> mTBins = {
  //     {0.5, 0.7}, {0.7, 1.0}, {1.0, 1.5}};
  std::vector<std::pair<double, double>> centBins = {
      {0, 10}, {10, 30}, {30, 50}, {50, 80}, {80, 100}};
  std::vector<std::pair<double, double>> mTBins = {{0.244131, 0.331059},
                                                   {0.331059, 0.423792},
                                                   {0.423792, 0.51923},
                                                   {0.51923, 0.713863},
                                                   {0.244131, 0.713863}};
  std::vector<std::string> pairphiBins = {"In_plane", "Out_of_plane"};

  std::vector<EPBin> epBins_se = {
      {0, TMath::Pi() / 4, 3 * TMath::Pi() / 4, TMath::Pi(), "In_plane"},
      {TMath::Pi() / 4, 3 * TMath::Pi() / 4, -1, -1, "Out_of_plane"}};
  std::vector<EPBin> epBins_me = {{{0, TMath::Pi(), -1, -1, "Min bias EP"}}};

  CFCalcWith_Cent_Mt_pairphi_full(rpath, rfilename, rtaskname, rsubtask_se,
                                  rsubtask_me, wpath, wfilename, centBins,
                                  mTBins, epBins_se, epBins_me, pairphiBins,
                                  {0.5, 0.8}, {0., 0.8});
}