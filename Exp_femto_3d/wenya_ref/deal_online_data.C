enum obsType { kOut = 0, kSide, kLong, kMt, kCent, kQn, kPhi };
enum evtType { kSAME = 0, kMIXED };
enum color { kRe = 0, kBl, kGr, kPur, kYl, kPk, kGry, kBlk };

double centBin[] = {30, 50};
double mTBin[] = {0.2, 0.3};
// double centBin[] = {10,30};
// double mTBin[]   = {0.2,0.3};

double qnBin[] = {0, 100};

const int nCentBins = sizeof(centBin) / sizeof(double) - 1;
const int nMtBins = sizeof(mTBin) / sizeof(double) - 1;
const int nPhiBins = 1;
const int nQnBins = 1;

TFile *fileInSame;
TFile *fileInMixed;

// SAME and MIXED histograms
vector<vector<vector<TH3F *>>> histoNum;
vector<vector<TH3F *>> histoDen;

// qi slices → 1D CF
vector<vector<vector<TH1F *>>> histoCF_SM;
vector<vector<vector<TF1 *>>> fnormbkg;
vector<vector<TH1F *>> histoCF_bkg;

// φ–resolved
vector<vector<vector<vector<TH3F *>>>> histoNum_phi;
vector<vector<vector<TH3F *>>> histoDen_phi;
vector<vector<vector<vector<TH3F *>>>> histoCF_SM_phi;

// function declarations
void readFile(TString filename, TFile *&file);
void readCF_generic(TFile *filein,
                    int evtType,
                    int phiBin,
                    bool usePhi,
                    vector<vector<vector<TH3F *>>> &numStore,
                    vector<vector<TH3F *>> &denStore);
void readCF(TFile *filein, int type);
void readCF_phi(TFile *filein, int type, int phiBin);

void projectCF1D(TH1F *&hCF1D,
                 TH1F *&hCF1DBkg,
                 TH3F *hNum3D,
                 TH3F *hDen3D,
                 TString target,
                 TF1 *&fNormBkg,
                 int qn,
                 int binsFromZero = 1);
TF1 *normBkg(TH1F *&hfCF1D, int iqn);
float norm(TH1F *&hfCF1D);
void drawAndCheck(TString drawOpt,
                  vector<int> color,
                  int canvaswidth,
                  int canvasHight,
                  bool isDivided,
                  bool isMultPlots,
                  int divideX,
                  int divideY,
                  TH1F *hDrawSingle,
                  vector<TH1F *> hDraw,
                  vector<TF1 *> fBkg = {});
void SetHistoSytle(TH1F *h, int color, int style);
void SetStyle(Bool_t graypalette = kFALSE);
int GetRootColor(int c);

int phibinCheck = 6;
// ===============================================================
// MAIN FUNCTION
// ===============================================================
void deal_online_data(TString fileNameSame, TString fileNameMixed, int analysisMode = 0) {
  switch (analysisMode) {
    case 0: {  // draw and test correlation functions

      // -----------------------------------------------------
      // read files
      // -----------------------------------------------------
      readFile(fileNameSame, fileInSame);
      readFile(fileNameMixed, fileInMixed);

      // resize storage
      histoNum.clear();
      histoNum.resize(nMtBins, vector<vector<TH3F *>>(nCentBins, vector<TH3F *>(nQnBins, nullptr)));

      histoDen.clear();
      histoDen.resize(nMtBins, vector<TH3F *>(nCentBins, nullptr));

      // read SAME as Num; MIXED as Den
      readCF(fileInSame, kSAME);
      readCF(fileInMixed, kMIXED);

      // prepare 1D hist storage
      histoCF_SM.clear();
      histoCF_SM.resize(nMtBins, vector<vector<TH1F *>>(nCentBins, vector<TH1F *>(nQnBins, nullptr)));

      fnormbkg.clear();
      fnormbkg.resize(nMtBins, vector<vector<TF1 *>>(nCentBins, vector<TF1 *>(nQnBins, nullptr)));

      histoCF_bkg.clear();
      histoCF_bkg.resize(nMtBins, vector<TH1F *>(nCentBins, nullptr));

      // -----------------------------------------------------
      // project to 1D (CORRECT METHOD!)
      // -----------------------------------------------------
      for (int iqn = 0; iqn < nQnBins; iqn++) {
        projectCF1D(histoCF_SM[0][0][iqn],
                    histoCF_bkg[0][0],
                    histoNum[0][0][iqn],
                    histoDen[0][0],
                    "out",
                    fnormbkg[0][0][iqn],
                    iqn);
      }

      // -----------------------------------------------------
      // draw
      // -----------------------------------------------------
      SetStyle();
      vector<TH1F *> outPlots;
      for (int iqn = 0; iqn < nQnBins; iqn++) outPlots.push_back(histoCF_SM[0][0][iqn]);

      vector<TF1 *> outfBkg;
      for (int iqn = 0; iqn < nQnBins; iqn++) outfBkg.push_back(fnormbkg[0][0][iqn]);

      drawAndCheck("pl", {kRe, kBl, kGr}, 900, 600, false, true, -1, -1, nullptr, outPlots);

      // // -----------------------------------------------------
      // // write
      // // -----------------------------------------------------
      // TFile* fout = new TFile(Form("output/CF/cf_out_3DCF_mT_cent_qn_phi%d_largeCPR.root",phibinCheck-1),"RECREATE");
      // fout->cd();
      // for (int iqn=0; iqn<nQnBins; iqn++){
      // 	outPlots[iqn]->Write(Form("cf_out_cent%d_mt%d_qn%d_phi%d",0,0,iqn,phibinCheck-1));
      // }

      break;
    }
    case 1: {
      // obtain correlation functions

      // -----------------------------------------------------
      // read files
      // -----------------------------------------------------
      readFile(fileNameSame, fileInSame);
      readFile(fileNameMixed, fileInMixed);

      // resize storage
      histoNum_phi.clear();
      histoNum_phi.resize(
          nPhiBins,
          vector<vector<vector<TH3F *>>>(nMtBins, vector<vector<TH3F *>>(nCentBins, vector<TH3F *>(nQnBins, nullptr))));

      histoDen_phi.clear();
      histoDen_phi.resize(nPhiBins, vector<vector<TH3F *>>(nMtBins, vector<TH3F *>(nCentBins, nullptr)));

      // read SAME as Num; MIXED as Den
      for (int iPhi(1); iPhi <= nPhiBins; ++iPhi) {
        readCF_phi(fileInSame, kSAME, iPhi);
        readCF_phi(fileInMixed, kMIXED, iPhi);
      }

      // -----------------------------------------------------
      // Same / Mixed
      // -----------------------------------------------------
      for (int iPhi = 0; iPhi < nPhiBins; iPhi++) {
        for (int iMt = 0; iMt < nMtBins; iMt++) {
          for (int iCent = 0; iCent < nCentBins; iCent++) {
            for (int iqn = 0; iqn < nQnBins; iqn++) {
              histoNum_phi[iPhi][iMt][iCent][iqn]->Scale(1. / histoNum_phi[iPhi][iMt][iCent][iqn]->Integral("width"));
            }
            histoDen_phi[iPhi][iMt][iCent]->Scale(1. / histoDen_phi[iPhi][iMt][iCent]->Integral("width"));
          }
        }
      }

      // -----------------------------------------------------
      // write
      // -----------------------------------------------------
      TFile *fout = new TFile("output/out_3DCF_mT_cent_qn_phi_noCPR.root", "RECREATE");
      fout->cd();
      for (int iCent = 0; iCent < nCentBins; iCent++) {
        for (int iMt = 0; iMt < nMtBins; iMt++) {
          for (int iPhi = 0; iPhi < nPhiBins; iPhi++) {
            histoDen_phi[iPhi][iMt][iCent]->Write();
            for (int iqn = 0; iqn < nQnBins; iqn++) {
              histoNum_phi[iPhi][iMt][iCent][iqn]->Write();
            }
          }
        }
      }

      break;
    }
    default: {
      break;
    }
  }

  return;
}

// ======================================================================
void readFile(TString filename, TFile *&file) {
  file = TFile::Open(filename);
}

// ======================================================================
// READ CF (SAME → NUM, MIXED → DEN)
// ======================================================================
void readCF_generic(TFile *filein,
                    int evtType,
                    int phiBin,
                    bool usePhi,
                    vector<vector<vector<TH3F *>>> &numStore,
                    vector<vector<TH3F *>> &denStore) {
  TDirectoryFile *base = (TDirectoryFile *)filein->Get("femto-dream-pair-task-track-track");
  TDirectoryFile *dir = (TDirectoryFile *)base->Get(evtType == kSAME ? "SameEvent_3Dqn" : "MixedEvent_3Dqn");

  THnSparseF *sp = (THnSparseF *)dir->Get("relPair3dRmTMultPercentileQnPairphi");

  for (int iMt = 0; iMt < nMtBins; iMt++) {
    sp->GetAxis(kMt)->SetRangeUser(mTBin[iMt] + 1e-3, mTBin[iMt + 1] - 1e-3);

    for (int iCent = 0; iCent < nCentBins; iCent++) {
      sp->GetAxis(kCent)->SetRangeUser(centBin[iCent], centBin[iCent + 1]);

      // φ-bin 选择
      if (usePhi) {
        sp->GetAxis(kPhi)->SetRange(phiBin, phiBin);
      } else {
        double phipair1 = 0. * TMath::Pi() / 180.;
        double phipair2 = 0. * TMath::Pi() / 180.;
        // int pb1 = sp->GetAxis(kPhi)->FindBin(phipair1);
        // int pb2 = sp->GetAxis(kPhi)->FindBin(phipair2);
        int pb1 = phibinCheck;
        int pb2 = phibinCheck;
        cout << "we constrain the phi(pair) at: " << pb1 << "  to  " << pb2 << endl;
        sp->GetAxis(kPhi)->SetRange(pb1, pb2);
      }

      if (evtType == kMIXED) {
        sp->GetAxis(kQn)->SetRange(1, 12);

        sp->GetAxis(kOut)->SetRangeUser(-2.0 + 0.00001, 2.0 - 0.00001);
        sp->GetAxis(kSide)->SetRangeUser(-2.0 + 0.00001, 2.0 - 0.00001);
        sp->GetAxis(kLong)->SetRangeUser(-2.0 + 0.00001, 2.0 - 0.00001);

        TH3F *h3 = (TH3F *)sp->Projection(kOut, kSide, kLong, "O");
        TString name = Form("H3_mixed_cent%d_mt%d_phi%d", iCent, iMt, phiBin - 1);
        h3->SetName(name);
        h3->SetDirectory(0);
        denStore[iMt][iCent] = h3;

      } else {
        for (int iqn = 0; iqn < nQnBins; iqn++) {
          if (evtType == kSAME) {
            if (iqn == 0)
              sp->GetAxis(kQn)->SetRange(1, 10);
            // if(iqn==1) sp->GetAxis(kQn)->SetRange(4,7);
            // if(iqn==2) sp->GetAxis(kQn)->SetRange(8,10);
            // sp->GetAxis(kQn)->SetRange(1,10);
          }
          sp->GetAxis(kOut)->SetRangeUser(-2.0 + 0.00001, 2.0 - 0.00001);
          sp->GetAxis(kSide)->SetRangeUser(-2.0 + 0.00001, 2.0 - 0.00001);
          sp->GetAxis(kLong)->SetRangeUser(-2.0 + 0.00001, 2.0 - 0.00001);

          TString name = Form("H3_same_cent%d_mt%d_qn%d_phi%d", iCent, iMt, iqn, phiBin - 1);
          TH3F *h3 = (TH3F *)sp->Projection(kOut, kSide, kLong, "O");
          h3->SetName(name);
          h3->SetDirectory(0);
          numStore[iMt][iCent][iqn] = h3;
        }
      }
    }
  }
  sp->GetAxis(kOut)->SetRange(0, -1);
  sp->GetAxis(kSide)->SetRange(0, -1);
  sp->GetAxis(kLong)->SetRange(0, -1);
  sp->GetAxis(kMt)->SetRange(0, -1);
  sp->GetAxis(kCent)->SetRange(0, -1);
  sp->GetAxis(kQn)->SetRange(0, -1);
  sp->GetAxis(kPhi)->SetRange(0, -1);
}

// ======================================================================
void readCF(TFile *filein, int type) {
  readCF_generic(filein, type, 0, false, histoNum, histoDen);
}

// ======================================================================
void readCF_phi(TFile *filein, int type, int phiBin) {
  readCF_generic(filein, type, phiBin, true, histoNum_phi[phiBin - 1], histoDen_phi[phiBin - 1]);
}

// ======================================================================
// CORRECT CF PROJECTION: CUT → PROJECT → RATIO
// ======================================================================
// binsFromZero 就是“从 0 bin 往两边取多少个 bin”，类似之前你用的 2
void projectCF1D(TH1F *&hCF1D,
                 TH1F *&hCF1DBkg,
                 TH3F *hNum3D,
                 TH3F *hDen3D,
                 TString target,
                 TF1 *&fNormBkg,
                 int qn,
                 int binsFromZero = 1) {
  TH1F *hNum1D = nullptr;
  TH1F *hDen1D = nullptr;

  hNum3D->Scale(1. / hNum3D->Integral("width"));
  if (qn == 0)
    hDen3D->Scale(1. / hDen3D->Integral("width"));

  // 开启误差
  // hNum1D->Sumw2();
  // hDen1D->Sumw2();

  hNum3D->Divide(hDen3D);

  if (target.EqualTo("out")) {
    // 投影到 X（q_out），在 Y/Z (side/long) 上切片

    int y1 = hNum3D->GetYaxis()->FindBin(-0.06);
    int y2 = hNum3D->GetYaxis()->FindBin(0.06);
    int z1 = hNum3D->GetZaxis()->FindBin(-0.06);
    int z2 = hNum3D->GetZaxis()->FindBin(0.06);

    // 和 CorrelationFunctionDataSimple::_apply_projection 完全同款调用
    hNum1D = (TH1F *)hNum3D->ProjectionX(Form("NumOut_%d", qn), y1, y2, z1, z2);
    hDen1D = (TH1F *)hDen3D->ProjectionX(Form("DenOut_%d", qn), y1, y2, z1, z2);

    int numBin = (y2 - y1 + 1) * (z2 - z1 + 1);
    hNum1D->Scale(1. / numBin);
  } else if (target.EqualTo("side")) {
    // 投影到 Y（q_side），在 X/Z 上切片

    int x1 = hNum3D->GetXaxis()->FindBin(-0.06);
    int x2 = hNum3D->GetXaxis()->FindBin(0.06);
    int z1 = hNum3D->GetZaxis()->FindBin(-0.06);
    int z2 = hNum3D->GetZaxis()->FindBin(0.06);

    hNum1D = (TH1F *)hNum3D->ProjectionY(Form("NumSide_%d", qn), x1, x2, z1, z2);
    hDen1D = (TH1F *)hDen3D->ProjectionY(Form("DenSide_%d", qn), x1, x2, z1, z2);

    int numBin = (x2 - x1 + 1) * (z2 - z1 + 1);
    hNum1D->Scale(1. / numBin);
  } else if (target.EqualTo("long")) {
    // 投影到 Z（q_long），在 X/Y 上切片

    int x1 = hNum3D->GetXaxis()->FindBin(-0.06);
    int x2 = hNum3D->GetXaxis()->FindBin(0.06);
    int y1 = hNum3D->GetYaxis()->FindBin(-0.06);
    int y2 = hNum3D->GetYaxis()->FindBin(0.06);

    hNum1D = (TH1F *)hNum3D->ProjectionZ(Form("NumLong_%d", qn), x1, x2, y1, y2);
    hDen1D = (TH1F *)hDen3D->ProjectionZ(Form("DenLong_%d", qn), x1, x2, y1, y2);

    int numBin = (y2 - y1 + 1) * (x2 - x1 + 1);
    hNum1D->Scale(1. / numBin);
  }

  hCF1D = hNum1D;
  hCF1DBkg = hDen1D;

  // fNormBkg = (TF1*)normBkg(hCF1D, qn);
  // for (int i=0; i<hCF1D->GetNbinsX(); i++){
  //    	hCF1D->SetBinContent(i+1, hCF1D->GetBinContent(i+1)/fNormBkg->Eval(hCF1D->GetBinCenter(i+1)) );
  // }
}

// ======================================================================
TF1 *normBkg(TH1F *&hfCF1D, int iqn) {
  TH1F *hCopy = (TH1F *)hfCF1D->Clone("clone");
  hCopy->Reset();
  for (int i = 0; i < hCopy->GetNbinsX(); i++) {
    if (fabs(hfCF1D->GetBinCenter(i + 1)) >= 0.06) {
      hCopy->SetBinContent(i + 1, hfCF1D->GetBinContent(i + 1));
      hCopy->SetBinError(i + 1, hfCF1D->GetBinError(i + 1));
    }
  }

  TF1 *fLine = new TF1(Form("fbkgline_qn%d", iqn), "[0]", -1, 1);

  hCopy->Fit(fLine, "NMR+0");

  new TCanvas;
  hCopy->Draw("pl");
  fLine->Draw("same l");

  return fLine;
}

float norm(TH1F *&hfCF1D) {
  TH1F *hCopy = (TH1F *)hfCF1D->Clone("clone");
  hCopy->Reset();
  for (int i = 0; i < hCopy->GetNbinsX(); i++) {
    if (fabs(hfCF1D->GetBinCenter(i + 1)) >= 0.15) {
      hCopy->SetBinContent(i + 1, hfCF1D->GetBinContent(i + 1));
      hCopy->SetBinError(i + 1, hfCF1D->GetBinError(i + 1));
    }
  }
  TF1 *fLine = new TF1("fline", "[0]+[1]*x+[2]*x*x", -2, 2);
  hCopy->Fit(fLine, "NMR+0");
  return fLine->GetParameter(0);
}

// ======================================================================
void drawAndCheck(TString drawOpt,
                  vector<int> color,
                  int canvaswidth,
                  int canvasHight,
                  bool isDivided,
                  bool isMultPlots,
                  int divideX,
                  int divideY,
                  TH1F *hDrawSingle,
                  vector<TH1F *> hDraw,
                  vector<TF1 *> fBkg = {}) {
  TCanvas *c1 = new TCanvas("c1", "c1", canvaswidth, canvasHight);

  if (!isMultPlots) {
    SetHistoSytle(hDrawSingle, color[0], 20);
    hDrawSingle->Draw(drawOpt);
    return;
  }

  TLegend *lg = new TLegend(0.6, 0.55, 0.9, 0.85, "Pb#minusPb, #pi^{+}-#pi^{+}");
  lg->SetFillStyle(0);
  lg->SetBorderSize(0);
  lg->SetTextSize(0.04);

  for (int i = 0; i < hDraw.size(); i++) {
    SetHistoSytle(hDraw[i], color[i], 29);
    hDraw[i]->GetXaxis()->SetTitle("q^{side}");
    hDraw[i]->GetYaxis()->SetTitle("C(q^{side})");
    hDraw[i]->GetYaxis()->SetRangeUser(0.95, 1.2);
    int rc = GetRootColor(color[i]);
    // fBkg[i]->SetLineColor(rc);
    // fBkg[i]->SetLineStyle(1);
    // fBkg[i]->SetLineWidth(4);
    if (i == 0) {
      hDraw[i]->Draw(drawOpt);
      // fBkg[i]->Draw("same l");
    } else {
      hDraw[i]->Draw(Form("same %s", drawOpt.Data()));
      // fBkg[i]->Draw("same l");
    }

    lg->AddEntry(hDraw[i], Form("qn %% %.f-%.f", qnBin[i], qnBin[i + 1]), "pl");
  }

  // lg->Draw();

  //     // -----------------------------------------------------
  //     // write
  //     // -----------------------------------------------------
  //     TFile* fout = new TFile("output/cf_side_cent30to50_mt0.2to0.3_qn_example.root","RECREATE");
  //     fout->cd();
  // for (int i=0; i<hDraw.size(); i++){
  // 	hDraw[i]->Write();
  // }
  //     fout->Close();
}

// ======================================================================
void SetHistoSytle(TH1F *h, int color, int style) {
  h->SetMarkerSize(1.3);
  h->SetLineWidth(1);
  h->SetMarkerStyle(style);
  int rc = GetRootColor(color);
  h->SetMarkerColor(rc);
  h->SetLineColor(rc);
}

// ======================================================================
void SetStyle(Bool_t graypalette) {
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  if (graypalette)
    gStyle->SetPalette(8, 0);
  else
    gStyle->SetPalette(1);

  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetLabelSize(0.04, "xyz");
  gStyle->SetTitleSize(0.05, "xyz");
  gStyle->SetTextFont(42);
}

int GetRootColor(int c) {
  if (c == kRe)
    return kRed - 4;
  if (c == kBl)
    return kAzure + 1;
  if (c == kGr)
    return kGreen - 2;
  if (c == kPk)
    return kPink + 8;
  if (c == kPur)
    return kViolet - 1;
  if (c == kYl)
    return kYellow - 3;
  if (c == kGry)
    return kGray + 1;
  if (c == kBlk)
    return kBlack;
  return kBlack;
}
