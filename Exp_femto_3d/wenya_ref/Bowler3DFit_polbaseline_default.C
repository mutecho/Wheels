// Bowler3DFit_polbaseline.C
// 3D Bowler–Sinyukov fit to C(q_out,q_side,q_long)
// 模型：C(q) = N(1+bq^2) [ λ G(q) K_coul(qinv) + (1 - λ) ]
// G(q) = 1 + exp( -Ro2*q_out^2 - Rs2*q_side^2 - Rl2*q_long^2
//                  -2*Ros2*q_out*q_side -2*Rol2*q_out*q_long -2*Rsl2*q_side*q_long )

#include <TH3F.h>
#include <TMinuit.h>
#include <TMath.h>
#include <iostream>
#include <map>

const double hbarc = 0.19732698;

// 全局指针：要 fit 的 3D CF (Num & Den)
static TH3F* gNum3D = nullptr;   // SAME
static TH3F* gDen3D = nullptr;   // MIXED
static double gQmax = 0.2;
static double gBestParForSlice[9];
static int gNPoints = 0;

TFile*  file_coul;
TSpline3* spl;
TFile*  file_KStar;
TH1F* hKStar;
TFile*  para_baseline;

static TFile* GetKStarFile()
{
  if (!file_KStar) {
    file_KStar = TFile::Open("./output/out_dist_kstar_same.root", "READ");
    if (!file_KStar || file_KStar->IsZombie()) {
      std::cerr << "ERROR: cannot open KStar file\n";
      if (file_KStar) { file_KStar->Close(); delete file_KStar; file_KStar = nullptr; }
    }
  }
  return file_KStar;
}

static void CloseKStarFile()
{
  if (file_KStar) {
    file_KStar->Close();
    delete file_KStar;
    file_KStar = nullptr;
  }
}

static std::map<int, TFile*>    gCoulFile;
static std::map<int, TSpline3*> gCoulSpl;

static TSpline3* GetCoulombSpline(int cent)
{
    if (cent != 0 && cent != 10 && cent != 30 && cent != 50) {
        std::cerr << "ERROR: GetCoulombSpline unknown cent=" << cent << std::endl;
        return nullptr;
    }

    // 已缓存 spline
    auto itS = gCoulSpl.find(cent);
    if (itS != gCoulSpl.end() && itS->second) return itS->second;

    const char* path = nullptr;
    if (cent == 0)  path = "./CRAB_version3/cf_kc_cent0to10.root";
    if (cent == 10) path = "./CRAB_version3/cf_kc_cent10to30.root";
    if (cent == 30) path = "./CRAB_version3/cf_kc_cent30to50.root";
    if (cent == 50) path = "./CRAB_version3/cf_kc_cent50to80.root";

    // 打开/复用文件
    TFile* f = nullptr;
    auto itF = gCoulFile.find(cent);
    if (itF != gCoulFile.end()) f = itF->second;

    if (!f) {
        f = TFile::Open(path, "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "ERROR: cannot open coulomb file: " << path << std::endl;
            if (f) { f->Close(); delete f; }
            gCoulFile[cent] = nullptr;
            return nullptr;
        }
        gCoulFile[cent] = f;
    }

    // 取 spline
    TSpline3* sp = (TSpline3*)f->Get("sp_cf_coulomb");
    if (!sp) {
        std::cerr << "ERROR: missing sp_cf_coulomb in " << f->GetName() << std::endl;
        return nullptr;
    }

    // Clone 一份放内存里（避免依赖文件内部对象生命周期）
    TSpline3* spClone = (TSpline3*)sp->Clone(Form("sp_cf_coulomb_cent%d", cent));
    gCoulSpl[cent] = spClone;
    return spClone;
}

static TSpline3* GetCoulombSpline1D(int cent)
{
  if (cent != 0 && cent != 10 && cent != 30 && cent != 50) {
    std::cerr << "ERROR: GetCoulombSpline1D unknown cent=" << cent << std::endl;
    return nullptr;
  }

  auto itS = gCoulSpl.find(cent);
  if (itS != gCoulSpl.end() && itS->second) return itS->second;

  const char* path = nullptr;
  // if (cent == 0)  path = "../CRAB_version3/cf_kc_cent0to10.root";
  // if (cent == 10) path = "../CRAB_version3/cf_kc_cent10to30.root";
  // if (cent == 30) path = "../CRAB_version3/cf_kc_cent30to50.root";
  // if (cent == 50) path = "../CRAB_version3/cf_kc_cent50to80.root";
  path = "./1DlamPara/coulomb.root";
  TFile* f = TFile::Open(path);
  f->ls();
  // auto itF = gCoulFile1D.find(cent);
  // if (itF != gCoulFile1D.end()) f = itF->second;

  // if (!f) {
  //   f = TFile::Open(path, "READ");
  //   if (!f || f->IsZombie()) {
  //     std::cerr << "ERROR: cannot open coulomb file: " << path << std::endl;
  //     if (f) { f->Close(); delete f; }
  //     gCoulFile1D[cent] = nullptr;
  //     return nullptr;
  //   }
  //   gCoulFile1D[cent] = f;
  // }

  TH1D* h_coulomb = (TH1D*)f->Get("h_coulomb");
  auto g = new TGraph(h_coulomb->GetNbinsX());
  int ip = 0;
  for (int i = 1; i <= h_coulomb->GetNbinsX(); ++i) {
    const double x = h_coulomb->GetBinCenter(i);
    const double y = h_coulomb->GetBinContent(i);
    g->SetPoint(ip++, x, y);
  }

  TSpline3* sp = new TSpline3("sp_coulomb", g);
  if (!sp) {
    std::cerr << "ERROR: missing sp_cf_coulomb in " << f->GetName() << std::endl;
    return nullptr;
  }

  // TSpline3* spClone = (TSpline3*)sp->Clone(Form("sp_cf_coulomb_cent%d_1D", cent));
  // gCoulSpl[cent] = spClone;
  return sp;
}

static void CloseCoulombFiles()
{
    // 先删 spline（clone 出来的归我们管）
    for (auto &kv : gCoulSpl) {
        delete kv.second;
        kv.second = nullptr;
    }
    gCoulSpl.clear();

    // 再关文件
    for (auto &kv : gCoulFile) {
        if (kv.second) {
            kv.second->Close();
            delete kv.second;
            kv.second = nullptr;
        }
    }
    gCoulFile.clear();
}

//--------------------------------------------------------
// Coulomb factor K_coul(q_inv) from Crab
//--------------------------------------------------------
double Kcoul_Crab(double qinv)
{
  if (qinv < 1e-6) qinv = 1e-6;    // 避免除零
  double K = spl->Eval(qinv);
  return K;
}

//--------------------------------------------------------
// Bowler–Sinyukov 3D 模型 (严格按 note 的 G(q) + C(q))
// par:
//   0-1: N(q)
//   2: lambda
//   3: Ro2   = R_o^2
//   4: Rs2   = R_s^2
//   5: Rl2   = R_l^2
//   6: Ros2  = R_os^2
//   7: Rol2  = R_ol^2
//   8: Rsl2  = R_sl^2
//--------------------------------------------------------
double CF_BowlerSinyukov_3D(double qout, double qside, double qlong,
                            const double *par)
{
    // baseline coefficients
    double a0 = par[0];
    double a1 = par[1];
    // femto parameters
    int startSourcePar = 2;
    double lambda = par[startSourcePar];
    double Ro2    = par[startSourcePar+1];
    double Rs2    = par[startSourcePar+2];
    double Rl2    = par[startSourcePar+3];
    double Ros2   = par[startSourcePar+4];
    double Rol2   = par[startSourcePar+5];
    double Rsl2   = par[startSourcePar+6];

    double qinvLow = TMath::Sqrt((fabs(qout)-0.01)*(fabs(qout)-0.01) + (fabs(qside)-0.01)*(fabs(qside)-0.01) + (fabs(qlong)-0.01)*(fabs(qlong)-0.01));
    double qinvUp = TMath::Sqrt((fabs(qout)+0.01)*(fabs(qout)+0.01) + (fabs(qside)+0.01)*(fabs(qside)+0.01) + (fabs(qlong)+0.01)*(fabs(qlong)+0.01));
    
    int kbinlow = hKStar->GetXaxis()->FindBin(0.5*qinvLow);
    int kbinup = hKStar->GetXaxis()->FindBin(0.5*qinvUp);
    
    hKStar->GetXaxis()->SetRange(kbinlow, kbinup);
    
    double kstarMean = hKStar->GetMean(1);
    
    double Kc   = Kcoul_Crab(kstarMean);

    double factor = 1.0 / (hbarc*hbarc);

    double exponent =
        - ( Ro2*qout*qout*factor
          + Rs2*qside*qside*factor
          + Rl2*qlong*qlong*factor
          + 2.0*Ros2*qout*qside*factor
          + 2.0*Rol2*qout*qlong*factor
          + 2.0*Rsl2*qside*qlong*factor );

    double G = 1.0 + TMath::Exp(exponent);

    double Nq = a0*(1 + a1*(2*kstarMean)*(2*kstarMean));

    double C = Nq * ( lambda * G * Kc + (1.0 - lambda) );
    return C;
}


//--------------------------------------------------------
// 给 TF1 用的 1D 切片函数（不能 capture lambda）
//--------------------------------------------------------
double fOut_1D(double* x, double*) {
    double qout = x[0];
    return (CF_BowlerSinyukov_3D(qout, -0.011, 0.0, gBestParForSlice) + 
            CF_BowlerSinyukov_3D(qout, 0.011, 0.0, gBestParForSlice) + 
            CF_BowlerSinyukov_3D(qout, 0.0, -0.011, gBestParForSlice) + 
            CF_BowlerSinyukov_3D(qout, 0.0, 0.011, gBestParForSlice))/4.;
    // return CF_BowlerSinyukov_3D(qout, -0.19, -0.19, gBestParForSlice);
}
double fSide_1D(double* x, double*) {
    double qside = x[0];
    return
        (CF_BowlerSinyukov_3D(-0.011, qside, 0.0, gBestParForSlice) + 
            CF_BowlerSinyukov_3D(0.011, qside, 0.0, gBestParForSlice) + 
            CF_BowlerSinyukov_3D(0.0, qside, -0.011, gBestParForSlice) + 
            CF_BowlerSinyukov_3D(0.0, qside, 0.011, gBestParForSlice))/4.;
        // CF_BowlerSinyukov_3D(-0.19, qside, -0.19, gBestParForSlice);
}

double fLong_1D(double* x, double*) {
    double qlong = x[0];
    return
        (CF_BowlerSinyukov_3D(-0.011, 0.0, qlong, gBestParForSlice) + 
            CF_BowlerSinyukov_3D(0.011, 0.0, qlong, gBestParForSlice) + 
            CF_BowlerSinyukov_3D(0.0, -0.011, qlong, gBestParForSlice) + 
            CF_BowlerSinyukov_3D(0.0, 0.011, qlong, gBestParForSlice))/4.;
        // CF_BowlerSinyukov_3D(-0.19, -0.19, qlong, gBestParForSlice);
}


struct Chi2EffResult {
    double chi2;
    double neff;
    double ndf_eff;
    double chi2ndf_eff;
    int    nPoints;   // raw counted bins (for reference)
};

void fcn_Bowler3D(Int_t& npar, Double_t* grad, Double_t& f, Double_t* par, Int_t flag)
{
    if (!gNum3D || !gDen3D) {
        f = 1e30;
        return;
    }

    double chi2pml = 0.0;

    const int Nx = gNum3D->GetNbinsX();
    const int Ny = gNum3D->GetNbinsY();
    const int Nz = gNum3D->GetNbinsZ();

    for (int ix=1; ix<=Nx; ix++){
        double qout = gNum3D->GetXaxis()->GetBinCenter(ix);
        if (fabs(qout) > gQmax) continue;

        for (int iy=1; iy<=Ny; iy++){
            double qside = gNum3D->GetYaxis()->GetBinCenter(iy);
            if (fabs(qside) > gQmax) continue;

            for (int iz=1; iz<=Nz; iz++){
                double qlong = gNum3D->GetZaxis()->GetBinCenter(iz);
                if (fabs(qlong) > gQmax) continue;

                double A = gNum3D->GetBinContent(ix,iy,iz);  // SAME
                double B = gDen3D->GetBinContent(ix,iy,iz);  // MIXED

                if (A <= 0 || B <= 0) continue;

                // 模型 SAME/MIX ratio
                double C = CF_BowlerSinyukov_3D(qout,qside,qlong,par);

                if (C <= 0) C = 1e-12;

                // =======  PML (Poisson Maximum Likelihood)  ========
                // 来自 femto note 的公式
                double term1 = A * log( C*(A+B) / (A*(C+1)) );
                double term2 = B * log( (A+B) / (B*(C+1)) );

                chi2pml += -2.0*(term1 + term2);
            }
        }
    }
    f = chi2pml;
}



//--------------------------------------------------------
// 拟合入口：
//
//   hNum :  3D same
//   hDen :  3D mixed
//   qmax : 在 |q_out|,|q_side|,|q_long|<qmax 区域里拟合
//--------------------------------------------------------
void Fit3D_BowlerSinyukov_Full(int cent, double mT, TH3F* hNum, TH3F* hDen, int iqn, int iphi,
                                double& lam, double& ro2, double& rs2, double& rl2, 
                                double& lamerr, double& ro2err, double& rs2err, double& rl2err, 
                                // double& n0, double& n1, 
                                // double& n0Err, double& n1Err,
                                int mod=0, double qmax = 0.2)
{
  if (!hNum || !hDen) {
    std::cerr << "[Fit3D_BowlerSinyukov_Full] ERROR: nullptr histogram\n";
    return;
  }

  gNum3D = hNum;
  gDen3D = hDen;
  gQmax  = qmax;

  // spl = GetCoulombSpline(cent);
  spl = GetCoulombSpline1D(cent);
  if (!spl) return;  // 或者跳过本次 fit

  TFile* fk = GetKStarFile();
  if (!fk) return;

  if (cent == 0) {
    if (mT == 0.2)hKStar = (TH1F*)fk->Get(Form("hkstar_same_cent0_mt0"));
    if (mT == 0.3)hKStar = (TH1F*)fk->Get(Form("hkstar_same_cent0_mt1"));
    if (mT == 0.5)hKStar = (TH1F*)fk->Get(Form("hkstar_same_cent0_mt2"));
  }
  if (cent == 10) {
    if (mT == 0.2)hKStar = (TH1F*)fk->Get(Form("hkstar_same_cent1_mt0"));
    if (mT == 0.3)hKStar = (TH1F*)fk->Get(Form("hkstar_same_cent1_mt1"));
    if (mT == 0.5)hKStar = (TH1F*)fk->Get(Form("hkstar_same_cent1_mt2"));
  }
  if (cent == 30) {
    if (mT == 0.2)hKStar = (TH1F*)fk->Get(Form("hkstar_same_cent2_mt0"));
    if (mT == 0.3)hKStar = (TH1F*)fk->Get(Form("hkstar_same_cent2_mt1"));
    if (mT == 0.5)hKStar = (TH1F*)fk->Get(Form("hkstar_same_cent2_mt2"));
  }
  if (cent == 50) {
    if (mT == 0.2)hKStar = (TH1F*)fk->Get(Form("hkstar_same_cent3_mt0"));
    if (mT == 0.3)hKStar = (TH1F*)fk->Get(Form("hkstar_same_cent3_mt1"));
    if (mT == 0.5)hKStar = (TH1F*)fk->Get(Form("hkstar_same_cent3_mt2"));
  }

  // ---- 创建 Minuit，? 个参数 ----
  const int npar = 9;
  TMinuit minuit(npar);
  minuit.SetFCN(fcn_Bowler3D);     // PML FCN
  minuit.SetPrintLevel(-1);
  
  Int_t ierr = 0;
  
  // ========== 参数初始化（更适合 PML 的设置） ==========
  
  // par[0]: n0 + a1*q*q
  minuit.mnparm(0, "n0",         1.0,   0.01,  0,        100,   ierr);
  minuit.mnparm(1, "a1",         0.0,   0.01,  -100,     100,   ierr);

  int startSourcePara = 2;
  // par[3]: lambda
  minuit.mnparm(startSourcePara, "lambda", 0.6,   0.01,  0.1,   4.0,   ierr);
  // minuit.FixParameter(startSourcePara);  // Ros2
  
  // par[4]: Ro2 >= 0
  minuit.mnparm(startSourcePara+1, "Ro2",    40.0,  0.1,    0.0,   100.0, ierr);
  
  // par[5]: Rs2 >= 0
  minuit.mnparm(startSourcePara+2, "Rs2",    35.0,  0.1,    0.0,   100.0, ierr);
  
  // par[6]: Rl2 >= 0
  minuit.mnparm(startSourcePara+3, "Rl2",    45.0,  0.1,    0.0,   100.0, ierr);
  
  // ===== 交叉项必须允许负值，并允许大范围 =====
  minuit.mnparm(startSourcePara+4, "Ros2",    0.0,  0.01,   0.0, 0.0, ierr);
  minuit.mnparm(startSourcePara+5, "Rol2",    0.0,  0.01,   0.0, 0.0, ierr);
  minuit.mnparm(startSourcePara+6, "Rsl2",    0.0,  0.01,   0.0, 0.0, ierr);
  minuit.FixParameter(startSourcePara+4);  // Ros2
  minuit.FixParameter(startSourcePara+5);  // Rol2
  minuit.FixParameter(startSourcePara+6);  // Rsl2

  // ------------------------------------------------------
  Double_t arglist[2];
  arglist[0] = 100000;   // max function calls
  arglist[1] = 0.1;     // tolerance
  
  minuit.mnexcm("MIGRAD", arglist, 2, ierr);
  
  // 再 HESSE 求误差
  arglist[0] = 0;
  minuit.mnexcm("HESSE", arglist, 1, ierr);


  // store best-fit parameters for slice plotting
  for(int i=0;i<9;i++){
      double v,e;
      minuit.GetParameter(i, v, e);
      gBestParForSlice[i] = v;
  }

  minuit.GetParameter(startSourcePara, lam, lamerr);
  minuit.GetParameter(startSourcePara+1, ro2, ro2err);
  minuit.GetParameter(startSourcePara+2, rs2, rs2err);
  minuit.GetParameter(startSourcePara+3, rl2, rl2err);

  // minuit.GetParameter(0, n0, n0Err);
  // minuit.GetParameter(1, n1, n1Err);

  if (mod ==0){
      // ---- 输出结果 ----
      double val, err;
    
      std::cout << "\n===== 3D Bowler–Sinyukov Fit (Full, with cross terms) =====\n";
    
      minuit.GetParameter(0, val, err);
      std::cout << "N0      = " << val << "  ±  " << err << "\n";
      minuit.GetParameter(1, val, err);
      std::cout << "A1      = " << val << "  ±  " << err << "\n";

      minuit.GetParameter(startSourcePara, val, err);
      std::cout << "lambda = " << val << "  ±  " << err << "\n";
    
      minuit.GetParameter(startSourcePara+1, val, err);
      std::cout << "Ro2    = " << val << "  ±  " << err
                << "   => R_out = " << TMath::Sqrt(TMath::Max(val,0.0)) << "\n";
    
      minuit.GetParameter(startSourcePara+2, val, err);
      std::cout << "Rs2    = " << val << "  ±  " << err
                << "   => R_side = " << TMath::Sqrt(TMath::Max(val,0.0)) << "\n";
    
      minuit.GetParameter(startSourcePara+3, val, err);
      std::cout << "Rl2    = " << val << "  ±  " << err
                << "   => R_long = " << TMath::Sqrt(TMath::Max(val,0.0)) << "\n";
    
      minuit.GetParameter(startSourcePara+4, val, err);
      std::cout << "Ros2   = " << val << "  ±  " << err << "\n";
    
      minuit.GetParameter(startSourcePara+5, val, err);
      std::cout << "Rol2   = " << val << "  ±  " << err << "\n";
    
      minuit.GetParameter(startSourcePara+6, val, err);
      std::cout << "Rsl2   = " << val << "  ±  " << err << "\n";
    
      std::cout << "===========================================================\n";       

      // chi^2
      Double_t fmin, fedm, errdef;
      Int_t npari, nparx, istat;
  
      minuit.mnstat(fmin, fedm, errdef, npari, nparx, istat);
  
      std::cout << "====================================\n";
      std::cout << "iqn: " << iqn << "  iphi: " << iphi << "\n";
      std::cout << "PML minimum (-2 ln L) = " << fmin << "\n";
      std::cout << "Estimated distance to minimum (EDM) = " << fedm << "\n";
      std::cout << "Number of free parameters           = " << npari << "\n";
      std::cout << "Fit status (istat)                  = " << istat << "\n";
      std::cout << "====================================\n";
    }
    else if (mod==1){
      TCanvas* c = new TCanvas("c_fit3d_slices","3D Fit Slices",1200,400);
      c->Divide(3,1);
  
      // slice range = [-0.01,0.01]
      int x1 = hNum->GetXaxis()->FindBin(-0.01);
      int x2 = hNum->GetXaxis()->FindBin( 0.01);
      int y1 = hNum->GetYaxis()->FindBin(-0.01);
      int y2 = hNum->GetYaxis()->FindBin( 0.01);
      int z1 = hNum->GetZaxis()->FindBin(-0.01);
      int z2 = hNum->GetZaxis()->FindBin( 0.01);
      
      hNum->Divide(hDen);
      // ---------------- OUT slice ----------------
      c->cd(1);
      TH1F* hOut = (TH1F*)hNum->ProjectionX("slice_out", y1,y2, z1,z2);
      hOut->Scale(1./((y2-y1+1)*(z2-z1+1)));

      hOut->SetMarkerStyle(20);
      hOut->GetYaxis()->SetRangeUser(0.97,1.4);
      hOut->Draw("P");
  
      TF1* fOut = new TF1("fOut", fOut_1D, -gQmax, gQmax, 0);
      fOut->SetLineColor(kRed);
      fOut->Draw("same");
  
      // ---------------- SIDE slice ----------------
      c->cd(2);
      TH1F* hSide = (TH1F*)hNum->ProjectionY("slice_side", x1,x2, z1,z2);
      hSide->Scale(1./((x2-x1+1)*(z2-z1+1)));

      hSide->SetMarkerStyle(20);
      hSide->GetYaxis()->SetRangeUser(0.97,1.4);
      hSide->Draw("P");
  
      TF1* fSide = new TF1("fSide", fSide_1D, -gQmax, gQmax, 0);
      fSide->SetLineColor(kRed);
      fSide->Draw("same");
  
      // ---------------- LONG slice ----------------
      c->cd(3);
      TH1F* hLong = (TH1F*)hNum->ProjectionZ("slice_long", x1,x2, y1,y2);
      hLong->Scale(1./((x2-x1+1)*(y2-y1+1)));

      hLong->SetMarkerStyle(20);
      hLong->GetYaxis()->SetRangeUser(0.97,1.4);
      hLong->Draw("P");
  
      TF1* fLong = new TF1("fLong", fLong_1D, -gQmax, gQmax, 0);
      fLong->SetLineColor(kRed);
      fLong->Draw("same");
  
      c->Update(); 

      // // -----------------------------------------------------
      // // write
      // // -----------------------------------------------------
      // TFile* fout = new TFile("output/fit_Npol_cent30to50_mt0.2to0.3_qn0_fit_example.root","RECREATE");
      // fout->cd();
      // hOut->Write("hOut");      
      // hSide->Write("hSide");      
      // hLong->Write("hLong");      
      // fOut->Write("fOut");
      // fSide->Write("fSide");
      // fLong->Write("fLong");
      // fout->Close();     
    }

}

