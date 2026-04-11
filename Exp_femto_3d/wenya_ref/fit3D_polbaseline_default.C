
#ifndef BOWLER3DFIT_POLBASELINE_C
#define BOWLER3DFIT_POLBASELINE_C

#include "Bowler3DFit_polbaseline.C"

#endif // BOWLER3DFIT_POLBASELINE_C

vector<vector<vector<vector<TH3F*>>>> histoSE;
vector<vector<vector<TH3F*>>> histoME;
vector<vector<vector<TH1F*>>> histoRoutVsPhi;
vector<vector<vector<TH1F*>>> histoRsideVsPhi;
vector<vector<vector<TH1F*>>> histoRlongVsPhi;
vector<vector<vector<TH1F*>>> histolamVsPhi;
vector<vector<vector<TH1F*>>> histoN0VsPhi;
vector<vector<vector<TH1F*>>> histoN1VsPhi;

double centBin[] = {0,10,30,50,80};
double mTBin[]   = {0.2,0.3,0.5,0.7};
double qnBin[]   = {0,30,70,100};

// const int nCentBins = sizeof(centBin)/sizeof(double)-1;
// const int nMtBins   = sizeof(mTBin)/sizeof(double)-1;
// const int nQnBins   = sizeof(qnBin)/sizeof(double)-1;
const int nCentBins = 3;
const int nMtBins   = 1;

const int nQnBins   = 1;
const int nPhiBins  = 6;

void fit3D_polbaseline(int mod=0){
	TFile* file = TFile::Open("output/subsamp/out_3DCF_mT_cent_qn_phi_samp13.root","READ");

	histoSE.clear();
	histoSE.resize(nCentBins, vector<vector<vector<TH3F*>>>(nMtBins, vector<vector<TH3F*>>(nQnBins, vector<TH3F*>(nPhiBins,nullptr))));

	histoME.clear();
	histoME.resize(nCentBins, vector<vector<TH3F*>>(nMtBins, vector<TH3F*>(nPhiBins, nullptr)));
							
	double lam[nCentBins][nMtBins][nQnBins][nPhiBins];
	double Ro2[nCentBins][nMtBins][nQnBins][nPhiBins];
	double Rs2[nCentBins][nMtBins][nQnBins][nPhiBins];
	double Rl2[nCentBins][nMtBins][nQnBins][nPhiBins];
	double n0[nCentBins][nMtBins][nQnBins][nPhiBins];
	double n1[nCentBins][nMtBins][nQnBins][nPhiBins];

	double lamerr[nCentBins][nMtBins][nQnBins][nPhiBins];
	double Ro2err[nCentBins][nMtBins][nQnBins][nPhiBins];
	double Rs2err[nCentBins][nMtBins][nQnBins][nPhiBins];
	double Rl2err[nCentBins][nMtBins][nQnBins][nPhiBins];
	double n0err[nCentBins][nMtBins][nQnBins][nPhiBins];
	double n1err[nCentBins][nMtBins][nQnBins][nPhiBins];


	for (int iCent=2; iCent<nCentBins; iCent++){
	for (int iMt=0; iMt<nMtBins; iMt++){
	for (int iqn=0; iqn<nQnBins; iqn++){				
	for (int iPhi=5; iPhi<nPhiBins; iPhi++){					
			histoSE[iCent][iMt][iqn][iPhi] = (TH3F*)file->Get(Form("H3_same_cent%d_mt%d_qn%d_phi%d", iCent, iMt, iqn, iPhi));
		}
	}}}


	for (int iCent=2; iCent<nCentBins; iCent++){
	for (int iMt=0; iMt<nMtBins; iMt++){
	for (int iPhi=5; iPhi<nPhiBins; iPhi++){					
			histoME[iCent][iMt][iPhi] = (TH3F*)file->Get(Form("H3_mixed_cent%d_mt%d_phi%d", iCent, iMt, iPhi));
		}
	}}


	for (int iCent=2; iCent<nCentBins; iCent++){
	for (int iMt=0; iMt<nMtBins; iMt++){
	for (int iqn=0; iqn<nQnBins; iqn++){				
	for (int iPhi=5; iPhi<nPhiBins; iPhi++){					
			Fit3D_BowlerSinyukov_Full((int)centBin[iCent], (double)mTBin[iMt], histoSE[iCent][iMt][iqn][iPhi], histoME[iCent][iMt][iPhi], iqn, iPhi,
									  lam[iCent][iMt][iqn][iPhi], Ro2[iCent][iMt][iqn][iPhi], Rs2[iCent][iMt][iqn][iPhi], Rl2[iCent][iMt][iqn][iPhi], 
									  lamerr[iCent][iMt][iqn][iPhi], Ro2err[iCent][iMt][iqn][iPhi], Rs2err[iCent][iMt][iqn][iPhi], Rl2err[iCent][iMt][iqn][iPhi],
									  n0[iCent][iMt][iqn][iPhi], n1[iCent][iMt][iqn][iPhi],
									  n0err[iCent][iMt][iqn][iPhi], n1err[iCent][iMt][iqn][iPhi],
									  mod);		}
	}}}

	// histoRoutVsPhi.clear();
	// histoRoutVsPhi.resize(nCentBins, vector<vector<TH1F*>>(nMtBins, vector<TH1F*>(nQnBins, nullptr)));

	// for (int iCent=0; iCent<nCentBins; iCent++){
	// for (int iMt=0; iMt<nMtBins; iMt++){
	// for (int iqn=0; iqn<nQnBins; iqn++){
	// 	histoRoutVsPhi[iCent][iMt][iqn] = new TH1F(Form("ro2_cent%d_mt%d_qn%d",iCent,iMt,iqn),";#phi^{pair} - #Psi_EP;R_{out}^{2}", 12, -0.2, 3.4);				
	// 	for (int iPhi=0; iPhi<nPhiBins; iPhi++){
	// 		histoRoutVsPhi[iCent][iMt][iqn]->SetBinContent(iPhi+1, Ro2[iCent][iMt][iqn][iPhi]);
	// 		histoRoutVsPhi[iCent][iMt][iqn]->SetBinError(iPhi+1, Ro2err[iCent][iMt][iqn][iPhi]);
	// 	}
	// }}}


	// histoRsideVsPhi.clear();
	// histoRsideVsPhi.resize(nCentBins, vector<vector<TH1F*>>(nMtBins, vector<TH1F*>(nQnBins, nullptr)));

	// for (int iCent=0; iCent<nCentBins; iCent++){
	// for (int iMt=0; iMt<nMtBins; iMt++){
	// for (int iqn=0; iqn<nQnBins; iqn++){
	// 	histoRsideVsPhi[iCent][iMt][iqn] = new TH1F(Form("rs2_cent%d_mt%d_qn%d",iCent,iMt,iqn),";#phi^{pair} - #Psi_EP;R_{out}^{2}", 12, -0.2, 3.4);				
	// 	for (int iPhi=0; iPhi<nPhiBins; iPhi++){
	// 		histoRsideVsPhi[iCent][iMt][iqn]->SetBinContent(iPhi+1, Rs2[iCent][iMt][iqn][iPhi]);
	// 		histoRsideVsPhi[iCent][iMt][iqn]->SetBinError(iPhi+1, Rs2err[iCent][iMt][iqn][iPhi]);
	// 	}
	// }}}


	// histoRlongVsPhi.clear();
	// histoRlongVsPhi.resize(nCentBins, vector<vector<TH1F*>>(nMtBins, vector<TH1F*>(nQnBins, nullptr)));

	// for (int iCent=0; iCent<nCentBins; iCent++){
	// for (int iMt=0; iMt<nMtBins; iMt++){
	// for (int iqn=0; iqn<nQnBins; iqn++){
	// 	histoRlongVsPhi[iCent][iMt][iqn] = new TH1F(Form("rl2_cent%d_mt%d_qn%d",iCent,iMt,iqn),";#phi^{pair} - #Psi_EP;R_{out}^{2}", 12, -0.2, 3.4);				
	// 	for (int iPhi=0; iPhi<nPhiBins; iPhi++){
	// 		histoRlongVsPhi[iCent][iMt][iqn]->SetBinContent(iPhi+1, Rl2[iCent][iMt][iqn][iPhi]);
	// 		histoRlongVsPhi[iCent][iMt][iqn]->SetBinError(iPhi+1, Rl2err[iCent][iMt][iqn][iPhi]);
	// 	}
	// }}}
		
	// // -----------------------------------------------------
	// // write
	// // -----------------------------------------------------
	// TFile* fout = new TFile("output/radii_cent_mt_qn_polbaseline_fitRange_2.root","RECREATE");
	// fout->cd();
	// for (int iCent=0; iCent<nCentBins; iCent++){
	// for (int iMt=0; iMt<nMtBins; iMt++){													
	// for (int iqn=0; iqn<nQnBins; iqn++){					
	// 	histoRoutVsPhi[iCent][iMt][iqn]->Write();
	// 	histoRsideVsPhi[iCent][iMt][iqn]->Write();
	// 	histoRlongVsPhi[iCent][iMt][iqn]->Write();
	// }}}	

	// histolamVsPhi.clear();
	// histolamVsPhi.resize(nCentBins, vector<vector<TH1F*>>(nMtBins, vector<TH1F*>(nQnBins, nullptr)));

	// for (int iCent=0; iCent<nCentBins; iCent++){
	// for (int iMt=0; iMt<nMtBins; iMt++){
	// for (int iqn=0; iqn<nQnBins; iqn++){
	// 	histolamVsPhi[iCent][iMt][iqn] = new TH1F(Form("lam_cent%d_mt%d_qn%d",iCent,iMt,iqn),";#phi^{pair} - #Psi_EP;R_{out}^{2}", 12, -0.2, 3.4);				
	// 	for (int iPhi=0; iPhi<nPhiBins; iPhi++){
	// 		histolamVsPhi[iCent][iMt][iqn]->SetBinContent(iPhi+1, lam[iCent][iMt][iqn][iPhi]);
	// 		histolamVsPhi[iCent][iMt][iqn]->SetBinError(iPhi+1, lamerr[iCent][iMt][iqn][iPhi]);
	// 	}
	// }}}
	// // -----------------------------------------------------
	// // write
	// // -----------------------------------------------------
	// TFile* fout = new TFile("output/lambda_mT0.3to0.5_cent10to30_polbaseline_allmomenta.root","RECREATE");
	// fout->cd();
	// for (int iCent=0; iCent<nCentBins; iCent++){
	// for (int iMt=0; iMt<nMtBins; iMt++){													
	// for (int iqn=0; iqn<nQnBins; iqn++){					
	// 	histolamVsPhi[iCent][iMt][iqn]->Write();
	// }}}


	// histoN0VsPhi.clear();
	// histoN0VsPhi.resize(nCentBins, vector<vector<TH1F*>>(nMtBins, vector<TH1F*>(nQnBins, nullptr)));

	// for (int iCent=0; iCent<nCentBins; iCent++){
	// for (int iMt=0; iMt<nMtBins; iMt++){
	// for (int iqn=0; iqn<nQnBins; iqn++){
	// 	histoN0VsPhi[iCent][iMt][iqn] = new TH1F(Form("n0_cent%d_mt%d_qn%d",iCent,iMt,iqn),";#phi^{pair} - #Psi_EP;R_{out}^{2}", 12, -0.2, 3.4);				
	// 	for (int iPhi=0; iPhi<nPhiBins; iPhi++){
	// 		histoN0VsPhi[iCent][iMt][iqn]->SetBinContent(iPhi+1, n0[iCent][iMt][iqn][iPhi]);
	// 		histoN0VsPhi[iCent][iMt][iqn]->SetBinError(iPhi+1, n0err[iCent][iMt][iqn][iPhi]);
	// 	}
	// }}}


	// histoN1VsPhi.clear();
	// histoN1VsPhi.resize(nCentBins, vector<vector<TH1F*>>(nMtBins, vector<TH1F*>(nQnBins, nullptr)));

	// for (int iCent=0; iCent<nCentBins; iCent++){
	// for (int iMt=0; iMt<nMtBins; iMt++){
	// for (int iqn=0; iqn<nQnBins; iqn++){
	// 	histoN1VsPhi[iCent][iMt][iqn] = new TH1F(Form("n1_cent%d_mt%d_qn%d",iCent,iMt,iqn),";#phi^{pair} - #Psi_EP;R_{out}^{2}", 12, -0.2, 3.4);				
	// 	for (int iPhi=0; iPhi<nPhiBins; iPhi++){
	// 		histoN1VsPhi[iCent][iMt][iqn]->SetBinContent(iPhi+1, n1[iCent][iMt][iqn][iPhi]);
	// 		histoN1VsPhi[iCent][iMt][iqn]->SetBinError(iPhi+1, n1err[iCent][iMt][iqn][iPhi]);
	// 	}
	// }}}
	// // -----------------------------------------------------
	// // write
	// // -----------------------------------------------------
	// TFile* foutn = new TFile("output/baseline_cent_mt_qn_basebasebase.root","RECREATE");
	// foutn->cd();
	// for (int iCent=0; iCent<nCentBins; iCent++){
	// for (int iMt=0; iMt<nMtBins; iMt++){													
	// for (int iqn=0; iqn<nQnBins; iqn++){					
	// 	histoN0VsPhi[iCent][iMt][iqn]->Write();
	// 	histoN1VsPhi[iCent][iMt][iqn]->Write();
	// }}}	
}