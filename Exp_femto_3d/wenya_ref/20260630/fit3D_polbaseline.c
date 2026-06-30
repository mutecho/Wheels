

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
vector<vector<vector<TH1F*>>> histoNoVsPhi;
vector<vector<vector<TH1F*>>> histoNsVsPhi;
vector<vector<vector<TH1F*>>> histoNlVsPhi;
vector<vector<vector<TH1F*>>> histoNosVsPhi;
vector<vector<vector<TH1F*>>> histoNolVsPhi;
vector<vector<vector<TH1F*>>> histoNslVsPhi;

// double centBin[] = {0,10,30,50};
// double mTBin[]   = {0.2,0.3};
double centBin[] = {50,80};
double mTBin[]   = {0.2,0.3};
double qnBin[]   = {0,30,70,100};

const int nCentBins = sizeof(centBin)/sizeof(double)-1;
const int nMtBins   = sizeof(mTBin)/sizeof(double)-1;
// const int nQnBins   = sizeof(qnBin)/sizeof(double)-1;
const int nQnBins   = 1;
const int nPhiBins  = 1;

void fit3D_polbaseline(int mod=0){
	TFile* file = TFile::Open("output/out_3DCF_mT_cent_qn_phi.root","READ");

	histoSE.clear();
	histoSE.resize(nCentBins, vector<vector<vector<TH3F*>>>(nMtBins, vector<vector<TH3F*>>(nQnBins, vector<TH3F*>(nPhiBins,nullptr))));

	histoME.clear();
	histoME.resize(nCentBins, vector<vector<TH3F*>>(nMtBins, vector<TH3F*>(nPhiBins, nullptr)));
							
	double lam[nCentBins][nMtBins][nQnBins][nPhiBins];
	double Ro2[nCentBins][nMtBins][nQnBins][nPhiBins];
	double Rs2[nCentBins][nMtBins][nQnBins][nPhiBins];
	double Rl2[nCentBins][nMtBins][nQnBins][nPhiBins];
	double n0[nCentBins][nMtBins][nQnBins][nPhiBins];
	double no[nCentBins][nMtBins][nQnBins][nPhiBins];
	double ns[nCentBins][nMtBins][nQnBins][nPhiBins];
	double nl[nCentBins][nMtBins][nQnBins][nPhiBins];
	double nos[nCentBins][nMtBins][nQnBins][nPhiBins];
	double nol[nCentBins][nMtBins][nQnBins][nPhiBins];
	double nsl[nCentBins][nMtBins][nQnBins][nPhiBins];

	double lamerr[nCentBins][nMtBins][nQnBins][nPhiBins];
	double Ro2err[nCentBins][nMtBins][nQnBins][nPhiBins];
	double Rs2err[nCentBins][nMtBins][nQnBins][nPhiBins];
	double Rl2err[nCentBins][nMtBins][nQnBins][nPhiBins];
	double n0err[nCentBins][nMtBins][nQnBins][nPhiBins];
	double noerr[nCentBins][nMtBins][nQnBins][nPhiBins];
	double nserr[nCentBins][nMtBins][nQnBins][nPhiBins];
	double nlerr[nCentBins][nMtBins][nQnBins][nPhiBins];
	double noserr[nCentBins][nMtBins][nQnBins][nPhiBins];
	double nolerr[nCentBins][nMtBins][nQnBins][nPhiBins];
	double nslerr[nCentBins][nMtBins][nQnBins][nPhiBins];


	for (int iCent=0; iCent<nCentBins; iCent++){
	for (int iMt=0; iMt<nMtBins; iMt++){
	for (int iqn=0; iqn<nQnBins; iqn++){				
	for (int iPhi=0; iPhi<nPhiBins; iPhi++){					
			histoSE[iCent][iMt][iqn][iPhi] = (TH3F*)file->Get(Form("H3_same_cent%d_mt%d_qn%d_phi%d", iCent, iMt, iqn, iPhi));
		}
	}}}


	for (int iCent=0; iCent<nCentBins; iCent++){
	for (int iMt=0; iMt<nMtBins; iMt++){
	for (int iPhi=0; iPhi<nPhiBins; iPhi++){					
			histoME[iCent][iMt][iPhi] = (TH3F*)file->Get(Form("H3_mixed_cent%d_mt%d_phi%d", iCent, iMt, iPhi));
		}
	}}


	for (int iCent=0; iCent<nCentBins; iCent++){
	for (int iMt=0; iMt<nMtBins; iMt++){
	for (int iqn=0; iqn<nQnBins; iqn++){				
	for (int iPhi=0; iPhi<nPhiBins; iPhi++){					
			Fit3D_BowlerSinyukov_Full((int)centBin[iCent], (double)mTBin[iMt], histoSE[iCent][iMt][iqn][iPhi], histoME[iCent][iMt][iPhi], iqn, iPhi,
									  lam[iCent][iMt][iqn][iPhi], Ro2[iCent][iMt][iqn][iPhi], Rs2[iCent][iMt][iqn][iPhi], Rl2[iCent][iMt][iqn][iPhi], 
									  lamerr[iCent][iMt][iqn][iPhi], Ro2err[iCent][iMt][iqn][iPhi], Rs2err[iCent][iMt][iqn][iPhi], Rl2err[iCent][iMt][iqn][iPhi],
									  // n0[iCent][iMt][iqn][iPhi], 
									  // no[iCent][iMt][iqn][iPhi], ns[iCent][iMt][iqn][iPhi], nl[iCent][iMt][iqn][iPhi],
									  // nos[iCent][iMt][iqn][iPhi], nol[iCent][iMt][iqn][iPhi], nsl[iCent][iMt][iqn][iPhi],
									  // n0err[iCent][iMt][iqn][iPhi], 
									  // noerr[iCent][iMt][iqn][iPhi], nserr[iCent][iMt][iqn][iPhi], nlerr[iCent][iMt][iqn][iPhi],
									  // noserr[iCent][iMt][iqn][iPhi], nolerr[iCent][iMt][iqn][iPhi], nslerr[iCent][iMt][iqn][iPhi],
									  mod);		}
	}}}

	// // -----------------------------------------------------
	// // read and write Radii
	// // -----------------------------------------------------
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
		
	// TFile* fout = new TFile("output/radii_mT0.3to0.5_cent30to50_polbaseline_OO_firstlook.root","RECREATE");
	// fout->cd();
	// for (int iCent=0; iCent<nCentBins; iCent++){
	// for (int iMt=0; iMt<nMtBins; iMt++){													
	// for (int iqn=0; iqn<nQnBins; iqn++){					
	// 	histoRoutVsPhi[iCent][iMt][iqn]->Write();
	// 	histoRsideVsPhi[iCent][iMt][iqn]->Write();
	// 	histoRlongVsPhi[iCent][iMt][iqn]->Write();
	// }}}	

	// // -----------------------------------------------------
	// // write lambda
	// // -----------------------------------------------------
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

	// TFile* foutla = new TFile("output/lambda_mT0.3to0.5_cent30to50_polbaseline_tocheck.root","RECREATE");
	// foutla->cd();
	// for (int iCent=0; iCent<nCentBins; iCent++){
	// for (int iMt=0; iMt<nMtBins; iMt++){													
	// for (int iqn=0; iqn<nQnBins; iqn++){					
	// 	histolamVsPhi[iCent][iMt][iqn]->Write();
	// }}}

	// // -----------------------------------------------------
	// // write baseline
	// // -----------------------------------------------------
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


	// histoNoVsPhi.clear();
	// histoNoVsPhi.resize(nCentBins, vector<vector<TH1F*>>(nMtBins, vector<TH1F*>(nQnBins, nullptr)));
	// for (int iCent=0; iCent<nCentBins; iCent++){
	// for (int iMt=0; iMt<nMtBins; iMt++){
	// for (int iqn=0; iqn<nQnBins; iqn++){
	// 	histoNoVsPhi[iCent][iMt][iqn] = new TH1F(Form("no_cent%d_mt%d_qn%d",iCent,iMt,iqn),";#phi^{pair} - #Psi_EP;R_{out}^{2}", 12, -0.2, 3.4);				
	// 	for (int iPhi=0; iPhi<nPhiBins; iPhi++){
	// 		histoNoVsPhi[iCent][iMt][iqn]->SetBinContent(iPhi+1, no[iCent][iMt][iqn][iPhi]);
	// 		histoNoVsPhi[iCent][iMt][iqn]->SetBinError(iPhi+1, noerr[iCent][iMt][iqn][iPhi]);
	// 	}
	// }}}

	// histoNsVsPhi.clear();
	// histoNsVsPhi.resize(nCentBins, vector<vector<TH1F*>>(nMtBins, vector<TH1F*>(nQnBins, nullptr)));
	// for (int iCent=0; iCent<nCentBins; iCent++){
	// for (int iMt=0; iMt<nMtBins; iMt++){
	// for (int iqn=0; iqn<nQnBins; iqn++){
	// 	histoNsVsPhi[iCent][iMt][iqn] = new TH1F(Form("ns_cent%d_mt%d_qn%d",iCent,iMt,iqn),";#phi^{pair} - #Psi_EP;R_{out}^{2}", 12, -0.2, 3.4);				
	// 	for (int iPhi=0; iPhi<nPhiBins; iPhi++){
	// 		histoNsVsPhi[iCent][iMt][iqn]->SetBinContent(iPhi+1, ns[iCent][iMt][iqn][iPhi]);
	// 		histoNsVsPhi[iCent][iMt][iqn]->SetBinError(iPhi+1, nserr[iCent][iMt][iqn][iPhi]);
	// 	}
	// }}}


	// histoNlVsPhi.clear();
	// histoNlVsPhi.resize(nCentBins, vector<vector<TH1F*>>(nMtBins, vector<TH1F*>(nQnBins, nullptr)));
	// for (int iCent=0; iCent<nCentBins; iCent++){
	// for (int iMt=0; iMt<nMtBins; iMt++){
	// for (int iqn=0; iqn<nQnBins; iqn++){
	// 	histoNlVsPhi[iCent][iMt][iqn] = new TH1F(Form("nl_cent%d_mt%d_qn%d",iCent,iMt,iqn),";#phi^{pair} - #Psi_EP;R_{out}^{2}", 12, -0.2, 3.4);				
	// 	for (int iPhi=0; iPhi<nPhiBins; iPhi++){
	// 		histoNlVsPhi[iCent][iMt][iqn]->SetBinContent(iPhi+1, nl[iCent][iMt][iqn][iPhi]);
	// 		histoNlVsPhi[iCent][iMt][iqn]->SetBinError(iPhi+1, nlerr[iCent][iMt][iqn][iPhi]);
	// 	}
	// }}}

	// histoNosVsPhi.clear();
	// histoNosVsPhi.resize(nCentBins, vector<vector<TH1F*>>(nMtBins, vector<TH1F*>(nQnBins, nullptr)));
	// for (int iCent=0; iCent<nCentBins; iCent++){
	// for (int iMt=0; iMt<nMtBins; iMt++){
	// for (int iqn=0; iqn<nQnBins; iqn++){
	// 	histoNosVsPhi[iCent][iMt][iqn] = new TH1F(Form("nos_cent%d_mt%d_qn%d",iCent,iMt,iqn),";#phi^{pair} - #Psi_EP;R_{out}^{2}", 12, -0.2, 3.4);				
	// 	for (int iPhi=0; iPhi<nPhiBins; iPhi++){
	// 		histoNosVsPhi[iCent][iMt][iqn]->SetBinContent(iPhi+1, nos[iCent][iMt][iqn][iPhi]);
	// 		histoNosVsPhi[iCent][iMt][iqn]->SetBinError(iPhi+1, noserr[iCent][iMt][iqn][iPhi]);
	// 	}
	// }}}

	// histoNolVsPhi.clear();
	// histoNolVsPhi.resize(nCentBins, vector<vector<TH1F*>>(nMtBins, vector<TH1F*>(nQnBins, nullptr)));
	// for (int iCent=0; iCent<nCentBins; iCent++){
	// for (int iMt=0; iMt<nMtBins; iMt++){
	// for (int iqn=0; iqn<nQnBins; iqn++){
	// 	histoNolVsPhi[iCent][iMt][iqn] = new TH1F(Form("nol_cent%d_mt%d_qn%d",iCent,iMt,iqn),";#phi^{pair} - #Psi_EP;R_{out}^{2}", 12, -0.2, 3.4);				
	// 	for (int iPhi=0; iPhi<nPhiBins; iPhi++){
	// 		histoNolVsPhi[iCent][iMt][iqn]->SetBinContent(iPhi+1, nol[iCent][iMt][iqn][iPhi]);
	// 		histoNolVsPhi[iCent][iMt][iqn]->SetBinError(iPhi+1, nolerr[iCent][iMt][iqn][iPhi]);
	// 	}
	// }}}

	// histoNslVsPhi.clear();
	// histoNslVsPhi.resize(nCentBins, vector<vector<TH1F*>>(nMtBins, vector<TH1F*>(nQnBins, nullptr)));
	// for (int iCent=0; iCent<nCentBins; iCent++){
	// for (int iMt=0; iMt<nMtBins; iMt++){
	// for (int iqn=0; iqn<nQnBins; iqn++){
	// 	histoNslVsPhi[iCent][iMt][iqn] = new TH1F(Form("nsl_cent%d_mt%d_qn%d",iCent,iMt,iqn),";#phi^{pair} - #Psi_EP;R_{out}^{2}", 12, -0.2, 3.4);				
	// 	for (int iPhi=0; iPhi<nPhiBins; iPhi++){
	// 		histoNslVsPhi[iCent][iMt][iqn]->SetBinContent(iPhi+1, nsl[iCent][iMt][iqn][iPhi]);
	// 		histoNslVsPhi[iCent][iMt][iqn]->SetBinError(iPhi+1, nslerr[iCent][iMt][iqn][iPhi]);
	// 	}
	// }}}

	// // -----------------------------------------------------
	// // write
	// // -----------------------------------------------------
	// TFile* foutn = new TFile("output/baselinefit_cent_mt_qn_tocheck.root","RECREATE");
	// foutn->cd();
	// for (int iCent=0; iCent<nCentBins; iCent++){
	// for (int iMt=0; iMt<nMtBins; iMt++){													
	// for (int iqn=0; iqn<nQnBins; iqn++){					
	// 	histoN0VsPhi[iCent][iMt][iqn]->Write();
	// 	histoNoVsPhi[iCent][iMt][iqn]->Write();
	// 	histoNsVsPhi[iCent][iMt][iqn]->Write();
	// 	histoNlVsPhi[iCent][iMt][iqn]->Write();
	// 	histoNosVsPhi[iCent][iMt][iqn]->Write();
	// 	histoNolVsPhi[iCent][iMt][iqn]->Write();
	// 	histoNslVsPhi[iCent][iMt][iqn]->Write();

	// }}}	
}