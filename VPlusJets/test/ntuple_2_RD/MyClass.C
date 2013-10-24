#define MyClass_cxx
#include "MyClass.h"
#include <TH2.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>
#include "TMath.h"
#include "iostream"
using namespace std;

bool equal(double a1, double a2, double delta=1e-3){
	if(TMath::Abs(a1-a2)<delta)return 1;
	else return 0;
}
bool match_dR(double a1, double a2, double b1, double b2, double delta=0.3){
	double dR= TMath::Sqrt( (a1-b1)*(a1-b1) + (a2-b2)*(a2-b2) ) ;
	if( dR<delta )return 1;
	else return 0;
}

void MyClass::Draw_and_Save(TH1D h1){
	h1.Write();
	TCanvas *c1 = new TCanvas(Form("c1_%s",h1.GetTitle()),Form("c1_%s",h1.GetTitle()),200,10,600,600);
	c1->cd();
	h1.Draw();
	c1->Print(Form("%s/%s_%s_%s.png",plot_Dir_DateTime.Data(), JetType.Data(), PfType.Data(), h1.GetTitle()));
	delete c1;
}
void MyClass::Draw_and_Save(TH1D h1, TH1D h2){
	TCanvas *c1 = new TCanvas(Form("c1_%s_and_%s",h1.GetTitle(),h2.GetTitle()),Form("c1_%s_%s",h1.GetTitle(),h2.GetTitle()),200,10,600,600);
	c1->cd();
	h1.GetYaxis()->SetRangeUser(0., TMath::Max(h1.GetMaximum(), h2.GetMaximum())*1.2);
	h1.Draw();
	h2.Draw("same");
	c1->Print(Form("%s/%s_%s_%s_and_%s.png",plot_Dir_DateTime.Data(), JetType.Data(), PfType.Data(),  h1.GetTitle(),h2.GetTitle()));
	delete c1;
}
void MyClass::Draw_and_Save(TH1D h1, TH1D h2, TH1D h3){
	TCanvas *c1 = new TCanvas(Form("c1_%s_and_%s_and_%s",h1.GetTitle(),h2.GetTitle(),h3.GetTitle()),Form("c1_%s_and_%s_and_%s",h1.GetTitle(),h2.GetTitle(),h3.GetTitle()),200,10,600,600);
	c1->cd();
	h1.GetYaxis()->SetRangeUser(0., TMath::Max(h1.GetMaximum(), h3.GetMaximum())*1.2);
	h1.Draw();
	h2.Draw("same");
	h3.Draw("same");
	c1->Print(Form("%s/%s_%s_%s_and_%s_and_%s.png",plot_Dir_DateTime.Data(), JetType.Data(), PfType.Data(),h1.GetTitle(),h2.GetTitle(),h3.GetTitle()));
	delete c1;
}
void MyClass::Draw_and_Save(TH1D h1, TH1D h2, TH1D h3, TH1D h4){
	TCanvas *c1 = new TCanvas(Form("c1_%s_and_%s_and_%s_and_%s",h1.GetTitle(),h2.GetTitle(),h3.GetTitle(),h4.GetTitle()), Form("c1_%s_and_%s_and_%s_and_%s",h1.GetTitle(),h2.GetTitle(),h3.GetTitle(),h4.GetTitle()),200,10,600,600);
	c1->cd();
	h1.GetYaxis()->SetRangeUser(0., TMath::Max(h1.GetMaximum(), h4.GetMaximum())*1.2);
	h1.Draw();
	h2.Draw("same");
	h3.Draw("same");
	h4.Draw("same");
	c1->Print(Form("%s/%s_%s_%s_and_%s_and_%s_and_%s.png",plot_Dir_DateTime.Data(), JetType.Data(), PfType.Data(),h1.GetTitle(),h2.GetTitle(),h3.GetTitle(),h4.GetTitle()));
	delete c1;
}
void MyClass::Draw_and_Save(TH2D h2, char* addtional_info){
	h2.Write();
	TCanvas *c1 = new TCanvas(Form("c1_%s",h2.GetTitle()),Form("c1_%s",h2.GetTitle()),200,10,600,600);
	c1->cd();
	h2.Draw("box");
	if( addtional_info){
		TLatex tl;
		tl.SetTextSize(0.04 ); tl.SetTextAlign(13);
		tl.DrawLatex(h2.GetXaxis()->GetXmin()*0.9+h2.GetXaxis()->GetXmax()*0.1,h2.GetYaxis()->GetXmin()*0.1+h2.GetYaxis()->GetXmax()*0.9,addtional_info);
	};
	c1->Print(Form("%s/%s_%s_%s.png",plot_Dir_DateTime.Data(), JetType.Data(), PfType.Data(),h2.GetTitle()));
	delete c1;
}

Bool_t MyClass::Select()
{

	if(!( Z_mass>65 && Z_mass<105 ))return 0;
	if(!( GenGroomedJet_pt[0]>20 && GenGroomedJet_pt[0]<5000 ))return 0;

	//alpha cut: z and jet are back to back
	if(GenGroomedJet_number_jet_central >1){
		//cout<<"alpha="<<GenGroomedJet_pt[1]/GenGroomedJet_pt[0]<<endl; 
		if( GenGroomedJet_pt[1]/GenGroomedJet_pt[0] >=0.3 ) return 0;
	}

	if(isBoosted){
		if(!( Z_pt>100))return 0;
		if(!( GenGroomedJet_pt[0]>100 ))return 0;
		//Double_t tmp_GEN_eta = GenGroomedJet_eta[0]; Double_t tmp_GEN_phi = GenGroomedJet_phi[0];
		//Double_t tmpDeltaR_Vj = TMath::Sqrt( (Z_eta-tmp_GEN_eta)*(Z_eta-tmp_GEN_eta) + (Z_phi-tmp_GEN_phi)*(Z_phi-tmp_GEN_phi) );
		//if(!( tmpDeltaR_Vj>2.0 ))return 0;
	}
	return 1;
}

void MyClass::Loop()
{
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();


	Double_t tmp_GEN_pt=0.;
	Double_t tmp_GEN_eta=0.;
	Double_t tmp_GEN_phi=0.;
	Double_t tmp_GEN_rhoSW=0.;
	Double_t tmp_GEN_rhoHand=0.;
	Double_t tmp_GEN_rhoHand2=0.;
	Double_t tmp_GEN_rhoGrid=0.;
	Double_t tmp_GEN_mass=0.;

	Double_t tmp_PFCor_pt=0.;
	Double_t tmp_PFCor_Pt_uncorr=0.;
	Double_t tmp_PFCor_Pt_afterL1=0.;
	Double_t tmp_PFCor_Pt_afterL2=0.;
	Double_t tmp_PFCor_mass=0.;

	Double_t tmp_RECO_eta=0.;
	Double_t tmp_RECO_phi=0.;
	Double_t tmp_RECO_pt_uncorr=0.;
	Double_t tmp_RECO_pt=0.; //after SW JEC
	Double_t tmp_RECO_pt_L1_rhoSW=0.;
	Double_t tmp_RECO_pt_L1_rhoHand=0.;
	Double_t tmp_RECO_pt_L1_rhoHand2=0.;
	Double_t tmp_RECO_pt_L1_rhoGrid=0.;
	Double_t tmp_RECO_pt_rho4A=0.;
	Double_t tmp_RECO_pt_rhom4A=0.;
	//Double_t tmp_RECO_pt_JetCleansing=0.;
	Double_t tmp_RECO_rhoSW=0.;
	Double_t tmp_RECO_rhoHand=0.;
	Double_t tmp_RECO_rhoHand2=0.;
	Double_t tmp_RECO_rhoGrid=0.;
	Double_t tmp_RECO_mass_uncorr=0.;
	Double_t tmp_RECO_mass_jec=0.;
	Double_t tmp_RECO_mass_rhoArea=0.;
	Double_t tmp_RECO_mass_rhoGArea=0.;
	Double_t tmp_RECO_mass_rho4Area=0.;
	Double_t tmp_RECO_mass_rhoG4Area=0.;
	Double_t tmp_RECO_mass_rhom4Area=0.;
	/*Double_t tmp_RECO_mass_JetCleansingATLASjvf=0.;
	  Double_t tmp_RECO_mass_JetCleansingATLASlin=0.;
	  Double_t tmp_RECO_mass_JetCleansingATLASgau=0.;
	  Double_t tmp_RECO_mass_JetCleansingCMSjvf=0.;
	  Double_t tmp_RECO_mass_JetCleansingCMSlin=0.;
	  Double_t tmp_RECO_mass_JetCleansingCMSgau=0.;*/
	Double_t tmp_RECO_tau2tau1=0.;
	Double_t tmp_RECO_tau2tau1_shapesubtract=0.;
	Double_t tmp_GEN_tau2tau1=0.;
	Double_t tmp_GEN_tau2tau1_shapesubtract=0.;


	Double_t tmp_event_nPV=0.;

	Double_t ratio=0.;// reco jet pt / gen jet pt
	Double_t dr=0.; // Delta R (Z, j)
	Double_t dphi=0.; // Delta Phi (Z, j) 


	Int_t nbin_rho=50; Double_t rhomin=0.; Double_t rhomax=50.;
	Int_t nbin_nPV=50; Double_t nPVmin=0.; Double_t nPVmax=50.;
	Int_t nbin_mass=60;Double_t jetmass_min=0;Double_t jetmass_max=300.;
	Int_t nbin_pt=40;Double_t jetpt_min=50;Double_t jetpt_max=450.;
	Int_t nbin_ratio=300; Double_t ratio_min=0; Double_t ratio_max=3.; 
	Double_t ratio_mrt_min=0.5; Double_t ratio_mrt_max=1.6; 
	Double_t ratio_mrt_uncorr_min=0.5; Double_t ratio_mrt_uncorr_max=1.6; 
	const Int_t nbin_eta=10;Double_t eta_min=-2.5;Double_t eta_max=2.5;
	//const Int_t nbin_tau2tau1=40;Double_t tau2tau1_min=-0.5;Double_t tau2tau1_max=1.5;
	const Int_t nbin_tau2tau1=40;Double_t tau2tau1_min=0.;Double_t tau2tau1_max=1.;

	if(!isBoosted){
		nbin_mass=40; jetmass_min=0; jetmass_max=80.;
		nbin_pt=30; jetpt_min=0; jetpt_max=300.;
	}

	TH1D h1_nPV("h1_nPV","h1_nPV;nPV", nbin_nPV, nPVmin, nPVmax);
	TH1D h1_z_mass("h1_z_mass","h1_z_mass;Z mass",60,60,120);
	TH1D h1_z_Pt("h1_z_Pt","h1_z_Pt;Z Pt",40,0,400);
	TH1D h1_muplus_Pt("h1_muplus_Pt","h1_muplus_Pt;#mu^{+} Pt",40,0,400);

	TH1D h1_GEN_eta("h1_GEN_eta","h1_GEN_eta; Gen Jet #eta", nbin_eta, eta_min, eta_max); h1_GEN_eta.SetLineColor(kRed);
	TH1D h1_GEN_phi("h1_GEN_phi","h1_GEN_phi; Gen Jet #phi",50,-4,4); h1_GEN_phi.SetLineColor(kRed);
	TH1D h1_GEN_zjet_dr("h1_GEN_zjet_dr","h1_GEN_zjet_dr;dR(Z,j)",50,0,10); h1_GEN_zjet_dr.SetLineColor(kRed);
	TH1D h1_GEN_zjet_dphi("h1_GEN_zjet_dphi","h1_GEN_zjet_dphi; d#phi(Z,j)",50,0,5); h1_GEN_zjet_dphi.SetLineColor(kRed);
	TH1D h1_GEN_rhoSW("h1_GEN_rhoSW","h1_GEN_rhoSW",nbin_rho,rhomin,rhomax); h1_GEN_rhoSW.SetLineColor(kRed);
	TH1D h1_GEN_rhoHand("h1_GEN_rhoHand","h1_GEN_rhoHand",nbin_rho,rhomin,rhomax); h1_GEN_rhoHand.SetLineColor(kRed);
	TH1D h1_GEN_rhoHand2("h1_GEN_rhoHand2","h1_GEN_rhoHand2",nbin_rho,rhomin,rhomax); h1_GEN_rhoHand2.SetLineColor(kRed);
	TH1D h1_GEN_rhoGrid("h1_GEN_rhoGrid","h1_GEN_rhoGrid",nbin_rho,rhomin,rhomax); h1_GEN_rhoGrid.SetLineColor(kRed);
	TH2D h2_GEN_rhoSW_vs_nPV("h2_GEN_rhoSW_vs_nPV","h2_GEN_rhoSW_vs_nPV",nbin_rho,rhomin,rhomax, nbin_nPV, nPVmin, nPVmax); h2_GEN_rhoSW_vs_nPV.SetMarkerColor(kRed);
	TH2D h2_GEN_rhoHand_vs_nPV("h2_GEN_rhoHand_vs_nPV","h2_GEN_rhoHand_vs_nPV",nbin_rho,rhomin,rhomax, nbin_nPV, nPVmin, nPVmax); h2_GEN_rhoHand_vs_nPV.SetMarkerColor(kRed);
	TH2D h2_GEN_rhoHand2_vs_nPV("h2_GEN_rhoHand2_vs_nPV","h2_GEN_rhoHand2_vs_nPV",nbin_rho,rhomin,rhomax, nbin_nPV, nPVmin, nPVmax); h2_GEN_rhoHand2_vs_nPV.SetMarkerColor(kRed);
	TH2D h2_GEN_rhoGrid_vs_nPV("h2_GEN_rhoGrid_vs_nPV","h2_GEN_rhoGrid_vs_nPV",nbin_rho,rhomin,rhomax, nbin_nPV, nPVmin, nPVmax); h2_GEN_rhoGrid_vs_nPV.SetMarkerColor(kRed);

	//PFCor: ak5PFCHS from SW
	TH1D h1_PFCor_Pt_afterL1("h1_PFCor_Pt_afterL1","h1_PFCor_Pt_afterL1;Jet Pt",nbin_pt, jetpt_min, jetpt_max);
	TH1D h1_PFCor_area("h1_PFCor_area","h1_PFCor_area",50,0.6,1.1);

	// RECO: clusted in the fly

	/*
	   TH1D h1_RECO_pt_uncorr("h1_RECO_pt_uncorr","h1_RECO_pt_uncorr;Jet Pt",nbin_pt, jetpt_min, jetpt_max);
	   TH1D h1_RECO_pt_rhoArea("h1_RECO_pt_rhoArea","h1_RECO_pt_rhoArea",nbin_pt, jetpt_min, jetpt_max);// rhoHand * area
	   TH1D h1_RECO_pt_rhoGArea("h1_RECO_pt_rhoGArea","h1_RECO_pt_rhoGArea",nbin_pt, jetpt_min, jetpt_max);
	   TH1D h1_RECO_pt_rho4Area("h1_RECO_pt_rho4Area","h1_RECO_pt_rho4Area",nbin_pt, jetpt_min, jetpt_max);
	   TH1D h1_RECO_pt_rhom4Area("h1_RECO_pt_rhom4Area","h1_RECO_pt_rhom4Area",nbin_pt, jetpt_min, jetpt_max);

	// reco pt VS gen pt
	TH2D h2_GEN_pt("h2_GEN_pt","h2_GEN_pt; gen jet pt; reco jet pt",nbin_pt,jetpt_min,jetpt_max, nbin_pt,jetpt_min,jetpt_max); 
	TH2D h2_RECO_pt_uncorr("h2_RECO_pt_uncorr","h2_RECO_pt_uncorr; gen jet pt; reco jet pt",nbin_pt,jetpt_min,jetpt_max, nbin_pt,jetpt_min,jetpt_max);
	TH2D h2_RECO_pt_rhoArea("h2_RECO_pt_rhoArea","h2_RECO_pt_rhoArea; gen jet pt; reco jet pt",nbin_pt,jetpt_min,jetpt_max, nbin_pt,jetpt_min,jetpt_max);
	TH2D h2_RECO_pt_rhoGArea("h2_RECO_pt_rhoGArea","h2_RECO_pt_rhoGArea; gen jet pt; reco jet pt",nbin_pt,jetpt_min,jetpt_max, nbin_pt,jetpt_min,jetpt_max);
	TH2D h2_RECO_pt_rho4Area("h2_RECO_pt_rho4Area","h2_RECO_pt_rho4Area; gen jet pt; reco jet pt",nbin_pt,jetpt_min,jetpt_max, nbin_pt,jetpt_min,jetpt_max);
	TH2D h2_RECO_pt_rhom4Area("h2_RECO_pt_rhom4Area","h2_RECO_pt_rhom4Area; gen jet pt; reco jet pt",nbin_pt,jetpt_min,jetpt_max, nbin_pt,jetpt_min,jetpt_max);

	//calculate correlationFactor;
	Int_t  num_points=0;
	TGraph gr_GEN_pt;
	TGraph gr_RECO_pt_uncorr;
	TGraph gr_RECO_pt_rhoArea;
	TGraph gr_RECO_pt_rhoGArea;
	TGraph gr_RECO_pt_rho4Area;
	TGraph gr_RECO_pt_rhom4Area;
	*/
	// reco pt VS gen pt
	correlation_tool ct_GEN_pt("ct_GEN_pt",nbin_pt,jetpt_min,jetpt_max, nbin_pt,jetpt_min,jetpt_max); 
	correlation_tool ct_RECO_pt_uncorr("ct_RECO_pt_uncorr",nbin_pt,jetpt_min,jetpt_max, nbin_pt,jetpt_min,jetpt_max);
	correlation_tool ct_RECO_pt_rhoArea("ct_RECO_pt_rhoArea",nbin_pt,jetpt_min,jetpt_max, nbin_pt,jetpt_min,jetpt_max);
	correlation_tool ct_RECO_pt_rhoGArea("ct_RECO_pt_rhoGArea",nbin_pt,jetpt_min,jetpt_max, nbin_pt,jetpt_min,jetpt_max);
	correlation_tool ct_RECO_pt_rho4Area("ct_RECO_pt_rho4Area",nbin_pt,jetpt_min,jetpt_max, nbin_pt,jetpt_min,jetpt_max);
	correlation_tool ct_RECO_pt_rhom4Area("ct_RECO_pt_rhom4Area",nbin_pt,jetpt_min,jetpt_max, nbin_pt,jetpt_min,jetpt_max);


	TH1D h1_RECO_eta("h1_RECO_eta","h1_RECO_eta", nbin_eta, eta_min, eta_max);
	TH1D h1_RECO_phi("h1_RECO_phi","h1_RECO_phi",50,-4,4);
	TH1D h1_RECO_zjet_dr("h1_RECO_zjet_dr","h1_RECO_zjet_dr;dR(Z,J)",50,0,7);
	TH1D h1_RECO_zjet_dphi("h1_RECO_zjet_dphi","h1_RECO_zjet_dphi",50,0,5);
	TH1D h1_RECO_rhoSW("h1_RECO_rhoSW","h1_RECO_rhoSW",nbin_rho,rhomin,rhomax);
	TH1D h1_RECO_rhoHand("h1_RECO_rhoHand","h1_RECO_rhoHand",nbin_rho,rhomin,rhomax);
	TH1D h1_RECO_rhoHand2("h1_RECO_rhoHand2","h1_RECO_rhoHand2",nbin_rho,rhomin,rhomax);
	TH1D h1_RECO_rhoGrid("h1_RECO_rhoGrid","h1_RECO_rhoGrid",nbin_rho,rhomin,rhomax);

	TH2D h2_RECO_rhoSW_vs_nPV("h2_RECO_rhoSW_vs_nPV","h2_RECO_rhoSW_vs_nPV",nbin_rho,rhomin,rhomax, nbin_nPV, nPVmin, nPVmax);
	TH2D h2_RECO_rhoHand_vs_nPV("h2_RECO_rhoHand_vs_nPV","h2_RECO_rhoHand_vs_nPV",nbin_rho,rhomin,rhomax, nbin_nPV, nPVmin, nPVmax);
	TH2D h2_RECO_rhoHand2_vs_nPV("h2_RECO_rhoHand2_vs_nPV","h2_RECO_rhoHand2_vs_nPV",nbin_rho,rhomin,rhomax, nbin_nPV, nPVmin, nPVmax);
	TH2D h2_RECO_rhoGrid_vs_nPV("h2_RECO_rhoGrid_vs_nPV","h2_RECO_rhoGrid_vs_nPV",nbin_rho,rhomin,rhomax, nbin_nPV, nPVmin, nPVmax);







	//tau2tau1

	//TH1D h1_RECO_tau2tau1("h1_RECO_tau2tau1","h1_RECO_tau2tau1", nbin_tau2tau1, tau2tau1_min, tau2tau1_max);
	//TH1D h1_RECO_tau2tau1_shapesubtract("h1_RECO_tau2tau1_shapesubtract","h1_RECO_tau2tau1_shapesubtract", nbin_tau2tau1, tau2tau1_min, tau2tau1_max);

	correlation_tool ct_RECO_tau2tau1("ct_RECO_tau2tau1", nbin_tau2tau1, tau2tau1_min, tau2tau1_max, nbin_tau2tau1, tau2tau1_min, tau2tau1_max);
	correlation_tool ct_RECO_tau2tau1_shapesubtract("ct_RECO_tau2tau1_shapesubtract", nbin_tau2tau1, tau2tau1_min, tau2tau1_max, nbin_tau2tau1, tau2tau1_min, tau2tau1_max);

	//TH1D h1_GEN_tau2tau1("h1_GEN_tau2tau1","h1_GEN_tau2tau1", nbin_tau2tau1, tau2tau1_min, tau2tau1_max);
	//TH1D h1_GEN_tau2tau1_shapesubtract("h1_GEN_tau2tau1_shapesubtract","h1_GEN_tau2tau1_shapesubtract", nbin_tau2tau1, tau2tau1_min, tau2tau1_max);

	correlation_tool ct_GEN_tau2tau1("ct_GEN_tau2tau1", nbin_tau2tau1, tau2tau1_min, tau2tau1_max, nbin_tau2tau1, tau2tau1_min, tau2tau1_max);
	correlation_tool ct_GEN_tau2tau1_shapesubtract("ct_GEN_tau2tau1_shapesubtract", nbin_tau2tau1, tau2tau1_min, tau2tau1_max, nbin_tau2tau1, tau2tau1_min, tau2tau1_max);

	TH1D h1_RECO_area("h1_RECO_area","h1_RECO_area",50,0.6,1.1);

	// gen and reco jet deltaR
	TH1D h1_RecoGen_matching("h1_RecoGen_matching","h1_RecoGen_matching; #delta R( reco j, gen j)",50,0,1.); h1_RecoGen_matching.SetLineColor(kRed);//matching with GEN

	// reco Jet Pt / Gen jet Pt
	TH1D h1_PFCor_jec_recogenptratio("h1_PFCor_jec_recogenptratio","h1_PFCor_jec_recogenptratio",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_PFCor_uncorr_recogenptratio("h1_PFCor_uncorr_recogenptratio","h1_PFCor_uncorr_recogenptratio",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_PFCor_afterL1_recogenptratio("h1_PFCor_afterL1_recogenptratio","h1_PFCor_afterL1_recogenptratio",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_PFCor_afterL2_recogenptratio("h1_PFCor_afterL2_recogenptratio","h1_PFCor_afterL2_recogenptratio",nbin_ratio, ratio_min, ratio_max);

	TH1D h1_RECO_jec_recogenptratio("h1_RECO_jec_recogenptratio","h1_RECO_jec_recogenptratio",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_RECO_uncorr_recogenptratio("h1_RECO_uncorr_recogenptratio","h1_RECO_uncorr_recogenptratio;Reco/Gen Pt ratio ",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_RECO_l1rhosw_recogenptratio("h1_RECO_l1rhosw_recogenptratio","h1_RECO_l1rhosw_recogenptratio",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_RECO_l1rhoHand_recogenptratio("h1_RECO_l1rhoHand_recogenptratio","h1_RECO_l1rhoHand_recogenptratio",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_RECO_l1rhoHand2_recogenptratio("h1_RECO_l1rhoHand2_recogenptratio","h1_RECO_l1rhoHand2_recogenptratio",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_RECO_l1rhoGrid_recogenptratio("h1_RECO_l1rhoGrid_recogenptratio","h1_RECO_l1rhoGrid_recogenptratio",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_RECO_rho4A_recogenptratio("h1_RECO_rho4A_recogenptratio","h1_RECO_rho4A_recogenptratio",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_RECO_rhom4A_recogenptratio("h1_RECO_rhom4A_recogenptratio","h1_RECO_rhom4A_recogenptratio",nbin_ratio, ratio_min, ratio_max);
	//TH1D h1_RECO_JetCleansing_recogenptratio("h1_RECO_JetCleansing_recogenptratio","h1_RECO_JetCleansing_recogenptratio",nbin_ratio, ratio_min, ratio_max);

	// reco/gen vs eta
	mean_rms_tool mrt_PFCor_jec_recogenptratio("mrt_PFCor_jec_recogenptratio_vs_eta",nbin_eta, eta_min, eta_max);
	mean_rms_tool mrt_PFCor_uncorr_recogenptratio("mrt_PFCor_uncorr_recogenptratio_vs_eta",nbin_eta, eta_min, eta_max);
	mean_rms_tool mrt_PFCor_afterL1_recogenptratio("mrt_PFCor_afterL1_recogenptratio_vs_eta",nbin_eta, eta_min, eta_max);
	mean_rms_tool mrt_PFCor_afterL2_recogenptratio("mrt_PFCor_afterL2_recogenptratio_vs_eta",nbin_eta, eta_min, eta_max);

	mean_rms_tool mrt_RECO_jec_recogenptratio("mrt_RECO_jec_recogenptratio_vs_eta",nbin_eta, eta_min, eta_max);
	mean_rms_tool mrt_RECO_uncorr_recogenptratio("mrt_RECO_uncorr_recogenptratio_vs_eta",nbin_eta, eta_min, eta_max);
	mean_rms_tool mrt_RECO_l1rhosw_recogenptratio("mrt_RECO_l1rhosw_recogenptratio_vs_eta",nbin_eta, eta_min, eta_max);
	mean_rms_tool mrt_RECO_l1rhoHand_recogenptratio("mrt_RECO_l1rhoHand_recogenptratio_vs_eta",nbin_eta, eta_min, eta_max);
	mean_rms_tool mrt_RECO_l1rhoHand2_recogenptratio("mrt_RECO_l1rhoHand2_recogenptratio_vs_eta",nbin_eta, eta_min, eta_max);
	mean_rms_tool mrt_RECO_l1rhoGrid_recogenptratio("mrt_RECO_l1rhoGrid_recogenptratio_vs_eta",nbin_eta, eta_min, eta_max);
	mean_rms_tool mrt_RECO_rho4A_recogenptratio("mrt_RECO_rho4A_recogenptratio_vs_eta",nbin_eta, eta_min, eta_max);
	mean_rms_tool mrt_RECO_rhom4A_recogenptratio("mrt_RECO_rhom4A_recogenptratio_vs_eta",nbin_eta, eta_min, eta_max);

	TH2D h2_RECO_l1rhoHand_recogenptratio_vs_nPV("h2_RECO_l1rhoHand_recogenptratio_vs_nPV","h2_RECO_l1rhoHand_recogenptratio_vs_nPV", nbin_nPV, nPVmin, nPVmax, 50, 0, 2.5);
	TH2D h2_RECO_l1rhoHand_recogenptratio_vs_ptHand("h2_RECO_l1rhoHand_recogenptratio_vs_ptHand","h2_RECO_l1rhoHand_recogenptratio_vs_ptHand",50,0,200, 50, 0, 2.5);
	TH2D h2_RECO_l1rhoHand_recogenptratio_vs_eta("h2_RECO_l1rhoHand_recogenptratio_vs_eta","h2_RECO_l1rhoHand_recogenptratio_vs_eta", nbin_eta, eta_min, eta_max, nbin_ratio, ratio_min, ratio_max);

	TH1D h1_PFCor_mass("h1_PFCor_mass","h1_PFCor_mass;jet mass;",nbin_mass,jetmass_min,jetmass_max); h1_PFCor_mass.SetLineColor(kBlack);

	TH2D h2_RECO_mass_jec_vs_PV("h1_RECO_mass_jec","h1_RECO_mass_jec; #PV;jet mass;", nbin_nPV, nPVmin, nPVmax, nbin_mass,jetmass_min,jetmass_max);
/*	//mass
	TH1D h1_GEN_mass("h1_GEN_mass","h1_GEN_mass;jet mass;",nbin_mass,jetmass_min,jetmass_max); h1_GEN_mass.SetLineColor(kRed);
	TH1D h1_RECO_mass_uncorr("h1_RECO_mass_uncorr","h1_RECO_mass_uncorr;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_RECO_mass_jec("h1_RECO_mass_jec","h1_RECO_mass_jec;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_RECO_mass_rhoArea("h1_RECO_mass_rhoArea","h1_RECO_mass_rhoArea;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_RECO_mass_rhoGArea("h1_RECO_mass_rhoGArea","h1_RECO_mass_rhoGArea;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_RECO_mass_rho4Area("h1_RECO_mass_rho4Area","h1_RECO_mass_rho4Area;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_RECO_mass_rhoG4Area("h1_RECO_mass_rhoG4Area","h1_RECO_mass_rhoG4Area;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_RECO_mass_rhom4Area("h1_RECO_mass_rhom4Area","h1_RECO_mass_rhom4Area;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	// reco mass VS gen mass
	TH2D h2_GEN_mass("h2_GEN_mass","h2_GEN_mass; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max); 
	TH2D h2_RECO_mass_uncorr("h2_RECO_mass_uncorr","h2_RECO_mass_uncorr; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_RECO_mass_jec("h2_RECO_mass_jec","h2_RECO_mass_jec; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_RECO_mass_rhoArea("h2_RECO_mass_rhoArea","h2_RECO_mass_rhoArea; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_RECO_mass_rhoGArea("h2_RECO_mass_rhoGArea","h2_RECO_mass_rhoGArea; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_RECO_mass_rho4Area("h2_RECO_mass_rho4Area","h2_RECO_mass_rho4Area; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_RECO_mass_rhoG4Area("h2_RECO_mass_rhoG4Area","h2_RECO_mass_rhoG4Area; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_RECO_mass_rhom4Area("h2_RECO_mass_rhom4Area","h2_RECO_mass_rhom4Area; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	//calculate correlationFactor;
	TGraph gr_GEN_mass;
	TGraph gr_RECO_mass_uncorr;
	TGraph gr_RECO_mass_jec;
	TGraph gr_RECO_mass_rhoArea;
	TGraph gr_RECO_mass_rhoGArea;
	TGraph gr_RECO_mass_rho4Area;
	TGraph gr_RECO_mass_rhoG4Area;
	TGraph gr_RECO_mass_rhom4Area;
*/

	// reco mass VS gen mass
	correlation_tool ct_GEN_mass("ct_GEN_mass","ct_GEN_mass; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max); 
	correlation_tool ct_RECO_mass_uncorr("ct_RECO_mass_uncorr","ct_RECO_mass_uncorr; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	correlation_tool ct_RECO_mass_jec("ct_RECO_mass_jec","ct_RECO_mass_jec; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	correlation_tool ct_RECO_mass_rhoArea("ct_RECO_mass_rhoArea","ct_RECO_mass_rhoArea; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	correlation_tool ct_RECO_mass_rhoGArea("ct_RECO_mass_rhoGArea","ct_RECO_mass_rhoGArea; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	correlation_tool ct_RECO_mass_rho4Area("ct_RECO_mass_rho4Area","ct_RECO_mass_rho4Area; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	correlation_tool ct_RECO_mass_rhoG4Area("ct_RECO_mass_rhoG4Area","ct_RECO_mass_rhoG4Area; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	correlation_tool ct_RECO_mass_rhom4Area("ct_RECO_mass_rhom4Area","ct_RECO_mass_rhom4Area; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);



	//const Double_t NUM_JETCLEANSING_DIFFMODE =200; 
	/*std::vector<TH1D>   vect_h1_RECO_mass_JetCleansing_DiffMode;
	std::vector<TH2D>   vect_h2_RECO_mass_JetCleansing_DiffMode;
	std::vector<TGraph> vect_gr_RECO_mass_JetCleansing_DiffMode;
	std::vector<TH1D>   vect_h1_RECO_pt_JetCleansing_DiffMode;
	std::vector<TH2D>   vect_h2_RECO_pt_JetCleansing_DiffMode;
	std::vector<TGraph> vect_gr_RECO_pt_JetCleansing_DiffMode;*/
	std::vector<mean_rms_tool> vect_mrt_RECO_pt_JetCleansing_DiffMode;
	std::vector<correlation_tool> vect_ct_RECO_mass_JetCleansing_DiffMode;
	//std::vector<correlation_tool> vect_ct_RECO_pt_JetCleansing_DiffMode;
	std::vector<correlation_tool> vect_ct_RECO_tau2tau1_JetCleansing_DiffMode;

	Int_t number_JetCleansing_DiffMode=0;
	for(Int_t i=0;i<NUM_JETCLEANSING_DIFFMODE;i++){
		/*TH1D h1_RECO_mass_JetCleansing_DiffMode(Form("h1_RECO_mass_JetCleansing_DiffMode%i",i),Form("h1_RECO_mass_JetCleansing_DiffMode%i;jet mass;",i),nbin_mass,jetmass_min,jetmass_max);
		vect_h1_RECO_mass_JetCleansing_DiffMode.push_back(h1_RECO_mass_JetCleansing_DiffMode);
		TH2D h2_RECO_mass_JetCleansing_DiffMode(Form("h2_RECO_mass_JetCleansing_DiffMode%i",i),Form("h2_RECO_mass_JetCleansing_DiffMode%i; gen jet mass; reco jet mass",i),nbin_mass,jetmass_min,jetmass_max,nbin_mass,jetmass_min,jetmass_max);
		vect_h2_RECO_mass_JetCleansing_DiffMode.push_back(h2_RECO_mass_JetCleansing_DiffMode);
		TGraph gr_RECO_mass_JetCleansing_DiffMode;
		vect_gr_RECO_mass_JetCleansing_DiffMode.push_back(gr_RECO_mass_JetCleansing_DiffMode);

		TH1D h1_RECO_pt_JetCleansing_DiffMode(Form("h1_RECO_pt_JetCleansing_DiffMode%i",i),Form("h1_RECO_pt_JetCleansing_DiffMode%i;jet pT",i),nbin_pt, jetpt_min, jetpt_max);
		vect_h1_RECO_pt_JetCleansing_DiffMode.push_back(h1_RECO_pt_JetCleansing_DiffMode);
		TH2D h2_RECO_pt_JetCleansing_DiffMode(Form("h2_RECO_pt_JetCleansing_DiffMode%i",i),Form("h2_RECO_pt_JetCleansing_DiffMode%i; gen jet pt; reco jet pt",i),nbin_pt,jetpt_min,jetpt_max,nbin_pt,jetpt_min,jetpt_max);
		vect_h2_RECO_pt_JetCleansing_DiffMode.push_back(h2_RECO_pt_JetCleansing_DiffMode);
		TGraph gr_RECO_pt_JetCleansing_DiffMode;
		vect_gr_RECO_pt_JetCleansing_DiffMode.push_back(gr_RECO_pt_JetCleansing_DiffMode);*/

		mean_rms_tool mrt_RECO_JetCleansing_recogenptratio_vs_eta_DiffMode(Form("mrt_RECO_JetCleansing_recogenptratio_vs_eta_DiffMode%i",i),nbin_eta, eta_min, eta_max, nbin_ratio, ratio_min, ratio_max);
		vect_mrt_RECO_pt_JetCleansing_DiffMode.push_back(mrt_RECO_JetCleansing_recogenptratio_vs_eta_DiffMode);

		correlation_tool ct_RECO_mass_JetCleansing_DiffMode(Form("ct_RECO_mass_JetCleansing_DiffMode%i",i), nbin_mass, jetmass_min, jetmass_max, nbin_mass, jetmass_min, jetmass_max);
		vect_ct_RECO_mass_JetCleansing_DiffMode.push_back(ct_RECO_mass_JetCleansing_DiffMode);

		//correlation_tool ct_RECO_pt_JetCleansing_DiffMode(Form("ct_RECO_pt_JetCleansing_DiffMode%i",i), nbin_pt, jetpt_min, jetpt_max, nbin_pt, jetpt_min, jetpt_max);
		//vect_ct_RECO_pt_JetCleansing_DiffMode.push_back(ct_RECO_pt_JetCleansing_DiffMode);

		correlation_tool ct_RECO_tau2tau1_JetCleansing_DiffMode(Form("ct_RECO_tau2tau1_JetCleansing_DiffMode%i",i), nbin_tau2tau1, tau2tau1_min, tau2tau1_max, nbin_tau2tau1, tau2tau1_min, tau2tau1_max);
		vect_ct_RECO_tau2tau1_JetCleansing_DiffMode.push_back(ct_RECO_tau2tau1_JetCleansing_DiffMode);
	}


	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		//cout<<"jentry="<<jentry<<endl;
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;
		if (!Select())continue;

		tmp_GEN_eta = GenGroomedJet_eta[0];
		tmp_GEN_phi = GenGroomedJet_phi[0];
		tmp_RECO_eta = GroomedJet_eta[0];
		tmp_RECO_phi = GroomedJet_phi[0];

		tmp_RECO_tau2tau1 = GroomedJet_tau2tau1[0];
		tmp_RECO_tau2tau1_shapesubtract = GroomedJet_tau2tau1_shapesubtract[0];
		//cout<<"tau2tau1={ "<<tmp_RECO_tau2tau1<<" , "<<tmp_RECO_tau2tau1_shapesubtract<<" }"<<endl;
		tmp_GEN_tau2tau1 = GenGroomedJet_tau2tau1[0];
		tmp_GEN_tau2tau1_shapesubtract = GenGroomedJet_tau2tau1_shapesubtract[0];
		//cout<<"tau2tau1 Gen={ "<<tmp_GEN_tau2tau1<<" , "<<tmp_GEN_tau2tau1_shapesubtract<<" }"<<endl;

		//RECO-GEN Jet matching efficiency
		dr = TMath::Sqrt( (tmp_GEN_eta-tmp_RECO_eta)*(tmp_GEN_eta-tmp_RECO_eta) +
					(tmp_GEN_phi-tmp_RECO_phi)*(tmp_GEN_phi-tmp_RECO_phi) ) ;
		h1_RecoGen_matching.Fill(dr);

		int i_PFCorJet_matching_PFCHS=0;//_matching PFCor with reco
		Bool_t isPFCor_matching_RECO=0;
		for(;i_PFCorJet_matching_PFCHS<numPFCorJets;i_PFCorJet_matching_PFCHS++){
			Double_t PFCor_jet_eta=JetPFCor_Eta[i_PFCorJet_matching_PFCHS];
			Double_t PFCor_jet_phi=JetPFCor_Phi[i_PFCorJet_matching_PFCHS];
			// phi and eta of PFCor and reco should be almost same, because all from PFCHS
			//if ( equal(tmp_RECO_eta, PFCor_jet_eta) && equal(tmp_RECO_phi, PFCor_jet_phi) ) 
			if ( match_dR(tmp_RECO_eta, tmp_RECO_phi, PFCor_jet_eta, PFCor_jet_phi) )
			{ 
				isPFCor_matching_RECO=1;
				break;
			}
		}
		//cout<<"isPFCor_matching_RECO="<<isPFCor_matching_RECO<<endl;
		//cout<<"i_PFCorJet_matching_PFCHS="<<i_PFCorJet_matching_PFCHS<<endl;

		if( dr <0.3 && isPFCor_matching_RECO) //Gen Jet matching with PF_uncorr
		{
			tmp_GEN_pt = GenGroomedJet_pt[0];

			tmp_PFCor_pt = JetPFCor_Pt[i_PFCorJet_matching_PFCHS];
			tmp_PFCor_Pt_uncorr = JetPFCor_Pt_uncorr[i_PFCorJet_matching_PFCHS];
			tmp_PFCor_Pt_afterL1 = JetPFCor_Pt_afterL1[i_PFCorJet_matching_PFCHS];
			tmp_PFCor_Pt_afterL2 = JetPFCor_Pt_afterL2[i_PFCorJet_matching_PFCHS];

			tmp_RECO_pt = GroomedJet_pt[0];
			tmp_RECO_pt_uncorr = GroomedJet_pt_uncorr[0];
			tmp_RECO_pt_L1_rhoSW = GroomedJet_pt_L1_rhoSW[0];
			tmp_RECO_pt_L1_rhoHand = GroomedJet_pt_L1_rhoHand[0];
			tmp_RECO_pt_L1_rhoHand2 = GroomedJet_pt_L1_rhoHand2[0];
			tmp_RECO_pt_L1_rhoGrid = GroomedJet_pt_L1_rhoGrid[0];

			tmp_RECO_pt_rho4A = GroomedJet_pt_rho4A[0];
			tmp_RECO_pt_rhom4A = GroomedJet_pt_rhom4A[0];
			//tmp_RECO_pt_JetCleansing = GroomedJet_pt_JetCleansing[0];

			tmp_event_nPV = event_nPV;

			tmp_GEN_rhoSW = GenGroomedJet_rhoSW;
			tmp_GEN_rhoHand = GenGroomedJet_rhohand;
			tmp_GEN_rhoHand2 = GenGroomedJet_rhohand2;
			tmp_GEN_rhoGrid = GenGroomedJet_rhogrid;

			tmp_RECO_rhoSW = GroomedJet_rhoSW;
			tmp_RECO_rhoHand = GroomedJet_rhohand;
			tmp_RECO_rhoHand2 = GroomedJet_rhohand2;
			tmp_RECO_rhoGrid = GroomedJet_rhogrid;

			//============= fill hist ==============
			//PFCor
			ratio = tmp_PFCor_pt/tmp_GEN_pt;
			h1_PFCor_jec_recogenptratio.Fill(ratio);
			mrt_PFCor_jec_recogenptratio.Fill(tmp_GEN_eta, ratio);
			//PFCor uncorr
			ratio = tmp_PFCor_Pt_uncorr/tmp_GEN_pt;
			h1_PFCor_uncorr_recogenptratio.Fill(ratio);
			mrt_PFCor_uncorr_recogenptratio.Fill(tmp_GEN_eta, ratio);
			//PFCor L1
			ratio = tmp_PFCor_Pt_afterL1/tmp_GEN_pt;
			h1_PFCor_afterL1_recogenptratio.Fill(ratio);
			mrt_PFCor_afterL1_recogenptratio.Fill(tmp_GEN_eta, ratio);
			h1_PFCor_Pt_afterL1.Fill(tmp_PFCor_Pt_afterL1);
			//PFCor L2
			ratio = tmp_PFCor_Pt_afterL2/tmp_GEN_pt;
			h1_PFCor_afterL2_recogenptratio.Fill(ratio);
			mrt_PFCor_afterL2_recogenptratio.Fill(tmp_GEN_eta, ratio);

			//RECO PF
			ratio = tmp_RECO_pt/tmp_GEN_pt;
			h1_RECO_jec_recogenptratio.Fill(ratio);
			mrt_RECO_jec_recogenptratio.Fill(tmp_GEN_eta, ratio);
			//RECO uncorr
			ratio = tmp_RECO_pt_uncorr/tmp_GEN_pt;
			h1_RECO_uncorr_recogenptratio.Fill(ratio);
			mrt_RECO_uncorr_recogenptratio.Fill(tmp_GEN_eta, ratio);
			//RECO rhoSW
			ratio = tmp_RECO_pt_L1_rhoSW/tmp_GEN_pt;
			h1_RECO_l1rhosw_recogenptratio.Fill(ratio);
			mrt_RECO_l1rhosw_recogenptratio.Fill(tmp_GEN_eta, ratio);
			//RECO rhohand
			ratio = tmp_RECO_pt_L1_rhoHand/tmp_GEN_pt;
			h1_RECO_l1rhoHand_recogenptratio.Fill(ratio);
			mrt_RECO_l1rhoHand_recogenptratio.Fill(tmp_GEN_eta, ratio);

			h2_RECO_l1rhoHand_recogenptratio_vs_nPV.Fill(tmp_event_nPV, ratio);
			h2_RECO_l1rhoHand_recogenptratio_vs_ptHand.Fill(tmp_RECO_pt_L1_rhoHand, ratio);
			h2_RECO_l1rhoHand_recogenptratio_vs_eta.Fill(tmp_RECO_eta, ratio );

			//RECO rhohand2
			ratio = tmp_RECO_pt_L1_rhoHand2/tmp_GEN_pt;
			h1_RECO_l1rhoHand2_recogenptratio.Fill(ratio);
			mrt_RECO_l1rhoHand2_recogenptratio.Fill(tmp_GEN_eta, ratio);
			//RECO rhogrid
			ratio = tmp_RECO_pt_L1_rhoGrid/tmp_GEN_pt;
			h1_RECO_l1rhoGrid_recogenptratio.Fill(ratio);
			mrt_RECO_l1rhoGrid_recogenptratio.Fill(tmp_GEN_eta, ratio);

			ratio = tmp_RECO_pt_rho4A/tmp_GEN_pt;
			h1_RECO_rho4A_recogenptratio.Fill(ratio);
			mrt_RECO_rho4A_recogenptratio.Fill(tmp_GEN_eta, ratio);

			ratio = tmp_RECO_pt_rhom4A/tmp_GEN_pt;
			h1_RECO_rhom4A_recogenptratio.Fill(ratio);
			mrt_RECO_rhom4A_recogenptratio.Fill(tmp_GEN_eta, ratio);

			//ratio = tmp_RECO_pt_JetCleansing/tmp_GEN_pt;
			//h1_RECO_JetCleansing_recogenptratio.Fill(ratio);

			//pt
			/*h1_GEN_pt.Fill(tmp_GEN_pt);
			h1_RECO_pt_uncorr.Fill(tmp_RECO_pt_uncorr);// pt of PF is with wrong JEC now
			h1_RECO_pt_rhoArea.Fill(tmp_RECO_pt_L1_rhoHand);// pt of PF is with wrong JEC now
			h1_RECO_pt_rhoGArea.Fill(tmp_RECO_pt_L1_rhoGrid);// pt of PF is with wrong JEC now
			h1_RECO_pt_rho4Area.Fill(tmp_RECO_pt_rho4A);// 
			h1_RECO_pt_rhom4Area.Fill(tmp_RECO_pt_rhom4A);// 
			//h1_RECO_pt_JetCleansing.Fill(tmp_RECO_pt_JetCleansing);// 

			h2_GEN_pt.Fill(tmp_GEN_pt, tmp_GEN_pt            ); 
			h2_RECO_pt_uncorr.Fill(tmp_GEN_pt, tmp_RECO_pt_uncorr  );
			h2_RECO_pt_rhoArea.Fill(tmp_GEN_pt, tmp_RECO_pt_L1_rhoHand );
			h2_RECO_pt_rhoGArea.Fill(tmp_GEN_pt, tmp_RECO_pt_L1_rhoGrid);
			h2_RECO_pt_rho4Area.Fill(tmp_GEN_pt, tmp_RECO_pt_rho4A);
			h2_RECO_pt_rhom4Area.Fill(tmp_GEN_pt, tmp_RECO_pt_rhom4A);

			gr_GEN_pt.SetPoint(num_points, tmp_GEN_pt, tmp_GEN_pt            ); 
			gr_RECO_pt_uncorr.SetPoint(num_points,tmp_GEN_pt, tmp_RECO_pt_uncorr  );
			gr_RECO_pt_rhoArea.SetPoint(num_points,tmp_GEN_pt, tmp_RECO_pt_L1_rhoHand );
			gr_RECO_pt_rhoGArea.SetPoint(num_points,tmp_GEN_pt, tmp_RECO_pt_L1_rhoGrid);
			gr_RECO_pt_rho4Area.SetPoint(num_points,tmp_GEN_pt, tmp_RECO_pt_rho4A);
			gr_RECO_pt_rhom4Area.SetPoint(num_points,tmp_GEN_pt, tmp_RECO_pt_rhom4A);*/

			ct_GEN_pt.Fill(tmp_GEN_pt, tmp_GEN_pt            ); 
			ct_RECO_pt_uncorr.Fill(tmp_GEN_pt, tmp_RECO_pt_uncorr  );
			ct_RECO_pt_rhoArea.Fill(tmp_GEN_pt, tmp_RECO_pt_L1_rhoHand );
			ct_RECO_pt_rhoGArea.Fill(tmp_GEN_pt, tmp_RECO_pt_L1_rhoGrid);
			ct_RECO_pt_rho4Area.Fill(tmp_GEN_pt, tmp_RECO_pt_rho4A);
			ct_RECO_pt_rhom4Area.Fill(tmp_GEN_pt, tmp_RECO_pt_rhom4A);

			if( GenGroomedJet_mass[0]>40 ){
				//tau2tau1
				//h1_RECO_tau2tau1.Fill(tmp_RECO_tau2tau1);
				//h1_RECO_tau2tau1_shapesubtract.Fill(tmp_RECO_tau2tau1_shapesubtract);
				//h1_GEN_tau2tau1.Fill(tmp_GEN_tau2tau1);
				//h1_GEN_tau2tau1_shapesubtract.Fill(tmp_GEN_tau2tau1_shapesubtract);

				ct_RECO_tau2tau1.Fill(tmp_GEN_tau2tau1, tmp_RECO_tau2tau1);
				ct_RECO_tau2tau1_shapesubtract.Fill(tmp_GEN_tau2tau1, tmp_RECO_tau2tau1_shapesubtract);
				ct_GEN_tau2tau1.Fill(tmp_GEN_tau2tau1, tmp_GEN_tau2tau1);
				ct_GEN_tau2tau1_shapesubtract.Fill(tmp_GEN_tau2tau1, tmp_GEN_tau2tau1_shapesubtract);
			}

			//eta
			h1_GEN_eta.Fill(tmp_GEN_eta);
			h1_RECO_eta.Fill(tmp_RECO_eta);
			//phi
			h1_GEN_phi.Fill(tmp_GEN_phi);
			h1_RECO_phi.Fill(tmp_RECO_phi);
			//dr
			dr = TMath::Sqrt( (Z_eta-tmp_GEN_eta)*(Z_eta-tmp_GEN_eta) + (Z_phi-tmp_GEN_phi)*(Z_phi-tmp_GEN_phi) );
			h1_GEN_zjet_dr.Fill(dr);
			dr = TMath::Sqrt( (Z_eta-tmp_RECO_eta)*(Z_eta-tmp_RECO_eta) + (Z_phi-tmp_RECO_phi)*(Z_phi-tmp_RECO_phi) );
			h1_RECO_zjet_dr.Fill(dr);
			//dphi
			dphi = TMath::Sqrt( (Z_phi-tmp_GEN_phi)*(Z_phi-tmp_GEN_phi) );
			h1_GEN_zjet_dphi.Fill(dphi);
			dphi = TMath::Sqrt( (Z_phi-tmp_RECO_phi)*(Z_phi-tmp_RECO_phi) );
			h1_RECO_zjet_dphi.Fill(dphi);
			//nPV
			h1_nPV.Fill(tmp_event_nPV);
			//z mass
			h1_z_mass.Fill(Z_mass);
			h1_z_Pt.Fill(Z_pt);
			//muplus pt
			h1_muplus_Pt.Fill(Z_muplus_pt);

			//rho
			h1_GEN_rhoSW.Fill(tmp_GEN_rhoSW);
			h1_RECO_rhoSW.Fill(tmp_RECO_rhoSW);

			h1_GEN_rhoHand.Fill(tmp_GEN_rhoHand);
			h1_RECO_rhoHand.Fill(tmp_RECO_rhoHand);

			h1_GEN_rhoHand2.Fill(tmp_GEN_rhoHand2);
			h1_RECO_rhoHand2.Fill(tmp_RECO_rhoHand2);

			h1_GEN_rhoGrid.Fill(tmp_GEN_rhoGrid);
			h1_RECO_rhoGrid.Fill(tmp_RECO_rhoGrid);

			h1_PFCor_area.Fill(JetPFCor_Area[i_PFCorJet_matching_PFCHS]);
			h1_RECO_area.Fill(GroomedJet_area[0]);  

			//rho vs nPV
			h2_GEN_rhoSW_vs_nPV.Fill(tmp_GEN_rhoSW,tmp_event_nPV);
			h2_GEN_rhoHand_vs_nPV.Fill(tmp_GEN_rhoHand,tmp_event_nPV);
			h2_GEN_rhoHand2_vs_nPV.Fill(tmp_GEN_rhoHand2,tmp_event_nPV);
			h2_GEN_rhoGrid_vs_nPV.Fill(tmp_GEN_rhoGrid,tmp_event_nPV);

			h2_RECO_rhoSW_vs_nPV.Fill(tmp_RECO_rhoSW,tmp_event_nPV);
			h2_RECO_rhoHand_vs_nPV.Fill(tmp_RECO_rhoHand,tmp_event_nPV);
			h2_RECO_rhoHand2_vs_nPV.Fill(tmp_RECO_rhoHand2,tmp_event_nPV);
			h2_RECO_rhoGrid_vs_nPV.Fill(tmp_RECO_rhoGrid,tmp_event_nPV);

			// jet mass
			tmp_GEN_mass=GenGroomedJet_mass_uncorr[0];
			tmp_PFCor_mass=JetPFCor_Mass[i_PFCorJet_matching_PFCHS];
			tmp_RECO_mass_uncorr=GroomedJet_mass_uncorr[0];
			tmp_RECO_mass_jec = GroomedJet_mass[0];
			tmp_RECO_mass_rhoArea=GroomedJet_mass_rhoArea[0];
			tmp_RECO_mass_rhoGArea=GroomedJet_mass_rhoGArea[0];
			tmp_RECO_mass_rho4Area=GroomedJet_mass_rho4Area[0];
			tmp_RECO_mass_rhoG4Area=GroomedJet_mass_rhoG4Area[0];
			tmp_RECO_mass_rhom4Area=GroomedJet_mass_rhom4Area[0];
			
			h1_PFCor_mass.Fill(tmp_PFCor_mass          );

			/*h1_GEN_mass.Fill(tmp_GEN_mass            ); 
			h1_RECO_mass_uncorr.Fill(tmp_RECO_mass_uncorr  );
			h1_RECO_mass_jec.Fill(tmp_RECO_mass_jec  );
			h1_RECO_mass_rhoArea.Fill(tmp_RECO_mass_rhoArea );
			h1_RECO_mass_rhoGArea.Fill(tmp_RECO_mass_rhoGArea);
			h1_RECO_mass_rho4Area.Fill(tmp_RECO_mass_rho4Area);
			h1_RECO_mass_rhoG4Area.Fill(tmp_RECO_mass_rhoG4Area);
			h1_RECO_mass_rhom4Area.Fill(tmp_RECO_mass_rhom4Area);

			h2_GEN_mass.Fill(tmp_GEN_mass, tmp_GEN_mass            ); 
			h2_RECO_mass_uncorr.Fill(tmp_GEN_mass, tmp_RECO_mass_uncorr  );
			h2_RECO_mass_jec.Fill(tmp_GEN_mass, tmp_RECO_mass_jec  );
			h2_RECO_mass_rhoArea.Fill(tmp_GEN_mass, tmp_RECO_mass_rhoArea );
			h2_RECO_mass_rhoGArea.Fill(tmp_GEN_mass, tmp_RECO_mass_rhoGArea);
			h2_RECO_mass_rho4Area.Fill(tmp_GEN_mass, tmp_RECO_mass_rho4Area);
			h2_RECO_mass_rhoG4Area.Fill(tmp_GEN_mass, tmp_RECO_mass_rhoG4Area);
			h2_RECO_mass_rhom4Area.Fill(tmp_GEN_mass, tmp_RECO_mass_rhom4Area);

			gr_GEN_mass.SetPoint(num_points, tmp_GEN_mass, tmp_GEN_mass            ); 
			gr_RECO_mass_uncorr.SetPoint(num_points,tmp_GEN_mass, tmp_RECO_mass_uncorr  );
			gr_RECO_mass_jec.SetPoint(num_points,tmp_GEN_mass, tmp_RECO_mass_jec  );
			gr_RECO_mass_rhoArea.SetPoint(num_points,tmp_GEN_mass, tmp_RECO_mass_rhoArea );
			gr_RECO_mass_rhoGArea.SetPoint(num_points,tmp_GEN_mass, tmp_RECO_mass_rhoGArea);
			gr_RECO_mass_rho4Area.SetPoint(num_points,tmp_GEN_mass, tmp_RECO_mass_rho4Area);
			gr_RECO_mass_rhoG4Area.SetPoint(num_points,tmp_GEN_mass, tmp_RECO_mass_rhoG4Area);
			gr_RECO_mass_rhom4Area.SetPoint(num_points,tmp_GEN_mass, tmp_RECO_mass_rhom4Area);*/
			ct_GEN_mass.Fill(tmp_GEN_mass, tmp_GEN_mass            ); 
			ct_RECO_mass_uncorr.Fill(tmp_GEN_mass, tmp_RECO_mass_uncorr  );
			ct_RECO_mass_jec.Fill(tmp_GEN_mass, tmp_RECO_mass_jec  );
			ct_RECO_mass_rhoArea.Fill(tmp_GEN_mass, tmp_RECO_mass_rhoArea );
			ct_RECO_mass_rhoGArea.Fill(tmp_GEN_mass, tmp_RECO_mass_rhoGArea);
			ct_RECO_mass_rho4Area.Fill(tmp_GEN_mass, tmp_RECO_mass_rho4Area);
			ct_RECO_mass_rhoG4Area.Fill(tmp_GEN_mass, tmp_RECO_mass_rhoG4Area);
			ct_RECO_mass_rhom4Area.Fill(tmp_GEN_mass, tmp_RECO_mass_rhom4Area);			

			h2_RECO_mass_jec_vs_PV.Fill(tmp_event_nPV, tmp_RECO_mass_jec  );

			Int_t tmp_number_JetCleansing_DiffMode=number_JetCleansing_DiffMode;
			if(tmp_number_JetCleansing_DiffMode==0)tmp_number_JetCleansing_DiffMode=NUM_JETCLEANSING_DIFFMODE;
			for(Int_t k=0;k<tmp_number_JetCleansing_DiffMode;k++){
				if (GroomedJet_mass_JetCleansing_DiffMode[k]>=0 && GroomedJet_pt_JetCleansing_DiffMode[k]>=0){
					/*vect_h1_RECO_mass_JetCleansing_DiffMode[k].Fill(GroomedJet_mass_JetCleansing_DiffMode[k]);
					vect_h2_RECO_mass_JetCleansing_DiffMode[k].Fill(tmp_GEN_mass, GroomedJet_mass_JetCleansing_DiffMode[k]);
					vect_gr_RECO_mass_JetCleansing_DiffMode[k].SetPoint(num_points,tmp_GEN_mass, GroomedJet_mass_JetCleansing_DiffMode[k]);

					vect_h1_RECO_pt_JetCleansing_DiffMode[k].Fill(GroomedJet_pt_JetCleansing_DiffMode[k]);
					vect_h2_RECO_pt_JetCleansing_DiffMode[k].Fill(tmp_GEN_pt, GroomedJet_pt_JetCleansing_DiffMode[k]);
					vect_gr_RECO_pt_JetCleansing_DiffMode[k].SetPoint(num_points,tmp_GEN_pt, GroomedJet_pt_JetCleansing_DiffMode[k]);
*/
					vect_ct_RECO_mass_JetCleansing_DiffMode[k].Fill(tmp_GEN_mass, GroomedJet_mass_JetCleansing_DiffMode[k]);
					//vect_ct_RECO_pt_JetCleansing_DiffMode[k].Fill(tmp_GEN_pt, GroomedJet_pt_JetCleansing_DiffMode[k]);

					//RECO rhohand
					ratio = GroomedJet_pt_JetCleansing_DiffMode[k]/tmp_GEN_pt;
					vect_mrt_RECO_pt_JetCleansing_DiffMode[k].Fill(tmp_GEN_eta, ratio);

					//tau2tau1
					if( GenGroomedJet_mass[0]>40 ){
						vect_ct_RECO_tau2tau1_JetCleansing_DiffMode[k].Fill(tmp_GEN_tau2tau1, GroomedJet_tau2tau1_JetCleansing_DiffMode[k]);
					}


					if(tmp_number_JetCleansing_DiffMode==NUM_JETCLEANSING_DIFFMODE)number_JetCleansing_DiffMode++;
				}else{ break;}
			}
			//num_points++;

		}
	}
	Draw_and_Save(h1_PFCor_jec_recogenptratio);
	Draw_and_Save(h1_PFCor_uncorr_recogenptratio);
	Draw_and_Save(h1_PFCor_afterL1_recogenptratio);
	Draw_and_Save(h1_PFCor_afterL2_recogenptratio);
	Draw_and_Save(h1_RECO_jec_recogenptratio);
	Draw_and_Save(h1_RECO_uncorr_recogenptratio);
	Draw_and_Save(h1_RECO_l1rhosw_recogenptratio);
	Draw_and_Save(h1_RECO_l1rhoHand_recogenptratio);
	Draw_and_Save(h1_RECO_l1rhoHand2_recogenptratio);
	Draw_and_Save(h1_RECO_l1rhoGrid_recogenptratio);
	Draw_and_Save(h1_RECO_rho4A_recogenptratio);
	Draw_and_Save(h1_RECO_rhom4A_recogenptratio);
	//Draw_and_Save(h1_RECO_JetCleansing_recogenptratio);

	Draw_and_Save(mrt_PFCor_jec_recogenptratio.get_mr_hist(ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save(mrt_PFCor_uncorr_recogenptratio.get_mr_hist(ratio_mrt_uncorr_min, ratio_mrt_uncorr_max));
	Draw_and_Save(mrt_PFCor_afterL1_recogenptratio.get_mr_hist(ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save(mrt_PFCor_afterL2_recogenptratio.get_mr_hist(ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save(mrt_RECO_jec_recogenptratio.get_mr_hist(ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save(mrt_RECO_uncorr_recogenptratio.get_mr_hist(ratio_mrt_uncorr_min, ratio_mrt_uncorr_max));
	Draw_and_Save(mrt_RECO_l1rhosw_recogenptratio.get_mr_hist(ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save(mrt_RECO_l1rhoHand_recogenptratio.get_mr_hist(ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save(mrt_RECO_l1rhoHand2_recogenptratio.get_mr_hist(ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save(mrt_RECO_l1rhoGrid_recogenptratio.get_mr_hist(ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save(mrt_RECO_rho4A_recogenptratio.get_mr_hist(ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save(mrt_RECO_rhom4A_recogenptratio.get_mr_hist(ratio_mrt_min, ratio_mrt_max));

	/*
	   h1_RECO_uncorr_recogenptratio.SetLineColor(1); h1_RECO_uncorr_recogenptratio.SetLineStyle(2); h1_RECO_uncorr_recogenptratio.SetLineWidth(2);
	   h1_RECO_l1rhoHand_recogenptratio.SetLineColor(2);h1_RECO_l1rhoHand_recogenptratio.SetLineStyle(2);h1_RECO_l1rhoHand_recogenptratio.SetLineWidth(2);
	   h1_RECO_rho4A_recogenptratio.SetLineColor(3);h1_RECO_rho4A_recogenptratio.SetLineStyle(1);h1_RECO_rho4A_recogenptratio.SetLineWidth(1);
	   h1_RECO_rhom4A_recogenptratio.SetLineColor(4);h1_RECO_rhom4A_recogenptratio.SetLineStyle(2);h1_RECO_rhom4A_recogenptratio.SetLineWidth(2);
	   h1_RECO_JetCleansing_recogenptratio.SetLineColor(6);h1_RECO_JetCleansing_recogenptratio.SetLineStyle(1);h1_RECO_JetCleansing_recogenptratio.SetLineWidth(1);
	   Draw_and_Save(h1_RECO_uncorr_recogenptratio, h1_RECO_l1rhoHand_recogenptratio, h1_RECO_rho4A_recogenptratio );
	   Draw_and_Save(h1_RECO_uncorr_recogenptratio, h1_RECO_rhom4A_recogenptratio, h1_RECO_JetCleansing_recogenptratio );

	   h1_GEN_pt.SetLineColor(1); 
	   h1_RECO_pt_uncorr.SetLineColor(1); h1_RECO_pt_uncorr.SetLineStyle(2); h1_RECO_pt_uncorr.SetLineWidth(2);
	   h1_RECO_pt_rhoArea.SetLineColor(2);h1_RECO_pt_rhoArea.SetLineStyle(2);h1_RECO_pt_rhoArea.SetLineWidth(2);
	   h1_RECO_pt_rho4Area.SetLineColor(3);h1_RECO_pt_rho4Area.SetLineStyle(1);
	   h1_RECO_pt_rhom4Area.SetLineColor(4);h1_RECO_pt_rhom4Area.SetLineStyle(2);h1_RECO_pt_rhom4Area.SetLineWidth(2);
	   h1_RECO_pt_JetCleansing.SetLineColor(6);h1_RECO_pt_JetCleansing.SetLineStyle(1);
	   Draw_and_Save(h1_GEN_pt, h1_RECO_pt_uncorr, h1_RECO_pt_rhoArea, h1_RECO_pt_rho4Area );
	   Draw_and_Save(h1_GEN_pt, h1_RECO_pt_uncorr, h1_RECO_pt_rhom4Area, h1_RECO_pt_JetCleansing );
	   */
	Draw_and_Save(h1_PFCor_area);
	Draw_and_Save(h1_RECO_area);

	h1_PFCor_area.SetLineColor(kBlue);
	h1_RECO_area.SetLineColor(kBlack); h1_RECO_area.SetLineStyle(2);
	Draw_and_Save(h1_PFCor_area, h1_RECO_area);

	//h1_rho_PFCor); 
	Draw_and_Save(h1_nPV);
	Draw_and_Save(h1_z_mass);
	Draw_and_Save(h1_z_Pt);
	Draw_and_Save(h1_muplus_Pt);

	Draw_and_Save(h1_GEN_eta);
	Draw_and_Save(h1_GEN_zjet_dr);
	Draw_and_Save(h1_GEN_zjet_dphi);
	/*Draw_and_Save(h1_GEN_rhoSW);
	  Draw_and_Save(h1_GEN_rhoHand);
	  Draw_and_Save(h1_GEN_rhoHand2);
	  Draw_and_Save(h1_GEN_rhoGrid);
	  Draw_and_Save(h2_GEN_rhoSW_vs_nPV);
	  Draw_and_Save(h2_GEN_rhoHand_vs_nPV);
	  Draw_and_Save(h2_GEN_rhoHand2_vs_nPV);
	  Draw_and_Save(h2_GEN_rhoGrid_vs_nPV);*/

	Draw_and_Save(ct_GEN_pt.get_hist1D());        
	Draw_and_Save(ct_RECO_pt_uncorr.get_hist1D());  
	Draw_and_Save(ct_RECO_pt_rhoArea.get_hist1D()); 
	Draw_and_Save(ct_RECO_pt_rhoGArea.get_hist1D()); 
	Draw_and_Save(ct_RECO_pt_rho4Area.get_hist1D()); 
	Draw_and_Save(ct_RECO_pt_rhom4Area.get_hist1D()); 

	Draw_and_Save(ct_GEN_pt.get_hist2D()           , Form("%g correlated", ct_GEN_pt.get_graph().GetCorrelationFactor())   );
	Draw_and_Save(ct_RECO_pt_uncorr.get_hist2D()   , Form("%g correlated", ct_RECO_pt_uncorr.get_graph().GetCorrelationFactor())				 );
	Draw_and_Save(ct_RECO_pt_rhoArea.get_hist2D()  , Form("%g correlated", ct_RECO_pt_rhoArea.get_graph().GetCorrelationFactor())				);
	Draw_and_Save(ct_RECO_pt_rhoGArea.get_hist2D() , Form("%g correlated", ct_RECO_pt_rhoGArea.get_graph().GetCorrelationFactor())				);
	Draw_and_Save(ct_RECO_pt_rho4Area.get_hist2D() , Form("%g correlated", ct_RECO_pt_rho4Area.get_graph().GetCorrelationFactor())				);
	Draw_and_Save(ct_RECO_pt_rhom4Area.get_hist2D(), Form("%g correlated", ct_RECO_pt_rhom4Area.get_graph().GetCorrelationFactor())				);

	Draw_and_Save(ct_RECO_tau2tau1.get_hist1D());
	Draw_and_Save(ct_RECO_tau2tau1_shapesubtract.get_hist1D());

	Draw_and_Save(ct_RECO_tau2tau1.get_hist1D_responce());
	Draw_and_Save(ct_RECO_tau2tau1_shapesubtract.get_hist1D_responce());
	//Draw_and_Save(ct_GEN_tau2tau1.get_hist1D());
	//Draw_and_Save(ct_GEN_tau2tau1_shapesubtract.get_hist1D());

	Draw_and_Save(ct_RECO_tau2tau1.get_hist2D(), Form("%g correlated",ct_RECO_tau2tau1.get_graph().GetCorrelationFactor()));
	Draw_and_Save(ct_RECO_tau2tau1_shapesubtract.get_hist2D(), Form("%g correlated",ct_RECO_tau2tau1_shapesubtract.get_graph().GetCorrelationFactor()));
	//Draw_and_Save(ct_GEN_tau2tau1.get_hist2D(), Form("%g correlated",ct_GEN_tau2tau1.get_graph().GetCorrelationFactor()));
	//Draw_and_Save(ct_GEN_tau2tau1_shapesubtract.get_hist2D(), Form("%g correlated",ct_GEN_tau2tau1_shapesubtract.get_graph().GetCorrelationFactor()));

	Draw_and_Save(h1_RECO_eta);
	Draw_and_Save(h1_RECO_zjet_dr);
	Draw_and_Save(h1_RECO_zjet_dphi);
	Draw_and_Save(h1_RECO_rhoSW);
	Draw_and_Save(h1_RECO_rhoHand);
	Draw_and_Save(h1_RECO_rhoHand2);
	Draw_and_Save(h1_RECO_rhoGrid);
	Draw_and_Save(h2_RECO_rhoSW_vs_nPV);
	Draw_and_Save(h2_RECO_rhoHand_vs_nPV);
	Draw_and_Save(h2_RECO_rhoHand2_vs_nPV);
	Draw_and_Save(h2_RECO_rhoGrid_vs_nPV);


	h1_RECO_rhoSW.SetLineColor(kBlue);
	h1_RECO_rhoHand.SetLineColor(kBlack);
	h1_RECO_rhoHand.SetLineStyle(2);
	Draw_and_Save(h1_RECO_rhoSW, h1_RECO_rhoHand);

	Draw_and_Save(h2_RECO_l1rhoHand_recogenptratio_vs_nPV);
	Draw_and_Save(h2_RECO_l1rhoHand_recogenptratio_vs_ptHand);
	Draw_and_Save(h2_RECO_l1rhoHand_recogenptratio_vs_eta);


	Draw_and_Save(h1_PFCor_mass              );
	Draw_and_Save(h2_RECO_mass_jec_vs_PV);
	Draw_and_Save(ct_GEN_mass.get_hist1D()                );
	Draw_and_Save(ct_RECO_mass_uncorr.get_hist1D()      );
	Draw_and_Save(ct_RECO_mass_jec.get_hist1D()      );
	Draw_and_Save(ct_RECO_mass_rhoArea.get_hist1D()     );
	Draw_and_Save(ct_RECO_mass_rhoGArea.get_hist1D()    );
	Draw_and_Save(ct_RECO_mass_rho4Area.get_hist1D()    );
	Draw_and_Save(ct_RECO_mass_rhoG4Area.get_hist1D()   );
	Draw_and_Save(ct_RECO_mass_rhom4Area.get_hist1D()   );

	Draw_and_Save(ct_GEN_mass.get_hist1D_responce()                );
	Draw_and_Save(ct_RECO_mass_uncorr.get_hist1D_responce()      );
	Draw_and_Save(ct_RECO_mass_jec.get_hist1D_responce()      );
	Draw_and_Save(ct_RECO_mass_rhoArea.get_hist1D_responce()     );
	Draw_and_Save(ct_RECO_mass_rhoGArea.get_hist1D_responce()    );
	Draw_and_Save(ct_RECO_mass_rho4Area.get_hist1D_responce()    );
	Draw_and_Save(ct_RECO_mass_rhoG4Area.get_hist1D_responce()   );
	Draw_and_Save(ct_RECO_mass_rhom4Area.get_hist1D_responce()   );


	Draw_and_Save(ct_GEN_mass.get_hist2D(), Form("%g correlated", ct_GEN_mass.get_graph().GetCorrelationFactor())   );
	Draw_and_Save(ct_RECO_mass_uncorr.get_hist2D()	, Form("%g correlated", ct_RECO_mass_uncorr.get_graph().GetCorrelationFactor())				 );
	Draw_and_Save(ct_RECO_mass_jec.get_hist2D()		, Form("%g correlated", ct_RECO_mass_jec.get_graph().GetCorrelationFactor())				 );
	Draw_and_Save(ct_RECO_mass_rhoArea.get_hist2D()	, Form("%g correlated", ct_RECO_mass_rhoArea.get_graph().GetCorrelationFactor())				);
	Draw_and_Save(ct_RECO_mass_rhoGArea.get_hist2D(), Form("%g correlated", ct_RECO_mass_rhoGArea.get_graph().GetCorrelationFactor())				);
	Draw_and_Save(ct_RECO_mass_rho4Area.get_hist2D(), Form("%g correlated", ct_RECO_mass_rho4Area.get_graph().GetCorrelationFactor())				);
	Draw_and_Save(ct_RECO_mass_rhoG4Area.get_hist2D(), Form("%g correlated", ct_RECO_mass_rhoG4Area.get_graph().GetCorrelationFactor())				);
	Draw_and_Save(ct_RECO_mass_rhom4Area.get_hist2D(), Form("%g correlated", ct_RECO_mass_rhom4Area.get_graph().GetCorrelationFactor())				);


	cout<<"number_JetCleansing_DiffMode="<<number_JetCleansing_DiffMode<<endl;
	if (number_JetCleansing_DiffMode>NUM_JETCLEANSING_DIFFMODE) number_JetCleansing_DiffMode=NUM_JETCLEANSING_DIFFMODE;
	for(Int_t k=0;k<number_JetCleansing_DiffMode;k++){
		Draw_and_Save(vect_ct_RECO_mass_JetCleansing_DiffMode[k].get_hist1D_responce() );
		Draw_and_Save(vect_ct_RECO_mass_JetCleansing_DiffMode[k].get_hist1D() );
		Draw_and_Save(vect_ct_RECO_mass_JetCleansing_DiffMode[k].get_hist2D(), Form("%g correlated", vect_ct_RECO_mass_JetCleansing_DiffMode[k].get_graph().GetCorrelationFactor()));

		//Draw_and_Save(vect_ct_RECO_pt_JetCleansing_DiffMode[k].get_hist1D_responce());
		//Draw_and_Save(vect_ct_RECO_pt_JetCleansing_DiffMode[k].get_hist1D());
		//Draw_and_Save(vect_ct_RECO_pt_JetCleansing_DiffMode[k].get_hist2D(), Form("%g correlated", vect_ct_RECO_pt_JetCleansing_DiffMode[k].get_graph().GetCorrelationFactor()));
		Draw_and_Save(vect_mrt_RECO_pt_JetCleansing_DiffMode[k].get_mr_hist(ratio_mrt_min, ratio_mrt_max)	);
		Draw_and_Save(vect_mrt_RECO_pt_JetCleansing_DiffMode[k].get_hist1D()	);

		Draw_and_Save(vect_ct_RECO_tau2tau1_JetCleansing_DiffMode[k].get_hist1D_responce());
		Draw_and_Save(vect_ct_RECO_tau2tau1_JetCleansing_DiffMode[k].get_hist2D(), Form("%g correlated",vect_ct_RECO_tau2tau1_JetCleansing_DiffMode[k].get_graph().GetCorrelationFactor()));
	}

}
