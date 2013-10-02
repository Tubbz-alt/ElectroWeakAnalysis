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

void Draw_and_Save(TH1D h1){
	h1.Write();
	TCanvas *c1 = new TCanvas(Form("c1_%s",h1.GetTitle()),Form("c1_%s",h1.GetTitle()),200,10,600,600);
	c1->cd();
	h1.Draw();
	//c1->Print(Form("%s.pdf",h1.GetTitle()));
	c1->Print(Form("%s.png",h1.GetTitle()));
	delete c1;
}
void Draw_and_Save(TH1D h1, TH1D h2){
	TCanvas *c1 = new TCanvas(Form("c1_%s_and_%s",h1.GetTitle(),h2.GetTitle()),Form("c1_%s_%s",h1.GetTitle(),h2.GetTitle()),200,10,600,600);
	c1->cd();
	h1.GetYaxis()->SetRangeUser(0., TMath::Max(h1.GetMaximum(), h2.GetMaximum())*1.2);
	h1.Draw();
	h2.Draw("same");
	//c1->Print(Form("%s.pdf",h1.GetTitle()));
	c1->Print(Form("%s_and_%s.png",h1.GetTitle(),h2.GetTitle()));
	delete c1;
}
void Draw_and_Save(TH1D h1, TH1D h2, TH1D h3){
	TCanvas *c1 = new TCanvas(Form("c1_%s_and_%s_and_%s",h1.GetTitle(),h2.GetTitle(),h3.GetTitle()),Form("c1_%s_and_%s_and_%s",h1.GetTitle(),h2.GetTitle(),h3.GetTitle()),200,10,600,600);
	c1->cd();
	h1.GetYaxis()->SetRangeUser(0., TMath::Max(h1.GetMaximum(), h3.GetMaximum())*1.2);
	h1.Draw();
	h2.Draw("same");
	h3.Draw("same");
	//c1->Print(Form("%s.pdf",h1.GetTitle()));
	c1->Print(Form("%s_and_%s_and_%s.png",h1.GetTitle(),h2.GetTitle(),h3.GetTitle()));
	delete c1;
}
void Draw_and_Save(TH1D h1, TH1D h2, TH1D h3, TH1D h4){
	TCanvas *c1 = new TCanvas(Form("c1_%s_and_%s_and_%s_and_%s",h1.GetTitle(),h2.GetTitle(),h3.GetTitle(),h4.GetTitle()), Form("c1_%s_and_%s_and_%s_and_%s",h1.GetTitle(),h2.GetTitle(),h3.GetTitle(),h4.GetTitle()),200,10,600,600);
	c1->cd();
	h1.GetYaxis()->SetRangeUser(0., TMath::Max(h1.GetMaximum(), h4.GetMaximum())*1.2);
	h1.Draw();
	h2.Draw("same");
	h3.Draw("same");
	h4.Draw("same");
	//c1->Print(Form("%s.pdf",h1.GetTitle()));
	c1->Print(Form("%s_and_%s_and_%s_and_%s.png",h1.GetTitle(),h2.GetTitle(),h3.GetTitle(),h4.GetTitle()));
	delete c1;
}
void Draw_and_Save(TH2D h2, char* addtional_info=""){
	h2.Write();
	TCanvas *c1 = new TCanvas(Form("c1_%s",h2.GetTitle()),Form("c1_%s",h2.GetTitle()),200,10,600,600);
	c1->cd();
	h2.Draw("box");
	if( addtional_info){
		TLatex tl;
		tl.SetTextSize(0.04 ); tl.SetTextAlign(13);
		tl.DrawLatex(h2.GetXaxis()->GetXmin()*0.9+h2.GetXaxis()->GetXmax()*0.1,h2.GetYaxis()->GetXmin()*0.1+h2.GetYaxis()->GetXmax()*0.9,addtional_info);
	};
	//c1->Print(Form("%s.pdf",h2.GetTitle()));
	c1->Print(Form("%s.png",h2.GetTitle()));
	delete c1;
}

Bool_t MyClass::Select()
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   if(!( Z_pt>100 && Z_mass>65 && Z_mass<105 ))return 0;
   if(!( GenGroomedJet_AK8_GEN_pt[0]>100 ))return 0;

   Double_t tmp_AK5_GEN_eta = GenGroomedJet_AK5_GEN_eta[0];
   Double_t tmp_AK5_GEN_phi = GenGroomedJet_AK5_GEN_phi[0];
   Double_t tmpDeltaR_Vj = TMath::Sqrt( (Z_eta-tmp_AK5_GEN_eta)*(Z_eta-tmp_AK5_GEN_eta) + (Z_phi-tmp_AK5_GEN_phi)*(Z_phi-tmp_AK5_GEN_phi) );

   if(!( tmpDeltaR_Vj>2.0 ))return 0;
   return 1;
}

void MyClass::Loop(){
	//LoopAK5();
	LoopAK8();
}

void MyClass::LoopAK5()
{
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();

	Double_t tmp_AK5_GEN_pt=0.;
	Double_t tmp_AK5_GEN_eta=0.;
	Double_t tmp_AK5_GEN_phi=0.;
	Double_t tmp_AK5_GEN_rhoSW=0.;
	Double_t tmp_AK5_GEN_rhoHand=0.;
	Double_t tmp_AK5_GEN_rhoHand2=0.;
	Double_t tmp_AK5_GEN_rhoGrid=0.;
	
	Double_t tmp_PFCor_pt=0.;
	Double_t tmp_PFCor_Pt_uncorr=0.;
	Double_t tmp_PFCor_Pt_afterL1=0.;
	Double_t tmp_PFCor_Pt_afterL2=0.;

	Double_t tmp_AK5_PF_pt=0.;
	Double_t tmp_AK5_PF_eta=0.;
	Double_t tmp_AK5_PF_phi=0.;
	Double_t tmp_AK5_PF_pt_uncorr=0.;
	Double_t tmp_AK5_PF_pt_L1_rhoSW=0.;
	Double_t tmp_AK5_PF_pt_L1_rhoHand=0.;
	Double_t tmp_AK5_PF_pt_L1_rhoHand2=0.;
	Double_t tmp_AK5_PF_pt_L1_rhoGrid=0.;
	Double_t tmp_AK5_PF_rhoSW=0.;
	Double_t tmp_AK5_PF_rhoHand=0.;
	Double_t tmp_AK5_PF_rhoHand2=0.;
	Double_t tmp_AK5_PF_rhoGrid=0.;

	Double_t tmp_AK5_PFCHS_pt=0.;
	Double_t tmp_AK5_PFCHS_eta=0.;
	Double_t tmp_AK5_PFCHS_phi=0.;
	Double_t tmp_AK5_PFCHS_pt_uncorr=0.;
	Double_t tmp_AK5_PFCHS_pt_L1_rhoSW=0.;
	Double_t tmp_AK5_PFCHS_pt_L1_rhoHand=0.;
	Double_t tmp_AK5_PFCHS_pt_L1_rhoHand2=0.;
	Double_t tmp_AK5_PFCHS_pt_L1_rhoGrid=0.;
	Double_t tmp_AK5_PFCHS_rhoSW=0.;
	Double_t tmp_AK5_PFCHS_rhoHand=0.;
	Double_t tmp_AK5_PFCHS_rhoHand2=0.;
	Double_t tmp_AK5_PFCHS_rhoGrid=0.;

	//Double_t tmp_Z_eta=0.; Double_t tmp_Z_phi=0.;

	Double_t tmp_event_nPV=0.;

	Double_t ratio=0.;
	Double_t dr=0.; // Delta R
	Double_t dphi=0.; // Delta Phi

	//jet mass
	Double_t tmp_AK5_GEN_mass=0.;
	Double_t tmp_PFCor_mass=0.;
	Double_t tmp_AK5_PF_mass_uncorr=0.;
	Double_t tmp_AK5_PF_mass_rhoArea=0.;
	Double_t tmp_AK5_PF_mass_rhoGArea=0.;
	Double_t tmp_AK5_PF_mass_rho4Area=0.;
	Double_t tmp_AK5_PF_mass_rhoG4Area=0.;
	Double_t tmp_AK5_PF_mass_rhom4Area=0.;
	Double_t tmp_AK5_PF_mass_JetCleansingATLASjvf=0.;
	Double_t tmp_AK5_PF_mass_JetCleansingATLASlin=0.;
	Double_t tmp_AK5_PF_mass_JetCleansingATLASgau=0.;
	Double_t tmp_AK5_PF_mass_JetCleansingCMSjvf=0.;
	Double_t tmp_AK5_PF_mass_JetCleansingCMSlin=0.;
	Double_t tmp_AK5_PF_mass_JetCleansingCMSgau=0.;


	Double_t tmp_AK5_PFCHS_mass_uncorr=0.;
	Double_t tmp_AK5_PFCHS_mass_rhoArea=0.;
	Double_t tmp_AK5_PFCHS_mass_rhoGArea=0.;
	Double_t tmp_AK5_PFCHS_mass_rho4Area=0.;
	Double_t tmp_AK5_PFCHS_mass_rhoG4Area=0.;
	Double_t tmp_AK5_PFCHS_mass_rhom4Area=0.;

	Double_t rhomin=0.; Double_t rhomax=50.;

	TH1D h1_nPV("h1_nPV","h1_nPV",50,0,50);
	TH1D h1_z_mass("h1_z_mass","h1_z_mass",60,60,120);
	TH1D h1_muplus_Pt("h1_muplus_Pt","h1_muplus_Pt",40,0,400);

	TH1D h1_AK5_GEN_pt("h1_AK5_GEN_pt","h1_AK5_GEN_pt",50,0,200); h1_AK5_GEN_pt.SetLineColor(kRed);
	TH1D h1_AK5_GEN_eta("h1_AK5_GEN_eta","h1_AK5_GEN_eta",50,-2.5,2.5); h1_AK5_GEN_eta.SetLineColor(kRed);
	TH1D h1_AK5_GEN_phi("h1_AK5_GEN_phi","h1_AK5_GEN_phi",50,-4,4); h1_AK5_GEN_phi.SetLineColor(kRed);
	TH1D h1_AK5_GEN_zjet_dr("h1_AK5_GEN_zjet_dr","h1_AK5_GEN_zjet_dr",50,0,10); h1_AK5_GEN_zjet_dr.SetLineColor(kRed);
	TH1D h1_AK5_GEN_zjet_dphi("h1_AK5_GEN_zjet_dphi","h1_AK5_GEN_zjet_dphi",50,0,5); h1_AK5_GEN_zjet_dphi.SetLineColor(kRed);
	TH1D h1_AK5_GEN_rhoSW("h1_AK5_GEN_rhoSW","h1_AK5_GEN_rhoSW",50,rhomin,rhomax); h1_AK5_GEN_rhoSW.SetLineColor(kRed);
	TH1D h1_AK5_GEN_rhoHand("h1_AK5_GEN_rhoHand","h1_AK5_GEN_rhoHand",50,rhomin,rhomax); h1_AK5_GEN_rhoHand.SetLineColor(kRed);
	TH1D h1_AK5_GEN_rhoHand2("h1_AK5_GEN_rhoHand2","h1_AK5_GEN_rhoHand2",50,rhomin,rhomax); h1_AK5_GEN_rhoHand2.SetLineColor(kRed);
	TH1D h1_AK5_GEN_rhoGrid("h1_AK5_GEN_rhoGrid","h1_AK5_GEN_rhoGrid",50,rhomin,rhomax); h1_AK5_GEN_rhoGrid.SetLineColor(kRed);
	TH2D h2_AK5_GEN_rhoSW_vs_nPV("h2_AK5_GEN_rhoSW_vs_nPV","h2_AK5_GEN_rhoSW_vs_nPV",50,rhomin,rhomax,50,0,50); h2_AK5_GEN_rhoSW_vs_nPV.SetMarkerColor(kRed);
	TH2D h2_AK5_GEN_rhoHand_vs_nPV("h2_AK5_GEN_rhoHand_vs_nPV","h2_AK5_GEN_rhoHand_vs_nPV",50,rhomin,rhomax,50,0,50); h2_AK5_GEN_rhoHand_vs_nPV.SetMarkerColor(kRed);
	TH2D h2_AK5_GEN_rhoHand2_vs_nPV("h2_AK5_GEN_rhoHand2_vs_nPV","h2_AK5_GEN_rhoHand2_vs_nPV",50,rhomin,rhomax,50,0,50); h2_AK5_GEN_rhoHand2_vs_nPV.SetMarkerColor(kRed);
	TH2D h2_AK5_GEN_rhoGrid_vs_nPV("h2_AK5_GEN_rhoGrid_vs_nPV","h2_AK5_GEN_rhoGrid_vs_nPV",50,rhomin,rhomax,50,0,50); h2_AK5_GEN_rhoGrid_vs_nPV.SetMarkerColor(kRed);

	TH1D h1_AK5_PF_pt_uncorr("h1_AK5_PF_pt_uncorr","h1_AK5_PF_pt_uncorr",50,0,200);
	TH1D h1_AK5_PF_eta("h1_AK5_PF_eta","h1_AK5_PF_eta",50,-2.5,2.5);
	TH1D h1_AK5_PF_phi("h1_AK5_PF_phi","h1_AK5_PF_phi",50,-4,4);
	TH1D h1_AK5_PF_zjet_dr("h1_AK5_PF_zjet_dr","h1_AK5_PF_zjet_dr",50,0,10);
	TH1D h1_AK5_PF_zjet_dphi("h1_AK5_PF_zjet_dphi","h1_AK5_PF_zjet_dphi",50,0,5);
	TH1D h1_AK5_PF_rhoSW("h1_AK5_PF_rhoSW","h1_AK5_PF_rhoSW",50,rhomin,rhomax);
	TH1D h1_AK5_PF_rhoHand("h1_AK5_PF_rhoHand","h1_AK5_PF_rhoHand",50,rhomin,rhomax);
	TH1D h1_AK5_PF_rhoHand2("h1_AK5_PF_rhoHand2","h1_AK5_PF_rhoHand2",50,rhomin,rhomax);
	TH1D h1_AK5_PF_rhoGrid("h1_AK5_PF_rhoGrid","h1_AK5_PF_rhoGrid",50,rhomin,rhomax);
	TH2D h2_AK5_PF_rhoSW_vs_nPV("h2_AK5_PF_rhoSW_vs_nPV","h2_AK5_PF_rhoSW_vs_nPV",50,rhomin,rhomax,50,0,50);
	TH2D h2_AK5_PF_rhoHand_vs_nPV("h2_AK5_PF_rhoHand_vs_nPV","h2_AK5_PF_rhoHand_vs_nPV",50,rhomin,rhomax,50,0,50);
	TH2D h2_AK5_PF_rhoHand2_vs_nPV("h2_AK5_PF_rhoHand2_vs_nPV","h2_AK5_PF_rhoHand2_vs_nPV",50,rhomin,rhomax,50,0,50);
	TH2D h2_AK5_PF_rhoGrid_vs_nPV("h2_AK5_PF_rhoGrid_vs_nPV","h2_AK5_PF_rhoGrid_vs_nPV",50,rhomin,rhomax,50,0,50);


	TH1D h1_AK5_PFCHS_pt_uncorr("h1_AK5_PFCHS_pt_uncorr","h1_AK5_PFCHS_pt_uncorr",50,0,200);
	TH1D h1_AK5_PFCHS_eta("h1_AK5_PFCHS_eta","h1_AK5_PFCHS_eta",50,-2.5,2.5);
	TH1D h1_AK5_PFCHS_phi("h1_AK5_PFCHS_phi","h1_AK5_PFCHS_phi",50,-4,4);
	TH1D h1_AK5_PFCHS_zjet_dr("h1_AK5_PFCHS_zjet_dr","h1_AK5_PFCHS_zjet_dr",50,0,10);
	TH1D h1_AK5_PFCHS_zjet_dphi("h1_AK5_PFCHS_zjet_dphi","h1_AK5_PFCHS_zjet_dphi",50,0,5);
	TH1D h1_AK5_PFCHS_rhoSW("h1_AK5_PFCHS_rhoSW","h1_AK5_PFCHS_rhoSW",50,rhomin,rhomax);
	TH1D h1_AK5_PFCHS_rhoHand("h1_AK5_PFCHS_rhoHand","h1_AK5_PFCHS_rhoHand",50,rhomin,rhomax);
	TH1D h1_AK5_PFCHS_rhoHand2("h1_AK5_PFCHS_rhoHand2","h1_AK5_PFCHS_rhoHand2",50,rhomin,rhomax);
	TH1D h1_AK5_PFCHS_rhoGrid("h1_AK5_PFCHS_rhoGrid","h1_AK5_PFCHS_rhoGrid",50,rhomin,rhomax);
	TH2D h2_AK5_PFCHS_rhoSW_vs_nPV("h2_AK5_PFCHS_rhoSW_vs_nPV","h2_AK5_PFCHS_rhoSW_vs_nPV",50,rhomin,rhomax,50,0,50);
	TH2D h2_AK5_PFCHS_rhoHand_vs_nPV("h2_AK5_PFCHS_rhoHand_vs_nPV","h2_AK5_PFCHS_rhoHand_vs_nPV",50,rhomin,rhomax,50,0,50);
	TH2D h2_AK5_PFCHS_rhoHand2_vs_nPV("h2_AK5_PFCHS_rhoHand2_vs_nPV","h2_AK5_PFCHS_rhoHand2_vs_nPV",50,rhomin,rhomax,50,0,50);
	TH2D h2_AK5_PFCHS_rhoGrid_vs_nPV("h2_AK5_PFCHS_rhoGrid_vs_nPV","h2_AK5_PFCHS_rhoGrid_vs_nPV",50,rhomin,rhomax,50,0,50);

	TH2D h2_AK5_PFCHS_ratioHand_vs_nPV("h2_AK5_PFCHS_ratioHand_vs_nPV","h2_AK5_PFCHS_ratioHand_vs_nPV",50,0,50, 50, 0, 2.5);
	TH2D h2_AK5_PFCHS_ratioHand_vs_ptHand("h2_AK5_PFCHS_ratioHand_vs_ptHand","h2_AK5_PFCHS_ratioHand_vs_ptHand",50,0,200, 50, 0, 2.5);
	TH2D h2_AK5_PFCHS_ratioHand_vs_eta("h2_AK5_PFCHS_ratioHand_vs_eta","h2_AK5_PFCHS_ratioHand_vs_eta", 50, -2.5, 2.5, 50, 0, 2.5);



	TH1D h1_PF_match("h1_PF_match","h1_PF_match",50,0,1.); h1_PF_match.SetLineColor(kRed);
	TH1D h1_PFCHS_match("h1_PFCHS_match","h1_PFCHS_match",50,0,1.);

	// Jet Pt
	Int_t nbin_pt=40;Double_t jetpt_min=50;Double_t jetpt_max=450.;
	TH1D h1_PFCor_Pt_afterL1("h1_PFCor_Pt_afterL1","h1_PFCor_Pt_afterL1;Jet Pt",nbin_pt, jetpt_min, jetpt_max);
	TH1D h1_AK5_PF_Pt_l1_rhoHand("h1_AK5_PF_Pt_l1_rhoHand","h1_AK5_PF_Pt_l1_rhoHand;Jet Pt",nbin_pt, jetpt_min, jetpt_max);
	TH1D h1_AK5_PFCHS_Pt_l1_rhoHand("h1_AK5_PFCHS_Pt_l1_rhoHand","h1_AK5_PFCHS_Pt_l1_rhoHand;Jet Pt",nbin_pt, jetpt_min, jetpt_max);
	// JetPt/GenPt
	Int_t nbin_ratio=20; Double_t ratio_min=0.3; Double_t ratio_max=1.7; 
	TH1D h1_PFCor("h1_PFCor","h1_PFCor",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_PFCor_uncorr("h1_PFCor_uncorr","h1_PFCor_uncorr",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_PFCor_afterL1("h1_PFCor_afterL1","h1_PFCor_afterL1",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_PFCor_afterL2("h1_PFCor_afterL2","h1_PFCor_afterL2",nbin_ratio, ratio_min, ratio_max);

	TH1D h1_AK5_PF("h1_AK5_PF","h1_AK5_PF",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_AK5_PF_uncorr("h1_AK5_PF_uncorr","h1_AK5_PF_uncorr",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_AK5_PF_l1_rhoSW("h1_AK5_PF_l1_rhoSW","h1_AK5_PF_l1_rhoSW",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_AK5_PF_l1_rhoHand("h1_AK5_PF_l1_rhoHand","h1_AK5_PF_l1_rhoHand",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_AK5_PF_l1_rhoHand2("h1_AK5_PF_l1_rhoHand2","h1_AK5_PF_l1_rhoHand2",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_AK5_PF_l1_rhoGrid("h1_AK5_PF_l1_rhoGrid","h1_AK5_PF_l1_rhoGrid",nbin_ratio, ratio_min, ratio_max);

	TH1D h1_AK5_PFCHS("h1_AK5_PFCHS","h1_AK5_PFCHS",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_AK5_PFCHS_uncorr("h1_AK5_PFCHS_uncorr","h1_AK5_PFCHS_uncorr",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_AK5_PFCHS_l1_rhoSW("h1_AK5_PFCHS_l1_rhoSW","h1_AK5_PFCHS_l1_rhoSW",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_AK5_PFCHS_l1_rhoHand("h1_AK5_PFCHS_l1_rhoHand","h1_AK5_PFCHS_l1_rhoHand",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_AK5_PFCHS_l1_rhoHand2("h1_AK5_PFCHS_l1_rhoHand2","h1_AK5_PFCHS_l1_rhoHand2",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_AK5_PFCHS_l1_rhoGrid("h1_AK5_PFCHS_l1_rhoGrid","h1_AK5_PFCHS_l1_rhoGrid",nbin_ratio, ratio_min, ratio_max);

	//area
	TH1D h1_PFCor_area("h1_PFCor_area","h1_PFCor_area",50,0.6,1.1);
	TH1D h1_AK5_PF_area("h1_AK5_PF_area","h1_AK5_PF_area",50,0.6,1.1);
	TH1D h1_AK5_PFCHS_area("h1_AK5_PFCHS_area","h1_AK5_PFCHS_area",50,0.6,1.1);

	//rho
	TH1D h1_rho_PFCor("h1_rho_PFCor","h1_rho_PFCor",30,0,30);


	//mass
	int nbin_mass=20;double jetmass_min=0;double jetmass_max=20.;
	TH1D h1_AK5_GEN_mass("h1_AK5_GEN_mass","h1_AK5_GEN_mass;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	h1_AK5_GEN_mass.SetLineColor(kRed);
	TH1D h1_PFCor_mass("h1_PFCor_mass","h1_PFCor_mass;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	h1_PFCor_mass.SetLineColor(kBlack);
	TH1D h1_AK5_PF_mass_uncorr("h1_AK5_PF_mass_uncorr","h1_AK5_PF_mass_uncorr;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK5_PF_mass_rhoArea("h1_AK5_PF_mass_rhoArea","h1_AK5_PF_mass_rhoArea;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK5_PF_mass_rhoGArea("h1_AK5_PF_mass_rhoGArea","h1_AK5_PF_mass_rhoGArea;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK5_PF_mass_rho4Area("h1_AK5_PF_mass_rho4Area","h1_AK5_PF_mass_rho4Area;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK5_PF_mass_rhoG4Area("h1_AK5_PF_mass_rhoG4Area","h1_AK5_PF_mass_rhoG4Area;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK5_PF_mass_rhom4Area("h1_AK5_PF_mass_rhom4Area","h1_AK5_PF_mass_rhom4Area;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK5_PF_mass_JetCleansingATLASjvf("h1_AK5_PF_mass_JetCleansingATLASjvf","h1_AK5_PF_mass_JetCleansingATLASjvf;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK5_PF_mass_JetCleansingATLASlin("h1_AK5_PF_mass_JetCleansingATLASlin","h1_AK5_PF_mass_JetCleansingATLASlin;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK5_PF_mass_JetCleansingATLASgau("h1_AK5_PF_mass_JetCleansingATLASgau","h1_AK5_PF_mass_JetCleansingATLASgau;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK5_PF_mass_JetCleansingCMSjvf("h1_AK5_PF_mass_JetCleansingCMSjvf","h1_AK5_PF_mass_JetCleansingCMSjvf;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK5_PF_mass_JetCleansingCMSlin("h1_AK5_PF_mass_JetCleansingCMSlin","h1_AK5_PF_mass_JetCleansingCMSlin;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK5_PF_mass_JetCleansingCMSgau("h1_AK5_PF_mass_JetCleansingCMSgau","h1_AK5_PF_mass_JetCleansingCMSgau;jet mass;",nbin_mass,jetmass_min,jetmass_max);

	TH1D h1_AK5_PFCHS_mass_uncorr(   "h1_AK5_PFCHS_mass_uncorr",   "h1_AK5_PFCHS_mass_uncorr;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK5_PFCHS_mass_rhoArea(  "h1_AK5_PFCHS_mass_rhoArea",  "h1_AK5_PFCHS_mass_rhoArea;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK5_PFCHS_mass_rhoGArea( "h1_AK5_PFCHS_mass_rhoGArea", "h1_AK5_PFCHS_mass_rhoGArea;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK5_PFCHS_mass_rho4Area( "h1_AK5_PFCHS_mass_rho4Area", "h1_AK5_PFCHS_mass_rho4Area;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK5_PFCHS_mass_rhoG4Area("h1_AK5_PFCHS_mass_rhoG4Area","h1_AK5_PFCHS_mass_rhoG4Area;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK5_PFCHS_mass_rhom4Area("h1_AK5_PFCHS_mass_rhom4Area","h1_AK5_PFCHS_mass_rhom4Area;jet mass;",nbin_mass,jetmass_min,jetmass_max);


	/*	TH1D h1_rhoSW_PF("h1_rhoSW_PF","h1_rhoSW_PF",30,0,30);
		TH1D h1_rhoHand_PF("h1_rhoHand_PF","h1_rhoHand_PF",30,0,30);
		TH1D h1_rhoHand2_PF("h1_rhoHand2_PF","h1_rhoHand2_PF",30,0,30);
		TH1D h1_rhoGrid_PF("h1_rhoGrid_PF","h1_rhoGrid_PF",30,0,30);

		TH1D h1_rhoSW_PFCHS("h1_rhoSW_PFCHS","h1_rhoSW_PFCHS",30,0,30);
		TH1D h1_rhoHand_PFCHS("h1_rhoHand_PFCHS","h1_rhoHand_PFCHS",30,0,30);
		TH1D h1_rhoHand2_PFCHS("h1_rhoHand2_PFCHS","h1_rhoHand2_PFCHS",30,0,30);
		TH1D h1_rhoGrid_PFCHS("h1_rhoGrid_PFCHS","h1_rhoGrid_PFCHS",30,0,30);
		*/
	// For GEN-RECO matching
	Double_t gen_jet_eta=0.; Double_t gen_jet_phi=0.;
	Double_t PFCHS_eta=0.; Double_t PFCHS_phi=0.;
	Double_t PFCor_jet_eta=0.; Double_t PFCor_jet_phi=0.;

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;
		if (!Select())continue;

		gen_jet_eta=GenGroomedJet_AK5_GEN_eta[0];
		gen_jet_phi=GenGroomedJet_AK5_GEN_phi[0];

		PFCHS_eta=GroomedJet_AK5_PFCHS_eta[0];
		PFCHS_phi=GroomedJet_AK5_PFCHS_phi[0];

		//Gen Jet matching with PFCHS_uncerr
		if( TMath::Sqrt( (gen_jet_eta-PFCHS_eta)*(gen_jet_eta-PFCHS_eta) + (gen_jet_phi-PFCHS_phi)*(gen_jet_phi-PFCHS_phi) ) <0.3 )
		{
			int i_PFCorJet_matching_PFCHS=0;//_matching with leading PFCHS_uncerr
			for(;i_PFCorJet_matching_PFCHS<numPFCorJets;i_PFCorJet_matching_PFCHS++){
				PFCor_jet_eta=JetPFCor_Eta[i_PFCorJet_matching_PFCHS]; PFCor_jet_phi=JetPFCor_Phi[i_PFCorJet_matching_PFCHS];
				if ( equal(PFCHS_eta, PFCor_jet_eta) && equal(PFCHS_phi, PFCor_jet_phi) )break;
			}

			tmp_AK5_GEN_pt = GenGroomedJet_AK5_GEN_pt[0];
			tmp_AK5_GEN_eta = GenGroomedJet_AK5_GEN_eta[0];
			tmp_AK5_GEN_phi = GenGroomedJet_AK5_GEN_phi[0];

			tmp_PFCor_pt = JetPFCor_Pt[i_PFCorJet_matching_PFCHS];
			tmp_PFCor_Pt_uncorr = JetPFCor_Pt_uncorr[i_PFCorJet_matching_PFCHS];
			tmp_PFCor_Pt_afterL1 = JetPFCor_Pt_afterL1[i_PFCorJet_matching_PFCHS];
			tmp_PFCor_Pt_afterL2 = JetPFCor_Pt_afterL2[i_PFCorJet_matching_PFCHS];

			tmp_AK5_PF_pt = GroomedJet_AK5_PF_pt[0];
			tmp_AK5_PF_eta = GroomedJet_AK5_PF_eta[0];
			tmp_AK5_PF_phi = GroomedJet_AK5_PF_phi[0];
			tmp_AK5_PF_pt_uncorr = GroomedJet_AK5_PF_pt_uncorr[0];
			tmp_AK5_PF_pt_L1_rhoSW = GroomedJet_AK5_PF_pt_L1_rhoSW[0];
			tmp_AK5_PF_pt_L1_rhoHand = GroomedJet_AK5_PF_pt_L1_rhoHand[0];
			tmp_AK5_PF_pt_L1_rhoHand2 = GroomedJet_AK5_PF_pt_L1_rhoHand2[0];
			tmp_AK5_PF_pt_L1_rhoGrid = GroomedJet_AK5_PF_pt_L1_rhoGrid[0];

			tmp_AK5_PFCHS_pt = GroomedJet_AK5_PFCHS_pt[0];
			tmp_AK5_PFCHS_eta = GroomedJet_AK5_PFCHS_eta[0];
			tmp_AK5_PFCHS_phi = GroomedJet_AK5_PFCHS_phi[0];
			tmp_AK5_PFCHS_pt_uncorr = GroomedJet_AK5_PFCHS_pt_uncorr[0];
			tmp_AK5_PFCHS_pt_L1_rhoSW = GroomedJet_AK5_PFCHS_pt_L1_rhoSW[0];
			tmp_AK5_PFCHS_pt_L1_rhoHand = GroomedJet_AK5_PFCHS_pt_L1_rhoHand[0];
			tmp_AK5_PFCHS_pt_L1_rhoHand2 = GroomedJet_AK5_PFCHS_pt_L1_rhoHand2[0];
			tmp_AK5_PFCHS_pt_L1_rhoGrid = GroomedJet_AK5_PFCHS_pt_L1_rhoGrid[0];

			//tmp_Z_eta = Z_eta; tmp_Z_phi = Z_phi;

			tmp_event_nPV = event_nPV;

			tmp_AK5_GEN_rhoSW = GenGroomedJet_AK5_GEN_rhoSW;
			tmp_AK5_GEN_rhoHand = GenGroomedJet_AK5_GEN_rhohand;
			tmp_AK5_GEN_rhoHand2 = GenGroomedJet_AK5_GEN_rhohand2;
			tmp_AK5_GEN_rhoGrid = GenGroomedJet_AK5_GEN_rhogrid;

			tmp_AK5_PF_rhoSW = GroomedJet_AK5_PF_rhoSW;
			tmp_AK5_PF_rhoHand = GroomedJet_AK5_PF_rhohand;
			tmp_AK5_PF_rhoHand2 = GroomedJet_AK5_PF_rhohand2;
			tmp_AK5_PF_rhoGrid = GroomedJet_AK5_PF_rhogrid;

			tmp_AK5_PFCHS_rhoSW = GroomedJet_AK5_PFCHS_rhoSW;
			tmp_AK5_PFCHS_rhoHand = GroomedJet_AK5_PFCHS_rhohand;
			tmp_AK5_PFCHS_rhoHand2 = GroomedJet_AK5_PFCHS_rhohand2;
			tmp_AK5_PFCHS_rhoGrid = GroomedJet_AK5_PFCHS_rhogrid;

			//============= begin to fill hist ==============
			//PFCor
			ratio = tmp_PFCor_pt/tmp_AK5_GEN_pt;
			h1_PFCor.Fill(ratio);
			//PFCor uncorr
			ratio = tmp_PFCor_Pt_uncorr/tmp_AK5_GEN_pt;
			h1_PFCor_uncorr.Fill(ratio);
			//PFCor L1
			ratio = tmp_PFCor_Pt_afterL1/tmp_AK5_GEN_pt;
			h1_PFCor_afterL1.Fill(ratio);
			h1_PFCor_Pt_afterL1.Fill(tmp_PFCor_Pt_afterL1);
			//PFCor L2
			ratio = tmp_PFCor_Pt_afterL2/tmp_AK5_GEN_pt;
			h1_PFCor_afterL2.Fill(ratio);

			//AK5 PF
			ratio = tmp_AK5_PF_pt/tmp_AK5_GEN_pt;
			h1_AK5_PF.Fill(ratio);
			//AK5 uncorr
			ratio = tmp_AK5_PF_pt_uncorr/tmp_AK5_GEN_pt;
			h1_AK5_PF_uncorr.Fill(ratio);
			//AK5 rhoSW
			ratio = tmp_AK5_PF_pt_L1_rhoSW/tmp_AK5_GEN_pt;
			h1_AK5_PF_l1_rhoSW.Fill(ratio);
			//AK5 rhohand
			ratio = tmp_AK5_PF_pt_L1_rhoHand/tmp_AK5_GEN_pt;
			h1_AK5_PF_l1_rhoHand.Fill(ratio);
			h1_AK5_PF_Pt_l1_rhoHand.Fill(tmp_AK5_PF_pt_L1_rhoHand);


			h2_AK5_PFCHS_ratioHand_vs_nPV.Fill(tmp_event_nPV, ratio);
			h2_AK5_PFCHS_ratioHand_vs_ptHand.Fill(tmp_AK5_PF_pt_L1_rhoHand, ratio);
			h2_AK5_PFCHS_ratioHand_vs_eta.Fill(tmp_AK5_PFCHS_eta, ratio );


			//AK5 rhohand2
			ratio = tmp_AK5_PF_pt_L1_rhoHand2/tmp_AK5_GEN_pt;
			h1_AK5_PF_l1_rhoHand2.Fill(ratio);
			//AK5 rhogrid
			ratio = tmp_AK5_PF_pt_L1_rhoGrid/tmp_AK5_GEN_pt;
			h1_AK5_PF_l1_rhoGrid.Fill(ratio);

			//AK5 PFCHS
			ratio = tmp_AK5_PFCHS_pt/tmp_AK5_GEN_pt;
			h1_AK5_PFCHS.Fill(ratio);

			ratio = tmp_AK5_PFCHS_pt_uncorr/tmp_AK5_GEN_pt;
			h1_AK5_PFCHS_uncorr.Fill(ratio);

			ratio = tmp_AK5_PFCHS_pt_L1_rhoSW/tmp_AK5_GEN_pt;
			h1_AK5_PFCHS_l1_rhoSW.Fill(ratio);

			ratio = tmp_AK5_PFCHS_pt_L1_rhoHand/tmp_AK5_GEN_pt;
			h1_AK5_PFCHS_l1_rhoHand.Fill(ratio);
			h1_AK5_PFCHS_Pt_l1_rhoHand.Fill(tmp_AK5_PFCHS_pt_L1_rhoHand);

			ratio = tmp_AK5_PFCHS_pt_L1_rhoHand2/tmp_AK5_GEN_pt;
			h1_AK5_PFCHS_l1_rhoHand2.Fill(ratio);

			ratio = tmp_AK5_PFCHS_pt_L1_rhoGrid/tmp_AK5_GEN_pt;
			h1_AK5_PFCHS_l1_rhoGrid.Fill(ratio);


			//pt
			h1_AK5_GEN_pt.Fill(tmp_AK5_GEN_pt);
			h1_AK5_PF_pt_uncorr.Fill(tmp_AK5_PF_pt_uncorr);// pt of PF and PFCHS are with wrong JEC now
			h1_AK5_PFCHS_pt_uncorr.Fill(tmp_AK5_PFCHS_pt_uncorr);
			//eta
			h1_AK5_GEN_eta.Fill(tmp_AK5_GEN_eta);
			h1_AK5_PF_eta.Fill(tmp_AK5_PF_eta);
			h1_AK5_PFCHS_eta.Fill(tmp_AK5_PFCHS_eta);
			//phi
			h1_AK5_GEN_phi.Fill(tmp_AK5_GEN_phi);
			h1_AK5_PF_phi.Fill(tmp_AK5_PF_phi);
			h1_AK5_PFCHS_phi.Fill(tmp_AK5_PFCHS_phi);
			//dr
			dr = TMath::Sqrt( (Z_eta-tmp_AK5_GEN_eta)*(Z_eta-tmp_AK5_GEN_eta) + (Z_phi-tmp_AK5_GEN_phi)*(Z_phi-tmp_AK5_GEN_phi) );
			h1_AK5_GEN_zjet_dr.Fill(dr);
			dr = TMath::Sqrt( (Z_eta-tmp_AK5_PF_eta)*(Z_eta-tmp_AK5_PF_eta) + (Z_phi-tmp_AK5_PF_phi)*(Z_phi-tmp_AK5_PF_phi) );
			h1_AK5_PF_zjet_dr.Fill(dr);
			dr = TMath::Sqrt( (Z_eta-tmp_AK5_PFCHS_eta)*(Z_eta-tmp_AK5_PFCHS_eta) + (Z_phi-tmp_AK5_PFCHS_phi)*(Z_phi-tmp_AK5_PFCHS_phi) );
			h1_AK5_PFCHS_zjet_dr.Fill(dr);
			//dphi
			dphi = TMath::Sqrt( (Z_phi-tmp_AK5_GEN_phi)*(Z_phi-tmp_AK5_GEN_phi) );
			h1_AK5_GEN_zjet_dphi.Fill(dphi);
			dphi = TMath::Sqrt( (Z_phi-tmp_AK5_PF_phi)*(Z_phi-tmp_AK5_PF_phi) );
			h1_AK5_PF_zjet_dphi.Fill(dphi);
			dphi = TMath::Sqrt( (Z_phi-tmp_AK5_PFCHS_phi)*(Z_phi-tmp_AK5_PFCHS_phi) );
			h1_AK5_PFCHS_zjet_dphi.Fill(dphi);
			//nPV
			h1_nPV.Fill(tmp_event_nPV);
			//z mass
			h1_z_mass.Fill(Z_mass);
			//muplus pt
			h1_muplus_Pt.Fill(Z_muplus_pt);
			//rho
			h1_AK5_GEN_rhoSW.Fill(tmp_AK5_GEN_rhoSW);
			h1_AK5_PF_rhoSW.Fill(tmp_AK5_PF_rhoSW);
			h1_AK5_PFCHS_rhoSW.Fill(tmp_AK5_PFCHS_rhoSW);

			h1_AK5_GEN_rhoHand.Fill(tmp_AK5_GEN_rhoHand);
			h1_AK5_PF_rhoHand.Fill(tmp_AK5_PF_rhoHand);
			h1_AK5_PFCHS_rhoHand.Fill(tmp_AK5_PFCHS_rhoHand);

			h1_AK5_GEN_rhoHand2.Fill(tmp_AK5_GEN_rhoHand2);
			h1_AK5_PF_rhoHand2.Fill(tmp_AK5_PF_rhoHand2);
			h1_AK5_PFCHS_rhoHand2.Fill(tmp_AK5_PFCHS_rhoHand2);

			h1_AK5_GEN_rhoGrid.Fill(tmp_AK5_GEN_rhoGrid);
			h1_AK5_PF_rhoGrid.Fill(tmp_AK5_PF_rhoGrid);
			h1_AK5_PFCHS_rhoGrid.Fill(tmp_AK5_PFCHS_rhoGrid);

			h1_PFCor_area.Fill(JetPFCor_Area[i_PFCorJet_matching_PFCHS]);
			h1_AK5_PF_area.Fill(GroomedJet_AK5_PF_area[0]);  
			h1_AK5_PFCHS_area.Fill(GroomedJet_AK5_PFCHS_area[0]); 

			/*//h1_rho_PFCor.Fill(); 
			  h1_rhoSW_PF.Fill(GroomedJet_AK5_PF_rhoSW); 
			  h1_rhoHand_PF.Fill(GroomedJet_AK5_PF_rhohand); 
			  h1_rhoHand2_PF.Fill(GroomedJet_AK5_PF_rhohand2); 
			  h1_rhoGrid_PF.Fill(GroomedJet_AK5_PF_rhogrid); 
			  h1_rhoSW_PFCHS.Fill(GroomedJet_AK5_PFCHS_rhoSW); 
			  h1_rhoHand_PFCHS.Fill(GroomedJet_AK5_PFCHS_rhohand); 
			  h1_rhoHand2_PFCHS.Fill(GroomedJet_AK5_PFCHS_rhohand2); 
			  h1_rhoGrid_PFCHS.Fill(GroomedJet_AK5_PFCHS_rhogrid); 
			  */

			//rho vs nPV
			h2_AK5_GEN_rhoSW_vs_nPV.Fill(tmp_AK5_GEN_rhoSW,tmp_event_nPV);
			h2_AK5_GEN_rhoHand_vs_nPV.Fill(tmp_AK5_GEN_rhoHand,tmp_event_nPV);
			h2_AK5_GEN_rhoHand2_vs_nPV.Fill(tmp_AK5_GEN_rhoHand2,tmp_event_nPV);
			h2_AK5_GEN_rhoGrid_vs_nPV.Fill(tmp_AK5_GEN_rhoGrid,tmp_event_nPV);

			h2_AK5_PF_rhoSW_vs_nPV.Fill(tmp_AK5_PF_rhoSW,tmp_event_nPV);
			h2_AK5_PF_rhoHand_vs_nPV.Fill(tmp_AK5_PF_rhoHand,tmp_event_nPV);
			h2_AK5_PF_rhoHand2_vs_nPV.Fill(tmp_AK5_PF_rhoHand2,tmp_event_nPV);
			h2_AK5_PF_rhoGrid_vs_nPV.Fill(tmp_AK5_PF_rhoGrid,tmp_event_nPV);
			//ratio vs pt
			//PF vs PFCHS matching efficiency
			dr = TMath::Sqrt( (gen_jet_eta-tmp_AK5_PF_eta)*(gen_jet_eta-tmp_AK5_PF_eta) + (gen_jet_phi-tmp_AK5_PF_phi)*(gen_jet_phi-tmp_AK5_PF_phi) ) ;
			h1_PF_match.Fill(dr);
			dr = TMath::Sqrt( (gen_jet_eta-tmp_AK5_PFCHS_eta)*(gen_jet_eta-tmp_AK5_PFCHS_eta) + (gen_jet_phi-tmp_AK5_PFCHS_phi)*(gen_jet_phi-tmp_AK5_PFCHS_phi) ) ;       
			h1_PFCHS_match.Fill(dr);

			tmp_AK5_GEN_mass=GenGroomedJet_AK5_GEN_mass_uncorr[0];
			tmp_PFCor_mass=JetPFCor_Mass[i_PFCorJet_matching_PFCHS];
			tmp_AK5_PF_mass_uncorr=GroomedJet_AK5_PF_mass_uncorr[0];
			tmp_AK5_PF_mass_rhoArea=GroomedJet_AK5_PF_mass_rhoArea[0];
			tmp_AK5_PF_mass_rhoGArea=GroomedJet_AK5_PF_mass_rhoGArea[0];
			tmp_AK5_PF_mass_rho4Area=GroomedJet_AK5_PF_mass_rho4Area[0];
			tmp_AK5_PF_mass_rhoG4Area=GroomedJet_AK5_PF_mass_rhoG4Area[0];
			tmp_AK5_PF_mass_rhom4Area=GroomedJet_AK5_PF_mass_rhom4Area[0];
			tmp_AK5_PF_mass_JetCleansingATLASjvf=GroomedJet_AK5_PF_mass_JetCleansingATLASjvf[0];
			tmp_AK5_PF_mass_JetCleansingATLASlin=GroomedJet_AK5_PF_mass_JetCleansingATLASlin[0];
			tmp_AK5_PF_mass_JetCleansingATLASgau=GroomedJet_AK5_PF_mass_JetCleansingATLASgau[0];
			tmp_AK5_PF_mass_JetCleansingCMSjvf=GroomedJet_AK5_PF_mass_JetCleansingCMSjvf[0];
			tmp_AK5_PF_mass_JetCleansingCMSlin=GroomedJet_AK5_PF_mass_JetCleansingCMSlin[0];
			tmp_AK5_PF_mass_JetCleansingCMSgau=GroomedJet_AK5_PF_mass_JetCleansingCMSgau[0];
			tmp_AK5_PFCHS_mass_uncorr=GroomedJet_AK5_PFCHS_mass_uncorr[0];
			tmp_AK5_PFCHS_mass_rhoArea=GroomedJet_AK5_PFCHS_mass_rhoArea[0];
			tmp_AK5_PFCHS_mass_rhoGArea=GroomedJet_AK5_PFCHS_mass_rhoG4Area[0];
			tmp_AK5_PFCHS_mass_rho4Area=GroomedJet_AK5_PFCHS_mass_rho4Area[0];
			tmp_AK5_PFCHS_mass_rhoG4Area=GroomedJet_AK5_PFCHS_mass_rhoG4Area[0];
			tmp_AK5_PFCHS_mass_rhom4Area=GroomedJet_AK5_PFCHS_mass_rhom4Area[0];

			h1_AK5_GEN_mass.Fill(tmp_AK5_GEN_mass            ); 
			h1_PFCor_mass.Fill(tmp_PFCor_mass          );
			h1_AK5_PF_mass_uncorr.Fill(tmp_AK5_PF_mass_uncorr  );
			h1_AK5_PF_mass_rhoArea.Fill(tmp_AK5_PF_mass_rhoArea );
			h1_AK5_PF_mass_rhoGArea.Fill(tmp_AK5_PF_mass_rhoGArea);
			h1_AK5_PF_mass_rho4Area.Fill(tmp_AK5_PF_mass_rho4Area);
			h1_AK5_PF_mass_rhoG4Area.Fill(tmp_AK5_PF_mass_rhoG4Area);
			h1_AK5_PF_mass_rhom4Area.Fill(tmp_AK5_PF_mass_rhom4Area);
			h1_AK5_PF_mass_JetCleansingATLASjvf.Fill(tmp_AK5_PF_mass_JetCleansingATLASjvf);
			h1_AK5_PF_mass_JetCleansingATLASlin.Fill(tmp_AK5_PF_mass_JetCleansingATLASlin);
			h1_AK5_PF_mass_JetCleansingATLASgau.Fill(tmp_AK5_PF_mass_JetCleansingATLASgau);
			h1_AK5_PF_mass_JetCleansingCMSjvf.Fill(tmp_AK5_PF_mass_JetCleansingCMSjvf);
			h1_AK5_PF_mass_JetCleansingCMSlin.Fill(tmp_AK5_PF_mass_JetCleansingCMSlin);
			h1_AK5_PF_mass_JetCleansingCMSgau.Fill(tmp_AK5_PF_mass_JetCleansingCMSgau);

			h1_AK5_PFCHS_mass_uncorr.Fill(tmp_AK5_PFCHS_mass_uncorr);  
			h1_AK5_PFCHS_mass_rhoArea.Fill(tmp_AK5_PFCHS_mass_rhoArea); 
			h1_AK5_PFCHS_mass_rhoGArea.Fill(tmp_AK5_PFCHS_mass_rhoGArea);
			h1_AK5_PFCHS_mass_rho4Area.Fill(tmp_AK5_PFCHS_mass_rho4Area);
			h1_AK5_PFCHS_mass_rhoG4Area.Fill(tmp_AK5_PFCHS_mass_rhoG4Area);
			h1_AK5_PFCHS_mass_rhom4Area.Fill(tmp_AK5_PFCHS_mass_rhom4Area);

		}
	}
	TCanvas *c1 = new TCanvas("c1","RECO_vs_GEN Pt",200,10,600,600);
	c1->cd();
	h1_AK5_GEN_pt.Draw();
	h1_AK5_PF_pt_uncorr.Draw("same");
	c1->Print("RECO vs_GEN_Pt.png");

	TCanvas *c2 = new TCanvas("c2","RECO_vs_GEN_Eta",200,10,600,600);
	c2->cd();
	h1_AK5_GEN_eta.Draw();
	h1_AK5_PF_eta.Draw("same");
	c2->Print("RECO_vs_GEN_Eta.png");

	TCanvas *c13 = new TCanvas("c13","RECO_vs_GEN_Phi",200,10,600,600);
	c13->cd();
	h1_AK5_GEN_phi.Draw();
	h1_AK5_PF_phi.Draw("same");
	c13->Print("RECO_vs_GEN_Phi.png");

	TCanvas *c3 = new TCanvas("c3","RECO_vs_GEN_dR",200,10,600,600);
	c3->cd();
	h1_AK5_GEN_zjet_dr.Draw();
	h1_AK5_PF_zjet_dr.Draw("same");
	c3->Print("RECO_vs_GEN_dR.png");

	TCanvas *c4 = new TCanvas("c4","RECO_vs_GEN_dphi",200,10,600,600);
	c4->cd();
	h1_AK5_GEN_zjet_dphi.Draw();
	h1_AK5_PF_zjet_dphi.Draw("same");
	c4->Print("RECO_vs_GEN_dphi.png");

	TCanvas *c5 = new TCanvas("c5","RECO_vs_GEN_rhoSW",200,10,600,600);
	c5->cd();
	h1_AK5_GEN_rhoSW.Draw();
	h1_AK5_PF_rhoSW.Draw("same");
	c5->Print("RECO_vs_GEN_rhoSW.png");

	TCanvas *c6 = new TCanvas("c6","RECO_vs_GEN_rhoHand",200,10,600,600);
	c6->cd();
	h1_AK5_GEN_rhoHand.Draw();
	h1_AK5_PF_rhoHand.Draw("same");
	c6->Print("RECO_vs_GEN_rhoHand.png");

	TCanvas *c7 = new TCanvas("c7","RECO_vs_GEN_rhoHand2",200,10,600,600);
	c7->cd();
	h1_AK5_GEN_rhoHand2.Draw();
	h1_AK5_PF_rhoHand2.Draw("same");
	c7->Print("RECO_vs_GEN_rhoHand2.png");

	TCanvas *c8 = new TCanvas("c8","RECO_vs_GEN_rhoGrid",200,10,600,600);
	c8->cd();
	h1_AK5_GEN_rhoGrid.Draw();
	h1_AK5_PF_rhoGrid.Draw("same");
	c8->Print("RECO_vs_GEN_rhoGrid.png");

	TCanvas *c9 = new TCanvas("c9","RECO_vs_GEN_rhoSW_vs_nPV",200,10,600,600);
	c9->cd();
	h2_AK5_GEN_rhoSW_vs_nPV.Draw();
	h2_AK5_PF_rhoSW_vs_nPV.Draw("same");
	c9->Print("RECO_vs_GEN_rhoSW_vs_nPV.png");

	TCanvas *c10 = new TCanvas("c10","RECO_vs_GEN_rhoHand_vs_nPV",200,10,600,600);
	c10->cd();
	h2_AK5_GEN_rhoHand_vs_nPV.Draw();
	h2_AK5_PF_rhoHand_vs_nPV.Draw("same");
	c10->Print("RECO_vs_GEN_rhoHand_vs_nPV.png");

	TCanvas *c11 = new TCanvas("c11","RECO_vs_GEN_rhoHand2_vs_nPV",200,10,600,600);
	c11->cd();
	h2_AK5_GEN_rhoHand2_vs_nPV.Draw();
	h2_AK5_PF_rhoHand2_vs_nPV.Draw("same");
	c11->Print("RECO_vs_GEN_rhoHand2_vs_nPV.png");

	TCanvas *c12 = new TCanvas("c12","RECO_vs_GEN_rhoGrid_vs_nPV",200,10,600,600);
	c12->cd();
	h2_AK5_GEN_rhoGrid_vs_nPV.Draw();
	h2_AK5_PF_rhoGrid_vs_nPV.Draw("same");
	c12->Print("RECO_vs_GEN_rhoGrid_vs_nPV.png");

	TCanvas *c14 = new TCanvas("c14","f",200,10,600,600);
	c14->cd();
	h1_PF_match.Draw();
	h1_PFCHS_match.Draw("same");
	c14->Print("PF_vs_PFCHS_matching_efficiency.png");
	//ratio

	Draw_and_Save(h1_PFCor);
	Draw_and_Save(h1_PFCor_uncorr);
	Draw_and_Save(h1_PFCor_afterL1);
	Draw_and_Save(h1_PFCor_Pt_afterL1);
	Draw_and_Save(h1_PFCor_afterL2);
	Draw_and_Save(h1_AK5_PF);
	Draw_and_Save(h1_AK5_PF_uncorr);
	Draw_and_Save(h1_AK5_PF_l1_rhoSW);
	Draw_and_Save(h1_AK5_PF_l1_rhoHand);
	Draw_and_Save(h1_AK5_PF_Pt_l1_rhoHand);
	Draw_and_Save(h1_AK5_PF_l1_rhoHand2);
	Draw_and_Save(h1_AK5_PF_l1_rhoGrid);
	Draw_and_Save(h1_AK5_PFCHS);
	Draw_and_Save(h1_AK5_PFCHS_uncorr);
	Draw_and_Save(h1_AK5_PFCHS_l1_rhoSW);
	Draw_and_Save(h1_AK5_PFCHS_l1_rhoHand);
	Draw_and_Save(h1_AK5_PFCHS_Pt_l1_rhoHand);
	Draw_and_Save(h1_AK5_PFCHS_l1_rhoHand2);
	Draw_and_Save(h1_AK5_PFCHS_l1_rhoGrid);

	//compare PFCor, AK5PF, AK5PFCHS
	h1_PFCor_Pt_afterL1.SetLineColor(kBlue);
	h1_AK5_PF_Pt_l1_rhoHand.SetLineColor(kBlack); h1_AK5_PF_Pt_l1_rhoHand.SetLineStyle(2);
	h1_AK5_PFCHS_Pt_l1_rhoHand.SetLineColor(kRed); h1_AK5_PFCHS_Pt_l1_rhoHand.SetLineStyle(2);
	Draw_and_Save(h1_PFCor_Pt_afterL1, h1_AK5_PF_Pt_l1_rhoHand, h1_AK5_PFCHS_Pt_l1_rhoHand);
	//compare PFCor, AK5PF, AK5PFCHS
	h1_PFCor_afterL1.SetLineColor(kBlue);
	h1_AK5_PF_l1_rhoHand.SetLineColor(kBlack); h1_AK5_PF_l1_rhoHand.SetLineStyle(2);
	h1_AK5_PFCHS_l1_rhoHand.SetLineColor(kRed); h1_AK5_PFCHS_l1_rhoHand.SetLineStyle(2);
	Draw_and_Save(h1_PFCor_afterL1, h1_AK5_PF_l1_rhoHand, h1_AK5_PFCHS_l1_rhoHand);

	Draw_and_Save(h1_PFCor_area);
	Draw_and_Save(h1_AK5_PF_area);
	Draw_and_Save(h1_AK5_PFCHS_area);

	h1_PFCor_area.SetLineColor(kBlue);
	h1_AK5_PFCHS_area.SetLineColor(kBlack);
	h1_AK5_PFCHS_area.SetLineStyle(2);
	Draw_and_Save(h1_PFCor_area, h1_AK5_PFCHS_area);

	//h1_rho_PFCor); 
	/*h1_rhoSW_PF); 
	  h1_rhoHand_PF);
	  h1_rhoHand2_PF);
	  h1_rhoGrid_PF);
	  h1_rhoSW_PFCHS);
	  h1_rhoHand_PFCHS);
	  h1_rhoHand2_PFCHS);
	  h1_rhoGrid_PFCHS);*/


	Draw_and_Save(h1_nPV);
	Draw_and_Save(h1_z_mass);
	Draw_and_Save(h1_muplus_Pt);

	Draw_and_Save(h1_AK5_GEN_pt);
	Draw_and_Save(h1_AK5_GEN_eta);
	Draw_and_Save(h1_AK5_GEN_zjet_dr);
	Draw_and_Save(h1_AK5_GEN_zjet_dphi);
	Draw_and_Save(h1_AK5_GEN_rhoSW);
	Draw_and_Save(h1_AK5_GEN_rhoHand);
	Draw_and_Save(h1_AK5_GEN_rhoHand2);
	Draw_and_Save(h1_AK5_GEN_rhoGrid);
	Draw_and_Save(h2_AK5_GEN_rhoSW_vs_nPV);
	Draw_and_Save(h2_AK5_GEN_rhoHand_vs_nPV);
	Draw_and_Save(h2_AK5_GEN_rhoHand2_vs_nPV);
	Draw_and_Save(h2_AK5_GEN_rhoGrid_vs_nPV);

	Draw_and_Save(h1_AK5_PF_pt_uncorr);
	Draw_and_Save(h1_AK5_PF_eta);
	Draw_and_Save(h1_AK5_PF_zjet_dr);
	Draw_and_Save(h1_AK5_PF_zjet_dphi);
	Draw_and_Save(h1_AK5_PF_rhoSW);
	Draw_and_Save(h1_AK5_PF_rhoHand);
	Draw_and_Save(h1_AK5_PF_rhoHand2);
	Draw_and_Save(h1_AK5_PF_rhoGrid);
	Draw_and_Save(h2_AK5_PF_rhoSW_vs_nPV);
	Draw_and_Save(h2_AK5_PF_rhoHand_vs_nPV);
	Draw_and_Save(h2_AK5_PF_rhoHand2_vs_nPV);
	Draw_and_Save(h2_AK5_PF_rhoGrid_vs_nPV);


	h1_AK5_PF_rhoSW.SetLineColor(kBlue);
	h1_AK5_PF_rhoHand.SetLineColor(kBlack);
	h1_AK5_PF_rhoHand.SetLineStyle(2);
	Draw_and_Save(h1_AK5_PF_rhoSW, h1_AK5_PF_rhoHand);



	Draw_and_Save(h1_AK5_PFCHS_pt_uncorr);
	Draw_and_Save(h1_AK5_PFCHS_eta);
	Draw_and_Save(h1_AK5_PFCHS_zjet_dr);
	Draw_and_Save(h1_AK5_PFCHS_zjet_dphi);
	Draw_and_Save(h1_AK5_PFCHS_rhoSW);
	Draw_and_Save(h1_AK5_PFCHS_rhoHand);
	Draw_and_Save(h1_AK5_PFCHS_rhoHand2);
	Draw_and_Save(h1_AK5_PFCHS_rhoGrid);
	Draw_and_Save(h2_AK5_PFCHS_rhoSW_vs_nPV);
	Draw_and_Save(h2_AK5_PFCHS_rhoHand_vs_nPV);
	Draw_and_Save(h2_AK5_PFCHS_rhoHand2_vs_nPV);
	Draw_and_Save(h2_AK5_PFCHS_rhoGrid_vs_nPV);


	Draw_and_Save(h2_AK5_PFCHS_ratioHand_vs_nPV);
	Draw_and_Save(h2_AK5_PFCHS_ratioHand_vs_ptHand);
	Draw_and_Save(h2_AK5_PFCHS_ratioHand_vs_eta);

	Draw_and_Save(h1_PF_match);
	Draw_and_Save(h1_PFCHS_match);


	Draw_and_Save(h1_AK5_GEN_mass                );
	Draw_and_Save(h1_PFCor_mass              );
	Draw_and_Save(h1_AK5_PF_mass_uncorr      );
	Draw_and_Save(h1_AK5_PF_mass_rhoArea     );
	Draw_and_Save(h1_AK5_PF_mass_rhoGArea    );
	Draw_and_Save(h1_AK5_PF_mass_rho4Area    );
	Draw_and_Save(h1_AK5_PF_mass_rhoG4Area   );
	Draw_and_Save(h1_AK5_PF_mass_rhom4Area   );
	Draw_and_Save(h1_AK5_PF_mass_JetCleansingATLASjvf   );
	Draw_and_Save(h1_AK5_PF_mass_JetCleansingATLASlin   );
	Draw_and_Save(h1_AK5_PF_mass_JetCleansingATLASgau   );
	Draw_and_Save(h1_AK5_PF_mass_JetCleansingCMSjvf   );
	Draw_and_Save(h1_AK5_PF_mass_JetCleansingCMSlin   );
	Draw_and_Save(h1_AK5_PF_mass_JetCleansingCMSgau   );
	Draw_and_Save(h1_AK5_PFCHS_mass_uncorr   );
	Draw_and_Save(h1_AK5_PFCHS_mass_rhoArea  );
	Draw_and_Save(h1_AK5_PFCHS_mass_rhoGArea );
	Draw_and_Save(h1_AK5_PFCHS_mass_rho4Area );
	Draw_and_Save(h1_AK5_PFCHS_mass_rhoG4Area);
	Draw_and_Save(h1_AK5_PFCHS_mass_rhom4Area);


	TCanvas *c15 = new TCanvas("c15","f",200,10,600,600);
	c15->cd();

	h1_AK5_GEN_mass.SetLineColor(2);
	h1_AK5_GEN_mass.Draw();
	//h1_PFCor_mass.Draw("same");
	h1_AK5_PF_mass_uncorr.SetLineColor(3);
	h1_AK5_PF_mass_uncorr.Draw("same");
	//h1_AK5_PF_mass_rhoArea.Draw("same");
	//h1_AK5_PF_mass_rhoGArea.Draw("same");
	h1_AK5_PF_mass_rho4Area.SetLineColor(4);
	h1_AK5_PF_mass_rho4Area.Draw("same");
	//h1_AK5_PF_mass_rhoG4Area.Draw("same");
	h1_AK5_PF_mass_rhom4Area.SetLineColor(5);
	h1_AK5_PF_mass_rhom4Area.Draw("same");
	//h1_AK5_PF_mass_JetCleansingATLASjvf.Draw("same");
	//h1_AK5_PF_mass_JetCleansingATLASlin.Draw("same");
	//h1_AK5_PF_mass_JetCleansingATLASgau.Draw("same");
	//h1_AK5_PF_mass_JetCleansingCMSjvf.Draw("same");
	h1_AK5_PF_mass_JetCleansingCMSlin.SetLineColor(6);
	h1_AK5_PF_mass_JetCleansingCMSlin.Draw("same");
	//h1_AK5_PF_mass_JetCleansingCMSgau.Draw("same");
	c15->Print("GEN_vs_PFCor_vs_PF_jetmass.png");


	TCanvas *c16 = new TCanvas("c16","f",200,10,600,600);
	c16->cd();
	h1_AK5_GEN_mass.Draw();
	//h1_PFCor_mass.Draw("same");
	h1_AK5_PFCHS_mass_uncorr.Draw("same");
	//h1_AK5_PFCHS_mass_rhoArea.Draw("same");
	//h1_AK5_PFCHS_mass_rhoGArea.Draw("same");
	h1_AK5_PFCHS_mass_rho4Area.Draw("same");
	//h1_AK5_PFCHS_mass_rhoG4Area.Draw("same");
	h1_AK5_PFCHS_mass_rhom4Area.Draw("same");
	c16->Print("GEN_vs_PFCor_vs_PFCHS_jetmass.png");

	c1->Write();
	c2->Write();
	c3->Write();
	c4->Write();
	c5->Write();
	c6->Write();
	c7->Write();
	c8->Write();
	c9->Write();
	c10->Write();
	c11->Write();
	c12->Write();
	c13->Write();
	c14->Write();
	c15->Write();
	c16->Write();

}


void MyClass::LoopAK8()
{
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();

	Double_t tmp_AK8_GEN_pt=0.;
	Double_t tmp_AK8_GEN_eta=0.;
	Double_t tmp_AK8_GEN_phi=0.;
	Double_t tmp_AK8_GEN_rhoSW=0.;
	Double_t tmp_AK8_GEN_rhoHand=0.;
	Double_t tmp_AK8_GEN_rhoHand2=0.;
	Double_t tmp_AK8_GEN_rhoGrid=0.;

	Double_t tmp_AK8_PF_pt=0.;
	Double_t tmp_AK8_PF_eta=0.;
	Double_t tmp_AK8_PF_phi=0.;
	Double_t tmp_AK8_PF_pt_uncorr=0.;
	Double_t tmp_AK8_PF_pt_L1_rhoSW=0.;
	Double_t tmp_AK8_PF_pt_L1_rhoHand=0.;
	Double_t tmp_AK8_PF_pt_L1_rhoHand2=0.;
	Double_t tmp_AK8_PF_pt_L1_rhoGrid=0.;
	Double_t tmp_AK8_PF_pt_rho4A=0.;
	Double_t tmp_AK8_PF_pt_rhom4A=0.;
	Double_t tmp_AK8_PF_pt_JetCleansing=0.;
	Double_t tmp_AK8_PF_rhoSW=0.;
	Double_t tmp_AK8_PF_rhoHand=0.;
	Double_t tmp_AK8_PF_rhoHand2=0.;
	Double_t tmp_AK8_PF_rhoGrid=0.;

	Double_t tmp_event_nPV=0.;

	Double_t ratio=0.;
	Double_t dr=0.; // Delta R
	Double_t dphi=0.; // Delta Phi

	//jet mass
	Double_t tmp_AK8_GEN_mass=0.;
	Double_t tmp_AK8_PF_mass_uncorr=0.;
	Double_t tmp_AK8_PF_mass_rhoArea=0.;
	Double_t tmp_AK8_PF_mass_rhoGArea=0.;
	Double_t tmp_AK8_PF_mass_rho4Area=0.;
	Double_t tmp_AK8_PF_mass_rhoG4Area=0.;
	Double_t tmp_AK8_PF_mass_rhom4Area=0.;
	Double_t tmp_AK8_PF_mass_JetCleansingATLASjvf=0.;
	Double_t tmp_AK8_PF_mass_JetCleansingATLASlin=0.;
	Double_t tmp_AK8_PF_mass_JetCleansingATLASgau=0.;
	Double_t tmp_AK8_PF_mass_JetCleansingCMSjvf=0.;
	Double_t tmp_AK8_PF_mass_JetCleansingCMSlin=0.;
	Double_t tmp_AK8_PF_mass_JetCleansingCMSgau=0.;


	Int_t nbin_rho=50; Double_t rhomin=0.; Double_t rhomax=50.;
	Int_t nbin_mass=60;Double_t jetmass_min=0;Double_t jetmass_max=300.;
	Int_t nbin_pt=40;Double_t jetpt_min=50;Double_t jetpt_max=450.;

	TH1D h1_nPV("h1_nPV","h1_nPV;nPV",50,0,50);
	TH1D h1_z_mass("h1_z_mass","h1_z_mass;Z mass",60,60,120);
	TH1D h1_z_Pt("h1_z_Pt","h1_z_Pt;Z Pt",60,60,360);
	TH1D h1_muplus_Pt("h1_muplus_Pt","h1_muplus_Pt;#mu^{+} Pt",40,0,400);

	TH1D h1_AK8_GEN_pt("h1_AK8_GEN_pt","h1_AK8_GEN_pt;Jet Pt",nbin_pt, jetpt_min, jetpt_max); h1_AK8_GEN_pt.SetLineColor(kRed);
	TH1D h1_AK8_GEN_eta("h1_AK8_GEN_eta","h1_AK8_GEN_eta",50,-2.5,2.5); h1_AK8_GEN_eta.SetLineColor(kRed);
	TH1D h1_AK8_GEN_phi("h1_AK8_GEN_phi","h1_AK8_GEN_phi",50,-4,4); h1_AK8_GEN_phi.SetLineColor(kRed);
	TH1D h1_AK8_GEN_zjet_dr("h1_AK8_GEN_zjet_dr","h1_AK8_GEN_zjet_dr;dR(Z,J)",50,0,10); h1_AK8_GEN_zjet_dr.SetLineColor(kRed);
	TH1D h1_AK8_GEN_zjet_dphi("h1_AK8_GEN_zjet_dphi","h1_AK8_GEN_zjet_dphi",50,0,5); h1_AK8_GEN_zjet_dphi.SetLineColor(kRed);

	TH1D h1_AK8_GEN_rhoSW("h1_AK8_GEN_rhoSW","h1_AK8_GEN_rhoSW",nbin_rho,rhomin,rhomax); h1_AK8_GEN_rhoSW.SetLineColor(kRed);
	TH1D h1_AK8_GEN_rhoHand("h1_AK8_GEN_rhoHand","h1_AK8_GEN_rhoHand",nbin_rho,rhomin,rhomax); h1_AK8_GEN_rhoHand.SetLineColor(kRed);
	TH1D h1_AK8_GEN_rhoHand2("h1_AK8_GEN_rhoHand2","h1_AK8_GEN_rhoHand2",nbin_rho,rhomin,rhomax); h1_AK8_GEN_rhoHand2.SetLineColor(kRed);
	TH1D h1_AK8_GEN_rhoGrid("h1_AK8_GEN_rhoGrid","h1_AK8_GEN_rhoGrid",nbin_rho,rhomin,rhomax); h1_AK8_GEN_rhoGrid.SetLineColor(kRed);
	TH2D h2_AK8_GEN_rhoSW_vs_nPV("h2_AK8_GEN_rhoSW_vs_nPV","h2_AK8_GEN_rhoSW_vs_nPV",nbin_rho,rhomin,rhomax,50,0,50); h2_AK8_GEN_rhoSW_vs_nPV.SetMarkerColor(kRed);
	TH2D h2_AK8_GEN_rhoHand_vs_nPV("h2_AK8_GEN_rhoHand_vs_nPV","h2_AK8_GEN_rhoHand_vs_nPV",nbin_rho,rhomin,rhomax,50,0,50); h2_AK8_GEN_rhoHand_vs_nPV.SetMarkerColor(kRed);
	TH2D h2_AK8_GEN_rhoHand2_vs_nPV("h2_AK8_GEN_rhoHand2_vs_nPV","h2_AK8_GEN_rhoHand2_vs_nPV",nbin_rho,rhomin,rhomax,50,0,50); h2_AK8_GEN_rhoHand2_vs_nPV.SetMarkerColor(kRed);
	TH2D h2_AK8_GEN_rhoGrid_vs_nPV("h2_AK8_GEN_rhoGrid_vs_nPV","h2_AK8_GEN_rhoGrid_vs_nPV",nbin_rho,rhomin,rhomax,50,0,50); h2_AK8_GEN_rhoGrid_vs_nPV.SetMarkerColor(kRed);

	TH1D h1_AK8_PF_pt_uncorr("h1_AK8_PF_pt_uncorr","h1_AK8_PF_pt_uncorr;Jet Pt",nbin_pt, jetpt_min, jetpt_max);
	TH1D h1_AK8_PF_pt_l1rhoHand("h1_AK8_PF_pt_l1rhoHand","h1_AK8_PF_pt_l1rhoHand",nbin_pt, jetpt_min, jetpt_max);
	TH1D h1_AK8_PF_pt_l1rhoGrid("h1_AK8_PF_pt_l1rhoGrid","h1_AK8_PF_pt_l1rhoGrid",nbin_pt, jetpt_min, jetpt_max);
	TH1D h1_AK8_PF_pt_rho4A("h1_AK8_PF_pt_rho4A","h1_AK8_PF_pt_rho4A",nbin_pt, jetpt_min, jetpt_max);
	TH1D h1_AK8_PF_pt_rhom4A("h1_AK8_PF_pt_rhom4A","h1_AK8_PF_pt_rhom4A",nbin_pt, jetpt_min, jetpt_max);
	TH1D h1_AK8_PF_pt_JetCleansing("h1_AK8_PF_pt_JetCleansing","h1_AK8_PF_pt_JetCleansing",nbin_pt, jetpt_min, jetpt_max);
	TH1D h1_AK8_PF_eta("h1_AK8_PF_eta","h1_AK8_PF_eta",50,-2.5,2.5);
	TH1D h1_AK8_PF_phi("h1_AK8_PF_phi","h1_AK8_PF_phi",50,-4,4);
	TH1D h1_AK8_PF_zjet_dr("h1_AK8_PF_zjet_dr","h1_AK8_PF_zjet_dr;dR(Z,J)",50,0,7);
	TH1D h1_AK8_PF_zjet_dphi("h1_AK8_PF_zjet_dphi","h1_AK8_PF_zjet_dphi",50,0,5);

	TH1D h1_AK8_PF_rhoSW("h1_AK8_PF_rhoSW","h1_AK8_PF_rhoSW",nbin_rho,rhomin,rhomax);
	TH1D h1_AK8_PF_rhoHand("h1_AK8_PF_rhoHand","h1_AK8_PF_rhoHand",nbin_rho,rhomin,rhomax);
	TH1D h1_AK8_PF_rhoHand2("h1_AK8_PF_rhoHand2","h1_AK8_PF_rhoHand2",nbin_rho,rhomin,rhomax);
	TH1D h1_AK8_PF_rhoGrid("h1_AK8_PF_rhoGrid","h1_AK8_PF_rhoGrid",nbin_rho,rhomin,rhomax);
	TH2D h2_AK8_PF_rhoSW_vs_nPV("h2_AK8_PF_rhoSW_vs_nPV","h2_AK8_PF_rhoSW_vs_nPV",nbin_rho,rhomin,rhomax,50,0,50);
	TH2D h2_AK8_PF_rhoHand_vs_nPV("h2_AK8_PF_rhoHand_vs_nPV","h2_AK8_PF_rhoHand_vs_nPV",nbin_rho,rhomin,rhomax,50,0,50);
	TH2D h2_AK8_PF_rhoHand2_vs_nPV("h2_AK8_PF_rhoHand2_vs_nPV","h2_AK8_PF_rhoHand2_vs_nPV",nbin_rho,rhomin,rhomax,50,0,50);
	TH2D h2_AK8_PF_rhoGrid_vs_nPV("h2_AK8_PF_rhoGrid_vs_nPV","h2_AK8_PF_rhoGrid_vs_nPV",nbin_rho,rhomin,rhomax,50,0,50);

	TH1D h1_PF_match("h1_PF_match","h1_PF_match",50,0,1.); h1_PF_match.SetLineColor(kRed);//matching with GEN

	// JetPt/GenPt
	Int_t nbin_ratio=20; Double_t ratio_min=0.3; Double_t ratio_max=1.7; 
	TH1D h1_AK8_PF_recogenptratio("h1_AK8_PF_recogenptratio","h1_AK8_PF_recogenptratio",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_AK8_PF_uncorr_recogenptratio("h1_AK8_PF_uncorr_recogenptratio","h1_AK8_PF_uncorr_recogenptratio;Reco/Gen Pt ratio ",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_AK8_PF_l1rhosw_recogenptratio("h1_AK8_PF_l1rhosw_recogenptratio","h1_AK8_PF_l1rhosw_recogenptratio",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_AK8_PF_l1rhoHand_recogenptratio("h1_AK8_PF_l1rhoHand_recogenptratio","h1_AK8_PF_l1rhoHand_recogenptratio",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_AK8_PF_l1rhoHand2_recogenptratio("h1_AK8_PF_l1rhoHand2_recogenptratio","h1_AK8_PF_l1rhoHand2_recogenptratio",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_AK8_PF_l1rhoGrid_recogenptratio("h1_AK8_PF_l1rhoGrid_recogenptratio","h1_AK8_PF_l1rhoGrid_recogenptratio",nbin_ratio, ratio_min, ratio_max);

	TH1D h1_AK8_PF_rho4A_recogenptratio("h1_AK8_PF_rho4A_recogenptratio","h1_AK8_PF_rho4A_recogenptratio",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_AK8_PF_rhom4A_recogenptratio("h1_AK8_PF_rhom4A_recogenptratio","h1_AK8_PF_rhom4A_recogenptratio",nbin_ratio, ratio_min, ratio_max);
	TH1D h1_AK8_PF_JetCleansing_recogenptratio("h1_AK8_PF_JetCleansing_recogenptratio","h1_AK8_PF_JetCleansing_recogenptratio",nbin_ratio, ratio_min, ratio_max);

	//area
	TH1D h1_AK8_PF_area("h1_AK8_PF_area","h1_AK8_PF_area",50,0.6,1.1);

	//mass
	TH1D h1_AK8_GEN_mass("h1_AK8_GEN_mass","h1_AK8_GEN_mass;jet mass;",nbin_mass,jetmass_min,jetmass_max); h1_AK8_GEN_mass.SetLineColor(kRed);
	TH1D h1_AK8_PF_mass_uncorr("h1_AK8_PF_mass_uncorr","h1_AK8_PF_mass_uncorr;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_rhoArea("h1_AK8_PF_mass_rhoArea","h1_AK8_PF_mass_rhoArea;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_rhoGArea("h1_AK8_PF_mass_rhoGArea","h1_AK8_PF_mass_rhoGArea;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_rho4Area("h1_AK8_PF_mass_rho4Area","h1_AK8_PF_mass_rho4Area;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_rhoG4Area("h1_AK8_PF_mass_rhoG4Area","h1_AK8_PF_mass_rhoG4Area;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_rhom4Area("h1_AK8_PF_mass_rhom4Area","h1_AK8_PF_mass_rhom4Area;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_JetCleansingATLASjvf("h1_AK8_PF_mass_JetCleansingATLASjvf","h1_AK8_PF_mass_JetCleansingATLASjvf;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_JetCleansingATLASlin("h1_AK8_PF_mass_JetCleansingATLASlin","h1_AK8_PF_mass_JetCleansingATLASlin;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_JetCleansingATLASgau("h1_AK8_PF_mass_JetCleansingATLASgau","h1_AK8_PF_mass_JetCleansingATLASgau;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_JetCleansingCMSjvf("h1_AK8_PF_mass_JetCleansingCMSjvf","h1_AK8_PF_mass_JetCleansingCMSjvf;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_JetCleansingCMSlin("h1_AK8_PF_mass_JetCleansingCMSlin","h1_AK8_PF_mass_JetCleansingCMSlin;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_JetCleansingCMSgau("h1_AK8_PF_mass_JetCleansingCMSgau","h1_AK8_PF_mass_JetCleansingCMSgau;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	// reco mass VS gen mass
	TH2D h2_AK8_GEN_mass("h2_AK8_GEN_mass","h2_AK8_GEN_mass; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max); 
	TH2D h2_AK8_PF_mass_uncorr("h2_AK8_PF_mass_uncorr","h2_AK8_PF_mass_uncorr; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_rhoArea("h2_AK8_PF_mass_rhoArea","h2_AK8_PF_mass_rhoArea; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_rhoGArea("h2_AK8_PF_mass_rhoGArea","h2_AK8_PF_mass_rhoGArea; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_rho4Area("h2_AK8_PF_mass_rho4Area","h2_AK8_PF_mass_rho4Area; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_rhoG4Area("h2_AK8_PF_mass_rhoG4Area","h2_AK8_PF_mass_rhoG4Area; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_rhom4Area("h2_AK8_PF_mass_rhom4Area","h2_AK8_PF_mass_rhom4Area; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_JetCleansingATLASjvf("h2_AK8_PF_mass_JetCleansingATLASjvf","h2_AK8_PF_mass_JetCleansingATLASjvf; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_JetCleansingATLASlin("h2_AK8_PF_mass_JetCleansingATLASlin","h2_AK8_PF_mass_JetCleansingATLASlin; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_JetCleansingATLASgau("h2_AK8_PF_mass_JetCleansingATLASgau","h2_AK8_PF_mass_JetCleansingATLASgau; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_JetCleansingCMSjvf("h2_AK8_PF_mass_JetCleansingCMSjvf","h2_AK8_PF_mass_JetCleansingCMSjvf; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_JetCleansingCMSlin("h2_AK8_PF_mass_JetCleansingCMSlin","h2_AK8_PF_mass_JetCleansingCMSlin; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_JetCleansingCMSgau("h2_AK8_PF_mass_JetCleansingCMSgau","h2_AK8_PF_mass_JetCleansingCMSgau; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);

	//calculate correlationFactor;
	Int_t num_points=0;
	TGraph gr_AK8_GEN_mass;
	TGraph gr_AK8_PF_mass_uncorr;
	TGraph gr_AK8_PF_mass_rhoArea;
	TGraph gr_AK8_PF_mass_rhoGArea;
	TGraph gr_AK8_PF_mass_rho4Area;
	TGraph gr_AK8_PF_mass_rhoG4Area;
	TGraph gr_AK8_PF_mass_rhom4Area;
	TGraph gr_AK8_PF_mass_JetCleansingATLASjvf;
	TGraph gr_AK8_PF_mass_JetCleansingATLASlin;
	TGraph gr_AK8_PF_mass_JetCleansingATLASgau;
	TGraph gr_AK8_PF_mass_JetCleansingCMSjvf;
	TGraph gr_AK8_PF_mass_JetCleansingCMSlin;
	TGraph gr_AK8_PF_mass_JetCleansingCMSgau;


	std::vector<TH1D> vect_h1_AK8_PF_mass_JetCleansing_DiffMode;
	std::vector<TH1D> vect_h1_AK8_PF_pt_JetCleansing_DiffMode;
	Int_t number_JetCleansing_DiffMode=0;
	for(Int_t i=0;i<50;i++){
		TH1D h1_AK8_PF_mass_JetCleansing_DiffMode(Form("h1_AK8_PF_mass_JetCleansing_DiffMode%i",i),Form("h1_AK8_PF_mass_JetCleansing_DiffMode%i;jet mass;",i),nbin_mass,jetmass_min,jetmass_max);
		TH1D h1_AK8_PF_pt_JetCleansing_DiffMode(Form("h1_AK8_PF_pt_JetCleansing_DiffMode%i",i),Form("h1_AK8_PF_pt_JetCleansing_DiffMode%i;jet pT",i),nbin_pt, jetpt_min, jetpt_max);
		vect_h1_AK8_PF_mass_JetCleansing_DiffMode.push_back(h1_AK8_PF_mass_JetCleansing_DiffMode);
		vect_h1_AK8_PF_pt_JetCleansing_DiffMode.push_back(h1_AK8_PF_pt_JetCleansing_DiffMode);
	}


	std::vector<TH2D> vect_h2_AK8_PF_mass_JetCleansing_DiffMode;
	std::vector<TGraph> vect_gr_AK8_PF_mass_JetCleansing_DiffMode;
	for(Int_t i=0;i<50;i++){
		TH2D h2_AK8_PF_mass_JetCleansing_DiffMode(Form("h2_AK8_PF_mass_JetCleansing_DiffMode%i",i),Form("h2_AK8_PF_mass_JetCleansing_DiffMode%i; gen jet mass; reco jet mass",i),nbin_mass,jetmass_min,jetmass_max,nbin_mass,jetmass_min,jetmass_max);
		vect_h2_AK8_PF_mass_JetCleansing_DiffMode.push_back(h2_AK8_PF_mass_JetCleansing_DiffMode);
		TGraph gr_AK8_PF_mass_JetCleansing_DiffMode;
		vect_gr_AK8_PF_mass_JetCleansing_DiffMode.push_back(gr_AK8_PF_mass_JetCleansing_DiffMode);
	}


	// For GEN-RECO matching
	Double_t gen_jet_eta=0.; Double_t gen_jet_phi=0.;
	Double_t pf_jet_eta=0.; Double_t pf_jet_phi=0.;

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;
		if (!Select())continue;

		gen_jet_eta=GenGroomedJet_AK8_GEN_eta[0];
		gen_jet_phi=GenGroomedJet_AK8_GEN_phi[0];

		pf_jet_eta=GroomedJet_AK8_PF_eta[0];
		pf_jet_phi=GroomedJet_AK8_PF_phi[0];

		//Gen Jet matching with PF_uncerr
		if( TMath::Sqrt( (gen_jet_eta-pf_jet_eta)*(gen_jet_eta-pf_jet_eta) + (gen_jet_phi-pf_jet_phi)*(gen_jet_phi-pf_jet_phi) ) <0.3 )
		{

			tmp_AK8_GEN_pt = GenGroomedJet_AK8_GEN_pt[0];
			tmp_AK8_GEN_eta = GenGroomedJet_AK8_GEN_eta[0];
			tmp_AK8_GEN_phi = GenGroomedJet_AK8_GEN_phi[0];

			tmp_AK8_PF_pt = GroomedJet_AK8_PF_pt[0];
			tmp_AK8_PF_eta = GroomedJet_AK8_PF_eta[0];
			tmp_AK8_PF_phi = GroomedJet_AK8_PF_phi[0];
			tmp_AK8_PF_pt_uncorr = GroomedJet_AK8_PF_pt_uncorr[0];
			tmp_AK8_PF_pt_L1_rhoSW = GroomedJet_AK8_PF_pt_L1_rhoSW[0];
			tmp_AK8_PF_pt_L1_rhoHand = GroomedJet_AK8_PF_pt_L1_rhoHand[0];
			tmp_AK8_PF_pt_L1_rhoHand2 = GroomedJet_AK8_PF_pt_L1_rhoHand2[0];
			tmp_AK8_PF_pt_L1_rhoGrid = GroomedJet_AK8_PF_pt_L1_rhoGrid[0];

			tmp_AK8_PF_pt_rho4A = GroomedJet_AK8_PF_pt_rho4A[0];
			tmp_AK8_PF_pt_rhom4A = GroomedJet_AK8_PF_pt_rhom4A[0];
			tmp_AK8_PF_pt_JetCleansing = GroomedJet_AK8_PF_pt_JetCleansing[0];

			tmp_event_nPV = event_nPV;

			tmp_AK8_GEN_rhoSW = GenGroomedJet_AK8_GEN_rhoSW;
			tmp_AK8_GEN_rhoHand = GenGroomedJet_AK8_GEN_rhohand;
			tmp_AK8_GEN_rhoHand2 = GenGroomedJet_AK8_GEN_rhohand2;
			tmp_AK8_GEN_rhoGrid = GenGroomedJet_AK8_GEN_rhogrid;

			tmp_AK8_PF_rhoSW = GroomedJet_AK8_PF_rhoSW;
			tmp_AK8_PF_rhoHand = GroomedJet_AK8_PF_rhohand;
			tmp_AK8_PF_rhoHand2 = GroomedJet_AK8_PF_rhohand2;
			tmp_AK8_PF_rhoGrid = GroomedJet_AK8_PF_rhogrid;

			//============= fill hist ==============
			//AK8 PF
			ratio = tmp_AK8_PF_pt/tmp_AK8_GEN_pt;
			h1_AK8_PF_recogenptratio.Fill(ratio);
			//AK8 uncorr
			ratio = tmp_AK8_PF_pt_uncorr/tmp_AK8_GEN_pt;
			h1_AK8_PF_uncorr_recogenptratio.Fill(ratio);
			//AK8 rhoSW
			ratio = tmp_AK8_PF_pt_L1_rhoSW/tmp_AK8_GEN_pt;
			h1_AK8_PF_l1rhosw_recogenptratio.Fill(ratio);
			//AK8 rhohand
			ratio = tmp_AK8_PF_pt_L1_rhoHand/tmp_AK8_GEN_pt;
			h1_AK8_PF_l1rhoHand_recogenptratio.Fill(ratio);
			//AK8 rhohand2
			ratio = tmp_AK8_PF_pt_L1_rhoHand2/tmp_AK8_GEN_pt;
			h1_AK8_PF_l1rhoHand2_recogenptratio.Fill(ratio);
			//AK8 rhogrid
			ratio = tmp_AK8_PF_pt_L1_rhoGrid/tmp_AK8_GEN_pt;
			h1_AK8_PF_l1rhoGrid_recogenptratio.Fill(ratio);

			ratio = tmp_AK8_PF_pt_rho4A/tmp_AK8_GEN_pt;
			h1_AK8_PF_rho4A_recogenptratio.Fill(ratio);

			ratio = tmp_AK8_PF_pt_rhom4A/tmp_AK8_GEN_pt;
			h1_AK8_PF_rhom4A_recogenptratio.Fill(ratio);

			ratio = tmp_AK8_PF_pt_JetCleansing/tmp_AK8_GEN_pt;
			h1_AK8_PF_JetCleansing_recogenptratio.Fill(ratio);

			//pt
			h1_AK8_GEN_pt.Fill(tmp_AK8_GEN_pt);
			h1_AK8_PF_pt_uncorr.Fill(tmp_AK8_PF_pt_uncorr);// pt of PF is with wrong JEC now
			h1_AK8_PF_pt_l1rhoHand.Fill(tmp_AK8_PF_pt_L1_rhoHand);// pt of PF is with wrong JEC now
			h1_AK8_PF_pt_l1rhoGrid.Fill(tmp_AK8_PF_pt_L1_rhoGrid);// pt of PF is with wrong JEC now
			h1_AK8_PF_pt_rho4A.Fill(tmp_AK8_PF_pt_rho4A);// 
			h1_AK8_PF_pt_rhom4A.Fill(tmp_AK8_PF_pt_rhom4A);// 
			h1_AK8_PF_pt_JetCleansing.Fill(tmp_AK8_PF_pt_JetCleansing);// 
			//eta
			h1_AK8_GEN_eta.Fill(tmp_AK8_GEN_eta);
			h1_AK8_PF_eta.Fill(tmp_AK8_PF_eta);
			//phi
			h1_AK8_GEN_phi.Fill(tmp_AK8_GEN_phi);
			h1_AK8_PF_phi.Fill(tmp_AK8_PF_phi);
			//dr
			dr = TMath::Sqrt( (Z_eta-tmp_AK8_GEN_eta)*(Z_eta-tmp_AK8_GEN_eta) + (Z_phi-tmp_AK8_GEN_phi)*(Z_phi-tmp_AK8_GEN_phi) );
			h1_AK8_GEN_zjet_dr.Fill(dr);
			dr = TMath::Sqrt( (Z_eta-tmp_AK8_PF_eta)*(Z_eta-tmp_AK8_PF_eta) + (Z_phi-tmp_AK8_PF_phi)*(Z_phi-tmp_AK8_PF_phi) );
			h1_AK8_PF_zjet_dr.Fill(dr);
			//dphi
			dphi = TMath::Sqrt( (Z_phi-tmp_AK8_GEN_phi)*(Z_phi-tmp_AK8_GEN_phi) );
			h1_AK8_GEN_zjet_dphi.Fill(dphi);
			dphi = TMath::Sqrt( (Z_phi-tmp_AK8_PF_phi)*(Z_phi-tmp_AK8_PF_phi) );
			h1_AK8_PF_zjet_dphi.Fill(dphi);
			//nPV
			h1_nPV.Fill(tmp_event_nPV);
			//z mass
			h1_z_mass.Fill(Z_mass);
			h1_z_Pt.Fill(Z_pt);
			//muplus pt
			h1_muplus_Pt.Fill(Z_muplus_pt);

			//rho
			h1_AK8_GEN_rhoSW.Fill(tmp_AK8_GEN_rhoSW);
			h1_AK8_PF_rhoSW.Fill(tmp_AK8_PF_rhoSW);

			h1_AK8_GEN_rhoHand.Fill(tmp_AK8_GEN_rhoHand);
			h1_AK8_PF_rhoHand.Fill(tmp_AK8_PF_rhoHand);

			h1_AK8_GEN_rhoHand2.Fill(tmp_AK8_GEN_rhoHand2);
			h1_AK8_PF_rhoHand2.Fill(tmp_AK8_PF_rhoHand2);

			h1_AK8_GEN_rhoGrid.Fill(tmp_AK8_GEN_rhoGrid);
			h1_AK8_PF_rhoGrid.Fill(tmp_AK8_PF_rhoGrid);

			h1_AK8_PF_area.Fill(GroomedJet_AK8_PF_area[0]);  

			//rho vs nPV
			h2_AK8_GEN_rhoSW_vs_nPV.Fill(tmp_AK8_GEN_rhoSW,tmp_event_nPV);
			h2_AK8_GEN_rhoHand_vs_nPV.Fill(tmp_AK8_GEN_rhoHand,tmp_event_nPV);
			h2_AK8_GEN_rhoHand2_vs_nPV.Fill(tmp_AK8_GEN_rhoHand2,tmp_event_nPV);
			h2_AK8_GEN_rhoGrid_vs_nPV.Fill(tmp_AK8_GEN_rhoGrid,tmp_event_nPV);

			h2_AK8_PF_rhoSW_vs_nPV.Fill(tmp_AK8_PF_rhoSW,tmp_event_nPV);
			h2_AK8_PF_rhoHand_vs_nPV.Fill(tmp_AK8_PF_rhoHand,tmp_event_nPV);
			h2_AK8_PF_rhoHand2_vs_nPV.Fill(tmp_AK8_PF_rhoHand2,tmp_event_nPV);
			h2_AK8_PF_rhoGrid_vs_nPV.Fill(tmp_AK8_PF_rhoGrid,tmp_event_nPV);


			// jet mass
			tmp_AK8_GEN_mass=GenGroomedJet_AK8_GEN_mass_uncorr[0];
			tmp_AK8_PF_mass_uncorr=GroomedJet_AK8_PF_mass_uncorr[0];
			tmp_AK8_PF_mass_rhoArea=GroomedJet_AK8_PF_mass_rhoArea[0];
			tmp_AK8_PF_mass_rhoGArea=GroomedJet_AK8_PF_mass_rhoGArea[0];
			tmp_AK8_PF_mass_rho4Area=GroomedJet_AK8_PF_mass_rho4Area[0];
			tmp_AK8_PF_mass_rhoG4Area=GroomedJet_AK8_PF_mass_rhoG4Area[0];
			tmp_AK8_PF_mass_rhom4Area=GroomedJet_AK8_PF_mass_rhom4Area[0];
			tmp_AK8_PF_mass_JetCleansingATLASjvf=GroomedJet_AK8_PF_mass_JetCleansingATLASjvf[0];
			tmp_AK8_PF_mass_JetCleansingATLASlin=GroomedJet_AK8_PF_mass_JetCleansingATLASlin[0];
			tmp_AK8_PF_mass_JetCleansingATLASgau=GroomedJet_AK8_PF_mass_JetCleansingATLASgau[0];
			tmp_AK8_PF_mass_JetCleansingCMSjvf=GroomedJet_AK8_PF_mass_JetCleansingCMSjvf[0];
			tmp_AK8_PF_mass_JetCleansingCMSlin=GroomedJet_AK8_PF_mass_JetCleansingCMSlin[0];
			tmp_AK8_PF_mass_JetCleansingCMSgau=GroomedJet_AK8_PF_mass_JetCleansingCMSgau[0];

			h1_AK8_GEN_mass.Fill(tmp_AK8_GEN_mass            ); 
			h1_AK8_PF_mass_uncorr.Fill(tmp_AK8_PF_mass_uncorr  );
			h1_AK8_PF_mass_rhoArea.Fill(tmp_AK8_PF_mass_rhoArea );
			h1_AK8_PF_mass_rhoGArea.Fill(tmp_AK8_PF_mass_rhoGArea);
			h1_AK8_PF_mass_rho4Area.Fill(tmp_AK8_PF_mass_rho4Area);
			h1_AK8_PF_mass_rhoG4Area.Fill(tmp_AK8_PF_mass_rhoG4Area);
			h1_AK8_PF_mass_rhom4Area.Fill(tmp_AK8_PF_mass_rhom4Area);
			h1_AK8_PF_mass_JetCleansingATLASjvf.Fill(tmp_AK8_PF_mass_JetCleansingATLASjvf);
			h1_AK8_PF_mass_JetCleansingATLASlin.Fill(tmp_AK8_PF_mass_JetCleansingATLASlin);
			h1_AK8_PF_mass_JetCleansingATLASgau.Fill(tmp_AK8_PF_mass_JetCleansingATLASgau);
			h1_AK8_PF_mass_JetCleansingCMSjvf.Fill(tmp_AK8_PF_mass_JetCleansingCMSjvf);
			h1_AK8_PF_mass_JetCleansingCMSlin.Fill(tmp_AK8_PF_mass_JetCleansingCMSlin);
			h1_AK8_PF_mass_JetCleansingCMSgau.Fill(tmp_AK8_PF_mass_JetCleansingCMSgau);

			h2_AK8_GEN_mass.Fill(tmp_AK8_GEN_mass, tmp_AK8_GEN_mass            ); 
			h2_AK8_PF_mass_uncorr.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_uncorr  );
			h2_AK8_PF_mass_rhoArea.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_rhoArea );
			h2_AK8_PF_mass_rhoGArea.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_rhoGArea);
			h2_AK8_PF_mass_rho4Area.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_rho4Area);
			h2_AK8_PF_mass_rhoG4Area.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_rhoG4Area);
			h2_AK8_PF_mass_rhom4Area.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_rhom4Area);
			h2_AK8_PF_mass_JetCleansingATLASjvf.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_JetCleansingATLASjvf);
			h2_AK8_PF_mass_JetCleansingATLASlin.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_JetCleansingATLASlin);
			h2_AK8_PF_mass_JetCleansingATLASgau.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_JetCleansingATLASgau);
			h2_AK8_PF_mass_JetCleansingCMSjvf.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_JetCleansingCMSjvf);
			h2_AK8_PF_mass_JetCleansingCMSlin.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_JetCleansingCMSlin);
			h2_AK8_PF_mass_JetCleansingCMSgau.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_JetCleansingCMSgau);


			gr_AK8_GEN_mass.SetPoint(num_points, tmp_AK8_GEN_mass, tmp_AK8_GEN_mass            ); 
			gr_AK8_PF_mass_uncorr.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_uncorr  );
			gr_AK8_PF_mass_rhoArea.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_rhoArea );
			gr_AK8_PF_mass_rhoGArea.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_rhoGArea);
			gr_AK8_PF_mass_rho4Area.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_rho4Area);
			gr_AK8_PF_mass_rhoG4Area.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_rhoG4Area);
			gr_AK8_PF_mass_rhom4Area.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_rhom4Area);
			gr_AK8_PF_mass_JetCleansingATLASjvf.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_JetCleansingATLASjvf);
			gr_AK8_PF_mass_JetCleansingATLASlin.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_JetCleansingATLASlin);
			gr_AK8_PF_mass_JetCleansingATLASgau.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_JetCleansingATLASgau);
			gr_AK8_PF_mass_JetCleansingCMSjvf.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_JetCleansingCMSjvf);
			gr_AK8_PF_mass_JetCleansingCMSlin.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_JetCleansingCMSlin);
			gr_AK8_PF_mass_JetCleansingCMSgau.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_JetCleansingCMSgau);

			Int_t tmp_number_JetCleansing_DiffMode=number_JetCleansing_DiffMode;
			if(tmp_number_JetCleansing_DiffMode==0)tmp_number_JetCleansing_DiffMode=50;
			for(Int_t k=0;k<tmp_number_JetCleansing_DiffMode;k++){
				if (GroomedJet_AK8_PF_mass_JetCleansing_DiffMode[k]>=0 && GroomedJet_AK8_PF_pt_JetCleansing_DiffMode[k]>=0){
					vect_h1_AK8_PF_mass_JetCleansing_DiffMode[k].Fill(GroomedJet_AK8_PF_mass_JetCleansing_DiffMode[k]);
					vect_h1_AK8_PF_pt_JetCleansing_DiffMode[k].Fill(GroomedJet_AK8_PF_pt_JetCleansing_DiffMode[k]);
					vect_h2_AK8_PF_mass_JetCleansing_DiffMode[k].Fill(tmp_AK8_GEN_mass, GroomedJet_AK8_PF_mass_JetCleansing_DiffMode[k]);
					vect_gr_AK8_PF_mass_JetCleansing_DiffMode[k].SetPoint(num_points,tmp_AK8_GEN_mass, GroomedJet_AK8_PF_mass_JetCleansing_DiffMode[k]);

					if(tmp_number_JetCleansing_DiffMode==50)number_JetCleansing_DiffMode++;
				}else{ break;}
			}
			num_points++;

		}
		//PF matching efficiency
		dr = TMath::Sqrt( (gen_jet_eta-tmp_AK8_PF_eta)*(gen_jet_eta-tmp_AK8_PF_eta) + (gen_jet_phi-tmp_AK8_PF_phi)*(gen_jet_phi-tmp_AK8_PF_phi) ) ;
		h1_PF_match.Fill(dr);
	}


	Draw_and_Save(h1_AK8_PF_recogenptratio);
	Draw_and_Save(h1_AK8_PF_uncorr_recogenptratio);
	Draw_and_Save(h1_AK8_PF_l1rhosw_recogenptratio);
	Draw_and_Save(h1_AK8_PF_l1rhoHand_recogenptratio);
	Draw_and_Save(h1_AK8_PF_l1rhoHand2_recogenptratio);
	Draw_and_Save(h1_AK8_PF_l1rhoGrid_recogenptratio);
	Draw_and_Save(h1_AK8_PF_rho4A_recogenptratio);
	Draw_and_Save(h1_AK8_PF_rhom4A_recogenptratio);
	Draw_and_Save(h1_AK8_PF_JetCleansing_recogenptratio);

	h1_AK8_PF_uncorr_recogenptratio.SetLineColor(1); h1_AK8_PF_uncorr_recogenptratio.SetLineStyle(2); h1_AK8_PF_uncorr_recogenptratio.SetLineWidth(2);
	h1_AK8_PF_l1rhoHand_recogenptratio.SetLineColor(2);h1_AK8_PF_l1rhoHand_recogenptratio.SetLineStyle(2);h1_AK8_PF_l1rhoHand_recogenptratio.SetLineWidth(2);
	h1_AK8_PF_rho4A_recogenptratio.SetLineColor(3);h1_AK8_PF_rho4A_recogenptratio.SetLineStyle(1);h1_AK8_PF_rho4A_recogenptratio.SetLineWidth(1);
	h1_AK8_PF_rhom4A_recogenptratio.SetLineColor(4);h1_AK8_PF_rhom4A_recogenptratio.SetLineStyle(2);h1_AK8_PF_rhom4A_recogenptratio.SetLineWidth(2);
	h1_AK8_PF_JetCleansing_recogenptratio.SetLineColor(6);h1_AK8_PF_JetCleansing_recogenptratio.SetLineStyle(1);h1_AK8_PF_JetCleansing_recogenptratio.SetLineWidth(1);
	Draw_and_Save(h1_AK8_PF_uncorr_recogenptratio, h1_AK8_PF_l1rhoHand_recogenptratio, h1_AK8_PF_rho4A_recogenptratio );
	Draw_and_Save(h1_AK8_PF_uncorr_recogenptratio, h1_AK8_PF_rhom4A_recogenptratio, h1_AK8_PF_JetCleansing_recogenptratio );

	h1_AK8_GEN_pt.SetLineColor(1); 
	h1_AK8_PF_pt_uncorr.SetLineColor(1); h1_AK8_PF_pt_uncorr.SetLineStyle(2); h1_AK8_PF_pt_uncorr.SetLineWidth(2);
	h1_AK8_PF_pt_l1rhoHand.SetLineColor(2);h1_AK8_PF_pt_l1rhoHand.SetLineStyle(2);h1_AK8_PF_pt_l1rhoHand.SetLineWidth(2);
	h1_AK8_PF_pt_rho4A.SetLineColor(3);h1_AK8_PF_pt_rho4A.SetLineStyle(1);
	h1_AK8_PF_pt_rhom4A.SetLineColor(4);h1_AK8_PF_pt_rhom4A.SetLineStyle(2);h1_AK8_PF_pt_rhom4A.SetLineWidth(2);
	h1_AK8_PF_pt_JetCleansing.SetLineColor(6);h1_AK8_PF_pt_JetCleansing.SetLineStyle(1);
	Draw_and_Save(h1_AK8_GEN_pt, h1_AK8_PF_pt_uncorr, h1_AK8_PF_pt_l1rhoHand, h1_AK8_PF_pt_rho4A );
	Draw_and_Save(h1_AK8_GEN_pt, h1_AK8_PF_pt_uncorr, h1_AK8_PF_pt_rhom4A, h1_AK8_PF_pt_JetCleansing );

	Draw_and_Save(h1_AK8_PF_area);

	Draw_and_Save(h1_nPV);
	Draw_and_Save(h1_z_mass);
	Draw_and_Save(h1_z_Pt);
	Draw_and_Save(h1_muplus_Pt);


	Draw_and_Save(h1_AK8_GEN_pt);
	Draw_and_Save(h1_AK8_GEN_eta);
	Draw_and_Save(h1_AK8_GEN_zjet_dr);
	Draw_and_Save(h1_AK8_GEN_zjet_dphi);
	/*Draw_and_Save(h1_AK8_GEN_rhoSW);
	  Draw_and_Save(h1_AK8_GEN_rhoHand);
	  Draw_and_Save(h1_AK8_GEN_rhoHand2);
	  Draw_and_Save(h1_AK8_GEN_rhoGrid);
	  Draw_and_Save(h2_AK8_GEN_rhoSW_vs_nPV);
	  Draw_and_Save(h2_AK8_GEN_rhoHand_vs_nPV);
	  Draw_and_Save(h2_AK8_GEN_rhoHand2_vs_nPV);
	  Draw_and_Save(h2_AK8_GEN_rhoGrid_vs_nPV);*/

	Draw_and_Save(h1_AK8_PF_pt_uncorr);
	Draw_and_Save(h1_AK8_PF_pt_l1rhoHand);
	Draw_and_Save(h1_AK8_PF_pt_l1rhoGrid);
	Draw_and_Save(h1_AK8_PF_pt_rho4A);
	Draw_and_Save(h1_AK8_PF_pt_rhom4A);
	Draw_and_Save(h1_AK8_PF_pt_JetCleansing);
	Draw_and_Save(h1_AK8_PF_eta);
	Draw_and_Save(h1_AK8_PF_zjet_dr);
	Draw_and_Save(h1_AK8_PF_zjet_dphi);
	Draw_and_Save(h1_AK8_PF_rhoSW);
	Draw_and_Save(h1_AK8_PF_rhoHand);
	Draw_and_Save(h1_AK8_PF_rhoHand2);
	Draw_and_Save(h1_AK8_PF_rhoGrid);
	Draw_and_Save(h2_AK8_PF_rhoSW_vs_nPV);
	Draw_and_Save(h2_AK8_PF_rhoHand_vs_nPV);
	Draw_and_Save(h2_AK8_PF_rhoHand2_vs_nPV);
	Draw_and_Save(h2_AK8_PF_rhoGrid_vs_nPV);

	Draw_and_Save(h1_PF_match);

	Draw_and_Save(h1_AK8_GEN_mass                );
	Draw_and_Save(h1_AK8_PF_mass_uncorr      );
	Draw_and_Save(h1_AK8_PF_mass_rhoArea     );
	Draw_and_Save(h1_AK8_PF_mass_rhoGArea    );
	Draw_and_Save(h1_AK8_PF_mass_rho4Area    );
	Draw_and_Save(h1_AK8_PF_mass_rhoG4Area   );
	Draw_and_Save(h1_AK8_PF_mass_rhom4Area   );
	Draw_and_Save(h1_AK8_PF_mass_JetCleansingATLASjvf   );
	Draw_and_Save(h1_AK8_PF_mass_JetCleansingATLASlin   );
	Draw_and_Save(h1_AK8_PF_mass_JetCleansingATLASgau   );
	Draw_and_Save(h1_AK8_PF_mass_JetCleansingCMSjvf   );
	Draw_and_Save(h1_AK8_PF_mass_JetCleansingCMSlin   );
	Draw_and_Save(h1_AK8_PF_mass_JetCleansingCMSgau   );

	cout<<"number_JetCleansing_DiffMode="<<number_JetCleansing_DiffMode<<endl;
	if (number_JetCleansing_DiffMode>50) number_JetCleansing_DiffMode=50;
	for(Int_t k=0;k<number_JetCleansing_DiffMode;k++){
		Draw_and_Save(vect_h1_AK8_PF_mass_JetCleansing_DiffMode[k]);
		Draw_and_Save(vect_h1_AK8_PF_pt_JetCleansing_DiffMode[k]);
		Draw_and_Save(vect_h2_AK8_PF_mass_JetCleansing_DiffMode[k]	, Form("%g correlated", vect_gr_AK8_PF_mass_JetCleansing_DiffMode[k].GetCorrelationFactor()));
	}

	Draw_and_Save(h2_AK8_GEN_mass, Form("%g correlated", gr_AK8_GEN_mass.GetCorrelationFactor())   );
	Draw_and_Save(h2_AK8_PF_mass_uncorr				, Form("%g correlated", gr_AK8_PF_mass_uncorr.GetCorrelationFactor())				 );
	Draw_and_Save(h2_AK8_PF_mass_rhoArea			, Form("%g correlated", gr_AK8_PF_mass_rhoArea.GetCorrelationFactor())				);
	Draw_and_Save(h2_AK8_PF_mass_rhoGArea			, Form("%g correlated", gr_AK8_PF_mass_rhoGArea.GetCorrelationFactor())				);
	Draw_and_Save(h2_AK8_PF_mass_rho4Area			, Form("%g correlated", gr_AK8_PF_mass_rho4Area.GetCorrelationFactor())				);
	Draw_and_Save(h2_AK8_PF_mass_rhoG4Area			, Form("%g correlated", gr_AK8_PF_mass_rhoG4Area.GetCorrelationFactor())				);
	Draw_and_Save(h2_AK8_PF_mass_rhom4Area			, Form("%g correlated", gr_AK8_PF_mass_rhom4Area.GetCorrelationFactor())				);
	Draw_and_Save(h2_AK8_PF_mass_JetCleansingATLASjvf	, Form("%g correlated", gr_AK8_PF_mass_JetCleansingATLASjvf.GetCorrelationFactor())		);
	Draw_and_Save(h2_AK8_PF_mass_JetCleansingATLASlin	, Form("%g correlated", gr_AK8_PF_mass_JetCleansingATLASlin.GetCorrelationFactor())		);
	Draw_and_Save(h2_AK8_PF_mass_JetCleansingATLASgau	, Form("%g correlated", gr_AK8_PF_mass_JetCleansingATLASgau.GetCorrelationFactor())		);
	Draw_and_Save(h2_AK8_PF_mass_JetCleansingCMSjvf	, Form("%g correlated", gr_AK8_PF_mass_JetCleansingCMSjvf.GetCorrelationFactor())		);
	Draw_and_Save(h2_AK8_PF_mass_JetCleansingCMSlin	, Form("%g correlated", gr_AK8_PF_mass_JetCleansingCMSlin.GetCorrelationFactor())		);
	Draw_and_Save(h2_AK8_PF_mass_JetCleansingCMSgau	, Form("%g correlated", gr_AK8_PF_mass_JetCleansingCMSgau.GetCorrelationFactor())		);
}
