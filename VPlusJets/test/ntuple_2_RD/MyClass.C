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
   if(!( GenGroomedJet_AK5_GEN_pt[0]>100 ))return 0;

   Double_t tmp_AK5_GEN_eta = GenGroomedJet_AK5_GEN_eta[0];
   Double_t tmp_AK5_GEN_phi = GenGroomedJet_AK5_GEN_phi[0];
   Double_t tmpDeltaPhi_Vj = TMath::Sqrt( (Z_phi-tmp_AK5_GEN_phi)*(Z_phi-tmp_AK5_GEN_phi) );
   Double_t tmpDeltaR_Vj = TMath::Sqrt( (Z_eta-tmp_AK5_GEN_eta)*(Z_eta-tmp_AK5_GEN_eta) + (Z_phi-tmp_AK5_GEN_phi)*(Z_phi-tmp_AK5_GEN_phi) );

   if(!( tmpDeltaPhi_Vj>2.0 && tmpDeltaR_Vj>1.0 ))return 0;
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
	Double_t tmp_AK5_PF_mass_cleansingATLASjvf=0.;
	Double_t tmp_AK5_PF_mass_cleansingATLASlin=0.;
	Double_t tmp_AK5_PF_mass_cleansingATLASgau=0.;
	Double_t tmp_AK5_PF_mass_cleansingCMSjvf=0.;
	Double_t tmp_AK5_PF_mass_cleansingCMSlin=0.;
	Double_t tmp_AK5_PF_mass_cleansingCMSgau=0.;


	Double_t tmp_AK5_PFCHS_mass_uncorr=0.;
	Double_t tmp_AK5_PFCHS_mass_rhoArea=0.;
	Double_t tmp_AK5_PFCHS_mass_rhoGArea=0.;
	Double_t tmp_AK5_PFCHS_mass_rho4Area=0.;
	Double_t tmp_AK5_PFCHS_mass_rhoG4Area=0.;
	Double_t tmp_AK5_PFCHS_mass_rhom4Area=0.;

	Double_t rhomin=0.; Double_t rhomax=50.;

	TH1D h1_nPV("h1_nPV","h1_nPV",50,0,50);

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

	// JetPt/GenPt
	TH1D h1_PFCor("h1_PFCor","h1_PFCor",50,0,2.5);
	TH1D h1_PFCor_uncorr("h1_PFCor_uncorr","h1_PFCor_uncorr",50,0,2.5);
	TH1D h1_PFCor_afterL1("h1_PFCor_afterL1","h1_PFCor_afterL1",50,0,2.5);
	TH1D h1_PFCor_afterL2("h1_PFCor_afterL2","h1_PFCor_afterL2",50,0,2.5);

	TH1D h1_AK5_PF("h1_AK5_PF","h1_AK5_PF",50,0,2.5);
	TH1D h1_AK5_PF_uncorr("h1_AK5_PF_uncorr","h1_AK5_PF_uncorr",50,0,2.5);
	TH1D h1_AK5_PF_l1_rhoSW("h1_AK5_PF_l1_rhoSW","h1_AK5_PF_l1_rhoSW",50,0,2.5);
	TH1D h1_AK5_PF_l1_rhoHand("h1_AK5_PF_l1_rhoHand","h1_AK5_PF_l1_rhoHand",50,0,2.5);
	TH1D h1_AK5_PF_l1_rhoHand2("h1_AK5_PF_l1_rhoHand2","h1_AK5_PF_l1_rhoHand2",50,0,2.5);
	TH1D h1_AK5_PF_l1_rhoGrid("h1_AK5_PF_l1_rhoGrid","h1_AK5_PF_l1_rhoGrid",50,0,2.5);

	TH1D h1_AK5_PFCHS("h1_AK5_PFCHS","h1_AK5_PFCHS",50,0,2.5);
	TH1D h1_AK5_PFCHS_uncorr("h1_AK5_PFCHS_uncorr","h1_AK5_PFCHS_uncorr",50,0,2.5);
	TH1D h1_AK5_PFCHS_l1_rhoSW("h1_AK5_PFCHS_l1_rhoSW","h1_AK5_PFCHS_l1_rhoSW",50,0,2.5);
	TH1D h1_AK5_PFCHS_l1_rhoHand("h1_AK5_PFCHS_l1_rhoHand","h1_AK5_PFCHS_l1_rhoHand",50,0,2.5);
	TH1D h1_AK5_PFCHS_l1_rhoHand2("h1_AK5_PFCHS_l1_rhoHand2","h1_AK5_PFCHS_l1_rhoHand2",50,0,2.5);
	TH1D h1_AK5_PFCHS_l1_rhoGrid("h1_AK5_PFCHS_l1_rhoGrid","h1_AK5_PFCHS_l1_rhoGrid",50,0,2.5);

	//area
	TH1D h1_area_PFCor("h1_area_PFCor","h1_area_PFCor",50,0.6,1.1);
	TH1D h1_area_PF("h1_area_PF","h1_area_PF",50,0.6,1.1);
	TH1D h1_area_PFCHS("h1_area_PFCHS","h1_area_PFCHS",50,0.6,1.1);

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
	TH1D h1_AK5_PF_mass_cleansingATLASjvf("h1_AK5_PF_mass_cleansingATLASjvf","h1_AK5_PF_mass_cleansingATLASjvf;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK5_PF_mass_cleansingATLASlin("h1_AK5_PF_mass_cleansingATLASlin","h1_AK5_PF_mass_cleansingATLASlin;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK5_PF_mass_cleansingATLASgau("h1_AK5_PF_mass_cleansingATLASgau","h1_AK5_PF_mass_cleansingATLASgau;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK5_PF_mass_cleansingCMSjvf("h1_AK5_PF_mass_cleansingCMSjvf","h1_AK5_PF_mass_cleansingCMSjvf;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK5_PF_mass_cleansingCMSlin("h1_AK5_PF_mass_cleansingCMSlin","h1_AK5_PF_mass_cleansingCMSlin;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK5_PF_mass_cleansingCMSgau("h1_AK5_PF_mass_cleansingCMSgau","h1_AK5_PF_mass_cleansingCMSgau;jet mass;",nbin_mass,jetmass_min,jetmass_max);

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

			h1_area_PFCor.Fill(JetPFCor_Area[i_PFCorJet_matching_PFCHS]);
			h1_area_PF.Fill(GroomedJet_AK5_PF_area[0]);  
			h1_area_PFCHS.Fill(GroomedJet_AK5_PFCHS_area[0]); 

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
			tmp_AK5_PF_mass_cleansingATLASjvf=GroomedJet_AK5_PF_mass_cleansingATLASjvf[0];
			tmp_AK5_PF_mass_cleansingATLASlin=GroomedJet_AK5_PF_mass_cleansingATLASlin[0];
			tmp_AK5_PF_mass_cleansingATLASgau=GroomedJet_AK5_PF_mass_cleansingATLASgau[0];
			tmp_AK5_PF_mass_cleansingCMSjvf=GroomedJet_AK5_PF_mass_cleansingCMSjvf[0];
			tmp_AK5_PF_mass_cleansingCMSlin=GroomedJet_AK5_PF_mass_cleansingCMSlin[0];
			tmp_AK5_PF_mass_cleansingCMSgau=GroomedJet_AK5_PF_mass_cleansingCMSgau[0];
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
			h1_AK5_PF_mass_cleansingATLASjvf.Fill(tmp_AK5_PF_mass_cleansingATLASjvf);
			h1_AK5_PF_mass_cleansingATLASlin.Fill(tmp_AK5_PF_mass_cleansingATLASlin);
			h1_AK5_PF_mass_cleansingATLASgau.Fill(tmp_AK5_PF_mass_cleansingATLASgau);
			h1_AK5_PF_mass_cleansingCMSjvf.Fill(tmp_AK5_PF_mass_cleansingCMSjvf);
			h1_AK5_PF_mass_cleansingCMSlin.Fill(tmp_AK5_PF_mass_cleansingCMSlin);
			h1_AK5_PF_mass_cleansingCMSgau.Fill(tmp_AK5_PF_mass_cleansingCMSgau);

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
	Draw_and_Save(h1_PFCor_afterL2);
	Draw_and_Save(h1_AK5_PF);
	Draw_and_Save(h1_AK5_PF_uncorr);
	Draw_and_Save(h1_AK5_PF_l1_rhoSW);
	Draw_and_Save(h1_AK5_PF_l1_rhoHand);
	Draw_and_Save(h1_AK5_PF_l1_rhoHand2);
	Draw_and_Save(h1_AK5_PF_l1_rhoGrid);
	Draw_and_Save(h1_AK5_PFCHS);
	Draw_and_Save(h1_AK5_PFCHS_uncorr);
	Draw_and_Save(h1_AK5_PFCHS_l1_rhoSW);
	Draw_and_Save(h1_AK5_PFCHS_l1_rhoHand);
	Draw_and_Save(h1_AK5_PFCHS_l1_rhoHand2);
	Draw_and_Save(h1_AK5_PFCHS_l1_rhoGrid);

	Draw_and_Save(h1_area_PFCor);
	Draw_and_Save(h1_area_PF);
	Draw_and_Save(h1_area_PFCHS);

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
	Draw_and_Save(h1_AK5_PF_mass_cleansingATLASjvf   );
	Draw_and_Save(h1_AK5_PF_mass_cleansingATLASlin   );
	Draw_and_Save(h1_AK5_PF_mass_cleansingATLASgau   );
	Draw_and_Save(h1_AK5_PF_mass_cleansingCMSjvf   );
	Draw_and_Save(h1_AK5_PF_mass_cleansingCMSlin   );
	Draw_and_Save(h1_AK5_PF_mass_cleansingCMSgau   );
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
	//h1_AK5_PF_mass_cleansingATLASjvf.Draw("same");
	//h1_AK5_PF_mass_cleansingATLASlin.Draw("same");
	//h1_AK5_PF_mass_cleansingATLASgau.Draw("same");
	//h1_AK5_PF_mass_cleansingCMSjvf.Draw("same");
	h1_AK5_PF_mass_cleansingCMSlin.SetLineColor(6);
	h1_AK5_PF_mass_cleansingCMSlin.Draw("same");
	//h1_AK5_PF_mass_cleansingCMSgau.Draw("same");
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
	Double_t tmp_AK8_PF_mass_cleansingATLASjvf=0.;
	Double_t tmp_AK8_PF_mass_cleansingATLASlin=0.;
	Double_t tmp_AK8_PF_mass_cleansingATLASgau=0.;
	Double_t tmp_AK8_PF_mass_cleansingCMSjvf=0.;
	Double_t tmp_AK8_PF_mass_cleansingCMSlin=0.;
	Double_t tmp_AK8_PF_mass_cleansingCMSgau=0.;


	Int_t nbin_rho=50; Double_t rhomin=0.; Double_t rhomax=50.;
	Int_t nbin_mass=60;Double_t jetmass_min=0;Double_t jetmass_max=300.;

	TH1D h1_nPV("h1_nPV","h1_nPV",50,0,50);

	TH1D h1_AK8_GEN_pt("h1_AK8_GEN_pt","h1_AK8_GEN_pt",50,0,200); h1_AK8_GEN_pt.SetLineColor(kRed);
	TH1D h1_AK8_GEN_eta("h1_AK8_GEN_eta","h1_AK8_GEN_eta",50,-2.5,2.5); h1_AK8_GEN_eta.SetLineColor(kRed);
	TH1D h1_AK8_GEN_phi("h1_AK8_GEN_phi","h1_AK8_GEN_phi",50,-4,4); h1_AK8_GEN_phi.SetLineColor(kRed);
	TH1D h1_AK8_GEN_zjet_dr("h1_AK8_GEN_zjet_dr","h1_AK8_GEN_zjet_dr",50,0,10); h1_AK8_GEN_zjet_dr.SetLineColor(kRed);
	TH1D h1_AK8_GEN_zjet_dphi("h1_AK8_GEN_zjet_dphi","h1_AK8_GEN_zjet_dphi",50,0,5); h1_AK8_GEN_zjet_dphi.SetLineColor(kRed);

	TH1D h1_AK8_GEN_rhoSW("h1_AK8_GEN_rhoSW","h1_AK8_GEN_rhoSW",nbin_rho,rhomin,rhomax); h1_AK8_GEN_rhoSW.SetLineColor(kRed);
	TH1D h1_AK8_GEN_rhoHand("h1_AK8_GEN_rhoHand","h1_AK8_GEN_rhoHand",nbin_rho,rhomin,rhomax); h1_AK8_GEN_rhoHand.SetLineColor(kRed);
	TH1D h1_AK8_GEN_rhoHand2("h1_AK8_GEN_rhoHand2","h1_AK8_GEN_rhoHand2",nbin_rho,rhomin,rhomax); h1_AK8_GEN_rhoHand2.SetLineColor(kRed);
	TH1D h1_AK8_GEN_rhoGrid("h1_AK8_GEN_rhoGrid","h1_AK8_GEN_rhoGrid",nbin_rho,rhomin,rhomax); h1_AK8_GEN_rhoGrid.SetLineColor(kRed);
	TH2D h2_AK8_GEN_rhoSW_vs_nPV("h2_AK8_GEN_rhoSW_vs_nPV","h2_AK8_GEN_rhoSW_vs_nPV",nbin_rho,rhomin,rhomax,50,0,50); h2_AK8_GEN_rhoSW_vs_nPV.SetMarkerColor(kRed);
	TH2D h2_AK8_GEN_rhoHand_vs_nPV("h2_AK8_GEN_rhoHand_vs_nPV","h2_AK8_GEN_rhoHand_vs_nPV",nbin_rho,rhomin,rhomax,50,0,50); h2_AK8_GEN_rhoHand_vs_nPV.SetMarkerColor(kRed);
	TH2D h2_AK8_GEN_rhoHand2_vs_nPV("h2_AK8_GEN_rhoHand2_vs_nPV","h2_AK8_GEN_rhoHand2_vs_nPV",nbin_rho,rhomin,rhomax,50,0,50); h2_AK8_GEN_rhoHand2_vs_nPV.SetMarkerColor(kRed);
	TH2D h2_AK8_GEN_rhoGrid_vs_nPV("h2_AK8_GEN_rhoGrid_vs_nPV","h2_AK8_GEN_rhoGrid_vs_nPV",nbin_rho,rhomin,rhomax,50,0,50); h2_AK8_GEN_rhoGrid_vs_nPV.SetMarkerColor(kRed);

	TH1D h1_AK8_PF_pt_uncorr("h1_AK8_PF_pt_uncorr","h1_AK8_PF_pt_uncorr",50,0,200);
	TH1D h1_AK8_PF_eta("h1_AK8_PF_eta","h1_AK8_PF_eta",50,-2.5,2.5);
	TH1D h1_AK8_PF_phi("h1_AK8_PF_phi","h1_AK8_PF_phi",50,-4,4);
	TH1D h1_AK8_PF_zjet_dr("h1_AK8_PF_zjet_dr","h1_AK8_PF_zjet_dr",50,0,10);
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
	TH1D h1_AK8_PF("h1_AK8_PF","h1_AK8_PF",50,0,2.5);
	TH1D h1_AK8_PF_uncorr("h1_AK8_PF_uncorr","h1_AK8_PF_uncorr",50,0,2.5);
	TH1D h1_AK8_PF_l1_rhoSW("h1_AK8_PF_l1_rhoSW","h1_AK8_PF_l1_rhoSW",50,0,2.5);
	TH1D h1_AK8_PF_l1_rhoHand("h1_AK8_PF_l1_rhoHand","h1_AK8_PF_l1_rhoHand",50,0,2.5);
	TH1D h1_AK8_PF_l1_rhoHand2("h1_AK8_PF_l1_rhoHand2","h1_AK8_PF_l1_rhoHand2",50,0,2.5);
	TH1D h1_AK8_PF_l1_rhoGrid("h1_AK8_PF_l1_rhoGrid","h1_AK8_PF_l1_rhoGrid",50,0,2.5);

	//area
	TH1D h1_area_PF("h1_area_PF","h1_area_PF",50,0.6,1.1);

	//mass
	TH1D h1_AK8_GEN_mass("h1_AK8_GEN_mass","h1_AK8_GEN_mass;jet mass;",nbin_mass,jetmass_min,jetmass_max); h1_AK8_GEN_mass.SetLineColor(kRed);
	TH1D h1_AK8_PF_mass_uncorr("h1_AK8_PF_mass_uncorr","h1_AK8_PF_mass_uncorr;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_rhoArea("h1_AK8_PF_mass_rhoArea","h1_AK8_PF_mass_rhoArea;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_rhoGArea("h1_AK8_PF_mass_rhoGArea","h1_AK8_PF_mass_rhoGArea;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_rho4Area("h1_AK8_PF_mass_rho4Area","h1_AK8_PF_mass_rho4Area;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_rhoG4Area("h1_AK8_PF_mass_rhoG4Area","h1_AK8_PF_mass_rhoG4Area;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_rhom4Area("h1_AK8_PF_mass_rhom4Area","h1_AK8_PF_mass_rhom4Area;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_cleansingATLASjvf("h1_AK8_PF_mass_cleansingATLASjvf","h1_AK8_PF_mass_cleansingATLASjvf;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_cleansingATLASlin("h1_AK8_PF_mass_cleansingATLASlin","h1_AK8_PF_mass_cleansingATLASlin;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_cleansingATLASgau("h1_AK8_PF_mass_cleansingATLASgau","h1_AK8_PF_mass_cleansingATLASgau;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_cleansingCMSjvf("h1_AK8_PF_mass_cleansingCMSjvf","h1_AK8_PF_mass_cleansingCMSjvf;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_cleansingCMSlin("h1_AK8_PF_mass_cleansingCMSlin","h1_AK8_PF_mass_cleansingCMSlin;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	TH1D h1_AK8_PF_mass_cleansingCMSgau("h1_AK8_PF_mass_cleansingCMSgau","h1_AK8_PF_mass_cleansingCMSgau;jet mass;",nbin_mass,jetmass_min,jetmass_max);

	// reco mass VS gen mass
	TH2D h2_AK8_GEN_mass("h2_AK8_GEN_mass","h2_AK8_GEN_mass; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max); 
	TH2D h2_AK8_PF_mass_uncorr("h2_AK8_PF_mass_uncorr","h2_AK8_PF_mass_uncorr; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_rhoArea("h2_AK8_PF_mass_rhoArea","h2_AK8_PF_mass_rhoArea; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_rhoGArea("h2_AK8_PF_mass_rhoGArea","h2_AK8_PF_mass_rhoGArea; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_rho4Area("h2_AK8_PF_mass_rho4Area","h2_AK8_PF_mass_rho4Area; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_rhoG4Area("h2_AK8_PF_mass_rhoG4Area","h2_AK8_PF_mass_rhoG4Area; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_rhom4Area("h2_AK8_PF_mass_rhom4Area","h2_AK8_PF_mass_rhom4Area; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_cleansingATLASjvf("h2_AK8_PF_mass_cleansingATLASjvf","h2_AK8_PF_mass_cleansingATLASjvf; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_cleansingATLASlin("h2_AK8_PF_mass_cleansingATLASlin","h2_AK8_PF_mass_cleansingATLASlin; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_cleansingATLASgau("h2_AK8_PF_mass_cleansingATLASgau","h2_AK8_PF_mass_cleansingATLASgau; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_cleansingCMSjvf("h2_AK8_PF_mass_cleansingCMSjvf","h2_AK8_PF_mass_cleansingCMSjvf; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_cleansingCMSlin("h2_AK8_PF_mass_cleansingCMSlin","h2_AK8_PF_mass_cleansingCMSlin; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);
	TH2D h2_AK8_PF_mass_cleansingCMSgau("h2_AK8_PF_mass_cleansingCMSgau","h2_AK8_PF_mass_cleansingCMSgau; gen jet mass; reco jet mass",nbin_mass,jetmass_min,jetmass_max, nbin_mass,jetmass_min,jetmass_max);

	//calculate correlationFactor;
	Int_t num_points=0;
	TGraph gr_AK8_GEN_mass;
	TGraph gr_AK8_PF_mass_uncorr;
	TGraph gr_AK8_PF_mass_rhoArea;
	TGraph gr_AK8_PF_mass_rhoGArea;
	TGraph gr_AK8_PF_mass_rho4Area;
	TGraph gr_AK8_PF_mass_rhoG4Area;
	TGraph gr_AK8_PF_mass_rhom4Area;
	TGraph gr_AK8_PF_mass_cleansingATLASjvf;
	TGraph gr_AK8_PF_mass_cleansingATLASlin;
	TGraph gr_AK8_PF_mass_cleansingATLASgau;
	TGraph gr_AK8_PF_mass_cleansingCMSjvf;
	TGraph gr_AK8_PF_mass_cleansingCMSlin;
	TGraph gr_AK8_PF_mass_cleansingCMSgau;


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
			h1_AK8_PF.Fill(ratio);
			//AK8 uncorr
			ratio = tmp_AK8_PF_pt_uncorr/tmp_AK8_GEN_pt;
			h1_AK8_PF_uncorr.Fill(ratio);
			//AK8 rhoSW
			ratio = tmp_AK8_PF_pt_L1_rhoSW/tmp_AK8_GEN_pt;
			h1_AK8_PF_l1_rhoSW.Fill(ratio);
			//AK8 rhohand
			ratio = tmp_AK8_PF_pt_L1_rhoHand/tmp_AK8_GEN_pt;
			h1_AK8_PF_l1_rhoHand.Fill(ratio);
			//AK8 rhohand2
			ratio = tmp_AK8_PF_pt_L1_rhoHand2/tmp_AK8_GEN_pt;
			h1_AK8_PF_l1_rhoHand2.Fill(ratio);
			//AK8 rhogrid
			ratio = tmp_AK8_PF_pt_L1_rhoGrid/tmp_AK8_GEN_pt;
			h1_AK8_PF_l1_rhoGrid.Fill(ratio);

			//pt
			h1_AK8_GEN_pt.Fill(tmp_AK8_GEN_pt);
			h1_AK8_PF_pt_uncorr.Fill(tmp_AK8_PF_pt_uncorr);// pt of PF is with wrong JEC now
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
			//rho
			h1_AK8_GEN_rhoSW.Fill(tmp_AK8_GEN_rhoSW);
			h1_AK8_PF_rhoSW.Fill(tmp_AK8_PF_rhoSW);

			h1_AK8_GEN_rhoHand.Fill(tmp_AK8_GEN_rhoHand);
			h1_AK8_PF_rhoHand.Fill(tmp_AK8_PF_rhoHand);

			h1_AK8_GEN_rhoHand2.Fill(tmp_AK8_GEN_rhoHand2);
			h1_AK8_PF_rhoHand2.Fill(tmp_AK8_PF_rhoHand2);

			h1_AK8_GEN_rhoGrid.Fill(tmp_AK8_GEN_rhoGrid);
			h1_AK8_PF_rhoGrid.Fill(tmp_AK8_PF_rhoGrid);

			h1_area_PF.Fill(GroomedJet_AK8_PF_area[0]);  

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
			tmp_AK8_PF_mass_cleansingATLASjvf=GroomedJet_AK8_PF_mass_cleansingATLASjvf[0];
			tmp_AK8_PF_mass_cleansingATLASlin=GroomedJet_AK8_PF_mass_cleansingATLASlin[0];
			tmp_AK8_PF_mass_cleansingATLASgau=GroomedJet_AK8_PF_mass_cleansingATLASgau[0];
			tmp_AK8_PF_mass_cleansingCMSjvf=GroomedJet_AK8_PF_mass_cleansingCMSjvf[0];
			tmp_AK8_PF_mass_cleansingCMSlin=GroomedJet_AK8_PF_mass_cleansingCMSlin[0];
			tmp_AK8_PF_mass_cleansingCMSgau=GroomedJet_AK8_PF_mass_cleansingCMSgau[0];

			h1_AK8_GEN_mass.Fill(tmp_AK8_GEN_mass            ); 
			h1_AK8_PF_mass_uncorr.Fill(tmp_AK8_PF_mass_uncorr  );
			h1_AK8_PF_mass_rhoArea.Fill(tmp_AK8_PF_mass_rhoArea );
			h1_AK8_PF_mass_rhoGArea.Fill(tmp_AK8_PF_mass_rhoGArea);
			h1_AK8_PF_mass_rho4Area.Fill(tmp_AK8_PF_mass_rho4Area);
			h1_AK8_PF_mass_rhoG4Area.Fill(tmp_AK8_PF_mass_rhoG4Area);
			h1_AK8_PF_mass_rhom4Area.Fill(tmp_AK8_PF_mass_rhom4Area);
			h1_AK8_PF_mass_cleansingATLASjvf.Fill(tmp_AK8_PF_mass_cleansingATLASjvf);
			h1_AK8_PF_mass_cleansingATLASlin.Fill(tmp_AK8_PF_mass_cleansingATLASlin);
			h1_AK8_PF_mass_cleansingATLASgau.Fill(tmp_AK8_PF_mass_cleansingATLASgau);
			h1_AK8_PF_mass_cleansingCMSjvf.Fill(tmp_AK8_PF_mass_cleansingCMSjvf);
			h1_AK8_PF_mass_cleansingCMSlin.Fill(tmp_AK8_PF_mass_cleansingCMSlin);
			h1_AK8_PF_mass_cleansingCMSgau.Fill(tmp_AK8_PF_mass_cleansingCMSgau);

			h2_AK8_GEN_mass.Fill(tmp_AK8_GEN_mass, tmp_AK8_GEN_mass            ); 
			h2_AK8_PF_mass_uncorr.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_uncorr  );
			h2_AK8_PF_mass_rhoArea.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_rhoArea );
			h2_AK8_PF_mass_rhoGArea.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_rhoGArea);
			h2_AK8_PF_mass_rho4Area.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_rho4Area);
			h2_AK8_PF_mass_rhoG4Area.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_rhoG4Area);
			h2_AK8_PF_mass_rhom4Area.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_rhom4Area);
			h2_AK8_PF_mass_cleansingATLASjvf.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_cleansingATLASjvf);
			h2_AK8_PF_mass_cleansingATLASlin.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_cleansingATLASlin);
			h2_AK8_PF_mass_cleansingATLASgau.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_cleansingATLASgau);
			h2_AK8_PF_mass_cleansingCMSjvf.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_cleansingCMSjvf);
			h2_AK8_PF_mass_cleansingCMSlin.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_cleansingCMSlin);
			h2_AK8_PF_mass_cleansingCMSgau.Fill(tmp_AK8_GEN_mass, tmp_AK8_PF_mass_cleansingCMSgau);


			gr_AK8_GEN_mass.SetPoint(num_points, tmp_AK8_GEN_mass, tmp_AK8_GEN_mass            ); 
			gr_AK8_PF_mass_uncorr.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_uncorr  );
			gr_AK8_PF_mass_rhoArea.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_rhoArea );
			gr_AK8_PF_mass_rhoGArea.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_rhoGArea);
			gr_AK8_PF_mass_rho4Area.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_rho4Area);
			gr_AK8_PF_mass_rhoG4Area.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_rhoG4Area);
			gr_AK8_PF_mass_rhom4Area.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_rhom4Area);
			gr_AK8_PF_mass_cleansingATLASjvf.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_cleansingATLASjvf);
			gr_AK8_PF_mass_cleansingATLASlin.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_cleansingATLASlin);
			gr_AK8_PF_mass_cleansingATLASgau.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_cleansingATLASgau);
			gr_AK8_PF_mass_cleansingCMSjvf.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_cleansingCMSjvf);
			gr_AK8_PF_mass_cleansingCMSlin.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_cleansingCMSlin);
			gr_AK8_PF_mass_cleansingCMSgau.SetPoint(num_points,tmp_AK8_GEN_mass, tmp_AK8_PF_mass_cleansingCMSgau);
			num_points++;

		}
		//PF matching efficiency
		dr = TMath::Sqrt( (gen_jet_eta-tmp_AK8_PF_eta)*(gen_jet_eta-tmp_AK8_PF_eta) + (gen_jet_phi-tmp_AK8_PF_phi)*(gen_jet_phi-tmp_AK8_PF_phi) ) ;
		h1_PF_match.Fill(dr);
	}


	/*Draw_and_Save(h1_AK8_PF);
	Draw_and_Save(h1_AK8_PF_uncorr);
	Draw_and_Save(h1_AK8_PF_l1_rhoSW);
	Draw_and_Save(h1_AK8_PF_l1_rhoHand);
	Draw_and_Save(h1_AK8_PF_l1_rhoHand2);
	Draw_and_Save(h1_AK8_PF_l1_rhoGrid);*/

	Draw_and_Save(h1_area_PF);

	Draw_and_Save(h1_nPV);

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
	Draw_and_Save(h1_AK8_PF_eta);
	Draw_and_Save(h1_AK8_PF_zjet_dr);
	Draw_and_Save(h1_AK8_PF_zjet_dphi);
	/*Draw_and_Save(h1_AK8_PF_rhoSW);
	Draw_and_Save(h1_AK8_PF_rhoHand);
	Draw_and_Save(h1_AK8_PF_rhoHand2);
	Draw_and_Save(h1_AK8_PF_rhoGrid);
	Draw_and_Save(h2_AK8_PF_rhoSW_vs_nPV);
	Draw_and_Save(h2_AK8_PF_rhoHand_vs_nPV);
	Draw_and_Save(h2_AK8_PF_rhoHand2_vs_nPV);
	Draw_and_Save(h2_AK8_PF_rhoGrid_vs_nPV);*/

	Draw_and_Save(h1_PF_match);

	Draw_and_Save(h1_AK8_GEN_mass                );
	Draw_and_Save(h1_AK8_PF_mass_uncorr      );
	Draw_and_Save(h1_AK8_PF_mass_rhoArea     );
	Draw_and_Save(h1_AK8_PF_mass_rhoGArea    );
	Draw_and_Save(h1_AK8_PF_mass_rho4Area    );
	Draw_and_Save(h1_AK8_PF_mass_rhoG4Area   );
	Draw_and_Save(h1_AK8_PF_mass_rhom4Area   );
	Draw_and_Save(h1_AK8_PF_mass_cleansingATLASjvf   );
	Draw_and_Save(h1_AK8_PF_mass_cleansingATLASlin   );
	Draw_and_Save(h1_AK8_PF_mass_cleansingATLASgau   );
	Draw_and_Save(h1_AK8_PF_mass_cleansingCMSjvf   );
	Draw_and_Save(h1_AK8_PF_mass_cleansingCMSlin   );
	Draw_and_Save(h1_AK8_PF_mass_cleansingCMSgau   );

	Draw_and_Save(h2_AK8_GEN_mass, Form("%g correlated", gr_AK8_GEN_mass.GetCorrelationFactor())   );
	Draw_and_Save(h2_AK8_PF_mass_uncorr				, Form("%g correlated", gr_AK8_PF_mass_uncorr.GetCorrelationFactor())				 );
	Draw_and_Save(h2_AK8_PF_mass_rhoArea			, Form("%g correlated", gr_AK8_PF_mass_rhoArea.GetCorrelationFactor())				);
	Draw_and_Save(h2_AK8_PF_mass_rhoGArea			, Form("%g correlated", gr_AK8_PF_mass_rhoGArea.GetCorrelationFactor())				);
	Draw_and_Save(h2_AK8_PF_mass_rho4Area			, Form("%g correlated", gr_AK8_PF_mass_rho4Area.GetCorrelationFactor())				);
	Draw_and_Save(h2_AK8_PF_mass_rhoG4Area			, Form("%g correlated", gr_AK8_PF_mass_rhoG4Area.GetCorrelationFactor())				);
	Draw_and_Save(h2_AK8_PF_mass_rhom4Area			, Form("%g correlated", gr_AK8_PF_mass_rhom4Area.GetCorrelationFactor())				);
	Draw_and_Save(h2_AK8_PF_mass_cleansingATLASjvf	, Form("%g correlated", gr_AK8_PF_mass_cleansingATLASjvf.GetCorrelationFactor())		);
	Draw_and_Save(h2_AK8_PF_mass_cleansingATLASlin	, Form("%g correlated", gr_AK8_PF_mass_cleansingATLASlin.GetCorrelationFactor())		);
	Draw_and_Save(h2_AK8_PF_mass_cleansingATLASgau	, Form("%g correlated", gr_AK8_PF_mass_cleansingATLASgau.GetCorrelationFactor())		);
	Draw_and_Save(h2_AK8_PF_mass_cleansingCMSjvf	, Form("%g correlated", gr_AK8_PF_mass_cleansingCMSjvf.GetCorrelationFactor())		);
	Draw_and_Save(h2_AK8_PF_mass_cleansingCMSlin	, Form("%g correlated", gr_AK8_PF_mass_cleansingCMSlin.GetCorrelationFactor())		);
	Draw_and_Save(h2_AK8_PF_mass_cleansingCMSgau	, Form("%g correlated", gr_AK8_PF_mass_cleansingCMSgau.GetCorrelationFactor())		);
}
