#define MyClass_cxx
#include "MyClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TMath.h"

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
void Draw_and_Save(TH2D h2){
	h2.Write();
	TCanvas *c1 = new TCanvas(Form("c1_%s",h2.GetTitle()),Form("c1_%s",h2.GetTitle()),200,10,600,600);
	c1->cd();
	h2.Draw();
	//c1->Print(Form("%s.pdf",h2.GetTitle()));
	c1->Print(Form("%s.png",h2.GetTitle()));
	delete c1;
}

void MyClass::Loop()
{
	//   In a ROOT session, you can do:
	//      Root > .L MyClass.C
	//      Root > MyClass t
	//      Root > t.GetEntry(12); // Fill t data members with entry number 12
	//      Root > t.Show();       // Show values of entry 12
	//      Root > t.Show(16);     // Read and show values of entry 16
	//      Root > t.Loop();       // Loop on all entries
	//

	//     This is the loop skeleton where:
	//    jentry is the global entry number in the chain
	//    ientry is the entry number in the current Tree
	//  Note that the argument to GetEntry must be:
	//    jentry for TChain::GetEntry
	//    ientry for TTree::GetEntry and TBranch::GetEntry
	//
	//       To read only selected branches, Insert statements like:
	// METHOD1:
	//    fChain->SetBranchStatus("*",0);  // disable all branches
	//    fChain->SetBranchStatus("branchname",1);  // activate branchname
	// METHOD2: replace line
	//    fChain->GetEntry(jentry);       //read all branches
	//by  b_branchname->GetEntry(ientry); //read only this branch
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();

	Double_t tmp_GEN_pt=0.;
	Double_t tmp_GEN_eta=0.;
	Double_t tmp_GEN_phi=0.;
	Double_t tmp_GEN_rhoSW=0,;
	Double_t tmp_GEN_rhoHand=0,;
	Double_t tmp_GEN_rhoHand2=0,;
	Double_t tmp_GEN_rhoGrid=0,;
	
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
	Double_t tmp_AK5_PF_rhoSW=0,;
	Double_t tmp_AK5_PF_rhoHand=0,;
	Double_t tmp_AK5_PF_rhoHand2=0,;
	Double_t tmp_AK5_PF_rhoGrid=0,;

	Double_t tmp_AK5_PFCHS_pt=0.;
	Double_t tmp_AK5_PFCHS_eta=0.;
	Double_t tmp_AK5_PFCHS_phi=0.;
	Double_t tmp_AK5_PFCHS_pt_uncorr=0.;
	Double_t tmp_AK5_PFCHS_pt_L1_rhoSW=0.;
	Double_t tmp_AK5_PFCHS_pt_L1_rhoHand=0.;
	Double_t tmp_AK5_PFCHS_pt_L1_rhoHand2=0.;
	Double_t tmp_AK5_PFCHS_pt_L1_rhoGrid=0.;
	Double_t tmp_AK5_PFCHS_rhoSW=0,;
	Double_t tmp_AK5_PFCHS_rhoHand=0,;
	Double_t tmp_AK5_PFCHS_rhoHand2=0,;
	Double_t tmp_AK5_PFCHS_rhoGrid=0,;

	Double_t tmp_Z_eta=0.;
	Double_t tmp_Z_phi=0.;

	Double_t tmp_event_nPV=0.;

	Double_t ratio=0.;
	Double_t dr=0,; // Delta R
	Double_t dphi=0.; // Delta Phi

	//jet mass
	Double_t tmp_GEN_mass=0.;
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

	TH1D h1_GEN_pt("h1_GEN_pt","h1_GEN_pt",50,0,200); h1_GEN_pt.SetLineColor(kRed);
	TH1D h1_GEN_eta("h1_GEN_eta","h1_GEN_eta",50,-2.5,2.5); h1_GEN_eta.SetLineColor(kRed);
	TH1D h1_GEN_phi("h1_GEN_phi","h1_GEN_phi",50,-4,4); h1_GEN_phi.SetLineColor(kRed);
	TH1D h1_GEN_zjet_dr("h1_GEN_zjet_dr","h1_GEN_zjet_dr",50,0,10); h1_GEN_zjet_dr.SetLineColor(kRed);
	TH1D h1_GEN_zjet_dphi("h1_GEN_zjet_dphi","h1_GEN_zjet_dphi",50,0,5); h1_GEN_zjet_dphi.SetLineColor(kRed);
	TH1D h1_GEN_rhoSW("h1_GEN_rhoSW","h1_GEN_rhoSW",50,rhomin,rhomax); h1_GEN_rhoSW.SetLineColor(kRed);
	TH1D h1_GEN_rhoHand("h1_GEN_rhoHand","h1_GEN_rhoHand",50,rhomin,rhomax); h1_GEN_rhoHand.SetLineColor(kRed);
	TH1D h1_GEN_rhoHand2("h1_GEN_rhoHand2","h1_GEN_rhoHand2",50,rhomin,rhomax); h1_GEN_rhoHand2.SetLineColor(kRed);
	TH1D h1_GEN_rhoGrid("h1_GEN_rhoGrid","h1_GEN_rhoGrid",50,rhomin,rhomax); h1_GEN_rhoGrid.SetLineColor(kRed);
	TH2D h2_GEN_rhoSW_vs_nPV("h2_GEN_rhoSW_vs_nPV","h2_GEN_rhoSW_vs_nPV",50,rhomin,rhomax,50,0,50); h2_GEN_rhoSW_vs_nPV.SetMarkerColor(kRed);
	TH2D h2_GEN_rhoHand_vs_nPV("h2_GEN_rhoHand_vs_nPV","h2_GEN_rhoHand_vs_nPV",50,rhomin,rhomax,50,0,50); h2_GEN_rhoHand_vs_nPV.SetMarkerColor(kRed);
	TH2D h2_GEN_rhoHand2_vs_nPV("h2_GEN_rhoHand2_vs_nPV","h2_GEN_rhoHand2_vs_nPV",50,rhomin,rhomax,50,0,50); h2_GEN_rhoHand2_vs_nPV.SetMarkerColor(kRed);
	TH2D h2_GEN_rhoGrid_vs_nPV("h2_GEN_rhoGrid_vs_nPV","h2_GEN_rhoGrid_vs_nPV",50,rhomin,rhomax,50,0,50); h2_GEN_rhoGrid_vs_nPV.SetMarkerColor(kRed);

	TH1D h1_ak5_pf_pt_uncorr("h1_ak5_pf_pt_uncorr","h1_ak5_pf_pt_uncorr",50,0,200);
	TH1D h1_ak5_pf_eta("h1_ak5_pf_eta","h1_ak5_pf_eta",50,-2.5,2.5);
	TH1D h1_ak5_pf_phi("h1_ak5_pf_phi","h1_ak5_pf_phi",50,-4,4);
	TH1D h1_ak5_pf_zjet_dr("h1_ak5_pf_zjet_dr","h1_ak5_pf_zjet_dr",50,0,10);
	TH1D h1_ak5_pf_zjet_dphi("h1_ak5_pf_zjet_dphi","h1_ak5_pf_zjet_dphi",50,0,5);
	TH1D h1_ak5_pf_rhoSW("h1_ak5_pf_rhoSW","h1_ak5_pf_rhoSW",50,rhomin,rhomax);
	TH1D h1_ak5_pf_rhoHand("h1_ak5_pf_rhoHand","h1_ak5_pf_rhoHand",50,rhomin,rhomax);
	TH1D h1_ak5_pf_rhoHand2("h1_ak5_pf_rhoHand2","h1_ak5_pf_rhoHand2",50,rhomin,rhomax);
	TH1D h1_ak5_pf_rhoGrid("h1_ak5_pf_rhoGrid","h1_ak5_pf_rhoGrid",50,rhomin,rhomax);
	TH2D h2_ak5_pf_rhoSW_vs_nPV("h2_ak5_pf_rhoSW_vs_nPV","h2_ak5_pf_rhoSW_vs_nPV",50,rhomin,rhomax,50,0,50);
	TH2D h2_ak5_pf_rhoHand_vs_nPV("h2_ak5_pf_rhoHand_vs_nPV","h2_ak5_pf_rhoHand_vs_nPV",50,rhomin,rhomax,50,0,50);
	TH2D h2_ak5_pf_rhoHand2_vs_nPV("h2_ak5_pf_rhoHand2_vs_nPV","h2_ak5_pf_rhoHand2_vs_nPV",50,rhomin,rhomax,50,0,50);
	TH2D h2_ak5_pf_rhoGrid_vs_nPV("h2_ak5_pf_rhoGrid_vs_nPV","h2_ak5_pf_rhoGrid_vs_nPV",50,rhomin,rhomax,50,0,50);


	TH1D h1_ak5_pfchs_pt_uncorr("h1_ak5_pfchs_pt_uncorr","h1_ak5_pfchs_pt_uncorr",50,0,200);
	TH1D h1_ak5_pfchs_eta("h1_ak5_pfchs_eta","h1_ak5_pfchs_eta",50,-2.5,2.5);
	TH1D h1_ak5_pfchs_phi("h1_ak5_pfchs_phi","h1_ak5_pfchs_phi",50,-4,4);
	TH1D h1_ak5_pfchs_zjet_dr("h1_ak5_pfchs_zjet_dr","h1_ak5_pfchs_zjet_dr",50,0,10);
	TH1D h1_ak5_pfchs_zjet_dphi("h1_ak5_pfchs_zjet_dphi","h1_ak5_pfchs_zjet_dphi",50,0,5);
	TH1D h1_ak5_pfchs_rhoSW("h1_ak5_pfchs_rhoSW","h1_ak5_pfchs_rhoSW",50,rhomin,rhomax);
	TH1D h1_ak5_pfchs_rhoHand("h1_ak5_pfchs_rhoHand","h1_ak5_pfchs_rhoHand",50,rhomin,rhomax);
	TH1D h1_ak5_pfchs_rhoHand2("h1_ak5_pfchs_rhoHand2","h1_ak5_pfchs_rhoHand2",50,rhomin,rhomax);
	TH1D h1_ak5_pfchs_rhoGrid("h1_ak5_pfchs_rhoGrid","h1_ak5_pfchs_rhoGrid",50,rhomin,rhomax);
	TH2D h2_ak5_pfchs_rhoSW_vs_nPV("h2_ak5_pfchs_rhoSW_vs_nPV","h2_ak5_pfchs_rhoSW_vs_nPV",50,rhomin,rhomax,50,0,50);
	TH2D h2_ak5_pfchs_rhoHand_vs_nPV("h2_ak5_pfchs_rhoHand_vs_nPV","h2_ak5_pfchs_rhoHand_vs_nPV",50,rhomin,rhomax,50,0,50);
	TH2D h2_ak5_pfchs_rhoHand2_vs_nPV("h2_ak5_pfchs_rhoHand2_vs_nPV","h2_ak5_pfchs_rhoHand2_vs_nPV",50,rhomin,rhomax,50,0,50);
	TH2D h2_ak5_pfchs_rhoGrid_vs_nPV("h2_ak5_pfchs_rhoGrid_vs_nPV","h2_ak5_pfchs_rhoGrid_vs_nPV",50,rhomin,rhomax,50,0,50);

	TH2D h2_ak5_pfchs_ratioHand_vs_nPV("h2_ak5_pfchs_ratioHand_vs_nPV","h2_ak5_pfchs_ratioHand_vs_nPV",50,0,50, 50, 0, 2.5);
	TH2D h2_ak5_pfchs_ratioHand_vs_ptHand("h2_ak5_pfchs_ratioHand_vs_ptHand","h2_ak5_pfchs_ratioHand_vs_ptHand",50,0,200, 50, 0, 2.5);
	TH2D h2_ak5_pfchs_ratioHand_vs_eta("h2_ak5_pfchs_ratioHand_vs_eta","h2_ak5_pfchs_ratioHand_vs_eta", 50, -2.5, 2.5, 50, 0, 2.5);



	TH1D h1_PF_match("h1_PF_match","h1_PF_match",50,0,1.); h1_PF_match.SetLineColor(kRed);
	TH1D h1_PFCHS_match("h1_PFCHS_match","h1_PFCHS_match",50,0,1.);

	// JetPt/GenPt
	TH1D h1_PFCor("h1_PFCor","h1_PFCor",50,0,2.5);
	TH1D h1_PFCor_uncorr("h1_PFCor_uncorr","h1_PFCor_uncorr",50,0,2.5);
	TH1D h1_PFCor_afterL1("h1_PFCor_afterL1","h1_PFCor_afterL1",50,0,2.5);
	TH1D h1_PFCor_afterL2("h1_PFCor_afterL2","h1_PFCor_afterL2",50,0,2.5);

	TH1D h1_ak5_pf("h1_ak5_pf","h1_ak5_pf",50,0,2.5);
	TH1D h1_ak5_pf_uncorr("h1_ak5_pf_uncorr","h1_ak5_pf_uncorr",50,0,2.5);
	TH1D h1_ak5_pf_l1_rhoSW("h1_ak5_pf_l1_rhoSW","h1_ak5_pf_l1_rhoSW",50,0,2.5);
	TH1D h1_ak5_pf_l1_rhoHand("h1_ak5_pf_l1_rhoHand","h1_ak5_pf_l1_rhoHand",50,0,2.5);
	TH1D h1_ak5_pf_l1_rhoHand2("h1_ak5_pf_l1_rhoHand2","h1_ak5_pf_l1_rhoHand2",50,0,2.5);
	TH1D h1_ak5_pf_l1_rhoGrid("h1_ak5_pf_l1_rhoGrid","h1_ak5_pf_l1_rhoGrid",50,0,2.5);

	TH1D h1_ak5_pfchs("h1_ak5_pfchs","h1_ak5_pfchs",50,0,2.5);
	TH1D h1_ak5_pfchs_uncorr("h1_ak5_pfchs_uncorr","h1_ak5_pfchs_uncorr",50,0,2.5);
	TH1D h1_ak5_pfchs_l1_rhoSW("h1_ak5_pfchs_l1_rhoSW","h1_ak5_pfchs_l1_rhoSW",50,0,2.5);
	TH1D h1_ak5_pfchs_l1_rhoHand("h1_ak5_pfchs_l1_rhoHand","h1_ak5_pfchs_l1_rhoHand",50,0,2.5);
	TH1D h1_ak5_pfchs_l1_rhoHand2("h1_ak5_pfchs_l1_rhoHand2","h1_ak5_pfchs_l1_rhoHand2",50,0,2.5);
	TH1D h1_ak5_pfchs_l1_rhoGrid("h1_ak5_pfchs_l1_rhoGrid","h1_ak5_pfchs_l1_rhoGrid",50,0,2.5);

	//area
	TH1D h1_area_PFCor("h1_area_PFCor","h1_area_PFCor",50,0.6,1.1);
	TH1D h1_area_pf("h1_area_pf","h1_area_pf",50,0.6,1.1);
	TH1D h1_area_pfchs("h1_area_pfchs","h1_area_pfchs",50,0.6,1.1);

	//rho
	TH1D h1_rho_PFCor("h1_rho_PFCor","h1_rho_PFCor",30,0,30);


	//mass
	int nbin_mass=20;double jetmass_min=0;double jetmass_max=20.;
	TH1D h1_GEN_mass("h1_GEN_mass","h1_GEN_mass;jet mass;",nbin_mass,jetmass_min,jetmass_max);
	h1_GEN_mass.SetLineColor(kRed);
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


	/*	TH1D h1_rhoSW_pf("h1_rhoSW_pf","h1_rhoSW_pf",30,0,30);
		TH1D h1_rhoHand_pf("h1_rhoHand_pf","h1_rhoHand_pf",30,0,30);
		TH1D h1_rhoHand2_pf("h1_rhoHand2_pf","h1_rhoHand2_pf",30,0,30);
		TH1D h1_rhoGrid_pf("h1_rhoGrid_pf","h1_rhoGrid_pf",30,0,30);

		TH1D h1_rhoSW_pfchs("h1_rhoSW_pfchs","h1_rhoSW_pfchs",30,0,30);
		TH1D h1_rhoHand_pfchs("h1_rhoHand_pfchs","h1_rhoHand_pfchs",30,0,30);
		TH1D h1_rhoHand2_pfchs("h1_rhoHand2_pfchs","h1_rhoHand2_pfchs",30,0,30);
		TH1D h1_rhoGrid_pfchs("h1_rhoGrid_pfchs","h1_rhoGrid_pfchs",30,0,30);
		*/
	// For GEN-RECO matching
	Double_t gen_jet_eta=gen_jet_phi=0.;
	Double_t pfchs_eta=pfchs_jet_phi=0.;
	Double_t PFCor_jet_eta=PFCor_jet_phi=0.;

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;

		gen_jet_eta=GenGroomedJet_AK5_GEN_eta[0];
		gen_jet_phi=GenGroomedJet_AK5_GEN_phi[0];

		pfchs_eta=GroomedJet_AK5_PFCHS_eta[0];
		pfchs_phi=GroomedJet_AK5_PFCHS_phi[0];

		//Gen Jet matching with pfchs_uncerr
		if( TMath::Sqrt( (gen_jet_eta-pfchs_eta)*(gen_jet_eta-pfchs_eta) + (gen_jet_phi-pfchs_phi)*(gen_jet_phi-pfchs_phi) ) <0.3 )
		{
			int i_PFCorJet_matching_pfchs=0;//_matching with leading pfchs_uncerr
			for(int i_PFCorJet_matching_pfchs=0;i_PFCorJet_matching_pfchs<numPFCorJets;i_PFCorJet_matching_pfchs++){
				PFCor_jet_eta=JetPFCor_Eta[i_PFCorJet_matching_pfchs]; PFCor_jet_phi=JetPFCor_Phi[i_PFCorJet_matching_pfchs];
				if ( equal(pfchs_eta, PFCor_jet_eta) && equal(pfchs_phi, PFCor_jet_phi) )break;
			}

			tmp_GEN_pt = GenGroomedJet_AK5_GEN_pt[0];
			tmp_GEN_eta = GenGroomedJet_AK5_GEN_eta[0];
			tmp_GEN_phi = GenGroomedJet_AK5_GEN_phi[0];

			tmp_PFCor_pt = JetPFCor_Pt[i_PFCorJet_matching_pfchs];
			tmp_PFCor_Pt_uncorr = JetPFCor_Pt_uncorr[i_PFCorJet_matching_pfchs];
			tmp_PFCor_Pt_afterL1 = JetPFCor_Pt_afterL1[i_PFCorJet_matching_pfchs];
			tmp_PFCor_Pt_afterL2 = JetPFCor_Pt_afterL2[i_PFCorJet_matching_pfchs];

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

			tmp_Z_eta = Z_eta;
			tmp_Z_phi = Z_phi;

			tmp_event_nPV = event_nPV[0];

			tmp_GEN_rhoSW = GenGroomedJet_AK5_GEN_rhoSW[0];
			tmp_GEN_rhoHand = GenGroomedJet_AK5_GEN_rhohand[0];
			tmp_GEN_rhoHand2 = GenGroomedJet_AK5_GEN_rhohand2[0];
			tmp_GEN_rhoGrid = GenGroomedJet_AK5_GEN_rhogrid[0];

			tmp_AK5_PF_rhoSW = GroomedJet_AK5_PF_rhoSW[0];
			tmp_AK5_PF_rhoHand = GroomedJet_AK5_PF_rhohand[0];
			tmp_AK5_PF_rhoHand2 = GroomedJet_AK5_PF_rhohand2[0];
			tmp_AK5_PF_rhoGrid = GroomedJet_AK5_PF_rhogrid[0];

			tmp_AK5_PFCHS_rhoSW = GroomedJet_AK5_PFCHS_rhoSW[0];
			tmp_AK5_PFCHS_rhoHand = GroomedJet_AK5_PFCHS_rhohand[0];
			tmp_AK5_PFCHS_rhoHand2 = GroomedJet_AK5_PFCHS_rhohand2[0];
			tmp_AK5_PFCHS_rhoGrid = GroomedJet_AK5_PFCHS_rhogrid[0];

			//============= begin to fill hist ==============
			//PFCor
			ratio = tmp_PFCor_pt/tmp_GEN_pt;
			h1_PFCor.Fill(ratio);
			//PFCor uncorr
			ratio = tmp_PFCor_Pt_uncorr/tmp_GEN_pt;
			h1_PFCor_uncorr.Fill(ratio);
			//PFCor L1
			ratio = tmp_PFCor_Pt_afterL1/tmp_GEN_pt;
			h1_PFCor_afterL1.Fill(ratio);
			//PFCor L2
			ratio = tmp_PFCor_Pt_afterL2/tmp_GEN_pt;
			h1_PFCor_afterL2.Fill(ratio);

			//ak5 PF
			ratio = tmp_AK5_PF_pt/tmp_GEN_pt;
			h1_ak5_pf.Fill(ratio);
			//ak5 uncorr
			ratio = tmp_AK5_PF_pt_uncorr/tmp_GEN_pt;
			h1_ak5_pf_uncorr.Fill(ratio);
			//ak5 rhoSW
			ratio = tmp_AK5_PF_pt_L1_rhoSW/tmp_GEN_pt;
			h1_ak5_pf_l1_rhoSW.Fill(ratio);
			//ak5 rhohand
			ratio = tmp_AK5_PF_pt_L1_rhoHand/tmp_GEN_pt;
			h1_ak5_pf_l1_rhoHand.Fill(ratio);


			h2_ak5_pfchs_ratioHand_vs_nPV.Fill(tmp_event_nPV, ratio);
			h2_ak5_pfchs_ratioHand_vs_ptHand.Fill(tmp_AK5_PF_pt_L1_rhoHand, ratio);
			h2_ak5_pfchs_ratioHand_vs_eta.Fill(tmp_AK5_PFCHS_eta, ratio );


			//ak5 rhohand2
			ratio = tmp_AK5_PF_pt_L1_rhoHand2/tmp_GEN_pt;
			h1_ak5_pf_l1_rhoHand2.Fill(ratio);
			//ak5 rhogrid
			ratio = tmp_AK5_PF_pt_L1_rhoGrid/tmp_GEN_pt;
			h1_ak5_pf_l1_rhoGrid.Fill(ratio);

			//ak5 PFCHS
			ratio = tmp_AK5_PFCHS_pt/tmp_GEN_pt;
			h1_ak5_pfchs.Fill(ratio);

			ratio = tmp_AK5_PFCHS_pt_uncorr/tmp_GEN_pt;
			h1_ak5_pfchs_uncorr.Fill(ratio);

			ratio = tmp_AK5_PFCHS_pt_L1_rhoSW/tmp_GEN_pt;
			h1_ak5_pfchs_l1_rhoSW.Fill(ratio);

			ratio = tmp_AK5_PFCHS_pt_L1_rhoHand/tmp_GEN_pt;
			h1_ak5_pfchs_l1_rhoHand.Fill(ratio);

			ratio = tmp_AK5_PFCHS_pt_L1_rhoHand2/tmp_GEN_pt;
			h1_ak5_pfchs_l1_rhoHand2.Fill(ratio);

			ratio = tmp_AK5_PFCHS_pt_L1_rhoGrid/tmp_GEN_pt;
			h1_ak5_pfchs_l1_rhoGrid.Fill(ratio);


			//pt
			h1_GEN_pt.Fill(tmp_GEN_pt);
			h1_ak5_pf_pt_uncorr.Fill(tmp_AK5_PF_pt_uncorr);// pt of PF and PFCHS are with wrong JEC now
			h1_ak5_pfchs_pt_uncorr.Fill(tmp_AK5_PFCHS_pt_uncorr);
			//eta
			h1_GEN_eta.Fill(tmp_GEN_eta);
			h1_ak5_pf_eta.Fill(tmp_AK5_PF_eta);
			h1_ak5_pfchs_eta.Fill(tmp_AK5_PFCHS_eta);
			//phi
			h1_GEN_phi.Fill(tmp_GEN_phi);
			h1_ak5_pf_phi.Fill(tmp_AK5_PF_phi);
			h1_ak5_pfchs_phi.Fill(tmp_AK5_PFCHS_phi);
			//dr
			dr = TMath::Sqrt( (Z_eta-tmp_GEN_eta)*(Z_eta-tmp_GEN_eta) + (Z_phi-tmp_GEN_phi)*(Z_phi-tmp_GEN_phi) );
			h1_GEN_zjet_dr.Fill(dr);
			dr = TMath::Sqrt( (Z_eta-tmp_AK5_PF_eta)*(Z_eta-tmp_AK5_PF_eta) + (Z_phi-tmp_AK5_PF_phi)*(Z_phi-tmp_AK5_PF_phi) );
			h1_ak5_pf_zjet_dr.Fill(dr);
			dr = TMath::Sqrt( (Z_eta-tmp_AK5_PFCHS_eta)*(Z_eta-tmp_AK5_PFCHS_eta) + (Z_phi-tmp_AK5_PFCHS_phi)*(Z_phi-tmp_AK5_PFCHS_phi) );
			h1_ak5_pfchs_zjet_dr.Fill(dr);
			//dphi
			dphi = TMath::Sqrt( (Z_phi-tmp_GEN_phi)*(Z_phi-tmp_GEN_phi) );
			h1_GEN_zjet_dphi.Fill(dphi);
			dphi = TMath::Sqrt( (Z_phi-tmp_AK5_PF_phi)*(Z_phi-tmp_AK5_PF_phi) );
			h1_ak5_pf_zjet_dphi.Fill(dphi);
			dphi = TMath::Sqrt( (Z_phi-tmp_AK5_PFCHS_phi)*(Z_phi-tmp_AK5_PFCHS_phi) );
			h1_ak5_pfchs_zjet_dphi.Fill(dphi);
			//nPV
			h1_nPV.Fill(tmp_event_nPV);
			//rho
			h1_GEN_rhoSW.Fill(tmp_GEN_rhoSW);
			h1_ak5_pf_rhoSW.Fill(tmp_AK5_PF_rhoSW);
			h1_ak5_pfchs_rhoSW.Fill(tmp_AK5_PFCHS_rhoSW);

			h1_GEN_rhoHand.Fill(tmp_GEN_rhoHand);
			h1_ak5_pf_rhoHand.Fill(tmp_AK5_PF_rhoHand);
			h1_ak5_pfchs_rhoHand.Fill(tmp_AK5_PFCHS_rhoHand);

			h1_GEN_rhoHand2.Fill(tmp_GEN_rhoHand2);
			h1_ak5_pf_rhoHand2.Fill(tmp_AK5_PF_rhoHand2);
			h1_ak5_pfchs_rhoHand2.Fill(tmp_AK5_PFCHS_rhoHand2);

			h1_GEN_rhoGrid.Fill(tmp_GEN_rhoGrid);
			h1_ak5_pf_rhoGrid.Fill(tmp_AK5_PF_rhoGrid);
			h1_ak5_pfchs_rhoGrid.Fill(tmp_AK5_PFCHS_rhoGrid);

			h1_area_PFCor.Fill(JetPFCor_Area[i_PFCorJet_matching_pfchs]);
			h1_area_pf.Fill(GroomedJet_AK5_PF_area[0]);  
			h1_area_pfchs.Fill(GroomedJet_AK5_PFCHS_area[0]); 

			/*//h1_rho_PFCor.Fill(); 
			  h1_rhoSW_pf.Fill(GroomedJet_AK5_PF_rhoSW); 
			  h1_rhoHand_pf.Fill(GroomedJet_AK5_PF_rhohand); 
			  h1_rhoHand2_pf.Fill(GroomedJet_AK5_PF_rhohand2); 
			  h1_rhoGrid_pf.Fill(GroomedJet_AK5_PF_rhogrid); 
			  h1_rhoSW_pfchs.Fill(GroomedJet_AK5_PFCHS_rhoSW); 
			  h1_rhoHand_pfchs.Fill(GroomedJet_AK5_PFCHS_rhohand); 
			  h1_rhoHand2_pfchs.Fill(GroomedJet_AK5_PFCHS_rhohand2); 
			  h1_rhoGrid_pfchs.Fill(GroomedJet_AK5_PFCHS_rhogrid); 
			  */

			//rho vs nPV
			h2_GEN_rhoSW_vs_nPV.Fill(tmp_GEN_rhoSW,tmp_event_nPV);
			h2_GEN_rhoHand_vs_nPV.Fill(tmp_GEN_rhoHand,tmp_event_nPV);
			h2_GEN_rhoHand2_vs_nPV.Fill(tmp_GEN_rhoHand2,tmp_event_nPV);
			h2_GEN_rhoGrid_vs_nPV.Fill(tmp_GEN_rhoGrid,tmp_event_nPV);

			h2_ak5_pf_rhoSW_vs_nPV.Fill(tmp_AK5_PF_rhoSW,tmp_event_nPV);
			h2_ak5_pf_rhoHand_vs_nPV.Fill(tmp_AK5_PF_rhoHand,tmp_event_nPV);
			h2_ak5_pf_rhoHand2_vs_nPV.Fill(tmp_AK5_PF_rhoHand2,tmp_event_nPV);
			h2_ak5_pf_rhoGrid_vs_nPV.Fill(tmp_AK5_PF_rhoGrid,tmp_event_nPV);
			//ratio vs pt
			//PF vs PFCHS matching efficiency
			dr = TMath::Sqrt( (gen_jet_eta-tmp_AK5_PF_eta)*(gen_jet_eta-tmp_AK5_PF_eta) + (gen_jet_phi-tmp_AK5_PF_phi)*(gen_jet_phi-tmp_AK5_PF_phi) ) ;
			h1_PF_match.Fill(dr);
			dr = TMath::Sqrt( (gen_jet_eta-tmp_AK5_PFCHS_eta)*(gen_jet_eta-tmp_AK5_PFCHS_eta) + (gen_jet_phi-tmp_AK5_PFCHS_phi)*(gen_jet_phi-tmp_AK5_PFCHS_phi) ) ;       
			h1_PFCHS_match.Fill(dr);

			tmp_GEN_mass=GenGroomedJet_AK5_GEN_mass_uncorr[0];
			tmp_PFCor_mass=JetPFCor_Mass[i_PFCorJet_matching_pfchs];
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

			h1_GEN_mass.Fill(tmp_GEN_mass            ); 
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
	h1_GEN_pt->Draw();
	h1_ak5_pf_pt_uncorr->Draw("same");
	c1->Print("RECO vs_GEN_Pt.png");

	TCanvas *c2 = new TCanvas("c2","RECO_vs_GEN_Eta",200,10,600,600);
	c2->cd();
	h1_GEN_eta->Draw();
	h1_ak5_pf_eta->Draw("same");
	c2->Print("RECO_vs_GEN_Eta.png");

	TCanvas *c13 = new TCanvas("c13","RECO_vs_GEN_Phi",200,10,600,600);
	c13->cd();
	h1_GEN_phi->Draw();
	h1_ak5_pf_phi->Draw("same");
	c13->Print("RECO_vs_GEN_Phi.png");

	TCanvas *c3 = new TCanvas("c3","RECO_vs_GEN_dR",200,10,600,600);
	c3->cd();
	h1_GEN_zjet_dr->Draw();
	h1_ak5_pf_zjet_dr->Draw("same");
	c3->Print("RECO_vs_GEN_dR.png");

	TCanvas *c4 = new TCanvas("c4","RECO_vs_GEN_dphi",200,10,600,600);
	c4->cd();
	h1_GEN_zjet_dphi->Draw();
	h1_ak5_pf_zjet_dphi->Draw("same");
	c4->Print("RECO_vs_GEN_dphi.png");

	TCanvas *c5 = new TCanvas("c5","RECO_vs_GEN_rhoSW",200,10,600,600);
	c5->cd();
	h1_GEN_rhoSW->Draw();
	h1_ak5_pf_rhoSW->Draw("same");
	c5->Print("RECO_vs_GEN_rhoSW.png");

	TCanvas *c6 = new TCanvas("c6","RECO_vs_GEN_rhoHand",200,10,600,600);
	c6->cd();
	h1_GEN_rhoHand->Draw();
	h1_ak5_pf_rhoHand->Draw("same");
	c6->Print("RECO_vs_GEN_rhoHand.png");

	TCanvas *c5 = new TCanvas("c7","RECO_vs_GEN_rhoHand2",200,10,600,600);
	c7->cd();
	h1_GEN_rhoHand2->Draw();
	h1_ak5_pf_rhoHand2->Draw("same");
	c7->Print("RECO_vs_GEN_rhoHand2.png");

	TCanvas *c8 = new TCanvas("c8","RECO_vs_GEN_rhoGrid",200,10,600,600);
	c8->cd();
	h1_GEN_rhoGrid->Draw();
	h1_ak5_pf_rhoGrid->Draw("same");
	c8->Print("RECO_vs_GEN_rhoGrid.png");

	TCanvas *c9 = new TCanvas("c9","RECO_vs_GEN_rhoSW_vs_nPV",200,10,600,600);
	c9->cd();
	h2_GEN_rhoSW_vs_nPV->Draw();
	h2_ak5_pf_rhoSW_vs_nPV->Draw("same");
	c9->Print("RECO_vs_GEN_rhoSW_vs_nPV.png");

	TCanvas *c10 = new TCanvas("c10","RECO_vs_GEN_rhoHand_vs_nPV",200,10,600,600);
	c10->cd();
	h2_GEN_rhoHand_vs_nPV->Draw();
	h2_ak5_pf_rhoHand_vs_nPV->Draw("same");
	c10->Print("RECO_vs_GEN_rhoHand_vs_nPV.png");

	TCanvas *c11 = new TCanvas("c11","RECO_vs_GEN_rhoHand2_vs_nPV",200,10,600,600);
	c11->cd();
	h2_GEN_rhoHand2_vs_nPV->Draw();
	h2_ak5_pf_rhoHand2_vs_nPV->Draw("same");
	c11->Print("RECO_vs_GEN_rhoHand2_vs_nPV.png");

	TCanvas *c12 = new TCanvas("c12","RECO_vs_GEN_rhoGrid_vs_nPV",200,10,600,600);
	c12->cd();
	h2_GEN_rhoGrid_vs_nPV->Draw();
	h2_ak5_pf_rhoGrid_vs_nPV->Draw("same");
	c12->Print("RECO_vs_GEN_rhoGrid_vs_nPV.png");

	TCanvas *c14 = new TCanvas("c14","f",200,10,600,600);
	c14->cd();
	h1_PF_match->Draw();
	h1_PFCHS_match->Draw("same");
	c14->Print("PF_vs_PFCHS_matching_efficiency.png");
	//ratio

	Draw_and_Save(h1_PFCor);
	Draw_and_Save(h1_PFCor_uncorr);
	Draw_and_Save(h1_PFCor_afterL1);
	Draw_and_Save(h1_PFCor_afterL2);
	Draw_and_Save(h1_ak5_pf);
	Draw_and_Save(h1_ak5_pf_uncorr);
	Draw_and_Save(h1_ak5_pf_l1_rhoSW);
	Draw_and_Save(h1_ak5_pf_l1_rhoHand);
	Draw_and_Save(h1_ak5_pf_l1_rhoHand2);
	Draw_and_Save(h1_ak5_pf_l1_rhoGrid);
	Draw_and_Save(h1_ak5_pfchs);
	Draw_and_Save(h1_ak5_pfchs_uncorr);
	Draw_and_Save(h1_ak5_pfchs_l1_rhoSW);
	Draw_and_Save(h1_ak5_pfchs_l1_rhoHand);
	Draw_and_Save(h1_ak5_pfchs_l1_rhoHand2);
	Draw_and_Save(h1_ak5_pfchs_l1_rhoGrid);

	Draw_and_Save(h1_area_PFCor);
	Draw_and_Save(h1_area_pf);
	Draw_and_Save(h1_area_pfchs);

	//h1_rho_PFCor); 
	/*h1_rhoSW_pf); 
	  h1_rhoHand_pf);
	  h1_rhoHand2_pf);
	  h1_rhoGrid_pf);
	  h1_rhoSW_pfchs);
	  h1_rhoHand_pfchs);
	  h1_rhoHand2_pfchs);
	  h1_rhoGrid_pfchs);*/


	Draw_and_Save(h1_nPV);

	Draw_and_Save(h1_GEN_pt);
	Draw_and_Save(h1_GEN_eta);
	Draw_and_Save(h1_GEN_zjet_dr);
	Draw_and_Save(h1_GEN_zjet_dphi);
	Draw_and_Save(h1_GEN_rhoSW);
	Draw_and_Save(h1_GEN_rhoHand);
	Draw_and_Save(h1_GEN_rhoHand2);
	Draw_and_Save(h1_GEN_rhoGrid);
	Draw_and_Save(h2_GEN_rhoSW_vs_nPV);
	Draw_and_Save(h2_GEN_rhoHand_vs_nPV);
	Draw_and_Save(h2_GEN_rhoHand2_vs_nPV);
	Draw_and_Save(h2_GEN_rhoGrid_vs_nPV);

	Draw_and_Save(h1_ak5_pf_pt_uncorr);
	Draw_and_Save(h1_ak5_pf_eta);
	Draw_and_Save(h1_ak5_pf_zjet_dr);
	Draw_and_Save(h1_ak5_pf_zjet_dphi);
	Draw_and_Save(h1_ak5_pf_rhoSW);
	Draw_and_Save(h1_ak5_pf_rhoHand);
	Draw_and_Save(h1_ak5_pf_rhoHand2);
	Draw_and_Save(h1_ak5_pf_rhoGrid);
	Draw_and_Save(h2_ak5_pf_rhoSW_vs_nPV);
	Draw_and_Save(h2_ak5_pf_rhoHand_vs_nPV);
	Draw_and_Save(h2_ak5_pf_rhoHand2_vs_nPV);
	Draw_and_Save(h2_ak5_pf_rhoGrid_vs_nPV);

	Draw_and_Save(h1_ak5_pfchs_pt_uncorr);
	Draw_and_Save(h1_ak5_pfchs_eta);
	Draw_and_Save(h1_ak5_pfchs_zjet_dr);
	Draw_and_Save(h1_ak5_pfchs_zjet_dphi);
	Draw_and_Save(h1_ak5_pfchs_rhoSW);
	Draw_and_Save(h1_ak5_pfchs_rhoHand);
	Draw_and_Save(h1_ak5_pfchs_rhoHand2);
	Draw_and_Save(h1_ak5_pfchs_rhoGrid);
	Draw_and_Save(h2_ak5_pfchs_rhoSW_vs_nPV);
	Draw_and_Save(h2_ak5_pfchs_rhoHand_vs_nPV);
	Draw_and_Save(h2_ak5_pfchs_rhoHand2_vs_nPV);
	Draw_and_Save(h2_ak5_pfchs_rhoGrid_vs_nPV);


	Draw_and_Save(h2_ak5_pfchs_ratioHand_vs_nPV);
	Draw_and_Save(h2_ak5_pfchs_ratioHand_vs_ptHand);
	Draw_and_Save(h2_ak5_pfchs_ratioHand_vs_eta);

	Draw_and_Save(h1_PF_match);
	Draw_and_Save(h1_PFCHS_match);


	Draw_and_Save(h1_GEN_mass                );
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

	h1_GEN_mass.SetLineColor(2);
	h1_GEN_mass.Draw();
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
	h1_GEN_mass.Draw();
	//h1_PFCor_mass.Draw("same");
	h1_AK5_PFCHS_mass_uncorr.Draw("same");
	//h1_AK5_PFCHS_mass_rhoArea.Draw("same");
	//h1_AK5_PFCHS_mass_rhoGArea.Draw("same");
	h1_AK5_PFCHS_mass_rho4Area.Draw("same");
	//h1_AK5_PFCHS_mass_rhoG4Area.Draw("same");
	h1_AK5_PFCHS_mass_rhom4Area.Draw("same");
	c16->Print("GEN_vs_PFCor_vs_PFCHS_jetmass.png");



	//c
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
