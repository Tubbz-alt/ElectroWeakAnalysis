#define MyClass_cxx
#include "MyClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

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
    
   Double_t tmp_PFCor_pt=0.;
   Double_t tmp_PFCor_Pt_uncorr=0.;
    
   Double_t tmp_AK5_PF_pt=0.;
   Double_t tmp_AK5_PF_eta=0.;
   Double_t tmp_AK5_PF_phi=0.;
   Double_t tmp_AK5_PF_pt_uncorr=0.;
   Double_t tmp_AK5_PF_pt_L1_rhoSW=0.;
   Double_t tmp_AK5_PF_pt_L1_rhoHand=0.;
   Double_t tmp_AK5_PF_pt_L1_rhoHand2=0.;
   Double_t tmp_AK5_PF_pt_L1_rhoGrid=0.;
    
   Double_t tmp_AK5_PFCHS_pt=0.;
    Double_t tmp_AK5_PFCHS_eta=0.;
    Double_t tmp_AK5_PFCHS_phi=0.;
   Double_t tmp_AK5_PFCHS_pt_uncorr=0.;
   Double_t tmp_AK5_PFCHS_pt_L1_rhoSW=0.;
   Double_t tmp_AK5_PFCHS_pt_L1_rhoHand=0.;
   Double_t tmp_AK5_PFCHS_pt_L1_rhoHand2=0.;
   Double_t tmp_AK5_PFCHS_pt_L1_rhoGrid=0.;
    
    Double_t tmp_Z_eta=0.;
    Double_t tmp_Z_phi=0.;
    
   Double_t tmp_event_nPV=0.;
    
   Double_t tmp_GEN_rhoSW=0,;
   Double_t tmp_GEN_rhoHand=0,;
   Double_t tmp_GEN_rhoHand2=0,;
   Double_t tmp_GEN_rhoGrid=0,;
   Double_t tmp_AK5_PF_rhoSW=0,;
   Double_t tmp_AK5_PF_rhoHand=0,;
   Double_t tmp_AK5_PF_rhoHand2=0,;
   Double_t tmp_AK5_PF_rhoGrid=0,;

   Double_t ratio=0.;
    Double_t dr=0,;
    Double_t dphi=0.;
   
    TH1D h1_GEN_pt("h1_GEN_pt","h1_GEN_pt",100,0,2.5);
    TH1D h1_ak5_pf_pt("h1_ak5_pf_pt","h1_ak5_pf_pt",100,0,2.5);

    TH1D h1_GEN_eta("h1_GEN_eta","h1_GEN_eta",100,0,2.5);
    TH1D h1_ak5_pf_eta("h1_ak5_pf_eta","h1_ak5_pf_eta",100,0,2.5);

    TH1D h1_GEN_phi("h1_GEN_phi","h1_GEN_phi",100,0,2.5);
    TH1D h1_ak5_pf_phi("h1_ak5_pf_phi","h1_ak5_pf_phi",100,0,2.5);
    
    TH1D h1_GEN_zjet_dr("h1_GEN_zjet_dr","h1_GEN_zjet_dr",100,0,2.5);
    TH1D h1_ak5_pf_zjet_dr("h1_ak5_pf_zjet_dr","h1_ak5_pf_zjet_dr",100,0,2.5);
    
    TH1D h1_GEN_zjet_dphi("h1_GEN_zjet_dphi","h1_GEN_zjet_dphi",100,0,2.5);
    TH1D h1_ak5_pf_zjet_dphi("h1_ak5_pf_zjet_dphi","h1_ak5_pf_zjet_dphi",100,0,2.5);
    
    
    TH1D h1_nPV("h1_nPV","h1_nPV",100,0,2.5);
    
    TH1D h1_GEN_rhoSW("h1_GEN_rhoSW","h1_GEN_rhoSW",100,0,2.5);
    TH1D h1_ak5_pf_rhoSW("h1_ak5_pf_rhoSW","h1_ak5_pf_rhoSW",100,0,2.5);
    
    TH1D h1_GEN_rhoHand("h1_GEN_rhoHand","h1_GEN_rhoHand",100,0,2.5);
    TH1D h1_ak5_pf_rhoHand("h1_ak5_pf_rhoHand","h1_ak5_pf_rhoHand",100,0,2.5);
    
    TH1D h1_GEN_rhoHand2("h1_GEN_rhoHand2","h1_GEN_rhoHand2",100,0,2.5);
    TH1D h1_ak5_pf_rhoHand2("h1_ak5_pf_rhoHand2","h1_ak5_pf_rhoHand2",100,0,2.5);
    
    TH1D h1_GEN_rhoGrid("h1_GEN_rhoGrid","h1_GEN_rhoGrid",100,0,2.5);
    TH1D h1_ak5_pf_rhoGrid("h1_ak5_pf_rhoGrid","h1_ak5_pf_rhoGrid",100,0,2.5);
    
    TH2D h1_GEN_rhoSW_vs_nPV("h1_GEN_rhoSW_vs_nPV","h1_GEN_rhoSW_vs_nPV",100,0,2.5,100,0,2.5);
    TH2D h1_ak5_pf_rhoSW_vs_nPV("h1_ak5_pf_rhoSW_vs_nPV","h1_ak5_pf_rhoSW_vs_nPV",100,0,2.5,100,0,2.5);
    
    TH2D h1_GEN_rhoHand_vs_nPV("h1_GEN_rhoHand_vs_nPV","h1_GEN_rhoHand_vs_nPV",100,0,2.5,100,0,2.5);
    TH2D h1_ak5_pf_rhoHand_vs_nPV("h1_ak5_pf_rhoHand_vs_nPV","h1_ak5_pf_rhoHand_vs_nPV",100,0,2.5,100,0,2.5);
    
    TH2D h1_GEN_rhoHand2_vs_nPV("h1_GEN_rhoHand2_vs_nPV","h1_GEN_rhoHand2_vs_nPV",100,0,2.5,100,0,2.5);
    TH2D h1_ak5_pf_rhoHand2_vs_nPV("h1_ak5_pf_rhoHand2_vs_nPV","h1_ak5_pf_rhoHand2_vs_nPV",100,0,2.5,100,0,2.5);
    
    TH2D h1_GEN_rhoGrid_vs_nPV("h1_GEN_rhoGrid_vs_nPV","h1_GEN_rhoGrid_vs_nPV",100,0,2.5,100,0,2.5);
    TH2D h1_ak5_pf_rhoGrid_vs_nPV("h1_ak5_pf_rhoGrid_vs_nPV","h1_ak5_pf_rhoGrid_vs_nPV",100,0,2.5,100,0,2.5);
    
    TH1D h1_PFCor_pt("h1_PFCor_pt","h1_PFCor_pt",100,0,2.5);
   TH1D h1_PFCor_pt_uncorr("h1_PFCor_pt_uncorr","h1_PFCor_pt_uncorr",100,0,2.5);
    
   TH1D h1_ak5_pf("h1_ak5_pf","h1_ak5_pf",100,0,2.5);
    
   TH1D h1_ak5_pf_uncorr("h1_ak5_pf_uncorr","h1_ak5_pf_uncorr",100,0,2.5);
   TH1D h1_ak5_pf_l1_rhoSW("h1_ak5_pf_l1_rhoSW","h1_ak5_pf_l1_rhoSW",100,0,2.5);
   TH1D h1_ak5_pf_l1_rhoHand("h1_ak5_pf_l1_rhoHand","h1_ak5_pf_l1_rhoHand",100,0,2.5);
   TH1D h1_ak5_pf_l1_rhoHand2("h1_ak5_pf_l1_rhoHand2","h1_ak5_pf_l1_rhoHand2",100,0,2.5);
   TH1D h1_ak5_pf_l1_rhoGrid("h1_ak5_pf_l1_rhoGrid","h1_ak5_pf_l1_rhoGrid",100,0,2.5);
   TH1D h1_ak5_pfchs("h1_ak5_pfchs","h1_ak5_pfchs",100,0,2.5);
   TH1D h1_ak5_pfchs_uncorr("h1_ak5_pfchs_uncorr","h1_ak5_pfchs_uncorr",100,0,2.5);
   TH1D h1_ak5_pfchs_l1_rhoSW("h1_ak5_pfchs_l1_rhoSW","h1_ak5_pfchs_l1_rhoSW",100,0,2.5);
   TH1D h1_ak5_pfchs_l1_rhoHand("h1_ak5_pfchs_l1_rhoHand","h1_ak5_pfchs_l1_rhoHand",100,0,2.5);
   TH1D h1_ak5_pfchs_l1_rhoHand2("h1_ak5_pfchs_l1_rhoHand2","h1_ak5_pfchs_l1_rhoHand2",100,0,2.5);
   TH1D h1_ak5_pfchs_l1_rhoGrid("h1_ak5_pfchs_l1_rhoGrid","h1_ak5_pfchs_l1_rhoGrid",100,0,2.5);


   // For GEN-RECO matching
   Double_t gen_jet_eta=gen_jet_phi=0.;
   Double_t PFCor_jet_eta=PFCor_jet_phi=0.;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

	  gen_jet_eta=GenGroomedJet_AK5_GEN_eta[0];
	  gen_jet_phi=GenGroomedJet_AK5_GEN_phi[0];
	  PFCor_jet_eta=JetPFCor_Eta[0];
	  PFCor_jet_phi=JetPFCor_Phi[0];

	  if( TMath::Sqrt( (gen_jet_eta-PFCor_jet_eta)*(gen_jet_eta-PFCor_jet_eta) + (gen_jet_phi-PFCor_jet_phi)*(gen_jet_phi-PFCor_jet_phi) ) <0.3 ){
		  tmp_GEN_pt = GenGroomedJet_AK5_GEN_pt[0];
          tmp_GEN_eta = GenGroomedJet_AK5_GEN_eta[0];
          tmp_GEN_phi = GenGroomedJet_AK5_GEN_phi[0];
          
          
		  tmp_PFCor_pt = JetPFCor_Pt[0];
		  tmp_PFCor_Pt_uncorr = JetPFCor_Pt_uncorr[0];
          
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
          tmp_AK5_PF_rhoSW = GroomedJet_AK5_PF_rhoSW[0];
          tmp_GEN_rhoHand = GenGroomedJet_AK5_GEN_rhohand[0];
          tmp_AK5_PF_rhoHand = GroomedJet_AK5_PF_rhohand[0];
          tmp_GEN_rhoHand2 = GenGroomedJet_AK5_GEN_rhohand1[0];
          tmp_AK5_PF_rhoHand2 = GroomedJet_AK5_PF_rhohand1[0];
          tmp_GEN_rhoGrid = GenGroomedJet_AK5_GEN_rhogrid[0];
          tmp_AK5_PF_rhoGrid = GroomedJet_AK5_PF_rhogrid[0];
          
          ratio = tmp_PFCor_pt/tmp_GEN_pt;
		  h1_PFCor_pt.Fill(ratio);
          
          ratio = tmp_PFCor_Pt_uncorr/tmp_GEN_pt;
		  h1_PFCor_pt_uncorr.Fill(ratio);
          
          ratio = tmp_AK5_PF_pt/tmp_GEN_pt;
		  h1_ak5_pf.Fill(ratio);
          
		  ratio = tmp_AK5_PF_pt_uncorr/tmp_GEN_pt;
		  h1_ak5_pf_uncorr.Fill(ratio);
          
          ratio = tmp_AK5_PF_pt_L1_rhoSW/tmp_GEN_pt;
		  h1_ak5_pf_l1_rhoSW.Fill(ratio);
          
          ratio = tmp_AK5_PF_pt_L1_rhoHand/tmp_GEN_pt;
		  h1_ak5_pf_l1_rhoHand.Fill(ratio);
          
          ratio = tmp_AK5_PF_pt_L1_rhoHand2/tmp_GEN_pt;
		  h1_ak5_pf_l1_rhoHand2.Fill(ratio);
          
          ratio = tmp_AK5_PF_pt_L1_rhoGrid/tmp_GEN_pt;
		  h1_ak5_pf_l1_rhoGrid.Fill(ratio);
          
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
          h1_ak5_pf_pt.Fill(tmp_AK5_PF_pt);
          //eta
          h1_GEN_eta.Fill(tmp_GEN_eta);
          h1_ak5_pf_eta.Fill(tmp_AK5_PF_eta);
          //dr
          dr = TMath::Sqrt( (Z_eta-tmp_GEN_eta)*(Z_eta-tmp_GEN_eta) + (Z_phi-tmp_GEN_phi)*(Z_phi-tmp_GEN_phi) );
          h1_GEN_zjet_dr.Fill(dr);
          dr = TMath::Sqrt( (Z_eta-tmp_AK5_PF_eta)*(Z_eta-tmp_AK5_PF_eta) + (Z_phi-tmp_AK5_PF_phi)*(Z_phi-tmp_AK5_PF_phi) );
          h1_ak5_pf_zjet_dr.Fill(dr);
          //dphi
          dphi = TMath::Sqrt( (Z_phi-tmp_GEN_phi)*(Z_phi-tmp_GEN_phi) );
          h1_GEN_zjet_dphi.Fill(dphi);
          dphi = TMath::Sqrt( (Z_phi-tmp_AK5_PF_phi)*(Z_phi-tmp_AK5_PF_phi) );
          h1_ak5_pf_zjet_dphi.Fill(dphi);
          
          h1_nPV.Fill(tmp_event_nPV);
          
          h1_GEN_rhoSW.Fill(tmp_GEN_rhoSW);
          h1_ak5_pf_rhoSW.Fill(tmp_AK5_PF_rhoSW);
          
          h1_GEN_rhoHand.Fill(tmp_GEN_rhoHand);
          h1_ak5_pf_rhoHand.Fill(tmp_AK5_PF_rhoHand);
          
          h1_GEN_rhoHand2.Fill(tmp_GEN_rhoHand2);
          h1_ak5_pf_rhoHand2.Fill(tmp_AK5_PF_rhoHand2);
          
          h1_GEN_rhoHand2.Fill(tmp_GEN_rhoGrid);
          h1_ak5_pf_rhoGrid.Fill(tmp_AK5_PF_rhoGrid);
          
          h1_GEN_rhoSW_vs_nPV.Fill(tmp_GEN_rhoSW,tmp_event_nPV);
          h1_ak5_pf_rhoSW_vs_nPV.Fill(tmp_AK5_PF_rhoSW,tmp_event_nPV);
          
          h1_GEN_rhoHand_vs_nPV.Fill(tmp_GEN_rhoHand,tmp_event_nPV);
          h1_ak5_pf_rhoHand_vs_nPV.Fill(tmp_AK5_PF_rhoHand,tmp_event_nPV);
          
          h1_GEN_rhoHand2_vs_nPV.Fill(tmp_GEN_rhoHand2,tmp_event_nPV);
          h1_ak5_pf_rhoHand2_vs_nPV.Fill(tmp_AK5_PF_rhoHand2,tmp_event_nPV);
          
          h1_GEN_rhoGrid_vs_nPV.Fill(tmp_GEN_rhoGrid,tmp_event_nPV);
          h1_ak5_pf_rhoGrid_vs_nPV.Fill(tmp_AK5_PF_rhoGrid,tmp_event_nPV);

          
          
		 
	  }
   }
    c1 = new TCanvas("c1","RECO vs. GEN Pt",200,10,700,500);
    c1.cd();
    h1_GEN_pt.Draw();
    h1_ak5_pf_pt.Draw("same");
    c1.Print("RECO vs. GEN Pt.pdf");
    
    c2 = new TCanvas("c2","RECO vs. GEN Eta",200,10,700,500);
    c2.cd();
    h1_GEN_eta.Draw();
    h1_ak5_pf_eta.Draw("same");
    c2.Print("RECO vs. GEN Eta.pdf");
    
    
    h1_PFCor_pt.Write();
    h1_PFCor_pt_uncorr.Write();
    h1_ak5_pf.Write();
    h1_ak5_pf_uncorr.Write();
    h1_ak5_pf_l1_rhoSW.Write();
    h1_ak5_pf_l1_rhoHand.Write();
    h1_ak5_pf_l1_rhoHand2.Write();
    h1_ak5_pf_l1_rhoGrid.Write();
    h1_ak5_pfchs.Write();
    h1_ak5_pfchs_uncorr.Write();
    h1_ak5_pfchs_l1_rhoSW.Write();
    h1_ak5_pfchs_l1_rhoHand.Write();
    h1_ak5_pfchs_l1_rhoHand2.Write();
    h1_ak5_pfchs_l1_rhoGrid.Write();
    
    h1_GEN_pt.Write();
    h1_ak5_pf_pt.Write();
    h1_GEN_eta.Write();
    h1_ak5_pf_eta.Write();
    h1_GEN_zjet_dr.Write();
    h1_ak5_pf_zjet_dr.Write();
    h1_GEN_zjet_dphi.Write();
    h1_ak5_pf_zjet_dphi.Write();
    h1_nPV.Write();
    h1_GEN_rhoSW.Write();
    h1_GEN_rhoHand.Write();
    h1_GEN_rhoHand2.Write();
    h1_GEN_rhoGrid.Write();
    h1_ak5_pf_rhoSW.Write();
    h1_ak5_pf_rhoHand.Write();
    h1_ak5_pf_rhoHand2.Write();
    h1_ak5_pf_rhoGrid.Write();
    h1_GEN_rhoSW_vs_nPV.Write();
    h1_GEN_rhoHand_vs_nPV.Write();
    h1_GEN_rhoHand2_vs_nPV.Write();
    h1_GEN_rhoGrid_vs_nPV.Write();
    h1_ak5_pf_rhoSW_vs_nPV.Write();
    h1_ak5_pf_rhoHand_vs_nPV.Write();
    h1_ak5_pf_rhoHand2_vs_nPV.Write();
    h1_ak5_pf_rhoGrid_vs_nPV.Write();
    
    
    
    
    
}
