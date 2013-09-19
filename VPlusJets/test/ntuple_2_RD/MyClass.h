//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Sep 13 01:06:43 2013 by ROOT version 5.34/05
// from TTree ZJet/V+jets Tree
// found on file: zmumujetsanalysisntuple.root
//////////////////////////////////////////////////////////

#ifndef MyClass_h
#define MyClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "iostream"
using namespace std;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class MyClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           numPFCorJets;
   Int_t           numPFCorJetBTags;
   Float_t         JetPFCor_Et[8];
   Float_t         JetPFCor_Pt[8];
   Float_t         JetPFCor_Pt_uncorr[8];
   Float_t         JetPFCor_Pt_afterL1[8];
   Float_t         JetPFCor_Pt_afterL2[8];
   Float_t         JetPFCor_Eta[8];
   Float_t         JetPFCor_Phi[8];
   Float_t         JetPFCor_Theta[8];
   Float_t         JetPFCor_Px[8];
   Float_t         JetPFCor_Py[8];
   Float_t         JetPFCor_Pz[8];
   Float_t         JetPFCor_E[8];
   Float_t         JetPFCor_Y[8];
   Float_t         JetPFCor_Mass[8];
   Float_t         JetPFCor_etaetaMoment[8];
   Float_t         JetPFCor_phiphiMoment[8];
   Float_t         JetPFCor_etaphiMoment[8];
   Float_t         JetPFCor_maxDistance[8];
   Int_t           JetPFCor_nConstituents[8];
   Float_t         JetPFCor_Area[8];
   Float_t         VplusPFCorJet_Mass[8];
   Float_t         JetPFCor_dphiBoson[8];
   Float_t         JetPFCor_detaBoson[8];
   Float_t         JetPFCor_dRBoson[8];
   Float_t         JetPFCor_dphiMET[8];
   Float_t         JetPFCor_bDiscriminator[8];
   Float_t         JetPFCor_bDiscriminatorSSVHE[8];
   Float_t         JetPFCor_bDiscriminatorTCHE[8];
   Float_t         JetPFCor_bDiscriminatorCSV[8];
   Float_t         JetPFCor_bDiscriminatorJP[8];
   Float_t         JetPFCor_bDiscriminatorSSVHP[8];
   Float_t         JetPFCor_bDiscriminatorTCHP[8];
   Float_t         JetPFCor_secVertexMass[8];
   Float_t         JetPFCor_ChargedHadronEnergy[8];
   Float_t         JetPFCor_ChargedHadronEnergyFrac[8];
   Float_t         JetPFCor_NeutralHadronEnergy[8];
   Float_t         JetPFCor_NeutralHadronEnergyFrac[8];
   Float_t         JetPFCor_ChargedEmEnergy[8];
   Float_t         JetPFCor_ChargedEmEnergyFrac[8];
   Float_t         JetPFCor_ChargedMuEnergy[8];
   Float_t         JetPFCor_ChargedMuEnergyFrac[8];
   Float_t         JetPFCor_NeutralEmEnergy[8];
   Float_t         JetPFCor_NeutralEmEnergyFrac[8];
   Int_t           JetPFCor_ChargedMultiplicity[8];
   Int_t           JetPFCor_NeutralMultiplicity[8];
   Int_t           JetPFCor_MuonMultiplicity[8];
   Float_t         JetPFCor_PhotonEnergy[8];
   Float_t         JetPFCor_PhotonEnergyFraction[8];
   Float_t         JetPFCor_ElectronEnergy[8];
   Float_t         JetPFCor_ElectronEnergyFraction[8];
   Float_t         JetPFCor_MuonEnergy[8];
   Float_t         JetPFCor_MuonEnergyFraction[8];
   Float_t         JetPFCor_HFHadronEnergy[8];
   Float_t         JetPFCor_HFHadronEnergyFraction[8];
   Float_t         JetPFCor_HFEMEnergy[8];
   Float_t         JetPFCor_HFEMEnergyFraction[8];
   Int_t           JetPFCor_ChargedHadronMultiplicity[8];
   Int_t           JetPFCor_NeutralHadronMultiplicity[8];
   Int_t           JetPFCor_PhotonMultiplicity[8];
   Int_t           JetPFCor_ElectronMultiplicity[8];
   Int_t           JetPFCor_HFHadronMultiplicity[8];
   Int_t           JetPFCor_HFEMMultiplicity[8];
   Float_t         JetPFCor_SumPtCands[8];
   Float_t         JetPFCor_SumPt2Cands[8];
   Float_t         JetPFCor_rmsCands[8];
   Float_t         JetPFCor_PtD[8];
   Float_t         JetPFCor_QGLikelihood[8];
   Float_t         GroomedJet_AK5_PF_pt_uncorr[6];
   Int_t           GroomedJet_AK5_PF_number_jet_central;
   Float_t         GroomedJet_AK5_PF_mass_uncorr[6];
   Float_t         GroomedJet_AK5_PF_mass_tr_uncorr[6];
   Float_t         GroomedJet_AK5_PF_mass_ft_uncorr[6];
   Float_t         GroomedJet_AK5_PF_mass_pr_uncorr[6];
   Float_t         GroomedJet_AK5_PF_tau2tau1[6];
   Float_t         GroomedJet_AK5_PF_tau1[6];
   Float_t         GroomedJet_AK5_PF_tau2[6];
   Float_t         GroomedJet_AK5_PF_tau3[6];
   Float_t         GroomedJet_AK5_PF_tau4[6];
   Float_t         GroomedJet_AK5_PF_massdrop_pr_uncorr[6];
   Float_t         GroomedJet_AK5_PF_pt[6];
   Float_t         GroomedJet_AK5_PF_eta[6];
   Float_t         GroomedJet_AK5_PF_phi[6];
   Float_t         GroomedJet_AK5_PF_e[6];
   Float_t         GroomedJet_AK5_PF_pt_L1_rhoSW[6];
   Float_t         GroomedJet_AK5_PF_pt_L1_rhoHand[6];
   Float_t         GroomedJet_AK5_PF_pt_L1_rhoHand2[6];
   Float_t         GroomedJet_AK5_PF_pt_L1_rhoGrid[6];
   Float_t         GroomedJet_AK5_PF_pt_tr_uncorr[6];
   Float_t         GroomedJet_AK5_PF_pt_tr[6];
   Float_t         GroomedJet_AK5_PF_eta_tr[6];
   Float_t         GroomedJet_AK5_PF_phi_tr[6];
   Float_t         GroomedJet_AK5_PF_e_tr[6];
   Float_t         GroomedJet_AK5_PF_pt_ft_uncorr[6];
   Float_t         GroomedJet_AK5_PF_pt_ft[6];
   Float_t         GroomedJet_AK5_PF_eta_ft[6];
   Float_t         GroomedJet_AK5_PF_phi_ft[6];
   Float_t         GroomedJet_AK5_PF_e_ft[6];
   Float_t         GroomedJet_AK5_PF_pt_pr_uncorr[6];
   Float_t         GroomedJet_AK5_PF_pt_pr[6];
   Float_t         GroomedJet_AK5_PF_eta_pr[6];
   Float_t         GroomedJet_AK5_PF_phi_pr[6];
   Float_t         GroomedJet_AK5_PF_e_pr[6];
   Float_t         GroomedJet_AK5_PF_prsubjet1_px[6];
   Float_t         GroomedJet_AK5_PF_prsubjet1_py[6];
   Float_t         GroomedJet_AK5_PF_prsubjet1_pz[6];
   Float_t         GroomedJet_AK5_PF_prsubjet1_e[6];
   Float_t         GroomedJet_AK5_PF_prsubjet2_px[6];
   Float_t         GroomedJet_AK5_PF_prsubjet2_py[6];
   Float_t         GroomedJet_AK5_PF_prsubjet2_pz[6];
   Float_t         GroomedJet_AK5_PF_prsubjet2_e[6];
   Float_t         GroomedJet_AK5_PF_mass[6];
   Float_t         GroomedJet_AK5_PF_mass_tr[6];
   Float_t         GroomedJet_AK5_PF_mass_ft[6];
   Float_t         GroomedJet_AK5_PF_mass_pr[6];
   Float_t         GroomedJet_AK5_PF_massdrop_pr[6];
   Float_t         GroomedJet_AK5_PF_area[6];
   Float_t         GroomedJet_AK5_PF_area_tr[6];
   Float_t         GroomedJet_AK5_PF_area_ft[6];
   Float_t         GroomedJet_AK5_PF_area_pr[6];
   Float_t         GroomedJet_AK5_PF_jetconstituents[6];
   Float_t         GroomedJet_AK5_PF_jetcharge[6];
   Float_t         GroomedJet_AK5_PF_rcores[11][6];
   Float_t         GroomedJet_AK5_PF_ptcores[11][6];
   Float_t         GroomedJet_AK5_PF_planarflow[11][6];
   Float_t         GroomedJet_AK5_PF_qjetmass[50];
   Float_t         GroomedJet_AK5_PF_qjetmassdrop[50];
   Float_t         GroomedJet_AK5_PF_constituents0_eta[100];
   Float_t         GroomedJet_AK5_PF_constituents0_phi[100];
   Float_t         GroomedJet_AK5_PF_constituents0_e[100];
   Int_t           GroomedJet_AK5_PF_nconstituents0;
   Float_t         GroomedJet_AK5_PF_constituents0pr_eta[100];
   Float_t         GroomedJet_AK5_PF_constituents0pr_phi[100];
   Float_t         GroomedJet_AK5_PF_constituents0pr_e[100];
   Int_t           GroomedJet_AK5_PF_nconstituents0pr;
   Double_t        GroomedJet_AK5_PF_rhoSW;
   Double_t        GroomedJet_AK5_PF_rhohand;
   Double_t        GroomedJet_AK5_PF_rhohand2;
   Double_t        GroomedJet_AK5_PF_rhogrid;
   Float_t         GroomedJet_AK5_PFCHS_pt_uncorr[6];
   Int_t           GroomedJet_AK5_PFCHS_number_jet_central;
   Float_t         GroomedJet_AK5_PFCHS_mass_uncorr[6];
   Float_t         GroomedJet_AK5_PFCHS_mass_tr_uncorr[6];
   Float_t         GroomedJet_AK5_PFCHS_mass_ft_uncorr[6];
   Float_t         GroomedJet_AK5_PFCHS_mass_pr_uncorr[6];
   Float_t         GroomedJet_AK5_PFCHS_tau2tau1[6];
   Float_t         GroomedJet_AK5_PFCHS_tau1[6];
   Float_t         GroomedJet_AK5_PFCHS_tau2[6];
   Float_t         GroomedJet_AK5_PFCHS_tau3[6];
   Float_t         GroomedJet_AK5_PFCHS_tau4[6];
   Float_t         GroomedJet_AK5_PFCHS_massdrop_pr_uncorr[6];
   Float_t         GroomedJet_AK5_PFCHS_pt[6];
   Float_t         GroomedJet_AK5_PFCHS_eta[6];
   Float_t         GroomedJet_AK5_PFCHS_phi[6];
   Float_t         GroomedJet_AK5_PFCHS_e[6];
   Float_t         GroomedJet_AK5_PFCHS_pt_L1_rhoSW[6];
   Float_t         GroomedJet_AK5_PFCHS_pt_L1_rhoHand[6];
   Float_t         GroomedJet_AK5_PFCHS_pt_L1_rhoHand2[6];
   Float_t         GroomedJet_AK5_PFCHS_pt_L1_rhoGrid[6];
   Float_t         GroomedJet_AK5_PFCHS_pt_tr_uncorr[6];
   Float_t         GroomedJet_AK5_PFCHS_pt_tr[6];
   Float_t         GroomedJet_AK5_PFCHS_eta_tr[6];
   Float_t         GroomedJet_AK5_PFCHS_phi_tr[6];
   Float_t         GroomedJet_AK5_PFCHS_e_tr[6];
   Float_t         GroomedJet_AK5_PFCHS_pt_ft_uncorr[6];
   Float_t         GroomedJet_AK5_PFCHS_pt_ft[6];
   Float_t         GroomedJet_AK5_PFCHS_eta_ft[6];
   Float_t         GroomedJet_AK5_PFCHS_phi_ft[6];
   Float_t         GroomedJet_AK5_PFCHS_e_ft[6];
   Float_t         GroomedJet_AK5_PFCHS_pt_pr_uncorr[6];
   Float_t         GroomedJet_AK5_PFCHS_pt_pr[6];
   Float_t         GroomedJet_AK5_PFCHS_eta_pr[6];
   Float_t         GroomedJet_AK5_PFCHS_phi_pr[6];
   Float_t         GroomedJet_AK5_PFCHS_e_pr[6];
   Float_t         GroomedJet_AK5_PFCHS_prsubjet1_px[6];
   Float_t         GroomedJet_AK5_PFCHS_prsubjet1_py[6];
   Float_t         GroomedJet_AK5_PFCHS_prsubjet1_pz[6];
   Float_t         GroomedJet_AK5_PFCHS_prsubjet1_e[6];
   Float_t         GroomedJet_AK5_PFCHS_prsubjet2_px[6];
   Float_t         GroomedJet_AK5_PFCHS_prsubjet2_py[6];
   Float_t         GroomedJet_AK5_PFCHS_prsubjet2_pz[6];
   Float_t         GroomedJet_AK5_PFCHS_prsubjet2_e[6];
   Float_t         GroomedJet_AK5_PFCHS_mass[6];
   Float_t         GroomedJet_AK5_PFCHS_mass_tr[6];
   Float_t         GroomedJet_AK5_PFCHS_mass_ft[6];
   Float_t         GroomedJet_AK5_PFCHS_mass_pr[6];
   Float_t         GroomedJet_AK5_PFCHS_massdrop_pr[6];
   Float_t         GroomedJet_AK5_PFCHS_area[6];
   Float_t         GroomedJet_AK5_PFCHS_area_tr[6];
   Float_t         GroomedJet_AK5_PFCHS_area_ft[6];
   Float_t         GroomedJet_AK5_PFCHS_area_pr[6];
   Float_t         GroomedJet_AK5_PFCHS_jetconstituents[6];
   Float_t         GroomedJet_AK5_PFCHS_jetcharge[6];
   Float_t         GroomedJet_AK5_PFCHS_rcores[11][6];
   Float_t         GroomedJet_AK5_PFCHS_ptcores[11][6];
   Float_t         GroomedJet_AK5_PFCHS_planarflow[11][6];
   Float_t         GroomedJet_AK5_PFCHS_qjetmass[50];
   Float_t         GroomedJet_AK5_PFCHS_qjetmassdrop[50];
   Float_t         GroomedJet_AK5_PFCHS_constituents0_eta[100];
   Float_t         GroomedJet_AK5_PFCHS_constituents0_phi[100];
   Float_t         GroomedJet_AK5_PFCHS_constituents0_e[100];
   Int_t           GroomedJet_AK5_PFCHS_nconstituents0;
   Float_t         GroomedJet_AK5_PFCHS_constituents0pr_eta[100];
   Float_t         GroomedJet_AK5_PFCHS_constituents0pr_phi[100];
   Float_t         GroomedJet_AK5_PFCHS_constituents0pr_e[100];
   Int_t           GroomedJet_AK5_PFCHS_nconstituents0pr;
   Double_t        GroomedJet_AK5_PFCHS_rhoSW;
   Double_t        GroomedJet_AK5_PFCHS_rhohand;
   Double_t        GroomedJet_AK5_PFCHS_rhohand2;
   Double_t        GroomedJet_AK5_PFCHS_rhogrid;
   Float_t         GenGroomedJet_AK5_GEN_pt_uncorr[6];
   Int_t           GenGroomedJet_AK5_GEN_number_jet_central;
   Float_t         GenGroomedJet_AK5_GEN_mass_uncorr[6];
   Float_t         GenGroomedJet_AK5_GEN_mass_tr_uncorr[6];
   Float_t         GenGroomedJet_AK5_GEN_mass_ft_uncorr[6];
   Float_t         GenGroomedJet_AK5_GEN_mass_pr_uncorr[6];
   Float_t         GenGroomedJet_AK5_GEN_tau2tau1[6];
   Float_t         GenGroomedJet_AK5_GEN_tau1[6];
   Float_t         GenGroomedJet_AK5_GEN_tau2[6];
   Float_t         GenGroomedJet_AK5_GEN_tau3[6];
   Float_t         GenGroomedJet_AK5_GEN_tau4[6];
   Float_t         GenGroomedJet_AK5_GEN_massdrop_pr_uncorr[6];
   Float_t         GenGroomedJet_AK5_GEN_pt[6];
   Float_t         GenGroomedJet_AK5_GEN_eta[6];
   Float_t         GenGroomedJet_AK5_GEN_phi[6];
   Float_t         GenGroomedJet_AK5_GEN_e[6];
   Float_t         GenGroomedJet_AK5_GEN_pt_L1_rhoSW[6];
   Float_t         GenGroomedJet_AK5_GEN_pt_L1_rhoHand[6];
   Float_t         GenGroomedJet_AK5_GEN_pt_L1_rhoHand2[6];
   Float_t         GenGroomedJet_AK5_GEN_pt_L1_rhoGrid[6];
   Float_t         GenGroomedJet_AK5_GEN_pt_tr_uncorr[6];
   Float_t         GenGroomedJet_AK5_GEN_pt_tr[6];
   Float_t         GenGroomedJet_AK5_GEN_eta_tr[6];
   Float_t         GenGroomedJet_AK5_GEN_phi_tr[6];
   Float_t         GenGroomedJet_AK5_GEN_e_tr[6];
   Float_t         GenGroomedJet_AK5_GEN_pt_ft_uncorr[6];
   Float_t         GenGroomedJet_AK5_GEN_pt_ft[6];
   Float_t         GenGroomedJet_AK5_GEN_eta_ft[6];
   Float_t         GenGroomedJet_AK5_GEN_phi_ft[6];
   Float_t         GenGroomedJet_AK5_GEN_e_ft[6];
   Float_t         GenGroomedJet_AK5_GEN_pt_pr_uncorr[6];
   Float_t         GenGroomedJet_AK5_GEN_pt_pr[6];
   Float_t         GenGroomedJet_AK5_GEN_eta_pr[6];
   Float_t         GenGroomedJet_AK5_GEN_phi_pr[6];
   Float_t         GenGroomedJet_AK5_GEN_e_pr[6];
   Float_t         GenGroomedJet_AK5_GEN_prsubjet1_px[6];
   Float_t         GenGroomedJet_AK5_GEN_prsubjet1_py[6];
   Float_t         GenGroomedJet_AK5_GEN_prsubjet1_pz[6];
   Float_t         GenGroomedJet_AK5_GEN_prsubjet1_e[6];
   Float_t         GenGroomedJet_AK5_GEN_prsubjet2_px[6];
   Float_t         GenGroomedJet_AK5_GEN_prsubjet2_py[6];
   Float_t         GenGroomedJet_AK5_GEN_prsubjet2_pz[6];
   Float_t         GenGroomedJet_AK5_GEN_prsubjet2_e[6];
   Float_t         GenGroomedJet_AK5_GEN_mass[6];
   Float_t         GenGroomedJet_AK5_GEN_mass_tr[6];
   Float_t         GenGroomedJet_AK5_GEN_mass_ft[6];
   Float_t         GenGroomedJet_AK5_GEN_mass_pr[6];
   Float_t         GenGroomedJet_AK5_GEN_massdrop_pr[6];
   Float_t         GenGroomedJet_AK5_GEN_area[6];
   Float_t         GenGroomedJet_AK5_GEN_area_tr[6];
   Float_t         GenGroomedJet_AK5_GEN_area_ft[6];
   Float_t         GenGroomedJet_AK5_GEN_area_pr[6];
   Float_t         GenGroomedJet_AK5_GEN_jetconstituents[6];
   Float_t         GenGroomedJet_AK5_GEN_jetcharge[6];
   Float_t         GenGroomedJet_AK5_GEN_rcores[11][6];
   Float_t         GenGroomedJet_AK5_GEN_ptcores[11][6];
   Float_t         GenGroomedJet_AK5_GEN_planarflow[11][6];
   Float_t         GenGroomedJet_AK5_GEN_qjetmass[50];
   Float_t         GenGroomedJet_AK5_GEN_qjetmassdrop[50];
   Float_t         GenGroomedJet_AK5_GEN_constituents0_eta[100];
   Float_t         GenGroomedJet_AK5_GEN_constituents0_phi[100];
   Float_t         GenGroomedJet_AK5_GEN_constituents0_e[100];
   Int_t           GenGroomedJet_AK5_GEN_nconstituents0;
   Float_t         GenGroomedJet_AK5_GEN_constituents0pr_eta[100];
   Float_t         GenGroomedJet_AK5_GEN_constituents0pr_phi[100];
   Float_t         GenGroomedJet_AK5_GEN_constituents0pr_e[100];
   Int_t           GenGroomedJet_AK5_GEN_nconstituents0pr;
   Double_t        GenGroomedJet_AK5_GEN_rhoSW;
   Double_t        GenGroomedJet_AK5_GEN_rhohand;
   Double_t        GenGroomedJet_AK5_GEN_rhohand2;
   Double_t        GenGroomedJet_AK5_GEN_rhogrid;
   Float_t         Z_mass;
   Float_t         Z_mt;
   Float_t         Z_mtMVA;
   Float_t         Z_px;
   Float_t         Z_py;
   Float_t         Z_pz;
   Float_t         Z_e;
   Float_t         Z_pt;
   Float_t         Z_et;
   Float_t         Z_eta;
   Float_t         Z_phi;
   Float_t         Z_vx;
   Float_t         Z_vy;
   Float_t         Z_vz;
   Float_t         Z_y;
   Float_t         Z_muplus_px;
   Float_t         Z_muplus_py;
   Float_t         Z_muplus_pz;
   Float_t         Z_muplus_e;
   Float_t         Z_muplus_pt;
   Float_t         Z_muplus_et;
   Float_t         Z_muplus_eta;
   Float_t         Z_muplus_theta;
   Float_t         Z_muplus_phi;
   Int_t           Z_muplus_charge;
   Float_t         Z_muplus_vx;
   Float_t         Z_muplus_vy;
   Float_t         Z_muplus_vz;
   Float_t         Z_muplus_y;
   Float_t         Z_muplus_trackiso;
   Float_t         Z_muplus_hcaliso;
   Float_t         Z_muplus_ecaliso;
   Int_t           Z_muplus_type;
   Int_t           Z_muplus_numberOfChambers;
   Int_t           Z_muplus_numberOfMatches;
   Float_t         Z_muplus_d0bsp;
   Float_t         Z_muplus_dz000;
   Float_t         Z_muplus_dzPV;
   Float_t         Z_muplus_pfiso_sumChargedHadronPt;
   Float_t         Z_muplus_pfiso_sumChargedParticlePt;
   Float_t         Z_muplus_pfiso_sumNeutralHadronEt;
   Float_t         Z_muplus_pfiso_sumPhotonEt;
   Float_t         Z_muplus_pfiso_sumPUPt;
   Float_t         Z_muminus_px;
   Float_t         Z_muminus_py;
   Float_t         Z_muminus_pz;
   Float_t         Z_muminus_e;
   Float_t         Z_muminus_pt;
   Float_t         Z_muminus_et;
   Float_t         Z_muminus_eta;
   Float_t         Z_muminus_theta;
   Float_t         Z_muminus_phi;
   Int_t           Z_muminus_charge;
   Float_t         Z_muminus_vx;
   Float_t         Z_muminus_vy;
   Float_t         Z_muminus_vz;
   Float_t         Z_muminus_y;
   Float_t         Z_muminus_trackiso;
   Float_t         Z_muminus_hcaliso;
   Float_t         Z_muminus_ecaliso;
   Int_t           Z_muminus_type;
   Int_t           Z_muminus_numberOfChambers;
   Int_t           Z_muminus_numberOfMatches;
   Float_t         Z_muminus_d0bsp;
   Float_t         Z_muminus_dz000;
   Float_t         Z_muminus_dzPV;
   Float_t         Z_muminus_pfiso_sumChargedHadronPt;
   Float_t         Z_muminus_pfiso_sumChargedParticlePt;
   Float_t         Z_muminus_pfiso_sumNeutralHadronEt;
   Float_t         Z_muminus_pfiso_sumPhotonEt;
   Float_t         Z_muminus_pfiso_sumPUPt;
   Float_t         Z_Photon_pt_gen;
   Float_t         Z_muplus_px_gen;
   Float_t         Z_muplus_py_gen;
   Float_t         Z_muplus_pz_gen;
   Float_t         Z_muplus_e_gen;
   Float_t         Z_muplus_pt_gen;
   Float_t         Z_muplus_et_gen;
   Float_t         Z_muplus_eta_gen;
   Float_t         Z_muplus_theta_gen;
   Float_t         Z_muplus_phi_gen;
   Int_t           Z_muplus_charge_gen;
   Float_t         Z_muplus_vx_gen;
   Float_t         Z_muplus_vy_gen;
   Float_t         Z_muplus_vz_gen;
   Float_t         Z_muplus_y_gen;
   Float_t         Z_muminus_px_gen;
   Float_t         Z_muminus_py_gen;
   Float_t         Z_muminus_pz_gen;
   Float_t         Z_muminus_e_gen;
   Float_t         Z_muminus_pt_gen;
   Float_t         Z_muminus_et_gen;
   Float_t         Z_muminus_eta_gen;
   Float_t         Z_muminus_theta_gen;
   Float_t         Z_muminus_phi_gen;
   Int_t           Z_muminus_charge_gen;
   Float_t         Z_muminus_vx_gen;
   Float_t         Z_muminus_vy_gen;
   Float_t         Z_muminus_vz_gen;
   Float_t         Z_muminus_y_gen;
   Int_t           event_runNo;
   Int_t           event_evtNo;
   Int_t           event_lumi;
   Int_t           event_bunch;
   Int_t           event_nPV;
   Float_t         event_met_pfmet;
   Float_t         event_met_pfsumet;
   Float_t         event_met_pfmetsignificance;
   Float_t         event_met_pfmetPhi;
   Float_t         event_metMVA_met;
   Float_t         event_metMVA_sumet;
   Float_t         event_metMVA_metsignificance;
   Float_t         event_metMVA_metPhi;
   Float_t         event_fastJetRho;
   Float_t         event_met_genmet;
   Float_t         event_met_gensumet;
   Float_t         event_met_genmetsignificance;
   Float_t         event_met_genmetPhi;
   Float_t         event_mcPU_totnvtx;
   Float_t         event_mcPU_trueInteractions;
   Float_t         event_mcPU_bx[3];
   Float_t         event_mcPU_nvtx[3];

   // List of branches
   TBranch        *b_numPFCorJets;   //!
   TBranch        *b_numPFCorJetBTags;   //!
   TBranch        *b_JetPFCor_Et;   //!
   TBranch        *b_JetPFCor_Pt;   //!
   TBranch        *b_JetPFCor_Pt_uncorr;   //!
   TBranch        *b_JetPFCor_Pt_afterL1;   //!
   TBranch        *b_JetPFCor_Pt_afterL2;   //!
   TBranch        *b_JetPFCor_Eta;   //!
   TBranch        *b_JetPFCor_Phi;   //!
   TBranch        *b_JetPFCor_Theta;   //!
   TBranch        *b_JetPFCor_Px;   //!
   TBranch        *b_JetPFCor_Py;   //!
   TBranch        *b_JetPFCor_Pz;   //!
   TBranch        *b_JetPFCor_E;   //!
   TBranch        *b_JetPFCor_Y;   //!
   TBranch        *b_JetPFCor_Mass;   //!
   TBranch        *b_JetPFCor_etaetaMoment;   //!
   TBranch        *b_JetPFCor_phiphiMoment;   //!
   TBranch        *b_JetPFCor_etaphiMoment;   //!
   TBranch        *b_JetPFCor_maxDistance;   //!
   TBranch        *b_JetPFCor_nConstituents;   //!
   TBranch        *b_JetPFCor_Area;   //!
   TBranch        *b_VplusPFCorJet_Mass;   //!
   TBranch        *b_JetPFCor_dphiBoson;   //!
   TBranch        *b_JetPFCor_detaBoson;   //!
   TBranch        *b_JetPFCor_dRBoson;   //!
   TBranch        *b_JetPFCor_dphiMET;   //!
   TBranch        *b_JetPFCor_bDiscriminator;   //!
   TBranch        *b_JetPFCor_bDiscriminatorSSVHE;   //!
   TBranch        *b_JetPFCor_bDiscriminatorTCHE;   //!
   TBranch        *b_JetPFCor_bDiscriminatorCSV;   //!
   TBranch        *b_JetPFCor_bDiscriminatorJP;   //!
   TBranch        *b_JetPFCor_bDiscriminatorSSVHP;   //!
   TBranch        *b_JetPFCor_bDiscriminatorTCHP;   //!
   TBranch        *b_JetPFCor_secVertexMass;   //!
   TBranch        *b_JetPFCor_ChargedHadronEnergy;   //!
   TBranch        *b_JetPFCor_ChargedHadronEnergyFrac;   //!
   TBranch        *b_JetPFCor_NeutralHadronEnergy;   //!
   TBranch        *b_JetPFCor_NeutralHadronEnergyFrac;   //!
   TBranch        *b_JetPFCor_ChargedEmEnergy;   //!
   TBranch        *b_JetPFCor_ChargedEmEnergyFrac;   //!
   TBranch        *b_JetPFCor_ChargedMuEnergy;   //!
   TBranch        *b_JetPFCor_ChargedMuEnergyFrac;   //!
   TBranch        *b_JetPFCor_NeutralEmEnergy;   //!
   TBranch        *b_JetPFCor_NeutralEmEnergyFrac;   //!
   TBranch        *b_JetPFCor_ChargedMultiplicity;   //!
   TBranch        *b_JetPFCor_NeutralMultiplicity;   //!
   TBranch        *b_JetPFCor_MuonMultiplicity;   //!
   TBranch        *b_JetPFCor_PhotonEnergy;   //!
   TBranch        *b_JetPFCor_PhotonEnergyFraction;   //!
   TBranch        *b_JetPFCor_ElectronEnergy;   //!
   TBranch        *b_JetPFCor_ElectronEnergyFraction;   //!
   TBranch        *b_JetPFCor_MuonEnergy;   //!
   TBranch        *b_JetPFCor_MuonEnergyFraction;   //!
   TBranch        *b_JetPFCor_HFHadronEnergy;   //!
   TBranch        *b_JetPFCor_HFHadronEnergyFraction;   //!
   TBranch        *b_JetPFCor_HFEMEnergy;   //!
   TBranch        *b_JetPFCor_HFEMEnergyFraction;   //!
   TBranch        *b_JetPFCor_ChargedHadronMultiplicity;   //!
   TBranch        *b_JetPFCor_NeutralHadronMultiplicity;   //!
   TBranch        *b_JetPFCor_PhotonMultiplicity;   //!
   TBranch        *b_JetPFCor_ElectronMultiplicity;   //!
   TBranch        *b_JetPFCor_HFHadronMultiplicity;   //!
   TBranch        *b_JetPFCor_HFEMMultiplicity;   //!
   TBranch        *b_JetPFCor_SumPtCands;   //!
   TBranch        *b_JetPFCor_SumPt2Cands;   //!
   TBranch        *b_JetPFCor_rmsCands;   //!
   TBranch        *b_JetPFCor_PtD;   //!
   TBranch        *b_JetPFCor_QGLikelihood;   //!
   TBranch        *b_GroomedJet_AK5_PF_pt_uncorr;   //!
   TBranch        *b_GroomedJet_AK5_PF_number_jet_central;   //!
   TBranch        *b_GroomedJet_AK5_PF_mass_uncorr;   //!
   TBranch        *b_GroomedJet_AK5_PF_mass_tr_uncorr;   //!
   TBranch        *b_GroomedJet_AK5_PF_mass_ft_uncorr;   //!
   TBranch        *b_GroomedJet_AK5_PF_mass_pr_uncorr;   //!
   TBranch        *b_GroomedJet_AK5_PF_tau2tau1;   //!
   TBranch        *b_GroomedJet_AK5_PF_tau1;   //!
   TBranch        *b_GroomedJet_AK5_PF_tau2;   //!
   TBranch        *b_GroomedJet_AK5_PF_tau3;   //!
   TBranch        *b_GroomedJet_AK5_PF_tau4;   //!
   TBranch        *b_GroomedJet_AK5_PF_massdrop_pr_uncorr;   //!
   TBranch        *b_GroomedJet_AK5_PF_pt;   //!
   TBranch        *b_GroomedJet_AK5_PF_eta;   //!
   TBranch        *b_GroomedJet_AK5_PF_phi;   //!
   TBranch        *b_GroomedJet_AK5_PF_e;   //!
   TBranch        *b_GroomedJet_AK5_PF_pt_L1_rhoSW;   //!
   TBranch        *b_GroomedJet_AK5_PF_pt_L1_rhoHand;   //!
   TBranch        *b_GroomedJet_AK5_PF_pt_L1_rhoHand2;   //!
   TBranch        *b_GroomedJet_AK5_PF_pt_L1_rhoGrid;   //!
   TBranch        *b_GroomedJet_AK5_PF_pt_tr_uncorr;   //!
   TBranch        *b_GroomedJet_AK5_PF_pt_tr;   //!
   TBranch        *b_GroomedJet_AK5_PF_eta_tr;   //!
   TBranch        *b_GroomedJet_AK5_PF_phi_tr;   //!
   TBranch        *b_GroomedJet_AK5_PF_e_tr;   //!
   TBranch        *b_GroomedJet_AK5_PF_pt_ft_uncorr;   //!
   TBranch        *b_GroomedJet_AK5_PF_pt_ft;   //!
   TBranch        *b_GroomedJet_AK5_PF_eta_ft;   //!
   TBranch        *b_GroomedJet_AK5_PF_phi_ft;   //!
   TBranch        *b_GroomedJet_AK5_PF_e_ft;   //!
   TBranch        *b_GroomedJet_AK5_PF_pt_pr_uncorr;   //!
   TBranch        *b_GroomedJet_AK5_PF_pt_pr;   //!
   TBranch        *b_GroomedJet_AK5_PF_eta_pr;   //!
   TBranch        *b_GroomedJet_AK5_PF_phi_pr;   //!
   TBranch        *b_GroomedJet_AK5_PF_e_pr;   //!
   TBranch        *b_GroomedJet_AK5_PF_prsubjet1_px;   //!
   TBranch        *b_GroomedJet_AK5_PF_prsubjet1_py;   //!
   TBranch        *b_GroomedJet_AK5_PF_prsubjet1_pz;   //!
   TBranch        *b_GroomedJet_AK5_PF_prsubjet1_e;   //!
   TBranch        *b_GroomedJet_AK5_PF_prsubjet2_px;   //!
   TBranch        *b_GroomedJet_AK5_PF_prsubjet2_py;   //!
   TBranch        *b_GroomedJet_AK5_PF_prsubjet2_pz;   //!
   TBranch        *b_GroomedJet_AK5_PF_prsubjet2_e;   //!
   TBranch        *b_GroomedJet_AK5_PF_mass;   //!
   TBranch        *b_GroomedJet_AK5_PF_mass_tr;   //!
   TBranch        *b_GroomedJet_AK5_PF_mass_ft;   //!
   TBranch        *b_GroomedJet_AK5_PF_mass_pr;   //!
   TBranch        *b_GroomedJet_AK5_PF_massdrop_pr;   //!
   TBranch        *b_GroomedJet_AK5_PF_area;   //!
   TBranch        *b_GroomedJet_AK5_PF_area_tr;   //!
   TBranch        *b_GroomedJet_AK5_PF_area_ft;   //!
   TBranch        *b_GroomedJet_AK5_PF_area_pr;   //!
   TBranch        *b_GroomedJet_AK5_PF_jetconstituents;   //!
   TBranch        *b_GroomedJet_AK5_PF_jetcharge;   //!
   TBranch        *b_GroomedJet_AK5_PF_rcores;   //!
   TBranch        *b_GroomedJet_AK5_PF_ptcores;   //!
   TBranch        *b_GroomedJet_AK5_PF_planarflow;   //!
   TBranch        *b_GroomedJet_AK5_PF_qjetmass;   //!
   TBranch        *b_GroomedJet_AK5_PF_qjetmassdrop;   //!
   TBranch        *b_GroomedJet_AK5_PF_constituents0_eta;   //!
   TBranch        *b_GroomedJet_AK5_PF_constituents0_phi;   //!
   TBranch        *b_GroomedJet_AK5_PF_constituents0_e;   //!
   TBranch        *b_GroomedJet_AK5_PF_nconstituents0;   //!
   TBranch        *b_GroomedJet_AK5_PF_constituents0pr_eta;   //!
   TBranch        *b_GroomedJet_AK5_PF_constituents0pr_phi;   //!
   TBranch        *b_GroomedJet_AK5_PF_constituents0pr_e;   //!
   TBranch        *b_GroomedJet_AK5_PF_nconstituents0pr;   //!
   TBranch        *b_GroomedJet_AK5_PF_rhoSW;   //!
   TBranch        *b_GroomedJet_AK5_PF_rhohand;   //!
   TBranch        *b_GroomedJet_AK5_PF_rhohand2;   //!
   TBranch        *b_GroomedJet_AK5_PF_rhogrid;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_pt_uncorr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_number_jet_central;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_mass_uncorr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_mass_tr_uncorr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_mass_ft_uncorr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_mass_pr_uncorr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_tau2tau1;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_tau1;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_tau2;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_tau3;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_tau4;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_massdrop_pr_uncorr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_pt;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_eta;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_phi;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_e;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_pt_L1_rhoSW;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_pt_L1_rhoHand;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_pt_L1_rhoHand2;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_pt_L1_rhoGrid;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_pt_tr_uncorr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_pt_tr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_eta_tr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_phi_tr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_e_tr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_pt_ft_uncorr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_pt_ft;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_eta_ft;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_phi_ft;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_e_ft;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_pt_pr_uncorr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_pt_pr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_eta_pr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_phi_pr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_e_pr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_prsubjet1_px;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_prsubjet1_py;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_prsubjet1_pz;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_prsubjet1_e;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_prsubjet2_px;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_prsubjet2_py;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_prsubjet2_pz;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_prsubjet2_e;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_mass;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_mass_tr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_mass_ft;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_mass_pr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_massdrop_pr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_area;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_area_tr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_area_ft;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_area_pr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_jetconstituents;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_jetcharge;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_rcores;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_ptcores;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_planarflow;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_qjetmass;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_qjetmassdrop;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_constituents0_eta;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_constituents0_phi;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_constituents0_e;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_nconstituents0;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_constituents0pr_eta;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_constituents0pr_phi;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_constituents0pr_e;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_nconstituents0pr;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_rhoSW;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_rhohand;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_rhohand2;   //!
   TBranch        *b_GroomedJet_AK5_PFCHS_rhogrid;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_pt_uncorr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_number_jet_central;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_mass_uncorr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_mass_tr_uncorr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_mass_ft_uncorr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_mass_pr_uncorr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_tau2tau1;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_tau1;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_tau2;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_tau3;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_tau4;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_massdrop_pr_uncorr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_pt;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_eta;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_phi;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_e;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_pt_L1_rhoSW;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_pt_L1_rhoHand;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_pt_L1_rhoHand2;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_pt_L1_rhoGrid;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_pt_tr_uncorr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_pt_tr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_eta_tr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_phi_tr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_e_tr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_pt_ft_uncorr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_pt_ft;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_eta_ft;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_phi_ft;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_e_ft;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_pt_pr_uncorr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_pt_pr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_eta_pr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_phi_pr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_e_pr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_prsubjet1_px;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_prsubjet1_py;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_prsubjet1_pz;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_prsubjet1_e;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_prsubjet2_px;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_prsubjet2_py;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_prsubjet2_pz;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_prsubjet2_e;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_mass;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_mass_tr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_mass_ft;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_mass_pr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_massdrop_pr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_area;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_area_tr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_area_ft;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_area_pr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_jetconstituents;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_jetcharge;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_rcores;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_ptcores;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_planarflow;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_qjetmass;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_qjetmassdrop;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_constituents0_eta;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_constituents0_phi;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_constituents0_e;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_nconstituents0;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_constituents0pr_eta;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_constituents0pr_phi;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_constituents0pr_e;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_nconstituents0pr;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_rhoSW;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_rhohand;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_rhohand2;   //!
   TBranch        *b_GenGroomedJet_AK5_GEN_rhogrid;   //!
   TBranch        *b_Z_mass;   //!
   TBranch        *b_Z_mt;   //!
   TBranch        *b_Z_mtMVA;   //!
   TBranch        *b_Z_px;   //!
   TBranch        *b_Z_py;   //!
   TBranch        *b_Z_pz;   //!
   TBranch        *b_Z_e;   //!
   TBranch        *b_Z_pt;   //!
   TBranch        *b_Z_et;   //!
   TBranch        *b_Z_eta;   //!
   TBranch        *b_Z_phi;   //!
   TBranch        *b_Z_vx;   //!
   TBranch        *b_Z_vy;   //!
   TBranch        *b_Z_vz;   //!
   TBranch        *b_Z_y;   //!
   TBranch        *b_Z_muplus_px;   //!
   TBranch        *b_Z_muplus_py;   //!
   TBranch        *b_Z_muplus_pz;   //!
   TBranch        *b_Z_muplus_e;   //!
   TBranch        *b_Z_muplus_pt;   //!
   TBranch        *b_Z_muplus_et;   //!
   TBranch        *b_Z_muplus_eta;   //!
   TBranch        *b_Z_muplus_theta;   //!
   TBranch        *b_Z_muplus_phi;   //!
   TBranch        *b_Z_muplus_charge;   //!
   TBranch        *b_Z_muplus_vx;   //!
   TBranch        *b_Z_muplus_vy;   //!
   TBranch        *b_Z_muplus_vz;   //!
   TBranch        *b_Z_muplus_y;   //!
   TBranch        *b_Z_muplus_trackiso;   //!
   TBranch        *b_Z_muplus_hcaliso;   //!
   TBranch        *b_Z_muplus_ecaliso;   //!
   TBranch        *b_Z_muplus_type;   //!
   TBranch        *b_Z_muplus_numberOfChambers;   //!
   TBranch        *b_Z_muplus_numberOfMatches;   //!
   TBranch        *b_Z_muplus_d0bsp;   //!
   TBranch        *b_Z_muplus_dz000;   //!
   TBranch        *b_Z_muplus_dzPV;   //!
   TBranch        *b_Z_muplus_pfiso_sumChargedHadronPt;   //!
   TBranch        *b_Z_muplus_pfiso_sumChargedParticlePt;   //!
   TBranch        *b_Z_muplus_pfiso_sumNeutralHadronEt;   //!
   TBranch        *b_Z_muplus_pfiso_sumPhotonEt;   //!
   TBranch        *b_Z_muplus_pfiso_sumPUPt;   //!
   TBranch        *b_Z_muminus_px;   //!
   TBranch        *b_Z_muminus_py;   //!
   TBranch        *b_Z_muminus_pz;   //!
   TBranch        *b_Z_muminus_e;   //!
   TBranch        *b_Z_muminus_pt;   //!
   TBranch        *b_Z_muminus_et;   //!
   TBranch        *b_Z_muminus_eta;   //!
   TBranch        *b_Z_muminus_theta;   //!
   TBranch        *b_Z_muminus_phi;   //!
   TBranch        *b_Z_muminus_charge;   //!
   TBranch        *b_Z_muminus_vx;   //!
   TBranch        *b_Z_muminus_vy;   //!
   TBranch        *b_Z_muminus_vz;   //!
   TBranch        *b_Z_muminus_y;   //!
   TBranch        *b_Z_muminus_trackiso;   //!
   TBranch        *b_Z_muminus_hcaliso;   //!
   TBranch        *b_Z_muminus_ecaliso;   //!
   TBranch        *b_Z_muminus_type;   //!
   TBranch        *b_Z_muminus_numberOfChambers;   //!
   TBranch        *b_Z_muminus_numberOfMatches;   //!
   TBranch        *b_Z_muminus_d0bsp;   //!
   TBranch        *b_Z_muminus_dz000;   //!
   TBranch        *b_Z_muminus_dzPV;   //!
   TBranch        *b_Z_muminus_pfiso_sumChargedHadronPt;   //!
   TBranch        *b_Z_muminus_pfiso_sumChargedParticlePt;   //!
   TBranch        *b_Z_muminus_pfiso_sumNeutralHadronEt;   //!
   TBranch        *b_Z_muminus_pfiso_sumPhotonEt;   //!
   TBranch        *b_Z_muminus_pfiso_sumPUPt;   //!
   TBranch        *b_Z_Photon_pt_gen;   //!
   TBranch        *b_Z_muplus_px_gen;   //!
   TBranch        *b_Z_muplus_py_gen;   //!
   TBranch        *b_Z_muplus_pz_gen;   //!
   TBranch        *b_Z_muplus_e_gen;   //!
   TBranch        *b_Z_muplus_pt_gen;   //!
   TBranch        *b_Z_muplus_et_gen;   //!
   TBranch        *b_Z_muplus_eta_gen;   //!
   TBranch        *b_Z_muplus_theta_gen;   //!
   TBranch        *b_Z_muplus_phi_gen;   //!
   TBranch        *b_Z_muplus_charge_gen;   //!
   TBranch        *b_Z_muplus_vx_gen;   //!
   TBranch        *b_Z_muplus_vy_gen;   //!
   TBranch        *b_Z_muplus_vz_gen;   //!
   TBranch        *b_Z_muplus_y_gen;   //!
   TBranch        *b_Z_muminus_px_gen;   //!
   TBranch        *b_Z_muminus_py_gen;   //!
   TBranch        *b_Z_muminus_pz_gen;   //!
   TBranch        *b_Z_muminus_e_gen;   //!
   TBranch        *b_Z_muminus_pt_gen;   //!
   TBranch        *b_Z_muminus_et_gen;   //!
   TBranch        *b_Z_muminus_eta_gen;   //!
   TBranch        *b_Z_muminus_theta_gen;   //!
   TBranch        *b_Z_muminus_phi_gen;   //!
   TBranch        *b_Z_muminus_charge_gen;   //!
   TBranch        *b_Z_muminus_vx_gen;   //!
   TBranch        *b_Z_muminus_vy_gen;   //!
   TBranch        *b_Z_muminus_vz_gen;   //!
   TBranch        *b_Z_muminus_y_gen;   //!
   TBranch        *b_event_runNo;   //!
   TBranch        *b_event_evtNo;   //!
   TBranch        *b_event_lumi;   //!
   TBranch        *b_event_bunch;   //!
   TBranch        *b_event_nPV;   //!
   TBranch        *b_event_met_pfmet;   //!
   TBranch        *b_event_met_pfsumet;   //!
   TBranch        *b_event_met_pfmetsignificance;   //!
   TBranch        *b_event_met_pfmetPhi;   //!
   TBranch        *b_event_metMVA_met;   //!
   TBranch        *b_event_metMVA_sumet;   //!
   TBranch        *b_event_metMVA_metsignificance;   //!
   TBranch        *b_event_metMVA_metPhi;   //!
   TBranch        *b_event_fastJetRho;   //!
   TBranch        *b_event_met_genmet;   //!
   TBranch        *b_event_met_gensumet;   //!
   TBranch        *b_event_met_genmetsignificance;   //!
   TBranch        *b_event_met_genmetPhi;   //!
   TBranch        *b_event_mcPU_totnvtx;   //!
   TBranch        *b_event_mcPU_trueInteractions;   //!
   TBranch        *b_event_mcPU_bx;   //!
   TBranch        *b_event_mcPU_nvtx;   //!

   MyClass(TTree *tree=0);
   virtual ~MyClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MyClass_cxx
MyClass::MyClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("zmumujetsanalysisntuple.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("zmumujetsanalysisntuple.root");
      }
      f->GetObject("ZJet",tree);

   }
   Init(tree);
}

MyClass::~MyClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MyClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("numPFCorJets", &numPFCorJets, &b_numPFCorJets);
   fChain->SetBranchAddress("numPFCorJetBTags", &numPFCorJetBTags, &b_numPFCorJetBTags);
   fChain->SetBranchAddress("JetPFCor_Et", JetPFCor_Et, &b_JetPFCor_Et);
   fChain->SetBranchAddress("JetPFCor_Pt", JetPFCor_Pt, &b_JetPFCor_Pt);
   fChain->SetBranchAddress("JetPFCor_Pt_uncorr", JetPFCor_Pt_uncorr, &b_JetPFCor_Pt_uncorr);
   fChain->SetBranchAddress("JetPFCor_Pt_afterL1", JetPFCor_Pt_afterL1, &b_JetPFCor_Pt_afterL1);
   fChain->SetBranchAddress("JetPFCor_Pt_afterL2", JetPFCor_Pt_afterL2, &b_JetPFCor_Pt_afterL2);
   fChain->SetBranchAddress("JetPFCor_Eta", JetPFCor_Eta, &b_JetPFCor_Eta);
   fChain->SetBranchAddress("JetPFCor_Phi", JetPFCor_Phi, &b_JetPFCor_Phi);
   fChain->SetBranchAddress("JetPFCor_Theta", JetPFCor_Theta, &b_JetPFCor_Theta);
   fChain->SetBranchAddress("JetPFCor_Px", JetPFCor_Px, &b_JetPFCor_Px);
   fChain->SetBranchAddress("JetPFCor_Py", JetPFCor_Py, &b_JetPFCor_Py);
   fChain->SetBranchAddress("JetPFCor_Pz", JetPFCor_Pz, &b_JetPFCor_Pz);
   fChain->SetBranchAddress("JetPFCor_E", JetPFCor_E, &b_JetPFCor_E);
   fChain->SetBranchAddress("JetPFCor_Y", JetPFCor_Y, &b_JetPFCor_Y);
   fChain->SetBranchAddress("JetPFCor_Mass", JetPFCor_Mass, &b_JetPFCor_Mass);
   fChain->SetBranchAddress("JetPFCor_etaetaMoment", JetPFCor_etaetaMoment, &b_JetPFCor_etaetaMoment);
   fChain->SetBranchAddress("JetPFCor_phiphiMoment", JetPFCor_phiphiMoment, &b_JetPFCor_phiphiMoment);
   fChain->SetBranchAddress("JetPFCor_etaphiMoment", JetPFCor_etaphiMoment, &b_JetPFCor_etaphiMoment);
   fChain->SetBranchAddress("JetPFCor_maxDistance", JetPFCor_maxDistance, &b_JetPFCor_maxDistance);
   fChain->SetBranchAddress("JetPFCor_nConstituents", JetPFCor_nConstituents, &b_JetPFCor_nConstituents);
   fChain->SetBranchAddress("JetPFCor_Area", JetPFCor_Area, &b_JetPFCor_Area);
   fChain->SetBranchAddress("VplusPFCorJet_Mass", VplusPFCorJet_Mass, &b_VplusPFCorJet_Mass);
   fChain->SetBranchAddress("JetPFCor_dphiBoson", JetPFCor_dphiBoson, &b_JetPFCor_dphiBoson);
   fChain->SetBranchAddress("JetPFCor_detaBoson", JetPFCor_detaBoson, &b_JetPFCor_detaBoson);
   fChain->SetBranchAddress("JetPFCor_dRBoson", JetPFCor_dRBoson, &b_JetPFCor_dRBoson);
   fChain->SetBranchAddress("JetPFCor_dphiMET", JetPFCor_dphiMET, &b_JetPFCor_dphiMET);
   fChain->SetBranchAddress("JetPFCor_bDiscriminator", JetPFCor_bDiscriminator, &b_JetPFCor_bDiscriminator);
   fChain->SetBranchAddress("JetPFCor_bDiscriminatorSSVHE", JetPFCor_bDiscriminatorSSVHE, &b_JetPFCor_bDiscriminatorSSVHE);
   fChain->SetBranchAddress("JetPFCor_bDiscriminatorTCHE", JetPFCor_bDiscriminatorTCHE, &b_JetPFCor_bDiscriminatorTCHE);
   fChain->SetBranchAddress("JetPFCor_bDiscriminatorCSV", JetPFCor_bDiscriminatorCSV, &b_JetPFCor_bDiscriminatorCSV);
   fChain->SetBranchAddress("JetPFCor_bDiscriminatorJP", JetPFCor_bDiscriminatorJP, &b_JetPFCor_bDiscriminatorJP);
   fChain->SetBranchAddress("JetPFCor_bDiscriminatorSSVHP", JetPFCor_bDiscriminatorSSVHP, &b_JetPFCor_bDiscriminatorSSVHP);
   fChain->SetBranchAddress("JetPFCor_bDiscriminatorTCHP", JetPFCor_bDiscriminatorTCHP, &b_JetPFCor_bDiscriminatorTCHP);
   fChain->SetBranchAddress("JetPFCor_secVertexMass", JetPFCor_secVertexMass, &b_JetPFCor_secVertexMass);
   fChain->SetBranchAddress("JetPFCor_ChargedHadronEnergy", JetPFCor_ChargedHadronEnergy, &b_JetPFCor_ChargedHadronEnergy);
   fChain->SetBranchAddress("JetPFCor_ChargedHadronEnergyFrac", JetPFCor_ChargedHadronEnergyFrac, &b_JetPFCor_ChargedHadronEnergyFrac);
   fChain->SetBranchAddress("JetPFCor_NeutralHadronEnergy", JetPFCor_NeutralHadronEnergy, &b_JetPFCor_NeutralHadronEnergy);
   fChain->SetBranchAddress("JetPFCor_NeutralHadronEnergyFrac", JetPFCor_NeutralHadronEnergyFrac, &b_JetPFCor_NeutralHadronEnergyFrac);
   fChain->SetBranchAddress("JetPFCor_ChargedEmEnergy", JetPFCor_ChargedEmEnergy, &b_JetPFCor_ChargedEmEnergy);
   fChain->SetBranchAddress("JetPFCor_ChargedEmEnergyFrac", JetPFCor_ChargedEmEnergyFrac, &b_JetPFCor_ChargedEmEnergyFrac);
   fChain->SetBranchAddress("JetPFCor_ChargedMuEnergy", JetPFCor_ChargedMuEnergy, &b_JetPFCor_ChargedMuEnergy);
   fChain->SetBranchAddress("JetPFCor_ChargedMuEnergyFrac", JetPFCor_ChargedMuEnergyFrac, &b_JetPFCor_ChargedMuEnergyFrac);
   fChain->SetBranchAddress("JetPFCor_NeutralEmEnergy", JetPFCor_NeutralEmEnergy, &b_JetPFCor_NeutralEmEnergy);
   fChain->SetBranchAddress("JetPFCor_NeutralEmEnergyFrac", JetPFCor_NeutralEmEnergyFrac, &b_JetPFCor_NeutralEmEnergyFrac);
   fChain->SetBranchAddress("JetPFCor_ChargedMultiplicity", JetPFCor_ChargedMultiplicity, &b_JetPFCor_ChargedMultiplicity);
   fChain->SetBranchAddress("JetPFCor_NeutralMultiplicity", JetPFCor_NeutralMultiplicity, &b_JetPFCor_NeutralMultiplicity);
   fChain->SetBranchAddress("JetPFCor_MuonMultiplicity", JetPFCor_MuonMultiplicity, &b_JetPFCor_MuonMultiplicity);
   fChain->SetBranchAddress("JetPFCor_PhotonEnergy", JetPFCor_PhotonEnergy, &b_JetPFCor_PhotonEnergy);
   fChain->SetBranchAddress("JetPFCor_PhotonEnergyFraction", JetPFCor_PhotonEnergyFraction, &b_JetPFCor_PhotonEnergyFraction);
   fChain->SetBranchAddress("JetPFCor_ElectronEnergy", JetPFCor_ElectronEnergy, &b_JetPFCor_ElectronEnergy);
   fChain->SetBranchAddress("JetPFCor_ElectronEnergyFraction", JetPFCor_ElectronEnergyFraction, &b_JetPFCor_ElectronEnergyFraction);
   fChain->SetBranchAddress("JetPFCor_MuonEnergy", JetPFCor_MuonEnergy, &b_JetPFCor_MuonEnergy);
   fChain->SetBranchAddress("JetPFCor_MuonEnergyFraction", JetPFCor_MuonEnergyFraction, &b_JetPFCor_MuonEnergyFraction);
   fChain->SetBranchAddress("JetPFCor_HFHadronEnergy", JetPFCor_HFHadronEnergy, &b_JetPFCor_HFHadronEnergy);
   fChain->SetBranchAddress("JetPFCor_HFHadronEnergyFraction", JetPFCor_HFHadronEnergyFraction, &b_JetPFCor_HFHadronEnergyFraction);
   fChain->SetBranchAddress("JetPFCor_HFEMEnergy", JetPFCor_HFEMEnergy, &b_JetPFCor_HFEMEnergy);
   fChain->SetBranchAddress("JetPFCor_HFEMEnergyFraction", JetPFCor_HFEMEnergyFraction, &b_JetPFCor_HFEMEnergyFraction);
   fChain->SetBranchAddress("JetPFCor_ChargedHadronMultiplicity", JetPFCor_ChargedHadronMultiplicity, &b_JetPFCor_ChargedHadronMultiplicity);
   fChain->SetBranchAddress("JetPFCor_NeutralHadronMultiplicity", JetPFCor_NeutralHadronMultiplicity, &b_JetPFCor_NeutralHadronMultiplicity);
   fChain->SetBranchAddress("JetPFCor_PhotonMultiplicity", JetPFCor_PhotonMultiplicity, &b_JetPFCor_PhotonMultiplicity);
   fChain->SetBranchAddress("JetPFCor_ElectronMultiplicity", JetPFCor_ElectronMultiplicity, &b_JetPFCor_ElectronMultiplicity);
   fChain->SetBranchAddress("JetPFCor_HFHadronMultiplicity", JetPFCor_HFHadronMultiplicity, &b_JetPFCor_HFHadronMultiplicity);
   fChain->SetBranchAddress("JetPFCor_HFEMMultiplicity", JetPFCor_HFEMMultiplicity, &b_JetPFCor_HFEMMultiplicity);
   fChain->SetBranchAddress("JetPFCor_SumPtCands", JetPFCor_SumPtCands, &b_JetPFCor_SumPtCands);
   fChain->SetBranchAddress("JetPFCor_SumPt2Cands", JetPFCor_SumPt2Cands, &b_JetPFCor_SumPt2Cands);
   fChain->SetBranchAddress("JetPFCor_rmsCands", JetPFCor_rmsCands, &b_JetPFCor_rmsCands);
   fChain->SetBranchAddress("JetPFCor_PtD", JetPFCor_PtD, &b_JetPFCor_PtD);
   fChain->SetBranchAddress("JetPFCor_QGLikelihood", JetPFCor_QGLikelihood, &b_JetPFCor_QGLikelihood);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_pt_uncorr", GroomedJet_AK5_PF_pt_uncorr, &b_GroomedJet_AK5_PF_pt_uncorr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_number_jet_central", &GroomedJet_AK5_PF_number_jet_central, &b_GroomedJet_AK5_PF_number_jet_central);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_mass_uncorr", GroomedJet_AK5_PF_mass_uncorr, &b_GroomedJet_AK5_PF_mass_uncorr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_mass_tr_uncorr", GroomedJet_AK5_PF_mass_tr_uncorr, &b_GroomedJet_AK5_PF_mass_tr_uncorr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_mass_ft_uncorr", GroomedJet_AK5_PF_mass_ft_uncorr, &b_GroomedJet_AK5_PF_mass_ft_uncorr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_mass_pr_uncorr", GroomedJet_AK5_PF_mass_pr_uncorr, &b_GroomedJet_AK5_PF_mass_pr_uncorr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_tau2tau1", GroomedJet_AK5_PF_tau2tau1, &b_GroomedJet_AK5_PF_tau2tau1);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_tau1", GroomedJet_AK5_PF_tau1, &b_GroomedJet_AK5_PF_tau1);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_tau2", GroomedJet_AK5_PF_tau2, &b_GroomedJet_AK5_PF_tau2);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_tau3", GroomedJet_AK5_PF_tau3, &b_GroomedJet_AK5_PF_tau3);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_tau4", GroomedJet_AK5_PF_tau4, &b_GroomedJet_AK5_PF_tau4);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_massdrop_pr_uncorr", GroomedJet_AK5_PF_massdrop_pr_uncorr, &b_GroomedJet_AK5_PF_massdrop_pr_uncorr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_pt", GroomedJet_AK5_PF_pt, &b_GroomedJet_AK5_PF_pt);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_eta", GroomedJet_AK5_PF_eta, &b_GroomedJet_AK5_PF_eta);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_phi", GroomedJet_AK5_PF_phi, &b_GroomedJet_AK5_PF_phi);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_e", GroomedJet_AK5_PF_e, &b_GroomedJet_AK5_PF_e);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_pt_L1_rhoSW", GroomedJet_AK5_PF_pt_L1_rhoSW, &b_GroomedJet_AK5_PF_pt_L1_rhoSW);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_pt_L1_rhoHand", GroomedJet_AK5_PF_pt_L1_rhoHand, &b_GroomedJet_AK5_PF_pt_L1_rhoHand);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_pt_L1_rhoHand2", GroomedJet_AK5_PF_pt_L1_rhoHand2, &b_GroomedJet_AK5_PF_pt_L1_rhoHand2);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_pt_L1_rhoGrid", GroomedJet_AK5_PF_pt_L1_rhoGrid, &b_GroomedJet_AK5_PF_pt_L1_rhoGrid);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_pt_tr_uncorr", GroomedJet_AK5_PF_pt_tr_uncorr, &b_GroomedJet_AK5_PF_pt_tr_uncorr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_pt_tr", GroomedJet_AK5_PF_pt_tr, &b_GroomedJet_AK5_PF_pt_tr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_eta_tr", GroomedJet_AK5_PF_eta_tr, &b_GroomedJet_AK5_PF_eta_tr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_phi_tr", GroomedJet_AK5_PF_phi_tr, &b_GroomedJet_AK5_PF_phi_tr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_e_tr", GroomedJet_AK5_PF_e_tr, &b_GroomedJet_AK5_PF_e_tr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_pt_ft_uncorr", GroomedJet_AK5_PF_pt_ft_uncorr, &b_GroomedJet_AK5_PF_pt_ft_uncorr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_pt_ft", GroomedJet_AK5_PF_pt_ft, &b_GroomedJet_AK5_PF_pt_ft);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_eta_ft", GroomedJet_AK5_PF_eta_ft, &b_GroomedJet_AK5_PF_eta_ft);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_phi_ft", GroomedJet_AK5_PF_phi_ft, &b_GroomedJet_AK5_PF_phi_ft);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_e_ft", GroomedJet_AK5_PF_e_ft, &b_GroomedJet_AK5_PF_e_ft);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_pt_pr_uncorr", GroomedJet_AK5_PF_pt_pr_uncorr, &b_GroomedJet_AK5_PF_pt_pr_uncorr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_pt_pr", GroomedJet_AK5_PF_pt_pr, &b_GroomedJet_AK5_PF_pt_pr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_eta_pr", GroomedJet_AK5_PF_eta_pr, &b_GroomedJet_AK5_PF_eta_pr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_phi_pr", GroomedJet_AK5_PF_phi_pr, &b_GroomedJet_AK5_PF_phi_pr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_e_pr", GroomedJet_AK5_PF_e_pr, &b_GroomedJet_AK5_PF_e_pr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_prsubjet1_px", GroomedJet_AK5_PF_prsubjet1_px, &b_GroomedJet_AK5_PF_prsubjet1_px);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_prsubjet1_py", GroomedJet_AK5_PF_prsubjet1_py, &b_GroomedJet_AK5_PF_prsubjet1_py);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_prsubjet1_pz", GroomedJet_AK5_PF_prsubjet1_pz, &b_GroomedJet_AK5_PF_prsubjet1_pz);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_prsubjet1_e", GroomedJet_AK5_PF_prsubjet1_e, &b_GroomedJet_AK5_PF_prsubjet1_e);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_prsubjet2_px", GroomedJet_AK5_PF_prsubjet2_px, &b_GroomedJet_AK5_PF_prsubjet2_px);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_prsubjet2_py", GroomedJet_AK5_PF_prsubjet2_py, &b_GroomedJet_AK5_PF_prsubjet2_py);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_prsubjet2_pz", GroomedJet_AK5_PF_prsubjet2_pz, &b_GroomedJet_AK5_PF_prsubjet2_pz);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_prsubjet2_e", GroomedJet_AK5_PF_prsubjet2_e, &b_GroomedJet_AK5_PF_prsubjet2_e);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_mass", GroomedJet_AK5_PF_mass, &b_GroomedJet_AK5_PF_mass);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_mass_tr", GroomedJet_AK5_PF_mass_tr, &b_GroomedJet_AK5_PF_mass_tr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_mass_ft", GroomedJet_AK5_PF_mass_ft, &b_GroomedJet_AK5_PF_mass_ft);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_mass_pr", GroomedJet_AK5_PF_mass_pr, &b_GroomedJet_AK5_PF_mass_pr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_massdrop_pr", GroomedJet_AK5_PF_massdrop_pr, &b_GroomedJet_AK5_PF_massdrop_pr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_area", GroomedJet_AK5_PF_area, &b_GroomedJet_AK5_PF_area);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_area_tr", GroomedJet_AK5_PF_area_tr, &b_GroomedJet_AK5_PF_area_tr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_area_ft", GroomedJet_AK5_PF_area_ft, &b_GroomedJet_AK5_PF_area_ft);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_area_pr", GroomedJet_AK5_PF_area_pr, &b_GroomedJet_AK5_PF_area_pr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_jetconstituents", GroomedJet_AK5_PF_jetconstituents, &b_GroomedJet_AK5_PF_jetconstituents);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_jetcharge", GroomedJet_AK5_PF_jetcharge, &b_GroomedJet_AK5_PF_jetcharge);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_rcores", GroomedJet_AK5_PF_rcores, &b_GroomedJet_AK5_PF_rcores);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_ptcores", GroomedJet_AK5_PF_ptcores, &b_GroomedJet_AK5_PF_ptcores);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_planarflow", GroomedJet_AK5_PF_planarflow, &b_GroomedJet_AK5_PF_planarflow);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_qjetmass", GroomedJet_AK5_PF_qjetmass, &b_GroomedJet_AK5_PF_qjetmass);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_qjetmassdrop", GroomedJet_AK5_PF_qjetmassdrop, &b_GroomedJet_AK5_PF_qjetmassdrop);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_constituents0_eta", GroomedJet_AK5_PF_constituents0_eta, &b_GroomedJet_AK5_PF_constituents0_eta);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_constituents0_phi", GroomedJet_AK5_PF_constituents0_phi, &b_GroomedJet_AK5_PF_constituents0_phi);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_constituents0_e", GroomedJet_AK5_PF_constituents0_e, &b_GroomedJet_AK5_PF_constituents0_e);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_nconstituents0", &GroomedJet_AK5_PF_nconstituents0, &b_GroomedJet_AK5_PF_nconstituents0);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_constituents0pr_eta", GroomedJet_AK5_PF_constituents0pr_eta, &b_GroomedJet_AK5_PF_constituents0pr_eta);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_constituents0pr_phi", GroomedJet_AK5_PF_constituents0pr_phi, &b_GroomedJet_AK5_PF_constituents0pr_phi);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_constituents0pr_e", GroomedJet_AK5_PF_constituents0pr_e, &b_GroomedJet_AK5_PF_constituents0pr_e);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_nconstituents0pr", &GroomedJet_AK5_PF_nconstituents0pr, &b_GroomedJet_AK5_PF_nconstituents0pr);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_rhoSW", &GroomedJet_AK5_PF_rhoSW, &b_GroomedJet_AK5_PF_rhoSW);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_rhohand", &GroomedJet_AK5_PF_rhohand, &b_GroomedJet_AK5_PF_rhohand);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_rhohand2", &GroomedJet_AK5_PF_rhohand2, &b_GroomedJet_AK5_PF_rhohand2);
   fChain->SetBranchAddress("GroomedJet_AK5_PF_rhogrid", &GroomedJet_AK5_PF_rhogrid, &b_GroomedJet_AK5_PF_rhogrid);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_pt_uncorr", GroomedJet_AK5_PFCHS_pt_uncorr, &b_GroomedJet_AK5_PFCHS_pt_uncorr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_number_jet_central", &GroomedJet_AK5_PFCHS_number_jet_central, &b_GroomedJet_AK5_PFCHS_number_jet_central);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_mass_uncorr", GroomedJet_AK5_PFCHS_mass_uncorr, &b_GroomedJet_AK5_PFCHS_mass_uncorr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_mass_tr_uncorr", GroomedJet_AK5_PFCHS_mass_tr_uncorr, &b_GroomedJet_AK5_PFCHS_mass_tr_uncorr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_mass_ft_uncorr", GroomedJet_AK5_PFCHS_mass_ft_uncorr, &b_GroomedJet_AK5_PFCHS_mass_ft_uncorr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_mass_pr_uncorr", GroomedJet_AK5_PFCHS_mass_pr_uncorr, &b_GroomedJet_AK5_PFCHS_mass_pr_uncorr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_tau2tau1", GroomedJet_AK5_PFCHS_tau2tau1, &b_GroomedJet_AK5_PFCHS_tau2tau1);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_tau1", GroomedJet_AK5_PFCHS_tau1, &b_GroomedJet_AK5_PFCHS_tau1);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_tau2", GroomedJet_AK5_PFCHS_tau2, &b_GroomedJet_AK5_PFCHS_tau2);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_tau3", GroomedJet_AK5_PFCHS_tau3, &b_GroomedJet_AK5_PFCHS_tau3);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_tau4", GroomedJet_AK5_PFCHS_tau4, &b_GroomedJet_AK5_PFCHS_tau4);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_massdrop_pr_uncorr", GroomedJet_AK5_PFCHS_massdrop_pr_uncorr, &b_GroomedJet_AK5_PFCHS_massdrop_pr_uncorr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_pt", GroomedJet_AK5_PFCHS_pt, &b_GroomedJet_AK5_PFCHS_pt);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_eta", GroomedJet_AK5_PFCHS_eta, &b_GroomedJet_AK5_PFCHS_eta);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_phi", GroomedJet_AK5_PFCHS_phi, &b_GroomedJet_AK5_PFCHS_phi);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_e", GroomedJet_AK5_PFCHS_e, &b_GroomedJet_AK5_PFCHS_e);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_pt_L1_rhoSW", GroomedJet_AK5_PFCHS_pt_L1_rhoSW, &b_GroomedJet_AK5_PFCHS_pt_L1_rhoSW);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_pt_L1_rhoHand", GroomedJet_AK5_PFCHS_pt_L1_rhoHand, &b_GroomedJet_AK5_PFCHS_pt_L1_rhoHand);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_pt_L1_rhoHand2", GroomedJet_AK5_PFCHS_pt_L1_rhoHand2, &b_GroomedJet_AK5_PFCHS_pt_L1_rhoHand2);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_pt_L1_rhoGrid", GroomedJet_AK5_PFCHS_pt_L1_rhoGrid, &b_GroomedJet_AK5_PFCHS_pt_L1_rhoGrid);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_pt_tr_uncorr", GroomedJet_AK5_PFCHS_pt_tr_uncorr, &b_GroomedJet_AK5_PFCHS_pt_tr_uncorr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_pt_tr", GroomedJet_AK5_PFCHS_pt_tr, &b_GroomedJet_AK5_PFCHS_pt_tr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_eta_tr", GroomedJet_AK5_PFCHS_eta_tr, &b_GroomedJet_AK5_PFCHS_eta_tr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_phi_tr", GroomedJet_AK5_PFCHS_phi_tr, &b_GroomedJet_AK5_PFCHS_phi_tr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_e_tr", GroomedJet_AK5_PFCHS_e_tr, &b_GroomedJet_AK5_PFCHS_e_tr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_pt_ft_uncorr", GroomedJet_AK5_PFCHS_pt_ft_uncorr, &b_GroomedJet_AK5_PFCHS_pt_ft_uncorr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_pt_ft", GroomedJet_AK5_PFCHS_pt_ft, &b_GroomedJet_AK5_PFCHS_pt_ft);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_eta_ft", GroomedJet_AK5_PFCHS_eta_ft, &b_GroomedJet_AK5_PFCHS_eta_ft);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_phi_ft", GroomedJet_AK5_PFCHS_phi_ft, &b_GroomedJet_AK5_PFCHS_phi_ft);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_e_ft", GroomedJet_AK5_PFCHS_e_ft, &b_GroomedJet_AK5_PFCHS_e_ft);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_pt_pr_uncorr", GroomedJet_AK5_PFCHS_pt_pr_uncorr, &b_GroomedJet_AK5_PFCHS_pt_pr_uncorr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_pt_pr", GroomedJet_AK5_PFCHS_pt_pr, &b_GroomedJet_AK5_PFCHS_pt_pr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_eta_pr", GroomedJet_AK5_PFCHS_eta_pr, &b_GroomedJet_AK5_PFCHS_eta_pr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_phi_pr", GroomedJet_AK5_PFCHS_phi_pr, &b_GroomedJet_AK5_PFCHS_phi_pr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_e_pr", GroomedJet_AK5_PFCHS_e_pr, &b_GroomedJet_AK5_PFCHS_e_pr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_prsubjet1_px", GroomedJet_AK5_PFCHS_prsubjet1_px, &b_GroomedJet_AK5_PFCHS_prsubjet1_px);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_prsubjet1_py", GroomedJet_AK5_PFCHS_prsubjet1_py, &b_GroomedJet_AK5_PFCHS_prsubjet1_py);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_prsubjet1_pz", GroomedJet_AK5_PFCHS_prsubjet1_pz, &b_GroomedJet_AK5_PFCHS_prsubjet1_pz);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_prsubjet1_e", GroomedJet_AK5_PFCHS_prsubjet1_e, &b_GroomedJet_AK5_PFCHS_prsubjet1_e);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_prsubjet2_px", GroomedJet_AK5_PFCHS_prsubjet2_px, &b_GroomedJet_AK5_PFCHS_prsubjet2_px);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_prsubjet2_py", GroomedJet_AK5_PFCHS_prsubjet2_py, &b_GroomedJet_AK5_PFCHS_prsubjet2_py);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_prsubjet2_pz", GroomedJet_AK5_PFCHS_prsubjet2_pz, &b_GroomedJet_AK5_PFCHS_prsubjet2_pz);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_prsubjet2_e", GroomedJet_AK5_PFCHS_prsubjet2_e, &b_GroomedJet_AK5_PFCHS_prsubjet2_e);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_mass", GroomedJet_AK5_PFCHS_mass, &b_GroomedJet_AK5_PFCHS_mass);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_mass_tr", GroomedJet_AK5_PFCHS_mass_tr, &b_GroomedJet_AK5_PFCHS_mass_tr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_mass_ft", GroomedJet_AK5_PFCHS_mass_ft, &b_GroomedJet_AK5_PFCHS_mass_ft);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_mass_pr", GroomedJet_AK5_PFCHS_mass_pr, &b_GroomedJet_AK5_PFCHS_mass_pr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_massdrop_pr", GroomedJet_AK5_PFCHS_massdrop_pr, &b_GroomedJet_AK5_PFCHS_massdrop_pr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_area", GroomedJet_AK5_PFCHS_area, &b_GroomedJet_AK5_PFCHS_area);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_area_tr", GroomedJet_AK5_PFCHS_area_tr, &b_GroomedJet_AK5_PFCHS_area_tr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_area_ft", GroomedJet_AK5_PFCHS_area_ft, &b_GroomedJet_AK5_PFCHS_area_ft);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_area_pr", GroomedJet_AK5_PFCHS_area_pr, &b_GroomedJet_AK5_PFCHS_area_pr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_jetconstituents", GroomedJet_AK5_PFCHS_jetconstituents, &b_GroomedJet_AK5_PFCHS_jetconstituents);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_jetcharge", GroomedJet_AK5_PFCHS_jetcharge, &b_GroomedJet_AK5_PFCHS_jetcharge);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_rcores", GroomedJet_AK5_PFCHS_rcores, &b_GroomedJet_AK5_PFCHS_rcores);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_ptcores", GroomedJet_AK5_PFCHS_ptcores, &b_GroomedJet_AK5_PFCHS_ptcores);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_planarflow", GroomedJet_AK5_PFCHS_planarflow, &b_GroomedJet_AK5_PFCHS_planarflow);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_qjetmass", GroomedJet_AK5_PFCHS_qjetmass, &b_GroomedJet_AK5_PFCHS_qjetmass);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_qjetmassdrop", GroomedJet_AK5_PFCHS_qjetmassdrop, &b_GroomedJet_AK5_PFCHS_qjetmassdrop);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_constituents0_eta", GroomedJet_AK5_PFCHS_constituents0_eta, &b_GroomedJet_AK5_PFCHS_constituents0_eta);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_constituents0_phi", GroomedJet_AK5_PFCHS_constituents0_phi, &b_GroomedJet_AK5_PFCHS_constituents0_phi);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_constituents0_e", GroomedJet_AK5_PFCHS_constituents0_e, &b_GroomedJet_AK5_PFCHS_constituents0_e);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_nconstituents0", &GroomedJet_AK5_PFCHS_nconstituents0, &b_GroomedJet_AK5_PFCHS_nconstituents0);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_constituents0pr_eta", GroomedJet_AK5_PFCHS_constituents0pr_eta, &b_GroomedJet_AK5_PFCHS_constituents0pr_eta);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_constituents0pr_phi", GroomedJet_AK5_PFCHS_constituents0pr_phi, &b_GroomedJet_AK5_PFCHS_constituents0pr_phi);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_constituents0pr_e", GroomedJet_AK5_PFCHS_constituents0pr_e, &b_GroomedJet_AK5_PFCHS_constituents0pr_e);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_nconstituents0pr", &GroomedJet_AK5_PFCHS_nconstituents0pr, &b_GroomedJet_AK5_PFCHS_nconstituents0pr);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_rhoSW", &GroomedJet_AK5_PFCHS_rhoSW, &b_GroomedJet_AK5_PFCHS_rhoSW);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_rhohand", &GroomedJet_AK5_PFCHS_rhohand, &b_GroomedJet_AK5_PFCHS_rhohand);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_rhohand2", &GroomedJet_AK5_PFCHS_rhohand2, &b_GroomedJet_AK5_PFCHS_rhohand2);
   fChain->SetBranchAddress("GroomedJet_AK5_PFCHS_rhogrid", &GroomedJet_AK5_PFCHS_rhogrid, &b_GroomedJet_AK5_PFCHS_rhogrid);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_pt_uncorr", GenGroomedJet_AK5_GEN_pt_uncorr, &b_GenGroomedJet_AK5_GEN_pt_uncorr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_number_jet_central", &GenGroomedJet_AK5_GEN_number_jet_central, &b_GenGroomedJet_AK5_GEN_number_jet_central);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_mass_uncorr", GenGroomedJet_AK5_GEN_mass_uncorr, &b_GenGroomedJet_AK5_GEN_mass_uncorr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_mass_tr_uncorr", GenGroomedJet_AK5_GEN_mass_tr_uncorr, &b_GenGroomedJet_AK5_GEN_mass_tr_uncorr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_mass_ft_uncorr", GenGroomedJet_AK5_GEN_mass_ft_uncorr, &b_GenGroomedJet_AK5_GEN_mass_ft_uncorr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_mass_pr_uncorr", GenGroomedJet_AK5_GEN_mass_pr_uncorr, &b_GenGroomedJet_AK5_GEN_mass_pr_uncorr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_tau2tau1", GenGroomedJet_AK5_GEN_tau2tau1, &b_GenGroomedJet_AK5_GEN_tau2tau1);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_tau1", GenGroomedJet_AK5_GEN_tau1, &b_GenGroomedJet_AK5_GEN_tau1);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_tau2", GenGroomedJet_AK5_GEN_tau2, &b_GenGroomedJet_AK5_GEN_tau2);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_tau3", GenGroomedJet_AK5_GEN_tau3, &b_GenGroomedJet_AK5_GEN_tau3);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_tau4", GenGroomedJet_AK5_GEN_tau4, &b_GenGroomedJet_AK5_GEN_tau4);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_massdrop_pr_uncorr", GenGroomedJet_AK5_GEN_massdrop_pr_uncorr, &b_GenGroomedJet_AK5_GEN_massdrop_pr_uncorr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_pt", GenGroomedJet_AK5_GEN_pt, &b_GenGroomedJet_AK5_GEN_pt);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_eta", GenGroomedJet_AK5_GEN_eta, &b_GenGroomedJet_AK5_GEN_eta);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_phi", GenGroomedJet_AK5_GEN_phi, &b_GenGroomedJet_AK5_GEN_phi);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_e", GenGroomedJet_AK5_GEN_e, &b_GenGroomedJet_AK5_GEN_e);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_pt_L1_rhoSW", GenGroomedJet_AK5_GEN_pt_L1_rhoSW, &b_GenGroomedJet_AK5_GEN_pt_L1_rhoSW);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_pt_L1_rhoHand", GenGroomedJet_AK5_GEN_pt_L1_rhoHand, &b_GenGroomedJet_AK5_GEN_pt_L1_rhoHand);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_pt_L1_rhoHand2", GenGroomedJet_AK5_GEN_pt_L1_rhoHand2, &b_GenGroomedJet_AK5_GEN_pt_L1_rhoHand2);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_pt_L1_rhoGrid", GenGroomedJet_AK5_GEN_pt_L1_rhoGrid, &b_GenGroomedJet_AK5_GEN_pt_L1_rhoGrid);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_pt_tr_uncorr", GenGroomedJet_AK5_GEN_pt_tr_uncorr, &b_GenGroomedJet_AK5_GEN_pt_tr_uncorr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_pt_tr", GenGroomedJet_AK5_GEN_pt_tr, &b_GenGroomedJet_AK5_GEN_pt_tr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_eta_tr", GenGroomedJet_AK5_GEN_eta_tr, &b_GenGroomedJet_AK5_GEN_eta_tr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_phi_tr", GenGroomedJet_AK5_GEN_phi_tr, &b_GenGroomedJet_AK5_GEN_phi_tr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_e_tr", GenGroomedJet_AK5_GEN_e_tr, &b_GenGroomedJet_AK5_GEN_e_tr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_pt_ft_uncorr", GenGroomedJet_AK5_GEN_pt_ft_uncorr, &b_GenGroomedJet_AK5_GEN_pt_ft_uncorr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_pt_ft", GenGroomedJet_AK5_GEN_pt_ft, &b_GenGroomedJet_AK5_GEN_pt_ft);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_eta_ft", GenGroomedJet_AK5_GEN_eta_ft, &b_GenGroomedJet_AK5_GEN_eta_ft);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_phi_ft", GenGroomedJet_AK5_GEN_phi_ft, &b_GenGroomedJet_AK5_GEN_phi_ft);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_e_ft", GenGroomedJet_AK5_GEN_e_ft, &b_GenGroomedJet_AK5_GEN_e_ft);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_pt_pr_uncorr", GenGroomedJet_AK5_GEN_pt_pr_uncorr, &b_GenGroomedJet_AK5_GEN_pt_pr_uncorr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_pt_pr", GenGroomedJet_AK5_GEN_pt_pr, &b_GenGroomedJet_AK5_GEN_pt_pr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_eta_pr", GenGroomedJet_AK5_GEN_eta_pr, &b_GenGroomedJet_AK5_GEN_eta_pr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_phi_pr", GenGroomedJet_AK5_GEN_phi_pr, &b_GenGroomedJet_AK5_GEN_phi_pr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_e_pr", GenGroomedJet_AK5_GEN_e_pr, &b_GenGroomedJet_AK5_GEN_e_pr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_prsubjet1_px", GenGroomedJet_AK5_GEN_prsubjet1_px, &b_GenGroomedJet_AK5_GEN_prsubjet1_px);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_prsubjet1_py", GenGroomedJet_AK5_GEN_prsubjet1_py, &b_GenGroomedJet_AK5_GEN_prsubjet1_py);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_prsubjet1_pz", GenGroomedJet_AK5_GEN_prsubjet1_pz, &b_GenGroomedJet_AK5_GEN_prsubjet1_pz);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_prsubjet1_e", GenGroomedJet_AK5_GEN_prsubjet1_e, &b_GenGroomedJet_AK5_GEN_prsubjet1_e);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_prsubjet2_px", GenGroomedJet_AK5_GEN_prsubjet2_px, &b_GenGroomedJet_AK5_GEN_prsubjet2_px);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_prsubjet2_py", GenGroomedJet_AK5_GEN_prsubjet2_py, &b_GenGroomedJet_AK5_GEN_prsubjet2_py);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_prsubjet2_pz", GenGroomedJet_AK5_GEN_prsubjet2_pz, &b_GenGroomedJet_AK5_GEN_prsubjet2_pz);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_prsubjet2_e", GenGroomedJet_AK5_GEN_prsubjet2_e, &b_GenGroomedJet_AK5_GEN_prsubjet2_e);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_mass", GenGroomedJet_AK5_GEN_mass, &b_GenGroomedJet_AK5_GEN_mass);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_mass_tr", GenGroomedJet_AK5_GEN_mass_tr, &b_GenGroomedJet_AK5_GEN_mass_tr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_mass_ft", GenGroomedJet_AK5_GEN_mass_ft, &b_GenGroomedJet_AK5_GEN_mass_ft);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_mass_pr", GenGroomedJet_AK5_GEN_mass_pr, &b_GenGroomedJet_AK5_GEN_mass_pr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_massdrop_pr", GenGroomedJet_AK5_GEN_massdrop_pr, &b_GenGroomedJet_AK5_GEN_massdrop_pr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_area", GenGroomedJet_AK5_GEN_area, &b_GenGroomedJet_AK5_GEN_area);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_area_tr", GenGroomedJet_AK5_GEN_area_tr, &b_GenGroomedJet_AK5_GEN_area_tr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_area_ft", GenGroomedJet_AK5_GEN_area_ft, &b_GenGroomedJet_AK5_GEN_area_ft);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_area_pr", GenGroomedJet_AK5_GEN_area_pr, &b_GenGroomedJet_AK5_GEN_area_pr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_jetconstituents", GenGroomedJet_AK5_GEN_jetconstituents, &b_GenGroomedJet_AK5_GEN_jetconstituents);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_jetcharge", GenGroomedJet_AK5_GEN_jetcharge, &b_GenGroomedJet_AK5_GEN_jetcharge);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_rcores", GenGroomedJet_AK5_GEN_rcores, &b_GenGroomedJet_AK5_GEN_rcores);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_ptcores", GenGroomedJet_AK5_GEN_ptcores, &b_GenGroomedJet_AK5_GEN_ptcores);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_planarflow", GenGroomedJet_AK5_GEN_planarflow, &b_GenGroomedJet_AK5_GEN_planarflow);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_qjetmass", GenGroomedJet_AK5_GEN_qjetmass, &b_GenGroomedJet_AK5_GEN_qjetmass);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_qjetmassdrop", GenGroomedJet_AK5_GEN_qjetmassdrop, &b_GenGroomedJet_AK5_GEN_qjetmassdrop);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_constituents0_eta", GenGroomedJet_AK5_GEN_constituents0_eta, &b_GenGroomedJet_AK5_GEN_constituents0_eta);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_constituents0_phi", GenGroomedJet_AK5_GEN_constituents0_phi, &b_GenGroomedJet_AK5_GEN_constituents0_phi);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_constituents0_e", GenGroomedJet_AK5_GEN_constituents0_e, &b_GenGroomedJet_AK5_GEN_constituents0_e);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_nconstituents0", &GenGroomedJet_AK5_GEN_nconstituents0, &b_GenGroomedJet_AK5_GEN_nconstituents0);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_constituents0pr_eta", GenGroomedJet_AK5_GEN_constituents0pr_eta, &b_GenGroomedJet_AK5_GEN_constituents0pr_eta);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_constituents0pr_phi", GenGroomedJet_AK5_GEN_constituents0pr_phi, &b_GenGroomedJet_AK5_GEN_constituents0pr_phi);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_constituents0pr_e", GenGroomedJet_AK5_GEN_constituents0pr_e, &b_GenGroomedJet_AK5_GEN_constituents0pr_e);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_nconstituents0pr", &GenGroomedJet_AK5_GEN_nconstituents0pr, &b_GenGroomedJet_AK5_GEN_nconstituents0pr);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_rhoSW", &GenGroomedJet_AK5_GEN_rhoSW, &b_GenGroomedJet_AK5_GEN_rhoSW);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_rhohand", &GenGroomedJet_AK5_GEN_rhohand, &b_GenGroomedJet_AK5_GEN_rhohand);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_rhohand2", &GenGroomedJet_AK5_GEN_rhohand2, &b_GenGroomedJet_AK5_GEN_rhohand2);
   fChain->SetBranchAddress("GenGroomedJet_AK5_GEN_rhogrid", &GenGroomedJet_AK5_GEN_rhogrid, &b_GenGroomedJet_AK5_GEN_rhogrid);
   fChain->SetBranchAddress("Z_mass", &Z_mass, &b_Z_mass);
   fChain->SetBranchAddress("Z_mt", &Z_mt, &b_Z_mt);
   fChain->SetBranchAddress("Z_mtMVA", &Z_mtMVA, &b_Z_mtMVA);
   fChain->SetBranchAddress("Z_px", &Z_px, &b_Z_px);
   fChain->SetBranchAddress("Z_py", &Z_py, &b_Z_py);
   fChain->SetBranchAddress("Z_pz", &Z_pz, &b_Z_pz);
   fChain->SetBranchAddress("Z_e", &Z_e, &b_Z_e);
   fChain->SetBranchAddress("Z_pt", &Z_pt, &b_Z_pt);
   fChain->SetBranchAddress("Z_et", &Z_et, &b_Z_et);
   fChain->SetBranchAddress("Z_eta", &Z_eta, &b_Z_eta);
   fChain->SetBranchAddress("Z_phi", &Z_phi, &b_Z_phi);
   fChain->SetBranchAddress("Z_vx", &Z_vx, &b_Z_vx);
   fChain->SetBranchAddress("Z_vy", &Z_vy, &b_Z_vy);
   fChain->SetBranchAddress("Z_vz", &Z_vz, &b_Z_vz);
   fChain->SetBranchAddress("Z_y", &Z_y, &b_Z_y);
   fChain->SetBranchAddress("Z_muplus_px", &Z_muplus_px, &b_Z_muplus_px);
   fChain->SetBranchAddress("Z_muplus_py", &Z_muplus_py, &b_Z_muplus_py);
   fChain->SetBranchAddress("Z_muplus_pz", &Z_muplus_pz, &b_Z_muplus_pz);
   fChain->SetBranchAddress("Z_muplus_e", &Z_muplus_e, &b_Z_muplus_e);
   fChain->SetBranchAddress("Z_muplus_pt", &Z_muplus_pt, &b_Z_muplus_pt);
   fChain->SetBranchAddress("Z_muplus_et", &Z_muplus_et, &b_Z_muplus_et);
   fChain->SetBranchAddress("Z_muplus_eta", &Z_muplus_eta, &b_Z_muplus_eta);
   fChain->SetBranchAddress("Z_muplus_theta", &Z_muplus_theta, &b_Z_muplus_theta);
   fChain->SetBranchAddress("Z_muplus_phi", &Z_muplus_phi, &b_Z_muplus_phi);
   fChain->SetBranchAddress("Z_muplus_charge", &Z_muplus_charge, &b_Z_muplus_charge);
   fChain->SetBranchAddress("Z_muplus_vx", &Z_muplus_vx, &b_Z_muplus_vx);
   fChain->SetBranchAddress("Z_muplus_vy", &Z_muplus_vy, &b_Z_muplus_vy);
   fChain->SetBranchAddress("Z_muplus_vz", &Z_muplus_vz, &b_Z_muplus_vz);
   fChain->SetBranchAddress("Z_muplus_y", &Z_muplus_y, &b_Z_muplus_y);
   fChain->SetBranchAddress("Z_muplus_trackiso", &Z_muplus_trackiso, &b_Z_muplus_trackiso);
   fChain->SetBranchAddress("Z_muplus_hcaliso", &Z_muplus_hcaliso, &b_Z_muplus_hcaliso);
   fChain->SetBranchAddress("Z_muplus_ecaliso", &Z_muplus_ecaliso, &b_Z_muplus_ecaliso);
   fChain->SetBranchAddress("Z_muplus_type", &Z_muplus_type, &b_Z_muplus_type);
   fChain->SetBranchAddress("Z_muplus_numberOfChambers", &Z_muplus_numberOfChambers, &b_Z_muplus_numberOfChambers);
   fChain->SetBranchAddress("Z_muplus_numberOfMatches", &Z_muplus_numberOfMatches, &b_Z_muplus_numberOfMatches);
   fChain->SetBranchAddress("Z_muplus_d0bsp", &Z_muplus_d0bsp, &b_Z_muplus_d0bsp);
   fChain->SetBranchAddress("Z_muplus_dz000", &Z_muplus_dz000, &b_Z_muplus_dz000);
   fChain->SetBranchAddress("Z_muplus_dzPV", &Z_muplus_dzPV, &b_Z_muplus_dzPV);
   fChain->SetBranchAddress("Z_muplus_pfiso_sumChargedHadronPt", &Z_muplus_pfiso_sumChargedHadronPt, &b_Z_muplus_pfiso_sumChargedHadronPt);
   fChain->SetBranchAddress("Z_muplus_pfiso_sumChargedParticlePt", &Z_muplus_pfiso_sumChargedParticlePt, &b_Z_muplus_pfiso_sumChargedParticlePt);
   fChain->SetBranchAddress("Z_muplus_pfiso_sumNeutralHadronEt", &Z_muplus_pfiso_sumNeutralHadronEt, &b_Z_muplus_pfiso_sumNeutralHadronEt);
   fChain->SetBranchAddress("Z_muplus_pfiso_sumPhotonEt", &Z_muplus_pfiso_sumPhotonEt, &b_Z_muplus_pfiso_sumPhotonEt);
   fChain->SetBranchAddress("Z_muplus_pfiso_sumPUPt", &Z_muplus_pfiso_sumPUPt, &b_Z_muplus_pfiso_sumPUPt);
   fChain->SetBranchAddress("Z_muminus_px", &Z_muminus_px, &b_Z_muminus_px);
   fChain->SetBranchAddress("Z_muminus_py", &Z_muminus_py, &b_Z_muminus_py);
   fChain->SetBranchAddress("Z_muminus_pz", &Z_muminus_pz, &b_Z_muminus_pz);
   fChain->SetBranchAddress("Z_muminus_e", &Z_muminus_e, &b_Z_muminus_e);
   fChain->SetBranchAddress("Z_muminus_pt", &Z_muminus_pt, &b_Z_muminus_pt);
   fChain->SetBranchAddress("Z_muminus_et", &Z_muminus_et, &b_Z_muminus_et);
   fChain->SetBranchAddress("Z_muminus_eta", &Z_muminus_eta, &b_Z_muminus_eta);
   fChain->SetBranchAddress("Z_muminus_theta", &Z_muminus_theta, &b_Z_muminus_theta);
   fChain->SetBranchAddress("Z_muminus_phi", &Z_muminus_phi, &b_Z_muminus_phi);
   fChain->SetBranchAddress("Z_muminus_charge", &Z_muminus_charge, &b_Z_muminus_charge);
   fChain->SetBranchAddress("Z_muminus_vx", &Z_muminus_vx, &b_Z_muminus_vx);
   fChain->SetBranchAddress("Z_muminus_vy", &Z_muminus_vy, &b_Z_muminus_vy);
   fChain->SetBranchAddress("Z_muminus_vz", &Z_muminus_vz, &b_Z_muminus_vz);
   fChain->SetBranchAddress("Z_muminus_y", &Z_muminus_y, &b_Z_muminus_y);
   fChain->SetBranchAddress("Z_muminus_trackiso", &Z_muminus_trackiso, &b_Z_muminus_trackiso);
   fChain->SetBranchAddress("Z_muminus_hcaliso", &Z_muminus_hcaliso, &b_Z_muminus_hcaliso);
   fChain->SetBranchAddress("Z_muminus_ecaliso", &Z_muminus_ecaliso, &b_Z_muminus_ecaliso);
   fChain->SetBranchAddress("Z_muminus_type", &Z_muminus_type, &b_Z_muminus_type);
   fChain->SetBranchAddress("Z_muminus_numberOfChambers", &Z_muminus_numberOfChambers, &b_Z_muminus_numberOfChambers);
   fChain->SetBranchAddress("Z_muminus_numberOfMatches", &Z_muminus_numberOfMatches, &b_Z_muminus_numberOfMatches);
   fChain->SetBranchAddress("Z_muminus_d0bsp", &Z_muminus_d0bsp, &b_Z_muminus_d0bsp);
   fChain->SetBranchAddress("Z_muminus_dz000", &Z_muminus_dz000, &b_Z_muminus_dz000);
   fChain->SetBranchAddress("Z_muminus_dzPV", &Z_muminus_dzPV, &b_Z_muminus_dzPV);
   fChain->SetBranchAddress("Z_muminus_pfiso_sumChargedHadronPt", &Z_muminus_pfiso_sumChargedHadronPt, &b_Z_muminus_pfiso_sumChargedHadronPt);
   fChain->SetBranchAddress("Z_muminus_pfiso_sumChargedParticlePt", &Z_muminus_pfiso_sumChargedParticlePt, &b_Z_muminus_pfiso_sumChargedParticlePt);
   fChain->SetBranchAddress("Z_muminus_pfiso_sumNeutralHadronEt", &Z_muminus_pfiso_sumNeutralHadronEt, &b_Z_muminus_pfiso_sumNeutralHadronEt);
   fChain->SetBranchAddress("Z_muminus_pfiso_sumPhotonEt", &Z_muminus_pfiso_sumPhotonEt, &b_Z_muminus_pfiso_sumPhotonEt);
   fChain->SetBranchAddress("Z_muminus_pfiso_sumPUPt", &Z_muminus_pfiso_sumPUPt, &b_Z_muminus_pfiso_sumPUPt);
   fChain->SetBranchAddress("Z_Photon_pt_gen", &Z_Photon_pt_gen, &b_Z_Photon_pt_gen);
   fChain->SetBranchAddress("Z_muplus_px_gen", &Z_muplus_px_gen, &b_Z_muplus_px_gen);
   fChain->SetBranchAddress("Z_muplus_py_gen", &Z_muplus_py_gen, &b_Z_muplus_py_gen);
   fChain->SetBranchAddress("Z_muplus_pz_gen", &Z_muplus_pz_gen, &b_Z_muplus_pz_gen);
   fChain->SetBranchAddress("Z_muplus_e_gen", &Z_muplus_e_gen, &b_Z_muplus_e_gen);
   fChain->SetBranchAddress("Z_muplus_pt_gen", &Z_muplus_pt_gen, &b_Z_muplus_pt_gen);
   fChain->SetBranchAddress("Z_muplus_et_gen", &Z_muplus_et_gen, &b_Z_muplus_et_gen);
   fChain->SetBranchAddress("Z_muplus_eta_gen", &Z_muplus_eta_gen, &b_Z_muplus_eta_gen);
   fChain->SetBranchAddress("Z_muplus_theta_gen", &Z_muplus_theta_gen, &b_Z_muplus_theta_gen);
   fChain->SetBranchAddress("Z_muplus_phi_gen", &Z_muplus_phi_gen, &b_Z_muplus_phi_gen);
   fChain->SetBranchAddress("Z_muplus_charge_gen", &Z_muplus_charge_gen, &b_Z_muplus_charge_gen);
   fChain->SetBranchAddress("Z_muplus_vx_gen", &Z_muplus_vx_gen, &b_Z_muplus_vx_gen);
   fChain->SetBranchAddress("Z_muplus_vy_gen", &Z_muplus_vy_gen, &b_Z_muplus_vy_gen);
   fChain->SetBranchAddress("Z_muplus_vz_gen", &Z_muplus_vz_gen, &b_Z_muplus_vz_gen);
   fChain->SetBranchAddress("Z_muplus_y_gen", &Z_muplus_y_gen, &b_Z_muplus_y_gen);
   fChain->SetBranchAddress("Z_muminus_px_gen", &Z_muminus_px_gen, &b_Z_muminus_px_gen);
   fChain->SetBranchAddress("Z_muminus_py_gen", &Z_muminus_py_gen, &b_Z_muminus_py_gen);
   fChain->SetBranchAddress("Z_muminus_pz_gen", &Z_muminus_pz_gen, &b_Z_muminus_pz_gen);
   fChain->SetBranchAddress("Z_muminus_e_gen", &Z_muminus_e_gen, &b_Z_muminus_e_gen);
   fChain->SetBranchAddress("Z_muminus_pt_gen", &Z_muminus_pt_gen, &b_Z_muminus_pt_gen);
   fChain->SetBranchAddress("Z_muminus_et_gen", &Z_muminus_et_gen, &b_Z_muminus_et_gen);
   fChain->SetBranchAddress("Z_muminus_eta_gen", &Z_muminus_eta_gen, &b_Z_muminus_eta_gen);
   fChain->SetBranchAddress("Z_muminus_theta_gen", &Z_muminus_theta_gen, &b_Z_muminus_theta_gen);
   fChain->SetBranchAddress("Z_muminus_phi_gen", &Z_muminus_phi_gen, &b_Z_muminus_phi_gen);
   fChain->SetBranchAddress("Z_muminus_charge_gen", &Z_muminus_charge_gen, &b_Z_muminus_charge_gen);
   fChain->SetBranchAddress("Z_muminus_vx_gen", &Z_muminus_vx_gen, &b_Z_muminus_vx_gen);
   fChain->SetBranchAddress("Z_muminus_vy_gen", &Z_muminus_vy_gen, &b_Z_muminus_vy_gen);
   fChain->SetBranchAddress("Z_muminus_vz_gen", &Z_muminus_vz_gen, &b_Z_muminus_vz_gen);
   fChain->SetBranchAddress("Z_muminus_y_gen", &Z_muminus_y_gen, &b_Z_muminus_y_gen);
   fChain->SetBranchAddress("event_runNo", &event_runNo, &b_event_runNo);
   fChain->SetBranchAddress("event_evtNo", &event_evtNo, &b_event_evtNo);
   fChain->SetBranchAddress("event_lumi", &event_lumi, &b_event_lumi);
   fChain->SetBranchAddress("event_bunch", &event_bunch, &b_event_bunch);
   fChain->SetBranchAddress("event_nPV", &event_nPV, &b_event_nPV);
   fChain->SetBranchAddress("event_met_pfmet", &event_met_pfmet, &b_event_met_pfmet);
   fChain->SetBranchAddress("event_met_pfsumet", &event_met_pfsumet, &b_event_met_pfsumet);
   fChain->SetBranchAddress("event_met_pfmetsignificance", &event_met_pfmetsignificance, &b_event_met_pfmetsignificance);
   fChain->SetBranchAddress("event_met_pfmetPhi", &event_met_pfmetPhi, &b_event_met_pfmetPhi);
   fChain->SetBranchAddress("event_metMVA_met", &event_metMVA_met, &b_event_metMVA_met);
   fChain->SetBranchAddress("event_metMVA_sumet", &event_metMVA_sumet, &b_event_metMVA_sumet);
   fChain->SetBranchAddress("event_metMVA_metsignificance", &event_metMVA_metsignificance, &b_event_metMVA_metsignificance);
   fChain->SetBranchAddress("event_metMVA_metPhi", &event_metMVA_metPhi, &b_event_metMVA_metPhi);
   fChain->SetBranchAddress("event_fastJetRho", &event_fastJetRho, &b_event_fastJetRho);
   fChain->SetBranchAddress("event_met_genmet", &event_met_genmet, &b_event_met_genmet);
   fChain->SetBranchAddress("event_met_gensumet", &event_met_gensumet, &b_event_met_gensumet);
   fChain->SetBranchAddress("event_met_genmetsignificance", &event_met_genmetsignificance, &b_event_met_genmetsignificance);
   fChain->SetBranchAddress("event_met_genmetPhi", &event_met_genmetPhi, &b_event_met_genmetPhi);
   fChain->SetBranchAddress("event_mcPU_totnvtx", &event_mcPU_totnvtx, &b_event_mcPU_totnvtx);
   fChain->SetBranchAddress("event_mcPU_trueInteractions", &event_mcPU_trueInteractions, &b_event_mcPU_trueInteractions);
   fChain->SetBranchAddress("event_mcPU_bx", event_mcPU_bx, &b_event_mcPU_bx);
   fChain->SetBranchAddress("event_mcPU_nvtx", event_mcPU_nvtx, &b_event_mcPU_nvtx);
   Notify();
}

Bool_t MyClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyClass_cxx
