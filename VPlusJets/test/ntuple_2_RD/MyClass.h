//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct  2 07:21:13 2013 by ROOT version 5.32/00
// from TTree ZJet/V+jets Tree
// found on fout: zmumujetsanalysisntuple.root
//////////////////////////////////////////////////////////

#ifndef MyClass_h
#define MyClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TROOT.h>
#include <iostream>
#include <TGraph.h>
#include <TMath.h>
#include <map>
#include "fstream"
using namespace std;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Efficiency_Tool{
	public:
		Int_t number_step;
		std::map< string, Int_t > map_stepname;
		std::map< Int_t, string > map_stepnumber;
		std::vector<Double_t> vect_step_event;
		Efficiency_Tool();
		~Efficiency_Tool(){}

		void Add_Step(string newstepname);
		void Add_Event(string newstepname, Double_t eventweight=1. );
		TH1D Get_Eff_hist();
};

class Table_Tool{ //print out tables into txt
	private:
		vector< TString> xaxis;
		vector< TString> yaxis;
		Int_t iter_x;
		Int_t iter_y;
		Double_t table_value[100][100];
	public:
		Table_Tool();
		~Table_Tool(){};
		void Insert(TString x_name, TString y_name, Double_t value);
		void PrintTable(ofstream &fout);
		//void PrintTable();
};


class RealVarArray{
	public:
		TString name;
		Int_t nbin;
		Double_t xmin;
		Double_t xmax;
		std::vector< Double_t > vect_x;
		std::vector< Double_t > vect_weight;

		RealVarArray(const char* in_name, Int_t in_nbin, Double_t in_xmin, Double_t in_xmax);
		RealVarArray();
		~RealVarArray(){};
};

typedef map<TString, RealVarArray> MAP_REALVARARRAY;


struct RESPONSE{
	TH1D hist;
	TH2D hist2D;
	Double_t mean;
	Double_t rms;
	Double_t correlationfactor;
};
typedef map<TString, RESPONSE> MAP_RESPONSE;

Bool_t find_var_in_map(MAP_REALVARARRAY &map_obs,  TString var_name){
	MAP_REALVARARRAY::iterator it=map_obs.find(var_name);
	if (it == map_obs.end() ) { return 0;}
	else { return 1; }
}
Bool_t find_var_in_map(MAP_RESPONSE &map_obs,  TString var_name){
	MAP_RESPONSE::iterator it=map_obs.find(var_name);
	if (it == map_obs.end() ) { return 0;}
	else { return 1; }
}
class JetCorrectionTool{
	private:
		TString name;
		MAP_REALVARARRAY map_obs; // reco pt, mass, jet shape, gen pt, eta, phi

	public:
		JetCorrectionTool(const char* in_name);
		~JetCorrectionTool(){};

		void addVar( RealVarArray x );
		Bool_t fill( TString var_name, Double_t in_val, Double_t eventweight=1.0, Bool_t debug=0);

		vector< TString > get_allvar_name();
		TH1D get_hist1D(TString var_name);
		TH1D get_hist1D_response( TString xdenominator_var_name,TString xnumerator_var_name , Int_t nbin=100, Double_t xmin=0, Double_t xmax=2. );
		RESPONSE get_response( TString xdenominator_var_name,TString xnumerator_var_name , Int_t nbin=100, Double_t xmin=0, Double_t xmax=2. );
		TH2D get_hist2D( TString x_var_name,TString y_var_name );
		TGraph get_graph( TString x_var_name,TString y_var_name );
		TH1D get_mean_rms_hist( TString x_var_name, TString ydenominator_var_name, TString ynumerator_var_name, Double_t ymin=0, Double_t ymax=2.0 );

		MAP_RESPONSE map_response;

};


class MyClass {
	public :
		TFile *fout;

		Efficiency_Tool efftool;//efficiency tool
		//Bool_t   isBoosted;//boosted Z or not
		Int_t   isBoosted;//boosted Z or not
		TTree          *fChain;   //!pointer to the analyzed TTree or TChain
		Int_t           fCurrent; //!current Tree number in a TChain

		TString      plot_Dir_DateTime;
		TString      FinalState;
		TString      JetType;
		TString      PfType;
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

		Float_t         GroomedJet_pt_uncorr[6];
		Int_t           GroomedJet_number_jet_central;
		Float_t         GroomedJet_mass_uncorr[6];
		Float_t         GroomedJet_mass_tr_uncorr[6];
		Float_t         GroomedJet_mass_ft_uncorr[6];
		Float_t         GroomedJet_mass_pr_uncorr[6];
		Float_t         GroomedJet_tau2tau1[6];
		Float_t         GroomedJet_tau2tau1_shapesubtraction[6];
		Float_t         GroomedJet_tau1[6];
		Float_t         GroomedJet_tau2[6];
		Float_t         GroomedJet_tau3[6];
		Float_t         GroomedJet_tau4[6];
		Float_t         GroomedJet_massdrop_pr_uncorr[6];
		Float_t         GroomedJet_pt[6];
		Float_t         GroomedJet_eta[6];
		Float_t         GroomedJet_phi[6];
		Float_t         GroomedJet_e[6];
		Float_t         GroomedJet_pt_JECL1[6];
		Float_t         GroomedJet_eta_JECL1[6];
		Float_t         GroomedJet_phi_JECL1[6];
		Float_t         GroomedJet_e_JECL1[6];
		Float_t         GroomedJet_pt_L1_rhoSW[6];
		Float_t         GroomedJet_pt_L1_rhoHand[6];
		Float_t         GroomedJet_pt_L1_rhoHand2[6];
		Float_t         GroomedJet_pt_L1_rhoGrid[6];
		Float_t         GroomedJet_pt_rho4A[6];
		Float_t         GroomedJet_pt_shapesubtraction[6];
		Float_t         GroomedJet_pt_tr_uncorr[6];
		Float_t         GroomedJet_pt_tr[6];
		Float_t         GroomedJet_eta_tr[6];
		Float_t         GroomedJet_phi_tr[6];
		Float_t         GroomedJet_e_tr[6];
		Float_t         GroomedJet_pt_ft_uncorr[6];
		Float_t         GroomedJet_pt_ft[6];
		Float_t         GroomedJet_eta_ft[6];
		Float_t         GroomedJet_phi_ft[6];
		Float_t         GroomedJet_e_ft[6];
		Float_t         GroomedJet_pt_pr_uncorr[6];
		Float_t         GroomedJet_pt_pr[6];
		Float_t         GroomedJet_eta_pr[6];
		Float_t         GroomedJet_phi_pr[6];
		Float_t         GroomedJet_e_pr[6];
		Float_t         GroomedJet_prsubjet1_px[6];
		Float_t         GroomedJet_prsubjet1_py[6];
		Float_t         GroomedJet_prsubjet1_pz[6];
		Float_t         GroomedJet_prsubjet1_e[6];
		Float_t         GroomedJet_prsubjet2_px[6];
		Float_t         GroomedJet_prsubjet2_py[6];
		Float_t         GroomedJet_prsubjet2_pz[6];
		Float_t         GroomedJet_prsubjet2_e[6];
		Float_t         GroomedJet_mass[6];
		Float_t         GroomedJet_mass_rhoArea[6];
		Float_t         GroomedJet_mass_rhoGArea[6];
		Float_t         GroomedJet_mass_rho4A[6];
		Float_t         GroomedJet_mass_rhoG4A[6];
		Float_t         GroomedJet_mass_rhom4Am[6];
		Float_t         GroomedJet_mass_shapesubtraction[6];
		Float_t         GroomedJet_mass_trimmingshapesubtraction[6];
		Float_t         GroomedJet_mass_tr[6];
		Float_t         GroomedJet_mass_ft[6];
		Float_t         GroomedJet_mass_pr[6];
		Float_t         GroomedJet_massdrop_pr[6];
		Float_t         GroomedJet_area[6];
		Float_t         GroomedJet_area_tr[6];
		Float_t         GroomedJet_area_ft[6];
		Float_t         GroomedJet_area_pr[6];
		Float_t         GroomedJet_jetconstituents[6];
		Float_t         GroomedJet_jetcharge[6];
		Float_t         GroomedJet_rcores[11][6];
		Float_t         GroomedJet_ptcores[11][6];
		Float_t         GroomedJet_planarflow[11][6];
		Float_t         GroomedJet_qjetmass[50];
		Float_t         GroomedJet_qjetmassdrop[50];
		static const Int_t NUM_JETCLEANSING_DIFFMODE =200; 
		Float_t         GroomedJet_mass_JetCleansing_DiffMode[NUM_JETCLEANSING_DIFFMODE];
		Float_t         GroomedJet_pt_JetCleansing_DiffMode[NUM_JETCLEANSING_DIFFMODE];
		Float_t         GroomedJet_eta_JetCleansing_DiffMode[NUM_JETCLEANSING_DIFFMODE];
		Float_t         GroomedJet_phi_JetCleansing_DiffMode[NUM_JETCLEANSING_DIFFMODE];
		Float_t         GroomedJet_tau2tau1_JetCleansing_DiffMode[NUM_JETCLEANSING_DIFFMODE];
		Float_t         GroomedJet_constituents0_eta[100];
		Float_t         GroomedJet_constituents0_phi[100];
		Float_t         GroomedJet_constituents0_e[100];
		Int_t           GroomedJet_nconstituents0;
		Float_t         GroomedJet_constituents0pr_eta[100];
		Float_t         GroomedJet_constituents0pr_phi[100];
		Float_t         GroomedJet_constituents0pr_e[100];
		Int_t           GroomedJet_nconstituents0pr;
		Double_t        GroomedJet_rhoSW;
		Double_t        GroomedJet_rhohand;
		Double_t        GroomedJet_rhohand2;
		Double_t        GroomedJet_rhogrid;

		Float_t         GenGroomedJet_pt_uncorr[6];
		Int_t           GenGroomedJet_number_jet_central;
		Float_t         GenGroomedJet_mass_uncorr[6];
		Float_t         GenGroomedJet_mass_tr_uncorr[6];
		Float_t         GenGroomedJet_mass_ft_uncorr[6];
		Float_t         GenGroomedJet_mass_pr_uncorr[6];
		Float_t         GenGroomedJet_tau2tau1[6];
		Float_t         GenGroomedJet_tau2tau1_shapesubtraction[6];
		Float_t         GenGroomedJet_tau1[6];
		Float_t         GenGroomedJet_tau2[6];
		Float_t         GenGroomedJet_tau3[6];
		Float_t         GenGroomedJet_tau4[6];
		Float_t         GenGroomedJet_massdrop_pr_uncorr[6];
		Float_t         GenGroomedJet_pt[6];
		Float_t         GenGroomedJet_eta[6];
		Float_t         GenGroomedJet_phi[6];
		Float_t         GenGroomedJet_e[6];
		Float_t         GenGroomedJet_pt_L1_rhoSW[6];
		Float_t         GenGroomedJet_pt_L1_rhoHand[6];
		Float_t         GenGroomedJet_pt_L1_rhoHand2[6];
		Float_t         GenGroomedJet_pt_L1_rhoGrid[6];
		Float_t         GenGroomedJet_pt_rho4A[6];
		Float_t         GenGroomedJet_pt_shapesubtraction[6];
		Float_t         GenGroomedJet_pt_tr_uncorr[6];
		Float_t         GenGroomedJet_pt_tr[6];
		Float_t         GenGroomedJet_eta_tr[6];
		Float_t         GenGroomedJet_phi_tr[6];
		Float_t         GenGroomedJet_e_tr[6];
		Float_t         GenGroomedJet_pt_ft_uncorr[6];
		Float_t         GenGroomedJet_pt_ft[6];
		Float_t         GenGroomedJet_eta_ft[6];
		Float_t         GenGroomedJet_phi_ft[6];
		Float_t         GenGroomedJet_e_ft[6];
		Float_t         GenGroomedJet_pt_pr_uncorr[6];
		Float_t         GenGroomedJet_pt_pr[6];
		Float_t         GenGroomedJet_eta_pr[6];
		Float_t         GenGroomedJet_phi_pr[6];
		Float_t         GenGroomedJet_e_pr[6];
		Float_t         GenGroomedJet_prsubjet1_px[6];
		Float_t         GenGroomedJet_prsubjet1_py[6];
		Float_t         GenGroomedJet_prsubjet1_pz[6];
		Float_t         GenGroomedJet_prsubjet1_e[6];
		Float_t         GenGroomedJet_prsubjet2_px[6];
		Float_t         GenGroomedJet_prsubjet2_py[6];
		Float_t         GenGroomedJet_prsubjet2_pz[6];
		Float_t         GenGroomedJet_prsubjet2_e[6];
		Float_t         GenGroomedJet_mass[6];
		Float_t         GenGroomedJet_mass_rhoArea[6];
		Float_t         GenGroomedJet_mass_rhoGArea[6];
		Float_t         GenGroomedJet_mass_rho4A[6];
		Float_t         GenGroomedJet_mass_rhoG4A[6];
		Float_t         GenGroomedJet_mass_rhom4Am[6];
		Float_t         GenGroomedJet_mass_shapesubtraction[6];
		Float_t         GenGroomedJet_mass_trimmingshapesubtraction[6];
		Float_t         GenGroomedJet_mass_tr[6];
		Float_t         GenGroomedJet_mass_ft[6];
		Float_t         GenGroomedJet_mass_pr[6];
		Float_t         GenGroomedJet_massdrop_pr[6];
		Float_t         GenGroomedJet_area[6];
		Float_t         GenGroomedJet_area_tr[6];
		Float_t         GenGroomedJet_area_ft[6];
		Float_t         GenGroomedJet_area_pr[6];
		Float_t         GenGroomedJet_jetconstituents[6];
		Float_t         GenGroomedJet_jetcharge[6];
		Float_t         GenGroomedJet_rcores[11][6];
		Float_t         GenGroomedJet_ptcores[11][6];
		Float_t         GenGroomedJet_planarflow[11][6];
		Float_t         GenGroomedJet_qjetmass[50];
		Float_t         GenGroomedJet_qjetmassdrop[50];
		Float_t         GenGroomedJet_mass_JetCleansing_DiffMode[NUM_JETCLEANSING_DIFFMODE];
		Float_t         GenGroomedJet_pt_JetCleansing_DiffMode[NUM_JETCLEANSING_DIFFMODE];
		Float_t         GenGroomedJet_eta_JetCleansing_DiffMode[NUM_JETCLEANSING_DIFFMODE];
		Float_t         GenGroomedJet_phi_JetCleansing_DiffMode[NUM_JETCLEANSING_DIFFMODE];
		Float_t         GenGroomedJet_tau2tau1_JetCleansing_DiffMode[NUM_JETCLEANSING_DIFFMODE];
		Float_t         GenGroomedJet_constituents0_eta[100];
		Float_t         GenGroomedJet_constituents0_phi[100];
		Float_t         GenGroomedJet_constituents0_e[100];
		Int_t           GenGroomedJet_nconstituents0;
		Float_t         GenGroomedJet_constituents0pr_eta[100];
		Float_t         GenGroomedJet_constituents0pr_phi[100];
		Float_t         GenGroomedJet_constituents0pr_e[100];
		Int_t           GenGroomedJet_nconstituents0pr;
		Double_t        GenGroomedJet_rhoSW;
		Double_t        GenGroomedJet_rhohand;
		Double_t        GenGroomedJet_rhohand2;
		Double_t        GenGroomedJet_rhogrid;


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
		Float_t         EventWeight;
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

		TBranch        *b_GroomedJet_pt_uncorr;   //!
		TBranch        *b_GroomedJet_number_jet_central;   //!
		TBranch        *b_GroomedJet_mass_uncorr;   //!
		TBranch        *b_GroomedJet_mass_tr_uncorr;   //!
		TBranch        *b_GroomedJet_mass_ft_uncorr;   //!
		TBranch        *b_GroomedJet_mass_pr_uncorr;   //!
		TBranch        *b_GroomedJet_tau2tau1;   //!
		TBranch        *b_GroomedJet_tau2tau1_shapesubtraction;   //!
		TBranch        *b_GroomedJet_tau1;   //!
		TBranch        *b_GroomedJet_tau2;   //!
		TBranch        *b_GroomedJet_tau3;   //!
		TBranch        *b_GroomedJet_tau4;   //!
		TBranch        *b_GroomedJet_massdrop_pr_uncorr;   //!
		TBranch        *b_GroomedJet_pt;   //!
		TBranch        *b_GroomedJet_eta;   //!
		TBranch        *b_GroomedJet_phi;   //!
		TBranch        *b_GroomedJet_e;   //!
		TBranch        *b_GroomedJet_pt_JECL1;   //!
		TBranch        *b_GroomedJet_eta_JECL1;   //!
		TBranch        *b_GroomedJet_phi_JECL1;   //!
		TBranch        *b_GroomedJet_e_JECL1;   //!
		TBranch        *b_GroomedJet_pt_L1_rhoSW;   //!
		TBranch        *b_GroomedJet_pt_L1_rhoHand;   //!
		TBranch        *b_GroomedJet_pt_L1_rhoHand2;   //!
		TBranch        *b_GroomedJet_pt_L1_rhoGrid;   //!
		TBranch        *b_GroomedJet_pt_rho4A;   //!
		TBranch        *b_GroomedJet_pt_shapesubtraction;   //!
		TBranch        *b_GroomedJet_pt_tr_uncorr;   //!
		TBranch        *b_GroomedJet_pt_tr;   //!
		TBranch        *b_GroomedJet_eta_tr;   //!
		TBranch        *b_GroomedJet_phi_tr;   //!
		TBranch        *b_GroomedJet_e_tr;   //!
		TBranch        *b_GroomedJet_pt_ft_uncorr;   //!
		TBranch        *b_GroomedJet_pt_ft;   //!
		TBranch        *b_GroomedJet_eta_ft;   //!
		TBranch        *b_GroomedJet_phi_ft;   //!
		TBranch        *b_GroomedJet_e_ft;   //!
		TBranch        *b_GroomedJet_pt_pr_uncorr;   //!
		TBranch        *b_GroomedJet_pt_pr;   //!
		TBranch        *b_GroomedJet_eta_pr;   //!
		TBranch        *b_GroomedJet_phi_pr;   //!
		TBranch        *b_GroomedJet_e_pr;   //!
		TBranch        *b_GroomedJet_prsubjet1_px;   //!
		TBranch        *b_GroomedJet_prsubjet1_py;   //!
		TBranch        *b_GroomedJet_prsubjet1_pz;   //!
		TBranch        *b_GroomedJet_prsubjet1_e;   //!
		TBranch        *b_GroomedJet_prsubjet2_px;   //!
		TBranch        *b_GroomedJet_prsubjet2_py;   //!
		TBranch        *b_GroomedJet_prsubjet2_pz;   //!
		TBranch        *b_GroomedJet_prsubjet2_e;   //!
		TBranch        *b_GroomedJet_mass;   //!
		TBranch        *b_GroomedJet_mass_rhoArea;   //!
		TBranch        *b_GroomedJet_mass_rhoGArea;   //!
		TBranch        *b_GroomedJet_mass_rho4A;   //!
		TBranch        *b_GroomedJet_mass_rhoG4A;   //!
		TBranch        *b_GroomedJet_mass_rhom4Am;   //!
		TBranch        *b_GroomedJet_mass_shapesubtraction;   //!
		TBranch        *b_GroomedJet_mass_trimmingshapesubtraction;   //!
		TBranch        *b_GroomedJet_mass_tr;   //!
		TBranch        *b_GroomedJet_mass_ft;   //!
		TBranch        *b_GroomedJet_mass_pr;   //!
		TBranch        *b_GroomedJet_massdrop_pr;   //!
		TBranch        *b_GroomedJet_area;   //!
		TBranch        *b_GroomedJet_area_tr;   //!
		TBranch        *b_GroomedJet_area_ft;   //!
		TBranch        *b_GroomedJet_area_pr;   //!
		TBranch        *b_GroomedJet_jetconstituents;   //!
		TBranch        *b_GroomedJet_jetcharge;   //!
		TBranch        *b_GroomedJet_rcores;   //!
		TBranch        *b_GroomedJet_ptcores;   //!
		TBranch        *b_GroomedJet_planarflow;   //!
		TBranch        *b_GroomedJet_qjetmass;   //!
		TBranch        *b_GroomedJet_qjetmassdrop;   //!
		TBranch        *b_GroomedJet_mass_JetCleansing_DiffMode;   //!
		TBranch        *b_GroomedJet_pt_JetCleansing_DiffMode;   //!
		TBranch        *b_GroomedJet_eta_JetCleansing_DiffMode;   //!
		TBranch        *b_GroomedJet_phi_JetCleansing_DiffMode;   //!
		TBranch        *b_GroomedJet_tau2tau1_JetCleansing_DiffMode;   //!
		TBranch        *b_GroomedJet_constituents0_eta;   //!
		TBranch        *b_GroomedJet_constituents0_phi;   //!
		TBranch        *b_GroomedJet_constituents0_e;   //!
		TBranch        *b_GroomedJet_nconstituents0;   //!
		TBranch        *b_GroomedJet_constituents0pr_eta;   //!
		TBranch        *b_GroomedJet_constituents0pr_phi;   //!
		TBranch        *b_GroomedJet_constituents0pr_e;   //!
		TBranch        *b_GroomedJet_nconstituents0pr;   //!
		TBranch        *b_GroomedJet_rhoSW;   //!
		TBranch        *b_GroomedJet_rhohand;   //!
		TBranch        *b_GroomedJet_rhohand2;   //!
		TBranch        *b_GroomedJet_rhogrid;   //!
		TBranch        *b_GenGroomedJet_pt_uncorr;   //!
		TBranch        *b_GenGroomedJet_number_jet_central;   //!
		TBranch        *b_GenGroomedJet_mass_uncorr;   //!
		TBranch        *b_GenGroomedJet_mass_tr_uncorr;   //!
		TBranch        *b_GenGroomedJet_mass_ft_uncorr;   //!
		TBranch        *b_GenGroomedJet_mass_pr_uncorr;   //!
		TBranch        *b_GenGroomedJet_tau2tau1;   //!
		TBranch        *b_GenGroomedJet_tau2tau1_shapesubtraction;   //!
		TBranch        *b_GenGroomedJet_tau1;   //!
		TBranch        *b_GenGroomedJet_tau2;   //!
		TBranch        *b_GenGroomedJet_tau3;   //!
		TBranch        *b_GenGroomedJet_tau4;   //!
		TBranch        *b_GenGroomedJet_massdrop_pr_uncorr;   //!
		TBranch        *b_GenGroomedJet_pt;   //!
		TBranch        *b_GenGroomedJet_eta;   //!
		TBranch        *b_GenGroomedJet_phi;   //!
		TBranch        *b_GenGroomedJet_e;   //!
		TBranch        *b_GenGroomedJet_pt_L1_rhoSW;   //!
		TBranch        *b_GenGroomedJet_pt_L1_rhoHand;   //!
		TBranch        *b_GenGroomedJet_pt_L1_rhoHand2;   //!
		TBranch        *b_GenGroomedJet_pt_L1_rhoGrid;   //!
		TBranch        *b_GenGroomedJet_pt_rho4A;   //!
		TBranch        *b_GenGroomedJet_pt_shapesubtraction;   //!
		TBranch        *b_GenGroomedJet_pt_tr_uncorr;   //!
		TBranch        *b_GenGroomedJet_pt_tr;   //!
		TBranch        *b_GenGroomedJet_eta_tr;   //!
		TBranch        *b_GenGroomedJet_phi_tr;   //!
		TBranch        *b_GenGroomedJet_e_tr;   //!
		TBranch        *b_GenGroomedJet_pt_ft_uncorr;   //!
		TBranch        *b_GenGroomedJet_pt_ft;   //!
		TBranch        *b_GenGroomedJet_eta_ft;   //!
		TBranch        *b_GenGroomedJet_phi_ft;   //!
		TBranch        *b_GenGroomedJet_e_ft;   //!
		TBranch        *b_GenGroomedJet_pt_pr_uncorr;   //!
		TBranch        *b_GenGroomedJet_pt_pr;   //!
		TBranch        *b_GenGroomedJet_eta_pr;   //!
		TBranch        *b_GenGroomedJet_phi_pr;   //!
		TBranch        *b_GenGroomedJet_e_pr;   //!
		TBranch        *b_GenGroomedJet_prsubjet1_px;   //!
		TBranch        *b_GenGroomedJet_prsubjet1_py;   //!
		TBranch        *b_GenGroomedJet_prsubjet1_pz;   //!
		TBranch        *b_GenGroomedJet_prsubjet1_e;   //!
		TBranch        *b_GenGroomedJet_prsubjet2_px;   //!
		TBranch        *b_GenGroomedJet_prsubjet2_py;   //!
		TBranch        *b_GenGroomedJet_prsubjet2_pz;   //!
		TBranch        *b_GenGroomedJet_prsubjet2_e;   //!
		TBranch        *b_GenGroomedJet_mass;   //!
		TBranch        *b_GenGroomedJet_mass_rhoArea;   //!
		TBranch        *b_GenGroomedJet_mass_rhoGArea;   //!
		TBranch        *b_GenGroomedJet_mass_rho4A;   //!
		TBranch        *b_GenGroomedJet_mass_rhoG4A;   //!
		TBranch        *b_GenGroomedJet_mass_rhom4Am;   //!
		TBranch        *b_GenGroomedJet_mass_shapesubtraction;   //!
		TBranch        *b_GenGroomedJet_mass_trimmingshapesubtraction;   //!
		TBranch        *b_GenGroomedJet_mass_tr;   //!
		TBranch        *b_GenGroomedJet_mass_ft;   //!
		TBranch        *b_GenGroomedJet_mass_pr;   //!
		TBranch        *b_GenGroomedJet_massdrop_pr;   //!
		TBranch        *b_GenGroomedJet_area;   //!
		TBranch        *b_GenGroomedJet_area_tr;   //!
		TBranch        *b_GenGroomedJet_area_ft;   //!
		TBranch        *b_GenGroomedJet_area_pr;   //!
		TBranch        *b_GenGroomedJet_jetconstituents;   //!
		TBranch        *b_GenGroomedJet_jetcharge;   //!
		TBranch        *b_GenGroomedJet_rcores;   //!
		TBranch        *b_GenGroomedJet_ptcores;   //!
		TBranch        *b_GenGroomedJet_planarflow;   //!
		TBranch        *b_GenGroomedJet_qjetmass;   //!
		TBranch        *b_GenGroomedJet_qjetmassdrop;   //!
		TBranch        *b_GenGroomedJet_mass_JetCleansing_DiffMode;   //!
		TBranch        *b_GenGroomedJet_pt_JetCleansing_DiffMode;   //!
		TBranch        *b_GenGroomedJet_eta_JetCleansing_DiffMode;   //!
		TBranch        *b_GenGroomedJet_phi_JetCleansing_DiffMode;   //!
		TBranch        *b_GenGroomedJet_tau2tau1_JetCleansing_DiffMode;   //!
		TBranch        *b_GenGroomedJet_constituents0_eta;   //!
		TBranch        *b_GenGroomedJet_constituents0_phi;   //!
		TBranch        *b_GenGroomedJet_constituents0_e;   //!
		TBranch        *b_GenGroomedJet_nconstituents0;   //!
		TBranch        *b_GenGroomedJet_constituents0pr_eta;   //!
		TBranch        *b_GenGroomedJet_constituents0pr_phi;   //!
		TBranch        *b_GenGroomedJet_constituents0pr_e;   //!
		TBranch        *b_GenGroomedJet_nconstituents0pr;   //!
		TBranch        *b_GenGroomedJet_rhoSW;   //!
		TBranch        *b_GenGroomedJet_rhohand;   //!
		TBranch        *b_GenGroomedJet_rhohand2;   //!
		TBranch        *b_GenGroomedJet_rhogrid;   //!

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
		TBranch        *b_EventWeight;   //!
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

		//MyClass(TTree *tree, char* inFinalState, char* inJetType, char* inPfType, Bool_t in_isBoosted,char* plot_dir); //JetType: AK5 AK8, PF PFCHS 
		MyClass(TTree *tree, char* inFinalState, char* inJetType, char* inPfType, Int_t in_isBoosted,char* plot_dir); //JetType: AK5 AK8, PF PFCHS 
		virtual ~MyClass();
		//virtual Int_t    Cut(Long64_t entry);
		virtual Bool_t   preSelect();
		virtual Int_t    GetEntry(Long64_t entry);
		virtual Long64_t LoadTree(Long64_t entry);
		virtual void     Init(TTree *tree);
		virtual void     Loop();
		virtual Bool_t   Notify();
		virtual void     Show(Long64_t entry = -1);

		void     Draw_and_Save(TH1D, char* addtional_info="");
		void     Draw_and_Save(TH1D, TH1D);
		void     Draw_and_Save(TH1D, TH1D, TH1D);
		void     Draw_and_Save(TH1D, TH1D, TH1D, TH1D);
		void     Draw_and_Save(TH2D, char* addtional_info="");

		void     DrawPlots(vector< TString >);

		void     Draw_and_Print_All();

		void     SaveResponse(JetCorrectionTool& jct, TString xdenominator_var_name,TString xnumerator_var_name , Int_t nbin=100, Double_t xmin=0, Double_t xmax=2. );

};

bool equal(double a1, double a2, double delta=1e-3);
bool match_dR(double a1, double a2, double b1, double b2, double delta=0.3, TH1D* h1=NULL);

#endif

#ifdef MyClass_cxx
//MyClass::MyClass(TTree *tree, char* inFinalState, char* inJetType, char* inPfType, Bool_t in_isBoosted, char* plot_dir) : fChain(0) 
MyClass::MyClass(TTree *tree, char* inFinalState, char* inJetType, char* inPfType, Int_t in_isBoosted, char* plot_dir) : fChain(0) 
{
	plot_Dir_DateTime=plot_dir;
	//fout = new TFile(Form("%s/out_%s.root", plot_Dir_DateTime , boostedlable.Data()), "RECREATE");
	fout = new TFile(Form("%s/out.root", plot_Dir_DateTime.Data() ), "RECREATE");

	// if parameter tree is not specified (or zero), connect the file
	// used to generate this class and read the Tree.

	if (tree == 0) {
		TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("zmumujetsanalysisntuple.root");
		if (!f || !f->IsOpen()) {
			f = new TFile("zmumujetsanalysisntuple.root");
		}
		f->GetObject("ZJet",tree);

	}
	FinalState=inFinalState;
	JetType=inJetType;
	PfType =inPfType;
	isBoosted=in_isBoosted;
	Init(tree);
	cout<<"FinalState="<<FinalState<<endl;

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



	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_pt_uncorr", JetType.Data(), PfType.Data()), GroomedJet_pt_uncorr, &b_GroomedJet_pt_uncorr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_number_jet_central", JetType.Data(), PfType.Data()), &GroomedJet_number_jet_central, &b_GroomedJet_number_jet_central);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_mass_uncorr", JetType.Data(), PfType.Data()), GroomedJet_mass_uncorr, &b_GroomedJet_mass_uncorr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_mass_tr_uncorr", JetType.Data(), PfType.Data()), GroomedJet_mass_tr_uncorr, &b_GroomedJet_mass_tr_uncorr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_mass_ft_uncorr", JetType.Data(), PfType.Data()), GroomedJet_mass_ft_uncorr, &b_GroomedJet_mass_ft_uncorr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_mass_pr_uncorr", JetType.Data(), PfType.Data()), GroomedJet_mass_pr_uncorr, &b_GroomedJet_mass_pr_uncorr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_tau2tau1", JetType.Data(), PfType.Data()), GroomedJet_tau2tau1, &b_GroomedJet_tau2tau1);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_tau2tau1_shapesubtraction", JetType.Data(), PfType.Data()), GroomedJet_tau2tau1_shapesubtraction, &b_GroomedJet_tau2tau1_shapesubtraction);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_tau1", JetType.Data(), PfType.Data()), GroomedJet_tau1, &b_GroomedJet_tau1);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_tau2", JetType.Data(), PfType.Data()), GroomedJet_tau2, &b_GroomedJet_tau2);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_tau3", JetType.Data(), PfType.Data()), GroomedJet_tau3, &b_GroomedJet_tau3);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_tau4", JetType.Data(), PfType.Data()), GroomedJet_tau4, &b_GroomedJet_tau4);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_massdrop_pr_uncorr", JetType.Data(), PfType.Data()), GroomedJet_massdrop_pr_uncorr, &b_GroomedJet_massdrop_pr_uncorr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_pt", JetType.Data(), PfType.Data()), GroomedJet_pt, &b_GroomedJet_pt);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_eta", JetType.Data(), PfType.Data()), GroomedJet_eta, &b_GroomedJet_eta);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_phi", JetType.Data(), PfType.Data()), GroomedJet_phi, &b_GroomedJet_phi);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_e", JetType.Data(), PfType.Data()), GroomedJet_e, &b_GroomedJet_e);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_pt_JECL1", JetType.Data(), PfType.Data()), GroomedJet_pt_JECL1, &b_GroomedJet_pt_JECL1);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_eta_JECL1", JetType.Data(), PfType.Data()), GroomedJet_eta_JECL1, &b_GroomedJet_eta_JECL1);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_phi_JECL1", JetType.Data(), PfType.Data()), GroomedJet_phi_JECL1, &b_GroomedJet_phi_JECL1);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_e_JECL1", JetType.Data(), PfType.Data()), GroomedJet_e_JECL1, &b_GroomedJet_e_JECL1);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_pt_L1_rhoSW", JetType.Data(), PfType.Data()), GroomedJet_pt_L1_rhoSW, &b_GroomedJet_pt_L1_rhoSW);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_pt_L1_rhoHand", JetType.Data(), PfType.Data()), GroomedJet_pt_L1_rhoHand, &b_GroomedJet_pt_L1_rhoHand);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_pt_L1_rhoHand2", JetType.Data(), PfType.Data()), GroomedJet_pt_L1_rhoHand2, &b_GroomedJet_pt_L1_rhoHand2);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_pt_L1_rhoGrid", JetType.Data(), PfType.Data()), GroomedJet_pt_L1_rhoGrid, &b_GroomedJet_pt_L1_rhoGrid);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_pt_rho4A", JetType.Data(), PfType.Data()), GroomedJet_pt_rho4A, &b_GroomedJet_pt_rho4A);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_pt_shapesubtraction", JetType.Data(), PfType.Data()), GroomedJet_pt_shapesubtraction, &b_GroomedJet_pt_shapesubtraction);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_pt_tr_uncorr", JetType.Data(), PfType.Data()), GroomedJet_pt_tr_uncorr, &b_GroomedJet_pt_tr_uncorr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_pt_tr", JetType.Data(), PfType.Data()), GroomedJet_pt_tr, &b_GroomedJet_pt_tr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_eta_tr", JetType.Data(), PfType.Data()), GroomedJet_eta_tr, &b_GroomedJet_eta_tr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_phi_tr", JetType.Data(), PfType.Data()), GroomedJet_phi_tr, &b_GroomedJet_phi_tr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_e_tr", JetType.Data(), PfType.Data()), GroomedJet_e_tr, &b_GroomedJet_e_tr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_pt_ft_uncorr", JetType.Data(), PfType.Data()), GroomedJet_pt_ft_uncorr, &b_GroomedJet_pt_ft_uncorr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_pt_ft", JetType.Data(), PfType.Data()), GroomedJet_pt_ft, &b_GroomedJet_pt_ft);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_eta_ft", JetType.Data(), PfType.Data()), GroomedJet_eta_ft, &b_GroomedJet_eta_ft);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_phi_ft", JetType.Data(), PfType.Data()), GroomedJet_phi_ft, &b_GroomedJet_phi_ft);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_e_ft", JetType.Data(), PfType.Data()), GroomedJet_e_ft, &b_GroomedJet_e_ft);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_pt_pr_uncorr", JetType.Data(), PfType.Data()), GroomedJet_pt_pr_uncorr, &b_GroomedJet_pt_pr_uncorr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_pt_pr", JetType.Data(), PfType.Data()), GroomedJet_pt_pr, &b_GroomedJet_pt_pr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_eta_pr", JetType.Data(), PfType.Data()), GroomedJet_eta_pr, &b_GroomedJet_eta_pr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_phi_pr", JetType.Data(), PfType.Data()), GroomedJet_phi_pr, &b_GroomedJet_phi_pr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_e_pr", JetType.Data(), PfType.Data()), GroomedJet_e_pr, &b_GroomedJet_e_pr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_prsubjet1_px", JetType.Data(), PfType.Data()), GroomedJet_prsubjet1_px, &b_GroomedJet_prsubjet1_px);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_prsubjet1_py", JetType.Data(), PfType.Data()), GroomedJet_prsubjet1_py, &b_GroomedJet_prsubjet1_py);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_prsubjet1_pz", JetType.Data(), PfType.Data()), GroomedJet_prsubjet1_pz, &b_GroomedJet_prsubjet1_pz);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_prsubjet1_e", JetType.Data(), PfType.Data()), GroomedJet_prsubjet1_e, &b_GroomedJet_prsubjet1_e);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_prsubjet2_px", JetType.Data(), PfType.Data()), GroomedJet_prsubjet2_px, &b_GroomedJet_prsubjet2_px);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_prsubjet2_py", JetType.Data(), PfType.Data()), GroomedJet_prsubjet2_py, &b_GroomedJet_prsubjet2_py);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_prsubjet2_pz", JetType.Data(), PfType.Data()), GroomedJet_prsubjet2_pz, &b_GroomedJet_prsubjet2_pz);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_prsubjet2_e", JetType.Data(), PfType.Data()), GroomedJet_prsubjet2_e, &b_GroomedJet_prsubjet2_e);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_mass", JetType.Data(), PfType.Data()), GroomedJet_mass, &b_GroomedJet_mass);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_mass_rhoArea", JetType.Data(), PfType.Data()), GroomedJet_mass_rhoArea, &b_GroomedJet_mass_rhoArea);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_mass_rhoGArea", JetType.Data(), PfType.Data()), GroomedJet_mass_rhoGArea, &b_GroomedJet_mass_rhoGArea);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_mass_rho4A", JetType.Data(), PfType.Data()), GroomedJet_mass_rho4A, &b_GroomedJet_mass_rho4A);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_mass_rhoG4A", JetType.Data(), PfType.Data()), GroomedJet_mass_rhoG4A, &b_GroomedJet_mass_rhoG4A);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_mass_rhom4Am", JetType.Data(), PfType.Data()), GroomedJet_mass_rhom4Am, &b_GroomedJet_mass_rhom4Am);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_mass_shapesubtraction", JetType.Data(), PfType.Data()), GroomedJet_mass_shapesubtraction, &b_GroomedJet_mass_shapesubtraction);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_mass_trimmingshapesubtraction", JetType.Data(), PfType.Data()), GroomedJet_mass_trimmingshapesubtraction, &b_GroomedJet_mass_trimmingshapesubtraction);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_mass_tr", JetType.Data(), PfType.Data()), GroomedJet_mass_tr, &b_GroomedJet_mass_tr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_mass_ft", JetType.Data(), PfType.Data()), GroomedJet_mass_ft, &b_GroomedJet_mass_ft);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_mass_pr", JetType.Data(), PfType.Data()), GroomedJet_mass_pr, &b_GroomedJet_mass_pr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_massdrop_pr", JetType.Data(), PfType.Data()), GroomedJet_massdrop_pr, &b_GroomedJet_massdrop_pr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_area", JetType.Data(), PfType.Data()), GroomedJet_area, &b_GroomedJet_area);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_area_tr", JetType.Data(), PfType.Data()), GroomedJet_area_tr, &b_GroomedJet_area_tr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_area_ft", JetType.Data(), PfType.Data()), GroomedJet_area_ft, &b_GroomedJet_area_ft);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_area_pr", JetType.Data(), PfType.Data()), GroomedJet_area_pr, &b_GroomedJet_area_pr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_jetconstituents", JetType.Data(), PfType.Data()), GroomedJet_jetconstituents, &b_GroomedJet_jetconstituents);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_jetcharge", JetType.Data(), PfType.Data()), GroomedJet_jetcharge, &b_GroomedJet_jetcharge);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_rcores", JetType.Data(), PfType.Data()), GroomedJet_rcores, &b_GroomedJet_rcores);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_ptcores", JetType.Data(), PfType.Data()), GroomedJet_ptcores, &b_GroomedJet_ptcores);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_planarflow", JetType.Data(), PfType.Data()), GroomedJet_planarflow, &b_GroomedJet_planarflow);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_qjetmass", JetType.Data(), PfType.Data()), GroomedJet_qjetmass, &b_GroomedJet_qjetmass);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_qjetmassdrop", JetType.Data(), PfType.Data()), GroomedJet_qjetmassdrop, &b_GroomedJet_qjetmassdrop);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_mass_JetCleansing_DiffMode", JetType.Data(), PfType.Data()), GroomedJet_mass_JetCleansing_DiffMode, &b_GroomedJet_mass_JetCleansing_DiffMode);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_pt_JetCleansing_DiffMode", JetType.Data(), PfType.Data()), GroomedJet_pt_JetCleansing_DiffMode, &b_GroomedJet_pt_JetCleansing_DiffMode);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_eta_JetCleansing_DiffMode", JetType.Data(), PfType.Data()), GroomedJet_eta_JetCleansing_DiffMode, &b_GroomedJet_eta_JetCleansing_DiffMode);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_phi_JetCleansing_DiffMode", JetType.Data(), PfType.Data()), GroomedJet_phi_JetCleansing_DiffMode, &b_GroomedJet_phi_JetCleansing_DiffMode);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_tau2tau1_JetCleansing_DiffMode", JetType.Data(), PfType.Data()), GroomedJet_tau2tau1_JetCleansing_DiffMode, &b_GroomedJet_tau2tau1_JetCleansing_DiffMode);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_constituents0_eta", JetType.Data(), PfType.Data()), GroomedJet_constituents0_eta, &b_GroomedJet_constituents0_eta);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_constituents0_phi", JetType.Data(), PfType.Data()), GroomedJet_constituents0_phi, &b_GroomedJet_constituents0_phi);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_constituents0_e", JetType.Data(), PfType.Data()), GroomedJet_constituents0_e, &b_GroomedJet_constituents0_e);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_nconstituents0", JetType.Data(), PfType.Data()), &GroomedJet_nconstituents0, &b_GroomedJet_nconstituents0);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_constituents0pr_eta", JetType.Data(), PfType.Data()), GroomedJet_constituents0pr_eta, &b_GroomedJet_constituents0pr_eta);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_constituents0pr_phi", JetType.Data(), PfType.Data()), GroomedJet_constituents0pr_phi, &b_GroomedJet_constituents0pr_phi);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_constituents0pr_e", JetType.Data(), PfType.Data()), GroomedJet_constituents0pr_e, &b_GroomedJet_constituents0pr_e);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_nconstituents0pr", JetType.Data(), PfType.Data()), &GroomedJet_nconstituents0pr, &b_GroomedJet_nconstituents0pr);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_rhoSW", JetType.Data(), PfType.Data()), &GroomedJet_rhoSW, &b_GroomedJet_rhoSW);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_rhohand", JetType.Data(), PfType.Data()), &GroomedJet_rhohand, &b_GroomedJet_rhohand);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_rhohand2", JetType.Data(), PfType.Data()), &GroomedJet_rhohand2, &b_GroomedJet_rhohand2);
	fChain->SetBranchAddress(Form("GroomedJet_%s_%s_rhogrid", JetType.Data(), PfType.Data()), &GroomedJet_rhogrid, &b_GroomedJet_rhogrid);

	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_pt_uncorr", JetType.Data() ), GenGroomedJet_pt_uncorr, &b_GenGroomedJet_pt_uncorr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_number_jet_central", JetType.Data() ), &GenGroomedJet_number_jet_central, &b_GenGroomedJet_number_jet_central);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_mass_uncorr", JetType.Data() ), GenGroomedJet_mass_uncorr, &b_GenGroomedJet_mass_uncorr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_mass_tr_uncorr", JetType.Data() ), GenGroomedJet_mass_tr_uncorr, &b_GenGroomedJet_mass_tr_uncorr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_mass_ft_uncorr", JetType.Data() ), GenGroomedJet_mass_ft_uncorr, &b_GenGroomedJet_mass_ft_uncorr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_mass_pr_uncorr", JetType.Data() ), GenGroomedJet_mass_pr_uncorr, &b_GenGroomedJet_mass_pr_uncorr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_tau2tau1", JetType.Data() ), GenGroomedJet_tau2tau1, &b_GenGroomedJet_tau2tau1);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_tau2tau1_shapesubtraction", JetType.Data() ), GenGroomedJet_tau2tau1_shapesubtraction, &b_GenGroomedJet_tau2tau1_shapesubtraction);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_tau1", JetType.Data() ), GenGroomedJet_tau1, &b_GenGroomedJet_tau1);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_tau2", JetType.Data() ), GenGroomedJet_tau2, &b_GenGroomedJet_tau2);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_tau3", JetType.Data() ), GenGroomedJet_tau3, &b_GenGroomedJet_tau3);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_tau4", JetType.Data() ), GenGroomedJet_tau4, &b_GenGroomedJet_tau4);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_massdrop_pr_uncorr", JetType.Data() ), GenGroomedJet_massdrop_pr_uncorr, &b_GenGroomedJet_massdrop_pr_uncorr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_pt", JetType.Data() ), GenGroomedJet_pt, &b_GenGroomedJet_pt);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_eta", JetType.Data() ), GenGroomedJet_eta, &b_GenGroomedJet_eta);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_phi", JetType.Data() ), GenGroomedJet_phi, &b_GenGroomedJet_phi);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_e", JetType.Data() ), GenGroomedJet_e, &b_GenGroomedJet_e);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_pt_L1_rhoSW", JetType.Data() ), GenGroomedJet_pt_L1_rhoSW, &b_GenGroomedJet_pt_L1_rhoSW);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_pt_L1_rhoHand", JetType.Data() ), GenGroomedJet_pt_L1_rhoHand, &b_GenGroomedJet_pt_L1_rhoHand);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_pt_L1_rhoHand2", JetType.Data() ), GenGroomedJet_pt_L1_rhoHand2, &b_GenGroomedJet_pt_L1_rhoHand2);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_pt_L1_rhoGrid", JetType.Data() ), GenGroomedJet_pt_L1_rhoGrid, &b_GenGroomedJet_pt_L1_rhoGrid);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_pt_rho4A", JetType.Data() ), GenGroomedJet_pt_rho4A, &b_GenGroomedJet_pt_rho4A);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_pt_shapesubtraction", JetType.Data() ), GenGroomedJet_pt_shapesubtraction, &b_GenGroomedJet_pt_shapesubtraction);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_pt_tr_uncorr", JetType.Data() ), GenGroomedJet_pt_tr_uncorr, &b_GenGroomedJet_pt_tr_uncorr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_pt_tr", JetType.Data() ), GenGroomedJet_pt_tr, &b_GenGroomedJet_pt_tr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_eta_tr", JetType.Data() ), GenGroomedJet_eta_tr, &b_GenGroomedJet_eta_tr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_phi_tr", JetType.Data() ), GenGroomedJet_phi_tr, &b_GenGroomedJet_phi_tr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_e_tr", JetType.Data() ), GenGroomedJet_e_tr, &b_GenGroomedJet_e_tr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_pt_ft_uncorr", JetType.Data() ), GenGroomedJet_pt_ft_uncorr, &b_GenGroomedJet_pt_ft_uncorr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_pt_ft", JetType.Data() ), GenGroomedJet_pt_ft, &b_GenGroomedJet_pt_ft);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_eta_ft", JetType.Data() ), GenGroomedJet_eta_ft, &b_GenGroomedJet_eta_ft);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_phi_ft", JetType.Data() ), GenGroomedJet_phi_ft, &b_GenGroomedJet_phi_ft);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_e_ft", JetType.Data() ), GenGroomedJet_e_ft, &b_GenGroomedJet_e_ft);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_pt_pr_uncorr", JetType.Data() ), GenGroomedJet_pt_pr_uncorr, &b_GenGroomedJet_pt_pr_uncorr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_pt_pr", JetType.Data() ), GenGroomedJet_pt_pr, &b_GenGroomedJet_pt_pr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_eta_pr", JetType.Data() ), GenGroomedJet_eta_pr, &b_GenGroomedJet_eta_pr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_phi_pr", JetType.Data() ), GenGroomedJet_phi_pr, &b_GenGroomedJet_phi_pr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_e_pr", JetType.Data() ), GenGroomedJet_e_pr, &b_GenGroomedJet_e_pr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_prsubjet1_px", JetType.Data() ), GenGroomedJet_prsubjet1_px, &b_GenGroomedJet_prsubjet1_px);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_prsubjet1_py", JetType.Data() ), GenGroomedJet_prsubjet1_py, &b_GenGroomedJet_prsubjet1_py);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_prsubjet1_pz", JetType.Data() ), GenGroomedJet_prsubjet1_pz, &b_GenGroomedJet_prsubjet1_pz);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_prsubjet1_e", JetType.Data() ), GenGroomedJet_prsubjet1_e, &b_GenGroomedJet_prsubjet1_e);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_prsubjet2_px", JetType.Data() ), GenGroomedJet_prsubjet2_px, &b_GenGroomedJet_prsubjet2_px);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_prsubjet2_py", JetType.Data() ), GenGroomedJet_prsubjet2_py, &b_GenGroomedJet_prsubjet2_py);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_prsubjet2_pz", JetType.Data() ), GenGroomedJet_prsubjet2_pz, &b_GenGroomedJet_prsubjet2_pz);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_prsubjet2_e", JetType.Data() ), GenGroomedJet_prsubjet2_e, &b_GenGroomedJet_prsubjet2_e);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_mass", JetType.Data() ), GenGroomedJet_mass, &b_GenGroomedJet_mass);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_mass_rhoArea", JetType.Data() ), GenGroomedJet_mass_rhoArea, &b_GenGroomedJet_mass_rhoArea);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_mass_rhoGArea", JetType.Data() ), GenGroomedJet_mass_rhoGArea, &b_GenGroomedJet_mass_rhoGArea);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_mass_rho4A", JetType.Data() ), GenGroomedJet_mass_rho4A, &b_GenGroomedJet_mass_rho4A);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_mass_rhoG4A", JetType.Data() ), GenGroomedJet_mass_rhoG4A, &b_GenGroomedJet_mass_rhoG4A);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_mass_rhom4Am", JetType.Data() ), GenGroomedJet_mass_rhom4Am, &b_GenGroomedJet_mass_rhom4Am);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_mass_shapesubtraction", JetType.Data() ), GenGroomedJet_mass_shapesubtraction, &b_GenGroomedJet_mass_shapesubtraction);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_mass_trimmingshapesubtraction", JetType.Data() ), GenGroomedJet_mass_trimmingshapesubtraction, &b_GenGroomedJet_mass_trimmingshapesubtraction);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_mass_tr", JetType.Data() ), GenGroomedJet_mass_tr, &b_GenGroomedJet_mass_tr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_mass_ft", JetType.Data() ), GenGroomedJet_mass_ft, &b_GenGroomedJet_mass_ft);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_mass_pr", JetType.Data() ), GenGroomedJet_mass_pr, &b_GenGroomedJet_mass_pr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_massdrop_pr", JetType.Data() ), GenGroomedJet_massdrop_pr, &b_GenGroomedJet_massdrop_pr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_area", JetType.Data() ), GenGroomedJet_area, &b_GenGroomedJet_area);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_area_tr", JetType.Data() ), GenGroomedJet_area_tr, &b_GenGroomedJet_area_tr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_area_ft", JetType.Data() ), GenGroomedJet_area_ft, &b_GenGroomedJet_area_ft);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_area_pr", JetType.Data() ), GenGroomedJet_area_pr, &b_GenGroomedJet_area_pr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_jetconstituents", JetType.Data() ), GenGroomedJet_jetconstituents, &b_GenGroomedJet_jetconstituents);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_jetcharge", JetType.Data() ), GenGroomedJet_jetcharge, &b_GenGroomedJet_jetcharge);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_rcores", JetType.Data() ), GenGroomedJet_rcores, &b_GenGroomedJet_rcores);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_ptcores", JetType.Data() ), GenGroomedJet_ptcores, &b_GenGroomedJet_ptcores);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_planarflow", JetType.Data() ), GenGroomedJet_planarflow, &b_GenGroomedJet_planarflow);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_qjetmass", JetType.Data() ), GenGroomedJet_qjetmass, &b_GenGroomedJet_qjetmass);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_qjetmassdrop", JetType.Data() ), GenGroomedJet_qjetmassdrop, &b_GenGroomedJet_qjetmassdrop);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_mass_JetCleansing_DiffMode", JetType.Data() ), GenGroomedJet_mass_JetCleansing_DiffMode, &b_GenGroomedJet_mass_JetCleansing_DiffMode);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_pt_JetCleansing_DiffMode", JetType.Data() ), GenGroomedJet_pt_JetCleansing_DiffMode, &b_GenGroomedJet_pt_JetCleansing_DiffMode);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_eta_JetCleansing_DiffMode", JetType.Data() ), GenGroomedJet_eta_JetCleansing_DiffMode, &b_GenGroomedJet_eta_JetCleansing_DiffMode);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_phi_JetCleansing_DiffMode", JetType.Data() ), GenGroomedJet_phi_JetCleansing_DiffMode, &b_GenGroomedJet_phi_JetCleansing_DiffMode);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_tau2tau1_JetCleansing_DiffMode", JetType.Data() ), GenGroomedJet_tau2tau1_JetCleansing_DiffMode, &b_GenGroomedJet_tau2tau1_JetCleansing_DiffMode);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_constituents0_eta", JetType.Data() ), GenGroomedJet_constituents0_eta, &b_GenGroomedJet_constituents0_eta);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_constituents0_phi", JetType.Data() ), GenGroomedJet_constituents0_phi, &b_GenGroomedJet_constituents0_phi);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_constituents0_e", JetType.Data() ), GenGroomedJet_constituents0_e, &b_GenGroomedJet_constituents0_e);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_nconstituents0", JetType.Data() ), &GenGroomedJet_nconstituents0, &b_GenGroomedJet_nconstituents0);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_constituents0pr_eta", JetType.Data() ), GenGroomedJet_constituents0pr_eta, &b_GenGroomedJet_constituents0pr_eta);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_constituents0pr_phi", JetType.Data() ), GenGroomedJet_constituents0pr_phi, &b_GenGroomedJet_constituents0pr_phi);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_constituents0pr_e", JetType.Data() ), GenGroomedJet_constituents0pr_e, &b_GenGroomedJet_constituents0pr_e);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_nconstituents0pr", JetType.Data() ), &GenGroomedJet_nconstituents0pr, &b_GenGroomedJet_nconstituents0pr);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_rhoSW", JetType.Data() ), &GenGroomedJet_rhoSW, &b_GenGroomedJet_rhoSW);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_rhohand", JetType.Data() ), &GenGroomedJet_rhohand, &b_GenGroomedJet_rhohand);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_rhohand2", JetType.Data() ), &GenGroomedJet_rhohand2, &b_GenGroomedJet_rhohand2);
	fChain->SetBranchAddress(Form("GenGroomedJet_%s_GEN_rhogrid", JetType.Data() ), &GenGroomedJet_rhogrid, &b_GenGroomedJet_rhogrid);



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
	fChain->SetBranchAddress("EventWeight", &EventWeight, &b_EventWeight);
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
/*Int_t MyClass::Cut(Long64_t entry)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}*/


Efficiency_Tool::Efficiency_Tool(){
	number_step=0;
	map_stepnumber.clear();
	map_stepname.clear();
	vect_step_event.clear();
}

void Efficiency_Tool::Add_Step(string newstepname){
	map_stepnumber.insert( map< Int_t, string >::value_type( number_step, newstepname) );
	map_stepname.insert( map< string, Int_t >::value_type( newstepname, number_step) );
	vect_step_event.push_back(0.);
	number_step++;
}
void Efficiency_Tool::Add_Event(string newstepname, Double_t eventweight){
	map<string, Int_t >::iterator iter;
	iter= map_stepname.find(newstepname);
	if(iter == map_stepname.end()){ 
		//cout<<"can't find this step"<<endl;return ;
		Add_Step(newstepname);
		iter= map_stepname.find(newstepname);
	}
	vect_step_event[iter->second]=vect_step_event[iter->second]+eventweight;
}
TH1D Efficiency_Tool::Get_Eff_hist(){
	TH1D hist_eff("hist_eff","hist_eff;;Event Number",number_step, 0, number_step);
	hist_eff.SetStats(0);
	hist_eff.SetBit(TH1::kCanRebin);
	map<Int_t, string>::iterator iter;
	for(Int_t i=0;i<number_step;i++){
		hist_eff.SetBinContent(i+1, vect_step_event[i]);
		iter=map_stepnumber.find(i);
		hist_eff.GetXaxis()->SetBinLabel(i+1, (iter->second).c_str() );
	}
	hist_eff.GetYaxis()->SetRangeUser( 0., 1.1*hist_eff.GetBinContent(1) );
	return hist_eff;
}


RealVarArray::RealVarArray(const char* in_name, Int_t in_nbin, Double_t in_xmin, Double_t in_xmax){
	name=TString(in_name);
	nbin=in_nbin;
	xmin=in_xmin;
	xmax=in_xmax;
}

RealVarArray::RealVarArray() { name=TString("tmp"); nbin=1; xmin=0; xmax=0; }



JetCorrectionTool::JetCorrectionTool(const char* in_name){ name=TString(in_name); }

void JetCorrectionTool::addVar( RealVarArray x) {
	map_obs.insert( MAP_REALVARARRAY::value_type(x.name, x) );
}

Bool_t JetCorrectionTool::fill( TString var_name, Double_t in_val, Double_t eventweight, Bool_t debug){
	//MAP_REALVARARRAY::iterator it=map_obs.find(var_name);
	if ( !find_var_in_map(map_obs, var_name) ) { cout<<"When fill, fail to find "<<var_name.Data()<<endl; return 0;}
	else {
		if(debug) cout<<var_name.Data()<<" = "<<in_val<<endl;
		map_obs[var_name].vect_x.push_back(in_val);
		map_obs[var_name].vect_weight.push_back(eventweight);
		return 1;
	}
}
vector< TString > JetCorrectionTool::get_allvar_name(){
	vector< TString > vect_allvar_name;
	MAP_REALVARARRAY::iterator iter;
	for( iter=map_obs.begin(); iter != map_obs.end(); iter++){
		vect_allvar_name.push_back( iter->first );
	}
	return vect_allvar_name;
}


TH1D JetCorrectionTool::get_hist1D(TString var_name){
	if ( !find_var_in_map(map_obs, var_name) ) { cout<<"When get_hist1D, fail to find "<<var_name.Data()<<endl;TH1D("h1","h1",1,0,1);}

	TH1D h1(Form("h1_JCT_%s_%s", name.Data(),var_name.Data()), Form("h1_JCT_%s_%s;%s;Event Number", name.Data(),var_name.Data(), var_name.Data()), map_obs[var_name].nbin, map_obs[var_name].xmin, map_obs[var_name].xmax );
	for(Int_t k=0;k<Int_t(map_obs[var_name].vect_x.size());k++){
		h1.Fill( map_obs[var_name].vect_x[k], map_obs[var_name].vect_weight[k] ); 
	}
	return h1;
}

Double_t RMS_weight(Int_t n, Double_t* array, Double_t* array_weight){

		Double_t mean=TMath::Mean(n, array, array_weight);
		Double_t rms=0.; //TMath::RMS(tmp_size, array_x);
		Double_t totalweight=0.;
		for(Int_t i=0;i<n;i++){
			rms+=(array[i]-mean)*(array[i]-mean)*array_weight[i]*array_weight[i];
			totalweight+=array_weight[i]*array_weight[i];
		}
		rms=TMath::Sqrt(rms/totalweight);
		return rms;
};

TH1D JetCorrectionTool::get_hist1D_response( TString xdenominator_var_name,TString xnumerator_var_name , Int_t nbin, Double_t xmin, Double_t xmax ){
	if ( !find_var_in_map(map_obs, xdenominator_var_name) ) { cout<<"When get_hist1D_response, fail to find "<<xdenominator_var_name.Data()<<endl;return TH1D("h1","h1",1,0,1);}
	if ( !find_var_in_map(map_obs, xnumerator_var_name) ) { cout<<"When get_hist1D_response, fail to find "<<xnumerator_var_name.Data()<<endl;return TH1D("h1","h1",1,0,1);}

	TH1D h1(Form("h1response_JCT_%s_%s_%s", name.Data(), xnumerator_var_name.Data(), xdenominator_var_name.Data()), Form("h1response_JCT_%s_%s_%s;%s/%s;Event Number", name.Data(),xnumerator_var_name.Data(), xdenominator_var_name.Data(), xnumerator_var_name.Data(), xdenominator_var_name.Data()), nbin, xmin, xmax);
	if( map_obs[xdenominator_var_name].vect_x.size() == map_obs[xnumerator_var_name].vect_x.size() ){


		Int_t tmp_size= Int_t(map_obs[xdenominator_var_name].vect_x.size());
		Double_t * array_x=new Double_t[tmp_size];
		Double_t * array_x_weight=new Double_t[tmp_size];

		for(Int_t k=0;k<tmp_size;k++){
			Double_t tmp = map_obs[xnumerator_var_name].vect_x[k]/ map_obs[xdenominator_var_name].vect_x[k]; 
			Double_t tmp_weight = map_obs[xnumerator_var_name].vect_weight[k] ;
			array_x[k]=tmp;
			array_x_weight[k]=tmp_weight;
			h1.Fill( tmp, tmp_weight); 
		}
		//Double_t mean=TMath::Mean(tmp_size, array_x);
		//Double_t rms=TMath::RMS(tmp_size, array_x);

		//Double_t mean_weight=TMath::Mean(tmp_size, array_x, array_x_weight);
		//Double_t rms_weight=RMS_weight(tmp_size, array_x, array_x_weight);

		//cout<<"======= Draw TH1-response: ========"<<h1.GetTitle()<<endl;
		//cout<<"     TH1 mean, rms, relative-rms= "<<h1.GetMean()<<", "<<h1.GetRMS()<<", "<<h1.GetRMS()/h1.GetMean()<<endl;
		//cout<<"         mean, rms, relative-rms= "<<mean<<", "<<rms<<", "<<rms/mean<<endl;
		//cout<<"Weighted mean, rms, relative-rms= "<<mean_weight<<", "<<rms_weight<<", "<<rms_weight/mean_weight<<endl;

		return h1;
	}else { cout<<"When get_hist1D_response, Two var didn't match: "<<xnumerator_var_name.Data()<<" , "<<xdenominator_var_name.Data()<<endl;return TH1D("h1","h1",1,0,1);}

}

RESPONSE JetCorrectionTool::get_response( TString xdenominator_var_name,TString xnumerator_var_name , Int_t nbin, Double_t xmin, Double_t xmax ){
	RESPONSE response0;
	if ( !find_var_in_map(map_obs, xdenominator_var_name) ) { cout<<"When get_response, fail to find "<<xdenominator_var_name.Data()<<endl;return response0;}
	if ( !find_var_in_map(map_obs, xnumerator_var_name) ) { cout<<"When get_response, fail to find "<<xnumerator_var_name.Data()<<endl;return response0;}

	if( map_obs[xdenominator_var_name].vect_x.size() == map_obs[xnumerator_var_name].vect_x.size() ){

		if( find_var_in_map( map_response, Form("h1response_JCT_%s_%s_%s", name.Data(), xnumerator_var_name.Data(), xdenominator_var_name.Data()) ) ){ 
			return map_response[Form("h1response_JCT_%s_%s_%s", name.Data(), xnumerator_var_name.Data(), xdenominator_var_name.Data())];
		}else{
			//1D plots
			TH1D h1(Form("h1response_JCT_%s_%s_%s", name.Data(), xnumerator_var_name.Data(), xdenominator_var_name.Data()), Form("h1response_JCT_%s_%s_%s;%s/%s;Event Number", name.Data(),xnumerator_var_name.Data(), xdenominator_var_name.Data(), xnumerator_var_name.Data(), xdenominator_var_name.Data()), nbin, xmin, xmax);

			Int_t tmp_size= Int_t(map_obs[xdenominator_var_name].vect_x.size());
			Double_t * array_x=new Double_t[tmp_size];
			Double_t * array_x_weight=new Double_t[tmp_size];

			for(Int_t k=0;k<tmp_size;k++){
				Double_t tmp = map_obs[xnumerator_var_name].vect_x[k]/ map_obs[xdenominator_var_name].vect_x[k]; 
				Double_t tmp_weight = map_obs[xnumerator_var_name].vect_weight[k] ;
				array_x[k]=tmp;
				array_x_weight[k]=tmp_weight;
				h1.Fill( tmp, tmp_weight); 
			}
			Double_t mean=TMath::Mean(tmp_size, array_x);
			Double_t rms=TMath::RMS(tmp_size, array_x);

			Double_t mean_weight=TMath::Mean(tmp_size, array_x, array_x_weight);
			Double_t rms_weight=RMS_weight(tmp_size, array_x, array_x_weight);

			cout<<"======= Draw TH1-response: "<<h1.GetTitle()<<endl;
			cout<<"     TH1 mean, rms, relative-rms= "<<h1.GetMean()<<", "<<h1.GetRMS()<<", "<<h1.GetRMS()/h1.GetMean()<<endl;
			cout<<"         mean, rms, relative-rms= "<<mean<<", "<<rms<<", "<<rms/mean<<endl;
			cout<<"Weighted mean, rms, relative-rms= "<<mean_weight<<", "<<rms_weight<<", "<<rms_weight/mean_weight<<endl;
			response0.hist=h1;
			response0.mean=mean_weight;
			response0.rms=rms_weight;

			//2D plots
			TH2D h2(Form("h2_JCT_%s_%s_%s", name.Data(),xnumerator_var_name.Data(),xdenominator_var_name.Data()), Form("h2_JCT_%s_%s_%s;%s;%s;", name.Data(),xnumerator_var_name.Data(),xdenominator_var_name.Data(), xdenominator_var_name.Data(), xnumerator_var_name.Data()), map_obs[xdenominator_var_name].nbin, map_obs[xdenominator_var_name].xmin, map_obs[xdenominator_var_name].xmax, map_obs[xnumerator_var_name].nbin, map_obs[xnumerator_var_name].xmin, map_obs[xnumerator_var_name].xmax  );
			TGraph gr;

			for(Int_t k=0;k<Int_t(map_obs[xdenominator_var_name].vect_x.size());k++){
				h2.Fill( map_obs[xdenominator_var_name].vect_x[k], map_obs[xnumerator_var_name].vect_x[k], map_obs[xnumerator_var_name].vect_weight[k]); 
				gr.SetPoint( k,  map_obs[xdenominator_var_name].vect_x[k], map_obs[xnumerator_var_name].vect_x[k]); 
			}
			response0.hist2D=h2;
			response0.correlationfactor=gr.GetCorrelationFactor();
			//correlation factor in 2D plots

			//add to map
			map_response.insert( MAP_RESPONSE::value_type(h1.GetTitle(), response0) );
			return response0;
		}
	}else { cout<<"When get_response, Two var didn't match: "<<xnumerator_var_name.Data()<<" , "<<xdenominator_var_name.Data()<<endl;return response0;}

}

TH2D JetCorrectionTool::get_hist2D( TString x_var_name,TString y_var_name ){
	if ( !find_var_in_map(map_obs, x_var_name) ) { cout<<"When get_hist2D, fail to find "<<x_var_name.Data()<<endl;return TH2D("h2","h2",1,0,1,1,0,1);}
	if ( !find_var_in_map(map_obs, y_var_name) ) { cout<<"When get_hist2D, fail to find "<<y_var_name.Data()<<endl;return TH2D("h2","h2",1,0,1,1,0,1);}

	TH2D h2(Form("h2_JCT_%s_%s_%s", name.Data(),y_var_name.Data(),x_var_name.Data()), Form("h2_JCT_%s_%s_%s;%s;%s;", name.Data(),y_var_name.Data(),x_var_name.Data(), x_var_name.Data(), y_var_name.Data()), map_obs[x_var_name].nbin, map_obs[x_var_name].xmin, map_obs[x_var_name].xmax, map_obs[y_var_name].nbin, map_obs[y_var_name].xmin, map_obs[y_var_name].xmax  );
	if( map_obs[x_var_name].vect_x.size() == map_obs[y_var_name].vect_x.size() ){
		for(Int_t k=0;k<Int_t(map_obs[x_var_name].vect_x.size());k++){
			h2.Fill( map_obs[x_var_name].vect_x[k], map_obs[y_var_name].vect_x[k], map_obs[y_var_name].vect_weight[k]); 
		}
		return h2;
	}else { cout<<"When get_hist2D, Two var didn't match: "<<y_var_name.Data()<<" , "<<x_var_name.Data()<<endl;return TH2D("h2","h2",1,0,1,1,0,1);}

}

TGraph JetCorrectionTool::get_graph( TString x_var_name,TString y_var_name ){
	TGraph gr;

	if ( !find_var_in_map(map_obs, x_var_name) ) { cout<<"When get_graph, fail to find "<<x_var_name.Data()<<endl; return gr;}
	if ( !find_var_in_map(map_obs, y_var_name) ) { cout<<"When get_graph, fail to find "<<y_var_name.Data()<<endl; return gr;}

	if( map_obs[x_var_name].vect_x.size() == map_obs[y_var_name].vect_x.size() ){
		for(Int_t k=0;k<Int_t(map_obs[x_var_name].vect_x.size());k++){
			gr.SetPoint( k,  map_obs[x_var_name].vect_x[k], map_obs[y_var_name].vect_x[k]); 
		}
		return gr;
	}else { cout<<"When get_graph, Two var didn't match: "<<y_var_name.Data()<<" , "<<x_var_name.Data()<<endl; return gr;}
}

TH1D JetCorrectionTool::get_mean_rms_hist( TString x_var_name, TString ydenominator_var_name, TString ynumerator_var_name, Double_t ymin, Double_t ymax )
{
	if ( !find_var_in_map(map_obs, x_var_name) ) { cout<<"When get_mean_rms_hist, fail to find "<<x_var_name.Data()<<endl;return TH1D("h1","h1",1,0,1);}
	if ( !find_var_in_map(map_obs, ydenominator_var_name) ) { cout<<"When get_mean_rms_hist, fail to find "<<ydenominator_var_name.Data()<<endl;return TH1D("h1","h1",1,0,1);}
	if ( !find_var_in_map(map_obs, ynumerator_var_name) ) { cout<<"When get_mean_rms_hist, fail to find "<<ynumerator_var_name.Data()<<endl;return TH1D("h1","h1",1,0,1);}



	if( map_obs[x_var_name].vect_x.size() == map_obs[ydenominator_var_name].vect_x.size() && map_obs[ynumerator_var_name].vect_x.size() == map_obs[ydenominator_var_name].vect_x.size() ){

		cout<<" get_mean_rms_hist: "<<x_var_name.Data()<<", "<<ydenominator_var_name.Data()<<", "<<ynumerator_var_name.Data()<<", "<<endl;

		TH1D h1(Form("h1_JCT_%s_mean_rms_%s_%s", name.Data(),ynumerator_var_name.Data(),x_var_name.Data()), Form("h1_JCT_%s_mean_rms_%s_%s_%s;%s;%s/%s;", name.Data(),ynumerator_var_name.Data(),ydenominator_var_name.Data(),x_var_name.Data(),x_var_name.Data(), ynumerator_var_name.Data(), ydenominator_var_name.Data()), map_obs[x_var_name].nbin, map_obs[x_var_name].xmin, map_obs[x_var_name].xmax);

		std::vector< std::vector<Double_t> > vectors_y;
		for(Int_t j=0; j<map_obs[x_var_name].nbin; j++){ 
			std::vector<Double_t> tmp_vect_y; 
			vectors_y.push_back(tmp_vect_y); 
		}


		for(Int_t k=0;k<Int_t(map_obs[x_var_name].vect_x.size());k++){
			if( map_obs[x_var_name].vect_x[k] < map_obs[x_var_name].xmin || map_obs[x_var_name].vect_x[k] > map_obs[x_var_name].xmax ) continue;

			for(Int_t i=1;i<=map_obs[x_var_name].nbin;i++){
				if( map_obs[x_var_name].vect_x[k] >=h1.GetBinLowEdge(i) && map_obs[x_var_name].vect_x[k] <h1.GetBinLowEdge(i+1)  ){
					vectors_y[i-1].push_back(map_obs[ynumerator_var_name].vect_x[k] / map_obs[ydenominator_var_name].vect_x[k]);
					//if(map_obs[ynumerator_var_name].vect_x[k] / map_obs[ydenominator_var_name].vect_x[k] >3 ||map_obs[ynumerator_var_name].vect_x[k] / map_obs[ydenominator_var_name].vect_x[k] <0  ) cout<<"x= "<<h1.GetBinLowEdge(i) <<" y= "<< map_obs[ynumerator_var_name].vect_x[k]<<" / "<< map_obs[ydenominator_var_name].vect_x[k] <<" ="<< map_obs[ynumerator_var_name].vect_x[k] / map_obs[ydenominator_var_name].vect_x[k]<<endl;
					break;
				}
			}

		}

		for(Int_t j=0; j<map_obs[x_var_name].nbin; j++){ 

			Int_t n_y=vectors_y[j].size(); 
			if(n_y==0) {
				h1.SetBinContent(j+1, 0.);
				h1.SetBinError(j+1, 0.);
			}else{
				Double_t * array_y=new Double_t[n_y];
				for(Int_t k=0;k<n_y;k++)array_y[k]=vectors_y[j][k];
				Double_t mean=TMath::Mean(n_y, array_y);
				Double_t rms=TMath::RMS(n_y, array_y);
				h1.SetBinContent(j+1,mean);
				h1.SetBinError(j+1, rms);
				cout<<"x="<<h1.GetBinLowEdge(j+1)<<" y_mean="<<mean<<", y_rms="<<rms<<endl;
			}
		}

		h1.GetYaxis()->SetRangeUser(ymin, ymax);

		return h1;
	}else { cout<<"When get_mean_rms_hist, Three var didn't match: "<<ynumerator_var_name.Data()<<" , "<<ydenominator_var_name.Data()<<" , "<<x_var_name.Data()<<endl;return TH1D("h1","h1",1,0,1);}

}


Table_Tool::Table_Tool(){
	iter_x=iter_y=0;
	xaxis.clear();
	yaxis.clear();
	for(Int_t m=0; m <100; m++){
		for(Int_t n=0; n <100; n++){
			table_value[m][n]=-10000;
		}
	}
}

void Table_Tool::Insert(TString x_name, TString y_name, Double_t value){
	Int_t x_pos=-1; Int_t y_pos=-1;
	for(Int_t m=0;m<iter_x;m++){
		if(xaxis[m]==x_name)x_pos=m;
	}
	for(Int_t m=0;m<iter_y;m++){
		if(yaxis[m]==y_name)y_pos=m;
	}

	if( x_pos <0) {
		xaxis.push_back(x_name);
		x_pos=iter_x;
		iter_x++; 
	}

	if( y_pos <0) {
		yaxis.push_back(y_name);
		y_pos=iter_y;
		iter_y++; 
	}

	table_value[x_pos][y_pos]=value;

}

void Table_Tool::PrintTable(ofstream &fout)
		//void Table_Tool::PrintTable()
{
	//first line
	fout<<"Table \t";
	for(Int_t i=0; i< iter_x; i++)
	  fout<<xaxis[i]<<" \t";
	fout<<endl;

	for(Int_t j=0; j< iter_y; j++){
		fout<<yaxis[j]<<" \t";
		for(Int_t k=0; k<iter_x; k++){
			fout<<table_value[k][j]<<" \t";
		}
		fout<<endl;
	}
	fout<<endl;
}

#endif // #ifdef MyClass_cxx
