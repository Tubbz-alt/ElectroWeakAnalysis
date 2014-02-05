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

bool equal(double a1, double a2, double delta){
	if(TMath::Abs(a1-a2)<delta)return 1;
	else return 0;
}
bool match_dR(double a1, double a2, double b1, double b2, double delta, TH1D* h1){
	double dR= TMath::Sqrt( (a1-b1)*(a1-b1) + (a2-b2)*(a2-b2) ) ;
	if(h1) h1->Fill(dR);
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

void MyClass::DrawPlots(vector< TString > plotnames)
{
	TString finalname("");
	//TCanvas *c1 = new TCanvas(Form("c1_%s",h2.GetTitle()),Form("c1_%s",h2.GetTitle()),200,10,600,600);
	TCanvas *c1 = new TCanvas(Form("c1_multiplots"),Form("c1_multiplots"),200,10,600,600);
	c1->cd();

	for(Int_t i=0; i< plotnames.size(); i++){
		finalname+=plotnames[i];
		TH1* h1;
		cout<<"plotnames[i]="<<plotnames[i]<<endl;
		fout->GetObject(plotnames[i].Data(), h1);
		if (i==0) h1->Draw();
		else h1->Draw("same");
	}

	c1->Print(Form("%s/multiplots_%s_%s_%s.png",plot_Dir_DateTime.Data(), JetType.Data(), PfType.Data(),finalname.Data()));
	delete c1;
}


/*
   Bool_t MyClass::preSelect()
   {
   efftool.Add_Event("Init");

   if( FinalState.Contains("Dijets") ){

   if(isBoosted){
   if(!( GenGroomedJet_pt[0]>100 ))return 0;
   efftool.Add_Event("Gen Jet Pt>100");
//Double_t tmp_GEN_eta = GenGroomedJet_eta[0]; Double_t tmp_GEN_phi = GenGroomedJet_phi[0];
//Double_t tmpDeltaR_Vj = TMath::Sqrt( (Z_eta-tmp_GEN_eta)*(Z_eta-tmp_GEN_eta) + (Z_phi-tmp_GEN_phi)*(Z_phi-tmp_GEN_phi) );
//if(!( tmpDeltaR_Vj>2.0 ))return 0;
}else{
if(!( GenGroomedJet_pt[0]>20 ))return 0;
efftool.Add_Event("Gen Jet Pt>20");
}

if(!( GenGroomedJet_mass[0]> 0.1 ))return 0;
efftool.Add_Event("Gen Jet mass> 0.1");

//alpha cut: z and jet are back to back
if(GroomedJet_number_jet_central >3){
//cout<<"alpha = "<<GenGroomedJet_pt[1]/GenGroomedJet_pt[0]<<endl; 
if( GroomedJet_pt[2]/GroomedJet_pt[1] >=0.3 ) return 0;
}
efftool.Add_Event("alpha<0.3");



}
else if ( FinalState.Contains("ZJet") ){

if(!( Z_mass>65 && Z_mass<105 ))return 0;
efftool.Add_Event("zmass: 65 105");

if(isBoosted){
if(!( Z_pt>100))return 0;
efftool.Add_Event("Z Pt>100");

if(!( GenGroomedJet_pt[0]>100 ))return 0;
efftool.Add_Event("Gen Jet Pt>100");
//Double_t tmp_GEN_eta = GenGroomedJet_eta[0]; Double_t tmp_GEN_phi = GenGroomedJet_phi[0];
//Double_t tmpDeltaR_Vj = TMath::Sqrt( (Z_eta-tmp_GEN_eta)*(Z_eta-tmp_GEN_eta) + (Z_phi-tmp_GEN_phi)*(Z_phi-tmp_GEN_phi) );
//if(!( tmpDeltaR_Vj>2.0 ))return 0;
}else{
if(!( GenGroomedJet_pt[0]>20 ))return 0;
efftool.Add_Event("Gen Jet Pt>20");
}

//alpha cut: z and jet are back to back
if(GroomedJet_number_jet_central >1){
//cout<<"alpha = "<<GenGroomedJet_pt[1]/GenGroomedJet_pt[0]<<endl; 
if( GroomedJet_pt[1]/GroomedJet_pt[0] >=0.3 ) return 0;
}
efftool.Add_Event("alpha<0.3");

}else{
cout<<"Wrong final state"<<endl;
return 0;
}

return 1;
} */
Bool_t MyClass::preSelect()
{
	efftool.Add_Event("Init");

	if( FinalState.Contains("Dijets") ){

		if(isBoosted==2){
			if(!( GroomedJet_pt[0]>300 && GroomedJet_pt[0]<500 ))return 0;
			efftool.Add_Event("recoJet Pt>300");
		}
		else if(isBoosted==1){
			//if(!( GroomedJet_pt[0]>100 ))return 0;
			if(!( GroomedJet_pt[0]>100 && GroomedJet_pt[0]<180 ))return 0;
			efftool.Add_Event("recoJet Pt>100");
		}else{
			//if(!( GroomedJet_pt[0]>20 ))return 0;
			if(!( GroomedJet_pt[0]>20 && GroomedJet_pt[0]<100 ))return 0;
			efftool.Add_Event("recoJet Pt>20");
		}

		if(!( GroomedJet_mass[0]> 0.1 ))return 0;
		efftool.Add_Event("recoJet mass> 0.1");

		//alpha cut: z and jet are back to back
		if(GroomedJet_number_jet_central >3){
			if( GroomedJet_pt[2]/GroomedJet_pt[1] >=0.3 ) return 0;
		}
		efftool.Add_Event("alpha<0.3");



	}
	else if ( FinalState.Contains("ZJet") ){

		if(!( Z_mass>65 && Z_mass<105 ))return 0;
		efftool.Add_Event("zmass: 65 105");

		if(isBoosted ==1){
			if(!( Z_pt>100))return 0;
			efftool.Add_Event("Z Pt>100");

			if(!( GroomedJet_pt[0]>100 ))return 0;
			efftool.Add_Event("recoJet Pt>100");
		}else{
			if(!( GroomedJet_pt[0]>20 ))return 0;
			efftool.Add_Event("recoJet Pt>20");
		}

		//alpha cut: z and jet are back to back
		if(GroomedJet_number_jet_central >1){
			if( GroomedJet_pt[1]/GroomedJet_pt[0] >=0.3 ) return 0;
		}
		efftool.Add_Event("alpha<0.3");

	}else{
		cout<<"Wrong final state"<<endl;
		return 0;
	}

	return 1;
}


/**/
void MyClass::Loop() {
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();

	/*Int_t nbin_rho=50; Double_t rho_min=0.; Double_t rho_max=50.;
	  Int_t nbin_nPV=50; Double_t nPV_min=0.; Double_t nPV_max=50.;
	  Int_t nbin_mass=60;Double_t jetmass_min=0;Double_t jetmass_max=300.;
	  Int_t nbin_pt=40;Double_t jetpt_min=50;Double_t jetpt_max=1050.;
	  Int_t nbin_ratio=100; Double_t ratio_min=-1; Double_t ratio_max=3.; 
	  Int_t nbin_eta=10;Double_t jeteta_min=-2.5;Double_t jeteta_max=2.5;
	  Int_t nbin_tau2tau1=40;Double_t jettau2tau1_min=0.;Double_t jettau2tau1_max=1.;
	  Double_t ratio_mrt_min=0.5; Double_t ratio_mrt_max=1.6; 
	  Double_t ratio_mrt_uncorr_min=0.5; Double_t ratio_mrt_uncorr_max=1.6; 
	  if(!isBoosted){
	  nbin_mass=40; jetmass_min=0; jetmass_max=300.;
	  nbin_pt=30; jetpt_min=0; jetpt_max=1050.;
	  }*/

	Int_t nbin_rho=60; Double_t rho_min=0.; Double_t rho_max=60.;
	Int_t nbin_nPV=60; Double_t nPV_min=0.; Double_t nPV_max=60.;
	Int_t nbin_mass=50;Double_t jetmass_min=-30;Double_t jetmass_max=120.;
	Int_t nbin_pt=30;Double_t jetpt_min=-100;Double_t jetpt_max= 500.;
	Int_t nbin_ratio=100; Double_t ratio_min=-1; Double_t ratio_max=3.; 
	Int_t nbin_ratio_mass=150; Double_t ratio_mass_min=-1; Double_t ratio_mass_max=5.; 
	Int_t nbin_eta=10;Double_t jeteta_min=-2.4;Double_t jeteta_max=2.4;
	Int_t nbin_phi=10;Double_t jetphi_min=-3.5;Double_t jetphi_max=7.0;
	Int_t nbin_tau2tau1=40;Double_t jettau2tau1_min=-0.2;Double_t jettau2tau1_max=1.2;
	Double_t ratio_mrt_min=0.5; Double_t ratio_mrt_max=1.6; 
	Double_t ratio_mass_mrt_min=0.; Double_t ratio_mass_mrt_max=3.5; 
	Double_t ratio_tau2tau1_mrt_min=0.; Double_t ratio_tau2tau1_mrt_max=3.0; 
	Double_t ratio_mrt_uncorr_min=0.5; Double_t ratio_mrt_uncorr_max=1.6; 

	JetCorrectionTool jct(FinalState.Data());

	RealVarArray rva_reco_pt_raw("reco_pt_raw",nbin_pt, jetpt_min, jetpt_max);
	RealVarArray rva_reco_pt_jecL1("reco_pt_jecL1",nbin_pt, jetpt_min, jetpt_max);
	RealVarArray rva_reco_pt_jecAll("reco_pt_jecAll",nbin_pt, jetpt_min, jetpt_max);
	RealVarArray rva_reco_pt_rhoswL1("reco_pt_rhoswL1",nbin_pt, jetpt_min, jetpt_max);
	RealVarArray rva_reco_pt_rhoHandL1("reco_pt_rhoHandL1",nbin_pt, jetpt_min, jetpt_max);
	RealVarArray rva_reco_pt_rhoHand2L1("reco_pt_rhoHand2L1",nbin_pt, jetpt_min, jetpt_max);
	RealVarArray rva_reco_pt_rhoGridL1("reco_pt_rhoGridL1",nbin_pt, jetpt_min, jetpt_max);
	RealVarArray rva_reco_pt_A4L1("reco_pt_A4L1",nbin_pt, jetpt_min, jetpt_max);
	RealVarArray rva_reco_pt_shapesubtraction("reco_pt_shapesubtraction",nbin_pt, jetpt_min, jetpt_max);
	RealVarArray rva_reco_pt_jetcleansing1("reco_pt_jetcleansing1",nbin_pt, jetpt_min, jetpt_max);
	RealVarArray rva_reco_pt_jetcleansing2("reco_pt_jetcleansing2",nbin_pt, jetpt_min, jetpt_max);
	RealVarArray rva_reco_pt_jetcleansing3("reco_pt_jetcleansing3",nbin_pt, jetpt_min, jetpt_max);
	RealVarArray rva_reco_pt_jetcleansing4("reco_pt_jetcleansing4",nbin_pt, jetpt_min, jetpt_max);
	RealVarArray rva_reco_pt_jetcleansing5("reco_pt_jetcleansing5",nbin_pt, jetpt_min, jetpt_max);
	RealVarArray rva_reco_pt_jetcleansing6("reco_pt_jetcleansing6",nbin_pt, jetpt_min, jetpt_max);
	RealVarArray rva_reco_pt_jetcleansing7("reco_pt_jetcleansing7",nbin_pt, jetpt_min, jetpt_max);
	RealVarArray rva_reco_pt_jetcleansing8("reco_pt_jetcleansing8",nbin_pt, jetpt_min, jetpt_max);
	RealVarArray rva_reco_pt_comb("reco_pt_comb",nbin_pt, jetpt_min, jetpt_max);

	RealVarArray rva_reco_mass_raw("reco_mass_raw",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_jecAll("reco_mass_jecAll",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_rhoHandL1("reco_mass_rhoHandL1",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_rhoGridL1("reco_mass_rhoGridL1",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_A4L1("reco_mass_A4L1",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_shapesubtraction("reco_mass_shapesubtraction",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_jetcleansing1("reco_mass_jetcleansing1",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_jetcleansing2("reco_mass_jetcleansing2",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_jetcleansing3("reco_mass_jetcleansing3",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_jetcleansing4("reco_mass_jetcleansing4",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_jetcleansing5("reco_mass_jetcleansing5",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_jetcleansing6("reco_mass_jetcleansing6",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_jetcleansing7("reco_mass_jetcleansing7",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_jetcleansing8("reco_mass_jetcleansing8",nbin_mass, jetmass_min, jetmass_max);

	RealVarArray rva_reco_tau2tau1_raw("reco_tau2tau1_raw",nbin_tau2tau1, jettau2tau1_min, jettau2tau1_max);
	RealVarArray rva_reco_tau2tau1_shapesubtraction("reco_tau2tau1_shapesubtraction",nbin_tau2tau1, jettau2tau1_min, jettau2tau1_max);
	RealVarArray rva_reco_tau2tau1_jetcleansing1("reco_tau2tau1_jetcleansing1",nbin_tau2tau1, jettau2tau1_min, jettau2tau1_max);
	RealVarArray rva_reco_tau2tau1_jetcleansing2("reco_tau2tau1_jetcleansing2",nbin_tau2tau1, jettau2tau1_min, jettau2tau1_max);
	RealVarArray rva_reco_tau2tau1_jetcleansing3("reco_tau2tau1_jetcleansing3",nbin_tau2tau1, jettau2tau1_min, jettau2tau1_max);
	RealVarArray rva_reco_tau2tau1_jetcleansing4("reco_tau2tau1_jetcleansing4",nbin_tau2tau1, jettau2tau1_min, jettau2tau1_max);
	RealVarArray rva_reco_tau2tau1_jetcleansing5("reco_tau2tau1_jetcleansing5",nbin_tau2tau1, jettau2tau1_min, jettau2tau1_max);
	RealVarArray rva_reco_tau2tau1_jetcleansing6("reco_tau2tau1_jetcleansing6",nbin_tau2tau1, jettau2tau1_min, jettau2tau1_max);
	RealVarArray rva_reco_tau2tau1_jetcleansing7("reco_tau2tau1_jetcleansing7",nbin_tau2tau1, jettau2tau1_min, jettau2tau1_max);
	RealVarArray rva_reco_tau2tau1_jetcleansing8("reco_tau2tau1_jetcleansing8",nbin_tau2tau1, jettau2tau1_min, jettau2tau1_max);
	RealVarArray rva_reco_eta4tau2tau1("reco_eta4tau2tau1",nbin_eta, jeteta_min, jeteta_max);

	jct.addVar(rva_reco_pt_raw); 
	jct.addVar(rva_reco_pt_jecL1);
	jct.addVar(rva_reco_pt_jecAll);
	jct.addVar(rva_reco_pt_rhoswL1);
	jct.addVar(rva_reco_pt_rhoHandL1);
	jct.addVar(rva_reco_pt_rhoHand2L1); 
	jct.addVar(rva_reco_pt_rhoGridL1);
	jct.addVar(rva_reco_pt_A4L1);
	jct.addVar(rva_reco_pt_shapesubtraction);
	jct.addVar(rva_reco_pt_jetcleansing1);
	jct.addVar(rva_reco_pt_jetcleansing2);
	jct.addVar(rva_reco_pt_jetcleansing3);
	jct.addVar(rva_reco_pt_jetcleansing4);
	jct.addVar(rva_reco_pt_jetcleansing5);
	jct.addVar(rva_reco_pt_jetcleansing6);
	jct.addVar(rva_reco_pt_jetcleansing7);
	jct.addVar(rva_reco_pt_jetcleansing8);
	jct.addVar(rva_reco_pt_comb);

	jct.addVar(rva_reco_mass_raw);
	jct.addVar(rva_reco_mass_jecAll);
	jct.addVar(rva_reco_mass_rhoHandL1);
	jct.addVar(rva_reco_mass_rhoGridL1);
	jct.addVar(rva_reco_mass_A4L1);
	jct.addVar(rva_reco_mass_shapesubtraction);
	jct.addVar(rva_reco_mass_jetcleansing1);
	jct.addVar(rva_reco_mass_jetcleansing2);
	jct.addVar(rva_reco_mass_jetcleansing3);
	jct.addVar(rva_reco_mass_jetcleansing4);
	jct.addVar(rva_reco_mass_jetcleansing5);
	jct.addVar(rva_reco_mass_jetcleansing6);
	jct.addVar(rva_reco_mass_jetcleansing7);
	jct.addVar(rva_reco_mass_jetcleansing8);

	jct.addVar(rva_reco_tau2tau1_raw);
	jct.addVar(rva_reco_tau2tau1_shapesubtraction);
	jct.addVar(rva_reco_tau2tau1_jetcleansing1);
	jct.addVar(rva_reco_tau2tau1_jetcleansing2);
	jct.addVar(rva_reco_tau2tau1_jetcleansing3);
	jct.addVar(rva_reco_tau2tau1_jetcleansing4);
	jct.addVar(rva_reco_tau2tau1_jetcleansing5);
	jct.addVar(rva_reco_tau2tau1_jetcleansing6);
	jct.addVar(rva_reco_tau2tau1_jetcleansing7);
	jct.addVar(rva_reco_tau2tau1_jetcleansing8);
	jct.addVar(rva_reco_eta4tau2tau1);

	RealVarArray rva_nPV("nPV",nbin_nPV, nPV_min, nPV_max);
	RealVarArray rva_reco_eta("reco_eta",nbin_eta, jeteta_min, jeteta_max);
	RealVarArray rva_reco_phi("reco_phi",nbin_phi, jetphi_min, jetphi_max);
	RealVarArray rva_gen_eta("gen_eta",nbin_eta, jeteta_min, jeteta_max);
	RealVarArray rva_gen_phi("gen_phi",nbin_phi, jetphi_min, jetphi_max);
	RealVarArray rva_gen_pt("gen_pt",nbin_pt, jetpt_min, jetpt_max);
	RealVarArray rva_gen_mass("gen_mass",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_gen_tau2tau1("gen_tau2tau1",nbin_tau2tau1, jettau2tau1_min, jettau2tau1_max);

	RealVarArray rva_reco_eta_jetcleansing1("reco_eta_jetcleansing1",nbin_eta, jeteta_min, jeteta_max);
	RealVarArray rva_reco_eta_jetcleansing2("reco_eta_jetcleansing2",nbin_eta, jeteta_min, jeteta_max);
	RealVarArray rva_reco_eta_jetcleansing3("reco_eta_jetcleansing3",nbin_eta, jeteta_min, jeteta_max);
	RealVarArray rva_reco_eta_jetcleansing4("reco_eta_jetcleansing4",nbin_eta, jeteta_min, jeteta_max);
	RealVarArray rva_reco_eta_jetcleansing5("reco_eta_jetcleansing5",nbin_eta, jeteta_min, jeteta_max);
	RealVarArray rva_reco_eta_jetcleansing6("reco_eta_jetcleansing6",nbin_eta, jeteta_min, jeteta_max);
	RealVarArray rva_reco_eta_jetcleansing7("reco_eta_jetcleansing7",nbin_eta, jeteta_min, jeteta_max);
	RealVarArray rva_reco_eta_jetcleansing8("reco_eta_jetcleansing8",nbin_eta, jeteta_min, jeteta_max);

	RealVarArray rva_reco_phi_jetcleansing1("reco_phi_jetcleansing1",nbin_phi, jetphi_min, jetphi_max);
	RealVarArray rva_reco_phi_jetcleansing2("reco_phi_jetcleansing2",nbin_phi, jetphi_min, jetphi_max);
	RealVarArray rva_reco_phi_jetcleansing3("reco_phi_jetcleansing3",nbin_phi, jetphi_min, jetphi_max);
	RealVarArray rva_reco_phi_jetcleansing4("reco_phi_jetcleansing4",nbin_phi, jetphi_min, jetphi_max);
	RealVarArray rva_reco_phi_jetcleansing5("reco_phi_jetcleansing5",nbin_phi, jetphi_min, jetphi_max);
	RealVarArray rva_reco_phi_jetcleansing6("reco_phi_jetcleansing6",nbin_phi, jetphi_min, jetphi_max);
	RealVarArray rva_reco_phi_jetcleansing7("reco_phi_jetcleansing7",nbin_phi, jetphi_min, jetphi_max);
	RealVarArray rva_reco_phi_jetcleansing8("reco_phi_jetcleansing8",nbin_phi, jetphi_min, jetphi_max);

	jct.addVar(rva_nPV);
	jct.addVar(rva_reco_eta);
	jct.addVar(rva_reco_phi);
	jct.addVar(rva_gen_eta);
	jct.addVar(rva_gen_phi);
	jct.addVar(rva_gen_pt);
	jct.addVar(rva_gen_mass);
	jct.addVar(rva_gen_tau2tau1);

	jct.addVar(rva_reco_eta_jetcleansing1);
	jct.addVar(rva_reco_eta_jetcleansing2);
	jct.addVar(rva_reco_eta_jetcleansing3);
	jct.addVar(rva_reco_eta_jetcleansing4);
	jct.addVar(rva_reco_eta_jetcleansing5);
	jct.addVar(rva_reco_eta_jetcleansing6);
	jct.addVar(rva_reco_eta_jetcleansing7);
	jct.addVar(rva_reco_eta_jetcleansing8);

	jct.addVar(rva_reco_phi_jetcleansing1);
	jct.addVar(rva_reco_phi_jetcleansing2);
	jct.addVar(rva_reco_phi_jetcleansing3);
	jct.addVar(rva_reco_phi_jetcleansing4);
	jct.addVar(rva_reco_phi_jetcleansing5);
	jct.addVar(rva_reco_phi_jetcleansing6);
	jct.addVar(rva_reco_phi_jetcleansing7);
	jct.addVar(rva_reco_phi_jetcleansing8);

	// gen and reco jet deltaR
	TH1D h1_RecoGen_matching("h1_RecoGen_matching","h1_RecoGen_matching; #delta R( reco j, gen j)",40,0,4.); h1_RecoGen_matching.SetLineColor(kRed);//matching with GEN

	Long64_t nbytes = 0, nb = 0;
	//for (Long64_t jentry=0; jentry<nentries;jentry++)
	//for (Long64_t jentry=0; jentry<nentries && jentry <100000;jentry++)
	for (Long64_t jentry=0; jentry<nentries && jentry <1000;jentry++)
	{
		//cout<<"jentry="<<jentry<<endl;
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if (!preSelect())continue;

		Double_t tmp_RECO_eta = GroomedJet_eta[0];
		Double_t tmp_RECO_phi = GroomedJet_phi[0];
		// gen jet matching with reco jet
		Int_t num_genreco_matching=-1;
		for(Int_t n=0;n<2;n++){
			Double_t tmp_GEN_eta = GenGroomedJet_eta[n];
			Double_t tmp_GEN_phi = GenGroomedJet_phi[n];
			//RECO-GEN Jet matching efficiency
			if( match_dR(tmp_RECO_eta, tmp_RECO_phi, tmp_GEN_eta, tmp_GEN_phi, 0.3, &h1_RecoGen_matching) ) {
				num_genreco_matching=n;
				break;
			}
		}
		if( num_genreco_matching==-1 ) { continue;}
		efftool.Add_Event("Reco-Gen jet matching");
		if( !(GenGroomedJet_mass[num_genreco_matching]>0.1) ){ continue;}
		efftool.Add_Event("Gen jet mass>0.1");
		//_matching PFCor with reco
		Bool_t isPFCor_matching_RECO=0;
		Int_t iposi_PFCorJet_matching_RECO=0;
		for(;iposi_PFCorJet_matching_RECO<numPFCorJets;iposi_PFCorJet_matching_RECO++){
			Double_t PFCor_jet_eta=JetPFCor_Eta[iposi_PFCorJet_matching_RECO];
			Double_t PFCor_jet_phi=JetPFCor_Phi[iposi_PFCorJet_matching_RECO];
			// phi and eta of PFCor(ak5PFCHS) and RECO should be matched
			if ( match_dR(tmp_RECO_eta, tmp_RECO_phi, PFCor_jet_eta, PFCor_jet_phi) ) { 
				isPFCor_matching_RECO=1;
				break;
			}
		}
		if( !isPFCor_matching_RECO) {continue;}
		efftool.Add_Event("PFCor-Reco jet matching");

		//cout<<"EventWeight="<<EventWeight<<endl;

		//cout<<"jentry="<<jentry<<"  pt_raw="<<GroomedJet_pt_uncorr[0]<<endl;
		//	cout<<"GroomedJet pt= "<<GenGroomedJet_pt[0]<<" , "<<GroomedJet_pt[0]<<" , "<<GroomedJet_pt_JetCleansing_DiffMode[0]<<" , "<<GroomedJet_pt_JetCleansing_DiffMode[1]<<" , "<<GroomedJet_pt_JetCleansing_DiffMode[2]<<endl;
		//	cout<<"GroomedJet eta= "<<GenGroomedJet_eta[0]<<" , "<<GroomedJet_eta[0]<<" , "<<GroomedJet_eta_JetCleansing_DiffMode[0]<<" , "<<GroomedJet_eta_JetCleansing_DiffMode[1]<<" , "<<GroomedJet_eta_JetCleansing_DiffMode[2]<<endl;
		//	cout<<"GroomedJet phi= "<<GenGroomedJet_phi[0]<<" , "<<GroomedJet_phi[0]<<" , "<<GroomedJet_phi_JetCleansing_DiffMode[0]<<" , "<<GroomedJet_phi_JetCleansing_DiffMode[1]<<" , "<<GroomedJet_phi_JetCleansing_DiffMode[2]<<endl;
		EventWeight=1.0;

		Bool_t debug =0.;
		//if( GroomedJet_mass_JetCleansing_DiffMode[2]/GenGroomedJet_mass[num_genreco_matching] < 0.2 ||GroomedJet_mass_JetCleansing_DiffMode[2]/GenGroomedJet_mass[num_genreco_matching] > 4 || GroomedJet_mass_shapesubtraction[0]/GenGroomedJet_mass[num_genreco_matching] < 0.2 ||GroomedJet_mass_shapesubtraction[0]/GenGroomedJet_mass[num_genreco_matching] > 4 )

		if( (GenGroomedJet_mass[num_genreco_matching]>40) && ( GroomedJet_tau2tau1_shapesubtraction[0] > 2 || GroomedJet_tau2tau1_shapesubtraction[0] < -2   
						||  GroomedJet_tau2tau1_JetCleansing_DiffMode[0]< -2 || GroomedJet_tau2tau1_JetCleansing_DiffMode[0]>2  ) )
		{
			debug =1;
			cout<<"--------------"<<endl<<"event_evtNo= "<<event_evtNo<<" event_lumi="<<event_lumi<<endl;
		}

		jct.fill( "reco_pt_raw", GroomedJet_pt_uncorr[0], EventWeight, debug);
		jct.fill( "reco_pt_jecL1", GroomedJet_pt_JECL1[0], EventWeight, debug); 
		jct.fill( "reco_pt_jecAll", GroomedJet_pt[0], EventWeight, debug);
		jct.fill( "reco_pt_rhoswL1", GroomedJet_pt_L1_rhoSW[0], EventWeight, debug);
		jct.fill( "reco_pt_rhoHandL1", GroomedJet_pt_L1_rhoHand[0], EventWeight, debug);
		jct.fill( "reco_pt_rhoHand2L1", GroomedJet_pt_L1_rhoHand2[0], EventWeight, debug);
		jct.fill( "reco_pt_rhoGridL1", GroomedJet_pt_L1_rhoGrid[0], EventWeight, debug);
		jct.fill( "reco_pt_A4L1", GroomedJet_pt_rho4A[0], EventWeight, debug);
		jct.fill( "reco_pt_shapesubtraction", GroomedJet_pt_shapesubtraction[0], EventWeight, debug);
		jct.fill( "reco_pt_jetcleansing1", GroomedJet_pt_JetCleansing_DiffMode[0], EventWeight, debug);//jvf r=0.3
		jct.fill( "reco_pt_jetcleansing2", GroomedJet_pt_JetCleansing_DiffMode[2], EventWeight, debug);;//jvf r=0.2 
		jct.fill( "reco_pt_jetcleansing3", GroomedJet_pt_JetCleansing_DiffMode[7], EventWeight, debug);//linear r=0.3, gamma0=0.55
		jct.fill( "reco_pt_jetcleansing4", GroomedJet_pt_JetCleansing_DiffMode[8], EventWeight, debug);//linear r=0.2, gamma0=0.55 
		jct.fill( "reco_pt_jetcleansing5", GroomedJet_pt_JetCleansing_DiffMode[9], EventWeight, debug);//jvf r=0.3
		jct.fill( "reco_pt_jetcleansing6", GroomedJet_pt_JetCleansing_DiffMode[25], EventWeight, debug);;//jvf r=0.2 
		jct.fill( "reco_pt_jetcleansing7", GroomedJet_pt_JetCleansing_DiffMode[26], EventWeight, debug);//linear r=0.3, gamma0=0.55
		jct.fill( "reco_pt_jetcleansing8", GroomedJet_pt_JetCleansing_DiffMode[27], EventWeight, debug);//linear r=0.2, gamma0=0.55 
		jct.fill( "reco_pt_comb", (GroomedJet_pt_JECL1[0]+GroomedJet_pt_JetCleansing_DiffMode[2])/2., EventWeight, debug);//linear r=0.2, gamma0=0.55 

		jct.fill( "reco_mass_raw", GroomedJet_mass_uncorr[0], EventWeight, debug);
		jct.fill( "reco_mass_jecAll", GroomedJet_mass[0], EventWeight, debug);
		jct.fill( "reco_mass_rhoHandL1", GroomedJet_mass_rhoArea[0], EventWeight, debug);
		jct.fill( "reco_mass_rhoGridL1", GroomedJet_mass_rhoGArea[0], EventWeight, debug);
		if( GroomedJet_mass_rho4Area[0]<0 ) GroomedJet_mass_rho4Area[0]=0.; 
		jct.fill( "reco_mass_A4L1", GroomedJet_mass_rho4Area[0], EventWeight, debug);
		if( GroomedJet_mass_shapesubtraction[0]<0 ) GroomedJet_mass_shapesubtraction[0]=0.; 
		jct.fill( "reco_mass_shapesubtraction", GroomedJet_mass_shapesubtraction[0], EventWeight, debug);
		jct.fill( "reco_mass_jetcleansing1", GroomedJet_mass_JetCleansing_DiffMode[0], EventWeight, debug);
		jct.fill( "reco_mass_jetcleansing2", GroomedJet_mass_JetCleansing_DiffMode[2], EventWeight, debug);
		jct.fill( "reco_mass_jetcleansing3", GroomedJet_mass_JetCleansing_DiffMode[7], EventWeight, debug);
		jct.fill( "reco_mass_jetcleansing4", GroomedJet_mass_JetCleansing_DiffMode[8], EventWeight, debug);
		jct.fill( "reco_mass_jetcleansing5", GroomedJet_mass_JetCleansing_DiffMode[9], EventWeight, debug);
		jct.fill( "reco_mass_jetcleansing6", GroomedJet_mass_JetCleansing_DiffMode[25], EventWeight, debug);
		jct.fill( "reco_mass_jetcleansing7", GroomedJet_mass_JetCleansing_DiffMode[26], EventWeight, debug);
		jct.fill( "reco_mass_jetcleansing8", GroomedJet_mass_JetCleansing_DiffMode[27], EventWeight, debug);

		//if( GenGroomedJet_mass[num_genreco_matching]>40 && GroomedJet_tau2tau1_JetCleansing_DiffMode[0]>-2000 && GroomedJet_tau2tau1_JetCleansing_DiffMode[0]<2000   )
		if( GenGroomedJet_mass[num_genreco_matching]>40 )
		{
			jct.fill( "reco_tau2tau1_raw", GroomedJet_tau2tau1[0], EventWeight, debug);
			if( GroomedJet_tau2tau1_shapesubtraction[0]<0. || GroomedJet_tau2tau1_shapesubtraction[0] >2 ) GroomedJet_tau2tau1_shapesubtraction[0]=0. ;
			jct.fill( "reco_tau2tau1_shapesubtraction", GroomedJet_tau2tau1_shapesubtraction[0], EventWeight, debug);
			jct.fill( "reco_tau2tau1_jetcleansing1", GroomedJet_tau2tau1_JetCleansing_DiffMode[0], EventWeight, debug);
			jct.fill( "reco_tau2tau1_jetcleansing2", GroomedJet_tau2tau1_JetCleansing_DiffMode[2], EventWeight, debug);
			jct.fill( "reco_tau2tau1_jetcleansing3", GroomedJet_tau2tau1_JetCleansing_DiffMode[7], EventWeight, debug);
			jct.fill( "reco_tau2tau1_jetcleansing4", GroomedJet_tau2tau1_JetCleansing_DiffMode[8], EventWeight, debug);
			jct.fill( "reco_tau2tau1_jetcleansing5", GroomedJet_tau2tau1_JetCleansing_DiffMode[9], EventWeight, debug);
			jct.fill( "reco_tau2tau1_jetcleansing6", GroomedJet_tau2tau1_JetCleansing_DiffMode[25], EventWeight, debug);
			jct.fill( "reco_tau2tau1_jetcleansing7", GroomedJet_tau2tau1_JetCleansing_DiffMode[26], EventWeight, debug);
			jct.fill( "reco_tau2tau1_jetcleansing8", GroomedJet_tau2tau1_JetCleansing_DiffMode[27], EventWeight, debug);
			jct.fill( "gen_tau2tau1", GenGroomedJet_tau2tau1[num_genreco_matching], EventWeight, debug);
			jct.fill( "reco_eta4tau2tau1", GroomedJet_eta[0], EventWeight, debug);
		}

		jct.fill( "gen_eta", GenGroomedJet_eta[num_genreco_matching], EventWeight, debug);
		//if( GenGroomedJet_phi[num_genreco_matching] <0 )GenGroomedJet_phi[num_genreco_matching] += 3.1415926*2; 
		jct.fill( "gen_phi", GenGroomedJet_phi[num_genreco_matching], EventWeight, debug);
		jct.fill( "gen_pt", GenGroomedJet_pt[num_genreco_matching], EventWeight, debug);
		jct.fill( "gen_mass", GenGroomedJet_mass[num_genreco_matching], EventWeight, debug);

		jct.fill( "nPV", event_nPV, EventWeight, debug);

		jct.fill( "reco_eta", GroomedJet_eta[0], EventWeight, debug);
		jct.fill( "reco_phi", GroomedJet_phi[0], EventWeight, debug);
		jct.fill( "reco_eta_jetcleansing1", GroomedJet_eta_JetCleansing_DiffMode[0], EventWeight, debug);//jvf r=0.3
		jct.fill( "reco_eta_jetcleansing2", GroomedJet_eta_JetCleansing_DiffMode[2], EventWeight, debug);;//jvf r=0.2 
		jct.fill( "reco_eta_jetcleansing3", GroomedJet_eta_JetCleansing_DiffMode[7], EventWeight, debug);//linear r=0.3, gamma0=0.55
		jct.fill( "reco_eta_jetcleansing4", GroomedJet_eta_JetCleansing_DiffMode[8], EventWeight, debug);//linear r=0.2, gamma0=0.55 
		jct.fill( "reco_eta_jetcleansing5", GroomedJet_eta_JetCleansing_DiffMode[9], EventWeight, debug);//jvf r=0.3
		jct.fill( "reco_eta_jetcleansing6", GroomedJet_eta_JetCleansing_DiffMode[25], EventWeight, debug);;//jvf r=0.2 
		jct.fill( "reco_eta_jetcleansing7", GroomedJet_eta_JetCleansing_DiffMode[26], EventWeight, debug);//linear r=0.3, gamma0=0.55
		jct.fill( "reco_eta_jetcleansing8", GroomedJet_eta_JetCleansing_DiffMode[27], EventWeight, debug);//linear r=0.2, gamma0=0.55 

		jct.fill( "reco_phi_jetcleansing1", GroomedJet_phi_JetCleansing_DiffMode[0], EventWeight, debug);//jvf r=0.3
		jct.fill( "reco_phi_jetcleansing2", GroomedJet_phi_JetCleansing_DiffMode[2], EventWeight, debug);;//jvf r=0.2 
		jct.fill( "reco_phi_jetcleansing3", GroomedJet_phi_JetCleansing_DiffMode[7], EventWeight, debug);//linear r=0.3, gamma0=0.55
		jct.fill( "reco_phi_jetcleansing4", GroomedJet_phi_JetCleansing_DiffMode[8], EventWeight, debug);//linear r=0.2, gamma0=0.55 
		jct.fill( "reco_phi_jetcleansing5", GroomedJet_phi_JetCleansing_DiffMode[9], EventWeight, debug);//jvf r=0.3
		jct.fill( "reco_phi_jetcleansing6", GroomedJet_phi_JetCleansing_DiffMode[25], EventWeight, debug);;//jvf r=0.2 
		jct.fill( "reco_phi_jetcleansing7", GroomedJet_phi_JetCleansing_DiffMode[26], EventWeight, debug);//linear r=0.3, gamma0=0.55
		jct.fill( "reco_phi_jetcleansing8", GroomedJet_phi_JetCleansing_DiffMode[27], EventWeight, debug);//linear r=0.2, gamma0=0.55 

	}

	//Draw all variable distri
	vector< TString > vect_allvar=jct.get_allvar_name();
	for( Int_t m=0;m <vect_allvar.size(); m++){
		cout<<vect_allvar[m].Data()<<endl; 
		Draw_and_Save(jct.get_hist1D(vect_allvar[m].Data()) );
	}


	vector< TString > vect_pt_corrected;
	vect_pt_corrected.clear();
	vect_pt_corrected.push_back( TString(Form("h1_JCT_%s_%s", FinalState.Data(), "gen_pt")) );
	vect_pt_corrected.push_back( TString(Form("h1_JCT_%s_%s", FinalState.Data(), "reco_pt_A4L")) );
	vect_pt_corrected.push_back( TString(Form("h1_JCT_%s_%s", FinalState.Data(), "reco_pt_jecL1")) );
	vect_pt_corrected.push_back( TString(Form("h1_JCT_%s_%s", FinalState.Data(), "reco_pt_comb")) );
	vect_pt_corrected.push_back( TString(Form("h1_JCT_%s_%s", FinalState.Data(), "reco_pt_jetcleansing2")) );
	DrawPlots(vect_pt_corrected);

	//	Draw_and_Save(h1_PFCor_jec_recogenptratio);
	//	Draw_and_Save(h1_PFCor_uncorr_recogenptratio);
	//	Draw_and_Save(h1_PFCor_afterL1_recogenptratio);
	//	Draw_and_Save(h1_PFCor_afterL2_recogenptratio);

	//	Draw_and_Save(mrt_PFCor_jec_recogenptratio.get_mr_hist(ratio_mrt_min, ratio_mrt_max));
	//	Draw_and_Save(mrt_PFCor_uncorr_recogenptratio.get_mr_hist(ratio_mrt_uncorr_min, ratio_mrt_uncorr_max));
	//	Draw_and_Save(mrt_PFCor_afterL1_recogenptratio.get_mr_hist(ratio_mrt_min, ratio_mrt_max));
	//	Draw_and_Save(mrt_PFCor_afterL2_recogenptratio.get_mr_hist(ratio_mrt_min, ratio_mrt_max));

	// pt, mass, tau2tau1 responce
	Draw_and_Save( jct.get_hist1D_response("gen_pt", "reco_pt_A4L1", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_pt", "reco_pt_jecAll", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_pt", "reco_pt_jecL1", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_pt", "reco_pt_jetcleansing1", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_pt", "reco_pt_jetcleansing2", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_pt", "reco_pt_jetcleansing3", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_pt", "reco_pt_jetcleansing4", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_pt", "reco_pt_jetcleansing5", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_pt", "reco_pt_jetcleansing6", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_pt", "reco_pt_jetcleansing7", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_pt", "reco_pt_jetcleansing8", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_pt", "reco_pt_comb", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_pt", "reco_pt_raw", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_pt", "reco_pt_rhoGridL1", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_pt", "reco_pt_rhoHand2L1", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_pt", "reco_pt_rhoHandL1", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_pt", "reco_pt_rhoswL1", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_pt", "reco_pt_shapesubtraction", nbin_ratio, ratio_min, ratio_max));

	Draw_and_Save( jct.get_hist1D_response("gen_mass", "reco_mass_A4L1", nbin_ratio_mass, ratio_mass_min, ratio_mass_max));
	Draw_and_Save( jct.get_hist1D_response("gen_mass", "reco_mass_jecAll", nbin_ratio_mass, ratio_mass_min, ratio_mass_max));
	Draw_and_Save( jct.get_hist1D_response("gen_mass", "reco_mass_jetcleansing1", nbin_ratio_mass, ratio_mass_min, ratio_mass_max));
	Draw_and_Save( jct.get_hist1D_response("gen_mass", "reco_mass_jetcleansing2", nbin_ratio_mass, ratio_mass_min, ratio_mass_max));
	Draw_and_Save( jct.get_hist1D_response("gen_mass", "reco_mass_jetcleansing3", nbin_ratio_mass, ratio_mass_min, ratio_mass_max));
	Draw_and_Save( jct.get_hist1D_response("gen_mass", "reco_mass_jetcleansing4", nbin_ratio_mass, ratio_mass_min, ratio_mass_max));
	Draw_and_Save( jct.get_hist1D_response("gen_mass", "reco_mass_jetcleansing5", nbin_ratio_mass, ratio_mass_min, ratio_mass_max));
	Draw_and_Save( jct.get_hist1D_response("gen_mass", "reco_mass_jetcleansing6", nbin_ratio_mass, ratio_mass_min, ratio_mass_max));
	Draw_and_Save( jct.get_hist1D_response("gen_mass", "reco_mass_jetcleansing7", nbin_ratio_mass, ratio_mass_min, ratio_mass_max));
	Draw_and_Save( jct.get_hist1D_response("gen_mass", "reco_mass_jetcleansing8", nbin_ratio_mass, ratio_mass_min, ratio_mass_max));
	Draw_and_Save( jct.get_hist1D_response("gen_mass", "reco_mass_raw", nbin_ratio_mass, ratio_mass_min, ratio_mass_max));
	Draw_and_Save( jct.get_hist1D_response("gen_mass", "reco_mass_rhoGridL1", nbin_ratio_mass, ratio_mass_min, ratio_mass_max));
	Draw_and_Save( jct.get_hist1D_response("gen_mass", "reco_mass_rhoHandL1", nbin_ratio_mass, ratio_mass_min, ratio_mass_max));
	Draw_and_Save( jct.get_hist1D_response("gen_mass", "reco_mass_shapesubtraction", nbin_ratio_mass, ratio_mass_min, ratio_mass_max));

	Draw_and_Save( jct.get_hist1D_response("gen_tau2tau1", "reco_tau2tau1_jetcleansing1", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_tau2tau1", "reco_tau2tau1_jetcleansing2", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_tau2tau1", "reco_tau2tau1_jetcleansing3", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_tau2tau1", "reco_tau2tau1_jetcleansing4", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_tau2tau1", "reco_tau2tau1_jetcleansing5", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_tau2tau1", "reco_tau2tau1_jetcleansing6", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_tau2tau1", "reco_tau2tau1_jetcleansing7", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_tau2tau1", "reco_tau2tau1_jetcleansing8", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_tau2tau1", "reco_tau2tau1_raw", nbin_ratio, ratio_min, ratio_max));
	Draw_and_Save( jct.get_hist1D_response("gen_tau2tau1", "reco_tau2tau1_shapesubtraction", nbin_ratio, ratio_min, ratio_max));

	// pt, mass, tau2tau1 2D: reco vs gen
	Draw_and_Save( jct.get_hist2D("gen_pt", "reco_pt_A4L1"  		), Form("%g correlated", jct.get_graph("gen_pt", "reco_pt_A4L1"  		).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_pt", "reco_pt_jecAll"		), Form("%g correlated", jct.get_graph("gen_pt", "reco_pt_jecAll"		).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_pt", "reco_pt_jecL1"			), Form("%g correlated", jct.get_graph("gen_pt", "reco_pt_jecL1"			).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_pt", "reco_pt_jetcleansing1" ), Form("%g correlated", jct.get_graph("gen_pt", "reco_pt_jetcleansing1" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_pt", "reco_pt_jetcleansing2" ), Form("%g correlated", jct.get_graph("gen_pt", "reco_pt_jetcleansing2" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_pt", "reco_pt_jetcleansing3" ), Form("%g correlated", jct.get_graph("gen_pt", "reco_pt_jetcleansing3" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_pt", "reco_pt_jetcleansing4" ), Form("%g correlated", jct.get_graph("gen_pt", "reco_pt_jetcleansing4" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_pt", "reco_pt_jetcleansing5" ), Form("%g correlated", jct.get_graph("gen_pt", "reco_pt_jetcleansing5" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_pt", "reco_pt_jetcleansing6" ), Form("%g correlated", jct.get_graph("gen_pt", "reco_pt_jetcleansing6" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_pt", "reco_pt_jetcleansing7" ), Form("%g correlated", jct.get_graph("gen_pt", "reco_pt_jetcleansing7" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_pt", "reco_pt_jetcleansing8" ), Form("%g correlated", jct.get_graph("gen_pt", "reco_pt_jetcleansing8" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_pt", "reco_pt_comb" ), Form("%g correlated", jct.get_graph("gen_pt", "reco_pt_comb" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_pt", "reco_pt_raw"			), Form("%g correlated", jct.get_graph("gen_pt", "reco_pt_raw"			).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_pt", "reco_pt_rhoGridL1"		), Form("%g correlated", jct.get_graph("gen_pt", "reco_pt_rhoGridL1"		).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_pt", "reco_pt_rhoHand2L1"	), Form("%g correlated", jct.get_graph("gen_pt", "reco_pt_rhoHand2L1"	).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_pt", "reco_pt_rhoHandL1"		), Form("%g correlated", jct.get_graph("gen_pt", "reco_pt_rhoHandL1"		).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_pt", "reco_pt_rhoswL1"		), Form("%g correlated", jct.get_graph("gen_pt", "reco_pt_rhoswL1"		).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_pt", "reco_pt_shapesubtraction" ), Form("%g correlated", jct.get_graph("gen_pt", "reco_pt_shapesubtraction" ).GetCorrelationFactor()) );

	Draw_and_Save( jct.get_hist2D("gen_mass", "reco_mass_A4L1"			), Form("%g correlated", jct.get_graph("gen_mass", "reco_mass_A4L1"			).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_mass", "reco_mass_jecAll"		), Form("%g correlated", jct.get_graph("gen_mass", "reco_mass_jecAll"		).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_mass", "reco_mass_jetcleansing1" ), Form("%g correlated", jct.get_graph("gen_mass", "reco_mass_jetcleansing1" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_mass", "reco_mass_jetcleansing2" ), Form("%g correlated", jct.get_graph("gen_mass", "reco_mass_jetcleansing2" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_mass", "reco_mass_jetcleansing3" ), Form("%g correlated", jct.get_graph("gen_mass", "reco_mass_jetcleansing3" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_mass", "reco_mass_jetcleansing4" ), Form("%g correlated", jct.get_graph("gen_mass", "reco_mass_jetcleansing4" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_mass", "reco_mass_jetcleansing5" ), Form("%g correlated", jct.get_graph("gen_mass", "reco_mass_jetcleansing5" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_mass", "reco_mass_jetcleansing6" ), Form("%g correlated", jct.get_graph("gen_mass", "reco_mass_jetcleansing6" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_mass", "reco_mass_jetcleansing7" ), Form("%g correlated", jct.get_graph("gen_mass", "reco_mass_jetcleansing7" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_mass", "reco_mass_jetcleansing8" ), Form("%g correlated", jct.get_graph("gen_mass", "reco_mass_jetcleansing8" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_mass", "reco_mass_raw"			), Form("%g correlated", jct.get_graph("gen_mass", "reco_mass_raw"			).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_mass", "reco_mass_rhoGridL1"		), Form("%g correlated", jct.get_graph("gen_mass", "reco_mass_rhoGridL1"		).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_mass", "reco_mass_rhoHandL1"		), Form("%g correlated", jct.get_graph("gen_mass", "reco_mass_rhoHandL1"		).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_mass", "reco_mass_shapesubtraction" ), Form("%g correlated", jct.get_graph("gen_mass", "reco_mass_shapesubtraction" ).GetCorrelationFactor()) );

	Draw_and_Save( jct.get_hist2D("gen_tau2tau1", "reco_tau2tau1_jetcleansing1" ), Form("%g correlated", jct.get_graph("gen_tau2tau1", "reco_tau2tau1_jetcleansing1" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_tau2tau1", "reco_tau2tau1_jetcleansing2" ), Form("%g correlated", jct.get_graph("gen_tau2tau1", "reco_tau2tau1_jetcleansing2" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_tau2tau1", "reco_tau2tau1_jetcleansing3" ), Form("%g correlated", jct.get_graph("gen_tau2tau1", "reco_tau2tau1_jetcleansing3" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_tau2tau1", "reco_tau2tau1_jetcleansing4" ), Form("%g correlated", jct.get_graph("gen_tau2tau1", "reco_tau2tau1_jetcleansing4" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_tau2tau1", "reco_tau2tau1_jetcleansing5" ), Form("%g correlated", jct.get_graph("gen_tau2tau1", "reco_tau2tau1_jetcleansing5" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_tau2tau1", "reco_tau2tau1_jetcleansing6" ), Form("%g correlated", jct.get_graph("gen_tau2tau1", "reco_tau2tau1_jetcleansing6" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_tau2tau1", "reco_tau2tau1_jetcleansing7" ), Form("%g correlated", jct.get_graph("gen_tau2tau1", "reco_tau2tau1_jetcleansing7" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_tau2tau1", "reco_tau2tau1_jetcleansing8" ), Form("%g correlated", jct.get_graph("gen_tau2tau1", "reco_tau2tau1_jetcleansing8" ).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_tau2tau1", "reco_tau2tau1_raw"			), Form("%g correlated", jct.get_graph("gen_tau2tau1", "reco_tau2tau1_raw"			).GetCorrelationFactor()) );
	Draw_and_Save( jct.get_hist2D("gen_tau2tau1", "reco_tau2tau1_shapesubtraction" ), Form("%g correlated", jct.get_graph("gen_tau2tau1", "reco_tau2tau1_shapesubtraction" ).GetCorrelationFactor()) );


	// pt, mass, tau2tau1 vs eta
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_A4L1", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_jecAll", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_jecL1", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_jetcleansing1", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_jetcleansing2", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_jetcleansing3", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_jetcleansing4", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_jetcleansing5", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_jetcleansing6", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_jetcleansing7", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_jetcleansing8", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_comb", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_raw", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_rhoGridL1", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_rhoHand2L1", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_rhoHandL1", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_rhoswL1", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_shapesubtraction", ratio_mrt_min, ratio_mrt_max));

	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_A4L1", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_jecAll", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_jetcleansing1", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_jetcleansing2", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_jetcleansing3", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_jetcleansing4", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_jetcleansing5", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_jetcleansing6", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_jetcleansing7", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_jetcleansing8", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_raw", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_rhoGridL1", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_rhoHandL1", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_shapesubtraction", ratio_mass_mrt_min, ratio_mass_mrt_max));

	/*// pt, mass, tau2tau1 vs pt 
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_A4L1", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_jecAll", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_jecL1", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_jetcleansing1", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_jetcleansing2", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_jetcleansing3", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_jetcleansing4", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_jetcleansing5", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_jetcleansing6", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_jetcleansing7", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_jetcleansing8", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_comb", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_raw", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_rhoGridL1", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_rhoHand2L1", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_rhoHandL1", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_rhoswL1", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_shapesubtraction", ratio_mrt_min, ratio_mrt_max));

	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_A4L1", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_jecAll", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_jetcleansing1", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_jetcleansing2", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_jetcleansing3", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_jetcleansing4", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_jetcleansing5", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_jetcleansing6", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_jetcleansing7", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_jetcleansing8", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_raw", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_rhoGridL1", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_rhoHandL1", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_shapesubtraction", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  */

	Draw_and_Save( jct.get_mean_rms_hist("reco_eta4tau2tau1","gen_tau2tau1", "reco_tau2tau1_jetcleansing1", ratio_tau2tau1_mrt_min, ratio_tau2tau1_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta4tau2tau1","gen_tau2tau1", "reco_tau2tau1_jetcleansing2", ratio_tau2tau1_mrt_min, ratio_tau2tau1_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta4tau2tau1","gen_tau2tau1", "reco_tau2tau1_jetcleansing3", ratio_tau2tau1_mrt_min, ratio_tau2tau1_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta4tau2tau1","gen_tau2tau1", "reco_tau2tau1_jetcleansing4", ratio_tau2tau1_mrt_min, ratio_tau2tau1_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta4tau2tau1","gen_tau2tau1", "reco_tau2tau1_jetcleansing5", ratio_tau2tau1_mrt_min, ratio_tau2tau1_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta4tau2tau1","gen_tau2tau1", "reco_tau2tau1_jetcleansing6", ratio_tau2tau1_mrt_min, ratio_tau2tau1_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta4tau2tau1","gen_tau2tau1", "reco_tau2tau1_jetcleansing7", ratio_tau2tau1_mrt_min, ratio_tau2tau1_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta4tau2tau1","gen_tau2tau1", "reco_tau2tau1_jetcleansing8", ratio_tau2tau1_mrt_min, ratio_tau2tau1_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta4tau2tau1","gen_tau2tau1", "reco_tau2tau1_raw", ratio_tau2tau1_mrt_min, ratio_tau2tau1_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta4tau2tau1","gen_tau2tau1", "reco_tau2tau1_shapesubtraction", ratio_tau2tau1_mrt_min, ratio_tau2tau1_mrt_max));

	Draw_and_Save( efftool.Get_Eff_hist() );



	//	Draw_and_Save(h1_PFCor_area);
	//	Draw_and_Save(h1_RECO_area);
	//	Draw_and_Save(h1_RecoGen_matching);
	//	Draw_and_Save(h1_RECO_eta);
	//	Draw_and_Save(h1_RECO_zjet_dr);
	//	Draw_and_Save(h1_RECO_zjet_dphi);
	//	Draw_and_Save(h1_RECO_rhoSW);
	//	Draw_and_Save(h1_RECO_rhoHand);
	//	Draw_and_Save(h1_RECO_rhoHand2);
	//	Draw_and_Save(h1_RECO_rhoGrid);
	//	Draw_and_Save(h2_RECO_rhoSW_vs_nPV);
	//	Draw_and_Save(h2_RECO_rhoHand_vs_nPV);
	//	Draw_and_Save(h2_RECO_rhoHand2_vs_nPV);
	//	Draw_and_Save(h2_RECO_rhoGrid_vs_nPV);
	//
	//	Draw_and_Save(h2_RECO_l1rhoHand_recogenptratio_vs_nPV);
	//	Draw_and_Save(h2_RECO_l1rhoHand_recogenptratio_vs_ptHand);
	//	Draw_and_Save(h2_RECO_l1rhoHand_recogenptratio_vs_eta);
	//
	//	Draw_and_Save(h1_PFCor_mass              );
	//	Draw_and_Save(h2_RECO_mass_jec_vs_PV);


	fout->Write();// fout->Close();

} /**/



void MyClass::Draw_and_Print_All() { }

