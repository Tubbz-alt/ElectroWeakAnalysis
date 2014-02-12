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

void MyClass::Draw_and_Save(TH1D h1, char* addtional_info){
	h1.Write();
	TCanvas *c1 = new TCanvas(Form("c1_%s",h1.GetTitle()),Form("c1_%s",h1.GetTitle()),200,10,600,600);
	c1->cd();
	h1.Draw();
	if( addtional_info){
		TLatex tl;
		tl.SetTextSize(0.04 ); tl.SetTextAlign(13);
		//tl.DrawLatex(h1.GetXaxis()->GetXmin()*0.9+h1.GetXaxis()->GetXmax()*0.1,h1.GetYaxis()->GetXmin()*0.1+h1.GetYaxis()->GetXmax()*0.9,addtional_info);
		tl.DrawLatex(h1.GetXaxis()->GetXmin()*0.9+h1.GetXaxis()->GetXmax()*0.1, h1.GetMinimum()*0.1 +h1.GetMaximum()*0.9,addtional_info);
		//cout<<"xmin="<<h1.GetXaxis()->GetXmin()<<endl;  cout<<"xmax="<<h1.GetXaxis()->GetXmax()<<endl;  cout<<"ymin="<<h1.GetYaxis()->GetXmin()<<endl;  cout<<"ymax="<<h1.GetYaxis()->GetXmax()<<endl;  
	};
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

void MyClass::DrawPlots(vector< PlotConfig > plotconfig, char* xaxis_title)
{
	TString finalname("");
	//TCanvas *c1 = new TCanvas(Form("c1_%s",h2.GetTitle()),Form("c1_%s",h2.GetTitle()),200,10,600,600);
	TCanvas *c1 = new TCanvas(Form("c1_multiplots"),Form("c1_multiplots"),200,10,600,600);
	c1->cd();

	vector<TH1*> vect_hist;
	TLegend *leg=new TLegend(0.65,0.6,0.9,0.9);
	Double_t max_yval=0.;

	leg->SetName("theLegend"); leg->SetBorderSize(0); leg->SetLineColor(0); leg->SetFillColor(0);
	leg->SetFillStyle(0); leg->SetLineWidth(0); leg->SetLineStyle(0); leg->SetTextFont(42);
	leg->SetTextSize(.045);

	for(Int_t i=0; i< Int_t(plotconfig.size()); i++){
		finalname+=plotconfig[i].name;
		TH1* h1;
		cout<<"plotconfig[i]="<<plotconfig[i].name.Data()<<endl;
		fout->GetObject(plotconfig[i].name.Data(), h1);
		if(h1){
			h1->SetLineColor( plotconfig[i].linecolor);
			h1->SetLineStyle( plotconfig[i].linestyle);
			vect_hist.push_back(h1);
			leg->AddEntry(h1,plotconfig[i].title,"lep");
			Double_t tmpmax=h1->GetMaximum();
			if(tmpmax>max_yval) max_yval=tmpmax;
			//if (i==0) h1->Draw(); else h1->Draw("same");
		}else{ cout<<"Can't find "<<plotconfig[i].name.Data()<<endl; }
	}

	vect_hist[0]->GetYaxis()->SetRangeUser(0., max_yval*1.2);
	if(xaxis_title)vect_hist[0]->GetXaxis()->SetTitle(xaxis_title);
	for(Int_t i=0; i< Int_t(vect_hist.size()); i++){
		if (i==0) vect_hist[i]->Draw(); 
		else vect_hist[i]->Draw("same");
	}
	leg->Draw();

	c1->Print(Form("%s/multiplots_%s_%s_%s.png",plot_Dir_DateTime.Data(), JetType.Data(), PfType.Data(),finalname.Data()));
	delete c1;
}


/* Bool_t MyClass::preSelect() {
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
		} else if(isBoosted==1){
			if(!( GroomedJet_pt[0]>100 && GroomedJet_pt[0]<180 ))return 0;
			efftool.Add_Event("recoJet Pt>100");
		}else{
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
	Int_t nbin_mass=50;Double_t jetmass_min=-10;Double_t jetmass_max=290.;
	Int_t nbin_pt=30;Double_t jetpt_min=0;Double_t jetpt_max= 900.;
	Int_t nbin_ratio=100; Double_t ratio_min=-1; Double_t ratio_max=3.; 
	Int_t nbin_ratio_mass=150; Double_t ratio_mass_min=-1; Double_t ratio_mass_max=5.; 
	Int_t nbin_eta=10;Double_t jeteta_min=-2.4;Double_t jeteta_max=2.4;
	Int_t nbin_phi=10;Double_t jetphi_min=-3.5;Double_t jetphi_max=7.0;
	Int_t nbin_tau2tau1=40;Double_t jettau2tau1_min=-0.2;Double_t jettau2tau1_max=1.2;
	Double_t ratio_mrt_min=0.5; Double_t ratio_mrt_max=2.0; 
	Double_t ratio_mass_mrt_min=0.; Double_t ratio_mass_mrt_max=4.0; 
	Double_t ratio_tau2tau1_mrt_min=0.; Double_t ratio_tau2tau1_mrt_max=3.0; 

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
	RealVarArray rva_reco_pt_jetcleansing9("reco_pt_jetcleansing9",nbin_pt, jetpt_min, jetpt_max);
	RealVarArray rva_reco_pt_comb("reco_pt_comb",nbin_pt, jetpt_min, jetpt_max);

	RealVarArray rva_reco_mass_raw("reco_mass_raw",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_jecAll("reco_mass_jecAll",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_rhoHandL1("reco_mass_rhoHandL1",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_rhoGridL1("reco_mass_rhoGridL1",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_A4L1("reco_mass_A4L1",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_A4mL1("reco_mass_A4mL1",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_shapesubtraction("reco_mass_shapesubtraction",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_jetcleansing1("reco_mass_jetcleansing1",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_jetcleansing2("reco_mass_jetcleansing2",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_jetcleansing3("reco_mass_jetcleansing3",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_jetcleansing4("reco_mass_jetcleansing4",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_jetcleansing5("reco_mass_jetcleansing5",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_jetcleansing6("reco_mass_jetcleansing6",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_jetcleansing7("reco_mass_jetcleansing7",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_jetcleansing8("reco_mass_jetcleansing8",nbin_mass, jetmass_min, jetmass_max);
	RealVarArray rva_reco_mass_jetcleansing9("reco_mass_jetcleansing9",nbin_mass, jetmass_min, jetmass_max);

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
	RealVarArray rva_reco_tau2tau1_jetcleansing9("reco_tau2tau1_jetcleansing9",nbin_tau2tau1, jettau2tau1_min, jettau2tau1_max);
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
	jct.addVar(rva_reco_pt_jetcleansing9);
	jct.addVar(rva_reco_pt_comb);

	jct.addVar(rva_reco_mass_raw);
	jct.addVar(rva_reco_mass_jecAll);
	jct.addVar(rva_reco_mass_rhoHandL1);
	jct.addVar(rva_reco_mass_rhoGridL1);
	jct.addVar(rva_reco_mass_A4L1);
	jct.addVar(rva_reco_mass_A4mL1);
	jct.addVar(rva_reco_mass_shapesubtraction);
	jct.addVar(rva_reco_mass_jetcleansing1);
	jct.addVar(rva_reco_mass_jetcleansing2);
	jct.addVar(rva_reco_mass_jetcleansing3);
	jct.addVar(rva_reco_mass_jetcleansing4);
	jct.addVar(rva_reco_mass_jetcleansing5);
	jct.addVar(rva_reco_mass_jetcleansing6);
	jct.addVar(rva_reco_mass_jetcleansing7);
	jct.addVar(rva_reco_mass_jetcleansing8);
	jct.addVar(rva_reco_mass_jetcleansing9);

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
	jct.addVar(rva_reco_tau2tau1_jetcleansing9);
	jct.addVar(rva_reco_eta4tau2tau1);

	RealVarArray rva_nPV("nPV",nbin_nPV, nPV_min, nPV_max);
	RealVarArray rva_rho("rho",nbin_rho, rho_min, rho_max);
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
	RealVarArray rva_reco_eta_jetcleansing9("reco_eta_jetcleansing9",nbin_eta, jeteta_min, jeteta_max);

	RealVarArray rva_reco_phi_jetcleansing1("reco_phi_jetcleansing1",nbin_phi, jetphi_min, jetphi_max);
	RealVarArray rva_reco_phi_jetcleansing2("reco_phi_jetcleansing2",nbin_phi, jetphi_min, jetphi_max);
	RealVarArray rva_reco_phi_jetcleansing3("reco_phi_jetcleansing3",nbin_phi, jetphi_min, jetphi_max);
	RealVarArray rva_reco_phi_jetcleansing4("reco_phi_jetcleansing4",nbin_phi, jetphi_min, jetphi_max);
	RealVarArray rva_reco_phi_jetcleansing5("reco_phi_jetcleansing5",nbin_phi, jetphi_min, jetphi_max);
	RealVarArray rva_reco_phi_jetcleansing6("reco_phi_jetcleansing6",nbin_phi, jetphi_min, jetphi_max);
	RealVarArray rva_reco_phi_jetcleansing7("reco_phi_jetcleansing7",nbin_phi, jetphi_min, jetphi_max);
	RealVarArray rva_reco_phi_jetcleansing8("reco_phi_jetcleansing8",nbin_phi, jetphi_min, jetphi_max);
	RealVarArray rva_reco_phi_jetcleansing9("reco_phi_jetcleansing9",nbin_phi, jetphi_min, jetphi_max);

	jct.addVar(rva_nPV);
	jct.addVar(rva_rho);
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
	jct.addVar(rva_reco_eta_jetcleansing9);

	jct.addVar(rva_reco_phi_jetcleansing1);
	jct.addVar(rva_reco_phi_jetcleansing2);
	jct.addVar(rva_reco_phi_jetcleansing3);
	jct.addVar(rva_reco_phi_jetcleansing4);
	jct.addVar(rva_reco_phi_jetcleansing5);
	jct.addVar(rva_reco_phi_jetcleansing6);
	jct.addVar(rva_reco_phi_jetcleansing7);
	jct.addVar(rva_reco_phi_jetcleansing8);
	jct.addVar(rva_reco_phi_jetcleansing9);

	// gen and reco jet deltaR
	TH1D h1_RecoGen_matching("h1_RecoGen_matching","h1_RecoGen_matching; #delta R( reco j, gen j)",40,0,4.); h1_RecoGen_matching.SetLineColor(kRed);//matching with GEN
	//Cleansing Diff-Mode array
	Int_t cleansing_diff_mode[9]={0, 2, 3, 12, 13, 14, 15, 16, 17};

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++)
	//for (Long64_t jentry=0; jentry<nentries && jentry <20000;jentry++)
	{
		//cout<<"jentry="<<jentry<<endl;
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if (!preSelect())continue;

		// gen jet matching with reco jet
		Double_t tmp_RECO_eta = GroomedJet_eta[0];
		Double_t tmp_RECO_phi = GroomedJet_phi[0];
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
		// Gen mass >0.1
		if( !(GenGroomedJet_mass[num_genreco_matching]>0.1) ){ continue;}
		efftool.Add_Event("Gen jet mass>0.1");
		// Matching PFCor(AK5PFChs) with reco
		Int_t num_pfcorreco_matching=-1;
		for(Int_t ip=0;ip<numPFCorJets;ip++){
			Double_t PFCor_jet_eta=JetPFCor_Eta[ip];
			Double_t PFCor_jet_phi=JetPFCor_Phi[ip];
			// phi and eta of PFCor(ak5PFCHS) and RECO should be matched
			if ( match_dR(tmp_RECO_eta, tmp_RECO_phi, PFCor_jet_eta, PFCor_jet_phi) ) { 
				num_pfcorreco_matching=ip;
				break;
			}
		}
		if( num_pfcorreco_matching==-1) {continue;}
		efftool.Add_Event("PFCor-Reco jet matching");


		//cout<<"jentry="<<jentry<<"  pt_raw="<<GroomedJet_pt_uncorr[0]<<endl;
		//	cout<<"GroomedJet pt= "<<GenGroomedJet_pt[0]<<" , "<<GroomedJet_pt[0]<<" , "<<GroomedJet_pt_JetCleansing_DiffMode[0]<<" , "<<GroomedJet_pt_JetCleansing_DiffMode[1]<<" , "<<GroomedJet_pt_JetCleansing_DiffMode[2]<<endl;
		//	cout<<"GroomedJet eta= "<<GenGroomedJet_eta[0]<<" , "<<GroomedJet_eta[0]<<" , "<<GroomedJet_eta_JetCleansing_DiffMode[0]<<" , "<<GroomedJet_eta_JetCleansing_DiffMode[1]<<" , "<<GroomedJet_eta_JetCleansing_DiffMode[2]<<endl;
		//	cout<<"GroomedJet phi= "<<GenGroomedJet_phi[0]<<" , "<<GroomedJet_phi[0]<<" , "<<GroomedJet_phi_JetCleansing_DiffMode[0]<<" , "<<GroomedJet_phi_JetCleansing_DiffMode[1]<<" , "<<GroomedJet_phi_JetCleansing_DiffMode[2]<<endl;

		//cout<<"EventWeight="<<EventWeight<<endl;
		EventWeight=1.0;

		Bool_t debug =0.;
		if( (GenGroomedJet_mass[num_genreco_matching]>40) && ( GroomedJet_tau2tau1_shapesubtraction[0] > 2 || GroomedJet_tau2tau1_shapesubtraction[0] < -2
						||  GroomedJet_tau2tau1_JetCleansing_DiffMode[0]< -2 || GroomedJet_tau2tau1_JetCleansing_DiffMode[0]>2  ) ) {
			debug =1;
			cout<<"--------------"<<endl<<"event_evtNo= "<<event_evtNo<<" event_lumi="<<event_lumi<<endl;
			cout<<"tau2tau1= "<<GroomedJet_tau2tau1_shapesubtraction[0] <<", "<< GroomedJet_tau2tau1_JetCleansing_DiffMode[0]<<endl;
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
		jct.fill( "reco_pt_jetcleansing1", GroomedJet_pt_JetCleansing_DiffMode[cleansing_diff_mode[0]], EventWeight, debug);//jvf r=0.3
		jct.fill( "reco_pt_jetcleansing2", GroomedJet_pt_JetCleansing_DiffMode[cleansing_diff_mode[1]], EventWeight, debug);;//jvf r=0.2 
		jct.fill( "reco_pt_jetcleansing3", GroomedJet_pt_JetCleansing_DiffMode[cleansing_diff_mode[2]], EventWeight, debug);//linear r=0.3, gamma0=0.55
		jct.fill( "reco_pt_jetcleansing4", GroomedJet_pt_JetCleansing_DiffMode[cleansing_diff_mode[3]], EventWeight, debug);//linear r=0.2, gamma0=0.55 
		jct.fill( "reco_pt_jetcleansing5", GroomedJet_pt_JetCleansing_DiffMode[cleansing_diff_mode[4]], EventWeight, debug);//jvf r=0.3
		jct.fill( "reco_pt_jetcleansing6", GroomedJet_pt_JetCleansing_DiffMode[cleansing_diff_mode[5]], EventWeight, debug);;//jvf r=0.2 
		jct.fill( "reco_pt_jetcleansing7", GroomedJet_pt_JetCleansing_DiffMode[cleansing_diff_mode[6]], EventWeight, debug);//linear r=0.3, gamma0=0.55
		jct.fill( "reco_pt_jetcleansing8", GroomedJet_pt_JetCleansing_DiffMode[cleansing_diff_mode[7]], EventWeight, debug);//linear r=0.2, gamma0=0.55 
		jct.fill( "reco_pt_jetcleansing9", GroomedJet_pt_JetCleansing_DiffMode[cleansing_diff_mode[8]], EventWeight, debug);//linear r=0.2, gamma0=0.55 
		jct.fill( "reco_pt_comb", (GroomedJet_pt_JECL1[0]+GroomedJet_pt_JetCleansing_DiffMode[2])/2., EventWeight, debug);//linear r=0.2, gamma0=0.55 

		jct.fill( "reco_mass_raw", GroomedJet_mass_uncorr[0], EventWeight, debug);
		jct.fill( "reco_mass_jecAll", GroomedJet_mass[0], EventWeight, debug);
		jct.fill( "reco_mass_rhoHandL1", GroomedJet_mass_rhoArea[0], EventWeight, debug);
		jct.fill( "reco_mass_rhoGridL1", GroomedJet_mass_rhoGArea[0], EventWeight, debug);
		if( GroomedJet_mass_rho4A[0]<0 ) GroomedJet_mass_rho4A[0]=0.; 
		jct.fill( "reco_mass_A4L1", GroomedJet_mass_rho4A[0], EventWeight, debug);
		if( GroomedJet_mass_rhom4Am[0]<0 ) GroomedJet_mass_rhom4Am[0]=0.; 
		jct.fill( "reco_mass_A4mL1", GroomedJet_mass_rhom4Am[0], EventWeight, debug);
		if( GroomedJet_mass_shapesubtraction[0]<0 ) GroomedJet_mass_shapesubtraction[0]=0.; 
		jct.fill( "reco_mass_shapesubtraction", GroomedJet_mass_shapesubtraction[0], EventWeight, debug);
		jct.fill( "reco_mass_jetcleansing1", GroomedJet_mass_JetCleansing_DiffMode[cleansing_diff_mode[0]], EventWeight, debug);
		jct.fill( "reco_mass_jetcleansing2", GroomedJet_mass_JetCleansing_DiffMode[cleansing_diff_mode[1]], EventWeight, debug);
		jct.fill( "reco_mass_jetcleansing3", GroomedJet_mass_JetCleansing_DiffMode[cleansing_diff_mode[2]], EventWeight, debug);
		jct.fill( "reco_mass_jetcleansing4", GroomedJet_mass_JetCleansing_DiffMode[cleansing_diff_mode[3]], EventWeight, debug);
		jct.fill( "reco_mass_jetcleansing5", GroomedJet_mass_JetCleansing_DiffMode[cleansing_diff_mode[4]], EventWeight, debug);
		jct.fill( "reco_mass_jetcleansing6", GroomedJet_mass_JetCleansing_DiffMode[cleansing_diff_mode[5]], EventWeight, debug);
		jct.fill( "reco_mass_jetcleansing7", GroomedJet_mass_JetCleansing_DiffMode[cleansing_diff_mode[6]], EventWeight, debug);
		jct.fill( "reco_mass_jetcleansing8", GroomedJet_mass_JetCleansing_DiffMode[cleansing_diff_mode[7]], EventWeight, debug);
		jct.fill( "reco_mass_jetcleansing9", GroomedJet_mass_JetCleansing_DiffMode[cleansing_diff_mode[8]], EventWeight, debug);

		//if( GenGroomedJet_mass[num_genreco_matching]>40 && GroomedJet_tau2tau1_JetCleansing_DiffMode[0]>-2000 && GroomedJet_tau2tau1_JetCleansing_DiffMode[0]<2000   )
		if( GenGroomedJet_mass[num_genreco_matching]>40 )
		{
			jct.fill( "reco_tau2tau1_raw", GroomedJet_tau2tau1[0], EventWeight, debug);
			//if( GroomedJet_tau2tau1_shapesubtraction[0]<0. || GroomedJet_tau2tau1_shapesubtraction[0] >2 ) GroomedJet_tau2tau1_shapesubtraction[0]=0. ;
			if( GroomedJet_tau2tau1_shapesubtraction[0]<0. ) GroomedJet_tau2tau1_shapesubtraction[0]=0. ;
			jct.fill( "reco_tau2tau1_shapesubtraction", GroomedJet_tau2tau1_shapesubtraction[0], EventWeight, debug);
			jct.fill( "reco_tau2tau1_jetcleansing1", GroomedJet_tau2tau1_JetCleansing_DiffMode[cleansing_diff_mode[0]], EventWeight, debug);
			jct.fill( "reco_tau2tau1_jetcleansing2", GroomedJet_tau2tau1_JetCleansing_DiffMode[cleansing_diff_mode[1]], EventWeight, debug);
			jct.fill( "reco_tau2tau1_jetcleansing3", GroomedJet_tau2tau1_JetCleansing_DiffMode[cleansing_diff_mode[2]], EventWeight, debug);
			jct.fill( "reco_tau2tau1_jetcleansing4", GroomedJet_tau2tau1_JetCleansing_DiffMode[cleansing_diff_mode[3]], EventWeight, debug);
			jct.fill( "reco_tau2tau1_jetcleansing5", GroomedJet_tau2tau1_JetCleansing_DiffMode[cleansing_diff_mode[4]], EventWeight, debug);
			jct.fill( "reco_tau2tau1_jetcleansing6", GroomedJet_tau2tau1_JetCleansing_DiffMode[cleansing_diff_mode[5]], EventWeight, debug);
			jct.fill( "reco_tau2tau1_jetcleansing7", GroomedJet_tau2tau1_JetCleansing_DiffMode[cleansing_diff_mode[6]], EventWeight, debug);
			jct.fill( "reco_tau2tau1_jetcleansing8", GroomedJet_tau2tau1_JetCleansing_DiffMode[cleansing_diff_mode[7]], EventWeight, debug);
			jct.fill( "reco_tau2tau1_jetcleansing9", GroomedJet_tau2tau1_JetCleansing_DiffMode[cleansing_diff_mode[8]], EventWeight, debug);
			jct.fill( "gen_tau2tau1", GenGroomedJet_tau2tau1[num_genreco_matching], EventWeight, debug);
			jct.fill( "reco_eta4tau2tau1", GroomedJet_eta[0], EventWeight, debug);
		}

		jct.fill( "gen_eta", GenGroomedJet_eta[num_genreco_matching], EventWeight, debug);
		//if( GenGroomedJet_phi[num_genreco_matching] <0 )GenGroomedJet_phi[num_genreco_matching] += 3.1415926*2; 
		jct.fill( "gen_phi", GenGroomedJet_phi[num_genreco_matching], EventWeight, debug);
		jct.fill( "gen_pt", GenGroomedJet_pt[num_genreco_matching], EventWeight, debug);
		jct.fill( "gen_mass", GenGroomedJet_mass[num_genreco_matching], EventWeight, debug);

		jct.fill( "nPV", event_nPV, EventWeight, debug);
		jct.fill( "rho", GroomedJet_rhohand, EventWeight, debug);

		jct.fill( "reco_eta", GroomedJet_eta[0], EventWeight, debug);
		jct.fill( "reco_phi", GroomedJet_phi[0], EventWeight, debug);
		jct.fill( "reco_eta_jetcleansing1", GroomedJet_eta_JetCleansing_DiffMode[cleansing_diff_mode[0]], EventWeight, debug);//jvf r=0.3
		jct.fill( "reco_eta_jetcleansing2", GroomedJet_eta_JetCleansing_DiffMode[cleansing_diff_mode[1]], EventWeight, debug);;//jvf r=0.2 
		jct.fill( "reco_eta_jetcleansing3", GroomedJet_eta_JetCleansing_DiffMode[cleansing_diff_mode[2]], EventWeight, debug);//linear r=0.3, gamma0=0.55
		jct.fill( "reco_eta_jetcleansing4", GroomedJet_eta_JetCleansing_DiffMode[cleansing_diff_mode[3]], EventWeight, debug);//linear r=0.2, gamma0=0.55 
		jct.fill( "reco_eta_jetcleansing5", GroomedJet_eta_JetCleansing_DiffMode[cleansing_diff_mode[4]], EventWeight, debug);//jvf r=0.3
		jct.fill( "reco_eta_jetcleansing6", GroomedJet_eta_JetCleansing_DiffMode[cleansing_diff_mode[5]], EventWeight, debug);;//jvf r=0.2 
		jct.fill( "reco_eta_jetcleansing7", GroomedJet_eta_JetCleansing_DiffMode[cleansing_diff_mode[6]], EventWeight, debug);//linear r=0.3, gamma0=0.55
		jct.fill( "reco_eta_jetcleansing8", GroomedJet_eta_JetCleansing_DiffMode[cleansing_diff_mode[7]], EventWeight, debug);//linear r=0.2, gamma0=0.55 
		jct.fill( "reco_eta_jetcleansing9", GroomedJet_eta_JetCleansing_DiffMode[cleansing_diff_mode[8]], EventWeight, debug);//linear r=0.2, gamma0=0.55 

		jct.fill( "reco_phi_jetcleansing1", GroomedJet_phi_JetCleansing_DiffMode[cleansing_diff_mode[0]], EventWeight, debug);//jvf r=0.3
		jct.fill( "reco_phi_jetcleansing2", GroomedJet_phi_JetCleansing_DiffMode[cleansing_diff_mode[1]], EventWeight, debug);;//jvf r=0.2 
		jct.fill( "reco_phi_jetcleansing3", GroomedJet_phi_JetCleansing_DiffMode[cleansing_diff_mode[2]], EventWeight, debug);//linear r=0.3, gamma0=0.55
		jct.fill( "reco_phi_jetcleansing4", GroomedJet_phi_JetCleansing_DiffMode[cleansing_diff_mode[3]], EventWeight, debug);//linear r=0.2, gamma0=0.55 
		jct.fill( "reco_phi_jetcleansing5", GroomedJet_phi_JetCleansing_DiffMode[cleansing_diff_mode[4]], EventWeight, debug);//jvf r=0.3
		jct.fill( "reco_phi_jetcleansing6", GroomedJet_phi_JetCleansing_DiffMode[cleansing_diff_mode[5]], EventWeight, debug);;//jvf r=0.2 
		jct.fill( "reco_phi_jetcleansing7", GroomedJet_phi_JetCleansing_DiffMode[cleansing_diff_mode[6]], EventWeight, debug);//linear r=0.3, gamma0=0.55
		jct.fill( "reco_phi_jetcleansing8", GroomedJet_phi_JetCleansing_DiffMode[cleansing_diff_mode[7]], EventWeight, debug);//linear r=0.2, gamma0=0.55 
		jct.fill( "reco_phi_jetcleansing9", GroomedJet_phi_JetCleansing_DiffMode[cleansing_diff_mode[8]], EventWeight, debug);//linear r=0.2, gamma0=0.55 

	}

	//Draw all variable distri
	cout<<"=========== Draw all variables distribution ============="<<endl;
	vector< TString > vect_allvar=jct.get_allvar_name();
	for( Int_t m=0;m <Int_t(vect_allvar.size()); m++){
		cout<<vect_allvar[m].Data()<<endl; 
		Draw_and_Save(jct.get_hist1D(vect_allvar[m].Data()) );
	}

	cout<<"=========== Draw Plots for comparing diff. correction ============="<<endl;
	Int_t colorlist[10]={1,2,3,4,6,7,8,10,13,15};
	vector< PlotConfig > vect_pt_corrected;
	vect_pt_corrected.clear();
	PlotConfig plot_pt_0(Form("h1_JCT_%s_%s", FinalState.Data(), "gen_pt"),                        "Gen",colorlist[0],1);
	PlotConfig plot_pt_1(Form("h1_JCT_%s_%s", FinalState.Data(), "reco_pt_raw"),              "RECO RAW",colorlist[1],2);
	PlotConfig plot_pt_2(Form("h1_JCT_%s_%s", FinalState.Data(), "reco_pt_A4L1"),          "4-vect area",colorlist[2],3);
	PlotConfig plot_pt_3(Form("h1_JCT_%s_%s", FinalState.Data(), "reco_pt_jecL1"),               "JECL1",colorlist[3],4);  
	PlotConfig plot_pt_4(Form("h1_JCT_%s_%s", FinalState.Data(), "reco_pt_shapesubtraction"), "Shape Subtraction",colorlist[4],5);
	PlotConfig plot_pt_5(Form("h1_JCT_%s_%s", FinalState.Data(), "reco_pt_jetcleansing2"),"JetCleansing",colorlist[5],6); 
	vect_pt_corrected.push_back( plot_pt_0);
	vect_pt_corrected.push_back( plot_pt_1);
	vect_pt_corrected.push_back( plot_pt_2);
	vect_pt_corrected.push_back( plot_pt_3);
	vect_pt_corrected.push_back( plot_pt_4);
	vect_pt_corrected.push_back( plot_pt_5);
	DrawPlots(vect_pt_corrected, "pT");

	vector< PlotConfig > vect_mass_corrected;
	vect_mass_corrected.clear();
	PlotConfig plot_mass_0(Form("h1_JCT_%s_%s", FinalState.Data(), "gen_mass"),                        "Gen",colorlist[0],1);
	PlotConfig plot_mass_1(Form("h1_JCT_%s_%s", FinalState.Data(), "reco_mass_raw"),              "RECO RAW",colorlist[1],2);
	PlotConfig plot_mass_2(Form("h1_JCT_%s_%s", FinalState.Data(), "reco_mass_A4L1"),          "4-vect area",colorlist[2],3);
	PlotConfig plot_mass_3(Form("h1_JCT_%s_%s", FinalState.Data(), "reco_mass_jecAll"),               "JEC",colorlist[3],4);  
	PlotConfig plot_mass_4(Form("h1_JCT_%s_%s", FinalState.Data(), "reco_mass_shapesubtraction"), "Shape Subtraction",colorlist[4],5);
	PlotConfig plot_mass_5(Form("h1_JCT_%s_%s", FinalState.Data(), "reco_mass_jetcleansing2"),"JetCleansing",colorlist[5],6); 
	vect_mass_corrected.push_back( plot_mass_0);
	vect_mass_corrected.push_back( plot_mass_1);
	vect_mass_corrected.push_back( plot_mass_2);
	vect_mass_corrected.push_back( plot_mass_3);
	vect_mass_corrected.push_back( plot_mass_4);
	vect_mass_corrected.push_back( plot_mass_5);
	DrawPlots(vect_mass_corrected, "mass");

	vector< PlotConfig > vect_tau2tau1_corrected;
	vect_tau2tau1_corrected.clear();
	PlotConfig plot_tau2tau1_0(Form("h1_JCT_%s_%s", FinalState.Data(), "gen_tau2tau1"),                        "Gen",colorlist[0],1);
	PlotConfig plot_tau2tau1_1(Form("h1_JCT_%s_%s", FinalState.Data(), "reco_tau2tau1_raw"),              "RECO RAW",colorlist[1],2);
	PlotConfig plot_tau2tau1_2(Form("h1_JCT_%s_%s", FinalState.Data(), "reco_tau2tau1_shapesubtraction"), "Shape Subtraction",colorlist[2],3);
	PlotConfig plot_tau2tau1_3(Form("h1_JCT_%s_%s", FinalState.Data(), "reco_tau2tau1_jetcleansing2"),"JetCleansing",colorlist[5],6); 
	vect_tau2tau1_corrected.push_back( plot_tau2tau1_0);
	vect_tau2tau1_corrected.push_back( plot_tau2tau1_1);
	vect_tau2tau1_corrected.push_back( plot_tau2tau1_2);
	vect_tau2tau1_corrected.push_back( plot_tau2tau1_3);
	DrawPlots(vect_tau2tau1_corrected, "tau2tau1");

	cout<<"=========== Draw Response Plots ============="<<endl; // pt, mass, tau2tau1 responce
	//Print out a table of different correction response
	Table_Tool table_ptreponse;
	SaveResponse(table_ptreponse, jct, "gen_pt", "reco_pt_A4L1", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_ptreponse, jct, "gen_pt", "reco_pt_jecAll", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_ptreponse, jct, "gen_pt", "reco_pt_jecL1", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_ptreponse, jct, "gen_pt", "reco_pt_jetcleansing1", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_ptreponse, jct, "gen_pt", "reco_pt_jetcleansing2", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_ptreponse, jct, "gen_pt", "reco_pt_jetcleansing3", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_ptreponse, jct, "gen_pt", "reco_pt_jetcleansing4", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_ptreponse, jct, "gen_pt", "reco_pt_jetcleansing5", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_ptreponse, jct, "gen_pt", "reco_pt_jetcleansing6", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_ptreponse, jct, "gen_pt", "reco_pt_jetcleansing7", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_ptreponse, jct, "gen_pt", "reco_pt_jetcleansing8", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_ptreponse, jct, "gen_pt", "reco_pt_jetcleansing9", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_ptreponse, jct, "gen_pt", "reco_pt_comb", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_ptreponse, jct, "gen_pt", "reco_pt_raw", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_ptreponse, jct, "gen_pt", "reco_pt_rhoGridL1", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_ptreponse, jct, "gen_pt", "reco_pt_rhoHand2L1", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_ptreponse, jct, "gen_pt", "reco_pt_rhoHandL1", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_ptreponse, jct, "gen_pt", "reco_pt_rhoswL1", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_ptreponse, jct, "gen_pt", "reco_pt_shapesubtraction", nbin_ratio, ratio_min, ratio_max);
	table_ptreponse.PrintTable();

	Table_Tool table_massreponse;
	SaveResponse(table_massreponse, jct, "gen_mass", "reco_mass_A4L1", nbin_ratio_mass, ratio_mass_min, ratio_mass_max);
	SaveResponse(table_massreponse, jct, "gen_mass", "reco_mass_A4mL1", nbin_ratio_mass, ratio_mass_min, ratio_mass_max);
	SaveResponse(table_massreponse, jct, "gen_mass", "reco_mass_jecAll", nbin_ratio_mass, ratio_mass_min, ratio_mass_max);
	SaveResponse(table_massreponse, jct, "gen_mass", "reco_mass_jetcleansing1", nbin_ratio_mass, ratio_mass_min, ratio_mass_max);
	SaveResponse(table_massreponse, jct, "gen_mass", "reco_mass_jetcleansing2", nbin_ratio_mass, ratio_mass_min, ratio_mass_max);
	SaveResponse(table_massreponse, jct, "gen_mass", "reco_mass_jetcleansing3", nbin_ratio_mass, ratio_mass_min, ratio_mass_max);
	SaveResponse(table_massreponse, jct, "gen_mass", "reco_mass_jetcleansing4", nbin_ratio_mass, ratio_mass_min, ratio_mass_max);
	SaveResponse(table_massreponse, jct, "gen_mass", "reco_mass_jetcleansing5", nbin_ratio_mass, ratio_mass_min, ratio_mass_max);
	SaveResponse(table_massreponse, jct, "gen_mass", "reco_mass_jetcleansing6", nbin_ratio_mass, ratio_mass_min, ratio_mass_max);
	SaveResponse(table_massreponse, jct, "gen_mass", "reco_mass_jetcleansing7", nbin_ratio_mass, ratio_mass_min, ratio_mass_max);
	SaveResponse(table_massreponse, jct, "gen_mass", "reco_mass_jetcleansing8", nbin_ratio_mass, ratio_mass_min, ratio_mass_max);
	SaveResponse(table_massreponse, jct, "gen_mass", "reco_mass_jetcleansing9", nbin_ratio_mass, ratio_mass_min, ratio_mass_max);
	SaveResponse(table_massreponse, jct, "gen_mass", "reco_mass_raw", nbin_ratio_mass, ratio_mass_min, ratio_mass_max);
	SaveResponse(table_massreponse, jct, "gen_mass", "reco_mass_rhoGridL1", nbin_ratio_mass, ratio_mass_min, ratio_mass_max);
	SaveResponse(table_massreponse, jct, "gen_mass", "reco_mass_rhoHandL1", nbin_ratio_mass, ratio_mass_min, ratio_mass_max);
	SaveResponse(table_massreponse, jct, "gen_mass", "reco_mass_shapesubtraction", nbin_ratio_mass, ratio_mass_min, ratio_mass_max);
	table_massreponse.PrintTable();

	Table_Tool table_tau2tau1reponse;
	SaveResponse(table_tau2tau1reponse, jct, "gen_tau2tau1", "reco_tau2tau1_jetcleansing1", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_tau2tau1reponse, jct, "gen_tau2tau1", "reco_tau2tau1_jetcleansing2", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_tau2tau1reponse, jct, "gen_tau2tau1", "reco_tau2tau1_jetcleansing3", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_tau2tau1reponse, jct, "gen_tau2tau1", "reco_tau2tau1_jetcleansing4", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_tau2tau1reponse, jct, "gen_tau2tau1", "reco_tau2tau1_jetcleansing5", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_tau2tau1reponse, jct, "gen_tau2tau1", "reco_tau2tau1_jetcleansing6", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_tau2tau1reponse, jct, "gen_tau2tau1", "reco_tau2tau1_jetcleansing7", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_tau2tau1reponse, jct, "gen_tau2tau1", "reco_tau2tau1_jetcleansing8", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_tau2tau1reponse, jct, "gen_tau2tau1", "reco_tau2tau1_jetcleansing9", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_tau2tau1reponse, jct, "gen_tau2tau1", "reco_tau2tau1_raw", nbin_ratio, ratio_min, ratio_max);
	SaveResponse(table_tau2tau1reponse, jct, "gen_tau2tau1", "reco_tau2tau1_shapesubtraction", nbin_ratio, ratio_min, ratio_max);
	table_tau2tau1reponse.PrintTable();





	cout<<"=========== Draw Eta-Depenence Plots ============="<<endl; // pt, mass, tau2tau1 vs eta
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
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_jetcleansing9", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_comb", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_raw", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_rhoGridL1", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_rhoHand2L1", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_rhoHandL1", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_rhoswL1", ratio_mrt_min, ratio_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_pt", "reco_pt_shapesubtraction", ratio_mrt_min, ratio_mrt_max));

	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_A4L1", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_A4mL1", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_jecAll", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_jetcleansing1", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_jetcleansing2", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_jetcleansing3", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_jetcleansing4", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_jetcleansing5", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_jetcleansing6", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_jetcleansing7", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_jetcleansing8", ratio_mass_mrt_min, ratio_mass_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta","gen_mass", "reco_mass_jetcleansing9", ratio_mass_mrt_min, ratio_mass_mrt_max));
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
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_jetcleansing9", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_comb", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_raw", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_rhoGridL1", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_rhoHand2L1", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_rhoHandL1", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_rhoswL1", ratio_mrt_min, ratio_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_pt", "reco_pt_shapesubtraction", ratio_mrt_min, ratio_mrt_max));

	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_A4L1", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_A4mL1", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_jecAll", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_jetcleansing1", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_jetcleansing2", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_jetcleansing3", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_jetcleansing4", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_jetcleansing5", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_jetcleansing6", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_jetcleansing7", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_jetcleansing8", ratio_mass_mrt_min, ratio_mass_mrt_max));
	  Draw_and_Save( jct.get_mean_rms_hist("gen_pt","gen_mass", "reco_mass_jetcleansing9", ratio_mass_mrt_min, ratio_mass_mrt_max));
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
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta4tau2tau1","gen_tau2tau1", "reco_tau2tau1_jetcleansing9", ratio_tau2tau1_mrt_min, ratio_tau2tau1_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta4tau2tau1","gen_tau2tau1", "reco_tau2tau1_raw", ratio_tau2tau1_mrt_min, ratio_tau2tau1_mrt_max));
	Draw_and_Save( jct.get_mean_rms_hist("reco_eta4tau2tau1","gen_tau2tau1", "reco_tau2tau1_shapesubtraction", ratio_tau2tau1_mrt_min, ratio_tau2tau1_mrt_max));

	cout<<"=========== Draw other suport Plots ============="<<endl; 
	Draw_and_Save( efftool.Get_Eff_hist() );
	Draw_and_Save(h1_RecoGen_matching);
	//	Draw_and_Save(h1_PFCor_area);
	//	Draw_and_Save(h1_RECO_area);
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

void MyClass::SaveResponse(Table_Tool& table, JetCorrectionTool &jct, TString xdenominator_var_name,TString xnumerator_var_name , Int_t nbin, Double_t xmin, Double_t xmax ){
	RESPONSE responce= jct.get_response(xdenominator_var_name.Data(), xnumerator_var_name.Data(), nbin, xmin, xmax);
	Draw_and_Save( responce.hist, Form("Mean=%g, RMS=%g", responce.mean, responce.rms ));
	Draw_and_Save( responce.hist2D, Form("%g correlated", responce.correlationfactor) );

	TString x_mean=xdenominator_var_name+"_mean";
	TString x_rms=xdenominator_var_name+"_rms";
	TString x_co=xdenominator_var_name+"_correlation";
	table.Insert(x_mean, xnumerator_var_name, responce.mean);
	table.Insert(x_rms, xnumerator_var_name, responce.rms/responce.mean);
	table.Insert(x_co, xnumerator_var_name, responce.correlationfactor);
}

