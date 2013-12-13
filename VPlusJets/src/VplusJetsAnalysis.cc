 /*****************************************************************************
 * Project: CMS detector at the CERN
 *
 * Package: ElectroWeakAnalysis/VPlusJets
 *
 *
 * Authors:
 *
 *   Kalanand Mishra, Fermilab - kalanand@fnal.gov
 *
 * Description:
 *   To fill jet related quantities into a specified TTree
 *   Can work with CaloJet, GenJet, JPT jet, PF jet.
 *   Can work with jets in RECO/AOD/PAT data formats.
 * History:
 *   
 *
 * Copyright (C) 2010 FNAL 
 *****************************************************************************/


// user include files
#include "ElectroWeakAnalysis/VPlusJets/interface/VplusJetsAnalysis.h" 

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TString.h"

//#include "JetSubstructure/SubstructureTools/interface/PseudoJetUserInfo.h"
//#include "JetSubstructure/SubstructureTools/interface/JetSubstructureTools.h" 

using namespace std;
using namespace edm;

/*
   ewk::VplusJetsAnalysis::VplusJetsAnalysis(const edm::ParameterSet& iConfig) :
   myTree ( fs -> mkdir("../").make<TTree>(iConfig.getParameter<std::string>("TreeName").c_str(),"V+jets Tree") ), 
   CorrectedPFJetFiller ( iConfig.existsAs<edm::InputTag>("srcPFCor") ?  new JetTreeFiller("CorrectedPFJetFiller", myTree, "PFCor", iConfig, 0) : 0),
   genAK5groomedJetFiller ((iConfig.existsAs<bool>("doGroomedAK5")&& iConfig.getParameter< bool >("doGroomedAK5")) ?  new GroomedJetFiller("genGroomedJetFiller", myTree, "AK5", "_GEN", iConfig, 1) : 0), //Gen
   AK5groomedJetFiller_PF ((iConfig.existsAs<bool>("doGroomedAK5")&& iConfig.getParameter< bool >("doGroomedAK5")) ?  new GroomedJetFiller("GroomedJetFiller", myTree, "AK5", "_PF", iConfig) : 0), 
   AK5groomedJetFiller_PFCHS ((iConfig.existsAs<bool>("doGroomedAK5")&& iConfig.getParameter< bool >("doGroomedAK5")) ?  new GroomedJetFiller("GroomedJetFiller", myTree, "AK5", "_PFCHS", iConfig) : 0), 
   genAK8groomedJetFiller ((iConfig.existsAs<bool>("doGroomedAK8")&& iConfig.getParameter< bool >("doGroomedAK8")) ?  new GroomedJetFiller("genGroomedJetFiller", myTree, "AK8", "_GEN", iConfig, 1) : 0), //Gen
   AK8groomedJetFiller_PF ((iConfig.existsAs<bool>("doGroomedAK8")&& iConfig.getParameter< bool >("doGroomedAK8")) ?  new GroomedJetFiller("GroomedJetFiller", myTree, "AK8", "_PF", iConfig) : 0),
   AK8groomedJetFiller_PFCHS ((iConfig.existsAs<bool>("doGroomedAK8")&& iConfig.getParameter< bool >("doGroomedAK8")) ?  new GroomedJetFiller("GroomedJetFiller", myTree, "AK8", "_PFCHS", iConfig) : 0),
   AK12groomedJetFiller_PF ((iConfig.existsAs<bool>("doGroomedAK12")&& iConfig.getParameter< bool >("doGroomedAK12")) ?  new GroomedJetFiller("GroomedJetFiller", myTree, "AK12", "_PF", iConfig) : 0), 
//PhotonFiller (  iConfig.existsAs<edm::InputTag>("srcPhoton") ?  new PhotonTreeFiller("PhotonFiller", myTree,  iConfig) : 0),
recoBosonFillerE( new VtoElectronTreeFiller( iConfig.getParameter<std::string>("VBosonType").c_str(), myTree, iConfig) ),
recoBosonFillerMu( new VtoMuonTreeFiller( iConfig.getParameter<std::string>("VBosonType").c_str(), myTree, iConfig) ), 
genBosonFiller( (iConfig.existsAs<bool>("runningOverMC") && iConfig.getParameter<bool>("runningOverMC")) ?  new MCTreeFiller(iConfig.getParameter<std::string>("VBosonType").c_str(), myTree, iConfig) : 0) 
*/
ewk::VplusJetsAnalysis::VplusJetsAnalysis(const edm::ParameterSet& iConfig) 
{

	myTree = fs -> mkdir("../").make<TTree>(iConfig.getParameter<std::string>("TreeName").c_str(),"V+jets Tree");
	if( iConfig.existsAs<edm::InputTag>("srcPFCor") )
	  CorrectedPFJetFiller.reset( new JetTreeFiller("CorrectedPFJetFiller", myTree, "PFCor", iConfig, 0) );
	//AK5
	if (iConfig.existsAs<bool>("doGroomedAK5")&& iConfig.getParameter< bool >("doGroomedAK5")){  
	  genAK5groomedJetFiller.reset(new GroomedJetFiller("genGroomedJetFiller", myTree, "AK5", "_GEN", iConfig, 1) );//Gen
	  AK5groomedJetFiller_PF.reset(new GroomedJetFiller("GroomedJetFiller", myTree, "AK5", "_PF", iConfig) );
	  AK5groomedJetFiller_PFCHS.reset(new GroomedJetFiller("GroomedJetFiller", myTree, "AK5", "_PFCHS", iConfig) ); 
	}
	//AK8
	if(iConfig.existsAs<bool>("doGroomedAK8")&& iConfig.getParameter< bool >("doGroomedAK8")){
	  genAK8groomedJetFiller.reset(new GroomedJetFiller("genGroomedJetFiller", myTree, "AK8", "_GEN", iConfig, 1) ); //Gen
	  AK8groomedJetFiller_PF.reset(new GroomedJetFiller("GroomedJetFiller", myTree, "AK8", "_PF", iConfig) );
	  AK8groomedJetFiller_PFCHS.reset(new GroomedJetFiller("GroomedJetFiller", myTree, "AK8", "_PFCHS", iConfig) );
	}
	//AK10
	if(iConfig.existsAs<bool>("doGroomedAK10")&& iConfig.getParameter< bool >("doGroomedAK10")){
	  genAK10groomedJetFiller.reset(new GroomedJetFiller("genGroomedJetFiller", myTree, "AK10", "_GEN", iConfig, 1) ); //Gen
	  AK10groomedJetFiller_PF.reset(new GroomedJetFiller("GroomedJetFiller", myTree, "AK10", "_PF", iConfig) );
	  AK10groomedJetFiller_PFCHS.reset(new GroomedJetFiller("GroomedJetFiller", myTree, "AK10", "_PFCHS", iConfig) );
	}
	//AK12
	if(iConfig.existsAs<bool>("doGroomedAK12")&& iConfig.getParameter< bool >("doGroomedAK12")){
	  genAK12groomedJetFiller.reset(new GroomedJetFiller("genGroomedJetFiller", myTree, "AK12", "_GEN", iConfig, 1) ); //Gen
	  AK12groomedJetFiller_PF.reset(new GroomedJetFiller("GroomedJetFiller", myTree, "AK12", "_PF", iConfig) );
	  AK12groomedJetFiller_PFCHS.reset(new GroomedJetFiller("GroomedJetFiller", myTree, "AK12", "_PFCHS", iConfig) );
	}
	if(iConfig.existsAs<bool>("runningOverMC") && iConfig.getParameter<bool>("runningOverMC"))
	  genBosonFiller.reset( new MCTreeFiller(iConfig.getParameter<std::string>("VBosonType").c_str(), myTree, iConfig) ); 
	recoBosonFillerE.reset( new VtoElectronTreeFiller( iConfig.getParameter<std::string>("VBosonType").c_str(), myTree, iConfig) );
	recoBosonFillerMu.reset( new VtoMuonTreeFiller( iConfig.getParameter<std::string>("VBosonType").c_str(), myTree, iConfig) );


	// final state; 0: ZJet; 1: Dijets
	if ( TString(iConfig.getParameter<std::string>("TreeName").c_str()).Contains("ZJet") ){ mFinalState=0;} 
	else if( TString(iConfig.getParameter<std::string>("TreeName").c_str()).Contains("Dijets") ){ mFinalState=1; }
	else{ cout<<"Don't support the final state: "<<iConfig.getParameter<std::string>("TreeName").c_str() <<endl; BREAK();}
	cout<<"Final State: "<<iConfig.getParameter<std::string>("TreeName").c_str() <<endl;

	// Are we running over Monte Carlo ?
	if( iConfig.existsAs<bool>("runningOverMC") ) runningOverMC_=iConfig.getParameter< bool >("runningOverMC");
	else runningOverMC_= false;

	if( mFinalState ==0){
		if(  iConfig.existsAs<edm::InputTag>("srcVectorBoson") ) mInputBoson = iConfig.getParameter<edm::InputTag>("srcVectorBoson");
		LeptonType_ = iConfig.getParameter<std::string>("LeptonType");
		VBosonType_ = iConfig.getParameter<std::string>("VBosonType");
	}

	if(  iConfig.existsAs<edm::InputTag>("srcPrimaryVertex") ) mPrimaryVertex = iConfig.getParameter<edm::InputTag>("srcPrimaryVertex"); 
	else mPrimaryVertex =  edm::InputTag("offlinePrimaryVertices");

	if(  iConfig.existsAs<edm::InputTag>("srcMet") ) mInputMet = iConfig.getParameter<edm::InputTag>("srcMet"); //if(  iConfig.existsAs<edm::InputTag>("srcMetMVA") ) mInputMetMVA = iConfig.getParameter<edm::InputTag>("srcMetMVA");

	//*********************  Run Over AOD or PAT  ***********//
	//if( iConfig.existsAs<bool>("runningOverAOD")) runoverAOD = iConfig.getParameter<bool>("runningOverAOD");

	//JetsFor_rho =  iConfig.getParameter<std::string>("srcJetsforRho") ; 

	//if(  iConfig.existsAs<edm::InputTag>("srcgenMet") ) mInputgenMet =  iConfig.getParameter<edm::InputTag>("srcgenMet") ; 

}



ewk::VplusJetsAnalysis::~VplusJetsAnalysis() {}


void ewk::VplusJetsAnalysis::beginJob() {
	declareTreeBranches(); // Declare all the branches of the tree
	std::cout<<"VplusJetsAnalysis begin work!"<<std::endl;
}



// ------------ method called to produce the data  ------------
void ewk::VplusJetsAnalysis::analyze(const edm::Event& iEvent, 
			const edm::EventSetup& iSetup) {
	// write event information: run, event, bunch crossing, ....
	run   = iEvent.id().run();
	event = iEvent.id().event();
	lumi  = iEvent.luminosityBlock();
	bunch = iEvent.bunchCrossing();
	//if( event != 5119780 && event!=5119805 && event!=5120829  ) return ;

	// primary/secondary vertices
	// edm::Handle<reco::VertexCollection > recVtxs;
	edm::Handle <edm::View<reco::Vertex> > recVtxs;
	iEvent.getByLabel( mPrimaryVertex, recVtxs);
	nPV = recVtxs->size();

	/////// PfMET information /////
	edm::Handle<edm::View<reco::MET> > pfmet;
	iEvent.getByLabel(mInputMet, pfmet);
	if (pfmet->size() == 0) {
		mpfMET   = -1;
		mpfSumET = -1;
		mpfMETSign = -1;
		mpfMETPhi   = -10.0;
	} else {
		//std::cout<<"pfmet size="<<pfmet->size()<<std::endl;
		mpfMET   = (*pfmet)[0].et();
		mpfSumET = (*pfmet)[0].sumEt();
		mpfMETSign = (*pfmet)[0].significance();
		mpfMETPhi   = (*pfmet)[0].phi();
	}

	/*/////// MVA MET information /////
	  edm::Handle<edm::View<reco::MET> > metMVA;
	  iEvent.getByLabel(mInputMetMVA, metMVA);
	  if (metMVA->size() == 0) {
	  mvaMET   = -1;
	  mvaSumET = -1;
	  mvaMETSign = -1;
	  mvaMETPhi   = -10.0;
	  } else {
	  mvaMET   = (*metMVA)[0].et();
	  mvaSumET = (*metMVA)[0].sumEt();
	  mvaMETSign = (*metMVA)[0].significance();
	  mvaMETPhi   = (*metMVA)[0].phi();
	  }*/

	/*/////// Pileup density "rho" in the event from fastJet pileup calculation /////
	  edm::Handle<double> rho;
	  const edm::InputTag eventrho(JetsFor_rho, "rho");
	  iEvent.getByLabel(eventrho,rho);
	  if( *rho == *rho) fastJetRho = *rho;
	  else  fastJetRho =  -999999.9;*/

	/////////// GenMET information & MC Pileup Summary Info  //////////
	mcPUtrueInteractions = 0;
	mcPUtotnvtx = 0;
	mcPUbx[0]   = -999; mcPUbx[1]   = -999; mcPUbx[2]   = -999;
	mcPUnvtx[0] = -999; mcPUnvtx[1] = -999; mcPUnvtx[2] = -999;
	if ( runningOverMC_ ){
		/*edm::Handle<reco::GenMETCollection> genMETs;
		  if(!runoverAOD){
		  iEvent.getByLabel(mInputgenMet,genMETs);
		  if ( genMETs->size() == 0) {
		  genMET   = -1.0;
		  genSumET = -1.0;
		  genMETSign  = -1.0;
		  genMETPhi   = -10.0;
		  } else {
		  genMET = (*genMETs)[0].et();
		  genSumET = (*genMETs)[0].sumEt();  
		  genMETSign = (*genMETs)[0].significance();  
		  genMETPhi = (*genMETs)[0].phi();
		  }
		  }*/
		// MC Pileup Summary Info
		const edm::InputTag PileupSrc("addPileupInfo");
		edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
		iEvent.getByLabel(PileupSrc, PupInfo);
		std::vector<PileupSummaryInfo>::const_iterator PVI;
		int ctid = 0;
		for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
			if (ctid>2) break;
			mcPUbx[ctid]   =  PVI->getBunchCrossing();
			mcPUnvtx[ctid] =  PVI->getPU_NumInteractions();
			mcPUtotnvtx   +=  PVI->getPU_NumInteractions();
			if(PVI->getBunchCrossing() == 0)mcPUtrueInteractions = PVI->getTrueNumInteractions();
			ctid++;
			//cout<<"bx="<<PVI->getBunchCrossing()<<" getTrueNumInteractions="<<PVI->getTrueNumInteractions()<<" getPU_NumInteractions()="<<PVI->getPU_NumInteractions()<<endl;
		}
	}

	if( mFinalState ==0){
		edm::Handle<edm::View< reco::Candidate> > boson;
		iEvent.getByLabel( mInputBoson, boson);
		mNVB = boson->size();
		if( mNVB<1 ) return; // Nothing to fill

		//if(GenJetFiller.get()) GenJetFiller->fill(iEvent);
		//if(PhotonFiller.get()) PhotonFiller->fill(iEvent);

		/**  Store reconstructed vector boson information */
		recoBosonFillerE->fill(iEvent, 0);
		// if(mNVB==2) recoBosonFillerE->fill(iEvent, 1);
		recoBosonFillerMu->fill(iEvent,0);
		// if(mNVB==2) recoBosonFillerMu->fill(iEvent, 1);
		if(genBosonFiller.get()) genBosonFiller->fill(iEvent); /**  Store generated vector boson information */
	}

	//Begin Jet substructure work
	// ---------------------------------------------------------------------
	// reading in different collections
	// ---------------------------------------------------------------------

	// pf inputs NoElectron
	edm::Handle< std::vector< reco::PFCandidate > >  pfs;
	iEvent.getByLabel("pfNoElectronPFlow", pfs);
	//std::cout << "pfs->size() = " << pfs->size() << std::endl;

	// PileUp
	edm::Handle< std::vector< reco::PFCandidate > >  pfs_PileUp;
	iEvent.getByLabel("pfPileUpPFlow", pfs_PileUp);
	//std::cout << "pfs_PileUp->size() = " << pfs_PileUp->size() << std::endl;

	edm::Handle< std::vector< reco::GenParticle > >  gens;
	iEvent.getByLabel("genParticles", gens);
	//std::cout << "gens->size() = " << gens->size() << std::endl;

	edm::Handle < reco::GenParticleRefVector > gens_nonu;
	iEvent.getByLabel("genParticlesForJetsNoNu", gens_nonu);
	//std::cout << "gens_nonu->size() = " << gens_nonu->size() << std::endl;

	//////// MC Simulation Event Weight ////////////////////////////
	evWeight = 1.0 ;
	edm::Handle< GenEventInfoProduct > genEventInfo;
	iEvent.getByLabel("generator", genEventInfo);

	if (genEventInfo.isValid()) {
		evWeight = genEventInfo->weight();
		//cout<<"evWeight="<<evWeight<<endl;
	}
	//////////////////////////////////////////////////////////////////

	// ---------------------------------------------------------------------
	// filling in pseudojets to be used as input for clustering
	// ---------------------------------------------------------------------

	/*std::vector<fastjet::PseudoJet> fjinputs_gens; fjinputs_gens.clear();
	  for (unsigned int i = 0; i < gens->size(); i++){
	  const reco::GenParticle P = (gens->at(i));
	  fastjet::PseudoJet tmp_psjet( P.px(), P.py(), P.pz(), P.energy() );
	// add user info about charge and pdgId
	tmp_psjet.set_user_info( new PseudoJetUserInfo(P.pdgId(), P.charge()) );
	print_p4(tmp_psjet,Form("gens_%i",i));
	fjinputs_gens.push_back( tmp_psjet );
	}*/


	std::vector<fastjet::PseudoJet> fjinputs_gens_nonu; fjinputs_gens_nonu.clear();
	for (unsigned int i = 0; i < gens_nonu->size(); i++){
		const reco::GenParticle& P = *(gens_nonu->at(i));
		fastjet::PseudoJet tmp_psjet( P.px(), P.py(), P.pz(), P.energy() );
		// add user info about charge and pdgId
		tmp_psjet.set_user_info( new PseudoJetUserInfo(P.pdgId(), P.charge()) );
		//print_p4(tmp_psjet,Form("gens_nonu_%i",i));
		fjinputs_gens_nonu.push_back( tmp_psjet );
	}	

	std::vector<fastjet::PseudoJet> fjinputs_pfs;  fjinputs_pfs.clear();
	std::vector<fastjet::PseudoJet> fjinputs_pfs_charge;  fjinputs_pfs_charge.clear();
	std::vector<fastjet::PseudoJet> fjinputs_pfs_neutral; fjinputs_pfs_neutral.clear();
	std::vector<fastjet::PseudoJet> fjinputs_pfs_noLep_noCHS; fjinputs_pfs_noLep_noCHS.clear();
	for (unsigned int i = 0; i < pfs->size(); i++){
		const reco::PFCandidate P = (pfs->at(i));
		fastjet::PseudoJet tmp_psjet( P.px(), P.py(), P.pz(), P.energy() );
		// add user info about charge and pdgId
		tmp_psjet.set_user_info( new PseudoJetUserInfo(P.pdgId(), P.charge()) );
		//print_p4(tmp_psjet,Form("pfs_%i",i),1);
		fjinputs_pfs.push_back( tmp_psjet );
		fjinputs_pfs_noLep_noCHS.push_back( tmp_psjet );
		if (P.charge()==0){
			fjinputs_pfs_neutral.push_back( tmp_psjet );
		}else{
			fjinputs_pfs_charge.push_back( tmp_psjet );
		}
	}

	std::vector<fastjet::PseudoJet> fjinputs_pfs_PileUp; fjinputs_pfs_PileUp.clear();
	for (unsigned int i = 0; i < pfs_PileUp->size(); i++){
		const reco::PFCandidate P = (pfs_PileUp->at(i));
		if (P.charge()==0){ std::cout<<"Error! there is a neutral particle in PileUp"<<std::endl; BREAK();}
		fastjet::PseudoJet tmp_psjet( P.px(), P.py(), P.pz(), P.energy() );
		// add user info about charge and pdgId
		tmp_psjet.set_user_info( new PseudoJetUserInfo(P.pdgId(), P.charge()) );
		//print_p4(tmp_psjet,Form("pfs_PileUp_%i",i));
		fjinputs_pfs_PileUp.push_back( tmp_psjet );
		fjinputs_pfs_noLep_noCHS.push_back( tmp_psjet );
	}
	fjinputs_pfs_noLep_noCHS = fastjet::sorted_by_pt(fjinputs_pfs_noLep_noCHS);
	// Now, we get 3 collections:
	// fjinputs_pfs_PileUp: PF charged, PileUp  
	// fjinputs_pfs_charge: PF charged, NoPileUp
	// fjinputs_pfs_neutral:PF neutral
	// fjinputs_pfs: PF neutral + charge, after CHS
	// fjinputs_pfs_noLep_noCHS: PF neutral + charge, no CHS




	/*
	// ---------------------------------------------------------------------
	// recluster on the fly
	// ---------------------------------------------------------------------
	double mJetRadius = 0.5;
	//fastjet::JetDefinition jetDef(fastjet::cambridge_algorithm, mJetRadius);
	fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, mJetRadius);
	std::cout<<"Clustered with "<< jetDef.description()<<std::endl;

	int activeAreaRepeats = 1;
	double ghostArea = 0.01;
	double ghostEtaMax = 5.0;

	fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
	fastjet::AreaDefinition fjAreaDefinition( fastjet::active_area, fjActiveArea );
	fastjet::AreaDefinition fjAreaDefinition_wGhosts( fastjet::active_area_explicit_ghosts, fjActiveArea );

	fastjet::ClusterSequenceArea* thisClustering_wGhosts_ = new fastjet::ClusterSequenceArea(fjinputs_gens_nonu, jetDef, fjAreaDefinition_wGhosts);
	fastjet::ClusterSequence* thisClustering_basic_ = new fastjet::ClusterSequence(fjinputs_gens_nonu, jetDef);

	std::vector<fastjet::PseudoJet> out_jets_wGhosts_ = sorted_by_pt(thisClustering_wGhosts_->inclusive_jets(15.0));    
	std::vector<fastjet::PseudoJet> out_jets_basic_ = sorted_by_pt(thisClustering_basic_->inclusive_jets(15.0));

	std::cout << "out_jets_wGhosts_->size() = " << out_jets_wGhosts_.size() << std::endl;
	std::cout << "out_jets_basic_->size() = " << out_jets_basic_.size() << std::endl;    

	// ---------------------------------------------------------------------
	// now pass on to a class which stores all the substructure information you would want...
	// inputs are jet definition, vector< PseudoJet >, area 
	//    - area, because you need to compute the area using all the constituents in an event
	//    - jet definition, because need to recluster the single jet to recover the cluster sequence
	// ---------------------------------------------------------------------

	//JetSubstructureTools *tmp = new JetSubstructureTools( jetDef, out_jets_basic_[0].constituents(), out_jets_wGhosts_[0].area() );

	for ( unsigned int i = 0; i < out_jets_basic_.size(); i++ ){

	std::cout << "out_jets_wGhosts_["<<i<<"].constituents().size() = " << out_jets_wGhosts_[i].constituents().size() << std::endl;
	std::cout << "out_jets_basic_["<<i<<"].constituents().size() = " << out_jets_basic_[i].constituents().size() << std::endl;        

	std::cout << "pt = " << out_jets_basic_[i].pt() << ", " << out_jets_wGhosts_[i].pt() << std::endl;

	JetSubstructureTools *tmp = new JetSubstructureTools( jetDef, out_jets_basic_[i].constituents(), out_jets_wGhosts_[i].area() );

	std::cout << "obj pt = " << tmp->getJet().pt() << ", jet charge = " << tmp->getJetCharge( 0.3 ) << std::endl;
	print_p4(tmp->getJet(),"fulljet");
	print_p4(tmp->getJet_wGhostes(),"fulljet_wGhost");
	print_p4(tmp->getJet_basic(),"basic jet");

	}*/


	bool doJetCleansing=1;
	//std::cout<<"====================== ak5 PF JET =========================="<<std::endl;
	if(genAK5groomedJetFiller.get() )genAK5groomedJetFiller->fill(iEvent, fjinputs_gens_nonu);
	if(AK5groomedJetFiller_PF.get() )AK5groomedJetFiller_PF->fill(iEvent, fjinputs_pfs_noLep_noCHS, doJetCleansing, fjinputs_pfs_charge,fjinputs_pfs_PileUp,fjinputs_pfs_neutral);//full event, bool, hard charge, pile up charge, full neutral
	if(AK5groomedJetFiller_PFCHS.get() )AK5groomedJetFiller_PFCHS->fill(iEvent, fjinputs_pfs);
	//std::cout<<"====================== ak8 PF JET =========================="<<std::endl;
	if(genAK8groomedJetFiller.get() )genAK8groomedJetFiller->fill(iEvent, fjinputs_gens_nonu);
	if(AK8groomedJetFiller_PF.get() )AK8groomedJetFiller_PF->fill(iEvent, fjinputs_pfs_noLep_noCHS, doJetCleansing, fjinputs_pfs_charge,fjinputs_pfs_PileUp,fjinputs_pfs_neutral);//full event, bool, hard charge, pile up charge, full neutral
	if(AK8groomedJetFiller_PFCHS.get() )AK8groomedJetFiller_PFCHS->fill(iEvent, fjinputs_pfs);
	//std::cout<<"====================== ak10PF JET =========================="<<std::endl;
	if(genAK10groomedJetFiller.get() )genAK10groomedJetFiller->fill(iEvent, fjinputs_gens_nonu);
	if(AK10groomedJetFiller_PF.get() )AK10groomedJetFiller_PF->fill(iEvent, fjinputs_pfs_noLep_noCHS, doJetCleansing, fjinputs_pfs_charge,fjinputs_pfs_PileUp,fjinputs_pfs_neutral);//full event, bool, hard charge, pile up charge, full neutral
	if(AK10groomedJetFiller_PFCHS.get() )AK10groomedJetFiller_PFCHS->fill(iEvent, fjinputs_pfs);
	//std::cout<<"====================== ak12PF JET =========================="<<std::endl;
	if(genAK12groomedJetFiller.get() )genAK12groomedJetFiller->fill(iEvent, fjinputs_gens_nonu);
	if(AK12groomedJetFiller_PF.get() )AK12groomedJetFiller_PF->fill(iEvent, fjinputs_pfs_noLep_noCHS, doJetCleansing, fjinputs_pfs_charge,fjinputs_pfs_PileUp,fjinputs_pfs_neutral);//full event, bool, hard charge, pile up charge, full neutral
	if(AK12groomedJetFiller_PFCHS.get() )AK12groomedJetFiller_PFCHS->fill(iEvent, fjinputs_pfs);

	//std::cout<<"====================== PAT::JET =========================="<<std::endl;
	if(CorrectedPFJetFiller.get()){ CorrectedPFJetFiller->fill(iEvent);}






	//BREAK();
	myTree->Fill();

} // analyze method






void ewk::VplusJetsAnalysis::endJob()
{
	std::cout<<"VplusJetsAnalysis end!"<<std::endl;
}




//  **** Utility: declare TTree branches for ntuple variables ***
void ewk::VplusJetsAnalysis::declareTreeBranches() {
	myTree->Branch("EventWeight", &evWeight, "EventWeight/F" );
	myTree->Branch("event_runNo",  &run,   "event_runNo/I");
	myTree->Branch("event_evtNo",  &event, "event_evtNo/I");
	myTree->Branch("event_lumi",   &lumi,  "event_lumi/I"); 
	myTree->Branch("event_bunch",  &bunch, "event_bunch/I"); 
	myTree->Branch("event_nPV",    &nPV,   "event_nPV/I"); 
	myTree->Branch("event_met_pfmet",             &mpfMET,      "event_met_pfmet/F"); 
	myTree->Branch("event_met_pfsumet",           &mpfSumET,    "event_met_pfsumet/F"); 
	myTree->Branch("event_met_pfmetsignificance", &mpfMETSign,  "event_met_pfmetsignificance/F"); 
	myTree->Branch("event_met_pfmetPhi",          &mpfMETPhi,   "event_met_pfmetPhi/F"); 
	myTree->Branch("event_metMVA_met",             &mvaMET,      "event_metMVA_met/F"); 
	myTree->Branch("event_metMVA_sumet",           &mvaSumET,    "event_metMVA_sumet/F"); 
	myTree->Branch("event_metMVA_metsignificance", &mvaMETSign,  "event_metMVA_metsignificance/F"); 
	myTree->Branch("event_metMVA_metPhi",          &mvaMETPhi,   "event_metMVA_metPhi/F"); 

	myTree->Branch("event_fastJetRho",                &fastJetRho,    "event_fastJetRho/F"); 

	if ( runningOverMC_ ){
		myTree->Branch("event_met_genmet",    &genMET,  "event_met_genmet/F"); 
		myTree->Branch("event_met_gensumet",  &genSumET,"event_met_gensumet/F"); 
		myTree->Branch("event_met_genmetsignificance", &genMETSign,  "event_met_genmetsignificance/F"); 
		myTree->Branch("event_met_genmetPhi",    &genMETPhi,  "event_met_genmetPhi/F"); 	  
		myTree->Branch("event_mcPU_totnvtx",    &mcPUtotnvtx,  "event_mcPU_totnvtx/F"); 
		myTree->Branch("event_mcPU_trueInteractions",    &mcPUtrueInteractions,  "event_mcPU_trueInteractions/F"); 
		myTree->Branch("event_mcPU_bx",         mcPUbx ,       "event_mcPU_bx[3]/F"); 
		myTree->Branch("event_mcPU_nvtx",       mcPUnvtx,      "event_mcPU_nvtx[3]/F"); 
	}

}  





// declare this class as a plugin
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
using ewk::VplusJetsAnalysis;
DEFINE_FWK_MODULE(VplusJetsAnalysis);
