#include "SelDiffZ/Selection/interface/MakeRootuplaForward.h" 
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HcalDetId/interface/HcalZDCDetId.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "DPGAnalysis/SiStripTools/interface/EventShape.h"
#include "DataFormats/HcalDetId/interface/HcalGenericDetId.h"
#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HcalDetId/interface/HcalZDCDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Candidate/interface/OverlapChecker.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

//PU reweight
#include "SelDiffZ/Selection/interface/Flat10.h"
#include "SelDiffZ/Selection/interface/ZSkim_v1.h"

MakeRootuplaForward::MakeRootuplaForward(const edm::ParameterSet& iConfig)
{
  electronCollectionTag_=iConfig.getUntrackedParameter<edm::InputTag>("electronCollectionTag");
  //  produces< pat::ElectronCollection > ("patElectrons").setBranchAlias("patElectrons");
  //std::string outputFile_D = iConfig.getUntrackedParameter<std::string>("filename");
  //outputFile_ = iConfig.getUntrackedParameter<std::string>("outputFile", outputFile_D);
  caloTowerTag_=iConfig.getParameter<edm::InputTag>("CaloTowerTag");
  if(electrons_)  zeeCollectionTag_ = iConfig.getUntrackedParameter<edm::InputTag>("zeeCollectionTag");
  zmumuCollectionTag_ = iConfig.getUntrackedParameter<edm::InputTag>("zmumuCollectionTag");
  ActivateMC_ = iConfig.getUntrackedParameter<Bool_t>("ActivateMC", true);
  useMB_ = iConfig.getUntrackedParameter<Bool_t>("useMinimumBias", false);
  JPsi_ = iConfig.getUntrackedParameter<Bool_t>("useJPsi", false);
  NoMassCuts_ = iConfig.getUntrackedParameter<Bool_t>("NoMassCuts", false);
  PVtxCollectionTag_      = iConfig.getParameter<edm::InputTag>("PVtxCollectionTag");
  TrackCollectionTag_      = iConfig.getParameter<edm::InputTag>("TrackCollectionTag");  
  VertexCollectionTag_      = iConfig.getParameter<edm::InputTag>("VertexCollectionTag"); 
  electrons_ = iConfig.getUntrackedParameter<Bool_t>("electrons");
  muons_ = iConfig.getUntrackedParameter<Bool_t>("muons");
  //  vertex_algo= iConfig.getParameter<std::string>("VERTEX_ALGO");
  //Flag: Vuoi riempire la entupla con gli elettroni

}

MakeRootuplaForward::~MakeRootuplaForward()
{ 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

bool
MakeRootuplaForward::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // *************************************************************************
  // General Run Info
  // *************************************************************************
  using namespace edm;
  using namespace std;
  using namespace pat;
  using namespace reco;
  using namespace edmNew;

  bool debug_deep=false;
  bool debug=false; 
  
  //if (NoMassCuts_) std::cout<<"ALERT: You are using the NoZCuts  analysis mode"<<std::endl; 
  //if (useMB_) std::cout<<"ALERT: You are using the MinimumBias analysis mode"<<std::endl; 
  //if (JPsi_) std::cout<<"ALERT: You are using the J/Psi analysis mode"<<std::endl; 
  
  double xi_Z_minus=0;
  double xi_Z_plus=0;
  edm::RunNumber_t runNumber   = iEvent.run();
  edm::EventNumber_t eventNumber = Long64_t( iEvent.eventAuxiliary().event() );
  edm::LuminosityBlockNumber_t lumiSection = iEvent.id().luminosityBlock();
  Rootuple->Zero();
  Rootuple->RunNumber=runNumber;
  Rootuple->EventNumber=eventNumber; 
  Rootuple->LumiSection=lumiSection; 

  if (debug_deep)
    std::cout<<"Event Number "<<eventNumber<<" Lumi "<<lumiSection<<endl;
  edm::Timestamp timestamp = iEvent.time();
  Rootuple->timestamp=timestamp;
  pat::Electron maxETelec2;
  pat::Electron maxETelec;
  
  
  //////////////////
  ////// Pile UP studies
  /////////////////
  
  if (ActivateMC_){
    try
      {
	edm::InputTag PileupSrc_ = (edm::InputTag) "addPileupInfo";
	Handle<std::vector< PileupSummaryInfo > >  PupInfo;
	iEvent.getByLabel(PileupSrc_, PupInfo);

  
	std::vector<PileupSummaryInfo>::const_iterator PVI;
	int npv = -1;
	
	for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
	  int BX = PVI->getBunchCrossing();
	  if(BX == 0) { 
	    npv = PVI->getPU_NumInteractions();
	    Rootuple->PU_NumInt       = PVI->getPU_NumInteractions() ;
	    Rootuple->PU_zpos         = PVI->getPU_zpositions()    ;
	    Rootuple->PU_ntrks_lowpT  = PVI->getPU_ntrks_lowpT()   ;
	    Rootuple->PU_ntrks_highpT = PVI->getPU_ntrks_highpT()  ; 
	    Rootuple->PU_sumpT_lowpT  = PVI->getPU_sumpT_lowpT()   ;  
	    Rootuple->PU_sumpT_highpT = PVI->getPU_sumpT_highpT()  ; 
	    std::cout << " Pileup Information: bunchXing, nvtx: " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << std::endl;
  
	    continue;
	  }      
	  
	  if (debug) std::cout << " Pileup Information: bunchXing, nvtx: " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << std::endl;
	}   
	
	double MyWeight = LumiWeights_.weight( npv );
	if (debug) cout<<"weight is "<<MyWeight<<endl;
	Rootuple->PUMCweight=MyWeight;
      }
    catch(cms::Exception& e)
      {
	std::cout<<"No PUinfoSummary block found"<<e.what();
      }
  }   
  else {
    Rootuple->PUMCweight=1;
  }
    


  if (electrons_){
    
    // *************************************************************************
    // TAKE ELECTRONS FROM EVENT and MAKE A PATTUPLA
    // *************************************************************************
    
    if (!useMB_){
      edm::Handle<reco::GsfElectronCollection> gsfElectrons;
      iEvent.getByLabel(electronCollectionTag_, gsfElectrons);
      if (!gsfElectrons.isValid()) {
	std::cout <<"MakeRootuplaForward: Could not get electron collection with label: "
		  <<electronCollectionTag_ << std::endl;
	//return false;
      }
      
      const reco::GsfElectronCollection *pElecs = gsfElectrons.product();
      reco::GsfElectronCollection::const_iterator elec;
      
      // *************************************************************************
      // GENERAL CUTS ON THE NUMBER OF ELECTRONS.
      // *************************************************************************
      
      //Do we have at least N electrons?
      int minNumbOfEle=2;
      const Int_t Nelecs = pElecs->size();
      if (Nelecs < minNumbOfEle) {
	if (debug_deep) std::cout << "No more than 1 electron found in this event" << std::endl;
	//return false; 
      } 
      
      // *************************************************************************
      // FILL AND CHECK ELECTRON PARAMETERS
      // *************************************************************************
         
      Int_t  counter = 0;
      std::vector<int> indices;
      std::vector<double> ETs;
      pat::ElectronCollection myElectrons;
      int count=0;
      double EtBasicCut=20;
      if (JPsi_) EtBasicCut=2;
      for (elec = pElecs->begin(); elec != pElecs->end(); ++elec) {
	Double_t sc_et=elec->caloEnergy()/TMath::CosH(elec->gsfTrack()->eta());
	indices.push_back(counter); ETs.push_back(sc_et);
	myElectrons.push_back(*elec);
	++counter;
      }
      const Int_t  event_elec_number = (int) indices.size();
      Int_t *sorted = new int[event_elec_number];
      Double_t *et = new double[event_elec_number];
      for (Int_t i=0; i<event_elec_number; i++) {
	et[i] = ETs[i];
	if (et[i]>EtBasicCut) count++; //Count the number of electrons over 20 GeV
      }
      // *************************************************************************
      // GENERAL CUTS ON THE ELECTRONS ENERGY
      // *************************************************************************
      Rootuple->numberOfLeptons=count;
      float ETCut=0;
      // array sorted now has the indices of the highest ET electrons
      TMath::Sort(event_elec_number, et, sorted, true);
      //
      // if the 2 highest electrons in the event has ET < ETCut_ return
      Int_t max_et_index = sorted[0];
      Int_t max_et_index2 = sorted[1];
      //Fill the electron properties (see ZeePlotsMatteo)
      if (ETs[max_et_index]<ETCut || ETs[max_et_index2]<ETCut) {
	if (debug) std::cout<<"The two most energetic electrons in the event have ET< ETcut"<<std::endl;
	//return false;
      } 
      delete [] sorted;
      delete [] et;
      
      // my electrons now:
      maxETelec = myElectrons[max_et_index];
      maxETelec2 = myElectrons[max_et_index2];   /// robi
      
      // *************************************************************************
      // FILL THE QUANTITY RELEATED TO THE ELECTRON
      // *************************************************************************
      
      for (elec = pElecs->begin(); elec != pElecs->end(); ++elec) {
	reco::GsfElectron mygsfelec = *elec;  
	//Rootuple->electrons.push_back(mygsfelec);  
	if (JPsi_ && mygsfelec.eEleClusterOverPout() <0.88) return false;
      }   
      //Rootuple->ElectronsWithEnoughEt=true;
      //if (JPsi_) Rootuple->maxETelec=maxETelec;
      //if (JPsi_) Rootuple->SecondETelec=maxETelec2;
      
      // *************************************************************************
      // Storing the parameters releated to the Zee
      // *************************************************************************
      
      //      if (!NoMassCuts_){
      //edm::Handle<pat::CompositeCandidateCollection> ZeeCands;
      //try{
      //  iEvent.getByLabel(zeeCollectionTag_, ZeeCands);
      //}
      //catch(cms::Exception& e)
      //  {
      //    std::cout<<"No Zee Composite Candidate found"<<e.what();
      //    return false;
      //  }
      //const pat::CompositeCandidateCollection *zcands = ZeeCands.product();
      //const pat::CompositeCandidateCollection::const_iterator zeeIter = zcands->begin();
      //const pat::CompositeCandidate zee = *zeeIter;
      //}
      //Rootuple->zee=zee;
      TLorentzVector e1;
      TLorentzVector e2;
      // Use directly the et,eta,phi from pat::Electron; assume e mass = 0.0
      e1.SetPtEtaPhiM(maxETelec.et(),maxETelec.eta(),maxETelec.phi(),0.0);
      e2.SetPtEtaPhiM(maxETelec2.et(),maxETelec2.eta(),maxETelec2.phi(),0.0);
      TLorentzVector Z = e1+e2;
      Double_t mee = Z.M();
      xi_Z_minus=(maxETelec.et()*pow(2.71,-maxETelec.eta()))/7000;
      xi_Z_minus+=(maxETelec2.et()*pow(2.71,-maxETelec2.eta()))/7000;
      xi_Z_plus=(maxETelec.et()*pow(2.71,maxETelec.eta()))/7000;
      xi_Z_plus+=(maxETelec2.et()*pow(2.71,maxETelec2.eta()))/7000;
      //xi_Z_minus=(Z.Et()*pow(2.71,-Z.Eta()))/7000;
      //xi_Z_plus=(Z.Et()*pow(2.71,Z.Eta()))/7000;
      //Double_t mee=zee.mass(); 
      double mee_down=60;
      double mee_up=120;
      //cout<<"mee "<<mee<<" and eta "<<Z.Eta()<<endl;
      
      Rootuple->ZMass=mee;  
      //Rootuple->etaZed=zee.eta();  
      Rootuple->etaZ=Z.Eta();   
      

      // *************************************************************************
      // Checking the other calorimetric parameters
      // *************************************************************************
    }  // if !useMB_
  } // if electrons_
 
  // *************************************************************************
  // Luminosity info:
  // ************************************************************************* 

  if (!ActivateMC_){
    edm::LuminosityBlock const& lumiBlock = iEvent.getLuminosityBlock();
    edm::Handle<LumiSummary> lumiSummary;
    edm::Handle<LumiDetails> d;
    lumiBlock.getByLabel("lumiProducer",d);
    lumiBlock.getByLabel("lumiProducer", lumiSummary);
    int bx=iEvent.bunchCrossing();
    const LumiSummary *lumi=lumiSummary.product();
    float istlumi=lumi->avgInsDelLumi();  // average (su LS) instant luminosity (delivered)
    // From 4_1
    const float istlumierr=d->lumiError(LumiDetails::kOCC1,bx);
    Rootuple->istlumi=istlumi;
    Rootuple->istlumierr=istlumierr;
    // From 4_1
    double  istlumiPerBX=d->lumiValue(LumiDetails::kOCC1,bx);
    Rootuple->istlumiPerBX=istlumiPerBX;
    Rootuple->bx=bx;
  }
  // *************************************************************************
  // Calo Towers
  // ************************************************************************* 


  /// New Thresholds 9 may 2011
  /// tuning should be done with unpaired Bunches (energy in HF)
  double energyThresholdHB = 1.25;
  double energyThresholdHE = 1.5; // 1.9;
  double energyThresholdHF = 6.0; // 4.0;
  double energyThresholdEB = 0.6;
  double energyThresholdEE = 1.5; //2.45;
 
  double etaMax=-999;
  double etaMin=999;
  double Epz_plus=0;  
  double Epz_minus=0;  

  int nTowersHF_plus = 0;
  int nTowersHF_minus = 0;
  int nTowersHE_plus = 0;
  int nTowersHE_minus = 0;
  int nTowersHB_plus = 0;
  int nTowersHB_minus = 0;
  int nTowersEE_plus = 0;
  int nTowersEE_minus = 0;
  int nTowersEB_plus = 0;
  int nTowersEB_minus = 0;
    
  //Sum(E)
  double sumEHF_plus = 0.;
  double sumEHF_minus = 0.;
  double sumEHE_plus = 0.;
  double sumEHE_minus = 0.;
  double sumEHB_plus = 0.;
  double sumEHB_minus = 0.;
  double sumEEE_plus = 0.;
  double sumEEE_minus = 0.;
  double sumEEB_plus = 0.;
  double sumEEB_minus = 0.;
  
  // Sum(ET)
  double sumETHF_plus = 0.;
  double sumETHF_minus = 0.;
  double sumETHE_plus = 0.;
  double sumETHE_minus = 0.;
  double sumETHB_plus = 0.;
  double sumETHB_minus = 0.;
  double sumETEB_plus = 0.;
  double sumETEB_minus = 0.;
  double sumETEE_plus = 0.;
  double sumETEE_minus = 0.;
  double xi_Calo_minus =0;
  double xi_Calo_plus =0;

  edm::Handle<CaloTowerCollection> towerCollectionH;
  iEvent.getByLabel(caloTowerTag_,towerCollectionH);
  const CaloTowerCollection& towerCollection = *towerCollectionH;

  CaloTowerCollection::const_iterator calotower;
  calotower = towerCollection.begin();
  CaloTowerCollection::const_iterator calotowers_end = towerCollection.end();
  //cout<< " before calotower loop "  <<endl;
  for(; calotower != calotowers_end; ++calotower) {
    
    if (fabs(calotower->eta())> 4.9) continue;
    
    bool hasHCAL = false;
    bool hasHF = false;
    bool hasHE = false;
    bool hasHB = false;
    bool hasHO = false;
    bool hasECAL = false;
    bool hasEE = false;
    bool hasEB = false;     
    for(size_t iconst = 0; iconst < calotower->constituentsSize(); iconst++){
      DetId adetId = calotower->constituent(iconst);
      if(adetId.det()==DetId::Hcal){
	hasHCAL = true;
	HcalDetId hcalDetId(adetId);
	if(hcalDetId.subdet()==HcalForward) hasHF = true;
	else if(hcalDetId.subdet()==HcalEndcap) hasHE = true;
	else if(hcalDetId.subdet()==HcalBarrel) hasHB = true;
	else if(hcalDetId.subdet()==HcalOuter) hasHO = true;  
      } 
      else if(adetId.det()==DetId::Ecal){
	hasECAL = true;
	EcalSubdetector ecalSubDet = (EcalSubdetector)adetId.subdetId();
	if(ecalSubDet == EcalEndcap) hasEE = true;
	else if(ecalSubDet == EcalBarrel) hasEB = true;
      }
    }
    int zside = calotower->zside();
    double caloTowerEnergy = calotower->energy();
    double caloTowerEmEnergy = calotower->emEnergy();
    double caloTowerHadEnergy = calotower->hadEnergy();
    double caloTowerPz = calotower->pz();

    // FIXME
    //double caloTowerET = calotower->et(primVtx.position());
    //double caloTowerET = calotower->et(primVtx.z());
    double caloTowerEt = calotower->et();
    double caloTowerEmEt = calotower->emEt();
    double caloTowerHadEt = calotower->hadEt();
    //double CaloTheta = calotower->p4(0.0).theta();
    //double pzck = calotower->p() * TMath::Cos(CaloTheta) ;
    ///// old wrong  Epz+= caloTowerEnergy + pt*zside;

    bool CalAboveTh = false;
    
    if( hasHF && !hasHE )
      {
	if( caloTowerEnergy > energyThresholdHF)
	  {
	    CalAboveTh = true;
	    //cout << "HF>> " <<  calotower->id() << " HAD energy "     << caloTowerHadEnergy << " EM energy " << caloTowerEmEnergy << " energy " << caloTowerEnergy << endl; 
	    if(zside >= 0)
	      {
		++nTowersHF_plus;
		sumEHF_plus += caloTowerEnergy;
		sumETHF_plus += caloTowerEt;
	      }
	    else
	      {
		++nTowersHF_minus;
		sumEHF_minus += caloTowerEnergy;
		sumETHF_minus += caloTowerEt;
	      }
	  }
      }
    else if( hasHE && !hasHF && !hasHB )
      {
	if( caloTowerHadEnergy > energyThresholdHE)
	  {
	    CalAboveTh = true;
	    //cout << "HE>> " <<  calotower->id() << "  HAD energy "     << caloTowerHadEnergy << " EM energy " << caloTowerEmEnergy << " energy " << caloTowerEnergy << endl;
	    if(zside >= 0)
	      {
		++nTowersHE_plus;
		sumEHE_plus += caloTowerHadEnergy;
		sumETHE_plus += caloTowerHadEt;
	      }
	    else
	      {
		++nTowersHE_minus;
		sumEHE_minus += caloTowerHadEnergy;
		sumETHE_minus += caloTowerHadEt;
	      }
	  }
      }
    else if( hasHB && !hasHE )
      {
	if( caloTowerHadEnergy > energyThresholdHB)
	  {
	    CalAboveTh = true;
	    //cout << "HB>> " <<  calotower->id() << "  HAD energy "     << caloTowerHadEnergy << " EM energy " << caloTowerEmEnergy << " energy " << caloTowerEnergy << endl;
	    if(zside >= 0)
	      {
		++nTowersHB_plus;
		sumEHB_plus += caloTowerHadEnergy;
		sumETHB_plus += caloTowerHadEt;
	      }
	    else
	      {
		++nTowersHB_minus;
		sumEHB_minus += caloTowerHadEnergy;
		sumETHB_minus += caloTowerHadEt;
	      }
	  }
      }
    
    if( hasEE && !hasEB )
      {
	if( caloTowerEmEnergy >= energyThresholdEE)
	  {
	    CalAboveTh = true;
	    //cout << "EE>> " <<  calotower->id() << "  HAD energy "     << caloTowerHadEnergy << " EM energy " << caloTowerEmEnergy << " energy " << caloTowerEnergy << endl;
	    if(zside >= 0)
	      {
		++nTowersEE_plus;
		sumEEE_plus += caloTowerEmEnergy;
		sumETEE_plus += caloTowerEmEt;
	      }
	    else
	      {
		++nTowersEE_minus;
		sumEEE_minus += caloTowerEmEnergy;
		sumETEE_minus += caloTowerEmEt;
	      }
	  }
      }
    else if( hasEB && !hasEE )
      {
	if( caloTowerEmEnergy >= energyThresholdEB)
	  {
	    CalAboveTh = true;
	    //cout << "EB>> " <<  calotower->id() << " HAD energy "     << caloTowerHadEnergy << " EM energy " << caloTowerEmEnergy << " energy " << caloTowerEnergy << endl; 
	    if(zside >= 0)
	      {
		++nTowersEB_plus;
		sumEEB_plus += caloTowerEmEnergy;
		sumETEB_plus += caloTowerEmEt;
	      }
	    else
	      {
		++nTowersEB_minus;
		sumEEB_minus += caloTowerEmEnergy;
		sumETEB_minus += caloTowerEmEt;
	      }
	  }
      }
    
    if(CalAboveTh)
      {
	if (calotower->eta() >= etaMax) etaMax=calotower->eta();
	if (calotower->eta() <= etaMin) etaMin=calotower->eta();
	xi_Calo_minus += caloTowerEt * pow(2.71,-calotower->eta()) / (7000);
	xi_Calo_plus += caloTowerEt * pow(2.71,calotower->eta()) / (7000);
	Epz_plus  += caloTowerEnergy + caloTowerPz;
	Epz_minus += caloTowerEnergy - caloTowerPz;

	//    cout << "NIK CALO " <<  calotower->id()  << " " <<  CaloTheta   << " - cos " << TMath::Cos(CaloTheta) ;
	//    cout << " -pz "     << caloTowerPz << " energy " << caloTowerEnergy << endl; 
	//    << " ck    " <<  pzck <<  endl;
	//    cout << " Epz_plus "     << Epz_plus << " Epz_minus   " <<  Epz_minus <<  endl;
      }
    
  }  ////has to close calotower loop

  
  //Fill the variables
  Rootuple->xi_Calo_minus=xi_Calo_minus;
  Rootuple->xi_Calo_plus=xi_Calo_plus;
  Rootuple->Epz_Calo_plus=Epz_plus;
  Rootuple->Epz_Calo_minus=Epz_minus;
  Rootuple->sumEHF_plus=sumEHF_plus;
  Rootuple->sumEHF_minus=sumEHF_minus;
  Rootuple->etaMax_Calo=etaMax;
  Rootuple->etaMin_Calo=etaMin;
  Rootuple->nTowersHF_plus=nTowersHF_plus;
  Rootuple->nTowersHF_minus=nTowersHF_minus;  
  if (sumEHF_plus>sumEHF_minus) Rootuple->minEHF=sumEHF_minus;
  if (sumEHF_plus<sumEHF_minus) Rootuple->minEHF=sumEHF_plus;

   edm::Handle<reco::VertexCollection> Vertexes;
  iEvent.getByLabel(VertexCollectionTag_, Vertexes); 
  Rootuple->numberOfVertexes = Vertexes->size(); 

  // *************************************************************************
  // Particle Flow Analysis
  // *************************************************************************
  
  double etaMin_PF,etaMax_PF;
  etaMin_PF=999;
  etaMax_PF=-999;
  double etaMin_PF_NOHF=999;
  double etaMax_PF_NOHF=-999;
  double etaMin_PF_Vertex_Selection=999;
  double etaMax_PF_Vertex_Selection=-999;
  edm::Handle <reco::PFCandidateCollection> PFCandidates;
  iEvent.getByLabel("particleFlow",PFCandidates); 
  reco::PFCandidateCollection::const_iterator iter;
  float etamin,etamax;
  etamin=etamax=0;
  float Epz_PF_plus=0;
  float Epz_PF_minus=0;
  float Epz_PF_plus_NOHF=0;
  float Epz_PF_minus_NOHF=0;
  double xi_PF_minus=0;
  double xi_PF_plus=0;
  double xi_PF_minus_charged_Vertex_Selection=0;
  double xi_PF_plus_charged_Vertex_Selection=0;  

  bool cold=false;
  bool MoreThanOneVertex=false;
  //std::vector<math::XYZPoint> vertexes;
  math::XYZPoint oldvertex;
  math::XYZPoint vertex;
  math::XYZPoint PrimaryVertex;
  double sumpx=0;
  double sumpy=0;
  double sumpz=0;
  double energyModule=0;
  double sumpxModule=0;
  double sumpyModule=0;
  double sumpzModule=0;
  double etaTimesEnergy=0;
  double energyTot_PF_Barrel_minus=0;
  double energyTot_PF_Barrel_plus=0;  
  double energyTot_PF_minus=0;
  double energyTot_PF_plus=0;
  std::vector<double> etas;
  int nPart_PF=0;
  double energyTot_PF_EE_minus=0;
  double energyTot_PF_EE_plus=0; 
  double sumEHF_minus_PF=0;
  double sumEHF_plus_PF=0;
  //double xi_PF_NOZED_minus=0;
  //double xi_PF_NOZED_plus=0;
  double xi_PF_NOHF_minus=0;
  double xi_PF_NOHF_plus=0; 
  TLorentzVector m1;
  TLorentzVector m2;
  bool first=false;
  bool second=false;
  double mmumu=0; 

  if (electrons_){
    // Detect the Z->ee vertex!
    for (reco::PFCandidateCollection::const_iterator iter = PFCandidates->begin(); iter != PFCandidates->end(); ++iter) {
      const reco::PFCandidate *particle = &(*iter);
      double et=particle->et();
      double charge=particle->charge();
      double phi=particle->phi();
      double eta=particle->eta();
      if (et>=20 && fabs(charge)>0){
	PrimaryVertex=particle->vertex();
	if (debug_deep) std::cout<<"The electron was in vertex "<<PrimaryVertex<<std::endl;
	if (charge<0) Rootuple->lepton2Phi=phi;
	if (charge<0) Rootuple->lepton2Eta=eta;
	if (charge>0) Rootuple->lepton1Phi=phi;
	if (charge>0) Rootuple->lepton1Eta=eta;
      }
    }
    
    //JPsi case
    if (JPsi_){
      for (reco::PFCandidateCollection::const_iterator iter = PFCandidates->begin(); iter != PFCandidates->end(); ++iter) {
	const reco::PFCandidate *particle = &(*iter);
	double et=particle->et();
	double charge=particle->charge();
	double phi=particle->phi();
	double eta=particle->eta();
	if (et>=4 && fabs(charge)>0){
	  PrimaryVertex=particle->vertex();
	  if (debug_deep) std::cout<<"The electron (JPsi case) was in vertex "<<PrimaryVertex<<std::endl;
	  if (charge<0) Rootuple->lepton2Phi=phi;
	  if (charge<0) Rootuple->lepton2Eta=eta;
	  if (charge>0) Rootuple->lepton1Phi=phi;
	  if (charge>0) Rootuple->lepton1Eta=eta;
	}
      }
    }
  }
  int count2=0;

  if (muons_){
    for (reco::PFCandidateCollection::const_iterator iter = PFCandidates->begin(); iter != PFCandidates->end(); ++iter) {
      const reco::PFCandidate *particle = &(*iter);
      double et=particle->et();
      double charge=particle->charge();
      double phi=particle->phi();
      double eta=particle->eta();
      double pt=particle->pt();
      if (et>=20 && (charge)>0){
	PrimaryVertex=particle->vertex();
	//std::cout<<"The muon was in vertex "<<PrimaryVertex<<std::endl;
	m1.SetPtEtaPhiM(pt,eta,phi,0.0);     
	//cout<<" the p4 is "<<particle->p4()<<endl;
	first=true;
	count2++;
	Rootuple->lepton1Phi=phi;
	Rootuple->lepton1Eta=eta;
        xi_Z_minus+=pt*pow(2.71,-eta)/7000;
	xi_Z_plus+=pt*pow(2.71,eta)/7000;	  
      }
      if (et>=20 && (charge)<0){
	PrimaryVertex=particle->vertex();
	//std::cout<<"The muon was in vertex "<<PrimaryVertex<<std::endl;
	m2.SetPtEtaPhiM(pt,eta,phi,0.0); 
	//cout<<" the p4 is "<<particle->p4()<<endl;
	second=true;
	count2++;
	Rootuple->lepton2Phi=phi;
	Rootuple->lepton2Eta=eta;
	xi_Z_minus+=pt*pow(2.71,-eta)/7000;
	xi_Z_plus+=pt*pow(2.71,eta)/7000;
      }
    }
    if (first==true && second==true){
      TLorentzVector Z = m1+m2;
      mmumu = Z.M();
      if (NoMassCuts_) Rootuple->ZMass=mmumu;     
      //cout<<"Invariant Mass is "<<mmumu<<endl;
      first=false;
      second=false;
    }
      
    if (!NoMassCuts_){
      edm::Handle<reco::CompositeCandidateCollection> ZmumuCands;
      iEvent.getByLabel(zmumuCollectionTag_, ZmumuCands);
      try{
	iEvent.getByLabel(zmumuCollectionTag_, ZmumuCands);
      }
      catch(cms::Exception& e)
	{
	  std::cout<<"No Zmumu Composite Candidate found"<<e.what();
	  return false;
	}
      const reco::CompositeCandidateCollection *zcands = ZmumuCands.product();
      const reco::CompositeCandidateCollection::const_iterator zmumuIter = zcands->begin();
      const reco::CompositeCandidate zmumu = *zmumuIter;
      //const pat::Muon * myElec1 = dynamic_cast<const pat::Muon*>( zmumu.daughter("muon1") );
      //const pat::Muon * myElec2 = dynamic_cast<const pat::Muon*>( zmumu.daughter("muon2") );
      mmumu=zmumu.mass(); 
      double eta=zmumu.eta();
      //xi_Z_minus=(zmumu.et()*pow(2.71,-eta))/7000;
      //xi_Z_plus=(zmumu.et()*pow(2.71,eta))/7000;
      double mmumu_down=60;
      double mmumu_up=120;
      if ((mmumu<mmumu_down || mmumu >mmumu_up) || NoMassCuts_){
	//std::cout<<"Zmumu invariant mass outside the limit "<<mmumu_down<<"-"<<mmumu_up<<" GeV"<<std::endl;
	return false;
      }
      Rootuple->ZMass=mmumu; 
      Rootuple->numberOfLeptons=count2;
      Rootuple->etaZ=eta;
    }
  }
  
  Rootuple->xi_Z_minus=xi_Z_minus;
  Rootuple->xi_Z_plus=xi_Z_plus;

  TLorentzVector dataMass(0.,0.,0.,0.);
  std::vector<double> electronEnergy;
  std::vector<double> muEnergy;

  /// PFlow loop 
  double PtThPFCharged = 0.1;  // it was 0.15
  double EnThPFBar = 1.5;
  double EnThPFEnd = 3.5; // 4.;
  double EnThPFFw  = 6.0; // 6; 

  for (reco::PFCandidateCollection::const_iterator iter = PFCandidates->begin(); iter != PFCandidates->end(); ++iter) {
    float MaxAllowableDistance=0.1; //Maximum difference between two vertixes
    const reco::PFCandidate *particle = &(*iter);
    double et=particle->et();
    double energy=particle->energy();
    double pt=particle->pt();
    double p=particle->p();
    double px=particle->px();
    double py=particle->py();
    double pz=particle->pz();
    double eta=particle->eta();
    //double phi=particle->phi();
    double charge=particle->charge();
    double theta=particle->theta();

    //eta cut - excluding last HF module 
    if (fabs(eta)>4.9) continue;

    if (cold==true) 
      {
	oldvertex=vertex;
      }
    vertex=particle->vertex();
    int type=particle->particleId();
    if (particle->particleId()==reco::PFCandidate::e) electronEnergy.push_back(et);
    if (particle->particleId()==reco::PFCandidate::mu) muEnergy.push_back(et);
    if (cold==true) {
      float dx=oldvertex.x()-vertex.x();
      float dy=oldvertex.y()-vertex.y();
      float dz=oldvertex.z()-vertex.z();
      if (fabs(dx) > 0.1 || fabs(dz) > 0.2 || fabs(dy) > 0.1  ){
	MoreThanOneVertex=true;
      }
    }
    cold=true;
    
    TLorentzVector tmp(px,py,pz,energy); 
    //Calculate the distance between the primary vertex and the relative vertex
    double distance= pow( pow(PrimaryVertex.x()-vertex.x(),2)+  pow(PrimaryVertex.y()-vertex.y(),2) +  pow(PrimaryVertex.z()-vertex.z(),2) , 0.5);

    //Rootuple->distance.push_back(distance);  
    //Rootuple->charge.push_back(charge);

    if (debug_deep)std::cout<<"The particle type "<<type<<" has eta "<<eta<<" and pt "<<pt<<" and p "<<p<<" and energy "<<energy<<" and charge "<<charge<<" and theta "<<theta << endl;
    //if (debug_deep) std::cout<<"The particle type "<<type<<" has eta "<<eta<<" and pt "<<pt<<" and energy "<<energy<<" and charge "<<charge<<" having vertex "<<vertex<<" (far from P.V "<<distance<<")"<<std::endl;


    if  (  (fabs(charge) >0 && pt >  PtThPFCharged ) ||
	   (fabs(charge) == 0  && ( (fabs(eta) <= 1.5 && energy > EnThPFBar)  ||
				     (fabs(eta) > 1.5 && fabs(eta) <= 3 && energy > EnThPFEnd) ||
				     (fabs(eta) > 3 && energy >EnThPFFw) ) )   )
      {        
	nPart_PF++;
	Epz_PF_plus+=p+p*TMath::Cos(theta);
	Epz_PF_minus+=p-p*TMath::Cos(theta);
	xi_PF_minus += et * pow(2.71,-eta) / (7000);
	xi_PF_plus += et * pow(2.71,eta) / (7000);

	etaTimesEnergy+=eta*energy;
	sumpxModule +=fabs(px);
	sumpyModule +=fabs(py);
	sumpzModule +=fabs(pz);
	sumpx +=px;
	sumpy +=py;
	sumpz +=pz;
	energyModule +=energy;

	if (fabs(charge) >0 && distance < MaxAllowableDistance ) {
	  if (eta<0) xi_PF_minus_charged_Vertex_Selection += et * pow(2.71,-eta) / (7000);
	  if (eta>0) xi_PF_plus_charged_Vertex_Selection += et * pow(2.71,eta) / (7000);
	  if (eta>etaMax_PF_Vertex_Selection ) etaMax_PF_Vertex_Selection=eta;
	  if (eta<etaMin_PF_Vertex_Selection ) etaMin_PF_Vertex_Selection=eta;
	    //// MaxAllow = 0.1 
	}
	
	if (eta>0) energyTot_PF_plus+=energy;
	if (eta<0) energyTot_PF_minus+=energy;
	if (fabs(eta)<3)
	  {
	    xi_PF_NOHF_minus += et * pow(2.71,-eta) / (7000);
	    xi_PF_NOHF_plus += et * pow(2.71,eta) / (7000);
	    Epz_PF_plus_NOHF+=p+p*TMath::Cos(theta);
	    Epz_PF_minus_NOHF+=p-p*TMath::Cos(theta);
	    if (eta>-3. && eta <-1.5 ) energyTot_PF_EE_minus+=energy;
	    else if (eta>=-1.5 && eta < 0) energyTot_PF_Barrel_minus+=energy;
	    else if (eta>0 && eta <=1.5) energyTot_PF_Barrel_plus+=energy;
	    else  energyTot_PF_EE_plus+=energy;
	  }
	
	if (eta>etaMax_PF) etaMax_PF=eta;
	if (eta<etaMin_PF) etaMin_PF=eta;
	if (eta>etaMax_PF_NOHF && eta<3) etaMax_PF_NOHF=eta;
	if (eta<etaMin_PF_NOHF && eta>-3) etaMin_PF_NOHF=eta;

	if (fabs(eta) > 3 ) {
	  if (eta<0) sumEHF_minus_PF += energy;
	  if (eta>0) sumEHF_plus_PF += energy;
	}
	etas.push_back(eta);
	dataMass+=tmp;
	HistoEtaEnergyW->Fill(eta,energy);
      }
  }  // PF loop



  Rootuple->Mx2=dataMass.M2();  /// massaquadro misurata
  Rootuple->P_x=dataMass.X();
  Rootuple->P_y=dataMass.Y();
  Rootuple->P_z=dataMass.Z();
  //  float * pippo =  HistoEtaEnergyW->GetArray();
  //cout << HistoEtaEnergyW->GetNbinsX() << endl;
  for (int i=1;i<=HistoEtaEnergyW->GetNbinsX();i++) {
    //cout << "ARRAY " <<  *(pippo+i) << "  " << HistoEtaEnergyW->GetBinContent(i) << endl;
    Rootuple->EnergyInEta.push_back(HistoEtaEnergyW->GetBinContent(i)) ;
  }

 
  

  if (debug_deep)    std::cout<<"DATA has "<<sumEHF_minus_PF<<" "<<sumEHF_plus_PF<<std::endl;  
  //if (debug_deep) std::cout<<"Eta Min PF is "<<etaMin_PF<<" while max is "<<etaMax_PF<<" and Eta Min No Vertex "<<etaMin_PF_Vertex_Selection<<" and Eta Max No Vertes "<<etaMax_PF_Vertex_Selection<<" and csi_PF_minus is "<<xi_PF_minus<<endl; 

 
  //// Computing GAPs
  //// adding two fake entries at +-4.9 in etas!!!

  etas.push_back(4.9);
  etas.push_back(-4.9);

  const int  size = (int) etas.size();
  int *sorted = new int[size];
  double *v = new double[size];
  double eta_gap_limplus = -10.0;
  double eta_gap_limminus = -10.0;
  
  for (int i=0; i<size; i++) {
    v[i] = etas[i];
    if (debug_deep) cout<<v[i]<<endl;
  }
  TMath::Sort(size, v, sorted, true);

  if (size > 1) {
    double *diff = new double[size-1];
    int *diffsorted = new int[size-1];
    for (int i=0; i<(size-1); i++) {
      diff[i] = fabs(etas[sorted[i+1]]-etas[sorted[i]]);
      //if (debug_deep) cout<<" etas " << i << " size " << size << " diff "<< diff[i]<<endl;
      //cout<<" etas "  << " = " << etas[sorted[i+1]] << " - " <<  etas[sorted[i]] <<  " diff "<< diff[i]<<endl;
    }
  
    TMath::Sort(size-1, diff, diffsorted, true);
    
    //checking the max gap
    double max_eta_gap=diff[diffsorted[0]];
    eta_gap_limminus = etas[sorted[diffsorted[0]+1]] ;
    eta_gap_limplus = etas[sorted[diffsorted[0]]] ;

    /// cout << "SUP " <<  eta_gap_limplus  << " " <<  eta_gap_limminus  << endl;

    Rootuple->max_eta_gap_PF=max_eta_gap;

    if (size>2) {
      double max_second_eta_gap=diff[diffsorted[1]];
      Rootuple->max_second_eta_gap_PF=max_second_eta_gap;
      if (debug_deep) cout<<" diff  " << diff[diffsorted[0]] << " sec " << diff[diffsorted[1]] << " diff size "<< diff[size-2] <<endl;
    }
    
    delete [] diff;
    delete [] diffsorted;
  }
  
  if (debug_deep) cout<<" GAPs  " << Rootuple->max_eta_gap_PF << "  " <<  Rootuple->max_second_eta_gap_PF <<endl;
  
  delete [] sorted;
  delete [] v;

  //sorting electron energy
  const int  size3 = (int) electronEnergy.size();
  int *sorted3 = new int[size3];
  double *v3 = new double[size3];
  
  for (int i=0; i<size3; i++) {
    v3[i] = electronEnergy[i];
  }
  TMath::Sort(size3, v3, sorted3, true);
  for (int i=0; i<size3; i++) {
    electronEnergy[i] = v3[sorted3[i]];
  }

  //sorting muon energy
  const int  size4 = (int) muEnergy.size();
  int *sorted4 = new int[size4];
  double *v4 = new double[size4];
  
  for (int i=0; i<size4; i++) {
    v4[i] = muEnergy[i];
  }
  TMath::Sort(size4, v4, sorted4, true);
  for (int i=0; i<size4; i++) {
    muEnergy[i] = v4[sorted4[i]];
  }
  delete [] sorted3;
  delete [] v3;
  delete [] sorted4;
  delete [] v4;

  /// Loop to compute Mx2 a destra e a sinistra del GAP
 
  TLorentzVector dataMass_plus(0.,0.,0.,0.);
  TLorentzVector dataMass_minus(0.,0.,0.,0.);
  int nplus =0;
  int nminus =0;

  for (reco::PFCandidateCollection::const_iterator iter = PFCandidates->begin(); iter != PFCandidates->end(); ++iter) {
    const reco::PFCandidate *particle = &(*iter);
    double energy=particle->energy();
    double pt=particle->pt();
    double px=particle->px();
    double py=particle->py();
    double pz=particle->pz();
    double eta=particle->eta();
    double charge=particle->charge();
    
    //eta cut - excluding last HF module 
    if (fabs(eta)>4.9) continue;
    
    TLorentzVector tmp(px,py,pz,energy); 
    
    if  (  (fabs(charge) >0 && pt >  PtThPFCharged ) ||
	   (fabs(charge) == 0  && ( (fabs(eta) <= 1.5 && energy > EnThPFBar)  ||
				    (fabs(eta) > 1.5 && fabs(eta) <= 3 && energy > EnThPFEnd) ||
				    (fabs(eta) > 3 && energy >EnThPFFw) ) )   )
      {        
	
	if ( eta >= eta_gap_limplus ) {
	  dataMass_plus+=tmp;
	  nplus++;
	}
	else {
	  dataMass_minus+=tmp;
	  nminus++;
	}
      }
  }  // PF loop


  Rootuple->Mx2_plus=dataMass_plus.M2();  /// massaquadro misurata
  Rootuple->Mx2_minus=dataMass_minus.M2();  /// massaquadro misurata
  Rootuple->N_mx2plus=nplus;  /// massaquadro misurata
  Rootuple->N_mx2minus=nminus;  /// massaquadro misurata
  Rootuple->eta_gap_limplus=eta_gap_limplus;  /// massaquadro misurata

  //  cout << "Mass2 "  << nPart_PF << " " << Rootuple->Mx2 << " " <<  nplus << " " << Rootuple->Mx2_plus << " " << nminus << " " <<  Rootuple->Mx2_minus << " " << "  eta  " <<  eta_gap_limplus <<   endl;
  Rootuple->electronEnergy=electronEnergy;
  Rootuple->muEnergy=muEnergy;


  if (debug_deep) std::cout<<"sumpxModule "<<sumpxModule<<" sumpyModule "<<sumpyModule<<" sumpzModule "<<sumpzModule<<" sumpx "<<sumpx<<" sumpy "<<sumpy<<" sumpz "<<sumpz<<std::endl;
  math::XYZVector ResultingAllTracks(sumpx,sumpy,sumpz);
  double etaAllTracks=ResultingAllTracks.eta();
  if (debug_deep) std::cout<<"etaAllTracks "<<etaAllTracks<<endl;
  Rootuple->etaAllTracks_PF=etaAllTracks;
  Rootuple->energyTot_PF=energyModule;
  double etaWeightedOnEnergy=etaTimesEnergy/energyModule;
  Rootuple->etaWeightedOnEnergy_PF=etaWeightedOnEnergy;
  Rootuple->Epz_PF_plus=Epz_PF_plus;
  Rootuple->Epz_PF_minus=Epz_PF_minus;
  Rootuple->etaMax_PF=etaMax_PF;
  Rootuple->etaMin_PF=etaMin_PF;
  Rootuple->Epz_NOHF_PF_plus=Epz_PF_plus_NOHF;
  Rootuple->Epz_NOHF_PF_minus=Epz_PF_minus_NOHF;
  Rootuple->etaMax_NOHF_PF=etaMax_PF_NOHF;
  Rootuple->etaMin_NOHF_PF=etaMin_PF_NOHF;
  Rootuple->etaMax_Charged_PV_PF=etaMax_PF_Vertex_Selection;
  Rootuple->etaMin_Charged_PV_PF=etaMin_PF_Vertex_Selection;
  Rootuple->energyTot_PF_minus=energyTot_PF_minus;
  Rootuple->energyTot_PF_plus=energyTot_PF_plus;
  Rootuple->xi_PF_minus=xi_PF_minus;  
  Rootuple->xi_PF_plus=xi_PF_plus; 
  Rootuple->xi_Z_minus=xi_Z_minus;  
  Rootuple->xi_Z_plus=xi_Z_plus; 
  Rootuple->xi_PF_NOHF_minus=xi_PF_NOHF_minus;  
  Rootuple->xi_PF_NOHF_plus=xi_PF_NOHF_plus; 
  Rootuple-> xi_PV_PF_charged_minus=xi_PF_minus_charged_Vertex_Selection;
  Rootuple-> xi_PV_PF_charged_plus=xi_PF_plus_charged_Vertex_Selection;
  Rootuple-> nPart_PF=nPart_PF;
  Rootuple-> energyTot_PF_EE_minus=energyTot_PF_EE_minus;
  Rootuple-> energyTot_PF_EE_plus=energyTot_PF_EE_plus;
  Rootuple-> energyTot_PF_Barrel_minus= energyTot_PF_Barrel_minus;
  Rootuple-> energyTot_PF_Barrel_plus= energyTot_PF_Barrel_plus; 
  Rootuple-> sumEHF_PF_minus=sumEHF_minus_PF;
  Rootuple-> sumEHF_PF_plus=sumEHF_plus_PF;



  // *************************************************************************
  // Montecarlo Quantities
  // ************************************************************************* 

  if (ActivateMC_){
    if (debug_deep) std::cout<<"You activate the MC variables analysis: getting the MC truth"<<std::endl;

    //Variable declaration
    int count=0;
    TLorentzVector part(0.,0.,0.,0.);
    TLorentzVector partVis(0.,0.,0.,0.);
    TLorentzVector partZ(0.,0.,0.,0.);
    double sumECastor_minus_gen=0;
    double sumECastor_plus_gen=0;
    double sumEZDC_minus_gen=0;
    double sumEZDC_plus_gen=0;
    double etaOutcomingProton=0;
    double energyOutcomingProton=0;
    double mostEnergeticXL=0;
    double mostEnergeticXLNum=0;
    vector<double> eta_gen_vec;
    double xi_Z_gen_minus=0;
    double xi_Z_gen_plus=0;
    double etaZ_gen=0;
    double energyZ_gen=0;
    double p_diss_mass=0;
    double p_diss=0;
    double xL_p_diss=0;

    double xL_etaGTP5=0;
    double xL_etaLTM5=0;
    int xL_LTM5Num=0;
    int xL_GTP5Num=0;

    std::vector<double> genpt;
    std::vector<double> tracks;
    std::vector<double> eta;
    std::vector<double> etaPT;
    std::vector<double> tracksPT;

    try
      {
	//////OLD BLOCK
	const HepMC::GenEvent *myGenEvent;
	// 3) generated electrons
	Handle<edm::HepMCProduct> hepMC;
	iEvent.getByLabel("generator",hepMC);
	myGenEvent = hepMC->GetEvent();
	//Loop on gen particles
	HepMC::GenEvent::particle_const_iterator mcIter;
	
	//edm::Handle<reco::GenParticleCollection> genParticles;
	//iEvent.getByLabel("generator", genParticles );
	
	//for (reco::GenParticleCollection::const_iterator mcIter=genParticles->begin();mcIter!=genParticles->end();++mcIter){
	
	for ( mcIter=myGenEvent->particles_begin(); mcIter != myGenEvent->particles_end(); mcIter++ ) {
	  bool electronFromZ=false;
	  int motherId=0;
	  int motherstatus=0;
	  int status=(*mcIter)->status(); 
	  int pdg=(*mcIter)->pdg_id();
	  HepMC::GenParticle* p = *(mcIter);
	  HepMC::GenVertex* productionVertex = p->production_vertex();
	  double part_pt = sqrt( p->momentum().px()* p->momentum().px()+ p->momentum().py()* p->momentum().py());
	  
	  //if (debug_deep) cout<<"The particle is "<<(*mcIter)->pdg_id()<<" in status "<<status<<endl;
	  if (fabs(pdg)==11 || fabs(pdg)==22){
	    
	    try{ 
	      motherId = (*(productionVertex->particles_in_const_begin()))->pdg_id();
	      motherstatus=(*(productionVertex->particles_in_const_begin()))->status();
	    }
	    catch(cms::Exception& e) {
	      std::cout << " Not possible to access to the motherId... !" << std::endl;
	    }
	    if (motherId==23) {
	      electronFromZ=true;
	      //cout<<"This particle (id "<<pdg<<" ) comes from the Z...the status is "<<status<<endl;
	      TLorentzVector tmp((p->momentum()).x(),(p->momentum()).y(),(p->momentum()).z(),(p->momentum()).e() );
	      partZ+=tmp;
	    }
	    if (fabs(motherId)==11) {
	      //cout<<"This particle (id "<<pdg<<" ) comes from an electron...the particle status is "<<status<<", while its mother's status is "<<motherstatus<<endl;
	    }
	  }
	  double eta_gen= p->momentum().eta();
	  //if (debug_deep) cout<<"While the mother is "<<motherId<<endl;
	  if (count==2) {
	    // cout <<  " 3rd MC " << (*mcIter)->pdg_id() << endl; 
	    p_diss_mass= p->momentum().m();
	    p_diss= p->momentum().z();
	    if ( pdg == 2212){
	      etaOutcomingProton= p->momentum().eta();
	      energyOutcomingProton= p->momentum().e();
	    }
	  }
	  if (( pdg == 23)){
	    xi_Z_gen_minus=( p->momentum().e() - p->momentum().z() )/7000;
	    xi_Z_gen_plus=( p->momentum().e() + p->momentum().z() )/7000;
	    etaZ_gen=eta_gen;
	    energyZ_gen= p->momentum().e();
	    //cout<<"Z generated main parameters: eta "<<etaZ_gen<<" energy "<<energyZ_gen<<endl;
	  }
	  
	  eta_gen_vec.push_back( p->momentum().eta());
	  if ( status == 1 && count>2) {   
	    TLorentzVector tmp(( p->momentum()).x(),( p->momentum()).y(),( p->momentum()).z(),( p->momentum()).e() );
	    part+=tmp;
	  }
	  if ( (status == 1) &&  (fabs(eta_gen) < 4.9) && (part_pt > 0.10) ) {   // if particle has a chance to reach the detector ...
	    TLorentzVector tmp(( p->momentum()).x(),( p->momentum()).y(),( p->momentum()).z(),( p->momentum()).e() );
	    partVis+=tmp;
	    //cout << " nel loop di Mx2_gen " <<  endl;
	  }
	  
	  //if (debug_deep) std::cout<<"Particle momentum along z is "<<setprecision(12)<<( p->momentum()).z()<<" and status "<<status<<std::endl;      
	  //if (debug_deep) cout<<"eta_gen "<<eta_gen<<" and status "<<status<<" pdg_id "<< p->pdg_id()<<" and energy "<< p->momentum().e()<<endl;
	  
	  
	  // new xL_gen definition (after Sasha)
	  if (count>=2 && status ==1)
	    {
	      if (eta_gen > 4.9)  
		{
		  xL_etaGTP5 += p->momentum().z();
		  xL_GTP5Num++;
		}
	      if (eta_gen < -4.9)  
		{
		  xL_etaLTM5 += p->momentum().z();
		  xL_LTM5Num++;
		}
	    }
	  
	  
	  if (count>=2 && status==1){
	    if (p_diss>0) {
	      if ( xL_p_diss < p->momentum().z() ){
		xL_p_diss= p->momentum().z();
	      }
	    }
	    if (p_diss<0) {
	      if ( xL_p_diss > p->momentum().z() ){
		xL_p_diss= p->momentum().z();
	      }
	    }
	  }
	  if ( fabs(eta_gen)>5.2 && fabs(eta_gen)<6.6 && status == 1){
	    //if (debug_deep) std::cout<<"Particle in Castor, having eta "<<eta_gen<<" and energy "<< p->momentum().e()<<endl;
	    if (eta_gen<0) sumECastor_minus_gen += p->momentum().e();
	    if (eta_gen>0) sumECastor_plus_gen += p->momentum().e();
	  }
	  
	  if ( fabs(eta_gen)>8.2 && status == 1 && ( pdg == 2112 || pdg == 22) ){
	    //if (debug_deep) std::cout<<"Particle in ZDC, having eta "<<eta_gen<<" and energy "<< p->momentum().e()<<endl;
	    if (eta_gen<0) sumEZDC_minus_gen += p->momentum().e();
	    if (eta_gen>0) sumEZDC_plus_gen += p->momentum().e();
	  }      
	  count++;
	} // loop over particles
	
	float Mx2_gen=partVis.M2(); /// massaquadro visibile generata
	TLorentzVector NOZ=partVis-partZ;
	float Mx2_NOZ_gen=NOZ.M2();
	//std::cout<<"Particle XL "<< mostEnergeticXL << " id "<< mostEnergeticXLType <<endl;
	if (debug_deep) cout<<"Mx2_gen is "<<Mx2_gen<<" while eta of the outcoming proton is "<<etaOutcomingProton<<" and the energy "<<energyOutcomingProton<<endl;
	
	mostEnergeticXL = xL_etaGTP5/3500.;
	mostEnergeticXLNum = xL_GTP5Num ;
	if (fabs(xL_etaGTP5)<fabs(xL_etaLTM5)) 
	  {
	    mostEnergeticXL = xL_etaLTM5/3500.;
	    mostEnergeticXLNum = xL_LTM5Num ;
	  }
	
	// cout << "* XLgen " << mostEnergeticXL << " num " << mostEnergeticXLNum << " + " << xL_etaGTP5 << " - " << xL_etaLTM5 <<  endl;
	
	Rootuple-> nParticles_gen= count;
	Rootuple-> Mx2_gen= Mx2_gen;
	Rootuple-> Mx2_NOZ_gen= Mx2_NOZ_gen;
	Rootuple-> sumECastor_gen_minus=sumECastor_minus_gen;
	Rootuple-> sumECastor_gen_plus=sumECastor_plus_gen;
	Rootuple-> sumEZDC_gen_minus=sumEZDC_minus_gen;
	Rootuple-> sumEZDC_gen_plus=sumEZDC_plus_gen;
	Rootuple-> etaOutcomingProton=etaOutcomingProton;
	Rootuple-> xL_gen=mostEnergeticXL;
	Rootuple-> xL_Num_gen=mostEnergeticXLNum;
	Rootuple-> xi_Z_gen_minus=xi_Z_gen_minus;
	Rootuple-> xi_Z_gen_plus=xi_Z_gen_plus;
	Rootuple-> etaZ_gen=etaZ_gen;
	Rootuple-> energyZ_gen=energyZ_gen;
	Rootuple-> p_diss_mass_gen=p_diss_mass;
	Rootuple-> xL_p_diss= xL_p_diss;
	
	//cout<<"Mx2_gen "<<Rootuple->Mx2_gen<<" Mx2 "<<Rootuple->Mx2 <<endl;	
	
      }
    catch(cms::Exception& e)
      {
	std::cout<<"No HEPMCProduct block found"<<e.what();
      }
    
    

    Handle<reco::GenParticleCollection> genParticles;     
    iEvent.getByLabel("genParticles",genParticles);  // standard PYTHIA collection
    //iEvent.getByLabel("genParticlePlusGEANT",genParticles);  // PYTHIA + GEANT collection
    int numseltracks =0;

    //cout << "XLGEN " << Rootuple->xL_gen  << endl;
    
    for(size_t i = 0; i < genParticles->size(); ++ i) {
      const GenParticle & p = (*genParticles)[i];
      int id = p.pdgId();
      int st = p.status();  
      int pz = p.pz(); 
      const Candidate * mom = p.mother();
      double pt = p.pt(), etagen = p.eta();
      //double phi = p.phi(), mass = p.mass();
      //double vx = p.vx(), vy = p.vy(), 
      //double vz = p.vz();
      int charge = p.charge();
      //int n = p.numberOfDaughters();
      // debug
      //      if (Rootuple->xL_gen>0 && log10(1-Rootuple->xL_gen)>-1 && Rootuple->etaMax_PF<1 && ( st ==1  ))
      //	{
      //	  cout <<" NICO XL "  <<  id << "  -st " << st << " pz " << pz << " eta " << etagen  ;   
      //	  if (mom) cout << " mother " << mom->pdgId(); 
      //	    cout << endl;
      //	}


      if ( st == 1 ) {
	   ///|| st == 8 ) { // when running over the  PYTHIA + GEANT collection
	//cout << "GEN PART " << i << " " <<  p.pdgId() << " " << p.status() << " daughters " << n << " zpos " << vz << endl;  
	if ( charge && fabs(etagen)<2.6 &&  pt >= 0.1 ) {  // !!  condition for xsec per 3 charged prompt particles
	  numseltracks++;
	  genpt.push_back(pt);
	  eta.push_back(etagen);
	}
      }
    }
    
    const int  size2 = (int) genpt.size();
    int *sorted = new int[size2];
    double *vv = new double[size2];
    for (int i=0; i<size2; i++) {
      vv[i] = genpt[i];
    }
    TMath::Sort(size2, vv, sorted, true);
    for (int i=0; i<size2; i++) {
      tracks.push_back(genpt[sorted[i]]);
      etaPT.push_back(eta[sorted[i]]);
      if (i>30) break;
    }  //  comes out size of 32!
    Rootuple->tracksPT_gen=tracks;
    Rootuple->etaOfTracksPT_gen=etaPT;   
    Rootuple->numberOfTracks_gen=tracks.size();  
    genpt.clear();
    eta.clear();
    delete [] sorted;
    delete [] vv;
  }

  // *************************************************************************
  // Zero Degree Calorimeter
  // *************************************************************************
  double ZDCMaxEnergy_minus=0;
  double ZDCMaxEnergy_plus=0;
  //ZDC_t ZDC;
  Handle <ZDCRecHitCollection> zdc_recHits_h;
  int counter=0;
  iEvent.getByLabel("zdcreco", zdc_recHits_h);
  //iEvent.getByType(zdc_recHits_h);
  const ZDCRecHitCollection *zdc_recHits = zdc_recHits_h.failedToGet()? 0 : &*zdc_recHits_h;
  
  if (zdc_recHits) { // object is available
    for (ZDCRecHitCollection::const_iterator zhit = zdc_recHits->begin(); zhit != zdc_recHits->end(); zhit++){		
      int ZDCSide      = (zhit->id()).zside();
      Float_t ZDCEnergy = zhit->energy();
      Float_t ZDCRecHitTime = zhit->time();
      int ZDCSection   = (zhit->id()).section();
      int ZDCChannel   = (zhit->id()).channel();
      //if (debug_deep) std::cout<<"ZDC Side "<<ZDCSide<<std::endl;
      //if (debug_deep) std::cout<<"ZDC Energy "<<ZDCEnergy<<std::endl;
      //if (debug_deep) std::cout<<"ZDC Section "<<ZDCSection<<std::endl;
      //if (debug_deep) std::cout<<"ZDC Channel "<<ZDCChannel<<std::endl;
      //if (debug_deep) std::cout<<"ZDC Timing "<<ZDCRecHitTime<<std::endl;

      ZDC.EventNumber=eventNumber;
      ZDC.side=ZDCSide;
      ZDC.section=ZDCSection;
      ZDC.channel=ZDCChannel;
      ZDC.energy=ZDCEnergy;
      ZDC.timing=ZDCRecHitTime;
      ZDCTree->Fill();
      //DetId CastorId= zhit->id();
      //HcalDetId CastorIdHcal(CastorId);
      if (ZDCEnergy>ZDCMaxEnergy_minus && ZDCSide==-1) ZDCMaxEnergy_minus=ZDCEnergy;
      if (ZDCEnergy>ZDCMaxEnergy_plus && ZDCSide==1) ZDCMaxEnergy_plus=ZDCEnergy;
    }
    counter++;
  }
  Rootuple->sumEZDC_minus=ZDCMaxEnergy_minus;
  Rootuple->sumEZDC_plus=ZDCMaxEnergy_plus;


  //=======================================================
  // Retrieve the Track information
  //=======================================================
  std::vector<Double_t> vertexNDOF;
  std::vector<Double_t> vertexChiNorm;
  std::vector<Double_t> vertexMolteplicity;
  double nhit=0;
  std::vector<double> V_x;
  std::vector<double> V_y;
  std::vector<double> V_z; 
  std::vector<double> pt;
  std::vector<double> tracks;
  std::vector<std::vector<double> > tracksPT;
 
  Handle<reco::VertexCollection>  vertexCollectionHandle;
  iEvent.getByLabel(PVtxCollectionTag_, vertexCollectionHandle);
  //cout <<" SIZE VTX " << vertexCollectionHandle->size() <<  endl;
  for(reco::VertexCollection::const_iterator vtx = vertexCollectionHandle->begin();vtx!=vertexCollectionHandle->end(); ++vtx)
    {
      reco::Vertex::trackRef_iterator it = vtx->tracks_begin();
      reco::Vertex::trackRef_iterator lastTrack = vtx->tracks_end();
      for(;it!=lastTrack;it++) {
	nhit+=(*it)->numberOfValidHits();
	pt.push_back((*it)->pt());
      }
      
      //cout << vtx->x() << endl;
      
      //Sorting the pt tracks, in order to take only the 31 most energetics
      const int  size = (int) pt.size();
      int *sorted = new int[size];
      double *v = new double[size];
      
      for (int i=0; i<size; i++) {
	v[i] = pt[i];
      }
      TMath::Sort(size, v, sorted, true);
      for (int i=0; i<size; i++) {
	tracks.push_back(pt[sorted[i]]);
	if (i>30) break;
      }
      tracksPT.push_back(tracks);
      tracks.clear();
      pt.clear();
      double ndof=vtx->ndof();
      double chiNorm=vtx->normalizedChi2();
      double NumbOfTracks=vtx->tracksSize();
      vertexNDOF.push_back(ndof);
      vertexChiNorm.push_back(chiNorm);
      vertexMolteplicity.push_back(NumbOfTracks);
      nhit=0;
      if ( ndof != 0 ) {
	V_x.push_back(vtx->x());
	V_y.push_back(vtx->y());
	V_z.push_back(vtx->z());
      } else {
	V_x.push_back(-999);
	V_y.push_back(-999);
	V_z.push_back(-999);
      }
      //if (ndof ==0) cout<<"VTX 0 "<<vtx->x() << " " << NumbOfTracks << " " << nhit << endl;
      //if (ndof >2) cout<<"VTX 2 "<<vtx->x() << " " << NumbOfTracks << " " << nhit << endl;
      delete [] sorted;
      delete [] v;
    } // loop over vtx
  Rootuple-> vertexMolteplicity = vertexMolteplicity;
  Rootuple-> vertexChiNorm = vertexChiNorm;
  Rootuple-> vertexNDOF = vertexNDOF;
  Rootuple-> V_z = V_z;
  Rootuple-> V_y = V_y;
  Rootuple-> V_x = V_x;
  Rootuple-> tracksPT = tracksPT;
  
  
  // *************************************************************************
  // Filling the rootuples
  // *************************************************************************
  tree_->Fill();
  return true;

}

// ------------ method called once each job just before starting event loop  -
void 
MakeRootuplaForward::beginJob() {
  //fOutputFile = new TFile(TString(outputFile_),"RECREATE"); 
  tree_= new TTree("tree_","tree_");
  tree_->Branch("Rootuple","DifNtuple",&Rootuple);
  HepHisto_PdgId = new TH1F("HepHisto_PdgId",   "Generated particles out ", 12000, -6000., 6000.);
  HistoEtaEnergyW = new TH1F("HistoEtaEnergyW",   "energy weighted dist ", 40, -5., 5.);

  //CastorTree = new TTree("CastorTree","CastorTree");
  //CastorTree->Branch("Castor",&Castor,"EventNumber/I:sector:module:timing/F:energy");
  ZDCTree = new TTree("ZDCTree","ZSCTree");
  ZDCTree->Branch("ZDC",&ZDC,"EventNumber/I:side:section:channel:energy/F:timing");

  std::vector<float> simulated;
  std::vector<float> trueD;
	
  //Calculate the distributions (our data and MC)
  for( int i=0; i<25; ++i) {
    trueD.push_back(ZSkim_2010[i]); // Name of the vector calculated with estimatedPU.py!
    simulated.push_back(probdistFlat10[i]); // Name of the vector included in Flat10.h !
  }
  
  LumiWeights_ = edm::LumiReWeighting(simulated, trueD);
}

// ------------ method called once each job just after ending the event loop  -
void 
MakeRootuplaForward::endJob() {
  //std::cout<<"End job"<<std::endl;
  //fOutputFile->cd();
  //tree_->Write();
  //CastorTree->Write();
  //HepHisto_PdgId->Write();
  //ZDCTree->Write();
  //if(fOutputFile){
    //fOutputFile->Write() ;
    //fOutputFile->Close() ;
  //}
}


//define this as a plug-in
DEFINE_FWK_MODULE(MakeRootuplaForward);
