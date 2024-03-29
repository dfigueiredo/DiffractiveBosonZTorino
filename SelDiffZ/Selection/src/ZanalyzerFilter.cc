// -*- C++ -*-
//
// Package:    Zanalyzer
// Class:      Zanalyzer
// 
/**\class Zanalyzer Zanalyzer.cc Zmonitoring/Zanalyzer/src/Zanalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Vieri Candelise & Matteo Marone
//         Created:  Wed May 11 14:53:26 CEST 2011
// $Id: ZanalyzerFilter.cc,v 1.3 2011/12/05 14:34:24 arcidiac Exp $
//
//


// system include files
#include <memory>

// user include files
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SelDiffZ/Selection/interface/ZanalyzerFilter.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include <string>
#include <cmath>
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <iostream>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"


using namespace std;
using namespace edm;
using namespace reco;





//
// member functions
//

// ------------ method called for each event  ------------
bool
ZanalyzerFilter::filter (edm::Event & iEvent, edm::EventSetup const & iSetup)
{
  if (debug) cout<<"------- NEW Event -----"<<endl;
  using namespace edm;
  Handle < GsfElectronCollection > electronCollection;
  iEvent.getByLabel (theElectronCollectionLabel, electronCollection);


  if(ActivateMC_) 
    {
      /// For every event, fill gnerated histograms

      Handle<reco::GenParticleCollection> genParticles;     
      iEvent.getByLabel("genParticles",genParticles);  // standard PYTHIA collection
      
      for(size_t i = 0; i < genParticles->size(); ++ i) {
	const GenParticle & p = (*genParticles)[i];
	int pdg = p.pdgId();
	int status = p.status();  	 
	//double eta_gen = p.eta();
	//double part_pt = p.pt();
	double ener_gen= p.energy();
	double px_gen  = p.px();
	double py_gen  = p.py();
	double pz_gen  = p.pz();
	double mass_gen= p.mass();
	//bool electronFromZ=false;
	int motherId=0;
	
	if (  p.mother() ) motherId =  p.mother()->pdgId();
	

	if (pdg==23 && status == 3 ) {
	  
	  h_invMassGEN->Fill(mass_gen);
	}
	if ( fabs(pdg) == 11 && status == 0 ) {
	  
	  if ( motherId == 23 ) {
	    //cout << " ele found " << pdg << endl;
	    h_elePxGEN->Fill(px_gen);
	    h_elePyGEN->Fill(py_gen);
	    h_elePzGEN->Fill(pz_gen);
	    h_eleEnergyGEN->Fill(ener_gen);	   
	  }
	}	
      }
    }   //  if Activate_MC

  if (!electronCollection.isValid ())
    return false;

  //////////////////////
  //Match The HLT Trigger
  //////////////////////

  edm::Handle<edm::TriggerResults> HLTResults;
  iEvent.getByLabel(triggerCollection_, HLTResults);
  if (!HLTResults.isValid ()) return false;
  
  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*HLTResults);   
  bool flag=false;
   
  if (HLTResults.isValid() && doTheHLTAnalysis_) {
    /// Storing the Prescale information: loop over the triggers and record prescale
    unsigned int minimalPrescale(10000);
    unsigned int prescale(0);
    bool bit(true);
    std::pair<int,int> prescalepair;
    std::vector<int>  triggerSubset;
    for(unsigned int itrig = 0; itrig < triggerNames_.size(); ++itrig) {
      if(triggerIndices_[itrig]!=2048) {
	// check trigger response
	bit = HLTResults->accept(triggerIndices_[itrig]);
	triggerSubset.push_back(bit);
	if(bit) {
	  flag=true;
	  if (debug) cout<<"Matched "<<triggerNames.triggerName(itrig)<<endl;
	  int prescaleset = hltConfig_.prescaleSet(iEvent,iSetup);
	  if(prescaleset!=-1) {
	    prescalepair = hltConfig_.prescaleValues(iEvent,iSetup,triggerNames_[itrig]);
	    if (debug) cout<<"prescale.first "<<prescalepair.first<<" prescalepair.second "<<prescalepair.second<<endl;
	    if((useCombinedPrescales_ && prescalepair.first<0) || prescalepair.second<0) {
	      edm::LogWarning("ZAnalyzerFilter") << " Unable to get prescale from event for trigger " << triggerNames.triggerName(itrig) << " :" 
						 << prescalepair.first << ", " << prescalepair.second;
	    }
	    prescale = useCombinedPrescales_ ? prescalepair.first*prescalepair.second : prescalepair.second;
	    minimalPrescale = minimalPrescale <  prescale ? minimalPrescale : prescale;
	    if (debug) cout<<"prescale "<<prescale<<" minimal Prescale "<<minimalPrescale<<" for trigger "<<triggerNames.triggerName(itrig)<<endl;
	  } 
	}
      }
      else {
	// that trigger is presently not in the menu
	triggerSubset.push_back(false);
      }
    }
  }
  
  if (!flag) 
    {
      if(!useAllTriggers_) return false;
    }
  
  bool isBarrelElectrons;
  bool isEndcapElectrons;
  bool isIsolatedBarrel;
  bool isIDBarrel;
  bool isConvertedBarrel;
  bool isIsolatedEndcap;
  bool isIDEndcap;
  bool isConvertedEndcap;
  int elIsAccepted=0;
  int elIsAcceptedEB=0;
  int elIsAcceptedEE=0;
  
  std::vector<TLorentzVector> LV;
  
  

  if (electronCollection->size()==1) return false;
  
  for (reco::GsfElectronCollection::const_iterator recoElectron = electronCollection->begin (); recoElectron != electronCollection->end (); recoElectron++) {

    if (recoElectron->et () <= 25)  continue;
    
    // Define Isolation variables
    double IsoTrk = (recoElectron->dr03TkSumPt () / recoElectron->et ());
    double IsoEcal = (recoElectron->dr03EcalRecHitSumEt () / recoElectron->et ());
    double IsoHcal = (recoElectron->dr03HcalTowerSumEt () / recoElectron->et ());
    double HE = recoElectron->hadronicOverEm();
    
    //Define ID variables
    
    float DeltaPhiTkClu = recoElectron->deltaPhiSuperClusterTrackAtVtx ();
    float DeltaEtaTkClu = recoElectron->deltaEtaSuperClusterTrackAtVtx ();
    float sigmaIeIe = recoElectron->sigmaIetaIeta ();
    
    //Define Conversion Rejection Variables

    float Dcot = recoElectron->convDcot ();
    float Dist = recoElectron->convDist ();
    int NumberOfExpectedInnerHits = recoElectron->gsfTrack ()->trackerExpectedHitsInner ().numberOfHits ();

    //quality flags

    isBarrelElectrons = false;
    isEndcapElectrons = false;
    isIsolatedBarrel = false;
    isIDBarrel = false;
    isConvertedBarrel = false;
    isIsolatedEndcap = false;
    isIDEndcap = false;
    isConvertedEndcap = false;
 
    /***** Barrel WP80 Cuts *****/

    if (fabs (recoElectron->eta ()) <= 1.4442) {

      /* Isolation */
      if (IsoTrk < 0.09 && IsoEcal < 0.07 && IsoHcal < 0.10) {
	isIsolatedBarrel = true;
      }

      /* Identification */
      if (fabs (DeltaEtaTkClu) < 0.004 && fabs (DeltaPhiTkClu) < 0.06
	  && sigmaIeIe < 0.01 && HE < 0.04) {
	isIDBarrel = true;
      }

      /* Conversion Rejection */
      if ((fabs (Dist) >= 0.02 || fabs (Dcot) >= 0.02)
	  && NumberOfExpectedInnerHits == 0) {
	isConvertedBarrel = true;
      }
      //cout<<"isIsolatedBarrel "<<isIsolatedBarrel<<" isIDBarrel "<< isIDBarrel<<" isConvertedBarrel "<<isConvertedBarrel<<endl;
    }

    if (isIsolatedBarrel && isIDBarrel && isConvertedBarrel) {
      elIsAccepted++;
      elIsAcceptedEB++;
      TLorentzVector b_e2;
      b_e2.SetPtEtaPhiM(recoElectron->pt(),recoElectron->eta(),recoElectron->phi(), 0.0);
      //      TLorentzVector b_e2(recoElectron->momentum ().x (),recoElectron->momentum ().y (),recoElectron->momentum ().z (), recoElectron->p ());
      LV.push_back(b_e2);
    }

    /***** Endcap WP80 Cuts *****/

    if (fabs (recoElectron->eta ()) >= 1.5660
	&& fabs (recoElectron->eta ()) <= 2.5000) {

      /* Isolation */
      if (IsoTrk < 0.04 && IsoEcal < 0.05 && IsoHcal < 0.025) {
	isIsolatedEndcap = true;
      }

      /* Identification */
      if (fabs (DeltaEtaTkClu) < 0.007 && fabs (DeltaPhiTkClu) < 0.03
	  && sigmaIeIe < 0.03 && HE < 0.15) {
	isIDEndcap = true;
      }

      /* Conversion Rejection */
      if ((fabs (Dcot) > 0.02 || fabs (Dist) > 0.02)
	  && NumberOfExpectedInnerHits == 0) {
	isConvertedEndcap = true;
      }
      //cout<<"isIsolatedEndcap "<<isIsolatedEndcap<<" isIDEndcap "<< isIDEndcap<<" isConvertedEndcap "<<isConvertedEndcap<<endl;
    }

    if (isIsolatedEndcap && isIDEndcap && isConvertedEndcap) {
      elIsAccepted++;
      elIsAcceptedEE++;
      TLorentzVector e_e2;
      e_e2.SetPtEtaPhiM(recoElectron->pt(),recoElectron->eta(),recoElectron->phi(), 0.0);
      LV.push_back(e_e2);
    }

  }
  if (elIsAccepted<=1)    return false;
   double e_ee_invMass=0; 
   if (elIsAccepted>2) cout<<"WARNING: In this events we have more than two electrons accpeted!!!!!!!"<<endl;
   if (LV.size()==2){
     TLorentzVector e_pair = LV[0] + LV[1];
     e_ee_invMass = e_pair.M ();
     h_invMass->Fill(e_ee_invMass);
   }  
   
   if (elIsAcceptedEB==2){
     h_invMassBB->Fill(e_ee_invMass);
   }
   if (elIsAcceptedEE==2){
     h_invMassEE->Fill(e_ee_invMass);
   }
   if (elIsAcceptedEB==1 && elIsAcceptedEE==1){
     h_invMassEB->Fill(e_ee_invMass);
   }
   
  eventAccept->Fill(elIsAccepted);
  LV.clear();
  
  return true;
 
}


// ------------ method called once each job just before starting event loop  ------------
void
ZanalyzerFilter::beginJob (){

//beginJob

}


// ------------ method called when starting to processes a run  ------------
bool
ZanalyzerFilter::beginRun(edm::Run &iRun, edm::EventSetup const& iSetup)
{
//HLT names
   std::vector<std::string>  hlNames;
   bool changed (true);
   if (hltConfig_.init(iRun,iSetup,triggerCollection_.process(),changed)) {
     if (changed) {
       hlNames = hltConfig_.triggerNames();
     }
   } else {
     edm::LogError("MyAnalyzer") << " HLT config extraction failure with process name " << triggerCollection_.process();
   }
   if (debug) cout<<"useAllTriggers?"<<useAllTriggers_<<endl;
   if(useAllTriggers_) triggerNames_ = hlNames;
   //triggerNames_ = hlNames;
   //HLT indices
   triggerIndices_.clear();
   for(unsigned int itrig = 0; itrig < triggerNames_.size(); ++itrig) {
     if(find(hlNames.begin(),hlNames.end(),triggerNames_[itrig])!=hlNames.end())
       triggerIndices_.push_back(hltConfig_.triggerIndex(triggerNames_[itrig]));
     else
       triggerIndices_.push_back(2048);
   }

   if (debug){
     // text (debug) output
     int i=0;
     for(std::vector<std::string>::const_iterator it = triggerNames_.begin(); it<triggerNames_.end();++it) {
       std::cout << (i++) << " = " << (*it) << std::endl;
     } 
   }
   return true;
}

// ------------ method called once each job just after ending the event loop  ------------
void
ZanalyzerFilter::endJob ()
{

//endJob

}

DEFINE_FWK_MODULE (ZanalyzerFilter);
