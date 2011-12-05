#ifndef MakeRootuplaForward_H
#define MakeRootuplaForward_H

/******************************************************************************
 *
 * Implementation Notes:
 *
 * contact:
 * Matteo.Marone@cern.ch
 *
 * 09 November 2010
 *
 *****************************************************************************/

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
#include <vector>
#include <iostream>
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/FWLite/interface/LuminosityBlock.h"
//
#include "TString.h"
#include "TMath.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "SelDiffZ/UtilitiesDATA/interface/DifNtuple.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "TFile.h"
#include "TTree.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

class DifNtuple;
DifNtuple *Rootuple;

struct Data {
  edm::EventNumber_t EventNumber;
  int           sector;
  int           module;
  float         timing;
  float         energy;
};

struct zdc {
  edm::EventNumber_t EventNumber;
  int side;
  int section;
  int channel;
  Float_t energy;
  Float_t timing;
};

class MakeRootuplaForward : public edm::EDFilter {
public:
  explicit MakeRootuplaForward(const edm::ParameterSet&);
  ~MakeRootuplaForward();
  TFile *fOutputFile;
  TTree *tree_;
  TTree *CastorTree;
  Data Castor; 
  zdc ZDC;
  TTree *ZDCTree;
  TH1F *HepHisto_PdgId;
  //TH1F HepHisto_PdgId("HepHisto_PdgId",   "Generated particles out ", 12000, -6000., 6000.);
  TH1F *HistoEtaEnergyW;

private:
  virtual void beginJob();
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  std::string   hltpath_;
  edm::InputTag triggerCollectionTag_; 
  edm::InputTag zeeCollectionTag_;
  edm::InputTag zmumuCollectionTag_;
  edm::InputTag electronCollectionTag_;
  edm::InputTag caloTowerTag_;
  edm::InputTag metCollectionTag_;
  edm::InputTag PVtxCollectionTag_;
  edm::InputTag TrackCollectionTag_;
  edm::InputTag VertexCollectionTag_;
  std::string outputFile_;
  edm::InputTag fPixelClusterLabel;
  bool ActivateMC_;
  bool useMB_;
  bool JPsi_;
  bool electrons_;
  bool muons_;
  bool NoMassCuts_;
  bool triggerAnalysis_;
  edm::LumiReWeighting LumiWeights_;

};
#endif
