// Ntuple used by ElectroWeakAnalysis
//
// Author: Matteo Marone
//
//

#include <vector>
#include <iostream>
#include "FWCore/Framework/interface/Event.h"
#include "TMath.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "TLorentzVector.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h" 
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "TTree.h"

class DiffractiveNtuple {

 public:
  
  DiffractiveNtuple(){};
  //std::vector<reco::GsfElectron> electrons;
  edm::RunNumber_t RunNumber;
  edm::EventNumber_t EventNumber;
  edm::LuminosityBlockNumber_t Lumi;
  bool HasTwoElectrons;
  bool ElectronsWithEnoughEt;
  //pat::Electron maxETelec;
  //pat::Electron SecondETelec;
  Double_t InvariantMass;
  //CaloTowerCollection tower;
  //pat::CompositeCandidate zee;
  float istlumi;
  float istlumierr;
  std::vector<double> ptTracks;
  std::vector<double> etaTracks;
  std::vector<double> ecalRawEnergyTracks;   
  std::vector<double> hcalRawEnergyTracks;   
  std::vector<int> TypeOfTracks; 
  std::vector<double> distance;   
  //std::vector<double> charge;   
  double Epz_Calo;
  double sumEHF_plus;
  double sumEHF_minus;
  double MinsumEHF_minus_plus;
  double etaMax_Calo;
  double etaMin_Calo;
  double Epz_PF_minus;
  double Epz_PF_plus;  
  double etaMax_PF;
  double etaMin_PF;
  double etaMax_PF_ChargedFromOneVertex;
  double etaMin_PF_ChargedFromOneVertex;
  double Epz_PF_minus_NOHF;
  double Epz_PF_plus_NOHF;  
  double etaMax_PF_NOHF;
  double etaMin_PF_NOHF;
  bool MoreThanOneVertex;
  double xi_gen;
  int nTowersHF_plus;
  int nTowersHF_minus;
  double etaAllTracks_PF;
  double energyTot_PF;
  double sumECastor_minus;
  double sumEZDC_minus;
  double sumEZDC_plus;
  //Montecarlo: we are "simulating" the forward calorimeters... 
  double sumECastor_plus_gen;
  double sumECastor_minus_gen;
  double sumEZDC_minus_gen;
  double sumEZDC_plus_gen;

  double etaWeightedOnEnergy_PF;
  double energyTot_PF_Barrel_minus;
  double energyTot_PF_Barrel_plus;
  double energyTot_PF_minus;
  double energyTot_PF_plus;
  double numberOfVertexesInEvent;
  double istlumiPerBX;
  double xi_PF_minus;
  double xi_PF_plus;
  double xi_PF_NOZED_minus;
  double xi_PF_NOZED_plus;
  double max_eta_gap_PF;
  double max_second_eta_gap_PF;
  double thrustValue;
  edm::Timestamp timestamp;
  double xi_PF_minus_charged;
  double xi_PF_plus_charged; 
  int nTracks_PF;
  double energyTot_PF_EE_minus;
  double energyTot_PF_EE_plus;  
  int bx;
  int numberOfVertexes;
  double sumEHF_minus_PF;
  double sumEHF_plus_PF;
  double etaOutcomingProton;
  double energyOutcomingProton;
  double ZDCRecHitTime;
  double mostEnergeticParticleEnergy_MC;
  double mostEnergeticParticleType_MC;
  double mostEnergeticParticleGap_MC;
  double mostEnergeticParticleEta_MC;
  double eMaxParticle_gen;
  double nTracks_gen;
  TTree *Castor;
  double CastoTotalEnergyTh;
  double sumECastor_Th_minus;
  int numberOfElectrons;
  Double_t etaZed;
  std::vector<int> vertexNDOF;
  double xi_PF_NOHF_minus;
  double xi_PF_NOHF_plus;
  std::vector<double> vertexChiNorm;
  std::vector<double> vertexMolteplicity;
  bool CastorActivity;
  double vertexNumberOfRecHits;
  double xL_gen;
  double xL_Type_gen;

  static const char* version(){return "$Revision: 1.1 $";}

  void Zero(){
    istlumi=0;
    bx=0;
    numberOfVertexes=0;
    //electrons.clear();
    //EventNumber=0;
    //RunNumber=0;
    Lumi=0;
    //HasTwoElectrons=false;
    //ElectronsWithEnoughEt=false;
    //InvariantMass=0;
  }
  

};

