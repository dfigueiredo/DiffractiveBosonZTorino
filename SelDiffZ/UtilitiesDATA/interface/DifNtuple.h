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
#include <TH1F.h>
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

class DifNtuple {

 public:
  
  DifNtuple(){};
  //std::vector<reco::GsfElectron> electrons;
  edm::RunNumber_t RunNumber;
  edm::EventNumber_t EventNumber;
  edm::LuminosityBlockNumber_t Lumi;
  //pat::Electron maxETelec;
  //pat::Electron SecondETelec;
  Double_t InvariantMass;
  //CaloTowerCollection tower;
  //pat::CompositeCandidate zee;
  float istlumi;
  float istlumierr;
  //std::vector<double> distance;   
  //std::vector<double> charge;   
  double Epz_Calo_plus;
  double Epz_Calo_minus;
  double sumEHF_plus;
  double sumEHF_minus;
  double minEHF;
  double etaMax_Calo;
  double etaMin_Calo;
  double Epz_PF_minus;
  double Epz_PF_plus;  
  double etaMax_PF;
  double etaMin_PF;
  double etaMax_Charged_PV_PF;
  double etaMin_Charged_PV_PF;
  double Epz_NOHF_PF_minus;
  double Epz_NOHF_PF_plus;  
  double etaMax_NOHF_PF;
  double etaMin_NOHF_PF;
  double xi_gen;
  double Mx2_gen;
  double xi_NOZ_gen;
  int nTowersHF_plus;
  int nTowersHF_minus;
  double etaAllTracks_PF;
  double energyTot_PF;
  double sumEZDC_minus;
  double sumEZDC_plus;
  //Montecarlo: we are "simulating" the forward calorimeters... 
  double sumECastor_gen_plus;
  double sumECastor_gen_minus;
  double sumEZDC_gen_minus;
  double sumEZDC_gen_plus;

  double etaWeightedOnEnergy_PF;
  double energyTot_PF_Barrel_minus;
  double energyTot_PF_Barrel_plus;
  double energyTot_PF_minus;
  double energyTot_PF_plus;
  float istlumiPerBX;
  double xi_PF_minus;
  double xi_PF_plus;
  double xi_Z_minus;
  double xi_Z_plus;
  double xi_Z_gen_minus;
  double xi_Z_gen_plus;
  double max_eta_gap_PF;
  double max_second_eta_gap_PF;
  //double thrustValue;
  //double thrustX;
  //double thrustY;
  //double thrustZ;
  //double sphericity;
  //double planarity;
  //double aplanarity;
  edm::Timestamp timestamp;
  double xi_PF_charged_minus;
  double xi_PF_charged_plus; 
  int nTracks_PF;
  double energyTot_PF_EE_minus;
  double energyTot_PF_EE_plus;  
  int bx;
  int numberOfVertexes;
  double sumEHF_PF_minus;
  double sumEHF_PF_plus;
  double etaOutcomingProton;
  double mostEnergeticParticleGap_MC;
  double nParticles_gen;
  double sumECastor_Th_minus;
  int numberOfLeptons;
  Double_t etaZ;
  std::vector<Double_t> vertexNDOF;
  double xi_PF_NOHF_minus;
  double xi_PF_NOHF_plus;
  std::vector<double> vertexChiNorm;
  std::vector<double> vertexMolteplicity;
  bool CastorActivity;
  std::vector<double> vertexNumberOfRH;
  double xL_gen;
  double xL_Num_gen;
  double sumECastorRaw_minus;
  double sumECastor_minus;
  double pixelNCluster;
  std::vector<double> V_x;
  std::vector<double> V_y;
  std::vector<double> V_z;
  double PV_x;
  double PV_y;
  double PV_z;
  double etaZ_gen;
  double lepton1Phi;
  double lepton1Eta;
  double lepton2Phi;
  double lepton2Eta;
  double Mx2;
  double M_x;
  double M_y;
  double M_z;
  double energyZ_gen;
  double xi_Calo_minus;
  double xi_Calo_plus;
  double p_diss_mass_gen;
  double  xL_p_diss;
  std::vector<std::vector<double> > tracksPT;
  std::vector<double> tracksPT_gen;
  int numberoOfTracks_gen; 
  std::vector<double> etaOfTracksPT_gen;
  //bool hltMatch;
  //int hlt9Match;
  //int hlt20Match;
  std::vector<double> electronEnergy;
  std::vector<double> muEnergy;

  int   PU_NumInt;
  std::vector<float> PU_zpos;
  std::vector<int>   PU_ntrks_lowpT;
  std::vector<int>   PU_ntrks_highpT;
  std::vector<float> PU_sumpT_lowpT;
  std::vector<float> PU_sumpT_highpT;



  static const char* version(){return "$Revision: 1.1 $";}

  void Zero(){
    istlumi=0;
    bx=0;
    numberOfVertexes=0;
    Lumi=0;
    InvariantMass=0;
    istlumi=-1;
    istlumierr=-1;
    Epz_Calo_plus=-1;
    Epz_Calo_minus=-1;
    sumEHF_plus=-1;
    sumEHF_minus=-1;
    minEHF=-1;
    etaMax_Calo=-999;
    etaMin_Calo=-999;
    Epz_PF_minus=-1;
    Epz_PF_plus=-1;  
    etaMax_PF=-999;
    etaMin_PF=-999;
    etaMax_Charged_PV_PF=-999;
    etaMin_Charged_PV_PF=-999;
    Epz_NOHF_PF_minus=-1;
    Epz_NOHF_PF_plus=-1;  
    etaMax_NOHF_PF=-999;
    etaMin_NOHF_PF=-999;
    xi_gen=-999;
    Mx2_gen=-999;
    xi_NOZ_gen=-999;
    nTowersHF_plus=-1;
    nTowersHF_minus=-1;
    etaAllTracks_PF=-999;
    energyTot_PF=-1;
    sumEZDC_minus=-1;
    sumEZDC_plus=-1;
    sumECastor_gen_plus=-1;
    sumECastor_gen_minus=-1;
    sumEZDC_gen_minus=-1;
    sumEZDC_gen_plus=-1;
    etaWeightedOnEnergy_PF=-999;
    energyTot_PF_Barrel_minus=-1;
    energyTot_PF_Barrel_plus=-1;
    energyTot_PF_minus=-1;
    energyTot_PF_plus=-1;
    istlumiPerBX=-999;
    xi_PF_minus=-1;
    xi_PF_plus=-1;
    xi_Z_minus=-1;
    xi_Z_plus=-1;
    xi_Z_gen_minus=-1;
    xi_Z_gen_plus=-1;
    max_eta_gap_PF=-999;
    max_second_eta_gap_PF=-999;
    //thrustValue=-999;
    //thrustX=-999;
    //thrustY=-999;
    //thrustZ=-999;
    //sphericity=-999;
    //planarity=-999;
    //aplanarity=-999;
    xi_PF_charged_minus=-1;
    xi_PF_charged_plus=-1; 
    nTracks_PF=-1;
    energyTot_PF_EE_minus=-1;
    energyTot_PF_EE_plus=-1;  
    bx=-999;
    numberOfVertexes=-1;
    sumEHF_PF_minus=-1;
    sumEHF_PF_plus=-1;
    etaOutcomingProton=-999;
    mostEnergeticParticleGap_MC=-1;
    nParticles_gen=-1;
    sumECastor_Th_minus=-1;
    numberOfLeptons=-1;
    etaZ=-999;
    vertexNDOF.clear();
    xi_PF_NOHF_minus=-1;
    xi_PF_NOHF_plus=-1;
    vertexChiNorm.clear();
    vertexMolteplicity.clear();
    CastorActivity=false;
    vertexNumberOfRH.clear();
    xL_gen=-999;
    xL_Num_gen=-999;
    sumECastorRaw_minus=-1;
    sumECastor_minus=-1;
    pixelNCluster=-1;
    V_x.clear();
    V_y.clear();
    V_z.clear();
    PV_x=-999;
    PV_y=-999;
    PV_z=-999;
    etaZ_gen=-999;
    lepton1Phi=-999;
    lepton1Eta=-999;
    lepton2Phi=-999;
    lepton2Eta=-999;
    Mx2=-1;
    M_x=-1;
    M_y=-1;
    M_z=-1;
    energyZ_gen=-1;
    xi_Calo_minus=-1;
    xi_Calo_plus=-1;
    p_diss_mass_gen=-1;
    xL_p_diss=-999;

    PU_NumInt = 0;
    PU_zpos.clear();
    PU_ntrks_lowpT.clear();
    PU_ntrks_highpT.clear();
    PU_sumpT_lowpT.clear();
    PU_sumpT_highpT.clear();
    
  }
  

};

