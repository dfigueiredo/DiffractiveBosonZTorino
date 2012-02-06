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
  edm::LuminosityBlockNumber_t LumiSection;
  edm::Timestamp timestamp;
  int bx;
  float istlumi;
  float istlumierr;
  float istlumiPerBX;

  double etaZ_gen;
  Double_t ZMass;
  Double_t etaZ;
  double energyZ_gen;
  double lepton1Phi;
  double lepton1Eta;
  double lepton2Phi;
  double lepton2Eta;
  std::vector<double> electronEnergy;
  std::vector<double> muEnergy;

  //CaloTowerCollection tower;
  //pat::CompositeCandidate zee;
  //std::vector<double> distance;   
  //std::vector<double> charge;   
  double Epz_Calo_minus;
  double Epz_Calo_plus;
  double Epz_PF_minus;
  double Epz_PF_plus;  
  double Epz_NOHF_PF_minus;
  double Epz_NOHF_PF_plus;  

  double etaMax_Calo;
  double etaMin_Calo;
  double etaMax_PF;
  double etaMin_PF;
  double etaMax_Charged_PV_PF;
  double etaMin_Charged_PV_PF;
  double etaMax_NOHF_PF;
  double etaMin_NOHF_PF;
  double etaWeightedOnEnergy_PF;
  double etaAllTracks_PF;
  double max_eta_gap_PF;
  double max_second_eta_gap_PF;
  double eta_gap_limplus;
  double max_eta_gap_gen;
  double max_second_eta_gap_gen;
  double eta_gap_limplus_gen;

  double Mx2;
  double Mx2_plus;
  double Mx2_minus;
  double P_x;
  double P_y;
  double P_z;
  double Mx2_gen;
  double Mx2_NOZ_gen;
  double Mx2_plus_gen;
  double Mx2_minus_gen;

  int nTowersHF_plus;
  int nTowersHF_minus;
  double sumEHF_plus;
  double sumEHF_minus;
  double sumEHF_L_plus;
  double sumEHF_L_minus;
  double sumEHF_S_plus;
  double sumEHF_S_minus;

  double sumEHF_PF_minus;
  double sumEHF_PF_plus;
  double minEHF;
  double sumEZDC_minus;
  double sumEZDC_plus;
  double sumECastor_gen_plus;
  double sumECastor_gen_minus;
  double sumEZDC_gen_minus;
  double sumEZDC_gen_plus;
  double sumECASTOR_minus;

  double energyTot_PF;
  double energyTot_PF_Barrel_minus;
  double energyTot_PF_Barrel_plus;
  double energyTot_PF_minus;
  double energyTot_PF_plus;
  double energyTot_PF_EE_minus;
  double energyTot_PF_EE_plus;  

  double xi_PF_minus;
  double xi_PF_plus;
  double xi_Z_minus;
  double xi_Z_plus;
  double xi_Z_gen_minus;
  double xi_Z_gen_plus;
  double xi_PV_PF_charged_minus;
  double xi_PV_PF_charged_plus; 
  double xi_PF_NOHF_minus;
  double xi_PF_NOHF_plus;
  double xi_Calo_minus;
  double xi_Calo_plus;
  double xL_gen;
  double xL_Num_gen;
  double p_diss_mass_gen;
  double xL_p_diss;

  int nPart_PF;
  int N_mx2plus;
  int N_mx2minus;
  int N_mx2plus_gen;
  int N_mx2minus_gen;
  double etaOutcomingProton;
  double nParticles_gen;
  int numberOfLeptons;

  int numberOfVertexes;
  std::vector<std::vector<double> > tracksPT;
  std::vector<Double_t> vertexNDOF;
  std::vector<double> vertexChiNorm;
  std::vector<double> vertexMolteplicity;
  std::vector<double> V_x;
  std::vector<double> V_y;
  std::vector<double> V_z;

  int numberOfTracks_gen; 
  std::vector<double> tracksPT_gen;
  std::vector<double> etaOfTracksPT_gen;

  double PUMCweight;
  int   PU_NumInt;
  std::vector<float> PU_zpos;
  std::vector<int>   PU_ntrks_lowpT;
  std::vector<int>   PU_ntrks_highpT;
  std::vector<float> PU_sumpT_lowpT;
  std::vector<float> PU_sumpT_highpT;

  std::vector<float> EnergyInEta;
  std::vector<float> EnergyInEtaHFL;
  std::vector<float> EnergyInEtaHFS;
  std::vector<float> EnergyCastorModule;

  static const char* version(){return "$Revision: 1.4 $";}

  void Zero(){
    istlumi=0;
    bx=0;
    numberOfVertexes=0;
    LumiSection=0;
    ZMass=0;
    istlumi=-1;
    istlumierr=-1;
    Epz_Calo_plus=-1;
    Epz_Calo_minus=-1;
    sumEHF_plus=-1;
    sumEHF_minus=-1;
    sumEHF_L_plus=-1;
    sumEHF_L_minus=-1;
    sumEHF_S_plus=-1;
    sumEHF_S_minus=-1;
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
    Mx2_gen=-999;
    Mx2_NOZ_gen=-999;
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
    eta_gap_limplus=-10;
    max_eta_gap_gen=-999;
    max_second_eta_gap_gen=-999;
    eta_gap_limplus_gen=-10;
    xi_PV_PF_charged_minus=-1;
    xi_PV_PF_charged_plus=-1; 
    nPart_PF=-1;
    N_mx2plus=-1;
    N_mx2minus=-1;
    N_mx2plus_gen=-1;
    N_mx2minus_gen=-1;
    energyTot_PF_EE_minus=-1;
    energyTot_PF_EE_plus=-1;  
    bx=-999;
    numberOfVertexes=-1;
    sumEHF_PF_minus=-1;
    sumEHF_PF_plus=-1;
    etaOutcomingProton=-999;
    nParticles_gen=-1;
    numberOfLeptons=-1;
    etaZ=-999;
    vertexNDOF.clear();
    xi_PF_NOHF_minus=-1;
    xi_PF_NOHF_plus=-1;
    vertexChiNorm.clear();
    vertexMolteplicity.clear();
    xL_gen=-999;
    xL_Num_gen=-999;
    sumECASTOR_minus=-1;
    V_x.clear();
    V_y.clear();
    V_z.clear();
    etaZ_gen=-999;
    lepton1Phi=-999;
    lepton1Eta=-999;
    lepton2Phi=-999;
    lepton2Eta=-999;
    Mx2=-1;
    Mx2_plus=-1;
    Mx2_minus=-1;
    Mx2_plus_gen=-1;
    Mx2_minus_gen=-1;
    P_x=-1;
    P_y=-1;
    P_z=-1;
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
    EnergyInEta.clear();
    EnergyInEtaHFL.clear();
    EnergyInEtaHFS.clear();
    EnergyCastorModule.clear();
 }
  

};

