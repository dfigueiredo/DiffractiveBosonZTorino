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

class MCGen {

 public:

  MCGen(){};
  double xL_gen;
  double xL_Num_gen;
  static const char* version(){return "$Revision: 1.1 $";}

 void Zero(){
    xL_gen=0;
    xL_Num_gen=0;
  }
  

};

