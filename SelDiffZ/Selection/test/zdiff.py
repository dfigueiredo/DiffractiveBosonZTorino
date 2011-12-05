
import FWCore.ParameterSet.Config as cms

import os

process = cms.Process("ZSelection")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#--- for data reprocessed 2010? 
process.GlobalTag.globaltag = 'GR_R_42_V8::All'

process.load("MagneticField.Engine.uniformMagneticField_cfi") 

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True) )

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True),
                                     makeTriggerResults=cms.untracked.bool(True),
                                     )


readFilesPompyt = cms.untracked.vstring()
readFilesPompyt.extend([
    "file:/tmp/marone/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/00CF200D-907C-E011-843E-003048D4363C.root",
#    "file:///tmp/marone/store/mc/Spring11/DYtoEE_M_20_TuneD6T_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/00CF200D-907C-E011-843E-003048D4363C.root",
#    "rfio:/castor/cern.ch/user/a/aproskur/z/pompyt_zg_pd_minus_reco_1.root",
#    "rfio:/castor/cern.ch/user/a/aproskur/z/pompyt_zg_pd_minus_reco_2.root",
#    "rfio:/castor/cern.ch/user/a/aproskur/z/pompyt_zg_pd_minus_reco_3.root",
#    "rfio:/castor/cern.ch/user/a/aproskur/z/pompyt_zg_pd_minus_reco_4.root",
#    "rfio:/castor/cern.ch/user/a/aproskur/z/pompyt_zg_pd_minus_reco_5.root",
#    "rfio:/castor/cern.ch/user/a/aproskur/z/pompyt_zg_pd_minus_reco_6.root",
    ]
                       )


process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
     reportEvery = cms.untracked.int32(500),
 )

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames=readFilesPompyt,
                            #fileNames =cms.untracked.vstring('file:/tmp/20467E47-D07B-E011-AF09-002618943861.root'),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            )

trigger2011v1  = cms.vstring("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3","HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v3","HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3","HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3","HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2","HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v3","HLT_Ele45_CaloIdVT_TrkIdT_v3","HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v4","HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v4","HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v4")

trigger2010    = cms.vstring("HLT_Ele17_CaloIdl_Ele8_CaloIsoIdL_CaloIsoVL_v3","HLT_Ele15_SW_L1R","HLT_Ele15_SW_CaloEleId_L1R","HLT_Ele17_SW_CaloEleId_L1R","HLT_Ele17_SW_TightEleId_L1R","HLT_Ele17_SW_TightEleId_L1R_v2","HLT_Ele17_SW_TightEleId_L1R_v3","HLT_Photon10_L1R","HLT_Photon15_L1R","HTL_Photon15_Cleaned_L1R")

alltriggers    = cms.vstring() # In this way, the HLT string is empty and it will trigger every event

process.Selection = cms.EDFilter('ZanalyzerFilter',
                                 electronCollection = cms.InputTag("gsfElectrons"),
                                 triggerCollectionTag = cms.untracked.InputTag("TriggerResults","","HLT"),
                                 filename=cms.untracked.string("ZAnalysisFilter.root"),
                                 UseCombinedPrescales = cms.bool(True),
                                 doTheHLTAnalysis = cms.bool(False),
                                 TriggerNames = alltriggers
                                 )

process.MakeRootuplaForward = cms.EDFilter('MakeRootuplaForward',
                                           electronCollectionTag = cms.untracked.InputTag("gsfElectrons","","RECO"),
                                           filename=cms.untracked.string("SelectedRootupla.root"),
                                           CaloTowerTag= cms.InputTag("towerMaker","","RECO"),
                                           zeeCollectionTag = cms.untracked.InputTag("zeeFilter","selectedZeeCandidates","testDiffractive"),
                                           ActivateMC=cms.untracked.bool(False),
                                           PVtxCollectionTag=cms.InputTag('offlinePrimaryVertices'),
                                           TrackCollectionTag = cms.InputTag("generalTracks"),
                                           VertexCollectionTag = cms.InputTag('offlinePrimaryVertices'),
                                           fPixelClusterLabel = cms.InputTag("siPixelClusters"),
                                           zmumuCollectionTag = cms.untracked.InputTag("zmmCands"),
                                           electrons =cms.untracked.bool(True),
                                           muons =cms.untracked.bool(False),
                                           )

process.load("JetCollections_cfi")


)

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('histo.root')
                                   )

process.p = cms.Path(
        process.PFJetPath*
        process.Selection*
        #process.demo*
        process.MakeRootuplaForward
                     )
