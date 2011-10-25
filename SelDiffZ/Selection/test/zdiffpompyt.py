
import FWCore.ParameterSet.Config as cms

import os

process = cms.Process("ZSelection")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_P_V20::All'

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
    "rfio:/castor/cern.ch/user/a/aproskur/z/pompyt_zg_pd_minus_reco_1.root",
    "rfio:/castor/cern.ch/user/a/aproskur/z/pompyt_zg_pd_minus_reco_2.root",
    "rfio:/castor/cern.ch/user/a/aproskur/z/pompyt_zg_pd_minus_reco_3.root",
    "rfio:/castor/cern.ch/user/a/aproskur/z/pompyt_zg_pd_minus_reco_4.root",
    "rfio:/castor/cern.ch/user/a/aproskur/z/pompyt_zg_pd_minus_reco_5.root",
    "rfio:/castor/cern.ch/user/a/aproskur/z/pompyt_zg_pd_minus_reco_6.root",
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
                                 doTheHLTAnalysis = cms.bool(True),
                                 TriggerNames = alltriggers
                                 )

process.MakeRootuplaForward = cms.EDFilter('MakeRootuplaForward',
                                           electronCollectionTag = cms.untracked.InputTag("gsfElectrons","","RECO"),
                                           filename=cms.untracked.string("SelectedRootupla.root"),
                                           CaloTowerTag= cms.InputTag("towerMaker","","RECO"),
                                           zeeCollectionTag = cms.untracked.InputTag("zeeFilter","selectedZeeCandidates","testDiffractive"),
                                           ActivateMC=cms.untracked.bool(True),
                                           PVtxCollectionTag=cms.InputTag('offlinePrimaryVertices'),
                                           TrackCollectionTag = cms.InputTag("generalTracks"),
                                           VertexCollectionTag = cms.InputTag('offlinePrimaryVertices'),
                                           fPixelClusterLabel = cms.InputTag("siPixelClusters"),
                                           #SiPixelClusteredmNewDetSetVector_siPixelClusters__RECO. 8221.09 2573.66
                                           zmumuCollectionTag = cms.untracked.InputTag("zmmCands"),
                                           electrons =cms.untracked.bool(True),
                                           muons =cms.untracked.bool(False),
                                           )


process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('histo.root')
                                   )

process.p = cms.Path(process.Selection
                     *process.MakeRootuplaForward
                     )
