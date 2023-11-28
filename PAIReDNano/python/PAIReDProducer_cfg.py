import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
       '/store/mc/Run3Summer22MiniAODv4/ZH_Hto2C_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/2530000/d848840d-0419-4863-93ff-d096472583e2.root'
    )
)

process.content = cms.EDAnalyzer("EventContentAnalyzer")

process.PAIReD = cms.EDProducer('PAIReDONNXProducer')

process.out = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAODSIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:testPAIReDTagger.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        "keep nanoaodFlatTable_*Table_*_*",     # event data
        "keep edmTriggerResults_*_*_*",  # event data
        "keep String_*_genModel_*",  # generator model data
        "keep nanoaodMergeableCounterTable_*Table_*_*", # accumulated per/run or per/lumi data
        "keep nanoaodUniqueString_nanoMetadata_*_*",   # basic metadata
    )
)

process.p = cms.Path(process.content*process.PAIReD)
process.e = cms.EndPath(process.out)