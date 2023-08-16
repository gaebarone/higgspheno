import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

def addPFCands(process, runOnMC=False):
    process.customizedPFCandsTask = cms.Task( )
    process.schedule.associate(process.customizedPFCandsTask)
    candInput = cms.InputTag("packedPFCandidates")
    '''process.customConstituentsExtTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                                                        src = candInput,
                                                        cut = cms.string(""), #we should not filter after pruning
                                                        name = cms.string("PAIReDCands"),
                                                        doc = cms.string("interesting particles from AK4 and AK8 jets"),
                                                        singleton = cms.bool(False), # the number of entries is variable
                                                        extension = cms.bool(False), # this is the extension table for the AK8 constituents
                                                        variables = cms.PSet(CandVars,
                                                            puppiWeight = Var("puppiWeight()", float, doc="Puppi weight",precision=10),
                                                            puppiWeightNoLep = Var("puppiWeightNoLep()", float, doc="Puppi weight removing leptons",precision=10),
                                                            vtxChi2 = Var("?hasTrackDetails()?vertexChi2():-1", float, doc="vertex chi2",precision=10),
                                                            trkChi2 = Var("?hasTrackDetails()?pseudoTrack().normalizedChi2():-1", float, doc="normalized trk chi2", precision=10),
                                                            dz = Var("?hasTrackDetails()?dz():-1", float, doc="pf dz", precision=10),
                                                            dzErr = Var("?hasTrackDetails()?dzError():-1", float, doc="pf dz err", precision=10),
                                                            d0 = Var("?hasTrackDetails()?dxy():-1", float, doc="pf d0", precision=10),
                                                            d0Err = Var("?hasTrackDetails()?dxyError():-1", float, doc="pf d0 err", precision=10),
                                                            pvAssocQuality = Var("pvAssociationQuality()", int, doc="primary vertex association quality"),
                                                            lostInnerHits = Var("lostInnerHits()", int, doc="lost inner hits"),
                                                            trkQuality = Var("?hasTrackDetails()?pseudoTrack().qualityMask():0", int, doc="track quality mask"),
                                                         )
                                    )'''
    process.PAIReDJetTable = cms.EDProducer("PAIReDJetProducer",
                                                        candidates = candInput,
                                                        jets = cms.InputTag("finalJets"),
                                                        name = cms.string("PAIReDJets"),
                                                        )
    process.customizedPFCandsTask.add(process.customConstituentsExtTable)
    process.customizedPFCandsTask.add(process.PAIReDJetTable)

    if runOnMC:
        process.genJetsConstituents = cms.EDProducer("GenJetPackedConstituentPtrSelector",
                                                    src = cms.InputTag("slimmedGenJets"),
                                                    cut = cms.string("pt > 20")
                                                    )
        genCandInput = cms.InputTag("packedGenParticles")
        '''process.genJetsParticleTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                                                         src = genCandInput,
                                                         cut = cms.string(""), #we should not filter after pruning
                                                         name= cms.string("PAIReDGenCands"),
                                                         doc = cms.string("interesting gen particles from AK4 and AK8 jets"),
                                                         singleton = cms.bool(False), # the number of entries is variable
                                                         extension = cms.bool(False), # this is the main table for the AK8 constituents
                                                         variables = cms.PSet(CandVars)
                                                         )'''
        process.genJetsPAIReD = cms.EDProducer("PAIReDGenJetProducer",
                                                         candidates = genCandInput,
                                                         jets = cms.InputTag("genJetsAK4Constituents"), # Note: The name has "Constituents" in it, but these are the jets
                                                         name = cms.string("PAIReDGenJets"),
                                                         )
        process.customizedPFCandsTask.add(process.genJetsPAIReD) #Note: For gen need to add jets to the process to keep pt cuts.
        #process.customizedPFCandsTask.add(process.genJetsParticleTable)
        process.customizedPFCandsTask.add(process.genJetsPAIReD)
        
    return process