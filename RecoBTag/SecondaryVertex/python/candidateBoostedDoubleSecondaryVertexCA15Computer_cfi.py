import FWCore.ParameterSet.Config as cms

from RecoBTag.SecondaryVertex.trackSelection_cff import *

candidateBoostedDoubleSecondaryVertexCA15Computer = cms.ESProducer("CandidateBoostedDoubleSecondaryVertexESProducer",
    trackSelectionBlock,
    beta = cms.double(1.0),
    R0 = cms.double(1.5),
    maxSVDeltaRToJet = cms.double(1.0),
    useCondDB = cms.bool(False),
    weightFile = cms.FileInPath('RecoBTag/SecondaryVertex/data/BoostedDoubleSV_CA15_BDT_v4_temp.weights.xml.gz'),
    useGBRForest = cms.bool(True),
    useAdaBoost = cms.bool(False),
    trackPairV0Filter = cms.PSet(k0sMassWindow = cms.double(0.03)),
    SubJetName = cms.untracked.string('AK4caPFJetsTrimmedCHS__SubJets'),
    CSVbtagSubJetName = cms.untracked.string('AK4CombinedInclusiveSecondaryVertexV2BJetTagsSJCHS')	
)

candidateBoostedDoubleSecondaryVertexCA15Computer.trackSelection.jetDeltaRMax = cms.double(1.5)
