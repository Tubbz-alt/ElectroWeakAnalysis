import FWCore.ParameterSet.Config as cms


from ElectroWeakAnalysis.VPlusJets.WenuCollectionsPAT_cfi import looseElectrons
from ElectroWeakAnalysis.VPlusJets.WenuCollectionsPAT_cfi import looseMuons


##########################################################################

## Apply loose PileUp PF jet ID
#ak5PFnoPUJets = cms.EDProducer("PATPuJetIdSelector",
#    src = cms.InputTag( "selectedPatJetsPFlow" ),
#    idLabel = cms.string("loose"),
#    valueMapLabel = cms.string("puJetMvaChs")
#)

# Apply loose PF jet ID
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
ak5PFGoodJets = cms.EDFilter("PFJetIDSelectionFunctorFilter",
     filterParams = pfJetIDSelector.clone(),
     src = cms.InputTag("goodPatJetsPFlow"),
     filter = cms.bool(True)
)

##-------- Remove electrons and muons from jet collection ----------------------
ak5PFJetsClean = cms.EDProducer("PFPATJetCleaner",
    srcJets = cms.InputTag("ak5PFGoodJets"),
    module_label = cms.string(""),
    srcObjects = cms.VInputTag(cms.InputTag("looseElectrons"),cms.InputTag("looseMuons")),
    deltaRMin = cms.double(0.3)
)


ak5PFJetsLooseIdAll = cms.EDFilter("PATJetRefSelector",
    src = cms.InputTag("ak5PFJetsClean"),
    cut = cms.string('pt > 20.0')
)

ak5PFJetsLooseId = cms.EDFilter("PATJetRefSelector",
    src = cms.InputTag("ak5PFJetsClean"),
    cut = cms.string('pt > 20.0 && abs(eta) < 2.4')
)

ak5PFJetsLooseIdVBFTag = cms.EDFilter("PATJetRefSelector",
    src = cms.InputTag("ak5PFJetsClean"),
    cut = cms.string('pt > 20.0 && abs(eta) > 2.4 && abs(eta) < 9.9')
)

##########################################
## Filter to require at least two jets in the event
RequireTwoJets = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(2),
    maxNumber = cms.uint32(100),
    src = cms.InputTag("ak5PFJetsLooseId"),
)
##########################################
## Filter to require at least one jets in the event
RequireOneJets = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(100),
    src = cms.InputTag("ak5PFJetsLooseId"),
)


############################################
#B2GPFJetPath = cms.Sequence( ak5PFnoPUJets + ak5PFGoodJets + ak5PFJetsClean + ak5PFJetsLooseId + ak5PFJetsLooseIdAll + ak5PFJetsLooseIdVBFTag + RequireOneJets )
B2GPFJetPath = cms.Sequence( ak5PFGoodJets + ak5PFJetsClean + ak5PFJetsLooseId + ak5PFJetsLooseIdAll + ak5PFJetsLooseIdVBFTag + RequireOneJets )

