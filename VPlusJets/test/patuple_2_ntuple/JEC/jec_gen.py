import FWCore.ParameterSet.Config as cms
process = cms.Process("jectxt")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# define your favorite global tag
#process.GlobalTag.globaltag = 'GR_R_53_V10::All'
#process.GlobalTag.globaltag = 'START53_V15::All'
#process.GlobalTag.globaltag = 'START53_V7E::All'
#process.GlobalTag.globaltag = 'FT_53_V10_AN3::All'
#process.GlobalTag.globaltag = 'GR_P_V42_AN4::All'
process.GlobalTag.globaltag = 'START53_V7G::All'
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
process.source = cms.Source("EmptySource")
process.readAK5PFchs    = cms.EDAnalyzer('JetCorrectorDBReader',  # below is the communication to the database 
        payloadName    = cms.untracked.string('AK5PFchs'),
        #payloadName    = cms.untracked.string('AK5PF'),
        # this is used ONLY for the name of the printed txt files. You can use any name that you like, 
        # but it is recommended to use the GT name that you retrieved the files from.
        #globalTag      = cms.untracked.string('GR_R_53_V10'),  
        #globalTag      = cms.untracked.string('START53_V15'),  
        #globalTag      = cms.untracked.string('START53_V7E'),  
        #globalTag      = cms.untracked.string('FT_53_V10_AN3'),  
        #globalTag      = cms.untracked.string('GR_P_V42_AN4'),  
        globalTag      = cms.untracked.string('START53_V7G'),  
        printScreen    = cms.untracked.bool(False),
        createTextFile = cms.untracked.bool(True)
        )
#process.readAK5Calo = process.readAK5PF.clone(payloadName = 'AK5Calo')
#process.readAK5JPT = process.readAK5PF.clone(payloadName = 'AK5JPT')
#process.p = cms.Path(process.readAK5PF * process.readAK5Calo * process.readAK5JPT)
process.p = cms.Path(process.readAK5PFchs)
