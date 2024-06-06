import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("ZeeDumper")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff") # gives deprecated message in 80X but still runs
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,'130X_mcRun3_2023_realistic_v14','')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1 )
                                                                       
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),                       
    fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch//store/mc/Run3Summer23MiniAODv4/DYto2L-4Jets_MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/MINIAODSIM/130X_mcRun3_2023_realistic_v14-v1/70002/a8a27830-f9cc-4e85-8c93-df6864ab0e32.root'
    ),
    secondaryFileNames = cms.untracked.vstring()
) 

process.load('ScaleAndSmearingTools.Dumper.Zee_dumper_MINIAOD_cfi')
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("output.root")
)

process.output = cms.OutputModule("PoolOutputModule",
                                   splitLevel = cms.untracked.int32(0),
                                   outputCommands = cms.untracked.vstring("keep *"),
                                   fileName = cms.untracked.string("miniAOD.root")
)


from Geometry.CaloEventSetup.CaloGeometryBuilder_cfi import *
CaloGeometryBuilder.SelectedCalos = ['HCAL', 'ZDC', 'EcalBarrel', 'EcalEndcap', 'EcalPreshower', 'TOWER']

process.eleNewEnergies_step = cms.Path(process.eleNewEnergiesProducer)
process.dumper_step = cms.Path(process.zeedumper)
process.output_step = cms.EndPath(process.output)

#process.schedule = cms.Schedule(process.eleNewEnergies_step,process.dumper_step,process.output_step)
process.schedule = cms.Schedule(process.eleNewEnergies_step,process.dumper_step)




