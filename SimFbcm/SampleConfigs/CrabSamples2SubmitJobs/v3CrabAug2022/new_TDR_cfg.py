import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing



options = VarParsing.VarParsing ('analysis')
options.register ('pu',
                  1.5, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "number of pile up events")
options.register ('aging',
                  0, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.float,          # string, int, or float
                  "number of pile up events")
options.parseArguments()

from Configuration.Eras.Era_Phase2_cff import Phase2
from Configuration.Eras.Modifier_fbcmDigi_cff import fbcmDigi
from Configuration.Eras.Modifier_OnlyfbcmDigi_cff import OnlyfbcmDigi

process = cms.Process("Demo",Phase2, fbcmDigi)

# import of standard configurations
#process.load('Configuration.StandardSequences.Services_cff')
#process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
# process.load("FWCore.MessageService.MessageLogger_cfi")

#process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D80Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D80_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic50ns13TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/FBCMV2/NuGun/FBCMNuGunPU50/210818_215438/0000/GEN_SIM_DIGI_99.root'
                )
                            )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('histodemonewTDR.root')
                                   )

process.demo = cms.EDAnalyzer('SimHitAnalyzer_TDR',
 RHU_InterestedHitBins = cms.vint32(0), # cms.vint32(0,1), first and last elements are included, higly depends on "BinOffset" in SiPadFrontEndParameters_cfi.py
#  InstanceNameTags = cms.vstring(
# 'Vtsh30mV', 
# 'Vtsh60mV',
# 'Vtsh90mV',
# 'Vtsh120mV'), # provide a list of valid instanse 
TreeName = cms.string( 'PU{0}'.format(options.pu)),
FbcmDigiTag = cms.string('simFbcmDigis'),
SubdetName = cms.string('FBCMHits'),
simHits = cms.InputTag("g4SimHits", "FBCMHits")
                                   )
  # tracks    = cms.untracked.InputTag('generalTracks'),
  # trackPtMin = cms.double(0.3),
  # trackEtaMin = cms.double(-2.4),
  # trackEtaMax = cms.double(2.4)

process.p = cms.Path(process.demo)
