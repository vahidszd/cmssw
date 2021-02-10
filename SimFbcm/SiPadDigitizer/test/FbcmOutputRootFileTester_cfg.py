import FWCore.ParameterSet.Config as cms

process = cms.Process("RootTester")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/FBCM/NuGun/FBCMNuGunPU100/210130_220345/0000/GEN_SIM_DIGI_1.root'
    )
)

process.Tester = cms.EDAnalyzer('FbcmOutputRootFileTester',
				FbcmDigiTag = cms.InputTag("simFbcmDigis", "SiPad"),
				#FbcmDigiTag = cms.InputTag("simFbcmDigis"),
)

process.p = cms.Path(process.Tester)
