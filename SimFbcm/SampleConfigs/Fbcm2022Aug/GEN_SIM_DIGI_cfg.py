# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: --evt_type SingleNuE10_cfi -s GEN,SIM,DIGI --mc --fileout file:GEN_SIM_DIGI.root --conditions auto:phase2_realistic --pileup_input file:MinBias_14TeV_pythia8_TuneCUETP8M1_GEN_SIM.root --pileup AVE_200_BX_25ns,{'B':(-3,3),'N':1.5} --era Phase2,fbcmDigi,OnlyfbcmDigi --datatier GEN-SIM-DIGI-RAW --geometry Extended2026D81 --eventcontent FEVTDEBUG --python_filename GEN_SIM_DIGI_cfg.py --customise SimFbcm/SiPadDigitizer/aging.no_aging,Configuration/DataProcessing/Utils.addMonitoring --nThreads 2 -n 2 --no_exe
import FWCore.ParameterSet.Config as cms
instanceName= "Vts25mVPreAmpModified"
import os
def getInputFileList(baseDirectory, beginsWithTheseChars ): # without "/" at the end
    # baseDirectory ='/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/FBCM/Aug2022Workshop/MinBias/FBCMV2MinBias/220820_202455'
    # baseDir = "/".join(baseDirectory.split('/')[:-1]) + "/"
    baseDir = baseDirectory + "/"
    storeDir = "/" + "/".join(baseDir.split('/')[3:])
    subDirs = os.listdir(baseDir)
    minBiasFiles = []
    for folder in subDirs:
        minBiasDirectory = baseDir + folder
        filesinDirectory = [storeDir + folder + "/" + f for f in os.listdir(minBiasDirectory) if f[:len(beginsWithTheseChars)] == beginsWithTheseChars]
        minBiasFiles = minBiasFiles + filesinDirectory
    return minBiasFiles


baseDirectory ='/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/FBCM/Aug2022Workshop/MinBias/FBCMV2MinBias/220820_202455'
minBiasFiles = getInputFileList(baseDirectory, "MinBias" )
minBiasFiles = [minBiasFiles[1]]


from Configuration.Eras.Era_Phase2_cff import Phase2
from Configuration.Eras.Modifier_fbcmDigi_cff import fbcmDigi
from Configuration.Eras.Modifier_OnlyfbcmDigi_cff import OnlyfbcmDigi

process = cms.Process('DIGI',Phase2,fbcmDigi,OnlyfbcmDigi)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D81Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D81_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic50ns13TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(1)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('SingleNuE10_cfi nevts:2'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# # Output definition
# process.TFileService = cms.Service("TFileService",
# fileName = cms.string('histodemo.root')
# )

# process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    # SelectEvents = cms.untracked.PSet(
        # SelectEvents = cms.vstring('generation_step')
    # ),
    # dataset = cms.untracked.PSet(
        # dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        # filterName = cms.untracked.string('')
    # ),
    # fileName = cms.untracked.string('file:GEN_SIM_DIGI.root'),
    # outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    # splitLevel = cms.untracked.int32(0)
# )

# Additional output definition

# Other statements
process.mix.input.nbPileupEvents.averageNumber = cms.double(1.5)
process.mix.bunchspace = cms.int32(25)
process.mix.minBunch = cms.int32(-3)
process.mix.maxBunch = cms.int32(3)
process.mix.input.fileNames = cms.untracked.vstring(['file:MinBias_14TeV_pythia8_TuneCUETP8M1_GEN_SIM.root'])
#process.mix.input.fileNames = cms.untracked.vstring(minBiasFiles)

process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.generator = cms.EDProducer("FlatRandomEGunProducer",
    AddAntiParticle = cms.bool(False),
    PGunParameters = cms.PSet(
        MaxE = cms.double(10.01),
        MaxEta = cms.double(2.5),
        MaxPhi = cms.double(3.14159265359),
        MinE = cms.double(9.99),
        MinEta = cms.double(-2.5),
        MinPhi = cms.double(-3.14159265359),
        PartID = cms.vint32(12)
    ),
    Verbosity = cms.untracked.int32(0),
    firstRun = cms.untracked.uint32(1),
    psethack = cms.string('single Nu E 10')
)
#########################


process.FbcmNtuple = cms.EDAnalyzer('FbcmNtuplizer_v4',
                                    # FbcmDigiTag = cms.InputTag("simFbcmDigis", instanceName),
                                    # FbcmDigiTag = cms.InputTag("simFbcmDigis", "SiPad"), # study another instance
                                    RHU_InterestedHitBins = cms.vint32(0), # cms.vint32(0,1), first and last elements are included, higly depends on "BinOffset" in SiPadFrontEndParameters_cfi.py
                                    TreeName = cms.string( 'PU{0}'.format(200)),
                                    InstanceNameTags = cms.vstring(
                                    'SiPadWithTimewalk','SiPadNoTimewalk',
                                    'SiPad25mVtsh', 'SiPad100mVtsh',
                                    'SiPadRfa38k154') # provide a list of valid instanse
                                    )

outFName = 'outFbcm2022_pu{0}.root'.format( '200' )
print("output is saved in {0}".format( outFName ) )
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outFName),
                                   closeFileFast = cms.untracked.bool(True) )

process.nTuple_step = cms.EndPath(process.FbcmNtuple)
#########################


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
# process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
#process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.endjob_step,process.FEVTDEBUGoutput_step)
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.endjob_step, process.nTuple_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(2)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path).insert(0, process.generator)

# customisation of the process.

process.mix.digitizers.SiPad.InstanceName = cms.string(instanceName)


# Automatic addition of the customisation function from SimFbcm.SiPadDigitizer.aging
from SimFbcm.SiPadDigitizer.aging import no_aging 

#call to customisation function no_aging imported from SimFbcm.SiPadDigitizer.aging
process = no_aging(process)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions


# Customisation from command line

SiPad1p=process.mix.digitizers.SiPad.clone()
SiPad2p=SiPad1p.clone()
SiPad3p=SiPad1p.clone()
SiPad4p=SiPad1p.clone()
SiPad5p=SiPad1p.clone()

SiPad1p.InstanceName = cms.string('SiPadWithTimewalk')
SiPad2p.InstanceName = cms.string('SiPadNoTimewalk')
SiPad3p.InstanceName = cms.string('SiPad25mVtsh')
SiPad4p.InstanceName = cms.string('SiPad100mVtsh')
SiPad5p.InstanceName = cms.string('SiPadRfa38k154')


SiPad1p.SiPadFrontEndParam[0].ApplyTimewalk =cms.bool(True)
SiPad2p.SiPadFrontEndParam[0].ApplyTimewalk =cms.bool(False)
SiPad3p.SiPadFrontEndParam[0].FE2022ASIC.ComparatorThreshold = cms.double(25.0)
SiPad4p.SiPadFrontEndParam[0].FE2022ASIC.ComparatorThreshold = cms.double(100.0)
SiPad5p.SiPadFrontEndParam[0].FE2022ASIC.R12 = cms.double(38.154)


process.theDigitizers = cms.PSet(
    SiPad1=SiPad1p,
    SiPad2=SiPad2p,
    SiPad3=SiPad3p,
    SiPad4=SiPad4p,
    SiPad5=SiPad5p,
)
process.theDigitizersValid = cms.PSet(
    SiPad1=SiPad1p,
    SiPad2=SiPad2p,
    SiPad3=SiPad3p,
    SiPad4=SiPad4p,
    SiPad5=SiPad5p,
)

process.mix.digitizers = cms.PSet(
    SiPad1=SiPad1p,
    SiPad2=SiPad2p,
    SiPad3=SiPad3p,
    SiPad4=SiPad4p,
    SiPad5=SiPad5p,
)



# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
