import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os


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
#process.load('Configuration.Geometry.GeometryExtended2026D81Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D81_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic50ns13TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )




##################################################################
#     Getting proper Directory for PU                            #
##################################################################
def get_files_in_directory(directory):
    file_list = []
    for root, dirs, files in os.walk(directory):
        for filename in files:
            file_list.append(os.path.join(root, filename))
    return file_list

# Provide the directory path you want to read files from

#directory_path = '/eos/home-v/vsedighz/vahidtest/NuGun/FbcmMultiInstanceNuGunPU0p5'

# Call the function to get the list of files in the directory and subdirectories


print(options.pu)
if options.pu == '0.5':
    directory_path = '/eos/home-v/vsedighz/vahidtest/FBCM_with_TDR_GEN_SIM_DIGI/NuGun/FbcmMultiInstanceNuGunPU0p5'
elif options.pu == '1':
    directory_path = '/eos/home-v/vsedighz/vahidtest/FBCM_with_TDR_GEN_SIM_DIGI/NuGun//FbcmMultiInstanceNuGunPU1'
elif options.pu == '1.5':
    directory_path = '/eos/home-v/vsedighz/vahidtest/FBCM_with_TDR_GEN_SIM_DIGI/NuGun/FbcmMultiInstanceNuGunPU1p5'
elif options.pu == 1.5:
    directory_path = '/eos/home-v/vsedighz/vahidtest/FBCM_with_TDR_GEN_SIM_DIGI/NuGun/FbcmMultiInstanceNuGunPU1p5'
elif options.pu == '2':
    directory_path = '/eos/home-v/vsedighz/vahidtest/FBCM_with_TDR_GEN_SIM_DIGI/NuGun/FbcmMultiInstanceNuGunPU2'
elif options.pu == '10':
    directory_path = '/eos/home-v/vsedighz/vahidtest/FBCM_with_TDR_GEN_SIM_DIGI/NuGun/FbcmMultiInstanceNuGunPU10'
elif options.pu == '30':
    directory_path = '/eos/home-v/vsedighz/vahidtest/FBCM_with_TDR_GEN_SIM_DIGI/NuGun/FbcmMultiInstanceNuGunPU30'
elif options.pu == '50':
    directory_path = '/eos/home-v/vsedighz/vahidtest/FBCM_with_TDR_GEN_SIM_DIGI/NuGun/FbcmMultiInstanceNuGunPU50'
elif options.pu == '100':
    directory_path = '/eos/home-v/vsedighz/vahidtest/FBCM_with_TDR_GEN_SIM_DIGI/NuGun/FbcmMultiInstanceNuGunPU100'
elif options.pu == '140':
    directory_path = '/eos/home-v/vsedighz/vahidtest/FBCM_with_TDR_GEN_SIM_DIGI/NuGun/FbcmMultiInstanceNuGunPU140'
elif options.pu == '200':
    directory_path = '/eos/home-v/vsedighz/vahidtest/FBCM_with_TDR_GEN_SIM_DIGI/NuGun/FbcmMultiInstanceNuGunPU200'
else:
    print("Invalid options.pu value!")

print (directory_path)
files =  get_files_in_directory(directory_path)
# for f in files:
#     print(f) 

file_paths = []
for root, directories, files in os.walk(directory_path):
    for filename in files:
        file_path = os.path.join(root, filename)
        file_paths.append('file:' + file_path)

print(file_paths)
#print(file)
# file =  ",".join(["'file:{}'".format(path) for path in files])
# print (file)    




#print ("type of root_files_200 is : ",type(root_files_200))
#print ("type of file_paths is : " , type (file_paths))

process.source = cms.Source("PoolSource",
                            fileNames =cms.untracked.vstring(file_paths
 #           'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/FBCM/NuGun/FBCMNuGunPU100/210130_220345/0000/GEN_SIM_DIGI_99.root'#(TDR files)
#'file:/eos/home-v/vsedighz/lasttest/NuGun/FbcmMultiInstanceNuGunPU50/230529_051020/0000/GEN_SIM_DIGI_1.root','file:/eos/home-v/vsedighz/lasttest/NuGun/FbcmMultiInstanceNuGunPU50/230529_051020/0000/GEN_SIM_DIGI_2.root','file:/eos/home-v/vsedighz/lasttest/NuGun/FbcmMultiInstanceNuGunPU50/230529_051020/0000/GEN_SIM_DIGI_3.root','file:/eos/home-v/vsedighz/lasttest/NuGun/FbcmMultiInstanceNuGunPU50/230529_051020/0000/GEN_SIM_DIGI_4.root','file:/eos/home-v/vsedighz/lasttest/NuGun/FbcmMultiInstanceNuGunPU50/230529_051020/0000/GEN_SIM_DIGI_5.root'
#'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/FBCM/Aug2022Workshop/GEN_SIM_DIGI_Validation/NuGun/FbcmMultiInstanceNuGunPU200/230611_161305/0000/GEN_SIM_DIGI_99.root'
                           )
                            )
output_folder = '/eos/home-v/vsedighz/vahidtest/FBCM_with_TDRFE_hists'
process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string(os.path.join(output_folder,'histodemonew_pu{0}.root'.format(options.pu)))
                                   )

process.demo = cms.EDAnalyzer('SimHitAnalyzer',
 InstanceNameTags = cms.vstring(
'Vtsh30mV', 
'Vtsh60mV',
'Vtsh90mV',
'Vtsh120mV',
#'Vtsh1000mV'
                                ), # provide a list of valid instanse 
TreeName = cms.vstring( 'PU{0}'.format(options.pu)),
hitsProducer = cms.string('g4SimHits'),
SubdetName = cms.string('FBCMHits'),
simHits = cms.InputTag("g4SimHits", "FBCMHits")
                                   )
  # tracks    = cms.untracked.InputTag('generalTracks'),
  # trackPtMin = cms.double(0.3),
  # trackEtaMin = cms.double(-2.4),
  # trackEtaMax = cms.double(2.4)

process.p = cms.Path(process.demo)
