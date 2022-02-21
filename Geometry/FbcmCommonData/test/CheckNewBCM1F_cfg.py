import FWCore.ParameterSet.Config as cms

#look towards the IP, countercoloclock wise, first is right lower
process = cms.Process("DUMP")
#process.load("Geometry.FbcmCommonData.fbcmGeometryXML_cfi")
#process.load("Geometry.CMSCommonData.cmsExtendedGeometry2026D60XML_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.XMLIdealGeometryESSource = cms.ESSource("XMLIdealGeometryESSource",
									geomXMLFiles = cms.vstring(     
                                        'Geometry/FbcmCommonData/test/pltbcm.xml',
                                        'Geometry/FbcmCommonData/test/bcm1f.xml',
                                        #'Geometry/FbcmCommonData/test/plt.xml',
                                        'Geometry/FbcmCommonData/test/cmsFbcm.xml',
										'Geometry/FbcmCommonData/test/fbcm.xml',
                                        'Geometry/FbcmCommonData/test/FbcmSpecs.xml',
										#'Geometry/CMSCommonData/data/materials.xml',
                                        'Geometry/CMSCommonData/data/materials/2021/v2/materials.xml',
										'Geometry/CMSCommonData/data/rotations.xml',
										'Geometry/FbcmCommonData/test/cmsextent.xml',
                                        #'Geometry/CMSCommonData/data/extend/v3/cmsextent.xml',
										'Geometry/CMSCommonData/data/cmsMother.xml',
										'Geometry/CMSCommonData/data/cmsTracker.xml',
										'Geometry/CMSCommonData/data/cavernData/2021/v1/cavernData.xml',
										#'Geometry/CMSCommonData/data/cms/2026/v3/cms.xml',
                                        'Geometry/CMSCommonData/data/cms/2021/v3/cms.xml',
										#'Geometry/TrackerCommonData/data/PhaseII/trackerParameters.xml',
										#'Geometry/TrackerCommonData/data/pixfwdCommon.xml',
										#'Geometry/TrackerCommonData/data/PhaseII/TiltedTracker613_MB_2019_04/pixfwd.xml',
										#'Geometry/TrackerCommonData/data/PhaseII/TiltedTracker613_MB_2019_04/pixbar.xml',
										#'Geometry/TrackerCommonData/data/trackermaterial.xml',
                                        'Geometry/TrackerCommonData/data/trackermaterial/2021/v1/trackermaterial.xml',
										#'Geometry/TrackerCommonData/data/PhaseII/TiltedTracker613_MB_2019_04/tracker.xml',
                                        'Geometry/TrackerCommonData/data/Run2/tracker.xml',
                                        'Geometry/TrackerCommonData/data/pixfwdMaterials/2021/v2/pixfwdMaterials.xml',
                                        'Geometry/TrackerCommonData/data/pixfwdCommon.xml',
                                        'Geometry/TrackerCommonData/data/PhaseI/v2/pixfwd.xml',
                                        'Geometry/TrackerCommonData/data/PhaseI/v2/pixfwdCylinder.xml',
                                        'Geometry/TrackerCommonData/data/PhaseI/v2/pixfwdDisks.xml',
                                        'Geometry/TrackerCommonData/data/PhaseI/v2/pixfwdDisks.xml',
                                        'Geometry/TrackerCommonData/data/PhaseI/v2/pixfwdSupportRingParameters.xml',
                                        'Geometry/TrackerCommonData/data/PhaseI/trackerpixfwd.xml',
										#'Geometry/TrackerCommonData/data/PhaseII/TiltedTracker613_MB_2019_04/pixel.xml',
										#'Geometry/TrackerCommonData/data/PhaseII/TiltedTracker404/trackerfwd.xml',
										#'Geometry/TrackerRecoData/data/PhaseII/TiltedTracker613_MB_2019_04/trackerRecoMaterial.xml'

									),
									rootNodeName = cms.string('cms:OCMS')
								)



process.MessageLogger = cms.Service("MessageLogger",
                                    debugModules = cms.untracked.vstring('*'),
                                    destinations = cms.untracked.vstring('cout'),
                                    cout = cms.untracked.PSet
                                    (
    threshold = cms.untracked.string('DEBUG'),
    noLineBreaks = cms.untracked.bool(True)
    )
                                    )

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.add_(cms.ESProducer("TGeoMgrFromDdd",
        verbose = cms.untracked.bool(False),
        level   = cms.untracked.int32(14)
))

process.dump = cms.EDAnalyzer("DumpSimGeometry")

process.p = cms.Path(process.dump)



