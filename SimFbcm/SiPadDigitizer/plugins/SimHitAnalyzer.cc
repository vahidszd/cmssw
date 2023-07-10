#include <memory>
#include <string>
#include <set>
#include <type_traits>
#include <typeinfo>
#include <vector>
#include <iostream>
#include <fstream>
//#include <format>

#include "SiPadDigitizer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/FbcmDigi/interface/SiPadDigiData.h"
//#include "DataFormats/FbcmDigi/interface/PSimHitInfo.h"
#include "DataFormats/FbcmDetId/interface/FbcmDetId.h"

#include "Geometry/Records/interface/FbcmGeometryRecord.h"
#include "Geometry/FbcmGeometry/interface/FbcmGeometry.h"
#include "DataFormats/FbcmDigi/interface/SiPadDigiData.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TTree.h"
#include <cmath>
//#include <matplot.h>



#define MAX_SIMHIT_SIZE 100
#define MAX_DIGHIT_SIZE 100
#define MAX_BXSlot_SIZE 100
#define AREA_GROUP_LEN 8

#define TEST_IND 7

using namespace std;

class SimHitAnalyzer : public edm::EDAnalyzer {
public:

  explicit SimHitAnalyzer (const edm::ParameterSet& iConfig);
  ~SimHitAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
 void accumulateSiPadHits(edm::Handle<std::vector<PSimHit> >, size_t globalSimHitIndex);
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

 virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  void endRun(edm::Run const&, edm::EventSetup const&) override;
  //----------------- Data Member---------------------------
    std::vector< edm::EDGetTokenT< edm::DetSetVector<SiPadDigiData> > > TokenTagsDigi_;	
  //const edm::ProducesCollector producesCollector  ;

  const std::string hitsProducer_;
  std::vector<std::string> Directory;
  std::vector<std::string> instanceNames;
 const std::string SubdetName_; 
  std::vector<edm::EDGetTokenT< edm::PSimHitContainer>> TokenTags_;  
  edm::EDGetTokenT<edm::PSimHitContainer> token_;
  edm::EDGetTokenT<std::vector<PSimHit>> token; //consumes<std::vector<PSimHit>>(edm::InputTag const&);
  std::unique_ptr<SiPadDigitizerAlgorithm> SiPadDigiAlgo;  
  edm::ESHandle<MagneticField> pSetup; 
  std::map<unsigned int, FbcmSiPadGeom const*> SiPadsIdGeomMap;
  std::map<unsigned int, float> rhoarea;
  edm::ESHandle<FbcmGeometry> theFbcmGeom;
  std::vector<TFileDirectory> tagDirectories;
std::vector<TFileDirectory> tagDirectoriesDigi;
  TFileDirectory Geom;
  edm::Handle<edm::PSimHitContainer> hsimHits;

  int numberOfSensorGroups;    
  
 
	std::map<std::string, size_t> crossingSimHitIndexOffset_;
  TH1* NofEvents;
  std::vector< TH1* > hist_NHits;
  std::vector< TH1* > hNEvents;
  std::vector< TH2* > RhoArea;
  std::vector< TH1* > PeakAmple;
  TH2*  RhoAreaGeom;
  std::vector< TH2* >  RhoAreaNormalizedNew;
std::vector< TH2* >  DigiHitsPerArea;
  std::vector < TH2* > DigiHitsPerSensor;
  TH2*  AreaGeom ;
  
  std::vector< TTree* > theTree;
  std::vector< TTree* > FixedValuesTree;
float SensorArea;
  float SensorRho;
  
  int SensorGroupIndex; // collected for each hit. 
  //int numberOfSensorGroups; // maximum number of SiGroupIndex used in the digi step. 
  int nbrOfDiesPerRing; 
    
  int SensorZ_end;
  
  //int NumOfSensors;
  
int nTotalRhuDigi;
  int nSimParticles;
  int BxSlotCnt; 
  int nValidDigiToAs;
  int nValidDigiToTs;
  int SimPdgIds[MAX_SIMHIT_SIZE];
  float SimPts[MAX_SIMHIT_SIZE];
  float SimCharges[MAX_SIMHIT_SIZE];
  float SimTofs[MAX_SIMHIT_SIZE];
  float SimTof_perBx[MAX_SIMHIT_SIZE];

  int DigiHitStatus[MAX_BXSlot_SIZE];
  int nDigiHits[MAX_BXSlot_SIZE]; 
float DigiToAs[MAX_DIGHIT_SIZE] ;
  float DigiToTs[MAX_DIGHIT_SIZE];
  float PeakAmpl[MAX_DIGHIT_SIZE];
  
  
  int nBinsX;
  int nBinsY; 
  int b;
  int val;
  float totArea ;
    int totval ;
  float hitPerArea;
  //FILE *fp;
  ofstream os;
  //char fi;
  //char fileName[] = "0.txt";

}; 
//--------------------------End of Class decleration------------------------------------
//--------------------------Constructor begins-----------------------------------------
SimHitAnalyzer::SimHitAnalyzer(const edm::ParameterSet& iConfig):

 

hitsProducer_(iConfig.getParameter<std::string>("hitsProducer")), //g4SimHits
Directory(iConfig.getParameter<std::vector<std::string>>("TreeName")),
instanceNames(iConfig.getParameter< std::vector<std::string> >("InstanceNameTags")),
SubdetName_(iConfig.getParameter<std::string>("SubdetName")),//FBCMHits
token_(consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("simHits"))),
numberOfSensorGroups(0)

{
    
   edm::InputTag tag(hitsProducer_,SubdetName_); 
  //token = consumes< std::vector<PSimHit> >(tag);
  numberOfSensorGroups=0;
  
  //consumes<edm::PSimHitContainer>(tag);
  token = consumes<edm::PSimHitContainer>(tag);
  // edm::Service<TFileService> fs;

  //  Directory.clear();

    //    hist_NHits = fs->make<TH1D>("nPSimHits","nPSimHits",100,0,100); 
  edm::Service<TFileService> fs;
  //TFileDirectory subDir;  
  tagDirectories.clear();
 hist_NHits.clear();
    theTree.clear();
  FixedValuesTree.clear();
  hNEvents.clear();
  //hNSensorGroups.clear();
 int m;
  for(auto it = std::begin(instanceNames); it != std::end(instanceNames); ++it) {

    //std::cout<<"stop for instannce name "<<it->c_str();
    //std::cin>> m;
    edm::InputTag tag("simFbcmDigis", it->c_str() );
    //cout<<"The tag is "<< tag;
    TokenTagsDigi_.emplace_back( consumes< edm::DetSetVector<SiPadDigiData> >( tag )  );
  }
//     tagDirectoriesDigi.emplace_back( fs->mkdir( it->c_str() ) ); 
// FixedValuesTree.emplace_back( tagDirectoriesDigi.back().make<TTree>( "GeometryInfo" , "FixedValues" )   );
//     hNEvents.emplace_back( tagDirectoriesDigi.back().make<TH1I>("hNEvents" , "" , 1 , 0 , 1 )  );
//     RhoAreaGeom.emplace_back( tagDirectoriesDigi.back().make<TH2F>("testGeom" , "" , 70 , 7.5 , 22,100,0,0.3 )  );

//   }

  //for (auto it = std::begin(Directory); it != std::end(Directory); it++)
  // {
       //tagDirectories.emplace_back( fs->mkdir( it->c_str() ) );
  //   cout<< "treename is : " << it->c_str() << endl;
  //   hist_NHits.emplace_back(tagDirectories.back().make<TH1I>("noPsimhits","",200,0,200));
  // }

  }
SimHitAnalyzer::~SimHitAnalyzer()
{
  //
}



void
SimHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
edm::Service<TFileService> fs;
    NofEvents->Fill(0.5);
   using namespace edm;
   using namespace std;
   int nBinsX();

   //  edm::Handle<std::vector<PSimHit> > hsimHits;
   




   //edm::Handle<edm::PSimHitContainer> hsimHits;

          edm::InputTag tag(hitsProducer_, SubdetName_);
       cout << "tag in SimHitAnalyzer is :" << tag << endl;
     
   
   cout<<"we are in analyze ------------------------------$$$$%%%%%%%%%%%%%"<<endl;
   //int n;
    for (auto it = std::begin(instanceNames); it != std::end(instanceNames); it++)
      {
 edm::InputTag tag("simFbcmDigis", it->c_str() );
    cout<<"The tag is "<< tag;
    // TokenTags_.emplace_back( consumes< edm::DetSetVector<SiPadDigiData> >( tag )  ); 
        tagDirectories.emplace_back( fs->mkdir( it->c_str() ) );
        cout<< "treename is : " << it->c_str() << endl;
        //cin >> n;
      }
   //iEvent.getByLabel(tag, hsimHits);
   iEvent.getByToken(token, hsimHits);
   // for (const auto& simHit : *hsimHits) {
   //   std::cout << "Track ID: " << simHit.trackId() << std::endl;
   //   std::cout << "Energy Deposit: " << simHit.energyLoss() << std::endl;
   //   std::cout << "Time of Flight: " << simHit.timeOfFlight() << std::endl;
   //   std::cout << "Local Position: " << simHit.localPosition() << std::endl;
   // }
   //const std::vector<PSimHit> & simHits = *hsimHits.product();

   //  for (std::vector<PSimHit>::const_iterator simHit = simHits.begin(); simHit != simHits.end(); ++simHit) {
      // Extract data from individual PSimHit objects
   // cout<<"we are in second loop"<<endl;
      // cout << "we are in PSimHits loop " << endl;
      // int particleType = simHit->particleType();
      // //GlobalPoint hitPosition = simHit->entryPoint();
      // double energyDeposit = simHit->energyLoss();
      // cout << "energyDeposit is: " << energyDeposit << endl;
      // // Store the extracted data
 
      //   }
   for(long unsigned int i = 0; i < tagDirectoriesDigi.size(); i++){
  if(hsimHits.isValid()){
    //cout << "SimHIt is Valid." << endl;
    // cout << "hsimhits->size is : " << hsimHits->size() << endl;
    //RhoArea[i]->Fill(5,6);
     hNEvents[i]->Fill(hsimHits->size());

  edm::Handle< edm::DetSetVector<SiPadDigiData> > handle; 
  iEvent.getByToken(TokenTagsDigi_[i] ,handle);
  
  
  for (edm::DetSetVector<SiPadDigiData>::const_iterator itDetSet = handle->begin() ;  itDetSet < handle->end() ; ++itDetSet) {
     cout << "firs loop" << endl;
    //cout << *itDetSet;
    // each Det in fbcmDets: multiple Data, but we have only one data
    if( itDetSet->size() != 1 )
       cout << "size == " << itDetSet->size() << endl; 

    for (std::vector<SiPadDigiData>::const_iterator Digidata_it = itDetSet->begin() ; Digidata_it < itDetSet->end() ; ++Digidata_it) {
      //const SiPadDigiData Digidata= *(Digidata_it);
      //SiDieId=Digidata_it->SiliconDieIndex() ;
      //SensorGroupIndex = SiDieId % nbrOfDiesPerRing;
      cout << "we are in second loop" << endl;
      SensorGroupIndex = Digidata_it->SiGroupIndex();
      
      if ( SensorGroupIndex > numberOfSensorGroups )
          numberOfSensorGroups=SensorGroupIndex;
      
      
       SensorZ_end = Digidata_it->SideIndex();
       cout << "SensorZ_end = " << SensorZ_end << endl;
       SensorRho=Digidata_it->Radius() ;
       cout << "SensorRho = " << SensorRho << endl;
       SensorArea=Digidata_it->Area() ;
        
             RhoArea[i]->Fill(SensorRho,SensorArea);
	    
            
       cout << "SensorArea = " << SensorArea << endl;
       nSimParticles = 0 ;
      // if ( SensorArea > 0.0290){
      // 	int tm ;
      //   cin >> tm ;
      // }
      for(auto sim : Digidata_it->CahrgePsimVector()) {
      	SimPdgIds[nSimParticles] = sim.second.ParticleType() ;
      	SimPts[nSimParticles] = sim.second.Pabs() ;
      	SimTofs[nSimParticles] = sim.second.Tof() ;
      	SimCharges[nSimParticles] = sim.first;
      	SimTof_perBx[nSimParticles] = sim.second.Tof() - 25.0 * sim.second.BunchCrossing() ;

      	nSimParticles++;
	//cout << " nSimParticles = " << nSimParticles << endl;
      }


     nValidDigiToAs=0;
      nValidDigiToTs=0;
      nTotalRhuDigi=0;
      // nInterstedRhuBins=0;
      
      BxSlotCnt=0; // just to make sure starting from non-zero index if in
   for(auto hit : Digidata_it->BxSlotHitAnalysisVector() ) {	  
            DigiHitStatus[BxSlotCnt] = hit.Bx_HitStatusInt(); 
            nDigiHits[BxSlotCnt] = hit.nbrOfRecognizedHitsInBx(); 
            BxSlotCnt++; 
                    
             for (auto ToaTot : hit.TotToaVectort()) {
                PeakAmple[i]->Fill(ToaTot.PeakAmplitude());
               if ( (ToaTot.ToAToAStatusInt() == ToAStatus::FullyWithinBx || ToaTot.ToAToAStatusInt() == ToAStatus::WithinBx_LastsAfter ) ) 
               {
                   {
		     //DigiToAs[nValidDigiToAs] = ToaTot.ToA();
                      PeakAmpl[nValidDigiToAs] = ToaTot.PeakAmplitude();
		      
                      nValidDigiToAs++;
                  }
                if (ToaTot.IsToTValid())
                  DigiToTs[nValidDigiToTs++] = ToaTot.ToT();
              
                //DigiRHUs[nTotalRhuDigi++] = ToaTot.SubBxBinNumber();
                //if (std::find(rhuInterestedHitBins_.begin(), rhuInterestedHitBins_.end(),ToaTot.SubBxBinNumber())!=rhuInterestedHitBins_.end())
		//  nInterstedRhuBins++;
              
              }

             }
	     //for (auto amp : PeakAmpl)
	     //cout << "Ample = " << PeakAmpl << endl;
		
      }
 
    }
    cout << "bxslotcont = " << BxSlotCnt << endl;
   
   RhoArea[i]->Draw("colz");

 } 
     std::set<unsigned int> detIds;

    }  
   }
}


void SimHitAnalyzer::beginJob()
{
//  int nBinsX;
//   int nBinsY;
//   int N;
//   float Rho;
//   int SiDieId=0;
//   int SiDieGroupIndex=0;
//   float padX, padY , padRho;
//   float minRho = 100, maxRho = 0, minArea = 10, maxArea= 0;
//   float area;
//   vector <Double_t> binedges;
//   //iSetup.get<FbcmGeometryRecord>().get(theFbcmGeom); 
//   const std::vector<const FbcmStationGeom*> AllStatitons = theFbcmGeom->Stations();
//   const std::vector<FbcmSiPadGeom const*> allSiPadGeoms = theFbcmGeom->SiPads();
//   nbrOfDiesPerRing = AllStatitons[0]->NumOfDiesPerRing(); 
//   cout<< "we are in beginjob simhitanalyzer ===============" << endl;
//   edm::Service<TFileService> fs;
 


//  for(const auto& siPad : allSiPadGeoms)
//     {
//       std::pair<float, float> SiPadDimension = siPad->SiPadTopology().pitch();

//       padX = SiPadDimension.first;
//       padY = SiPadDimension.second;
//       area = padX*padY;
//       bool exists = find(binedges.begin(), binedges.end(), area) != binedges.end();
//       if ( !exists ){
// 	cout << " area = " << area << endl;
// 	binedges.push_back(area);
// 	cout << "biedges content  = " << binedges[binedges.size()-1] << endl;
//         cout<< "binages size = "  << binedges.size() << endl ;
// 	//	cin >> n;
//       }
 
//       // cout<< "we are in beginjob simhitanalyzer ===============" << endl;
//       // edm::Service<TFileService> fs;
//  for(auto it = std::begin(instanceNames); it != std::end(instanceNames); ++it) {

//     std::cout<<"stop for instannce name "<<it->c_str();
//     //std::cin>> m;
//     edm::InputTag tag("simFbcmDigis", it->c_str() );
//     //cout<<"The tag is "<< tag;
//     //   TokenTagsDigi_.emplace_back( consumes< edm::DetSetVector<SiPadDigiData> >( tag )  );
//     tagDirectoriesDigi.emplace_back( fs->mkdir( it->c_str() ) ); 
// FixedValuesTree.emplace_back( tagDirectoriesDigi.back().make<TTree>( "GeometryInfo" , "FixedValues" )   );
//     hNEvents.emplace_back( tagDirectoriesDigi.back().make<TH1I>("hNEvents" , "" , 1 , 0 , 1 )  );
//     RhoAreaGeom.emplace_back( tagDirectoriesDigi.back().make<TH2F>("testGeom" , "" , 70 , 7.5 , 22,100,0,0.3 )  );
// RhoAreaNormalizedNew.emplace_back( tagDirectoriesDigi.back().make<TH2F>("testGeomNew" , "" , 70 , 7.5 , 22,100,0,0.3 )  );
//  RhoArea.emplace_back( tagDirectoriesDigi.back().make<TH2F>("test" , "" , 70 , 7.5 , 22,100,0,0.3 )  );
//   }
//     }
 }
 void SimHitAnalyzer::endJob()
 {

}

void SimHitAnalyzer::beginRun(edm::Run const&, edm::EventSetup const& iSetup)
{

 

  vector <Double_t> binedges;
  float binwidthRho;
  int nBinsX;
  int nBinsY;
  int N;
  float Rho;
  int SiDieId=0;
  int SiDieGroupIndex=0;
  float padX, padY , padRho;
  float minRho = 100, maxRho = 0, minArea = 10, maxArea= 0;
  float area;

 edm::Service<TFileService> fs;

  // RhoAreaGeom.emplace_back( tagDirectoriesDigi.back().make<TH2F>("testGeom" , "" , 70 , 7.5 , 22,100,0,0.3 )  );
 cout << "we are in beginRun SimHItanalyzer&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& " << endl <<"tagDirectoriesDigi = "  << tagDirectoriesDigi.size() << endl;

 iSetup.get<FbcmGeometryRecord>().get(theFbcmGeom); 
  const std::vector<const FbcmStationGeom*> AllStatitons = theFbcmGeom->Stations();
  const std::vector<FbcmSiPadGeom const*> allSiPadGeoms = theFbcmGeom->SiPads();
 

 for(const auto& siPad : allSiPadGeoms)
    {
      cout << "we are in sipadGeometry loop" << endl;
      std::pair<float, float> SiPadDimension = siPad->SiPadTopology().pitch();
     
      padX = SiPadDimension.first;
      padY = SiPadDimension.second;
      area = padX*padY;
      // we could also store only the Area instead of padX and padY
      padRho = siPad->surface().position().perp(); 
bool exists = find(binedges.begin(), binedges.end(), area) != binedges.end();
      if ( !exists ){
		int m;
	cout << " area = " << area << endl;
	binedges.push_back(area);
	cout << "biedges content  = " << binedges[binedges.size()-1] << endl;
        cout<< "binages size = "  << binedges.size() << endl ;
		cin >> m;
      }
      //   SiDieId=siPad->id().SiliconDie();
      // SiDieGroupIndex = SiDieId % nbrOfDiesPerRing;
 
 // cout << "SiDieGroupIndex = " << SiDieGroupIndex << endl;
      cout << "padx is " << padX << endl;	  
      cout << "padY is " << padY << endl;	  
      cout << "padRho is " << padRho << endl;
      //RhoAreaGeom[i]->Fill(padRho,padX*padY);	  
        //std::cout << SiDieId <<", " << SiDieGroupIndex <<"\n"; 
      minRho = min(minRho, padRho);
      maxRho = max(maxRho,padRho);
      minArea = min(minArea, padX*padY);
      maxArea = max(maxArea, padX*padY);
      binwidthRho = (maxRho -minRho)/sqrt(minArea);
      cout << " minRho = " << minRho << endl << "maxRho = " << maxRho << endl;
      cout << "minAre = " << minArea << endl << "maxArea = " << maxArea << endl;
    }
 binedges.push_back(0.04);
 const Int_t  n = binedges.size();
  const Int_t binNo = n-1;
  Double_t BinEdges[n]  ;
 for (int i = 0 ; i<n ; i++)
  BinEdges[i] = binedges[i] ;
 double* BinEdgesNew  = &binedges[0]; 
 //cout << "size of BinEdges = " << BinEdges.size() << endl;
 // copy(binedges.begin(), binedges.end(), BinEdges); 
  for (double it :BinEdges )
    cout << " array  = " << it << endl;
  int m;
  cin >> m;
  const Int_t NBINS = 1;
  Double_t edges[NBINS+1 ] = {BinEdges[0],BinEdges[0]+0.03};

  Geom = fs->mkdir("Geom");
AreaGeom = Geom.make<TH2F>(" SensorArea" , " Sensor Area" ,Int_t(binwidthRho)  , minRho ,maxRho,binNo,BinEdgesNew   );   
NofEvents = fs->make<TH1I>("No of Events ","NofEvents",1,0,1); 
 RhoAreaGeom = Geom.make<TH2F>("RhoAreaGeom" , "Number of sensors in a Radius" , Int_t(binwidthRho)  , minRho ,maxRho,binNo,BinEdgesNew   );
 for(auto it = std::begin(instanceNames); it != std::end(instanceNames); ++it) {

   string label = "DigiHitPerSensor_" + string(it->c_str()) ;
   string label1 ="nOfDigiHits" + string(it->c_str());
string title = " Number of DigiHits Per Sensor Per Event_" + string(it->c_str()) ;
string title1 = "Total Number of DigiHits_" + string(it->c_str()) ;
    std::cout<<"stop for instannce name "<<it->c_str();
    //std::cin>> m;
    edm::InputTag tag("simFbcmDigis", it->c_str() );
    //cout<<"The tag is "<< tag;
    //   TokenTagsDigi_.emplace_back( consumes< edm::DetSetVector<SiPadDigiData> >( tag )  );
    tagDirectoriesDigi.emplace_back( fs->mkdir( it->c_str() ) );
    

    hNEvents.emplace_back( tagDirectoriesDigi.back().make<TH1I>("hNEvents" , "" , 1 , 0 , 1 )  );
   
    RhoAreaNormalizedNew.emplace_back( tagDirectoriesDigi.back().make<TH2F>("RhoAreaNormalizedNew" , "RhoAreaNormalizedNew" , Int_t(binwidthRho)  , minRho ,maxRho,binNo,BinEdgesNew  )  );
    RhoArea.emplace_back( tagDirectoriesDigi.back().make<TH2F>(label1.c_str() ,title1.c_str() ,Int_t(binwidthRho)  , minRho ,maxRho,binNo,BinEdgesNew )  );
    DigiHitsPerSensor.emplace_back( tagDirectoriesDigi.back().make<TH2F>(label.c_str() , label.c_str() ,Int_t(binwidthRho)  , minRho ,maxRho,binNo,BinEdgesNew )  );
 DigiHitsPerArea.emplace_back( tagDirectoriesDigi.back().make<TH2F>(" DigiHitsPerArea" , " Number of DigiHits Per mm^{2}" ,Int_t(binwidthRho)  , minRho ,maxRho,binNo,BinEdgesNew )  );
 PeakAmple.emplace_back( tagDirectoriesDigi.back().make<TH1F>(" PeakAmple" , " PeakAmple",100,0,1000 )  );
 
}

  // nbrOfDiesPerRing = AllStatitons[0]->NumOfDiesPerRing(); 
  
  cout << " Number of stations = " << nbrOfDiesPerRing << endl;
  for(long unsigned int i=0; i < tagDirectoriesDigi.size(); i++){
    //    int n;
    //int SiDieId=0;
    //int SiDieGroupIndex=0;

 
  for(const auto& siPad : allSiPadGeoms)
    {
     
      std::pair<float, float> SiPadDimension = siPad->SiPadTopology().pitch();
     
      padX = SiPadDimension.first;
      padY = SiPadDimension.second;
     
     

      // we could also store only the Area instead of padX and padY
     
      // we could also store only the Area instead of padX and padY
      padRho = siPad->surface().position().perp(); 
bool exists = find(binedges.begin(), binedges.end(), area) != binedges.end();
     
      //RhoAreaGeom[i]->Fill(padRho,padX*padY);	  
        //std::cout << SiDieId <<", " << SiDieGroupIndex <<"\n"; 
      minRho = min(minRho, padRho);
      maxRho = max(maxRho,padRho);
      minArea = min(minArea, padX*padY);
      maxArea = max(maxArea, padX*padY);
      cout << " minRho = " << minRho << endl << "maxRho = " << maxRho << endl;
      cout << "minAre = " << minArea << endl << "maxArea = " << maxArea << endl;
     
     
      RhoAreaNormalizedNew[i]->Fill(padRho,padX*padY);
    }

 
  }
for(const auto& siPad : allSiPadGeoms)
    {
      cout << "we are in sipadGeometry loop" << endl;
      std::pair<float, float> SiPadDimension = siPad->SiPadTopology().pitch();
     
      padX = SiPadDimension.first;
      padY = SiPadDimension.second;
      area = padX*padY;
      // we could also store only the Area instead of padX and padY
      padRho = siPad->surface().position().perp();
      RhoAreaGeom->Fill(padRho,padX*padY); 
    }
nBinsX = RhoAreaGeom->GetNbinsX();
 nBinsY = RhoAreaGeom->GetNbinsY();
       for (int binx  = 0 ; binx <=nBinsX-1 ; binx ++)
          for ( int biny = 0 ; biny <=nBinsY-1; biny ++){
            int  bin = RhoAreaGeom->GetBin(binx,biny);
	    cout << "BinEdgesNew is " << BinEdgesNew[biny] << endl;
	    int newm;
	    // cin >> newm;
	    AreaGeom->SetBinContent(binx,biny,BinEdgesNew[biny]);
		   }

}

  //  edm::InputTag tag(hitsProducer_);
  //consumes<edm::PSimHitContainer>(tag);

  //   std::cout << "vahid" <<endl;
//   edm::Handle<std::vector<PSimHit> > hSimHit;
//    if(hSimHit.isValid()){
//   std::vector<PSimHit> const& simHits = *hSimHit.product();
//   std::set<unsigned int> detIds;
//   int i = 0;
//    for (auto it = simHits.begin(); it!= simHits.end(); ++it)
//    	 {
// 	   std::cout << "i is " << ++i;
//            std::cout << "vahid"<<endl; 
// 	   	 }
// } 

void SimHitAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
  int N = NofEvents->GetBinContent(1);
  cout << " number of events = " << N; 
  for(long unsigned int i=0; i < tagDirectoriesDigi.size(); i++){
 DigiHitsPerSensor[i]->Divide(RhoArea[i],RhoAreaGeom);
 DigiHitsPerArea[i]->Divide(DigiHitsPerSensor[i],AreaGeom);
 DigiHitsPerSensor[i]->Scale(1./N);

 }

  
}

void
SimHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  //  edm::ParameterSetDescription desc;
  //desc.setUnknown();
  //descriptions.addDefault(desc);
}

    // edm::ESHandle<FbcmGeometry> theFbcmGeom;
// using cms::SiPadDigitizer;
//  DEFINE_DIGI_ACCUMULATOR(SiPadDigitizer);
DEFINE_FWK_MODULE(SimHitAnalyzer);
