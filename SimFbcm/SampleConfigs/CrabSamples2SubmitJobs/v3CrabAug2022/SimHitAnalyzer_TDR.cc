
// system include files
#include <typeinfo>
#include <vector>
#include<array>
#include <memory>
#include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
//#include "DataFormats/...

#include "DataFormats/FbcmDigi/interface/SiPadDigiData.h"
#include "DataFormats/FbcmDetId/interface/FbcmDetId.h"

#include "Geometry/Records/interface/FbcmGeometryRecord.h"
#include "Geometry/FbcmGeometry/interface/FbcmGeometry.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include <string> 
#include <TH2.h>
#include <TLatex.h>
#define MAX_SIMHIT_SIZE 100
#define MAX_DIGHIT_SIZE 100
#define MAX_BXSlot_SIZE 100
#define AREA_GROUP_LEN 8

#define TEST_IND 7

using namespace std;

class SimHitAnalyzer_TDR : public edm::EDAnalyzer {
public:
  explicit SimHitAnalyzer_TDR(const edm::ParameterSet&);
  ~SimHitAnalyzer_TDR();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;


  // ----------member data ---------------------------
//  edm::InputTag tag(hitsProducer_, SubdetName_);
  edm::EDGetTokenT< edm::DetSetVector<SiPadDigiData> > TokenTag_;
  edm::ESHandle<FbcmGeometry> theFbcmGeom;
  edm::ESHandle<FbcmGeometryRecord> theFbcmGeomTDR;

  //const std::vector<edm::ParameterSet> sectors_;
  
  
  TH1* nSensorsPerRho;  
 TH1* nSensorsPerRho2;  
  TH1* hNEvents;
  TH1* hNSensorGroups;
  
  TTree* theTree;
  TTree* FixedValuesTree;
    
  float SensorArea;
  float SensorRho;
  
  int SensorGroupIndex; // collected for each hit. 
  int numberOfSensorGroups; // maximum number of SiGroupIndex used in the digi step. 
  int nbrOfDiesPerRing; 

  int SensorZ_end;
  //int NumOfSensors;
  
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
  
  
  int nTotalRhuDigi;
  int nInterstedRhuBins; 
  int DigiRHUs[MAX_DIGHIT_SIZE] ;
  
  std::vector < int > tmpVect;
  std::vector < int > rhuInterestedHitBins_;

  TH2* DigiHitsPerArea;
  TH2*  AreaGeomTDR;
 TH2*  RhoArea;
 TH2*  RhoAreaGeom;
 TH2*   RhoAreaNormalized;
TH2*   RhoAreaNormalizedNew;

  //std::vector < int > rhuFullBxBinRange_;
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SimHitAnalyzer_TDR::SimHitAnalyzer_TDR(const edm::ParameterSet& iConfig) :
 
  numberOfSensorGroups(0)
{
  edm::InputTag tag("simFbcmDigis","SiPad");
  TokenTag_ = (consumes< edm::DetSetVector<SiPadDigiData> >(tag));
  numberOfSensorGroups=0;
  edm::Service<TFileService> fs;
  TFileDirectory subDir = fs->mkdir( "mySubDirectory" );
  
  
  hNEvents = subDir.make<TH1I>("hNEvents" , "" , 1 , 0 , 1 );
 nSensorsPerRho = subDir.make<TH1I>("nSensorsPerRho", "nSensorsPerRho", 70, 7.5, 22);
 // nSensorsPerRho2 = subDir.make<TH2I>("nSensorsPerRho", "nSensorsPerRho", 92, 8.14624, 21.9207, 1500, 0, 1500);


 nSensorsPerRho2 = subDir.make<TH2I>("nSensorsPerRho", "nSensorsPerRho", 92, 8.14624, 21.9207, 1500, 0, 1500);




  hNSensorGroups = subDir.make<TH1I>("hNSensorGroups", "", 1, 0, 1);
  //RhoArea = subDir.make<TH2F>("TDRRhoArea","" , 92 , 8.14624 ,21.9207 , 60 , 0.02, 0.09);
  //  RhoAreaNormalized = subDir.make<TH2F>("normalized","" , 92 , 8.14624 ,21.9207 , 100 , 0, 0.09);

 
  // RhoAreaGeom = subDir.make<TH2F>("TDRRhoAreaGeom","" , 92 , 8.14624 ,21.9207 , 60 , 0.02, 0.09);
  Int_t bsize = 32000; //default

  //theTree->Branch("SensorRho" , &SensorRho , "SensorRho/f[8,22]", bsize );
  //theTree->Branch("SensorArea" , &SensorArea , "SensorArea/f[0.01,1.0]" , bsize  );
  

  
  tmpVect = iConfig.getParameter< std::vector<int> >("RHU_InterestedHitBins");
  rhuInterestedHitBins_.clear();
  for (int i = tmpVect.front(); i <= tmpVect.back(); i++) 
      rhuInterestedHitBins_.emplace_back(i);
  
  //tmpVect = iConfig.getUntrackedParameter< std::vector<int> >("RHU_FullBxBinRange", {-2,1});
  //rhuFullBxBinRange_.clear();
  //for (int i = tmpVect.front(); i <= tmpVect.back(); i++) 
      //rhuFullBxBinRange_.emplace_back(i);

}


SimHitAnalyzer_TDR::~SimHitAnalyzer_TDR()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SimHitAnalyzer_TDR::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  using namespace std;

  hNEvents->Fill( 0.5 );

  edm::Handle< edm::DetSetVector<SiPadDigiData> > handle;
  iEvent.getByToken(TokenTag_,handle);

  
  //int SiDieId;
  for (edm::DetSetVector<SiPadDigiData>::const_iterator itDetSet = handle->begin() ;  itDetSet < handle->end() ; ++itDetSet) {
    // each Det in fbcmDets: multiple Data, but we have only one data
    if( itDetSet->size() != 1 )
      cout << "size == " << itDetSet->size() << endl; 
    for (std::vector<SiPadDigiData>::const_iterator Digidata_it = itDetSet->begin() ; Digidata_it < itDetSet->end() ; ++Digidata_it) {
      //const SiPadDigiData Digidata= *(Digidata_it);
      //SiDieId=Digidata_it->SiliconDieIndex() ;
      //SensorGroupIndex = SiDieId % nbrOfDiesPerRing;
      SensorGroupIndex = Digidata_it->SiGroupIndex();
      
      if ( SensorGroupIndex > numberOfSensorGroups )
          numberOfSensorGroups=SensorGroupIndex;
      
      SensorZ_end = Digidata_it->SideIndex();
      
      SensorRho=Digidata_it->Radius() ;
      SensorArea=Digidata_it->Area() ;
      //cout << "Sensor Area in TDR = " << SensorArea << endl;
      #RhoArea->Fill(SensorRho, SensorArea);
      nSimParticles = 0 ; 
      for(auto sim : Digidata_it->CahrgePsimVector()) {
	SimPdgIds[nSimParticles] = sim.second.ParticleType() ;
	SimPts[nSimParticles] = sim.second.Pabs() ;
	SimTofs[nSimParticles] = sim.second.Tof() ;
	SimCharges[nSimParticles] = sim.first;
	SimTof_perBx[nSimParticles] = sim.second.Tof() - 25.0 * sim.second.BunchCrossing() ;
	//SimHitTime = sim.second.time() - 25.0 * sim.second.BunchCrossing() ;
	nSimParticles++;
      }
 
      nValidDigiToAs=0;
      nValidDigiToTs=0;
	  nTotalRhuDigi=0;
      nInterstedRhuBins=0;
      
      BxSlotCnt=0; // just to make sure starting from non-zero index if in the case of using BxSlotNo < 0. 
      for(auto hit : Digidata_it->BxSlotHitAnalysisVector() ) {	  
            DigiHitStatus[BxSlotCnt] = hit.Bx_HitStatusInt(); 
            nDigiHits[BxSlotCnt] = hit.nbrOfRecognizedHitsInBx(); 
            BxSlotCnt++; 
                    
            for (auto ToaTot : hit.TotToaVectort()) {
                
              if ( (ToaTot.ToAToAStatusInt() == ToAStatus::FullyWithinBx || ToaTot.ToAToAStatusInt() == ToAStatus::WithinBx_LastsAfter ) ) 
              {
                  {
                      DigiToAs[nValidDigiToAs] = ToaTot.ToA();
                      PeakAmpl[nValidDigiToAs] = ToaTot.PeakAmplitude();
                      nValidDigiToAs++;
		      RhoArea->Fill(SensorRho,SensorArea);
                  }
                if (ToaTot.IsToTValid())
                  DigiToTs[nValidDigiToTs++] = ToaTot.ToT();
              
                DigiRHUs[nTotalRhuDigi++] = ToaTot.SubBxBinNumber();
                if (std::find(rhuInterestedHitBins_.begin(), rhuInterestedHitBins_.end(),ToaTot.SubBxBinNumber())!=rhuInterestedHitBins_.end())
                    nInterstedRhuBins++;
              
              }

            }
	
      }
  cout<< "BxSlotCnt = " << BxSlotCnt << endl;	
      //
      //      theTree->Fill();

    }

  }
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
SimHitAnalyzer_TDR::beginJob()
{
  
}

// ------------ method called once each job just after ending the event loop  ------------
void SimHitAnalyzer_TDR::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------

void SimHitAnalyzer_TDR::beginRun(edm::Run const&, edm::EventSetup const& iSetup) {
  
  edm::Service<TFileService> fs;
  TFileDirectory subDir = fs->mkdir( "mySubDirectory" );
 

  using namespace std;
  //std::cout << rhuInterestedHitBins_.front() << ", " << rhuInterestedHitBins_.back() << "\n" ; 
  //std::cout << rhuFullBxBinRange_.front() << ", " << rhuFullBxBinRange_.back() << "\n" ; 
  
    //for (auto v : rhuInterestedHitBins_) 
      //std::cout << v << ", ";
    //std::cout << "\n"; 
    //for (auto v : rhuFullBxBinRange_) 
      //std::cout << v << ", ";
 vector <Double_t> binedges;

  int nBinsX;
  int nBinsY;
  int N;
  float Rho;
  int SiDieId=0;
  int SiDieGroupIndex=0;
  float padX, padY , padRho;
  float minRho = 100, maxRho = 0, minArea = 10, maxArea= 0;
  float area;
  
  // const Int_t n;

  // ----------------------------------making necessary Histograms--------------------------------//

 
  iSetup.get<FbcmGeometryRecord>().get(theFbcmGeom); 
  const std::vector<const FbcmStationGeom*> AllStatitons = theFbcmGeom->Stations();
  const std::vector<FbcmSiPadGeom const*> allSiPadGeoms = theFbcmGeom->SiPads();
  cout << "type of AllStations is " << typeid(AllStatitons).name() << endl;
 
  unsigned int AllStatSize = AllStatitons.size();
  for (unsigned int i = 0 ; i <= AllStatSize ; i++)
    cout << " Number of Dies Per Ring = " <<  AllStatitons[i]->NumOfDiesPerRing() << " and ***** Number of Rings = " << AllStatitons[i]->NumOfRings() << endl;
  // for (const auto& station : AllStatitons) 
  //for(const auto& i :  station)
      //for (const auto& j :  
  //  cout << "Number of Dies per ring = " << i->NumOfDiesPerRing() << endl;
  
  // cin >> n ;
  binedges.clear();
  nBinsX = 0;
float  binwidthRho; 
  N = 0;
  Rho = 0;
  SiDieId = 0;
  SiDieGroupIndex = 0;
  //padX = 0;
  //padY = 0;
  //padRho = 0;

  minRho = 100; maxRho = 0; minArea = 10; maxArea= 0 ;

  //binedges.push_back(0);
  for(const auto& siPad : allSiPadGeoms)
    {
      std::pair<float, float> SiPadDimension = siPad->SiPadTopology().pitch();

      padX = SiPadDimension.first;
      padY = SiPadDimension.second;
      area = padX*padY;
      bool exists = find(binedges.begin(), binedges.end(), area) != binedges.end();
      if ( !exists ){
	cout << " area = " << area << endl;
	binedges.push_back(area);
	cout << "biedges content  = " << binedges[binedges.size()-1] << endl;
        cout<< "binages size = "  << binedges.size() << endl ;
	//	cin >> n;
      }
      // we could also store only the Area instead of padX and padY
      padRho = siPad->surface().position().perp(); 
	  
      SiDieId=siPad->id().SiliconDie();
      SiDieGroupIndex = SiDieId % nbrOfDiesPerRing;
	  
        //std::cout << SiDieId <<", " << SiDieGroupIndex <<"\n"; 
        cout << "SiDieGroupIndex = " << SiDieGroupIndex << endl;
      cout << "padx is " << padX << endl;	  
      cout << "padY is " << padY << endl;	  
      cout << "padRho is " << padRho << endl;
      minRho = min(minRho, padRho);
      maxRho = max(maxRho,padRho);
      
      minArea = min(minArea, padX*padY);
      maxArea = max(maxArea, padX*padY);
      binwidthRho = (maxRho -minRho)/sqrt(minArea);
      cout << " binwidth = "  << binwidthRho << endl;
      cout << " minRho = " << minRho << endl << "maxRho = " << maxRho << endl;
      cout << "minAre = " << minArea << endl << "maxArea = " << maxArea << endl;
      
      nSensorsPerRho->Fill(padRho);
      // FixedValuesTree->Fill();
    }

  cout<< "binages size = " << binedges.size() << endl ;
  const Int_t  n = binedges.size();
  binedges.push_back(0.09);
    int m  = binedges.size();
  const Int_t binNo = m-1;
  Double_t BinEdges[m]  ;
 for (int i = 0 ; i<m ; i++)
   BinEdges[i] =Double_t( binedges[i]) ;
 double* BinEdgesNew = &binedges[0]; 
 //cout << "size of BinEdges = " << BinEdges.size() << endl;
 // copy(binedges.begin(), binedges.end(), BinEdges); 
  for (double it :BinEdges )
    cout << " array  = " << it << endl;
  //int m;
  //cin >> m;
  const Int_t NBINS = 8;
  Double_t edges[NBINS + 1] = {BinEdges[0],BinEdges[1],BinEdges[2],BinEdges[3],BinEdges[4],BinEdges[5],BinEdges[6],BinEdges[7],BinEdges[7]+0.01};
  //  Double_t edges[NBINS + 1] = { 0.022, 0.028, 0.033, 0.038, 0.046, 0.053, 0.065, 0.080, 0.09};
  RhoAreaNormalizedNew = subDir.make<TH2F>("Newnormalized","Normalized " , Int_t(binwidthRho) , minRho ,maxRho , binNo,BinEdgesNew  ); 
 RhoAreaGeom = subDir.make<TH2F>("TDRRhoAreaGeom","Number of Sensors" , Int_t(binwidthRho) , minRho ,maxRho , binNo, BinEdgesNew);
 AreaGeomTDR = subDir.make<TH2F>("TDRAreaOfSensors","Area of Sensor in Rho" , Int_t(binwidthRho) , minRho ,maxRho , binNo, BinEdgesNew);
  RhoArea = subDir.make<TH2F>("DigisTDR","Number of digitized hits in TDR" , Int_t(binwidthRho) , minRho ,maxRho , binNo, BinEdgesNew);
DigiHitsPerArea = subDir.make<TH2F>("DigiHitPerAreaTDR","Number of digitized hits Per mm^{2} in TDR" , Int_t(binwidthRho) , minRho ,maxRho , binNo, BinEdgesNew);
 RhoAreaNormalized = subDir.make<TH2F>("normalized"," Digi Hits Per Sensor Per Event" , Int_t(binwidthRho) , minRho ,maxRho , binNo, BinEdgesNew );
  nBinsX = nSensorsPerRho->GetNbinsX();
  //  int nBinsY = nSensorsPerRho->GetNBinsY();

  for(const auto& siPad : allSiPadGeoms)
    {
      std::pair<float, float> SiPadDimension = siPad->SiPadTopology().pitch();

      padX = SiPadDimension.first;
      padY = SiPadDimension.second;
      area = padX*padY;
      // we could also store only the Area instead of padX and padY
      padRho = siPad->surface().position().perp(); 
	  
      SiDieId=siPad->id().SiliconDie();
      SiDieGroupIndex = SiDieId % nbrOfDiesPerRing;
      RhoAreaGeom->Fill(padRho,padX*padY); 
      RhoAreaNormalizedNew->Fill(padRho,padX*padY);
    }

 nBinsX = RhoAreaGeom->GetNbinsX();
 nBinsY = RhoAreaGeom->GetNbinsY();
       for (int binx  = 0 ; binx <=nBinsX ; binx ++)
          for ( int biny = 0 ; biny <=nBinsY; biny ++){
	    cout << "edges biny = " << edges[biny];
	    AreaGeomTDR->SetBinContent(binx,biny,edges[biny]);
		   }

  

  for (int binx = 0 ; binx <=nBinsX ; binx++)
   {
     N = nSensorsPerRho->GetBinContent(binx);
     Rho = nSensorsPerRho->GetBinCenter(binx);
     cout << "N = " << N << endl << "Rho = " << Rho << endl;
     nSensorsPerRho2->Fill(Rho,N);

   }
 }



// ------------ method called when ending the processing of a run  ------------

  void 
  SimHitAnalyzer_TDR::endRun(edm::Run const&, edm::EventSetup const&)
  {
    int myint, nBinsX, nBinsY, N, NEvents  ;
    float Rho;
    cout << "we are in endRun " << endl;
    NEvents=hNEvents->GetBinContent(1);
    RhoAreaNormalized->Divide(RhoArea,RhoAreaGeom);
    RhoAreaNormalized->Scale(1./NEvents);
    DigiHitsPerArea->Divide(RhoAreaNormalized, AreaGeomTDR);
    DigiHitsPerArea->Scale(1./100);
    RhoAreaGeom->GetXaxis()->SetTitle("Radius [cm]");
    RhoAreaGeom->GetYaxis()->SetTitle("Area[cm^{2}]");
    RhoArea->GetXaxis()->SetTitle("Radius[cm]");
    RhoArea->GetYaxis()->SetTitle("Area[cm^{2}]");
    RhoAreaNormalized->GetXaxis()->SetTitle("Radius [cm]");
 RhoAreaNormalized->GetYaxis()->SetTitle("Sensor Area [cm^{2}]");

nBinsX = nSensorsPerRho->GetNbinsX();
  //  int nBinsY = nSensorsPerRho->GetNBinsY();
  for (int binx = 0 ; binx <=nBinsX ; binx++)
   {
     N = nSensorsPerRho->GetBinContent(binx);
     Rho = nSensorsPerRho->GetBinCenter(binx);
     cout << "N = " << N << endl << "Rho = " << Rho << endl;
     //nSensorsPerRho2->Fill(Rho,N);
   }
 
  cout << "hNEvents = " << NEvents;
  cin >> myint ;
      numberOfSensorGroups++; // becuse indecies beginfrom zero
      
      //hNSensorGroups->Fill(0.5, nbrOfDiesPerRing); 
      hNSensorGroups->Fill(0.5, numberOfSensorGroups); 
     
  }


// ------------ method called when starting to processes a luminosity block  ------------
/*
  void 
  SimHitAnalyzer_TDR::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void 
  SimHitAnalyzer_TDR::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SimHitAnalyzer_TDR::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  //  edm::ParameterSetDescription desc;
  //desc.setUnknown();
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimHitAnalyzer_TDR);
