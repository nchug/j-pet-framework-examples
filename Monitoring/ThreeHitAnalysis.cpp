/**
 *  @copyright Copyright 2019 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *  @file ThreeHitAnalysis.cpp
 */

#include <JPetAnalysisTools/JPetAnalysisTools.h>
#include <JPetWriter/JPetWriter.h>
#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetStatistics/JPetStatistics.h>
#include "ThreeHitAnalysis.h"
#include "../LargeBarrelAnalysis/EventCategorizerTools.h"
#include "../LargeBarrelAnalysis/HitFinderTools.h"
#include <iostream>
#include <TMath.h>
#include <fstream>
#include <string>
#include <TVector3.h>
#include <TRandom3.h>

#include "reconstructor.h"
#include "./JPetHit/JPetHit.h"
#include <JPetMCHit/JPetMCHit.h>


using namespace jpet_options_tools;
using namespace std;

ThreeHitAnalysis::ThreeHitAnalysis(const char *name) : JPetUserTask(name) {}

bool ThreeHitAnalysis::init() {
  INFO("Event analysis started.");

 
 if (isOptionSet(fParams.getOptions(), kMC)) {
    fIsMC = getOptionAsBool(fParams.getOptions(), kMC);
 } 
     
  fOutputEvents = new JPetTimeWindow("JPetEvent");

  fReconstructor = Reconstructor::GetInstance();
  double scin_length = getParamBank().getScintillators().begin()->second->getScinSize(JPetScin::kLength);
  scin_length *= 0.1; // in cm
  fReconstructor->setBarrelLength( scin_length );
  fReconstructor->setChamberRadius( 10.0 );

  std::map<EvtType, std::string> EventTypeNames;
  EvtType event_type;

  initialiseHistograms();
  
  return true;
}


bool ThreeHitAnalysis::exec() {
  
  JPetTimeWindowMC *time_window_mc = nullptr;
  
  if (time_window_mc = dynamic_cast<JPetTimeWindowMC *const>(fEvent)) {
    fIsMC = true;
    EventTypeNames = EventTypeNames_MC;
  }
  else {
    fIsMC = false;
    EventTypeNames = EventTypeNames_data;
  }


   if (auto timeWindow = dynamic_cast<const JPetTimeWindow *const>(fEvent)) {
     
    for (uint i = 0; i < timeWindow->getNumberOfEvents(); i++) {
      
      const auto &event = dynamic_cast<const JPetEvent &>(timeWindow->operator[](i));

      int barrel_size = getParamBank().getBarrelSlotsSize();
     

      std::vector<JPetHit> hits = event.getHits();

      if( fIsMC == false){

      if (event.getHits().size() > 0 ) {
	 
	for (int j = 0; j <hits.size(); j++){

	  JPetHit hit = hits.at(j);
	  double tot_allHits = HitFinderTools::calculateTOT(hit, HitFinderTools::TOTCalculationType::kThresholdTrapeze)/1000;
	   
	  getStatistics().fillHistogram("tot_allHits", tot_allHits);	  
	  
	 }
       }
      }
      
      // LifeTime -- 
      
      if (event.getHits().size() == 1 ) {
	for (int j = 0; j <hits.size(); j++){
	  JPetHit hit = hits.at(j);

	  if (fIsMC == false){

	    double tot = HitFinderTools::calculateTOT(hit, HitFinderTools::TOTCalculationType::kThresholdTrapeze)/1000;
	    getStatistics().fillHistogram("tot_1HitEvt", tot);
	  
	    if (tot > 65){
	      pmtEmmTime = hit.getTime()/1000 - hit.getPos().Mag()/29.9792458;
	   
	    }
	  }

	  else if (fIsMC == true){
	    std::vector<double> eng = {hit.getEnergy(), 0, 0};
	    double tot = TOT(eng).at(0);
	    getStatistics().fillHistogram("tot_1HitEvt", tot);
	  
	    if (tot > 20){
	      pmtEmmTime = hit.getTime()/1000 - hit.getPos().Mag()/29.9792458;
	   
	    }
	  }
	    	
	}
      }
      
      
      EvtType event_type = UNKNOWN;
      
      getStatistics().fillHistogram( ("HitMult_"+ EventTypeNames.at(UNKNOWN)).c_str(), event.getHits().size());

      if(event.getHits().size() == 3 ) {            

	event_type = UNKNOWN;
      
        getStatistics().fillHistogram( "HitMult",  event_type);
	 
	vector<JPetHit> hits = event.getHits();
	  
	ReconHitsCalculation( &event, event_type = UNKNOWN);
	
	if ( fIsMC == true ){ //continue;
	
	int vtxIndex1, vtxIndex2, vtxIndex3;
	int VtxIndex1, VtxIndex2, VtxIndex3;
	Int_t comb[6][3] = { {0,1,2}, {0,2,1}, {1,2,0}};  // {1,0,2}, {2,0,1}, {2,1,0} };
	Int_t comb1[6][3] = { {0,1,2}, {0,2,1}, {1,2,0}, {1,0,2}, {2,0,1}, {2,1,0} };

       	genOPSCalculation( &event, *time_window_mc, event_type );
	pPs3G(&event, *time_window_mc, event_type);

       	if( event_type != OPS && event_type != PPS ){ 
	  
	  for ( int i = 0; i < 3; i++){
	  	     
	    JPetMCHit mchit1 = time_window_mc->getMCHit<JPetMCHit>( event.getHits().at(comb[i][0]).getMCindex() );
	    JPetMCHit mchit2 = time_window_mc->getMCHit<JPetMCHit>( event.getHits().at(comb[i][1]).getMCindex() );
	    JPetMCHit mchit3 = time_window_mc->getMCHit<JPetMCHit>( event.getHits().at(comb[i][2]).getMCindex() );

	    vtxIndex1 = mchit1.getMCVtxIndex();   
	    vtxIndex2 = mchit2.getMCVtxIndex();
	    vtxIndex3 = mchit3.getMCVtxIndex();
	    
	    
	    if ( mchit1.getGenGammaMultiplicity() == 3 && mchit2.getGenGammaMultiplicity() == 3 && mchit3.getGenGammaMultiplicity() == 103){

	       if (vtxIndex1 == vtxIndex2) {

		 if( event_type != UNKNOWN){
		   WARNING("Event type is assigned twice");
		 }
	        event_type = OPS_AND_1SCATTER;
		
	        
		getStatistics().fillHistogram( ("HitMult_"+ EventTypeNames.at(OPS_AND_1SCATTER)).c_str(), event.getHits().size() );
		getStatistics().fillHistogram( "HitMult",  event_type);
	         }
	    }
  
	    else if ( mchit1.getGenGammaMultiplicity() == 2 && mchit2.getGenGammaMultiplicity() == 2 && mchit3.getGenGammaMultiplicity() == 102){
	        if (vtxIndex1 == vtxIndex2) {

		  if( event_type != UNKNOWN){
		   WARNING("Event type is assigned twice");
		 }
	          
		  event_type = B2B_AND_SCATTER;
	    
		  getStatistics().fillHistogram( ("HitMult_"+ EventTypeNames.at(B2B_AND_SCATTER)).c_str(), event.getHits().size() );
		  getStatistics().fillHistogram( "HitMult", event_type );
	        }
	    }


	    else if ( mchit1.getGenGammaMultiplicity() == 2 && mchit2.getGenGammaMultiplicity() == 2 && mchit3.getGenGammaMultiplicity() == 1){
	        if (vtxIndex1 == vtxIndex2) {

		  if( event_type != UNKNOWN){
		   WARNING("Event type is assigned twice");
		 }
	          event_type = B2B_AND_PROMPT;
		 
		  getStatistics().fillHistogram( ("HitMult_"+ EventTypeNames.at(B2B_AND_PROMPT)).c_str(), event.getHits().size() );
		  getStatistics().fillHistogram( "HitMult",  event_type);

	        }
	    }

	    else if ( mchit1.getGenGammaMultiplicity() == 3 && mchit2.getGenGammaMultiplicity() == 3 && mchit3.getGenGammaMultiplicity() == 1){

	       if (vtxIndex1 == vtxIndex2) {

		 if( event_type != UNKNOWN){
		   WARNING("Event type is assigned twice");
		 }
	        
	        event_type = OPS_AND_PROMPT;
	     
		getStatistics().fillHistogram( ("HitMult_"+ EventTypeNames.at(OPS_AND_PROMPT)).c_str(), event.getHits().size() );
		getStatistics().fillHistogram( "HitMult",  event_type);
	      
	         }
	       }
	  }

	  // from here 6 different combination of mcHits start --- IMP
	  for ( int i = 0; i < 6; i++){

	    JPetMCHit mcHit1 = time_window_mc->getMCHit<JPetMCHit>( event.getHits().at(comb1[i][0]).getMCindex() );
	    JPetMCHit mcHit2 = time_window_mc->getMCHit<JPetMCHit>( event.getHits().at(comb1[i][1]).getMCindex() );
	    JPetMCHit mcHit3 = time_window_mc->getMCHit<JPetMCHit>( event.getHits().at(comb1[i][2]).getMCindex() );

	    VtxIndex1 = mcHit1.getMCVtxIndex();   
	    VtxIndex2 = mcHit2.getMCVtxIndex();
	    VtxIndex3 = mcHit3.getMCVtxIndex();

	    if ( (mcHit1.getGenGammaMultiplicity() == 3 && mcHit2.getGenGammaMultiplicity() == 103 && mcHit3.getGenGammaMultiplicity() == 203) )
	      {

	       if( event_type != UNKNOWN){
		   WARNING("Event type is assigned twice");
		 }
	       event_type = OPS_AND_2SCATTER;
	       
	       getStatistics().fillHistogram( ("HitMult_"+ EventTypeNames.at(OPS_AND_2SCATTER)).c_str() , event.getHits().size() );
	       getStatistics().fillHistogram( "HitMult",  event_type);
	    }
	    
	    else if ( mcHit1.getGenGammaMultiplicity() == 3 && mcHit2.getGenGammaMultiplicity() == 1 && mcHit3.getGenGammaMultiplicity() == 101){

		 if( event_type != UNKNOWN){
		   WARNING("Event type is assigned twice");
		 }
	        
	        event_type = OPS_AND_SCATTER_PROMPT;
	
		getStatistics().fillHistogram( ("HitMult_"+ EventTypeNames.at(OPS_AND_SCATTER_PROMPT)).c_str(), event.getHits().size() );
		getStatistics().fillHistogram( "HitMult",  event_type);
	      
	       }
	    
	    else if ( (mcHit1.getGenGammaMultiplicity() == 3 && mcHit2.getGenGammaMultiplicity() == 103 && mcHit3.getGenGammaMultiplicity() == 1) )
	      {

	       if( event_type != UNKNOWN){
		   WARNING("Event type is assigned twice");
		 }
	       event_type = OPS_SCATTER_AND_PROMPT ;
	       
	       getStatistics().fillHistogram( ("HitMult_"+ EventTypeNames.at(OPS_SCATTER_AND_PROMPT)).c_str() , event.getHits().size() );
	       getStatistics().fillHistogram( "HitMult",  event_type);
	    }
	    

	    else if ( (mcHit1.getGenGammaMultiplicity() == 2 && mcHit2.getGenGammaMultiplicity() == 102 && mcHit3.getGenGammaMultiplicity() == 202) )
	      {

	       if( event_type != UNKNOWN){
		   WARNING("Event type is assigned twice");
		 }
	       event_type = B2B_AND_2SCATTER;
	       
	       getStatistics().fillHistogram( ("HitMult_"+ EventTypeNames.at(B2B_AND_2SCATTER)).c_str() , event.getHits().size() );
	       getStatistics().fillHistogram( "HitMult",  event_type);
	    }

	    else if ( (mcHit1.getGenGammaMultiplicity() == 1 && mcHit2.getGenGammaMultiplicity() == 101 && mcHit3.getGenGammaMultiplicity() == 201) )
	      {

	       if( event_type != UNKNOWN){
		   WARNING("Event type is assigned twice");
		 }
	       event_type = PROMPT_SCATTER;
	       
	       
	       getStatistics().fillHistogram( ("HitMult_"+ EventTypeNames.at(PROMPT_SCATTER)).c_str() , event.getHits().size() );
	       getStatistics().fillHistogram( "HitMult",  event_type);
	    }

	    else if ( (mcHit1.getGenGammaMultiplicity() == 2 && mcHit2.getGenGammaMultiplicity() == 102 && mcHit3.getGenGammaMultiplicity() == 1) )
	      {

	       if( event_type != UNKNOWN){
		   WARNING("Event type is assigned twice");
		 }
	       event_type = B2B_SCATTER_AND_PROMPT;
	       
	       
	       getStatistics().fillHistogram( ("HitMult_"+ EventTypeNames.at(B2B_SCATTER_AND_PROMPT)).c_str() , event.getHits().size() );
	       getStatistics().fillHistogram( "HitMult",  event_type);
	    }

	     else if ( (mcHit1.getGenGammaMultiplicity() == 2 && mcHit2.getGenGammaMultiplicity() == 1 && mcHit3.getGenGammaMultiplicity() == 101) )
	      {

	       if( event_type != UNKNOWN){
		   WARNING("Event type is assigned twice");
		 }
	       event_type = B2B_PROMPT_AND_SCATTER;
	       
	       getStatistics().fillHistogram( ("HitMult_"+ EventTypeNames.at(B2B_PROMPT_AND_SCATTER)).c_str() , event.getHits().size() );
	       getStatistics().fillHistogram( "HitMult",  event_type);
	    }
	    
	     else if ( (mcHit1.getGenGammaMultiplicity() == 1 && mcHit2.getGenGammaMultiplicity() == 2 && mcHit3.getGenGammaMultiplicity() == 202) )
	      {

	       if( event_type != UNKNOWN){
		   WARNING("Event type is assigned twice");
		 }
	       event_type = PRMT_B2B_2SCAT;
	       
	       
	       getStatistics().fillHistogram( ("HitMult_"+ EventTypeNames.at(PRMT_B2B_2SCAT)).c_str() , event.getHits().size() );
	       getStatistics().fillHistogram( "HitMult",  event_type);
	      }

	    else if ( (mcHit1.getGenGammaMultiplicity() == 1 && mcHit2.getGenGammaMultiplicity() == 2 && mcHit3.getGenGammaMultiplicity() == 201) )
	      {

	       if( event_type != UNKNOWN){
		   WARNING("Event type is assigned twice");
		 }
	       event_type = B2B_PRMT_2SCAT;
	       
	       
	       getStatistics().fillHistogram( ("HitMult_"+ EventTypeNames.at(B2B_PRMT_2SCAT)).c_str() , event.getHits().size() );
	       getStatistics().fillHistogram( "HitMult",  event_type);
	      }


	  } 

	if ( event_type == UNKNOWN) {

	  JPetMCHit MCHit1 = time_window_mc->getMCHit<JPetMCHit>( event.getHits().at(0).getMCindex() );
          JPetMCHit MCHit2 = time_window_mc->getMCHit<JPetMCHit>( event.getHits().at(1).getMCindex() );
          JPetMCHit MCHit3 = time_window_mc->getMCHit<JPetMCHit>( event.getHits().at(2).getMCindex() );

	  vtxIndex1 = MCHit1.getMCVtxIndex();   
	  vtxIndex2 = MCHit2.getMCVtxIndex();
	  vtxIndex3 = MCHit3.getMCVtxIndex();
 
	  getStatistics().fillHistogram( "genGammaMult_g1",  MCHit1.getGenGammaMultiplicity(), (MCHit2.getGenGammaMultiplicity()+ MCHit3.getGenGammaMultiplicity()));
	  getStatistics().fillHistogram( "genGammaMult_g2",  MCHit2.getGenGammaMultiplicity(), (MCHit1.getGenGammaMultiplicity()+ MCHit3.getGenGammaMultiplicity()));
	  getStatistics().fillHistogram( "genGammaMult_g3",  MCHit3.getGenGammaMultiplicity(), (MCHit2.getGenGammaMultiplicity()+ MCHit1.getGenGammaMultiplicity()));

	  getStatistics().fillHistogram( "hit_genMult", 1  , MCHit1.getGenGammaMultiplicity() );
	  getStatistics().fillHistogram( "hit_genMult", 2  , MCHit2.getGenGammaMultiplicity() );
	  getStatistics().fillHistogram( "hit_genMult", 3  , MCHit3.getGenGammaMultiplicity() );

	  getStatistics().fillHistogram( "hit_vtx", 1  , vtxIndex1 );
	  getStatistics().fillHistogram( "hit_vtx", 2  , vtxIndex2 );
	  getStatistics().fillHistogram( "hit_vtx", 3  , vtxIndex3 );

	  
	   if (vtxIndex1 == vtxIndex2 && vtxIndex2 == vtxIndex3){
	     event_type = ALL_THREE_FROM_SAME_VTX;

	     getStatistics().fillHistogram( ("HitMult_"+ EventTypeNames.at(ALL_THREE_FROM_SAME_VTX)).c_str(), event.getHits().size() );
	     getStatistics().fillHistogram( "HitMult", event_type);
	     //getStatistics().fillHistogram( "HitMult_vtx",  1);

	     getStatistics().fillHistogram( "hit_genMult_sameVtx", 1  , MCHit1.getGenGammaMultiplicity() );
	     getStatistics().fillHistogram( "hit_genMult_sameVtx", 2  , MCHit2.getGenGammaMultiplicity() );
	     getStatistics().fillHistogram( "hit_genMult_sameVtx", 3  , MCHit3.getGenGammaMultiplicity() );

	     if ( scatterTest (&event) > 20.){
	     
	     if( MCHit1.getGenGammaMultiplicity() == 1 ){
 
	       if( MCHit2.getGenGammaMultiplicity() < MCHit3.getGenGammaMultiplicity()) {
		 
		 getStatistics().fillHistogram( "genMult1_sameVtx_g2g3", MCHit2.getGenGammaMultiplicity()  , MCHit3.getGenGammaMultiplicity() );
	       }
	       else {
		 getStatistics().fillHistogram( "genMult1_sameVtx_g2g3", MCHit3.getGenGammaMultiplicity()  , MCHit2.getGenGammaMultiplicity() );
	       }    
	     }

	     if( MCHit1.getGenGammaMultiplicity() == 2 ){
 
	       if( MCHit2.getGenGammaMultiplicity() < MCHit3.getGenGammaMultiplicity()) {
		 
		 getStatistics().fillHistogram( "genMult2_sameVtx_g2g3", MCHit2.getGenGammaMultiplicity()  , MCHit3.getGenGammaMultiplicity() );
	       }
	       else {
		 getStatistics().fillHistogram( "genMult2_sameVtx_g2g3", MCHit3.getGenGammaMultiplicity()  , MCHit2.getGenGammaMultiplicity() );
	       }    
	     }

	     if( MCHit1.getGenGammaMultiplicity() == 3 ){
 
	       if( MCHit2.getGenGammaMultiplicity() < MCHit3.getGenGammaMultiplicity()) {
		 
		 getStatistics().fillHistogram( "genMult3_sameVtx_g2g3", MCHit2.getGenGammaMultiplicity()  , MCHit3.getGenGammaMultiplicity() );
	       }
	       else {
		 getStatistics().fillHistogram( "genMult3_sameVtx_g2g3", MCHit3.getGenGammaMultiplicity()  , MCHit2.getGenGammaMultiplicity() );
	       }    
	     }
	     
	     if ( calculatedLOR( &event) > 5.){

	       if( MCHit1.getGenGammaMultiplicity() == 1 ){
 
	       if( MCHit2.getGenGammaMultiplicity() < MCHit3.getGenGammaMultiplicity()) {
		 
		 getStatistics().fillHistogram( "genMult1_sameVtx_g2g3_afterDLOR", MCHit2.getGenGammaMultiplicity()  , MCHit3.getGenGammaMultiplicity() );
	       }
	       else{
		 getStatistics().fillHistogram( "genMult1_sameVtx_g2g3_afterDLOR", MCHit3.getGenGammaMultiplicity()  , MCHit2.getGenGammaMultiplicity() );
	       }
	       }
	     }
	     }
	     
	  }

	  else if (vtxIndex1 == vtxIndex2 ||  vtxIndex2 == vtxIndex3 || vtxIndex1 == vtxIndex3 ){
	    event_type = TWO_FROM_SAME_VTX;
	    
	    getStatistics().fillHistogram( ("HitMult_"+ EventTypeNames.at(TWO_FROM_SAME_VTX)).c_str(), event.getHits().size() );
	    getStatistics().fillHistogram( "HitMult", event_type);
	    
	    // getStatistics().fillHistogram( "HitMult_vtx", 2);

	     getStatistics().fillHistogram( "hit_genMult_2gSameVtx", 1  , MCHit1.getGenGammaMultiplicity() );
	     getStatistics().fillHistogram( "hit_genMult_2gSameVtx", 2  , MCHit2.getGenGammaMultiplicity() );
	     getStatistics().fillHistogram( "hit_genMult_2gSameVtx", 3  , MCHit3.getGenGammaMultiplicity() );
	     
	  }

	  else if(vtxIndex1 != vtxIndex2 && vtxIndex2 != vtxIndex3 ){
	    event_type = NONE_FROM_SAME_VTX;

	    getStatistics().fillHistogram( ("HitMult_"+ EventTypeNames.at(NONE_FROM_SAME_VTX)).c_str(), event.getHits().size() );
	    getStatistics().fillHistogram( "HitMult", event_type);
	    
	  }
	   
	  else {

	    event_type = OTHER;
	    getStatistics().fillHistogram( ("HitMult_"+ EventTypeNames.at(OTHER)).c_str(), event.getHits().size() );
	    getStatistics().fillHistogram( "HitMult", event_type);
	    
	  }
	    }
	}
      
	ReconHitsCalculation( &event, event_type);
	getStatistics().fillHistogram(("lifeTime_"+ EventTypeNames.at(event_type)).c_str(), (decayTime_3hits(&event,event_type) - pmtEmmTime)/1000);

	} //end of if bool for data

	if (fIsMC == false){
	  
	  decayTime = decayTime_3hits(&event,event_type );	  
	  lifeTime = decayTime - pmtEmmTime;
	  getStatistics().fillHistogram(("lifeTime_"+ EventTypeNames.at(event_type)).c_str(), lifeTime/1000);
	  
	}
      }
      
    }
   }
   
 else {
    return false;
  }

  return true;
}

void ThreeHitAnalysis::pPs3G( const JPetEvent *event, const JPetTimeWindowMC& time_window_mc, EvtType& event_type ) {

  int vtxIndex1, vtxIndex2, vtxIndex3;
  
  JPetMCHit mchit1 = time_window_mc.getMCHit<JPetMCHit>( event->getHits().at(0).getMCindex() );
  JPetMCHit mchit2 = time_window_mc.getMCHit<JPetMCHit>( event->getHits().at(1).getMCindex() );
  JPetMCHit mchit3 = time_window_mc.getMCHit<JPetMCHit>( event->getHits().at(2).getMCindex() );

  if( mchit1.getGenGammaMultiplicity() == 2 && mchit2.getGenGammaMultiplicity() == 2 && mchit3.getGenGammaMultiplicity() == 2){
                                          	      
    vtxIndex1 = mchit1.getMCVtxIndex();
    vtxIndex2 = mchit2.getMCVtxIndex();
    vtxIndex3 = mchit3.getMCVtxIndex();
         
    if ( vtxIndex1 == vtxIndex2 && vtxIndex2 == vtxIndex3) {

      event_type = PPS;
      
      getStatistics().fillHistogram(("HitMult_"+ EventTypeNames.at(PPS)).c_str(), event->getHits().size());
      getStatistics().fillHistogram( "HitMult", event_type);
    }
  }
}


void ThreeHitAnalysis::genOPSCalculation( const JPetEvent *event, const JPetTimeWindowMC& time_window_mc, EvtType& event_type ) {
  	    
  JPetHit hit1 = event->getHits().at(0);
  JPetHit hit2 = event->getHits().at(1);
  JPetHit hit3 = event->getHits().at(2);

  JPetMCHit mchit1 = time_window_mc.getMCHit<JPetMCHit>( hit1.getMCindex() );
  JPetMCHit mchit2 = time_window_mc.getMCHit<JPetMCHit>( hit2.getMCindex() );
  JPetMCHit mchit3 = time_window_mc.getMCHit<JPetMCHit>( hit3.getMCindex() );

  //genHitsCalculation_3G( &event, mchit1, mchit2, mchit3, hit1, hit2, hit3);
  Int_t combs[][2] = { {0,1}, {1,0}, {2,1}, {1,2}, {2,0}, {0,2} };
  TRandom3 rndm;	
  double mass_e = 510.99;
  int vtxIndex1, vtxIndex2, vtxIndex3;
  
    //--------------------Generated gamma multiplicity condition --------------//

  if( mchit1.getGenGammaMultiplicity() == 3 && mchit2.getGenGammaMultiplicity() == 3 && mchit3.getGenGammaMultiplicity() == 3){
                                          	      
    vtxIndex1 = mchit1.getMCVtxIndex();
    vtxIndex2 = mchit2.getMCVtxIndex();
    vtxIndex3 = mchit3.getMCVtxIndex();
          
    //-----Hits from same vertex ------//
    if ( vtxIndex1 == vtxIndex2 && vtxIndex2 == vtxIndex3) {  

       
       event_type = OPS;
      
       getStatistics().fillHistogram(("HitMult_"+ EventTypeNames.at(OPS)).c_str(), event->getHits().size());
       getStatistics().fillHistogram( "HitMult", event_type);

       TVector3 sum_mom;
       TVector3 p1 = mchit1.getMomentum();
       TVector3 p2 = mchit2.getMomentum();
       TVector3 p3 = mchit3.getMomentum();

       sum_mom = p1 + p2 + p3;
	      
       TVector3 Pos1 = mchit1.getPos();
       TVector3 Pos2 = mchit2.getPos();
       TVector3 Pos3 = mchit3.getPos();
             
       // Reconstructed hits
       // ReconHitsCalculation( hit1, hit2, hit3, event_type);
       
       //---------------Annihilation Point Calculation -------------------------//
       
       double k1, k2;             // Constants
	     
       TVector3 vec1( Pos2.X() - Pos1.X(), Pos2.Y() - Pos1.Y(), Pos2.Z() - Pos1.Z());
	     
       double a1 = p1.X()/p1.Mag();
       double a2 = p2.X()/p2.Mag();

       double b1 = p1.Y()/p1.Mag();
       double b2 = p2.Y()/p2.Mag();
		
       k1 = (vec1.X() * b2 - vec1.Y() * a2)/(a1 * b2 - b1 * a2);
       k2 = (vec1.X() * b1 - vec1.Y() * a1)/(b2 *a1 - a2*b1);

       double constValue1 = k1/p1.Mag();
       double constValue2 = k2/p2.Mag();

       TVector3 annhPoint = Pos1 + k1*p1.Unit();
             
       //double constValue = k1/Pos1.Mag();

       // TVector3 annhPoint = Pos1 + constValue1*Pos1;
	    
       double annh_radius = annhPoint.Perp();
	     
       // ----------------------------Energy Deposited-------------------------//
       vector<double> EngDep;
       
       EngDep.push_back( mchit1.getEnergy());	
       EngDep.push_back( mchit2.getEnergy());
       EngDep.push_back( mchit3.getEnergy());
       sort( EngDep.begin(), EngDep.end());


       //-----------------Energies from mom magnitude
       double e1, e2, e3;
       e1 = p1.Mag();
       e2 = p2.Mag();
       e3 = p3.Mag();

       //--------------------- Relative Angles (in radian) -------------------------//
       std::vector<double> angles(3);
       std::vector<double> angles_ordered;
       angles[0] = p1.Angle(p2);//*TMath::RadToDeg();
       angles[1] = p2.Angle(p3);//*TMath::RadToDeg();
       angles[2] = p1.Angle(p3);//*TMath::RadToDeg();
            

       // Randomly 
       int wc = rndm.Integer(3);

       angles_ordered.push_back(angles[0]);
       angles_ordered.push_back(angles[1]);
       angles_ordered.push_back(angles[2]);
       sort( angles_ordered.begin(), angles_ordered.end());
	
       //-------------------Energies from angles
       std::vector<double> energy(3);
       energy[0] = -2*mass_e*(-cos(angles[2])+cos(angles[0])*cos(angles[1]))/((-1+cos(angles[0]))*(1+cos(angles[0])-cos(angles[2])-cos(angles[1])));
       energy[1] = -2*mass_e*(cos(angles[0])*cos(angles[2])-cos(angles[1]))/((-1+cos(angles[0]))*(1+cos(angles[0])-cos(angles[2])-cos(angles[1])));
       energy[2] = 2*mass_e*(1+cos(angles[0]))/(1+cos(angles[0])-cos(angles[2])-cos(angles[1]));

       if (energy[0] > 0 && energy[1] > 0 && energy[2] > 0){
	  
          // --------------------- Spin ---------------//
	  vector <pair <double, int>> momentum;
	  momentum.push_back({p1.Mag(),0});
	  momentum.push_back({p2.Mag(),1});
	  momentum.push_back({p3.Mag(),2});	  
	  sort(momentum.begin(), momentum.end());

	  // hits arranged in decreasing order of mangnitudes of momenta
	  JPetHit hit1_sorted = event->getHits().at(momentum.at(2).second);
          JPetHit hit2_sorted = event->getHits().at(momentum.at(1).second);
	  JPetHit hit3_sorted = event->getHits().at(momentum.at(0).second);

	  JPetMCHit mchit1_sorted = time_window_mc.getMCHit<JPetMCHit>( hit1_sorted.getMCindex() );
	  JPetMCHit mchit2_sorted = time_window_mc.getMCHit<JPetMCHit>( hit2_sorted.getMCindex() );
	  JPetMCHit mchit3_sorted = time_window_mc.getMCHit<JPetMCHit>( hit3_sorted.getMCindex() );   
  
	  TVector3 p1_sorted = mchit1_sorted.getMomentum();
	  TVector3 p2_sorted = mchit2_sorted.getMomentum();
	  TVector3 p3_sorted = mchit3_sorted.getMomentum();

	  // p1_gen = p1_sorted.Unit();
	  TVector3 planeNormal = p1_sorted.Cross(p2_sorted).Unit();
	  TVector3 s = annhPoint.Unit();
          double operators[2];
	  operators[0] = s.Dot(p1_sorted.Unit());
	  operators[1] = s.Dot(planeNormal);
	  
       }

    }
  }
  
}


void ThreeHitAnalysis::ReconHitsCalculation( const JPetEvent *event,  EvtType& event_type) {         
  // Reconstructed hits -- Annihilation Point ( Trilateration Method ) (Alek's code)

  JPetHit hit1 = event->getHits().at(0);
  JPetHit hit2 = event->getHits().at(1);
  JPetHit hit3 = event->getHits().at(2);
  
  fReconstructor->setGamma(0, event->getHits().at(0));
  fReconstructor->setGamma(1, event->getHits().at(1));
  fReconstructor->setGamma(2, event->getHits().at(2));
  
  int error = 0;
  TVector3 sol[2];
  double t[2];
  error = fReconstructor->getSolution(sol[1], t[1], 1);
  getStatistics().fillHistogram("error_tri", error);

  TVector3 decayPoint = sol[1];
  std::vector<TVector3> momenta(3);	
  std::vector<pair<double, int>> k;
  std::vector<double> theta(3);
  std::vector<double> energies(3);
  TVector3 spin, normaltoPlane;
  double OperatorStudy[2];
  
  momenta.at(0) = hit1.getPos() - decayPoint;
  momenta.at(1) = hit2.getPos() - decayPoint;
  momenta.at(2) = hit3.getPos() - decayPoint;
  double radius = decayPoint.Perp();

  double transRad = pow( (pow(decayPoint.X(),2) + pow(decayPoint.Y(),2)),0.5);

  //p1_recs = hit1.getPos().Unit();
  //p1_uncer = p1_recs.Angle(p1_gen);
 

  std::vector<double> time;
  double vel =  29.9792458;
  time.push_back(momenta[0].Mag()/vel);
  time.push_back(momenta[1].Mag()/vel);
  time.push_back(momenta[2].Mag()/vel);
  meanTime_3Hits = (time[0]+time[1]+time[2])/3;
  
  Int_t combs[][2] = { {0,1}, {1,0}, {2,1}, {1,2}, {2,0}, {0,2} };
  TRandom3 rndm;
  int wc = rndm.Integer(3);

  theta[0] = momenta[0].Angle(momenta[1]);
  theta[1] = momenta[1].Angle(momenta[2]);
  theta[2] = momenta[0].Angle(momenta[2]);

  std::vector<double> theta_rndm;
  theta_rndm.push_back(theta[0]);
  theta_rndm.push_back(theta[1]);
  theta_rndm.push_back(theta[2]);
  
  double mass_e = 510.99;
  energies[0] = -2*mass_e*(-cos(theta[2])+cos(theta[0])*cos(theta[1]))/((-1+cos(theta[0]))*(1+cos(theta[0])-cos(theta[2])-cos(theta[1])));
  energies[1] = -2*mass_e*(cos(theta[0])*cos(theta[2])-cos(theta[1]))/((-1+cos(theta[0]))*(1+cos(theta[0])-cos(theta[2])-cos(theta[1])));
  energies[2] = 2*mass_e*(1+cos(theta[0]))/(1+cos(theta[0])-cos(theta[2])-cos(theta[1]));
 
  std::vector<double> eng;       // deposited energy by gammas
  eng.push_back(hit1.getEnergy());
  eng.push_back(hit2.getEnergy());
  eng.push_back(hit3.getEnergy());
    
  k.push_back( {momenta[0].Mag(),0});
  k.push_back( {momenta[1].Mag(),1});
  k.push_back( {momenta[2].Mag(),2});
  std::sort( k.begin(), k.end());

  // using x and y  cord of hits
  std::vector<double> thetas;
  thetas.push_back(calculateAngle (hit1, hit2));
  thetas.push_back(calculateAngle (hit2, hit3));
  thetas.push_back(calculateAngle (hit3, hit1));
  std::sort (thetas.begin(), thetas.end());

  // using x y and z pos of hits
  std::vector<double> angle_3D;
  angle_3D.push_back(calculateAngle_3D (hit1, hit2));
  angle_3D.push_back(calculateAngle_3D (hit2, hit3));
  angle_3D.push_back(calculateAngle_3D (hit3, hit1));
  
  if ( error == 0){

    getStatistics().fillHistogram( "HitMult_beforeEng", event_type);

    getStatistics().fillHistogram(("engDep_g1_"+ EventTypeNames.at(event_type)).c_str(), hit1.getEnergy());   //eng.at(0)); //ns
    getStatistics().fillHistogram(("engDep_g2_"+ EventTypeNames.at(event_type)).c_str(), eng.at(1));
    getStatistics().fillHistogram(("engDep_g3_"+ EventTypeNames.at(event_type)).c_str(), eng.at(2));
    
    if (energies[0] > 0 && energies[1] > 0 && energies[2] > 0){
      getStatistics().fillHistogram( "HitMult_afterEng", event_type);

      getStatistics().fillHistogram(("totHit1_"+ EventTypeNames.at(event_type)).c_str(), TOT( eng ).at(0));
      getStatistics().fillHistogram(("totHit2_"+ EventTypeNames.at(event_type)).c_str(), TOT( eng ).at(1));
      getStatistics().fillHistogram(("totHit3_"+ EventTypeNames.at(event_type)).c_str(), TOT( eng ).at(2));

      getStatistics().fillHistogram(("scattTest_"+ EventTypeNames.at(event_type)).c_str(), scatterTest( event) );
     
      std::sort(theta.begin(), theta.end()); // from photons' momenta
      getStatistics().fillHistogram(("sum_diff_angles_"+ EventTypeNames.at(event_type)).c_str(), (theta.at(0)+ theta.at(1))*TMath::RadToDeg(), (theta.at(1)- theta.at(0))*TMath::RadToDeg());

      getStatistics().fillHistogram(("dLOR_"+ EventTypeNames.at(event_type)).c_str(), calculatedLOR( event));

      getStatistics().fillHistogram(("sumDiff_azmTheta_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), (thetas[1]- thetas[0]));

      getStatistics().fillHistogram(("angle_dLOR_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calculatedLOR( event));

      getStatistics().fillHistogram(("angle_dLOR_fromSph_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calcLOR_sphRadius(event) );

      getStatistics().fillHistogram(("angle_dLOR_fromCyd_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calcLOR_cydRadius(event) );

      getStatistics().fillHistogram(("angle_dLOR_fromTril_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calcLOR_recsAnnhPt( event, decayPoint ).at(0));
      
      getStatistics().fillHistogram(("dLOR_fromAnnhPt_"+ EventTypeNames.at(event_type)).c_str(), calcLOR_recsAnnhPt( event, decayPoint ).at(0), calcLOR_recsAnnhPt( event, decayPoint ).at(2));
      getStatistics().fillHistogram(("dLOR_fromAnnhPt_total_"+ EventTypeNames.at(event_type)).c_str(), calcLOR_recsAnnhPt( event, decayPoint ).at(0), (calcLOR_recsAnnhPt( event, decayPoint ).at(0)+calcLOR_recsAnnhPt( event, decayPoint ).at(1)+calcLOR_recsAnnhPt( event, decayPoint ).at(2)));
				 
      getStatistics().fillHistogram(("shortDlor_comp_"+ EventTypeNames.at(event_type)).c_str(), calcLOR_sphRadius(event), calcLOR_recsAnnhPt( event, decayPoint ).at(0));

      getStatistics().fillHistogram(("distance_hits_"+ EventTypeNames.at(event_type)).c_str(), distance_hits(event).at(0) );
      getStatistics().fillHistogram(("distance_hits_"+ EventTypeNames.at(event_type)).c_str(), distance_hits(event).at(1) );
      getStatistics().fillHistogram(("distance_hits_"+ EventTypeNames.at(event_type)).c_str(), distance_hits(event).at(2) );

      getStatistics().fillHistogram(("hitsDis_diff_"+ EventTypeNames.at(event_type)).c_str(), hitsDis_diff(event).at(0), hitsDis_diff(event).at(2) );
	   
      getStatistics().fillHistogram(("tot_data_"+ EventTypeNames.at(event_type)).c_str(), tot_data(event).at(0) );
      getStatistics().fillHistogram(("tot_data_"+ EventTypeNames.at(event_type)).c_str(), tot_data(event).at(1) );
      getStatistics().fillHistogram(("tot_data_"+ EventTypeNames.at(event_type)).c_str(), tot_data(event).at(2) );

      getStatistics().fillHistogram(("AnnhPointXYComp_"+ EventTypeNames.at(event_type)).c_str(), sol[1].X(), sol[1].Y());
      getStatistics().fillHistogram(("AnnhPointYZComp_"+ EventTypeNames.at(event_type)).c_str(), sol[1].Y(), sol[1].Z());

      getStatistics().fillHistogram(("tot_adj_"+ EventTypeNames.at(event_type)).c_str(), tot_adj(event).at(0) );
      getStatistics().fillHistogram(("tot_adj_"+ EventTypeNames.at(event_type)).c_str(), tot_adj(event).at(1) );
      getStatistics().fillHistogram(("tot_adj_"+ EventTypeNames.at(event_type)).c_str(), tot_adj(event).at(2) );
   
      
      if( (tot_data(event).at(0) < 40) && (tot_data(event).at(1) < 40) && (tot_data(event).at(2) < 40) ){

	   getStatistics().fillHistogram( "HitMult_afterTOT", event_type);

	   getStatistics().fillHistogram(("scattTest_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), scatterTest( event) );

	   getStatistics().fillHistogram(("sum_diff_angles_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), (theta.at(0)+ theta.at(1))*TMath::RadToDeg(), (theta.at(1)- theta.at(0))*TMath::RadToDeg());

	   getStatistics().fillHistogram(("dLOR_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), calculatedLOR( event));

	   getStatistics().fillHistogram(("sumDiff_azmTheta_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), (thetas[1]- thetas[0]) );

	   getStatistics().fillHistogram(("angle_dLOR_afterTOT_"+ EventTypeNames.at(event_type)).c_str(),(thetas[0]+ thetas[1]) , calculatedLOR( event));

	   getStatistics().fillHistogram(("angle_dLOR_fromSph_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calcLOR_sphRadius(event) );

	   getStatistics().fillHistogram(("angle_dLOR_fromCyd_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calcLOR_cydRadius(event) );

	   getStatistics().fillHistogram(("angle_dLOR_fromTril_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calcLOR_recsAnnhPt( event, decayPoint ).at(0) );

	   getStatistics().fillHistogram(("dLOR_fromAnnhPt_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), calcLOR_recsAnnhPt( event, decayPoint ).at(0), calcLOR_recsAnnhPt( event, decayPoint ).at(2));
	   getStatistics().fillHistogram(("dLOR_fromAnnhPt_total_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), calcLOR_recsAnnhPt( event, decayPoint ).at(0), (calcLOR_recsAnnhPt( event, decayPoint ).at(0)+calcLOR_recsAnnhPt( event, decayPoint ).at(1)+calcLOR_recsAnnhPt( event, decayPoint ).at(2)));

	   getStatistics().fillHistogram(("shortDlor_comp_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), calcLOR_sphRadius(event), calcLOR_recsAnnhPt( event, decayPoint ).at(0));

	   getStatistics().fillHistogram(("distance_hits_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), distance_hits(event).at(0) );
	   getStatistics().fillHistogram(("distance_hits_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), distance_hits(event).at(1) );
	   getStatistics().fillHistogram(("distance_hits_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), distance_hits(event).at(2) );

	   getStatistics().fillHistogram(("hitsDis_diff_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), hitsDis_diff(event).at(0), hitsDis_diff(event).at(2) );

	    getStatistics().fillHistogram(("tot_data_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), tot_data(event).at(0) );
	    getStatistics().fillHistogram(("tot_data_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), tot_data(event).at(1) );
	    getStatistics().fillHistogram(("tot_data_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), tot_data(event).at(2) );

	   getStatistics().fillHistogram(("AnnhPointXYComp_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), sol[1].X(), sol[1].Y());
	   getStatistics().fillHistogram(("AnnhPointYZComp_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), sol[1].Y(), sol[1].Z());
      
	   
	   /*   
	if ( (TOT( eng).at(0) < 15) &&  (TOT( eng).at(1) < 15 ) && (TOT( eng ).at(2) < 15) ) {
	  
	  // if ( fIsMC == true ){
      
	    getStatistics().fillHistogram( "HitMult_afterTOT", event_type);

	    getStatistics().fillHistogram(("totHit1_afterTOTcut_"+ EventTypeNames.at(event_type)).c_str(), TOT( eng ).at(0));
	    getStatistics().fillHistogram(("totHit2_afterTOTcut_"+ EventTypeNames.at(event_type)).c_str(), TOT( eng ).at(1));
	    getStatistics().fillHistogram(("totHit3_afterTOTcut_"+ EventTypeNames.at(event_type)).c_str(), TOT( eng ).at(2));

	    getStatistics().fillHistogram(("scattTest_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), scatterTest( event) );

	    getStatistics().fillHistogram(("sum_diff_angles_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), (theta.at(0)+ theta.at(1))*TMath::RadToDeg(), (theta.at(1)- theta.at(0))*TMath::RadToDeg());

	    getStatistics().fillHistogram(("dLOR_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), calculatedLOR( event));

	    getStatistics().fillHistogram(("sumDiff_azmTheta_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), (thetas[1]- thetas[0]) );

	    getStatistics().fillHistogram(("angle_dLOR_afterTOT_"+ EventTypeNames.at(event_type)).c_str(),(thetas[0]+ thetas[1]) , calculatedLOR( event));

	    getStatistics().fillHistogram(("angle_dLOR_fromSph_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calcLOR_sphRadius(event) );

	    getStatistics().fillHistogram(("angle_dLOR_fromCyd_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calcLOR_cydRadius(event) );

	    getStatistics().fillHistogram(("angle_dLOR_fromTril_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calcLOR_recsAnnhPt( event, decayPoint ).at(0) );

            getStatistics().fillHistogram(("dLOR_fromAnnhPt_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), calcLOR_recsAnnhPt( event, decayPoint ).at(0), calcLOR_recsAnnhPt( event, decayPoint ).at(2));
            getStatistics().fillHistogram(("dLOR_fromAnnhPt_total_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), calcLOR_recsAnnhPt( event, decayPoint ).at(0), (calcLOR_recsAnnhPt( event, decayPoint ).at(0)+calcLOR_recsAnnhPt( event, decayPoint ).at(1)+calcLOR_recsAnnhPt( event, decayPoint ).at(2)));

	    getStatistics().fillHistogram(("shortDlor_comp_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), calcLOR_sphRadius(event), calcLOR_recsAnnhPt( event, decayPoint ).at(0));

	    getStatistics().fillHistogram(("distance_hits_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), distance_hits(event).at(0) );
	    getStatistics().fillHistogram(("distance_hits_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), distance_hits(event).at(1) );
	    getStatistics().fillHistogram(("distance_hits_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), distance_hits(event).at(2) );

	    getStatistics().fillHistogram(("hitsDis_diff_afterTOT_"+ EventTypeNames.at(event_type)).c_str(), hitsDis_diff(event).at(0), hitsDis_diff(event).at(2) );
	    //}
	    */
	
	   if ( scatterTest( event ) > 20 ){
	  
	     getStatistics().fillHistogram( "HitMult_afterSTest", event_type);

	 //!!! toCheckEvtType ( event, time_window_mc,  event_type = ALL_THREE_FROM_SAME_VTX );

	     getStatistics().fillHistogram(("totHit1_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), TOT( eng ).at(0));
	     getStatistics().fillHistogram(("totHit2_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), TOT( eng ).at(1));
	     getStatistics().fillHistogram(("totHit3_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), TOT( eng ).at(2));

	     getStatistics().fillHistogram(("scattTest_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), scatterTest( event) );

	     getStatistics().fillHistogram(("dLOR_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), calculatedLOR( event));

	     getStatistics().fillHistogram(("sumDiff_azmTheta_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), (thetas[1]- thetas[0]));

	     getStatistics().fillHistogram(("angle_dLOR_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calculatedLOR( event));

	     getStatistics().fillHistogram(("angle_dLOR_fromSph_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calcLOR_sphRadius(event) );

	     getStatistics().fillHistogram(("angle_dLOR_fromCyd_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calcLOR_cydRadius(event) );

	     getStatistics().fillHistogram(("angle_dLOR_fromTril_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calcLOR_recsAnnhPt( event, decayPoint ).at(0));

	     getStatistics().fillHistogram(("dLOR_fromAnnhPt_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), calcLOR_recsAnnhPt( event, decayPoint ).at(0), calcLOR_recsAnnhPt( event, decayPoint ).at(2));
	     getStatistics().fillHistogram(("dLOR_fromAnnhPt_total_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), calcLOR_recsAnnhPt( event, decayPoint ).at(0), (calcLOR_recsAnnhPt( event, decayPoint ).at(0)+calcLOR_recsAnnhPt( event, decayPoint ).at(1)+calcLOR_recsAnnhPt( event, decayPoint ).at(2)));

	     getStatistics().fillHistogram(("shortDlor_comp_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), calcLOR_sphRadius(event), calcLOR_recsAnnhPt( event, decayPoint ).at(0));

	     getStatistics().fillHistogram(("distance_hits_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), distance_hits(event).at(0) );
	     getStatistics().fillHistogram(("distance_hits_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), distance_hits(event).at(1) );
	     getStatistics().fillHistogram(("distance_hits_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), distance_hits(event).at(2) );

	     getStatistics().fillHistogram(("hitsDis_diff_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), hitsDis_diff(event).at(0), hitsDis_diff(event).at(2) );

	     getStatistics().fillHistogram(("tot_data_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), tot_data(event).at(0) );
	     getStatistics().fillHistogram(("tot_data_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), tot_data(event).at(1) );
	     getStatistics().fillHistogram(("tot_data_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), tot_data(event).at(2) );

	     getStatistics().fillHistogram(("AnnhPointXYComp_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), sol[1].X(), sol[1].Y());
	     getStatistics().fillHistogram(("AnnhPointYZComp_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), sol[1].Y(), sol[1].Z());

	     if ( sol[1].Z() < 1. && sol[1].Z() > -1.){
	       getStatistics().fillHistogram(("AnnhPointXY_Zto0_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), sol[1].X(), sol[1].Y());
	     }

	     getStatistics().fillHistogram(("sum_diff_angles_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), (theta.at(0)+ theta.at(1))*TMath::RadToDeg(), (theta.at(1)- theta.at(0))*TMath::RadToDeg());

	     getStatistics().fillHistogram(("lifeTime_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), (decayTime_3hits(event,event_type) - pmtEmmTime)/1000);
	     getStatistics().fillHistogram(("lifeTime_ns_afterSTest_"+ EventTypeNames.at(event_type)).c_str(), (decayTime_3hits(event,event_type) - pmtEmmTime));


	     if ( calculatedLOR( event) > 5 ){

	       getStatistics().fillHistogram( "HitMult_afterDLOR", event_type);

	       getStatistics().fillHistogram(("tot_data_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), tot_data(event).at(0) );
	       getStatistics().fillHistogram(("tot_data_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), tot_data(event).at(1) );
	       getStatistics().fillHistogram(("tot_data_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), tot_data(event).at(2) );

	       getStatistics().fillHistogram(("sum_diff_angles_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), (theta.at(0)+ theta.at(1))*TMath::RadToDeg(), (theta.at(1)- theta.at(0))*TMath::RadToDeg());

	       getStatistics().fillHistogram(("scattTest_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), scatterTest( event ) );

	       getStatistics().fillHistogram(("dLOR_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), calculatedLOR( event));

	       getStatistics().fillHistogram(("sumDiff_azmTheta_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), (thetas[1]- thetas[0]));

	       getStatistics().fillHistogram(("angle_dLOR_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calculatedLOR( event));

	       getStatistics().fillHistogram(("angle_dLOR_fromSph_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calcLOR_sphRadius(event) );

	       getStatistics().fillHistogram(("angle_dLOR_fromCyd_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calcLOR_cydRadius(event) );

	       getStatistics().fillHistogram(("angle_dLOR_fromTril_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calcLOR_recsAnnhPt( event, decayPoint ).at(0));

	       getStatistics().fillHistogram(("shortDlor_comp_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), calcLOR_sphRadius(event), calcLOR_recsAnnhPt( event, decayPoint ).at(0));

	       getStatistics().fillHistogram(("distance_hits_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), distance_hits(event).at(0) );
	       getStatistics().fillHistogram(("distance_hits_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), distance_hits(event).at(1) );
	       getStatistics().fillHistogram(("distance_hits_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), distance_hits(event).at(2) );

	       getStatistics().fillHistogram(("hitsDis_diff_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), hitsDis_diff(event).at(0), hitsDis_diff(event).at(2) );

	       getStatistics().fillHistogram(("dLOR_fromAnnhPt_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), calcLOR_recsAnnhPt( event, decayPoint ).at(0), calcLOR_recsAnnhPt( event, decayPoint ).at(2));
	       getStatistics().fillHistogram(("dLOR_fromAnnhPt_total_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), calcLOR_recsAnnhPt( event, decayPoint ).at(0), (calcLOR_recsAnnhPt( event, decayPoint ).at(0)+calcLOR_recsAnnhPt( event, decayPoint ).at(1)+calcLOR_recsAnnhPt( event, decayPoint ).at(2)));

	       getStatistics().fillHistogram(("AnnhPointXYComp_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), sol[1].X(), sol[1].Y());
	       getStatistics().fillHistogram(("AnnhPointYZComp_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), sol[1].Y(), sol[1].Z());

	       if ( sol[1].Z() < 1. && sol[1].Z() > -1.){
		 getStatistics().fillHistogram(("AnnhPointXY_Zto0_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), sol[1].X(), sol[1].Y());
	       }

	       getStatistics().fillHistogram(("lifeTime_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), (decayTime_3hits(event,event_type) - pmtEmmTime)/1000);
	       getStatistics().fillHistogram(("lifeTime_ns_afterDLOR_"+ EventTypeNames.at(event_type)).c_str(), (decayTime_3hits(event,event_type) - pmtEmmTime));


	       
	       if ( (thetas[0]+ thetas[1]) > 180){

		 getStatistics().fillHistogram( "HitMult_after2DAngle", event_type);

		 getStatistics().fillHistogram(("tot_data_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), tot_data(event).at(0) );
		 getStatistics().fillHistogram(("tot_data_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), tot_data(event).at(1) );
		 getStatistics().fillHistogram(("tot_data_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), tot_data(event).at(2) );

		 getStatistics().fillHistogram(("sum_diff_angles_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), (theta.at(0)+ theta.at(1))*TMath::RadToDeg(), (theta.at(1)- theta.at(0))*TMath::RadToDeg());

		 getStatistics().fillHistogram(("scattTest_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), scatterTest( event ) );

		 getStatistics().fillHistogram(("dLOR_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), calculatedLOR( event));

		 getStatistics().fillHistogram(("sumDiff_azmTheta_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), (thetas[1]- thetas[0]));

		 getStatistics().fillHistogram(("angle_dLOR_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calculatedLOR( event));

		 getStatistics().fillHistogram(("angle_dLOR_fromSph_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calcLOR_sphRadius(event) );

		 getStatistics().fillHistogram(("angle_dLOR_fromCyd_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calcLOR_cydRadius(event) );

		 getStatistics().fillHistogram(("angle_dLOR_fromTril_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), (thetas[0]+ thetas[1]), calcLOR_recsAnnhPt( event, decayPoint ).at(0));

		 getStatistics().fillHistogram(("shortDlor_comp_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), calcLOR_sphRadius(event), calcLOR_recsAnnhPt( event, decayPoint ).at(0));

		 getStatistics().fillHistogram(("distance_hits_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), distance_hits(event).at(0) );
		 getStatistics().fillHistogram(("distance_hits_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), distance_hits(event).at(1) );
		 getStatistics().fillHistogram(("distance_hits_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), distance_hits(event).at(2) );

		 getStatistics().fillHistogram(("hitsDis_diff_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), hitsDis_diff(event).at(0), hitsDis_diff(event).at(2) );
		 getStatistics().fillHistogram(("dLOR_fromAnnhPt_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), calcLOR_recsAnnhPt( event, decayPoint ).at(0), calcLOR_recsAnnhPt( event, decayPoint ).at(2));
	       getStatistics().fillHistogram(("dLOR_fromAnnhPt_total_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), calcLOR_recsAnnhPt( event, decayPoint ).at(0), (calcLOR_recsAnnhPt( event, decayPoint ).at(0)+calcLOR_recsAnnhPt( event, decayPoint ).at(1)+calcLOR_recsAnnhPt( event, decayPoint ).at(2)));
	       
		 getStatistics().fillHistogram(("AnnhPointXYComp_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), sol[1].X(), sol[1].Y());
		 getStatistics().fillHistogram(("AnnhPointYZComp_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), sol[1].Y(), sol[1].Z());

		 if ( sol[1].Z() < 1. && sol[1].Z() > -1.){
		   getStatistics().fillHistogram(("AnnhPointXY_Zto0_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), sol[1].X(), sol[1].Y());
		 }

		 getStatistics().fillHistogram(("lifeTime_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), (decayTime_3hits(event,event_type) - pmtEmmTime)/1000);
		 getStatistics().fillHistogram(("lifeTime_ns_after2DAngle_"+ EventTypeNames.at(event_type)).c_str(), (decayTime_3hits(event,event_type) - pmtEmmTime));


		 getStatistics().fillHistogram(("AnnhPointXComp_"+ EventTypeNames.at(event_type)).c_str(), sol[1].X());
		 getStatistics().fillHistogram(("AnnhPointYComp_"+ EventTypeNames.at(event_type)).c_str(), sol[1].Y());
		 getStatistics().fillHistogram(("AnnhPointZComp_"+ EventTypeNames.at(event_type)).c_str(), sol[1].Z());

		 if ( sol[1].Z() < 10.1 && sol[1].Z() > -10.1){
		   getStatistics().fillHistogram(("AnnhPointXY_Zto10_"+ EventTypeNames.at(event_type)).c_str(), sol[1].X(), sol[1].Y());
		 }
      
	     
		 getStatistics().fillHistogram(("annh_radius_"+ EventTypeNames.at(event_type)).c_str(), radius);
		 getStatistics().fillHistogram(("annh_radius_trans_"+ EventTypeNames.at(event_type)).c_str(), transRad);
	    
		 getStatistics().getHisto1D(("annh_radius_jacobian_"+ EventTypeNames.at(event_type)).c_str())->Fill(radius, 1./radius);
		 getStatistics().fillHistogram(("annh_radius_3D_"+ EventTypeNames.at(event_type)).c_str(), decayPoint.Mag());
		 getStatistics().getHisto1D(("annh_radius3D_jacobian_"+ EventTypeNames.at(event_type)).c_str())->Fill(decayPoint.Mag() , 1./decayPoint.Mag());

	         // getStatistics().fillHistogram(("momUncer_"+ EventTypeNames.at(event_type)).c_str(), p1_uncer);

		 getStatistics().fillHistogram(("ang12_rndm_"+ EventTypeNames.at(event_type)).c_str(), theta_rndm[combs[wc][0]]*TMath::RadToDeg(), theta_rndm[combs[wc][1]]*TMath::RadToDeg());
		 getStatistics().fillHistogram(("ang23_rndm_"+ EventTypeNames.at(event_type)).c_str(), theta_rndm[combs[wc][1]]*TMath::RadToDeg(), theta_rndm[combs[wc][2]]*TMath::RadToDeg());
		 getStatistics().fillHistogram(("ang13_rndm_"+ EventTypeNames.at(event_type)).c_str(), theta[combs[wc][0]]*TMath::RadToDeg(), theta[combs[wc][1]]*TMath::RadToDeg());

		 getStatistics().fillHistogram(("th1th2_"+ EventTypeNames.at(event_type)).c_str(), angle_3D[0], angle_3D[1] );
		 getStatistics().fillHistogram(("th2th3_"+ EventTypeNames.at(event_type)).c_str(), angle_3D[1], angle_3D[2] );
		 getStatistics().fillHistogram(("th1th3_"+ EventTypeNames.at(event_type)).c_str(), angle_3D[0], angle_3D[2] );
      
		 getStatistics().fillHistogram(("E1E2_"+ EventTypeNames.at(event_type)).c_str(), energies[0], energies[1]);
		 getStatistics().fillHistogram(("E2E3_"+ EventTypeNames.at(event_type)).c_str(), energies[1], energies[2]);
		 getStatistics().fillHistogram(("E1E3_"+ EventTypeNames.at(event_type)).c_str(), energies[0], energies[2]);
                 //getStatistics().fillHistogram("eng_sum", energies[0]+energies[1]+energies[2]);

		 getStatistics().fillHistogram(("E1E2_rndm_"+ EventTypeNames.at(event_type)).c_str(), energies[combs[wc][0]], energies[combs[wc][1]]);
		 getStatistics().fillHistogram(("E2E3_rndm_"+ EventTypeNames.at(event_type)).c_str(), energies[combs[wc][1]], energies[combs[wc][2]]);
		 getStatistics().fillHistogram(("E1E3_rndm_"+ EventTypeNames.at(event_type)).c_str(), energies[combs[wc][0]], energies[combs[wc][2]]);

		 
		 getStatistics().fillHistogram(("angles12_"+ EventTypeNames.at(event_type)).c_str(), theta[0]*TMath::RadToDeg(), theta[1]*TMath::RadToDeg());
		 getStatistics().fillHistogram(("angles23_"+ EventTypeNames.at(event_type)).c_str(), theta[1]*TMath::RadToDeg(), theta[2]*TMath::RadToDeg());
		 getStatistics().fillHistogram(("angles13_"+ EventTypeNames.at(event_type)).c_str(), theta[0]*TMath::RadToDeg(), theta[2]*TMath::RadToDeg());
         

		 normaltoPlane = momenta.at( k.at(2).second).Cross( momenta.at( k.at(1).second));
		 spin = decayPoint.Unit();

		 OperatorStudy[0]  = spin.Dot(momenta.at(k.at(2).second).Unit());
		 OperatorStudy[1]  = spin.Dot(normaltoPlane.Unit());

		 getStatistics().fillHistogram(("S_k1_"+ EventTypeNames.at(event_type)).c_str(), OperatorStudy[0] );
		 getStatistics().fillHistogram(("S_k1_k2_"+ EventTypeNames.at(event_type)).c_str(), OperatorStudy[1]);
	
	         //getStatistics().fillHistogram("decayTime", t[1]/1000);
	
	       }
	     }

	}
      }
    }  
  }
}


double ThreeHitAnalysis::scatterTest( const JPetEvent *event  ) {

  double tDiff;
  double  c = 29.979246;   // cm/ns
  std::vector<double> deltaScatt;

  for(int i=0; i < 3; ++i){
      for(int j=i+1; j < 3; ++j){

	if (event->getHits().at(i).getTime() < event->getHits().at(j).getTime()) {
	  tDiff = (event->getHits().at(j).getTime() - event->getHits().at(i).getTime())/1000.0;
	} else {
	  tDiff = (event->getHits().at(i).getTime() - event->getHits().at(j).getTime())/1000.0;
	}

	float distance = sqrt(pow(event->getHits().at(j).getPosX()- event->getHits().at(i).getPosX(), 2) + pow(event->getHits().at(j).getPosY()- event->getHits().at(i).getPosY(), 2)
			      + pow(event->getHits().at(j).getPosZ()- event->getHits().at(i).getPosZ(), 2));

	deltaScatt.push_back(  fabs(distance - tDiff*c ));
        
      }
  }
  
  std::sort(deltaScatt.begin(), deltaScatt.end());
  return deltaScatt.at(0);
  
}

std::vector<double> ThreeHitAnalysis::TOT( std::vector<double> eng){ // EvtType& event_type ){
  
  double param[4], sumTOT = 0;
  std::vector<double> tot;
  
  param[0] =  29.84;  // ns
  param[1] =  2.446; // ns/keV
  param[2] =  310; // keV
  param[3] =  2.41; // ns/keV^2
  
  for ( int i = 0; i < 3; i++){
    tot.push_back( param[0] +  (param[1] - param[0])/ (1 + pow(eng.at(i)/param[2], param[3]))); 							      
  }
  
    return tot;
    tot.clear(); 
}

double ThreeHitAnalysis:: calculateAngle(JPetHit hit1, JPetHit hit2 ){
  
  double scalerProd = hit1.getPosX()*hit2.getPosX() + hit1.getPosY()*hit2.getPosY() ;
  double magnitude = sqrt((pow(hit1.getPosX(),2) + pow(hit1.getPosY(),2) )*(pow(hit2.getPosX(),2) + pow(hit2.getPosY(),2)));

  return acos(scalerProd/magnitude)*180/3.14159265;
}

double ThreeHitAnalysis:: calculatedLOR( const JPetEvent *event ){

  std::vector<Double_t> d_LOR;
  double dLOR_min;
  for(int i=0; i < 3; ++i){
      for(int j=i+1; j < 3; ++j){
	TVector3 vtx2g = EventCategorizerTools::calculateAnnihilationPoint( event->getHits().at(i), event->getHits().at(j));
        d_LOR.push_back( vtx2g.Mag());
      }
  }
  
  std::sort (d_LOR.begin(), d_LOR.end());
  dLOR_min = d_LOR[0];
  return dLOR_min;
    
}

double ThreeHitAnalysis:: calcLOR_sphRadius( const JPetEvent *event ){

  std::vector<Double_t> distance;
  double distance_min;
  for(int i=0; i < 3; ++i){
      for(int j=i+1; j < 3; ++j){
	TVector3 vtx2g = EventCategorizerTools::calculateAnnihilationPoint( event->getHits().at(i), event->getHits().at(j));
        distance.push_back (vtx2g.Mag() - 10.);
        
      }
  }
  
  std::sort (distance.begin(), distance.end());
  distance_min = distance.at(0);
  return distance_min;

}

double ThreeHitAnalysis:: calcLOR_cydRadius( const JPetEvent *event ){

  std::vector<Double_t> distance;
  double distance_min;
  for(int i=0; i < 3; ++i){
      for(int j=i+1; j < 3; ++j){
	TVector3 vtx2g = EventCategorizerTools::calculateAnnihilationPoint( event->getHits().at(i), event->getHits().at(j));
        distance.push_back (vtx2g.Mag() - 12.);
        
      }
  }
  
  std::sort (distance.begin(), distance.end());
  distance_min = distance.at(0);
  return distance_min;

}


std::vector<double> ThreeHitAnalysis:: calcLOR_recsAnnhPt( const JPetEvent *event, TVector3 decayPt ){

  std::vector<Double_t> distance;
  double distance_min;
  for(int i=0; i < 3; ++i){
      for(int j=i+1; j < 3; ++j){
	TVector3 vtx2g = EventCategorizerTools::calculateAnnihilationPoint( event->getHits().at(i), event->getHits().at(j));
        distance.push_back ((vtx2g - decayPt).Mag());
        
      }
  }
  
  std::sort (distance.begin(), distance.end());
  //distance_min = distance.at(0);
  return distance;

}

double ThreeHitAnalysis::calculateAngleSum (JPetHit hit1, JPetHit hit2, JPetHit hit3){
  
  std::vector<double> theta;
  theta.push_back(hit1.getBarrelSlot().getTheta());
  theta.push_back(hit2.getBarrelSlot().getTheta());
  theta.push_back(hit3.getBarrelSlot().getTheta());
  std::sort (theta.begin(), theta.end());

  double theta12 = theta[1] - theta[0];
  double theta23 = theta[2] - theta[1];
  double theta31 = 360 - theta12 - theta23;
  theta.clear();
  theta.push_back (theta12);
  theta.push_back (theta23);
  theta.push_back (theta31);
  std::sort (theta.begin(), theta.end());

  return (theta[0] + theta[1]);

}

std::vector<double> ThreeHitAnalysis::distance_hits (const JPetEvent *event) {

  std::vector<double> distance;
  
  for(int i=0; i < 3; ++i){
      for(int j=i+1; j < 3; ++j){
	
        distance.push_back( sqrt( pow(event->getHits().at(j).getPosX() - event->getHits().at(i).getPosX(),2) + pow( event->getHits().at(j).getPosY() - event->getHits().at(i).getPosY(),2) + pow( event->getHits().at(j).getPosZ() - event->getHits().at(i).getPosZ(),2) ));
	
      }
  }

  std::sort (distance.begin(), distance.end());
  
  return distance;
}


std::vector<double> ThreeHitAnalysis::hitsDis_diff (const JPetEvent *event ) {

  std::vector<double> distance_diff;
  
  for(int i=0; i < 3; ++i){
      for(int j=i+1; j < 3; ++j){
 
	distance_diff.push_back( fabs (distance_hits(event).at(i) - distance_hits(event).at(j)) );
      }
  }

  std::sort (distance_diff.begin(), distance_diff.end());
  
  return distance_diff;
}

double ThreeHitAnalysis:: calculateAngle_3D(JPetHit hit1, JPetHit hit2 ){
  
  double scalerProd = hit1.getPosX()*hit2.getPosY() + hit1.getPosY()*hit2.getPosY() ;
  double magnitude = sqrt((pow(hit1.getPosX(),2) + pow(hit1.getPosY(),2) )*(pow(hit2.getPosX(),2) + pow(hit2.getPosY(),2))*
			  (pow(hit1.getPosZ(),2) + pow(hit1.getPosZ(),2) ));

  return acos(scalerProd/magnitude)*180/3.14159265;
}

std::vector<double> ThreeHitAnalysis::tot_data ( const JPetEvent *event ) {

  std::vector<double> tot;
  
  for(int i=0; i < 3; ++i){
    
    tot.push_back( HitFinderTools::calculateTOT(event->getHits().at(i), HitFinderTools::TOTCalculationType::kThresholdTrapeze)/1000);
  }
  
  return tot;
  tot.clear();
}

std::vector<double> ThreeHitAnalysis::tot_adj ( const JPetEvent *event ) {

  double tot = 0;
  std::vector<double> totAdj;
  std::vector<JPetHit> hit;
  
  for(int i=0; i < 3; ++i){
    hit = event->getHits();

    auto sigALead = hit.at(i).getSignalA().getRecoSignal().getRawSignal().getPoints(JPetSigCh::Leading, JPetRawSignal::ByThrNum);
    auto sigBLead = hit.at(i).getSignalB().getRecoSignal().getRawSignal().getPoints(JPetSigCh::Leading, JPetRawSignal::ByThrNum);

    auto sigATrail = hit.at(i).getSignalA().getRecoSignal().getRawSignal().getPoints(JPetSigCh::Trailing, JPetRawSignal::ByThrNum);
    auto sigBTrail = hit.at(i).getSignalB().getRecoSignal().getRawSignal().getPoints(JPetSigCh::Trailing, JPetRawSignal::ByThrNum);

    if (sigALead.size()>0 && sigATrail.size()>0) {
      for( unsigned i = 0; i < sigALead.size() && i < sigATrail.size(); i++){
	switch(i)
	  {
	  case 0:
	    tot += (sigATrail.at(i).getValue() - sigALead.at(i).getValue()) * 3.0/11.0;
	    break;
	  case 1:
	    tot += (sigATrail.at(i).getValue() - sigALead.at(i).getValue()) * 5.0/11.0;
	    break;
	  case 2:
	    tot += (sigATrail.at(i).getValue() - sigALead.at(i).getValue());
	    break;
	  case 3:
	    tot += (sigATrail.at(i).getValue() - sigALead.at(i).getValue());
	    break;
	  }
      }
    }
    
    if (sigBLead.size() > 0 && sigBTrail.size() > 0) {
      for (unsigned i = 0; i < sigBLead.size() && i < sigBTrail.size(); i++) {
	switch( i )
	  {
	  case 0:
	    tot += (sigBTrail.at(i).getValue() - sigBLead.at(i).getValue()) * 3.0/11.0;
	    break;
	  case 1:
	    tot += (sigBTrail.at(i).getValue() - sigBLead.at(i).getValue()) * 5.0/11.0;
	    break;
	  case 2:
	    tot += (sigBTrail.at(i).getValue() - sigBLead.at(i).getValue());
	    break;
	  case 3:
	    tot += (sigBTrail.at(i).getValue() - sigBLead.at(i).getValue());
	    break;
	  }
      }
    }
    
    totAdj.push_back(tot);
  }
  return totAdj;
 
}


double ThreeHitAnalysis::decayTime_3hits (const JPetEvent *event, EvtType& event_type){

  double time;
  double vel = 29.9792458 ;
 
  fReconstructor->setGamma(0, event->getHits().at(0));
  fReconstructor->setGamma(1, event->getHits().at(1));
  fReconstructor->setGamma(2, event->getHits().at(2));
  
  int error = 0;
  TVector3 sol[2];
  double t[2];
  error = fReconstructor->getSolution(sol[1], t[1], 1);
  time = t[1];

  return time;

}


bool ThreeHitAnalysis::terminate() {
  INFO("Event analysis completed.");
  return true;
}
	
void ThreeHitAnalysis::initialiseHistograms()
{

  if (fIsMC == true){
  	EventTypeNames = EventTypeNames_MC;
    }
  else {
    EventTypeNames = EventTypeNames_data;
  }
   
  std::vector<std::pair<unsigned, std::string>> binLabels;
  for (auto& i: EventTypeNames){
    
    getStatistics().createHistogramWithAxes(new TH1D( ("HitMult_"+ i.second).c_str(), "Number of Hits", 11, -0.5, 10.5), "Number of hits", "Counts");
    binLabels.push_back(std::make_pair(i.first, i.second));

    getStatistics().createHistogramWithAxes(new TH1D( ("engDep_g1_"+ i.second).c_str(), "", 1290, -10.5, 1279.5), "Eng [keV]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("engDep_g2_"+ i.second).c_str(), "", 1290, -10.5, 1279.5), "Eng [keV]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("engDep_g3_"+ i.second).c_str(), "", 1290, -10.5, 1279.5), "Eng [keV]", "Counts");

    getStatistics().createHistogramWithAxes(new TH2D( ("distance_scat_"+ i.second).c_str(), "", 161, -0.5, 159.5, 21, -10.5, 10.5), "distance b/w two hits [cm]", "ScattTest [ns]");

    getStatistics().createHistogramWithAxes(new TH1D( ("totHit1_"+ i.second).c_str(), "", 41, -0.5, 40.5), "Time Over Threshold [ns]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("totHit2_"+ i.second).c_str(), "", 41, -0.5, 40.5), "Time Over Threshold [ns]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("totHit3_"+ i.second).c_str(), "", 41, -0.5, 40.5), "Time Over Threshold [ns]", "Counts");

    getStatistics().createHistogramWithAxes(new TH1D( ("totHit1_afterTOTcut_"+ i.second).c_str(), "", 41, -0.5, 40.5), "Time Over Threshold [ns]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("totHit2_afterTOTcut_"+ i.second).c_str(), "", 41, -0.5, 40.5), "Time Over Threshold [ns]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("totHit3_afterTOTcut_"+ i.second).c_str(), "", 41, -0.5, 40.5), "Time Over Threshold [ns]", "Counts");

    getStatistics().createHistogramWithAxes(new TH1D( ("totHit1_afterSTest_"+ i.second).c_str(), "", 41, -0.5, 40.5), "Time Over Threshold [ns]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("totHit2_afterSTest_"+ i.second).c_str(), "", 41, -0.5, 40.5), "Time Over Threshold [ns]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("totHit3_afterSTest_"+ i.second).c_str(), "", 41, -0.5, 40.5), "Time Over Threshold [ns]", "Counts");

    getStatistics().createHistogramWithAxes(new TH1D( ("scattTest_"+ i.second).c_str(), "Scatter Test", 120, -10.5, 119.5), "#delta d [cm]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("scattTest_afterTOT_"+ i.second).c_str(), "Scatter Test", 120, -10.5, 119.5), "#delta d [cm]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("scattTest_afterSTest_"+ i.second).c_str(), "Scatter Test", 120, -10.5, 119.5), "#delta d [cm]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("scattTest_afterDLOR_"+ i.second).c_str(), "Scatter Test", 120, -10.5, 119.5), "#delta d [cm]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("scattTest_after2DAngle_"+ i.second).c_str(), "Scatter Test", 120, -10.5, 119.5), "#delta d [cm]", "Counts");

    getStatistics().createHistogramWithAxes(new TH1D( ("dLOR_"+ i.second).c_str(), "", 120, -0.5, 59.5), "LOR [cm]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("dLOR_afterTOT_"+ i.second).c_str(), "", 120, -0.5, 59.5), " LOR [cm]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("dLOR_afterSTest_"+ i.second).c_str(), "", 120, -0.5, 59.5), " LOR [cm]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("dLOR_afterDLOR_"+ i.second).c_str(), "", 120, -0.5, 59.5), " LOR [cm]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("dLOR_after2DAngle_"+ i.second).c_str(), "", 120, -0.5, 59.5), " LOR [cm]", "Counts");

    getStatistics().createHistogramWithAxes(new TH1D( ("AnnhPointXComp_"+ i.second).c_str(), " Annihilation Point (X) ", 120, -60, 60), " AnnihilationPoint_X [cm]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("AnnhPointYComp_"+ i.second).c_str(), " Annihilation Point (Y) ", 120, -60, 60), " AnnihilationPoint_Y [cm]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("AnnhPointZComp_"+ i.second).c_str(), " Annihilation Point (Z) ", 120, -60, 60), " AnnihilationPoint_Z [cm]", "Counts");

    getStatistics().createHistogramWithAxes(new TH2D( ("AnnhPointYZComp_"+ i.second).c_str(), "Annihilation Point YZ ", 120, -60, 60, 120, -60, 60), " AnnihilationPoint_Y [cm]", "AnnihilationPoint_Z [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("AnnhPointXYComp_"+ i.second).c_str(), " Annihilation Point XY ", 120, -60, 60, 120, -60, 60), " AnnihilationPoint_X [cm]", "AnnihilationPoint_Y [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("AnnhPointXY_Zto0_"+ i.second).c_str(), "Annihilation Point XY ", 120, -60, 60, 120, -60, 60), " AnnihilationPoint_X [cm]", "AnnihilationPoint_Y [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("AnnhPointXY_Zto10_"+ i.second).c_str(), "Annihilation Point XY ", 120, -60, 60, 120, -60, 60), " AnnihilationPoint_X [cm]", "AnnihilationPoint_Y [cm]");

    getStatistics().createHistogramWithAxes(new TH2D( ("AnnhPointYZComp_afterTOT_"+ i.second).c_str(), "Annihilation Point YZ ", 120, -60, 60, 120, -60, 60), " AnnihilationPoint_Y [cm]", "AnnihilationPoint_Z [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("AnnhPointXYComp_afterTOT_"+ i.second).c_str(), " Annihilation Point XY ", 120, -60, 60, 120, -60, 60), " AnnihilationPoint_X [cm]", "AnnihilationPoint_Y [cm]");
    

    getStatistics().createHistogramWithAxes(new TH2D( ("AnnhPointYZComp_afterSTest_"+ i.second).c_str(), "Annihilation Point YZ ", 120, -60, 60, 120, -60, 60), " AnnihilationPoint_Y [cm]", "AnnihilationPoint_Z [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("AnnhPointXYComp_afterSTest_"+ i.second).c_str(), " Annihilation Point XY ", 120, -60, 60, 120, -60, 60), " AnnihilationPoint_X [cm]", "AnnihilationPoint_Y [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("AnnhPointXY_Zto0_afterSTest_"+ i.second).c_str(), "Annihilation Point XY ", 120, -60, 60, 120, -60, 60), " AnnihilationPoint_X [cm]", "AnnihilationPoint_Y [cm]");
    

    getStatistics().createHistogramWithAxes(new TH2D( ("AnnhPointYZComp_afterDLOR_"+ i.second).c_str(), "Annihilation Point YZ ", 120, -60, 60, 120, -60, 60), " AnnihilationPoint_Y [cm]", "AnnihilationPoint_Z [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("AnnhPointXYComp_afterDLOR_"+ i.second).c_str(), " Annihilation Point XY ", 120, -60, 60, 120, -60, 60), " AnnihilationPoint_X [cm]", "AnnihilationPoint_Y [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("AnnhPointXY_Zto0_afterDLOR_"+ i.second).c_str(), "Annihilation Point XY ", 120, -60, 60, 120, -60, 60), " AnnihilationPoint_X [cm]", "AnnihilationPoint_Y [cm]");

    getStatistics().createHistogramWithAxes(new TH2D( ("AnnhPointYZComp_after2DAngle_"+ i.second).c_str(), "Annihilation Point YZ ", 120, -60, 60, 120, -60, 60), " AnnihilationPoint_Y [cm]", "AnnihilationPoint_Z [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("AnnhPointXYComp_after2DAngle_"+ i.second).c_str(), " Annihilation Point XY ", 120, -60, 60, 120, -60, 60), " AnnihilationPoint_X [cm]", "AnnihilationPoint_Y [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("AnnhPointXY_Zto0_after2DAngle_"+ i.second).c_str(), "Annihilation Point XY ", 120, -60, 60, 120, -60, 60), " AnnihilationPoint_X [cm]", "AnnihilationPoint_Y [cm]");
    

    getStatistics().createHistogramWithAxes( new TH1F(("annh_radius_"+ i.second).c_str(), " ", 120, -60, 60), " radius [cm]", "Counts");
    getStatistics().createHistogramWithAxes( new TH1F(("annh_radius_trans_"+ i.second).c_str(), "Transverse radius of decay point ", 120, -60, 60), " radius [cm]", "Counts");
    getStatistics().createHistogram( new TH1F(("annh_radius_jacobian_"+ i.second).c_str(), "Transverse radius of decay point (weighted);" " radius [cm]; Counts", 120, -60, 60));
    getStatistics().createHistogramWithAxes( new TH1F(("annh_radius_3D_"+ i.second).c_str(), " ", 120, -60, 60), " radius [cm]", "Counts");
    getStatistics().createHistogram( new TH1F(("annh_radius3D_jacobian_"+ i.second).c_str(), "Radius of Annh point (weighted);" " radius [cm]; Counts", 120, -60, 60));

    getStatistics().createHistogramWithAxes( new TH2D (("angles12_"+ i.second).c_str(), "Angle12_vs_Angle23 ", 200, 0., 180., 200, 0., 180), "#theta_{12} [deg]", "#theta_{23} [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("angles23_"+ i.second).c_str(), "Angle23_vs_Angle13 ", 200, 0., 180., 200, 0., 180), "#theta_{23} [deg]", "#theta_{13} [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("angles13_"+ i.second).c_str(), "Angle12_vs_Angle13 ", 200, 0., 180., 200, 0., 180), "#theta_{12} [deg]", "#theta_{13} [deg]");

    getStatistics().createHistogramWithAxes( new TH2D (("ang12_rndm_"+ i.second).c_str(), "Angle12_vs_Angle23 ", 200, 0., 200., 200, 0., 200.), "#theta_{12} [deg]", "#theta_{23} [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("ang23_rndm_"+ i.second).c_str(), "Angle23_vs_Angle13 ", 200, 0., 200., 200, 0., 200.), "#theta_{23} [deg]", "#theta_{13} [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("ang13_rndm_"+ i.second).c_str(), "Angle12_vs_Angle13 ", 200, 0., 200., 200, 0., 200.), "#theta_{12} [deg]", "#theta_{13} [deg]");

    getStatistics().createHistogramWithAxes( new TH2D (("th1th2_"+ i.second).c_str(), "Angle12_vs_Angle23 ", 200, 0., 180., 200, 0., 180), "#theta_{12} [deg]", "#theta_{23} [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("th2th3_"+ i.second).c_str(), "Angle23_vs_Angle13 ", 200, 0., 180., 200, 0., 180), "#theta_{23} [deg]", "#theta_{13} [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("th1th3_"+ i.second).c_str(), "Angle12_vs_Angle13 ", 200, 0., 180., 200, 0., 180), "#theta_{12} [deg]", "#theta_{13} [deg]");

    getStatistics().createHistogramWithAxes( new TH2D (("sum_diff_angles_"+ i.second).c_str()," ", 301, -0.5, 299.5, 301, -0.5, 299.5), "#theta_{1} + #theta_{2} [deg]", "#theta_{1} - #theta_{2} [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("sum_diff_angles_afterTOT_"+ i.second).c_str()," ", 301, -0.5, 299.5, 301, -0.5, 299.5), "#theta_{1} + #theta_{2} [deg]", "#theta_{1} - #theta_{2} [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("sum_diff_angles_afterSTest_"+ i.second).c_str()," ", 301, -0.5, 299.5, 301, -0.5, 299.5), "#theta_{1} + #theta_{2} [deg]", "#theta_{1} - #theta_{2} [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("sum_diff_angles_afterDLOR_"+ i.second).c_str()," ", 301, -0.5, 299.5, 301, -0.5, 299.5), "#theta_{1} + #theta_{2} [deg]", "#theta_{1} - #theta_{2} [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("sum_diff_angles_after2DAngle_"+ i.second).c_str()," ", 301, -0.5, 299.5, 301, -0.5, 299.5), "#theta_{1}+#theta_{2} [deg]", "#theta_{1} - #theta_{2} [deg]");

    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_"+ i.second).c_str()," ", 301, -0.5, 299.5, 120, -0.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "min d_{LOR} [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_afterTOT_"+ i.second).c_str()," ", 301, -0.5, 299.5, 120, -0.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "min d_{LOR} [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_afterSTest_"+ i.second).c_str()," ", 301, -0.5, 299.5, 120, -0.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "min d_{LOR} [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_afterDLOR_"+ i.second).c_str()," ", 301, -0.5, 299.5, 120, -0.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "min d_{LOR} [cm]");
     getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_after2DAngle_"+ i.second).c_str()," ", 301, -0.5, 299.5, 120, -0.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "min d_{LOR} [cm]");

    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_fromSph_"+ i.second).c_str()," ", 301, -0.5, 299.5, 140, -10.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "d_LOR [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_fromSph_afterTOT_"+ i.second).c_str()," ", 301, -0.5, 299.5, 140, -10.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "d_LOR [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_fromSph_afterSTest_"+ i.second).c_str()," ", 301, -0.5, 299.5, 140, -10.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "d_LOR [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_fromSph_afterDLOR_"+ i.second).c_str()," ", 301, -0.5, 299.5, 140, -10.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "d_LOR [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_fromSph_after2DAngle_"+ i.second).c_str()," ", 301, -0.5, 299.5, 140, -10.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "d_LOR [cm]");

    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_fromCyd_"+ i.second).c_str()," ", 301, -0.5, 299.5, 140, -12.5, 59.5), "Sum of 2 smallest angles [deg]", "d_LOR [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_fromCyd_afterTOT_"+ i.second).c_str()," ", 301, -0.5, 299.5, 140, -12.5, 59.5), "Sum of 2 smallest angles [deg]", "d_LOR [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_fromCyd_afterSTest_"+ i.second).c_str()," ", 301, -0.5, 299.5, 140, -12.5, 59.5), "Sum of 2 smallest angles [deg]", "d_LOR [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_fromCyd_afterDLOR_"+ i.second).c_str()," ", 301, -0.5, 299.5, 140, -12.5, 59.5), "Sum of 2 smallest angles [deg]", "d_LOR [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_fromCyd_after2DAngle_"+ i.second).c_str()," ", 301, -0.5, 299.5, 140, -12.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "min d_{LOR} [cm]");

    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_fromTril_"+ i.second).c_str()," ", 301, -0.5, 299.5, 140, -10.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "min d_{LOR} [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_fromTril_afterTOT_"+ i.second).c_str()," ", 301, -0.5, 299.5, 140, -10.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "min d_{LOR} [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_fromTril_afterSTest_"+ i.second).c_str()," ", 301, -0.5, 299.5, 140, -10.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "min d_{LOR} [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_fromTril_afterDLOR_"+ i.second).c_str()," ", 301, -0.5, 299.5, 140, -10.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "min d_{LOR} [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("angle_dLOR_fromTril_after2DAngle_"+ i.second).c_str()," ", 301, -0.5, 299.5, 140, -10.5, 59.5), "#theta_{1} + #theta_{2} [deg]", "min d_{LOR} [cm]");

    getStatistics().createHistogramWithAxes( new TH2D (("dLOR_fromAnnhPt_"+ i.second).c_str()," ", 140, -10.5, 59.5, 340, -10.5, 159.5), "d_min [cm]", "d_max [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("dLOR_fromAnnhPt_afterTOT_"+ i.second).c_str()," ", 140, -10.5, 59.5, 340, -10.5, 159.5), "d_min [cm]", "d_max [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("dLOR_fromAnnhPt_afterSTest_"+ i.second).c_str()," ", 140, -10.5, 59.5, 340, -10.5, 159.5), "d_min [cm]", "d_max [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("dLOR_fromAnnhPt_afterDLOR_"+ i.second).c_str()," ", 140, -10.5, 59.5, 340, -10.5, 159.5), "d_min [cm]", "d_max [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("dLOR_fromAnnhPt_after2DAngle_"+ i.second).c_str()," ", 140, -10.5, 59.5, 340, -10.5, 159.5), "d_min [cm]", "d_max [cm]");

    getStatistics().createHistogramWithAxes( new TH2D (("dLOR_fromAnnhPt_total_"+ i.second).c_str()," ",  140, -10.5, 59.5, 600,-10.5, 289.5), "d_min [cm]", "d_total [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("dLOR_fromAnnhPt_total_afterTOT_"+ i.second).c_str()," ", 140, -10.5, 59.5, 600, -10.5, 289.5), "d_min [cm]", "d_total [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("dLOR_fromAnnhPt_total_afterSTest_"+ i.second).c_str()," ", 140, -10.5, 59.5, 600, -10.5, 289.5), "d_min [cm]", "d_total [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("dLOR_fromAnnhPt_total_afterDLOR_"+ i.second).c_str()," ", 140, -10.5, 59.5, 600, -10.5, 289.5), "d_min [cm]", "d_total [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("dLOR_fromAnnhPt_total_after2DAngle_"+ i.second).c_str()," ", 140, -10.5, 59.5, 600, -10.5, 289.5), "d_min [cm]", "d_total [cm]");

    getStatistics().createHistogramWithAxes( new TH2D (("shortDlor_comp_"+ i.second).c_str()," ",140, -10.5, 59.5 , 140, -10.5, 59.5), "d_LOR_SphRadius [cm]", "d_LOR_trilAnnhPt [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("shortDlor_comp_afterTOT_"+ i.second).c_str()," ",140, -10.5, 59.5 , 140, -10.5, 59.5), "d_LOR_SphRadius [cm]", "d_LOR_trilAnnhPt [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("shortDlor_comp_afterSTest_"+ i.second).c_str()," ",140, -10.5, 59.5 , 140, -10.5, 59.5), "d_LOR_SphRadius [cm]", "d_LOR_trilAnnhPt [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("shortDlor_comp_afterDLOR_"+ i.second).c_str()," ",140, -10.5, 59.5 , 140, -10.5, 59.5), "d_LOR_SphRadius [cm]", "d_LOR_trilAnnhPt [cm]");
    getStatistics().createHistogramWithAxes( new TH2D (("shortDlor_comp_after2DAngle_"+ i.second).c_str()," ",140, -10.5, 59.5 , 140, -10.5, 59.5), "d_LOR_SphRadius [cm]", "d_LOR_trilAnnhPt [cm]");

    getStatistics().createHistogramWithAxes(new TH1D( ("distance_hits_"+ i.second).c_str(), "", 160, -10.5, 149.5), "distance bw hits [cm]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("distance_hits_afterTOT_"+ i.second).c_str(), "", 160, -10.5, 149.5), "distance bw hits [cm]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("distance_hits_afterSTest_"+ i.second).c_str(), "", 160, -10.5, 149.5), "distance bw hits [cm]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("distance_hits_afterDLOR_"+ i.second).c_str(), "", 160, -10.5, 149.5), "distance bw hits [cm]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("distance_hits_after2DAngle_"+ i.second).c_str(), "", 160, -10.5, 149.5), "distance bw hits [cm]", "Counts");

    getStatistics().createHistogramWithAxes(new TH2D( ("hitsDis_diff_"+ i.second).c_str(), "", 160, -10.5, 149.5, 160, -10.5, 149.5), "min diff. b/w 2 hits distance [cm]", "max diff. b/w 2 hits distance [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("hitsDis_diff_afterTOT_"+ i.second).c_str(), "", 160, -10.5, 149.5, 160, -10.5, 149.5), "min diff. b/w 2 hits distance [cm]", "max diff. b/w 2 hits distance [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("hitsDis_diff_afterSTest_"+ i.second).c_str(), "", 160, -10.5, 149.5, 160, -10.5, 149.5), "min diff. b/w 2 hits distance [cm]", "max diff. b/w 2 hits distance [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("hitsDis_diff_afterDLOR_"+ i.second).c_str(), "", 160, -10.5, 149.5, 160, -10.5, 149.5), "min diff. b/w 2 hits distance [cm]", "max diff. b/w 2 hits distance [cm]");
    getStatistics().createHistogramWithAxes(new TH2D( ("hitsDis_diff_after2DAngle_"+ i.second).c_str(), "", 160, -10.5, 149.5, 160, -10.5, 149.5), "min diff. b/w 2 hits distance [cm]", "max diff. b/w 2 hits distance [cm]");

    getStatistics().createHistogramWithAxes( new TH2D (("sumDiff_azmTheta_"+ i.second).c_str()," ", 401, -0.5, 299.5, 401, -0.5, 299.5), "Sum of 2 smallest angles [deg]", "Diff between 2 smallest angles [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("sumDiff_azmTheta_afterTOT_"+ i.second).c_str()," ", 401, -0.5, 299.5, 401, -0.5, 299.5), "Sum of 2 smallest angles [deg]", "Diff between 2 smallest angles [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("sumDiff_azmTheta_afterSTest_"+ i.second).c_str()," ", 401, -0.5, 299.5, 401, -0.5, 299.5), "Sum of 2 smallest angles [deg]", "Diff between 2 smallest angles [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("sumDiff_azmTheta_afterDLOR_"+ i.second).c_str()," ", 401, -0.5, 299.5, 401, -0.5, 299.5), "Sum of 2 smallest angles [deg]", "Diff between 2 smallest angles [deg]");
    getStatistics().createHistogramWithAxes( new TH2D (("sumDiff_azmTheta_after2DAngle_"+ i.second).c_str()," ", 401, -0.5, 299.5, 401, -0.5, 299.5), "Sum of 2 smallest angles [deg]", "Diff between 2 smallest angles [deg]");

    getStatistics().createHistogramWithAxes(new TH1D( ("lifeTime_"+ i.second).c_str(), "LifeTime", 101, -50.5, 50.5), "#delta t [mu sec]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("lifeTime_afterSTest_"+ i.second).c_str(), "LifeTime", 101, -50.5, 50.5), "#delta t [mu sec]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("lifeTime_ns_afterSTest_"+ i.second).c_str(), "LifeTime", 1001, -500.5, 500.5), "#delta t [ns]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("lifeTime_afterDLOR_"+ i.second).c_str(), "LifeTime", 101, -50.5, 50.5), "#delta t [mu sec]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("lifeTime_ns_afterDLOR_"+ i.second).c_str(), "LifeTime", 1001, -500.5, 500.5), "#delta t [ns]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("lifeTime_after2DAngle_"+ i.second).c_str(), "LifeTime", 101, -50.5, 50.5), "#delta t [mu sec]", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D( ("lifeTime_ns_after2DAngle_"+ i.second).c_str(), "LifeTime", 1001, -500.5, 500.5), "#delta t [ns]", "Counts");

    getStatistics().createHistogramWithAxes(new TH1D(("tot_data_"+ i.second).c_str(), "TOT", 341, -10.25, 160.25), "TOT [ns] ", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D(("tot_adj_"+ i.second).c_str(), "TOT", 341, -10.25, 160.25), "TOT [ns] ", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D(("tot_data_afterTOT_"+ i.second).c_str(), "TOT", 341, -10.25, 160.25), "TOT [ns] ", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D(("tot_data_afterSTest_"+ i.second).c_str(), "TOT", 341, -10.25, 160.25), "TOT [ns] ", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D(("tot_data_afterDLOR_"+ i.second).c_str(), "TOT", 341, -10.25, 160.25), "TOT [ns] ", "Counts");
    getStatistics().createHistogramWithAxes(new TH1D(("tot_data_after2DAngle_"+ i.second).c_str(), "TOT", 341, -10.25, 160.25), "TOT [ns] ", "Counts");
         
     
    getStatistics().createHistogramWithAxes( new TH2D (("E1E2_"+ i.second).c_str(), "E1 vs E2 ", 301, -0.5, 599.5, 301, -0.5, 599.5), "E_{1} [keV]", "E_{2} [keV]");
    getStatistics().createHistogramWithAxes( new TH2D (("E2E3_"+ i.second).c_str(), "E2 vs E3 ", 301, -0.5, 599.5, 301, -0.5, 599.5 ), "E_{2} [keV]", "E_{3} [keV]");
    getStatistics().createHistogramWithAxes( new TH2D (("E1E3_"+ i.second).c_str(), "E1 vs E3 ", 301, -0.5, 599.5, 301, -0.5, 599.5 ), "E_{1} [keV]", "E_{3} [keV]");

    getStatistics().createHistogramWithAxes( new TH2D (("E1E2_rndm_"+ i.second).c_str(), "E1 vs E2 ", 301, -0.5, 599.5, 301, -0.5, 599.5), "E_{1} [keV]", "E_{2} [keV]");
    getStatistics().createHistogramWithAxes( new TH2D (("E2E3_rndm_"+ i.second).c_str(), "E2 vs E3 ", 301, -0.5, 599.5, 301, -0.5, 599.5 ), "E_{2} [keV]", "E_{3} [keV]");
    getStatistics().createHistogramWithAxes( new TH2D (("E1E3_rndm_"+ i.second).c_str(), "E1 vs E3 ", 301, -0.5, 599.5, 301, -0.5, 599.5 ), "E_{1} [keV]", "E_{3} [keV]");

    getStatistics().createHistogramWithAxes( new TH1D(("S_k1_"+ i.second).c_str(), "S.k1 ", 100, -1, 1), "S.k1", "Counts");
    getStatistics().createHistogramWithAxes( new TH1D(("S_k1_k2_"+ i.second).c_str(), "S.(k1*k2) ", 100, -1, 1), "S.(k1*k2) ", "Counts");

    getStatistics().createHistogramWithAxes( new TH1D(("momUncer_"+ i.second).c_str(), " ", 40, -20, 20), " mom_uncertainity", "Counts");
					     
  }

  getStatistics().createHistogramWithAxes(new TH1D("tot_allHits", "TOT", 851, -10.25, 160.25), "TOT [ns] ", "Counts");
  getStatistics().createHistogramWithAxes(new TH1D("tot_1HitEvt", "TOT", 851, -10.25, 160.25), "TOT [ns] ", "Counts");
  
  getStatistics().createHistogramWithAxes(new TH2D( "HitMult_ScinID", "Number of Hits v/s scinID", 200, -0.5, 199.5,  11, -0.5, 10.5), "Scin_ID", "Number of Hits");
  
  
  getStatistics().createHistogramWithAxes(new TH1D( "HitMult", "Number of hits", 24, 0.5, 24.5), "Event Type", "Counts");
  getStatistics().setHistogramBinLabel("HitMult", getStatistics().AxisLabel::kXaxis, binLabels);
  
  getStatistics().createHistogramWithAxes(new TH1D( "HitMult_beforeEng", "Number of hits", 24, 0.5, 24.5), "Event Type", "Counts");
  getStatistics().setHistogramBinLabel("HitMult_beforeEng", getStatistics().AxisLabel::kXaxis, binLabels);
					     
  getStatistics().createHistogramWithAxes(new TH1D( "HitMult_afterEng", "Number of hits", 24, 0.5, 24.5), "Event Type", "Counts");
  getStatistics().setHistogramBinLabel("HitMult_afterEng", getStatistics().AxisLabel::kXaxis, binLabels);

  getStatistics().createHistogramWithAxes(new TH1D( "HitMult_afterTOT", "Number of hits", 24, 0.5, 24.5), "Event Type", "Counts");
  getStatistics().setHistogramBinLabel("HitMult_afterTOT", getStatistics().AxisLabel::kXaxis, binLabels);

  getStatistics().createHistogramWithAxes(new TH1D( "HitMult_afterSTest", "Number of hits", 24, 0.5, 24.5), "Event Type", "Counts");
  getStatistics().setHistogramBinLabel("HitMult_afterSTest", getStatistics().AxisLabel::kXaxis, binLabels);

  getStatistics().createHistogramWithAxes(new TH1D( "HitMult_afterDLOR", "Number of hits", 24, 0.5, 24.5), "Event Type", "Counts");
  getStatistics().setHistogramBinLabel("HitMult_afterDLOR", getStatistics().AxisLabel::kXaxis, binLabels);

  getStatistics().createHistogramWithAxes(new TH1D( "HitMult_after2DAngle", "Number of hits", 24, 0.5, 24.5), "Event Type", "Counts");
  getStatistics().setHistogramBinLabel("HitMult_after2DAngle", getStatistics().AxisLabel::kXaxis, binLabels);

  //getStatistics().createHistogramWithAxes(new TH1D( "HitMult_afterScaTest", "Number of hits", 10, 0.5, 10.5), "Event Type", "Counts");
  //getStatistics().setHistogramBinLabel("HitMult_afterScaTest", getStatistics().AxisLabel::kXaxis, binLabels);

  getStatistics().createHistogramWithAxes( new TH1D("error_tri", "error_tri", 100, -0.5, 9.5), "error ", "Counts");

  getStatistics().createHistogramWithAxes( new TH2D ("genGammaMult_g1", "GenGammaMult", 250, -1, 249, 250, -1, 249), "GenGammaMult_g1", " GenGammaMult_g2g3" );
  getStatistics().createHistogramWithAxes( new TH2D ("genGammaMult_g2", "GenGammaMult", 250, -1, 249, 250, -1, 249), "GenGammaMult_g2", " GenGammaMult_g1g3" );
  getStatistics().createHistogramWithAxes( new TH2D ("genGammaMult_g3", "GenGammaMult", 250, -1, 249, 250, -1, 249), "GenGammaMult_g3", " GenGammaMult_g1g2" );

  getStatistics().createHistogramWithAxes( new TH2D ("hit_genMult", " ", 8, -2, 6 , 250, -1, 249), "Hit", " genGammaMult" );
  getStatistics().createHistogramWithAxes( new TH2D ("hit_genMult_sameVtx", " ", 8, -2, 6 , 250, -1, 249), "Hit", " genGammaMult" );
  getStatistics().createHistogramWithAxes( new TH2D ("hit_genMult_2gSameVtx", " ", 8, -2, 6 , 250, -1, 249), "Hit", " genGammaMult" );

  //getStatistics().createHistogramWithAxes( new TH2D ("genMult1_sameVtx_g1g2", " ", 250, -1, 249 , 250, -1, 249), "genGammaMult_g1 (small)", " genGammaMult_g2" );
  //getStatistics().createHistogramWithAxes( new TH2D ("genMult2_sameVtx_g1g2", " ", 250, -1, 249 , 250, -1, 249), "genGammaMult_g1 (small)", " genGammaMult_g2" );
  //getStatistics().createHistogramWithAxes( new TH2D ("genMult3_sameVtx_g1g2", " ", 250, -1, 249 , 250, -1, 249), "genGammaMult_g1 (small)", " genGammaMult_g2" );

  getStatistics().createHistogramWithAxes( new TH2D ("genMult1_sameVtx_g2g3", " ", 250, -1, 249 , 250, -1, 249), "genGammaMult_g2 (small)", " genGammaMult_g3" );
  getStatistics().createHistogramWithAxes( new TH2D ("genMult2_sameVtx_g2g3", " ", 250, -1, 249 , 250, -1, 249), "genGammaMult_g2 (small)", " genGammaMult_g3" );
  getStatistics().createHistogramWithAxes( new TH2D ("genMult3_sameVtx_g2g3", " ", 250, -1, 249 , 250, -1, 249), "genGammaMult_g2 (small)", " genGammaMult_g3" );

  getStatistics().createHistogramWithAxes( new TH2D ("genMult1_sameVtx_g2g3_afterDLOR", " ", 250, -1, 249 , 250, -1, 249), "genGammaMult_g2 (small)", " genGammaMult_g3" );
  
  getStatistics().createHistogramWithAxes( new TH1D ("HitMult_vtx", " ", 3, 0.5, 3.5), "", " Counts" );
  
  getStatistics().createHistogramWithAxes( new TH2D ("hit_vtx", " ", 8, -2, 6, 2500, 0, 5000), "", " Counts" );

  getStatistics().createHistogramWithAxes(new TH1D( "lifeTime_3Hits", "", 101, -50.5, 50.5), "LifeTime [micro sec]", "Counts");
  
  std::vector<std::pair<unsigned, std::string>> bin_label;
  bin_label.push_back(std::make_pair( 1, "vtx1 = vtx2 = vtx3" ));
  bin_label.push_back(std::make_pair( 2, "vtx1 = vtx2"));
  bin_label.push_back(std::make_pair( 3,  " no common vtx"));
  getStatistics().setHistogramBinLabel("HitMult_vtx", getStatistics().AxisLabel::kXaxis, bin_label);  

      
}
