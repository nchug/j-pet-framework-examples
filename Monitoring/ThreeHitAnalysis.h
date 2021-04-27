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
 *  @file ThreeHitAnalysis.h
 */

#ifndef THREEHITANALYSIS_H
#define THREEHITANALYSIS_H

#include <JPetCommonTools/JPetCommonTools.h>
#include <JPetUserTask/JPetUserTask.h>
#include "../LargeBarrelAnalysis/EventCategorizerTools.h"
#include "reconstructor.h"
#include <JPetEvent/JPetEvent.h>
#include <JPetHit/JPetHit.h>
#include <JPetMCHit/JPetMCHit.h>
#include <JPetTimer/JPetTimer.h>
#include <vector>
#include <map>

using namespace std;
class JPetWriter;

#ifdef __CINT__
#define override
#endif

static const double kLightVelocity_cm_ns = 29.9792458;

enum EvtType{

  UNKNOWN = 1,
  OPS,
  PPS,
  OPS_AND_1SCATTER,
  B2B_AND_SCATTER,
  B2B_AND_PROMPT,
  OPS_AND_PROMPT,
  OPS_AND_SCATTER_PROMPT,
  OPS_SCATTER_AND_PROMPT,
  B2B_AND_2SCATTER,
  OPS_AND_2SCATTER,
  B2B_SCATTER_AND_PROMPT,
  B2B_PROMPT_AND_SCATTER,
  PROMPT_SCATTER,
  B2B_PRMT_2SCAT,
  PRMT_B2B_2SCAT,
  ALL_THREE_FROM_SAME_VTX,
  TWO_FROM_SAME_VTX,
  NONE_FROM_SAME_VTX,
  OTHER  
  
};



class ThreeHitAnalysis : public JPetUserTask{
public:
	ThreeHitAnalysis(const char * name);
	virtual ~ThreeHitAnalysis(){};
	virtual bool init() override;
	virtual bool exec() override;
	virtual bool terminate() override;
		

	double meanTime_3Hits;
	
	double Magnitude(TVector3 p1);
	void saveEvents(const std::vector<JPetEvent>& event);
	
	void initialiseHistograms();
	void genOPSCalculation( const JPetEvent *event, const JPetTimeWindowMC& time_window_mc, EvtType& event_type);
	void pPs3G( const JPetEvent *event, const JPetTimeWindowMC& time_window_mc, EvtType& event_type);
	//void ReconHitsCalculation( const JPetEvent *event);
        void ReconHitsCalculation( const JPetEvent *event, EvtType& event_type);
	double scatterTest( const JPetEvent *event );
	std::vector<double> TOT(std::vector<double> eng);
	double calculateAngle(JPetHit hit1, JPetHit hit2 );
	double calculatedLOR( const JPetEvent *event);
	double calcLOR_sphRadius( const JPetEvent *event);
	double calcLOR_cydRadius( const JPetEvent *event);
	std::vector<double> calcLOR_recsAnnhPt( const JPetEvent *event, TVector3 decayPt);
	double calculateAngleSum (JPetHit hit1, JPetHit hit2, JPetHit hit3);
	std::vector<double> distance_hits (const JPetEvent *event);
	std::vector<double> hitsDis_diff (const JPetEvent *event);
	double calculateAngle_3D(JPetHit hit1, JPetHit hit2 );
	std::vector<double> tot_data (const JPetEvent *event);
	double decayTime_3hits(const JPetEvent *event , EvtType& event_type);
	std::vector<double> tot_adj (const JPetEvent *event);
	
	//void toCheckEvtType (const JPetEvent *event, const JPetTimeWindowMC& time_window_mc, EvtType& event_type);
        

private:
	double decayTime, lifeTime, pmtEmmTime;
	
   
protected:

        Reconstructor * fReconstructor;
	
	//vector<string> oPsType = {"gen", "recs"};

	const std::string kMC = "Save_MC_bool";
	bool fIsMC = false;

	std::map<EvtType, std::string> EventTypeNames_MC= {
	  {UNKNOWN, "3 hit evts"},
	  {OPS, "oPs"},
	  {PPS, "pPs"},
	  {OPS_AND_1SCATTER, "oPs and scattered photon"},
	  {B2B_AND_SCATTER, "back-to-back and scattered photon"},
	  {B2B_AND_PROMPT, "back-to-back and prompt photon"},
	  {OPS_AND_PROMPT, "oPs and prompt photon"},
	  {OPS_AND_SCATTER_PROMPT, "o-Ps, prompt, sactter prompt"},
	  {OPS_SCATTER_AND_PROMPT, "o-Ps, scatter, prompt"},
	  {OPS_AND_2SCATTER, "oPs and 2 scattered photon"},
	  {B2B_AND_2SCATTER, "back-to-back and 2 scattered photon"},
	  {B2B_SCATTER_AND_PROMPT, "b2b, b2b scattered and prompt"},
	  {B2B_PROMPT_AND_SCATTER, "b2b, prompt and prompt scatter"},
	  {PROMPT_SCATTER, "prompt and 2 scattered prompt"},
	  {PRMT_B2B_2SCAT, "prmt, b2b, b2b-doubleScat"},
	  {B2B_PRMT_2SCAT, "b2b, prmt, prmt-doubleScat"},
	  {ALL_THREE_FROM_SAME_VTX, "3g from same vtx"},
	  {TWO_FROM_SAME_VTX, "Any two from same vtx"},
	  {NONE_FROM_SAME_VTX, "all from different vtx"},
	  {OTHER, "other evts"}	  
	  
	};

	std::map<EvtType, std::string> EventTypeNames_data = {
	  {UNKNOWN, "3 hit evts"}
	};
	
	EvtType event_type;
	std::map<EvtType, std::string> EventTypeNames;

	//TVector3 p1_gen , p1_recs ;
	//double p1_uncer;
       	
	  
};

#endif /* !THREEHITANALYSIS_H */
