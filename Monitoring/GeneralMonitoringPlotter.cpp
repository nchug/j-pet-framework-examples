/**
 *  @copyright Copyright 2020 The J-PET Framework Authors. All rights reserved.
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
 *  @file GeneralMonitoringPlotter.cpp
 */

#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetWriter/JPetWriter.h>
#include "../LargeBarrelAnalysis/EventCategorizerTools.h"
#include "../LargeBarrelAnalysis/HitFinderTools.h"
#include "GeneralMonitoringPlotter.h"
#include <iostream>

using namespace jpet_options_tools;
using namespace std;

GeneralMonitoringPlotter::GeneralMonitoringPlotter(const char* name): JPetUserTask(name) {}

GeneralMonitoringPlotter::~GeneralMonitoringPlotter() {}

bool GeneralMonitoringPlotter::init()
{
  INFO("Event categorization started.");

  // Input events type
  fOutputEvents = new JPetTimeWindow("JPetEvent");
  // Initialise hisotgrams
  initialiseHistograms();
  return true;
}

bool GeneralMonitoringPlotter::exec()
{
  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {
    vector<JPetEvent> events;
    for (uint i = 0; i < timeWindow->getNumberOfEvents(); i++) {
      const auto& event = dynamic_cast<const JPetEvent&>(timeWindow->operator[](i));

      const std::vector<JPetHit>& hits = event.getHits();    

      analyzeSingleHits(hits);

      
      
      events.push_back(event);
    }
    saveEvents(events);
  } else { return false; }
  return true;
}

bool GeneralMonitoringPlotter::terminate()
{
  INFO("Event categorization completed.");
  return true;
}

void GeneralMonitoringPlotter::saveEvents(const vector<JPetEvent>& events)
{
  for (const auto& event : events) { fOutputEvents->add<JPetEvent>(event); }
}

void GeneralMonitoringPlotter::initialiseHistograms(){

  getStatistics().createHistogramWithAxes(new TH1D("hitCount_vs_ID", "HitCounts_vs_ID", 200, -0.5, 199.5), "Scin_ID", "Number of hits");
  getStatistics().createHistogramWithAxes(new TH1D("TOT_singleHit", "TOT of all hits", 340, -10.5, 159.5), "TOT [ns] ", "Counts");
  getStatistics().createHistogramWithAxes(new TH2D("TOT_vs_ID", "TOT of all hits vs Scin ID",  200, -0.5, 199.5, 340, -10.5, 159.5), "Scin_ID","TOT [ns]");
  getStatistics().createHistogramWithAxes(new TH2D("XY_plane_singleHit", "XY plane of scintillators", 400, -100, 100, 400, -100, 100), "Y axis", "X axis");

  // General histograms
  // getStatistics().createHistogramWithAxes(
  //   new TH2D("All_XYpos", "Hit position XY", 240, -60.25, 59.75, 240, -60.25, 59.75),
  //   "Hit X position [cm]", "Hit Y position [cm]"
  // );


}

void GeneralMonitoringPlotter::analyzeSingleHits(const std::vector<JPetHit>& hits) {

   if (hits.size() > 0 ) {
     
     for (int j = 0; j < hits.size(); j++){

       JPetHit singleHit = hits.at(j);
       double tot = HitFinderTools::calculateTOT(singleHit, HitFinderTools::TOTCalculationType::kThresholdTrapeze)/1000;
       
       getStatistics().fillHistogram("hitCount_vs_ID", singleHit.getScintillator().getID());
       getStatistics().fillHistogram( "TOT_singleHit", tot );
       getStatistics().fillHistogram( "TOT_vs_ID",  singleHit.getScintillator().getID(), tot );
       getStatistics().fillHistogram( "XY_plane_singleHit", singleHit.getPosY(), singleHit.getPosX() );
     }
  }
   
}
