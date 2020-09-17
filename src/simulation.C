// @(#)root/simulation
// Author: Alberto Perro 15/11/19

/* COLLIDER EVENT MONTECARLO SIMULATION
 *  Generating particles collisions and intersecting with two
 *  layers of vertex detectors and saving hits to file
 */

#include <Riostream.h>
#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "lib/detector.h"
#include "lib/collision.h"
#include "lib/hit.h"

#ifndef N_EVENTS
#define N_EVENTS 1e6
#endif

int simulation(bool multScat = false, int randomNoise = 0)
{
   if(multScat)printf("MULTIPLE SCATTERING ENABLED!\n");
   int        debug = 1;
   TStopwatch watch;
   // creating TTree for simulation
   TFile *simulation = new TFile("simulation.root", "RECREATE");
   TTree *hitTree    = new TTree("hits", "Hits tree");

   if (debug > 0) printf("instantiating detector class\n");
   detector *det = new detector();
   if (debug > 0) printf("instantiating collision class\n");
   collision *vtx = new collision();
   double     ptc[4]; // vertex coordinates array
   // instantiating TClonesArrays to store hits
   TClonesArray *hitsL1 = new TClonesArray("hit", 100);
   TClonesArray *hitsL2 = new TClonesArray("hit", 100);
   TClonesArray &L1HITS = *hitsL1;
   TClonesArray &L2HITS = *hitsL2;
   // referring tree branches to arrays
   hitTree->Branch("Vertex", &ptc[0], "x0/D:y0:z0:m");
   hitTree->Branch("L1hit", &hitsL1);
   hitTree->Branch("L2hit", &hitsL2);
   //----------------------------- MONTECARLO ---------------------------
   for (int i = 0; i < N_EVENTS; i++) {
      if (debug > 3) printf("\ngenerating Collision %d\n", i);
      vtx->generateCollision(KINEMATIC);
      vtx->getCoordinates(ptc);
      if (debug > 2) printf("\rx: %f\t y: %f\t z: %f\t m: %f", ptc[0], ptc[1], ptc[2], ptc[3]);
      // generating particles with a given multiplicity distribution
      particle *part = new particle(); // direction holder
      uint8_t   L1 = 0, L2 = 0;        // number of hits
      for (uint8_t mult = 0; mult < ptc[3]; mult++) {
         vtx->getDir(part); // update the direction
         double beampt[3];  // hit coordinates in cartesian
         det->intersect(ptc, part, beampt, 1000., 8.0, 0.12, 0.03, multScat, 0.008);
         double l1pt[3];
         hit *  l1hit = det->intersect(beampt, part, l1pt, 270., 40., 0.012, 0.003, multScat, 0.008);
         if (l1hit != NULL) {
            L1HITS[L1++] = l1hit;
         }
         double l2pt[3];
         hit *  l2hit = det->intersect(l1pt, part, l2pt, 270., 70., 0.012, 0.003, false, 0);
         if (l2hit != NULL) {
            L2HITS[L2++] = l2hit;
         }
         if (l2hit != NULL) {
         //printf("%f,%f  %f,%f,%f\n",l2hit->getTheta(), TMath::ATan(beampt[1] / beampt[0]), (beampt[1] / beampt[0]), (l1pt[1] / l1pt[0]), l2pt[1] / l2pt[0]);
         }
      }
      hitTree->Fill();
      hitsL1->Clear();
      hitsL2->Clear();
      //printf("\nL1 hits: %2d\t L2 hits: %2d\n", L1, L2);
      printf("\rSimulation Progress: %u%%", (uint8_t)(i * 100 / N_EVENTS + 1));
   }
   if (debug > 0) printf("\nmontecarlo finished\n");
   // saving simulation file
   simulation->Write("", TObject::kOverwrite);
   simulation->Close();
   hitsL1->Delete();
   hitsL2->Delete();
   delete vtx;
   delete det;
   if (debug > 0) printf("simulation run succesful, cpu time: %.2f s\n", watch.CpuTime());
   return 0;
}