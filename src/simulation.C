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

#include "lib/detector.h"
#include "lib/collision.h"
#include "lib/hit.h"

#ifndef N_EVENTS
#define N_EVENTS 1e6
#endif

int simulation(bool multScat = false, int randomNoise = 0)
{
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
   TClonesArray *hitsL1    = new TClonesArray("hit", 100);
   TClonesArray *hitsL2    = new TClonesArray("hit", 100);
   TClonesArray *particles = new TClonesArray("particle", 100);
   // referring tree branches to arrays
   hitTree->Branch("Vertex", &ptc[0], "x0/D:y0:z0:m");
   hitTree->Branch("L1hit", &hitsL1);
   hitTree->Branch("L2hit", &hitsL2);
   //----------------------------- MONTECARLO ---------------------------
   for (int i = 0; i < N_EVENTS; i++) {
      if (debug > 3) printf("\ngenerating Collision %d\n", i);
      vtx->generateVertex(ptc);
      // generating particles with a given multiplicity distribution
      vtx->generateParticles(ptc, particles, KINEMATIC);
      if (debug > 2) printf("\rx: %f\t y: %f\t z: %f\t m: %f", ptc[0], ptc[1], ptc[2], ptc[3]);
      // beampipe scattering
      if (multScat) det->multipleScattering(particles, 0.001);
      if (debug > 3) printf("intersecting layer 1\n");
      // layer 1 intersection
      det->intersect(ptc, particles, hitsL1, 270., 40., 0.12, 0.03, randomNoise);
      // layer 1 scattering
      if (multScat) det->multipleScattering(particles, 0.001);
      if (debug > 3) printf("intersecting layer 2\n");
      // layer 2 intersection
      det->intersect(ptc, particles, hitsL2, 270., 70., 0.12, 0.03, randomNoise);
      // laier 2 scattering
      if (multScat) det->multipleScattering(particles, 0.001);
      // filling TTree
      hitTree->Fill();
      hitsL1->Clear();
      hitsL2->Clear();
      printf("\rSimulation Progress: %u%%", (uint8_t)(i * 100 / N_EVENTS + 1));
   }
   if (debug > 0) printf("\nmontecarlo finished\n");
   // saving simulation file
   simulation->Write("", TObject::kOverwrite);
   simulation->Close();
   delete vtx;
   delete det;
   if (debug > 0) printf("simulation run succesful, cpu time: %.2f s\n", watch.CpuTime());
   return 0;
}
