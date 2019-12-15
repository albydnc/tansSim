// @(#)root/simulation
// Author: Alberto Perro 15/11/19
#include "lib/detector.h"
#include "lib/collision.h"
#include "lib/hit.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include <Riostream.h>
#include <vector>
#include <array>
#define N_EVENTS 1e6
int simulation(bool multScat = false, int randomNoise = 0)
{
  int debug = 1;
  TStopwatch watch;
  TFile * simulation = new TFile("simulation.root","RECREATE");
  TTree * hitTree =  new TTree("hits", "Hits tree");
   if(debug>0)printf("instantiating detector\n");
   detector *det = new detector();
   if(debug>0)printf("instantiating collision\n");
   collision *vtx = new collision();
   double     ptc[4];
   TClonesArray * hitsL1 = new TClonesArray("hit",50);
   TClonesArray * hitsL2 = new TClonesArray("hit",50);
   TClonesArray * particles = new TClonesArray("particle",100);
   hitTree->Branch("Vertex",&ptc[0],"x0/D:y0:z0:m");
   hitTree->Branch("L1hit",&hitsL1);
   hitTree->Branch("L2hit",&hitsL2);

   for (int i = 0; i < N_EVENTS; i++) {
      if(debug >3)printf("\ngenerating Collision %d\n",i);
      vtx->generateVertex(ptc);
      vtx->generateParticles(ptc,particles);
      if(debug>2)printf("\rx: %f\t y: %f\t z: %f\t m: %f",ptc[0],ptc[1],ptc[2],ptc[3]);
      //beampipe scattering
      if(multScat) det->multipleScattering(particles, 0.001);
      if(debug >3)printf("intersecting layer 1\n");
      det->intersect(ptc, particles, hitsL1, 270., 40.,0.12,0.03,randomNoise);
      if(multScat) det->multipleScattering(particles, 0.001);
      if(debug >3)printf("intersecting layer 2\n");
      det->intersect(ptc, particles, hitsL2, 270., 70.,0.12,0.03,randomNoise);
      if(multScat) det->multipleScattering(particles, 0.001);
      hitTree->Fill();
      hitsL1->Clear();
      hitsL2->Clear();
      printf("\rSimulation Progress: %u%%",(uint8_t)(i*100/N_EVENTS+1));
   }
   if(debug>0)printf("\nmontecarlo finished\n");
   simulation->Write("", TObject::kOverwrite);
   simulation->Close();
   delete vtx;
   delete det;
   if(debug>0)printf("simulation run succesful, cpu time: %.2f s\n", watch.CpuTime());
   return 0;
}
