// @(#)root/simulation
// Author: Alberto Perro 15/11/19
#include "lib/std_impl/detector_std.h"
#include "lib/std_impl/collision_std.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include <Riostream.h>
#include <vector>
#include <array>

#ifdef __MAKECINT__
#pragma link C++ class std::array<double,2> +;
#pragma link C++ class std::vector<std::array<double, 2> > +;
#endif

#define N_EVENTS 1e6
int simulation_std()
{
  int debug = 1;
  TStopwatch watch;
  TFile * simulation = new TFile("simulation_std.root","RECREATE");
  TTree * hitTree =  new TTree("hits", "Hits tree");
   if(debug>0)printf("instantiating detector\n");
   detector *det = new detector();
   if(debug>0)printf("instantiating collision\n");
   collision *vtx = new collision();
   double     ptc[4];
   bool multScat = false;
   std::vector<std::array<double,2>> particles;
   std::vector<std::array<double,2>> hitsL1;
   std::vector<std::array<double,2>> hitsL2;
   hitTree->Branch("Vertex",&ptc[0],"x0/D:y0:z0:m/I");
   hitTree->Branch("L1hit",&hitsL1);
   hitTree->Branch("L2hit",&hitsL2);

   for (int i = 0; i < N_EVENTS; i++) {
      if(debug >3)printf("\ngenerating Collision %d\n",i);
      vtx->generateVertex(ptc);
      vtx->generateParticles(particles);
      //beampipe scattering
      if(multScat) det->multipleScattering(particles, 0.001);
      if(debug >3)printf("intersecting layer 1\n");
      det->intersect(ptc, particles, hitsL1, 270., 40.);
      if(multScat) det->multipleScattering(particles, 0.001);
      if(debug >3)printf("intersecting layer 2\n");
      det->intersect(ptc, particles, hitsL2, 270., 70.);
      if(multScat) det->multipleScattering(particles, 0.001);
      hitTree->Fill();
      particles.clear();
      hitsL1.clear();
      hitsL2.clear();
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
