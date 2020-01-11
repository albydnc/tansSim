// @(#)root/simulationclass
// Author: Alberto Perro 15/11/19

#ifndef DETECTOR_H
#define DETECTOR_H

#include "TObject.h"
#include "TFile.h"
#include "TString.h"
#include "TClonesArray.h"
#include "hit.h"
#include "particle.h"

class detector : public TObject {
public:
   detector();
   // intersecting layer of detector given parameters
   hit *intersect(double *ptc, particle *part, double *new_pt, double width, double radius, double z_sm, double rho_sm,
                  bool multScat, double rmsTh); // multiple scattering generator
   hit *randomHits(int width);
   void multipleScattering(particle *part, double rmsTh);

   ClassDef(detector, 0)
};

#endif
