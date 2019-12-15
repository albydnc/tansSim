// @(#)root/simulationclass
// Author: Alberto Perro 15/11/19

#ifndef DETECTOR_H
#define DETECTOR_H

#include "TObject.h"
#include "TFile.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include <vector>
#include <array>

class detector : public TObject {
public:
   detector();
   detector(const detector &source);
   virtual ~detector();
   detector &operator=(const detector &source); // assignment operator

   // intersecting layer of detector given parameters
   void intersect(double *ptc, TClonesArray *particles, TClonesArray *Hits, double width, double radius, double z_sm,
                  double rho_sm, int randomNoise);
   // multiple scattering generator
   void multipleScattering(TClonesArray *particles, double rmsTh);

   ClassDef(detector, 0)
};

#endif
