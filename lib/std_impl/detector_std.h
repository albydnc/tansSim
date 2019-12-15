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
   detector(const detector &source);            // copia oggetto già dichiarato SENZA modificarlo
   virtual ~detector();                         // distruttore virtuale perché ereditato da TObject
   detector &operator=(const detector &source); // assignment operator
   void intersect(double *ptc, std::vector<std::array<double, 2>> &particles, std::vector<std::array<double, 2>> &hits,
                  double width, double radius);
   void multipleScattering(std::vector<std::array<double, 2>> &particles, double rmsTh);

   ClassDef(detector, 1) private:
};

#endif
