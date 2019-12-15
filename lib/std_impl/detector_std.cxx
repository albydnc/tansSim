// @(#)root/simulationclass
// Author: Alberto Perro 15/11/19
#include <Riostream.h>
#include "detector_std.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TString.h"
#include "TClonesArray.h"
#include <vector>
#include <array>

ClassImp(detector)

   detector::detector()
   : TObject()
{
}

detector::detector(const detector &source) : TObject(source) {}

detector::~detector() {}

detector &detector::operator=(const detector &source)
{
   this->~detector();
   new (this) detector(source);
   return *this;
}

void detector::intersect(double * ptc, std::vector<std::array<double, 2>> &particles,
                         std::vector<std::array<double, 2>> &hits, double width, double radius)
{
  uint16_t success = 0;
   for (uint32_t i = 0; i < particles.size(); i++) {
      double theta = particles.at(i).at(0);
      double phi   = particles.at(i).at(1);
      if (theta > TMath::Pi()) theta = theta - TMath::Pi();
      double c1    = TMath::Sin(theta) * TMath::Cos(phi);
      double c2    = TMath::Sin(theta) * TMath::Sin(phi);
      double c3    = TMath::Cos(theta);
      double a     = ptc[0] * c1 + ptc[1] * c2;
      double b     = c1 * c1 + c2 * c2;
      double c     = ptc[0] * ptc[0] + ptc[1] * ptc[1] - radius * radius;
      double delta = a * a - b * c;
      double t;
      double t_temp = -((ptc[0] * c1 + ptc[1] * c2) + TMath::Sqrt(delta)) / b;
      if (t_temp > 0)
         t = t_temp;
      else {
         t = -((ptc[0] * c1 + ptc[1] * c2) - TMath::Sqrt(delta)) / b;
      }
      double z = ptc[2] + c3 * t;
      double rho;
      if (z <= width / 2. && z >= -width / 2) {
         double x = ptc[0] + c1 * t;
         double y = ptc[1] + c2 * t;
         if (x > 0) {
            rho = TMath::ATan(y / x);
         } else if (x == 0) {
            rho = (y > 0 ? 1. : -1.) * TMath::Pi() / 2;
         } else {
            rho = (y > 0 ? 1. : -1.) * TMath::Pi() + TMath::ATan(y / x);
         }
         std::array<double,2> hit = {rho,z};
         hits.push_back(hit);
         success++;
      }else{
        particles.erase(particles.begin()+i);
      }
   }
}

void detector::multipleScattering(std::vector<std::array<double, 2>> &particles, double rmsTh)
{
   for (uint32_t k = 0; k < particles.size(); k++) {
      double theta = particles.at(k).at(0);
      double phi   = particles.at(k).at(1);

      double dTh  = gRandom->Gaus(0, rmsTh); // angoli di cui viene ruotata la direzione entrante
      double dPhi = 2 * TMath::Pi() * gRandom->Rndm();
      if (theta > TMath::Pi()) theta = theta - TMath::Pi();
      double rotMatx[3][3];
      rotMatx[0][0] = TMath::Sin(phi);
      rotMatx[1][0] = TMath::Cos(phi);
      rotMatx[2][0] = 0.;
      rotMatx[0][1] = TMath::Cos(theta) * TMath::Cos(phi);
      rotMatx[1][1] = TMath::Cos(theta) * TMath::Sin(phi);
      rotMatx[2][1] = TMath::Sin(theta);
      rotMatx[0][2] = TMath::Sin(theta) * TMath::Cos(phi);
      rotMatx[1][2] = TMath::Sin(theta) * TMath::Sin(phi);
      rotMatx[2][2] = TMath::Cos(theta);

      double rotVec[3];
      rotVec[0] = TMath::Sin(dTh) * TMath::Cos(dPhi);
      rotVec[1] = TMath::Sin(dTh) * TMath::Sin(dPhi);
      rotVec[2] = TMath::Cos(dTh);

      double res[3];
      for (int i = 0; i < 3; i++) {
         res[i] = 0.;
         for (int j = 0; j < 3; j++) {
            res[i] += rotMatx[i][j] * rotVec[j];
         }
      }

      particles.at(k).at(0) = TMath::ACos(res[2]);
      particles.at(k).at(1) = TMath::ATan2(res[1], res[2]);
   }
}
