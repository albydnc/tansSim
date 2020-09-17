// @(#)root/simulationclass
// Author: Alberto Perro 15/11/19

#include <Riostream.h>
#include "detector.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TString.h"
#include "TClonesArray.h"
#include "hit.h"
#include "particle.h"

ClassImp(detector)

   detector::detector()
   : TObject()
{
}

hit *detector::intersect(double *ptc, particle *part, double *new_pt, double width, double radius, double z_sm,
                         double rho_sm, bool multScat, double rmsTh)
{
   uint16_t success = 0;
   double   theta   = part->getTheta(); // beam axis angle
   double   phi     = part->getPhi();   // anomaly
   if (theta > TMath::Pi()) theta = theta - TMath::Pi();
   double c1    = TMath::Sin(theta) * TMath::Cos(phi);
   double c2    = TMath::Sin(theta) * TMath::Sin(phi);
   double c3    = TMath::Cos(theta);
   double a     = ptc[0] * c1 + ptc[1] * c2;
   double b     = c1 * c1 + c2 * c2;
   double c     = ptc[0] * ptc[0] + ptc[1] * ptc[1] - radius * radius;
   double delta = a * a - (b * c);
   double t;
   double t_temp = -((ptc[0] * c1 + ptc[1] * c2) + TMath::Sqrt(delta)) / b;
   if (t_temp > 0) {
      t = t_temp;
   } else {
      t = -((ptc[0] * c1 + ptc[1] * c2) - TMath::Sqrt(delta)) / b;
   }
   double z = ptc[2] + c3 * t;
   double rho;
   if (z <= width / 2. && z >= -width / 2) {
      double x = ptc[0] + c1 * t;
      double y = ptc[1] + c2 * t;

      if (x > 0 && y > 0) {
         rho = TMath::ATan(y / x); //ok
      } else if (x == 0) {
         rho = (y > 0 ? 1. : -1.) * TMath::Pi() / 2;
      } else if (x > 0 && y < 0) {
         rho = 2 * TMath::Pi() + TMath::ATan(y / x);
      } else if (x < 0 && y > 0) {
         rho = TMath::Pi() + TMath::ATan(y / x);
      } else if (x < 0 && y < 0) {
         rho = TMath::Pi() + TMath::ATan(y / x);
      }

      rho_sm /= radius; // normalize smearing
      double z_smear   = gRandom->Gaus(z, z_sm);     // smearing z
      double rho_smear = gRandom->Gaus(rho, rho_sm); // smearing polar angle
      if (rho_smear < 0) {
         rho_smear = rho_smear + 2 * TMath::Pi();
      } else if (rho_smear > 2 * TMath::Pi()) {
         rho_smear = rho_smear - 2 * TMath::Pi();
      }
      if (multScat) multipleScattering(part, rmsTh);
      new_pt[0] = x;
      new_pt[1] = y;
      new_pt[2] = z;
      return new hit(rho_smear, z_smear);
   }
   return NULL;
}

hit *detector::randomHits(int width)
{
   double z   = width - (gRandom->Rndm() * width / 2);
   double rho = 2 * TMath::Pi() * gRandom->Rndm();
   return new hit(rho, z);
}

void detector::multipleScattering(particle *part, double rmsTh)
{
   double theta = part->getTheta();
   double phi   = part->getPhi();

   double dTh  = gRandom->Gaus(0, rmsTh); // angoli di cui viene ruotata la direzione entrante
   double dPhi = 2 * TMath::Pi() * gRandom->Rndm();
   if (theta > TMath::Pi()) theta = theta - TMath::Pi();
   double rotMatx[3][3];
   rotMatx[0][0] = -1*TMath::Sin(phi);
   rotMatx[1][0] = TMath::Cos(phi);
   rotMatx[2][0] = 0.;
   rotMatx[0][1] = -1* TMath::Cos(theta) * TMath::Cos(phi);
   rotMatx[1][1] = -1* TMath::Cos(theta) * TMath::Sin(phi);
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

   part->setTheta(TMath::ACos(res[2]));
   part->setPhi(TMath::ATan2(res[1], res[0]));
}
