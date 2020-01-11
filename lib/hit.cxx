// @(#)root/simulationclass
// Author: Alberto Perro 15/11/19

#include <Riostream.h>
#include "hit.h"
#include "TMath.h"

ClassImp(hit)

   hit::hit()
   : TObject()
{
}
hit::hit(double theta, double z) : TObject(), _theta(theta), _z(z) {}

hit::hit(const hit &source) : TObject(source), _theta(source._theta), _z(source._z) {}

hit::~hit() {}

hit &hit::operator=(const hit &source)
{
   this->~hit();
   new (this) hit(source);
   return *this;
}

double hit::getTheta()
{
   return _theta;
}

double hit::getZ()
{
   return _z;
}

double hit::deltaTheta(hit &hitL1, hit &hitL2)
{
   return deltaTheta(&hitL1, &hitL2);
}

double hit::deltaTheta(hit *hitL1, hit *hitL2)
{
   double th1 = hitL1->getTheta();
   double th2 = hitL2->getTheta();

   double result = th1 - th2;
   while (result > TMath::Pi()) result -= 2 * TMath::Pi();
   while (result <= -TMath::Pi()) result += 2 * TMath::Pi();
   return fabs(result);
}
