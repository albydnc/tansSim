// @(#)root/simulationclass
// Author: Alberto Perro 15/11/19

#include <Riostream.h>
#include "tracklet.h"

ClassImp(tracklet)

   tracklet::tracklet()
   : TObject()
{
}
tracklet::tracklet(const hit &hitL1, const hit &hitL2, const double l1rad, const double l2rad)
   : TObject(), _L1(hitL1), _L2(hitL2), _l1R(l1rad), _l2R(l2rad)
{
}
tracklet::tracklet(const hit *const phitL1, const hit *const phitL2, const double l1rad, const double l2rad)
   : TObject(), _L1(*phitL1), _L2(*phitL2), _l1R(l1rad), _l2R(l2rad)
{
}
tracklet::~tracklet() {}

tracklet &tracklet::operator=(const tracklet &source)
{
   this->~tracklet();
   new (this) tracklet(source);
   return *this;
}

double tracklet::findZVertex()
{
   // r = a*z + b
   double a = (_l2R - _l1R) / (_L2.getZ() - _L1.getZ()); // a=(r2-r1)/(z2-z1)
   double b = _l2R - a * _L2.getZ();                     // b = r2 - a*z2
   return -b / a;                                        // r==0 --> z = -b/a
}
