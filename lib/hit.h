// @(#)root/simulationclass
// Author: Alberto Perro 15/11/19

#ifndef hit_H
#define hit_H

#include "TObject.h"

class hit : public TObject {
public:
   hit();
   hit(double theta, double z);
   hit(const hit &source);            // copia oggetto già dichiarato SENZA modificarlo
   virtual ~hit();                          // distruttore virtuale perché ereditato da TObject
   hit &operator=(const hit &source); // assignment operator
   double getTheta();
   double getZ();
   double deltaTheta( hit * hitL1, hit * hitL2);
   double deltaTheta( hit &hitL1, hit &hitL2);
   ClassDef(hit,1)
private:
  double _theta, _z;
};
#endif
