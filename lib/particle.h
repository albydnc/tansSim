// @(#)root/simulationclass
// Author: Alberto Perro 15/11/19

#ifndef particle_H
#define particle_H

#include "TObject.h"

class particle : public TObject {
public:
   particle();
   particle(double theta, double phi);
   particle(const particle &source);
   virtual ~particle();
   particle &operator=(const particle &source); // assignment operator
   // getters
   double getTheta();
   double getPhi();
   // setters
   void setTheta(double theta);
   void setPhi(double phi);
   ClassDef(particle, 1)
private:
   double _theta, _phi;
};
#endif
