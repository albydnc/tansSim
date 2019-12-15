#ifndef particle_H
#define particle_H

#include "TObject.h"

class particle : public TObject {
public:
   particle();
   particle(double theta, double phi);
   particle(const particle &source);            // copia oggetto già dichiarato SENZA modificarlo
   virtual ~particle();                          // distruttore virtuale perché ereditato da TObject
   particle &operator=(const particle &source); // assignment operator
   double getTheta();
   double getPhi();
   void setTheta(double theta);
   void setPhi(double phi);
   ClassDef(particle,1)
private:
  double _theta, _phi;
};
#endif
