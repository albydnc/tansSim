// @(#)root/simulationclass
// Author: Alberto Perro 15/11/19

#ifndef tracklet_H
#define tracklet_H

#include "TObject.h"
#include "hit.h"

class tracklet : public TObject {
public:
   tracklet();
   tracklet(const hit& hitL1, const hit& hitL2);
	 tracklet(const hit* const phitL1, const hit* const phitL2);
   virtual ~tracklet();                          // distruttore virtuale perché ereditato da TObject
   tracklet &operator=(const tracklet &source); // assignment operator
   double findZVertex();
   ClassDef(tracklet,1)
private:
  hit _L1,_L2;
  double _l1R = 40.0, _l2R = 70.0;
};
#endif
