// @(#)root/simulationclass
// Author: Alberto Perro 15/11/19

#ifndef collision_H
#define collision_H

#include "TObject.h"
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "particle.h"

class collision : public TObject {
public:
   collision();
   collision(TFile *sourceFile);
   collision(const collision &source);            // copia oggetto già dichiarato SENZA modificarlo
   virtual ~collision();                          // distruttore virtuale perché ereditato da TObject
   collision &operator=(const collision &source); // assignment operator
   void       generateVertex();
   void       generateVertex(double *ptc);
   // getters
   void getCoordinates(double *ptc);
   void generateParticles(double * ptc, TClonesArray * particles);
   // setters
   void ImportKinem(TString path);
   ClassDef(collision, 1)

private :
      // point coordinates
   double _ptc[4]; //[4]
   // distrib params
   double _mx = 0.;
   double _my = 0.;
   double _mz = 0.;
   double _sx = 0.1;
   double _sy = 0.1;
   double _sz = 53.0;
   // kinematic file
   TH1F * _distMult;
   TH1F * _distEta;
   TFile *_kinemFile;
   void   getDir(particle &dir);
};

#endif
