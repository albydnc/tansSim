// @(#)root/simulationclass
// Author: Alberto Perro 15/11/19

#ifndef collision_H
#define collision_H

#include "TObject.h"
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include <vector>
#include <array>

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
   void generateParticles(std::vector<std::array<double, 2>> &particles);
   void generateCollision(std::vector<std::array<double, 2>> &particles);
   // setters
   void ImportKinem(TString path);
   ClassDef(collision, 1)

private :
      // point coordinates
   double _ptc[3]; //[3]
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
   void   getDir(std::array<double, 2> &dir);
};

#endif
