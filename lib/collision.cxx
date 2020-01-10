// @(#)root/simulationclass
// Author: Alberto Perro 15/11/19

#include <Riostream.h>
#include "collision.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TString.h"
#include "hit.h"
#include "particle.h"

ClassImp(collision)

   collision::collision()
   : TObject()
{
   _kinemFile = new TFile("./data/kinem.root");
   _distMult  = (TH1F *)_kinemFile->Get("hmul");
   _distEta   = (TH1F *)_kinemFile->Get("heta");
}

collision::collision(TFile *sourceFile) : TObject(), _kinemFile(sourceFile)
{
   _distMult = (TH1F *)_kinemFile->Get("hmul");
   _distEta  = (TH1F *)_kinemFile->Get("heta");
}

collision::collision(const collision &source) : TObject(source), _kinemFile(source._kinemFile) {}

collision::~collision()
{
   delete _distEta;
   delete _distMult;
   _kinemFile->Close();
   delete _kinemFile;
}

collision &collision::operator=(const collision &source)
{
   this->~collision();
   new (this) collision(source);
   return *this;
}

void collision::generateCollision(const uint8_t distType)
{
   // generate coordinates
   _ptc[0] = gRandom->Gaus(_mx, _sx);
   _ptc[1] = gRandom->Gaus(_my, _sy);
   _ptc[2] = gRandom->Gaus(_mz, _sz);
   uint16_t      mult  = 0;
   switch (distType) {
   case KINEMATIC: mult = _distMult->GetRandom(); break;
   case FIXED: mult = 25; break;
   case UNIFORM: mult = (uint16_t)(1 + 99 * gRandom->Rndm()); break;
   case GAUSSIAN: mult = gRandom->Gaus(20, 10); break;
   }
   _ptc[3] = mult;
}

// getters
void collision::getCoordinates(double *ptc)
{
   ptc[0] = _ptc[0];
   ptc[1] = _ptc[1];
   ptc[2] = _ptc[2];
   ptc[3] = _ptc[3];
}

// setters
void collision::ImportKinem(TString path)
{
   if (!_kinemFile->IsOpen()) _kinemFile = new TFile(path);
   _distMult = (TH1F *)_kinemFile->Get("hmul");
   _distEta  = (TH1F *)_kinemFile->Get("heta");
}

void collision::getDir(particle *dir)
{
   double phi = 2. * TMath::Pi() * gRandom->Rndm();
   double theta   = 2. * TMath::ATan(TMath::Exp(-(_distEta->GetRandom())));
   dir->setTheta(theta);
   dir->setPhi(phi);
}
