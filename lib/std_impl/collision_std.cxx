// @(#)root/simulationclass
// Author: Alberto Perro 15/11/19
#include <Riostream.h>
#include "collision_std.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TString.h"
#include <vector>
#include <array>

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

void collision::generateVertex()
{
   // generate coordinates
   _ptc[0] = gRandom->Gaus(_mx, _sx);
   _ptc[1] = gRandom->Gaus(_my, _sy);
   _ptc[2] = gRandom->Gaus(_mz, _sz);
}

void collision::generateVertex(double *ptc)
{
   // generate coordinates
   _ptc[0] = gRandom->Gaus(_mx, _sx);
   _ptc[1] = gRandom->Gaus(_my, _sy);
   _ptc[2] = gRandom->Gaus(_mz, _sz);
   ptc[0]  = _ptc[0];
   ptc[1]  = _ptc[1];
   ptc[2]  = _ptc[2];
}

// getters
void collision::getCoordinates(double *ptc)
{
   ptc[0] = _ptc[0];
   ptc[1] = _ptc[1];
   ptc[2] = _ptc[2];
}

void collision::generateParticles(std::vector<std::array<double, 2>> &particles)
{
   // generate multiplicity
   uint16_t mult = _distMult->GetRandom();
   for (int i = 0; i < mult; i++) {
      std::array<double, 2> part;
      getDir(part);
      particles.push_back(part);
   }
}

// setters
void collision::ImportKinem(TString path)
{
   if (!_kinemFile->IsOpen()) _kinemFile = new TFile(path);
   _distMult = (TH1F *)_kinemFile->Get("hmul");
   _distEta  = (TH1F *)_kinemFile->Get("heta");
}

void collision::getDir(std::array<double, 2> &dir)
{
   double theta = 2. * TMath::Pi() * gRandom->Rndm();
   double phi   = 2. * TMath::ATan(TMath::Exp(-(_distEta->GetRandom())));
   dir          = {theta, phi};

}
