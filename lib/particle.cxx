// @(#)root/simulationclass
// Author: Alberto Perro 15/11/19

#include <Riostream.h>
#include "particle.h"

ClassImp(particle)

   particle::particle()
   : TObject()
{
}
particle::particle(double theta, double phi) : TObject(), _theta(theta), _phi(phi) {}

particle::particle(const particle &source) : TObject(source), _theta(source._theta), _phi(source._phi) {}

particle::~particle() {}

particle &particle::operator=(const particle &source)
{
   this->~particle();
   new (this) particle(source);
   return *this;
}

double particle::getTheta()
{
   return _theta;
}

double particle::getPhi()
{
   return _phi;
}
void particle::setTheta(double theta)
{
   _theta = theta;
}

void particle::setPhi(double phi)
{
   _phi = phi;
}
