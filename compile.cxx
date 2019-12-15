// @(#)root/simulation
// Author: Alberto Perro 15/12/19

#include <Riostream.h>
#include <sys/stat.h>
#include "TSystem.h"

void compile_simulation(){

	gSystem->CleanCompiledMacros();
  gSystem->SetBuildDir("bin", kTRUE);
  gSystem->CompileMacro("lib/hit.cxx","kg");
  gSystem->CompileMacro("lib/particle.cxx","kg");
  gSystem->CompileMacro("lib/collision.cxx","kg");
  gSystem->CompileMacro("lib/detector.cxx","kg");
  gSystem->CompileMacro("lib/tracklet.cxx","kg");
  gSystem->CompileMacro("src/simulation.C","kg");
  gSystem->CompileMacro("src/reconstruction.C","kg");
  printf("done!\n");
}
