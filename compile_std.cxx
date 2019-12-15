#include <Riostream.h>
#include <sys/stat.h>
#include "TSystem.h"

void compile_simulation_std(){

	gSystem->CleanCompiledMacros();
  gSystem->SetBuildDir("bin_std", kTRUE);
  gSystem->CompileMacro("lib/std_impl/detector_std.cxx","kg");
  gSystem->CompileMacro("lib/std_impl/collision_std.cxx","kg");
  gSystem->CompileMacro("src/simulation_std.C","kg");
  //gSystem->CompileMacro("src/reconstruction_std.C","kg");
  printf("done!\n");
}
