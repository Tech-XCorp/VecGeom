/// \file OrbBenchmark.cpp
/// \author Raman Sehgal (raman.sehgal@cern.ch)

#include "volumes/LogicalVolume.h"
#include "volumes/Box.h"
#include "volumes/Orb.h"
#include "benchmarking/Benchmarker.h"
#include "management/GeoManager.h"
#include "ArgParser.h"

using namespace vecgeom;

int main(int argc, char* argv[]) {
  OPTION_INT(npoints);
  OPTION_INT(nrep);
  OPTION_DOUBLE(r);

  UnplacedBox worldUnplaced = UnplacedBox(r*4, r*4, r*4);
  UnplacedOrb orbUnplaced = UnplacedOrb(r);
  LogicalVolume world = LogicalVolume("w0rld", &worldUnplaced);
  LogicalVolume orb = LogicalVolume("p4r4", &orbUnplaced);
  Transformation3D placement = Transformation3D(5, 5, 5);
  world.PlaceDaughter(&orb, &placement);
  //world.PlaceDaughter(&orb, &Transformation3D::kIdentity);

  VPlacedVolume *worldPlaced = world.Place();

  GeoManager::Instance().set_world(worldPlaced);

  Benchmarker tester(GeoManager::Instance().world());
  tester.SetVerbosity(3);
  tester.SetPoolMultiplier(1);
  tester.SetPointCount(npoints);
  tester.SetRepetitions(nrep);
  //tester.RunInsideBenchmark();
  //tester.RunToOutBenchmark();
  tester.RunBenchmark();

  return 0;
}
