/// @file ParaboloidBenchmark.cpp
/// @author Marilena Bandieramonte (marilena.bandieramonte@cern.ch)

#include "volumes/LogicalVolume.h"
#include "volumes/Box.h"
#include "volumes/Paraboloid.h"
#include "benchmarking/Benchmarker.h"
#include "management/GeoManager.h"
#include <iostream>

using namespace vecgeom;

int main()
{

  std::cout << "Paraboloid Benchmark\n";
  UnplacedBox worldUnplaced             = UnplacedBox(10., 10., 10.);
  UnplacedParaboloid paraboloidUnplaced = UnplacedParaboloid(3., 5., 7.); // rlo=3. - rhi=5. dz=7
  std::cout << "Paraboloid created\n";
  LogicalVolume world("MBworld", &worldUnplaced);
  LogicalVolume paraboloid("paraboloid", &paraboloidUnplaced);
  world.PlaceDaughter(&paraboloid, &Transformation3D::kIdentity);
  VPlacedVolume *worldPlaced = world.Place();
  GeoManager::Instance().SetWorldAndClose(worldPlaced);
  std::cout << "World set\n";

  Benchmarker tester(GeoManager::Instance().GetWorld());
  tester.SetPoolMultiplier(1);
  tester.SetVerbosity(3);
  tester.SetPointCount(1 << 10);
  tester.SetRepetitions(4);
  return tester.RunBenchmark();
}
