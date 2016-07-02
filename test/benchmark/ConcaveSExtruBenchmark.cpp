#include "volumes/LogicalVolume.h"
#include "volumes/SExtru.h"
#include "volumes/Box.h"
#include "benchmarking/Benchmarker.h"
#include "management/GeoManager.h"
#include "ArgParser.h"
#include "base/Stopwatch.h"
#include <iostream>

using namespace vecgeom;

int main(int argc, char *argv[])
{
  OPTION_INT(npoints, 1024);
  OPTION_INT(nrep, 4);
  OPTION_INT(N, 12);

  double dx = 5;
  double dy = 5;
  double dz = 5;
  UnplacedBox worldUnplaced(dx * 4, dy * 4, dz * 4);

  double x[N], y[N];
  for (size_t i = 0; i < (size_t)N; ++i) {
    x[i] = dx * std::sin(i * (2. * M_PI) / N);
    y[i] = dy * std::cos(i * (2. * M_PI) / N);
  }
  // make it concave by overwriting on of the points
  // to be at origin
  x[0] = 0.;
  y[1] = 0.;

  UnplacedSExtruVolume extru(N, x, y, -dz, dz);

  LogicalVolume world("world", &worldUnplaced);
  LogicalVolume extrul("extrul", &extru);

  Transformation3D placement(0., 0, 0);
  world.PlaceDaughter("extrul", &extrul, &placement);
  VPlacedVolume *worldPlaced = world.Place();
  GeoManager::Instance().SetWorldAndClose(worldPlaced);

  Benchmarker tester(GeoManager::Instance().GetWorld());
  tester.SetVerbosity(3);
  tester.SetRepetitions(nrep);
  tester.SetPointCount(npoints);
  tester.SetPoolMultiplier(1);
  tester.RunBenchmark();
  tester.RunToOutFromBoundaryBenchmark();
  tester.RunToOutFromBoundaryExitingBenchmark();
  tester.RunToInFromBoundaryBenchmark();
}
