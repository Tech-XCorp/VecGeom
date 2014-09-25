#include "volumes/LogicalVolume.h"
#include "volumes/Polyhedron.h"
#include "benchmarking/Benchmarker.h"
#include "management/GeoManager.h"

using namespace vecgeom;

int main() {

  UnplacedBox worldUnplaced = UnplacedBox(3., 3., 3.);
  constexpr int nPlanes = 2;
  Precision zPlanes[nPlanes] = {-1, 1};
  Precision rInner[nPlanes] = {0.5, 1};
  Precision rOuter[nPlanes] = {1, 2};
  UnplacedPolyhedron polyhedronUnplaced(5, nPlanes, zPlanes, rInner, rOuter);

  LogicalVolume world = LogicalVolume("world", &worldUnplaced);
  LogicalVolume polyhedron = LogicalVolume("polyhedron", &polyhedronUnplaced);

  Transformation3D placement;
  world.PlaceDaughter("polyhedron", &polyhedron, &placement);

  VPlacedVolume *worldPlaced = world.Place();

  GeoManager::Instance().set_world(worldPlaced);

  Benchmarker tester(GeoManager::Instance().world());
  tester.SetVerbosity(3);
  tester.SetPoolMultiplier(1);
  tester.SetRepetitions(1);
  tester.SetPointCount(16);
  tester.RunBenchmark();

  return 0;
}
