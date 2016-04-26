#include "volumes/LogicalVolume.h"
#include "volumes/Box.h"
#include "volumes/Parallelepiped.h"
#include "benchmarking/Benchmarker.h"
#include "management/GeoManager.h"
#include "ArgParser.h"

using namespace vecgeom;

int main(int argc, char *argv[]) {
  OPTION_INT(npoints, 16384);
  OPTION_INT(nrep, 1024);
  OPTION_DOUBLE(dx, 20.);
  OPTION_DOUBLE(dy, 30.);
  OPTION_DOUBLE(dz, 40.);
  OPTION_DOUBLE(alpha, 30. / 180. * kPi);
  OPTION_DOUBLE(theta, 15. / 180. * kPi);
  OPTION_DOUBLE(phi, 30. / 180. * kPi);

  UnplacedBox worldUnplaced = UnplacedBox(dx * 4, dy * 4, dz * 4);
  UnplacedParallelepiped paraUnplaced(dx, dy, dz, alpha, theta, phi);
  LogicalVolume world("w0rld", &worldUnplaced);
  LogicalVolume para("p4r4", &paraUnplaced);
  Transformation3D placement = Transformation3D(5, 5, 5);
  world.PlaceDaughter(&para, &placement);

  VPlacedVolume *worldPlaced = world.Place();

  GeoManager::Instance().SetWorldAndClose(worldPlaced);

  Benchmarker tester(GeoManager::Instance().GetWorld());
  tester.SetVerbosity(3);
  tester.SetRepetitions(nrep);
  tester.SetPointCount(npoints);
  tester.RunBenchmark();

  return 0;
}
