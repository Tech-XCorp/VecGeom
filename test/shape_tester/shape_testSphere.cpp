#include "../benchmark/ArgParser.h"
#include "ShapeTester.h"
#include "VecGeom/volumes/PlacedVolume.h"
#include "VecGeom/volumes/Sphere.h"

using VPlacedVolume = vecgeom::VPlacedVolume;
using VGSphere      = vecgeom::PlacedSphere;

template <typename ImplT>
int runTester(ImplT const *shape, int npoints, bool debug, bool stat);

int main(int argc, char *argv[])
{
  OPTION_INT(npoints, 10000);
  OPTION_BOOL(debug, false);
  OPTION_BOOL(stat, false);

  OPTION_DOUBLE(rmin, 15.);
  OPTION_DOUBLE(rmax, 20.);
  OPTION_DOUBLE(sphi, 0.);
  OPTION_DOUBLE(dphi, vecgeom::kTwoPi / 3.);
  OPTION_DOUBLE(stheta, 0.);
  OPTION_DOUBLE(dtheta, vecgeom::kTwoPi);

  auto sphere = new VGSphere("vecgeomSphere", rmin, rmax, sphi, dphi, stheta, dtheta);
  sphere->Print();
  return runTester<VPlacedVolume>(sphere, npoints, debug, stat);
}

template <typename ImplT>
int runTester(ImplT const *shape, int npoints, bool debug, bool stat)
{

  ShapeTester<ImplT> tester;
  tester.setDebug(debug);
  tester.setStat(stat);
  tester.SetMaxPoints(npoints);
  int errcode = tester.Run(shape);

  std::cout << "Final Error count for Shape *** " << shape->GetName() << "*** = " << errcode << "\n";
  std::cout << "=========================================================\n";

  if (shape) delete shape;
  return errcode;
}
