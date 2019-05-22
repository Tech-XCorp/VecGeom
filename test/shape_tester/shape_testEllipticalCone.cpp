#include "../benchmark/ArgParser.h"
#include "ShapeTester.h"
#include "volumes/PlacedVolume.h"
#include "volumes/EllipticalCone.h"

using VPlacedVolume = vecgeom::VPlacedVolume;
using VGCone        = vecgeom::SimpleEllipticalCone;

template <typename ImplT>
int runTester(ImplT const *shape, int npoints, bool debug, bool stat);

int main(int argc, char *argv[])
{
  OPTION_INT(npoints, 10000);
  OPTION_BOOL(debug, false);
  OPTION_BOOL(stat, false);

  OPTION_DOUBLE(a, 0.5);
  OPTION_DOUBLE(b, 0.4);
  OPTION_DOUBLE(h, 10.);
  OPTION_DOUBLE(zcut, 5.);

  auto cone = new VGCone("vecgeomEllipticalCone", a, b, h, zcut);
  cone->Print();
  return runTester<VPlacedVolume>(cone, npoints, debug, stat);
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