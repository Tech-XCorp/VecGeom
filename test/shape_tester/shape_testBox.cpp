#include "../benchmark/ArgParser.h"
#include "ShapeTester.h"
#include "volumes/PlacedVolume.h"
#include "volumes/Box.h"

using VPlacedVolume = vecgeom::VPlacedVolume;
using VGBox         = vecgeom::SimpleBox;

template <typename ImplT>
int runTester(ImplT const *shape, int npoints, bool debug, bool stat);

int main(int argc, char *argv[])
{
  OPTION_INT(npoints, 10000);
  OPTION_BOOL(debug, false);
  OPTION_BOOL(stat, false);

  OPTION_DOUBLE(dx, 10.);
  OPTION_DOUBLE(dy, 15.);
  OPTION_DOUBLE(dz, 20.);

  auto box = new VGBox("vecgeomBox", dx, dy, dz);
  box->Print();
  return runTester<VPlacedVolume>(box, npoints, debug, stat);
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
