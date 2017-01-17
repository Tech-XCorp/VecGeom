#include "../benchmark/ArgParser.h"
#include "ShapeTester.h"
#include "volumes/Hype.h"
typedef vecgeom::SimpleHype Hype_t;

int main(int argc, char *argv[])
{
  OPTION_INT(npoints, 10000);
  OPTION_BOOL(debug, false);
  OPTION_BOOL(stat, false);
  OPTION_BOOL(usolids, false);

  if (usolids) {
    std::cerr << "\n*** ERROR: '-usolids true' is not valid for SExtru shape!\n Aborting...\n\n";
    return 1;
  }

  using vecgeom::kPi;
  auto hype = new Hype_t("test_VecGeomHype", 5., 20, kPi / 6, kPi / 3, 50);
  hype->Print();

  ShapeTester<vecgeom::VPlacedVolume> tester;
  tester.setConventionsMode(usolids);
  tester.setDebug(debug);
  tester.setStat(stat);
  tester.SetMaxPoints(npoints);
  tester.SetSolidTolerance(1.e-9);
  tester.SetTestBoundaryErrors(true);
  int errCode = tester.Run(hype);

  std::cout << "Final Error count for Shape *** " << hype->GetName() << "*** = " << errCode << " ("
            << (tester.getConventionsMode() ? "USolids" : "VecGeom") << " conventions)\n";
  std::cout << "=========================================================" << std::endl;

  if (hype) delete hype;
  return 0;
}
