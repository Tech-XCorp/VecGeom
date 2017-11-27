#include "../benchmark/ArgParser.h"
#include "ShapeTester.h"
#include "volumes/PlacedVolume.h"

#ifdef VECGEOM_USOLIDS
#include "UCons.hh"
#endif
#include "volumes/Cone.h"

using VPlacedVolume = vecgeom::VPlacedVolume;
using VGCone        = vecgeom::SimpleCone;

template <typename ImplT>
int runTester(ImplT const *shape, int npoints, bool usolids, bool debug, bool stat);

int main(int argc, char *argv[])
{
  OPTION_INT(npoints, 10000);
  OPTION_BOOL(debug, false);
  OPTION_BOOL(stat, false);
  OPTION_BOOL(usolids, false);

  OPTION_DOUBLE(rmin1, 4.);
  OPTION_DOUBLE(rmax1, 6.);
  OPTION_DOUBLE(rmin2, 5.);
  OPTION_DOUBLE(rmax2, 7.);
  OPTION_DOUBLE(dz, 2.);
  OPTION_DOUBLE(sphi, 0.);
  OPTION_DOUBLE(dphi, vecgeom::kTwoPi);

  if (usolids) {
#ifndef VECGEOM_USOLIDS
    std::cerr << "\n*** ERROR: library built with -DUSOLIDS=OFF and user selected '-usolids true'!\n Aborting...\n\n";
    return 1;
#else
    // USolids conventions
    auto cone = new UCons("usolidsCone", rmin1, rmax1, rmin2, rmax2, dz, sphi, dphi);
    cone->StreamInfo(std::cout);
    return runTester<VUSolid>(cone, npoints, usolids, debug, stat);
#endif
  }

  else {
    // VecGeom conventions
    auto cone = new VGCone("vecgeomCone", rmin1, rmax1, rmin2, rmax2, dz, sphi, dphi);
    cone->Print();
    return runTester<VPlacedVolume>(cone, npoints, usolids, debug, stat);
  }
}

template <typename ImplT>
int runTester(ImplT const *shape, int npoints, bool usolids, bool debug, bool stat)
{

  ShapeTester<ImplT> tester;
  tester.setConventionsMode(usolids);
  tester.SetSolidTolerance(1.e-7);
  tester.setDebug(debug);
  tester.setStat(stat);
  tester.SetMaxPoints(npoints);
  int errcode = tester.Run(shape);

  std::cout << "Final Error count for Shape *** " << shape->GetName() << "*** = " << errcode << " ("
            << (tester.getConventionsMode() ? "USolids" : "VecGeom") << " conventions)\n";
  std::cout << "=========================================================\n";
  if (shape) delete shape;
  return errcode;
}
