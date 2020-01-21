//
// TestHype
//

//.. Ensure asserts are compiled in
#undef NDEBUG
#include "VecGeom/base/FpeEnable.h"

#include "VecGeom/base/Global.h"
#include "VecGeom/base/Vector3D.h"
#include "VecGeom/volumes/Box.h"
#include "VecGeom/volumes/Hype.h"
#include "ApproxEqual.h"
#include <iomanip>
#include <cassert>
#include <cmath>
#include <string>

using vecgeom::kInfLength;
using vecgeom::kPi;

template <class Hype_t, class Vec_t = vecgeom::Vector3D<vecgeom::Precision>>
bool TestHype()
{

  std::cout << std::setprecision(15);
  // int verbose=0;
  double fRmin = 10, fRmax = 20, stIn = kPi / 4, stOut = kPi / 3, halfZ = 50;
  vecgeom::Precision fR = 9.;
  Vec_t pzero(0, 0, 0);
  Vec_t pbigx(100, 0, 0), pbigy(0, 100, 0), pbigz(0, 0, 100);
  Vec_t pbigmx(-100, 0, 0), pbigmy(0, -100, 0), pbigmz(0, 0, -100);
  Vec_t ponx(fR, 0., 0.);   // point on surface on X axis
  Vec_t ponmx(-fR, 0., 0.); // point on surface on minus X axis
  Vec_t pony(0., fR, 0.);   // point on surface on Y axis
  Vec_t ponmy(0., -fR, 0.); // point on surface on minus Y axis
  Vec_t ponz(0., 0., fR);   // point on surface on Z axis
  Vec_t ponmz(0., 0., -fR); // point on surface on minus Z axis

  Vec_t ponxsideO(fRmax, 0, 0), ponysideO(0, fRmax, 0); //, ponzsideO(fRmax,0,halfZ);
  Vec_t ponxsideI(fRmin, 0, 0), ponysideI(0, fRmin, 0); //, ponzsideI(fRmin,0,halfZ);;

  Vec_t ponmxside(-fR, 0, 0), ponmyside(0, -fR, 0), ponmzside(0, 0, -fR);

  Vec_t vx(1, 0, 0), vy(0, 1, 0), vz(0, 0, 1);
  Vec_t vmx(-1, 0, 0), vmy(0, -1, 0), vmz(0, 0, -1);
  Vec_t vxy(1 / std::sqrt(2.0), 1 / std::sqrt(2.0), 0);
  Vec_t vmxy(-1 / std::sqrt(2.0), 1 / std::sqrt(2.0), 0);
  Vec_t vmxmy(-1 / std::sqrt(2.0), -1 / std::sqrt(2.0), 0);
  Vec_t vxmy(1 / std::sqrt(2.0), -1 / std::sqrt(2.0), 0);
  Vec_t vxmz(1 / std::sqrt(2.0), 0, -1 / std::sqrt(2.0));

  Hype_t b1("Solid VecGeomHype #1", fRmin, fRmax, kPi / 4, kPi / 3, 50);
  Hype_t b2("Solid VecGeomHype #2", 10, 20, kPi / 4, kPi / 4, 50);

  // Check name
  std::cout << "Name : " << b1.GetName() << std::endl;
  // assert(b1.GetName()=="Solid VecGeomHype #1");
  assert(!strcmp(b1.GetName(), "Solid VecGeomHype #1"));
  assert(!strcmp(b2.GetName(), "Solid VecGeomHype #2"));

  Vec_t minExtent, maxExtent;
  b1.Extent(minExtent, maxExtent);
  double tSTout = std::tan(stOut);
  double tSTin  = std::tan(stIn);

  double xy = std::sqrt(fRmax * fRmax + tSTout * tSTout * halfZ * halfZ);
  assert(ApproxEqual(minExtent, Vec_t(-xy, -xy, -halfZ)));
  assert(ApproxEqual(maxExtent, Vec_t(xy, xy, halfZ)));

  // Check Surface Normal
  Vec_t normal;
  bool valid;

  double Dist;

  // COMMENTING NORMAL FOR THE TIME BEING
  // These Normal tests Needs a relook ---

  valid = b1.Normal(pbigx, normal);
  assert(ApproxEqual(normal, Vec_t(0, 0, 1)));

  valid = b1.Normal(pbigy, normal);
  assert(ApproxEqual(normal, Vec_t(0, 0, 1)));

  valid = b1.Normal(pbigz, normal);
  assert(ApproxEqual(normal, Vec_t(0, 0, 1)));

  valid = b1.Normal(pbigmx, normal);
  assert(ApproxEqual(normal, Vec_t(0, 0, 1)));

  valid = b1.Normal(pbigmy, normal);
  assert(ApproxEqual(normal, Vec_t(0, 0, 1)));

  valid = b1.Normal(pbigmz, normal);
  assert(ApproxEqual(normal, Vec_t(0, 0, -1)));

  // valid = b1.Normal(pony,normal);
  // assert(ApproxEqual(normal,Vec_t(0,1,0)));

  valid = b1.Normal(ponz, normal);
  assert(ApproxEqual(normal, Vec_t(0, 0, 1)));

  valid = b1.Normal(ponmx, normal);
  // assert(ApproxEqual(normal,Vec_t(-1,0,0)));

  // valid = b1.Normal(ponmy,normal);
  // assert(ApproxEqual(normal,Vec_t(0,-1,0)));

  valid = b1.Normal(ponmz, normal);
  assert(ApproxEqual(normal, Vec_t(0, 0, -1)));

  // DistanceToOut(P,V) with asserts for norm and convex

  Vec_t pmidx(fRmin + ((fRmax - fRmin) / 2), 0, 0);
  Vec_t pmidy(0, fRmin + ((fRmax - fRmin) / 2), 0);

  Dist  = b1.DistanceToOut(pmidx, vx);
  valid = b1.Normal(pmidx + Dist * vx, normal);
  assert(ApproxEqual(Dist, ((fRmax - fRmin) / 2)) && ApproxEqual(normal, vx));

  Dist  = b1.DistanceToOut(-pmidx, vmx);
  valid = b1.Normal(-pmidx + Dist * vmx, normal);
  assert(ApproxEqual(Dist, ((fRmax - fRmin) / 2)) && ApproxEqual(normal, vmx));

  Dist  = b1.DistanceToOut(pmidy, vy);
  valid = b1.Normal(pmidy + Dist * vy, normal);
  assert(ApproxEqual(Dist, ((fRmax - fRmin) / 2)) && ApproxEqual(normal, vy));

  Dist  = b1.DistanceToOut(-pmidy, vmy);
  valid = b1.Normal(-pmidy + Dist * vmy, normal);
  assert(ApproxEqual(Dist, ((fRmax - fRmin) / 2)) && ApproxEqual(normal, vmy));

  double distZ = std::sqrt((pmidx.Mag2() - (fRmin * fRmin)) / (tSTin * tSTin));

  Dist  = b1.DistanceToOut(pmidx, vz);
  valid = b1.Normal(pmidx + Dist * vz, normal);
  assert(ApproxEqual(Dist, distZ)); //&& ApproxEqual(normal,vz));

  Dist = b1.DistanceToOut(pmidx, vmz);
  assert(ApproxEqual(Dist, distZ)); //&& ApproxEqual(normal,vmz));

  Dist = b1.DistanceToOut(-pmidx, vz);
  assert(ApproxEqual(Dist, distZ)); //&& ApproxEqual(normal,vz));

  Dist = b1.DistanceToOut(-pmidx, vmz);
  assert(ApproxEqual(Dist, distZ)); //&& ApproxEqual(normal,vmz));

  Dist = b1.DistanceToOut(pmidy, vz);
  assert(ApproxEqual(Dist, distZ)); //&& ApproxEqual(normal,vz));

  Dist = b1.DistanceToOut(pmidy, vmz);
  assert(ApproxEqual(Dist, distZ)); //&& ApproxEqual(normal,vmz));

  Dist = b1.DistanceToOut(-pmidy, vz);
  assert(ApproxEqual(Dist, distZ)); //&& ApproxEqual(normal,vz));

  Dist = b1.DistanceToOut(-pmidy, vmz);
  assert(ApproxEqual(Dist, distZ)); //&& ApproxEqual(normal,vmz));

  // Point is already outside and checking DistanceToOut. In this case distance is shoule be set to -1.

  Dist = b1.DistanceToOut(pbigx, vx);
  std::cout << "Dist : " << Dist << std::endl;
  assert(ApproxEqual(Dist, -1.));

  Dist = b1.DistanceToOut(pbigy, vy);
  assert(ApproxEqual(Dist, -1.));

  Dist = b1.DistanceToOut(pbigz, vz);
  assert(ApproxEqual(Dist, -1.));

  Dist = b1.DistanceToOut(pbigx, vmx);
  assert(ApproxEqual(Dist, -1.));

  Dist = b1.DistanceToOut(pbigy, vmy);
  assert(ApproxEqual(Dist, -1.));

  Dist = b1.DistanceToOut(pbigz, vmz);
  assert(ApproxEqual(Dist, -1.));

  // Check Inside
  assert(b1.Inside(pzero) == vecgeom::EInside::kOutside);
  assert(b1.Inside(pbigx) == vecgeom::EInside::kOutside);
  assert(b1.Inside(pbigy) == vecgeom::EInside::kOutside);
  assert(b1.Inside(pbigz) == vecgeom::EInside::kOutside);

  assert(b1.Inside(pmidx) == vecgeom::EInside::kInside);
  assert(b1.Inside(-pmidx) == vecgeom::EInside::kInside);
  assert(b1.Inside(pmidy) == vecgeom::EInside::kInside);
  assert(b1.Inside(-pmidy) == vecgeom::EInside::kInside);

  assert(b1.Inside(ponxsideO) == vecgeom::EInside::kSurface);
  assert(b1.Inside(ponysideO) == vecgeom::EInside::kSurface);
  assert(b1.Inside(ponxsideI) == vecgeom::EInside::kSurface);
  assert(b1.Inside(ponysideI) == vecgeom::EInside::kSurface);

  Dist = b1.DistanceToIn(pbigx, vmx);
  assert(ApproxEqual(Dist, 100 - fRmax));
  Dist = b1.DistanceToIn(pbigmx, vx);
  assert(ApproxEqual(Dist, 100 - fRmax));
  Dist = b1.DistanceToIn(pbigy, vmy);
  assert(ApproxEqual(Dist, 100 - fRmax));
  Dist = b1.DistanceToIn(pbigmy, vy);
  assert(ApproxEqual(Dist, 100 - fRmax));

  Dist = b1.DistanceToIn(pbigz, vmz);
  if (Dist >= kInfLength) Dist = kInfLength;
  assert(ApproxEqual(Dist, kInfLength));
  Dist = b1.DistanceToIn(pbigmz, vz);
  if (Dist >= kInfLength) Dist = kInfLength;
  assert(ApproxEqual(Dist, kInfLength));

  Dist = b1.DistanceToIn(pbigx, vxy);
  if (Dist >= kInfLength) Dist = kInfLength;
  assert(ApproxEqual(Dist, kInfLength));

  Dist = b1.DistanceToIn(pbigmx, vxy);
  if (Dist >= kInfLength) Dist = kInfLength;
  assert(ApproxEqual(Dist, kInfLength));

  std::cout << "------------------------------------------------" << std::endl;
  Vec_t pointOTol(20 + vecgeom::cxx::kTolerance, 0., 0.);
  Dist = b1.DistanceToIn(pointOTol, vx);
  std::cout << "DistToIn : " << Dist << std::endl;
  if (Dist >= kInfLength) Dist = kInfLength;
  assert(ApproxEqual(Dist, kInfLength));

  Dist = b1.DistanceToOut(pointOTol, vx); // This case fails
  std::cout << "DistToOut : " << Dist << std::endl;
  assert(ApproxEqual(Dist, 0.));

  // Point inside outer tolerance of outer hyperboloid  and directing in
  Dist = b1.DistanceToIn(pointOTol, vmx);
  std::cout << "DistToIn : " << Dist << std::endl;
  assert(ApproxEqual(Dist, 0.));

  Dist = b1.DistanceToOut(pointOTol, vmx); // This case fails
  std::cout << "DistToOut : " << Dist << std::endl;
  assert(ApproxEqual(Dist, 10.));

  // Point outside inner tolerance of outer hyperboloid and directing out
  std::cout << "------------------------------------------------" << std::endl;
  Vec_t pointOTolI(20 - vecgeom::cxx::kTolerance, 0., 0.);
  Dist = b1.DistanceToIn(pointOTolI, vx);
  std::cout << "DistToIn : " << Dist << std::endl;
  if (Dist >= kInfLength) Dist = kInfLength;
  assert(ApproxEqual(Dist, kInfLength));

  Dist = b1.DistanceToOut(pointOTolI, vx); // This case fails
  std::cout << "DistToOut : " << Dist << std::endl;
  assert(ApproxEqual(Dist, 0.));

  // Point Outside inner tolerance of outer hyperboloid and directing in

  Dist = b1.DistanceToIn(pointOTolI, vmx);
  std::cout << "DistToIn : " << Dist << std::endl;
  assert(ApproxEqual(Dist, 0.));

  Dist = b1.DistanceToOut(pointOTolI, vmx);
  std::cout << "DistToOut : " << Dist << std::endl;
  assert(ApproxEqual(Dist, 10.));

  std::cout << "------------------------------------------------" << std::endl;
  Vec_t pointOTol_IH(10 + vecgeom::cxx::kTolerance, 0., 0.);
  Dist = b1.DistanceToIn(pointOTol_IH, vx);
  std::cout << "DistToIn : " << Dist << std::endl;

  Dist = b1.DistanceToOut(pointOTol_IH, vx); // This case fails
  std::cout << "DistToOut : " << Dist << std::endl;
  assert(ApproxEqual(Dist, 10.));

  Dist = b1.DistanceToIn(pointOTol_IH, vmx); // May fail
  std::cout << "DistToIn : " << Dist << std::endl;
  assert(ApproxEqual(Dist, 20.));

  Dist = b1.DistanceToOut(pointOTol_IH, vmx); // This case fails
  std::cout << "DistToOut : " << Dist << std::endl;
  assert(ApproxEqual(Dist, 0.));

  std::cout << "------------------------------------------------" << std::endl;
  Vec_t pointOTol_IHm(10 - vecgeom::cxx::kTolerance, 0., 0.);
  Dist = b1.DistanceToIn(pointOTol_IHm, vx);
  std::cout << "DistToIn : " << Dist << std::endl;
  assert(ApproxEqual(Dist, 0.));

  Dist = b1.DistanceToOut(pointOTol_IHm, vx); // This case fails
  std::cout << "DistToOut : " << Dist << std::endl;
  assert(ApproxEqual(Dist, 10.));
  Dist = b1.DistanceToIn(pointOTol_IHm, vmx); // May fail
  std::cout << "DistToIn : " << Dist << std::endl;
  assert(ApproxEqual(Dist, 20.));

  Dist = b1.DistanceToOut(pointOTol_IH, vmx); // This case fails
  std::cout << "DistToOut : " << Dist << std::endl;
  assert(ApproxEqual(Dist, 0.));

  std::cout << "------------------------------------------------" << std::endl;
  Vec_t pointOTol_Z(60, 0., 50. + vecgeom::cxx::kTolerance);

  Dist = b1.DistanceToIn(pointOTol_Z, vz);
  if (Dist >= kInfLength) Dist = kInfLength;
  std::cout << "DistToIn : " << Dist << std::endl;
  assert(ApproxEqual(Dist, kInfLength));

  Dist = b1.DistanceToOut(pointOTol_Z, vz); // This case fails
  std::cout << "DistToOut : " << Dist << std::endl;
  assert(ApproxEqual(Dist, 0.));

  Dist = b1.DistanceToIn(pointOTol_Z, vmz);
  std::cout << "DistToIn : " << Dist << std::endl;
  assert(ApproxEqual(Dist, 0.));

  double Dist3 = b1.DistanceToOut(pointOTol_Z, vmz); // This case fails
  std::cout << "DistToOut : " << Dist3 << std::endl;
  assert(ApproxEqual(Dist, 0.));

  std::cout << "------------------------------------------------" << std::endl;

  Vec_t pointOTol_ZNeg(60, 0., -50. - vecgeom::cxx::kTolerance);

  Dist = b1.DistanceToIn(pointOTol_ZNeg, vmz);
  if (Dist >= kInfLength) Dist = kInfLength;
  std::cout << "DistToIn : " << Dist << std::endl;
  assert(ApproxEqual(Dist, kInfLength));

  Dist = b1.DistanceToOut(pointOTol_ZNeg, vmz); // This case fails
  std::cout << "DistToOut : " << Dist << std::endl;
  assert(ApproxEqual(Dist, 0.));

  Dist = b1.DistanceToIn(pointOTol_ZNeg, vz);
  std::cout << "DistToIn : " << Dist << std::endl;
  assert(ApproxEqual(Dist, 0.));

  double Dist4 = b1.DistanceToOut(pointOTol_ZNeg, vz); // This case fails
  std::cout << "DistToOut : " << Dist4 << std::endl;
  assert(ApproxEqual(Dist3, Dist4));

  // UNIT TESTS FOR WRONG SIDE POINTS

  // Testing DistanceToOut for outside points
  Dist = b1.DistanceToOut(pbigx, vmx);
  assert(ApproxEqual(Dist, -1.));
  Dist = b1.DistanceToOut(pbigy, vmy);
  assert(ApproxEqual(Dist, -1.));
  Dist = b1.DistanceToOut(pbigz, vmz);
  assert(ApproxEqual(Dist, -1.));
  Dist = b1.DistanceToOut(pbigmx, vmx);
  assert(ApproxEqual(Dist, -1.));
  Dist = b1.DistanceToOut(pbigmy, vmy);
  assert(ApproxEqual(Dist, -1.));
  Dist = b1.DistanceToOut(pbigmz, vmz);
  assert(ApproxEqual(Dist, -1.));

  // Testing SafetFromInside for outside points
  Dist = b1.SafetyToOut(pbigx);
  assert(ApproxEqual(Dist, -1));
  Dist = b1.SafetyToOut(pbigy);
  assert(ApproxEqual(Dist, -1));
  Dist = b1.SafetyToOut(pbigz);
  assert(ApproxEqual(Dist, -1));
  Dist = b1.SafetyToOut(pbigmx);
  assert(ApproxEqual(Dist, -1));
  Dist = b1.SafetyToOut(pbigmy);
  assert(ApproxEqual(Dist, -1));
  Dist = b1.SafetyToOut(pbigmz);
  assert(ApproxEqual(Dist, -1));

  // Testing DistanceToIn for inside points
  Dist = b1.DistanceToIn(Vec_t(15., 0., 0.), vx);
  assert(ApproxEqual(Dist, -1.));
  Dist = b1.DistanceToIn(Vec_t(0., 15., 0.), vy);
  assert(ApproxEqual(Dist, -1.));
  Dist = b1.DistanceToIn(Vec_t(-15., 0., 0.), vmx);
  assert(ApproxEqual(Dist, -1.));
  Dist = b1.DistanceToIn(Vec_t(0., -15., 0.), vmy);
  assert(ApproxEqual(Dist, -1.));

  // Testing SafetyFromOut for inside points
  Dist = b1.SafetyToIn(Vec_t(15., 0., 0.));
  assert(ApproxEqual(Dist, -1.));
  Dist = b1.SafetyToIn(Vec_t(0., 15., 0.));
  assert(ApproxEqual(Dist, -1.));
  Dist = b1.SafetyToIn(Vec_t(-15., 0., 0.));
  assert(ApproxEqual(Dist, -1.));
  Dist = b1.SafetyToIn(Vec_t(0., -15., 0.));
  assert(ApproxEqual(Dist, -1.));
  std::cout << valid << std::endl;
  return true;
}

int main(int argc, char *argv[])
{
  assert(TestHype<vecgeom::SimpleHype>());
  std::cout << "VecGeomHype passed\n";

  return 0;
}
