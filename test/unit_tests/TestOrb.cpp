// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.
// For the full list of authors see CONTRIBUTORS.txt and `git log`.

/// \brief Unit test for the Orb.
/// \file test/unit_tests/TestOrb.cpp
/// \author Raman Sehgal

// ensure asserts are compiled in
#undef NDEBUG
#include "VecGeom/base/FpeEnable.h"

#include "VecGeom/base/Global.h"
#include "VecGeom/base/Vector3D.h"
#include "VecGeom/volumes/Orb.h"
#include "ApproxEqual.h"

#include <cmath>

bool testvecgeom = false;

using vecgeom::kInfLength;
using vecgeom::kPi;
using vecgeom::Precision;

template <class Orb_t, class Vec_t = vecgeom::Vector3D<vecgeom::Precision>>
bool TestOrb()
{

  int verbose = 0;

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

  Vec_t ponxside(fR, 0, 0), ponyside(0, fR, 0), ponzside(0, 0, fR);
  Vec_t ponmxside(-fR, 0, 0), ponmyside(0, -fR, 0), ponmzside(0, 0, -fR);

  Vec_t vx(1, 0, 0), vy(0, 1, 0), vz(0, 0, 1);
  Vec_t vmx(-1, 0, 0), vmy(0, -1, 0), vmz(0, 0, -1);
  Vec_t vxy(1 / std::sqrt(2.0), 1 / std::sqrt(2.0), 0);
  Vec_t vmxy(-1 / std::sqrt(2.0), 1 / std::sqrt(2.0), 0);
  Vec_t vmxmy(-1 / std::sqrt(2.0), -1 / std::sqrt(2.0), 0);
  Vec_t vxmy(1 / std::sqrt(2.0), -1 / std::sqrt(2.0), 0);
  Vec_t vxmz(1 / std::sqrt(2.0), 0, -1 / std::sqrt(2.0));

  Orb_t b1("Solid VecGeomOrb #1", fR);
  Orb_t b2("Solid VecGeomOrb #2", 6);

  // Check cubic volume
  assert(ApproxEqual<Precision>(b1.Capacity(), (4 * kPi / 3) * fR * fR * fR));
  assert(ApproxEqual<Precision>(b2.Capacity(), (4 * kPi / 3) * 6 * 6 * 6));

  // Check Surface area
  assert(ApproxEqual<Precision>(b1.SurfaceArea(), ((4 * kPi) * fR * fR)));
  assert(ApproxEqual<Precision>(b2.SurfaceArea(), ((4 * kPi) * 6 * 6)));

  // Check Extent and cached BBox
  Vec_t minExtent, maxExtent;
  Vec_t minBBox, maxBBox;
  b1.Extent(minExtent, maxExtent);
  b1.GetUnplacedVolume()->GetBBox(minBBox, maxBBox);
  assert(ApproxEqual(minExtent, Vec_t(-fR, -fR, -fR)));
  assert(ApproxEqual(maxExtent, Vec_t(fR, fR, fR)));
  assert(ApproxEqual(minExtent, minBBox));
  assert(ApproxEqual(maxExtent, maxBBox));
  b2.Extent(minExtent, maxExtent);
  b2.GetUnplacedVolume()->GetBBox(minBBox, maxBBox);
  assert(ApproxEqual(minExtent, Vec_t(-6, -6, -6)));
  assert(ApproxEqual(maxExtent, Vec_t(6, 6, 6)));
  assert(ApproxEqual(minExtent, minBBox));
  assert(ApproxEqual(maxExtent, maxBBox));

  // Check Surface Normal
  Vec_t normal;
  bool valid;

  double Dist;

  valid = b1.Normal(ponx, normal);
  assert(ApproxEqual(normal, Vec_t(1, 0, 0)));
  assert(valid);

  valid = b1.Normal(pony, normal);
  assert(ApproxEqual(normal, Vec_t(0, 1, 0)));
  assert(valid);

  valid = b1.Normal(ponz, normal);
  assert(ApproxEqual(normal, Vec_t(0, 0, 1)));
  assert(valid);

  valid = b1.Normal(ponmx, normal);
  assert(ApproxEqual(normal, Vec_t(-1, 0, 0)));
  assert(valid);

  valid = b1.Normal(ponmy, normal);
  assert(ApproxEqual(normal, Vec_t(0, -1, 0)));
  assert(valid);

  valid = b1.Normal(ponmz, normal);
  assert(ApproxEqual(normal, Vec_t(0, 0, -1)));
  assert(valid);
  // DistanceToOut(P,V) with asserts for norm and convex
  Dist  = b1.DistanceToOut(pzero, vx);
  valid = b1.Normal(pzero + Dist * vx, normal);
  assert(ApproxEqual<Precision>(Dist, fR) && ApproxEqual(normal, vx));

  Dist  = b1.DistanceToOut(pzero, vmx);
  valid = b1.Normal(pzero + Dist * vmx, normal);
  assert(ApproxEqual<Precision>(Dist, fR) && ApproxEqual(normal, vmx));

  Dist  = b1.DistanceToOut(pzero, vy);
  valid = b1.Normal(pzero + Dist * vy, normal);
  assert(ApproxEqual<Precision>(Dist, fR) && ApproxEqual(normal, vy));

  Dist  = b1.DistanceToOut(pzero, vmy);
  valid = b1.Normal(pzero + Dist * vmy, normal);
  assert(ApproxEqual<Precision>(Dist, fR) && ApproxEqual(normal, vmy));

  Dist  = b1.DistanceToOut(pzero, vz);
  valid = b1.Normal(pzero + Dist * vz, normal);
  assert(ApproxEqual<Precision>(Dist, fR) && ApproxEqual(normal, vz));

  Dist  = b1.DistanceToOut(pzero, vmz);
  valid = b1.Normal(pzero + Dist * vmz, normal);
  assert(ApproxEqual<Precision>(Dist, fR) && ApproxEqual(normal, vmz));

  Dist  = b1.DistanceToOut(ponxside, vx);
  valid = b1.Normal(ponxside + Dist * vx, normal);
  assert(ApproxEqual<Precision>(Dist, 0) && ApproxEqual(normal, vx));

  Dist  = b1.DistanceToOut(ponxside, vmx);
  valid = b1.Normal(ponxside + Dist * vmx, normal);
  assert(ApproxEqual<Precision>(Dist, 2 * fR) && ApproxEqual(normal, vmx));

  Dist  = b1.DistanceToOut(ponmxside, vx);
  valid = b1.Normal(ponmxside + Dist * vx, normal);
  assert(ApproxEqual<Precision>(Dist, 2 * fR) && ApproxEqual(normal, vx));

  Dist  = b1.DistanceToOut(ponmxside, vmx);
  valid = b1.Normal(ponmxside + Dist * vmx, normal);
  assert(ApproxEqual<Precision>(Dist, 0) && ApproxEqual(normal, vmx));

  Dist  = b1.DistanceToOut(ponyside, vy);
  valid = b1.Normal(ponyside + Dist * vy, normal);
  assert(ApproxEqual<Precision>(Dist, 0) && ApproxEqual(normal, vy));

  Dist  = b1.DistanceToOut(ponyside, vmy);
  valid = b1.Normal(ponyside + Dist * vmy, normal);
  assert(ApproxEqual<Precision>(Dist, 2 * fR) && ApproxEqual(normal, vmy));

  Dist  = b1.DistanceToOut(ponmyside, vy);
  valid = b1.Normal(ponmyside + Dist * vy, normal);
  assert(ApproxEqual<Precision>(Dist, 2 * fR) && ApproxEqual(normal, vy));

  Dist  = b1.DistanceToOut(ponmyside, vmy);
  valid = b1.Normal(ponmyside + Dist * vmy, normal);
  assert(ApproxEqual<Precision>(Dist, 0) && ApproxEqual(normal, vmy));

  Dist  = b1.DistanceToOut(ponzside, vz);
  valid = b1.Normal(ponzside + Dist * vz, normal);
  assert(ApproxEqual<Precision>(Dist, 0) && ApproxEqual(normal, vz));

  Dist  = b1.DistanceToOut(ponzside, vmz);
  valid = b1.Normal(ponzside + Dist * vmz, normal);
  assert(ApproxEqual<Precision>(Dist, 2 * fR) && ApproxEqual(normal, vmz));

  Dist  = b1.DistanceToOut(ponmzside, vz);
  valid = b1.Normal(ponmzside + Dist * vz, normal);
  assert(ApproxEqual<Precision>(Dist, 2 * fR) && ApproxEqual(normal, vz));

  Dist  = b1.DistanceToOut(ponmzside, vmz);
  valid = b1.Normal(ponmzside + Dist * vmz, normal);
  assert(ApproxEqual<Precision>(Dist, 0) && ApproxEqual(normal, vmz));

  // Check Inside
  assert(b1.Inside(pzero) == vecgeom::EInside::kInside);
  assert(b1.Inside(pbigx) == vecgeom::EInside::kOutside);
  assert(b1.Inside(pbigy) == vecgeom::EInside::kOutside);
  assert(b1.Inside(pbigz) == vecgeom::EInside::kOutside);

  assert(b1.Inside(ponxside) == vecgeom::EInside::kSurface);
  assert(b1.Inside(ponyside) == vecgeom::EInside::kSurface);
  assert(b1.Inside(ponzside) == vecgeom::EInside::kSurface);

  assert(b2.Inside(pzero) == vecgeom::EInside::kInside);
  assert(b2.Inside(ponxside) == vecgeom::EInside::kOutside);
  assert(b2.Inside(ponyside) == vecgeom::EInside::kOutside);
  assert(b2.Inside(ponzside) == vecgeom::EInside::kOutside);
  assert(b2.Inside(Vec_t(6, 0, 0)) == vecgeom::EInside::kSurface);
  assert(b2.Inside(Vec_t(0, 6, 0)) == vecgeom::EInside::kSurface);
  assert(b2.Inside(Vec_t(0, 0, 6)) == vecgeom::EInside::kSurface);

  // SafetyToOut(P)
  Dist = b1.SafetyToOut(pzero);
  assert(ApproxEqual<Precision>(Dist, fR));
  Dist = b1.SafetyToOut(vx);
  assert(ApproxEqual<Precision>(Dist, fR - 1));
  Dist = b1.SafetyToOut(vy);
  assert(ApproxEqual<Precision>(Dist, fR - 1));
  Dist = b1.SafetyToOut(vz);
  assert(ApproxEqual<Precision>(Dist, fR - 1));

  // SafetyToIn(P)
  Dist = b1.SafetyToIn(pbigx);
  assert(ApproxEqual<Precision>(Dist, 100 - fR));
  Dist = b1.SafetyToIn(pbigmx);
  assert(ApproxEqual<Precision>(Dist, 100 - fR));
  Dist = b1.SafetyToIn(pbigy);
  assert(ApproxEqual<Precision>(Dist, 100 - fR));
  Dist = b1.SafetyToIn(pbigmy);
  assert(ApproxEqual<Precision>(Dist, 100 - fR));
  Dist = b1.SafetyToIn(pbigz);
  assert(ApproxEqual<Precision>(Dist, 100 - fR));
  Dist = b1.SafetyToIn(pbigmz);
  assert(ApproxEqual<Precision>(Dist, 100 - fR));

  // DistanceToIn(P,V)
  Dist = b1.DistanceToIn(pbigx, vmx);
  assert(ApproxEqual<Precision>(Dist, 100 - fR));
  Dist = b1.DistanceToIn(pbigmx, vx);
  assert(ApproxEqual<Precision>(Dist, 100 - fR));
  Dist = b1.DistanceToIn(pbigy, vmy);
  assert(ApproxEqual<Precision>(Dist, 100 - fR));
  Dist = b1.DistanceToIn(pbigmy, vy);
  assert(ApproxEqual<Precision>(Dist, 100 - fR));
  Dist = b1.DistanceToIn(pbigz, vmz);
  assert(ApproxEqual<Precision>(Dist, 100 - fR));
  Dist = b1.DistanceToIn(pbigmz, vz);
  assert(ApproxEqual<Precision>(Dist, 100 - fR));

  Dist = b1.DistanceToIn(pbigx, vxy);
  if (Dist >= kInfLength) Dist = kInfLength;
  assert(ApproxEqual<Precision>(Dist, kInfLength));

  Dist = b1.DistanceToIn(pbigmx, vxy);
  if (Dist >= kInfLength) Dist = kInfLength;
  assert(ApproxEqual<Precision>(Dist, kInfLength));

  Vec_t pJohnXZ(9, 0, 12);
  Dist = b2.DistanceToIn(pJohnXZ, vxmz);
  if (Dist >= kInfLength) Dist = kInfLength;
  assert(ApproxEqual<Precision>(Dist, kInfLength));

  Vec_t pJohnXY(12, 9, 0);
  Dist = b2.DistanceToIn(pJohnXY, vmxy);
  if (Dist >= kInfLength) Dist = kInfLength;
  assert(ApproxEqual<Precision>(Dist, kInfLength));

  Dist = b2.DistanceToIn(pJohnXY, vmx);
  if (Dist >= kInfLength) Dist = kInfLength;
  assert(ApproxEqual<Precision>(Dist, kInfLength));

  Vec_t pJohnX(8, 0, 0);
  Dist = b2.DistanceToIn(pJohnX, vmx);
  assert(ApproxEqual<Precision>(Dist, 2));

  Vec_t p2JohnXY(7, 5, 0);
  Dist = b2.DistanceToIn(p2JohnXY, vmx);
  assert(ApproxEqual<Precision>(Dist, 3.6833752));

  Dist = b1.DistanceToIn(Vec_t(-25, -35, 0), vx);
  if (Dist >= kInfLength) Dist = kInfLength;
  assert(ApproxEqual<Precision>(Dist, kInfLength));

  Dist = b1.DistanceToIn(Vec_t(-25, -35, 0), vy);
  if (Dist >= kInfLength) Dist = kInfLength;
  assert(ApproxEqual<Precision>(Dist, kInfLength));

  Vec_t pointO(-8.363470934547895, 2.754420966126675, -2.665617952433236);
  Vec_t dirO(-8.363470934547895 / 9.2, 2.754420966126675 / 9.2, -2.665617952433236 / 9.2);
  Dist = b1.DistanceToIn(pointO, dirO);
  if (Dist >= kInfLength) Dist = kInfLength;
  assert(ApproxEqual<Precision>(Dist, kInfLength));

  Dist = b1.DistanceToOut(pointO, dirO);
  if (verbose) std::cout << "DistanceToOut is : " << Dist << std::endl;
  assert(Dist <= 0.);

  if (verbose) std::cout << " Now testing point out directing in " << std::endl;
  Vec_t dirI(8.363470934547895 / 9.2, -2.754420966126675 / 9.2, 2.665617952433236 / 9.2);
  Dist = b1.DistanceToIn(pointO, dirI);
  assert(ApproxEqual<Precision>(Dist, 0.2));
  if (verbose) std::cout << "DistanceToIn PODI is : " << Dist << std::endl;

  Dist = b1.DistanceToOut(pointO, dirO);
  assert(Dist <= 0.);
  if (verbose) std::cout << "DistanceToOut PODI is : " << Dist << std::endl;

  // Point inside outer tolerance and directing out
  if (verbose) std::cout << "Testing point inside outer tolerance and directing out" << std::endl;
  Vec_t pointOTol(8.884242447222299, 0.134875592787852, -1.432495973274375);
  Vec_t dirOTol(8.884242447222299 / 9, 0.134875592787852 / 9, -1.432495973274375 / 9);
  Dist = b1.DistanceToIn(pointOTol, dirOTol);
  if (Dist >= kInfLength) Dist = kInfLength;
  assert(ApproxEqual<Precision>(Dist, kInfLength));
  if (verbose) std::cout << "DistanceToIn for point inside outer tolerance and directing out: " << Dist << std::endl;

  // Point inside outer tolerance and directing in
  if (verbose) std::cout << "Testing point inside outer tolerance and directing IN" << std::endl;
  Vec_t dirOTolI(-8.884242447222299 / 9, -0.134875592787852 / 9, 1.432495973274375 / 9);
  Dist = b1.DistanceToIn(pointOTol, dirOTolI);
  assert(ApproxEqual<Precision>(Dist, 0));
  if (verbose) std::cout << "DistanceToIn for point inside outer tolerance and directing IN: " << Dist << std::endl;

  // Point inside  and directing out
  if (verbose) std::cout << "Testing point inside and directing OUT" << std::endl;
  Vec_t pointI(-3.618498437781364, 2.401810108299175, -6.718465394675017);
  Vec_t dirIO(-3.618498437781364 / 8, 2.401810108299175 / 8, -6.718465394675017 / 8);
  Dist = b1.DistanceToIn(pointI, dirIO);
  if (Dist >= kInfLength) Dist = kInfLength;
  // assert(ApproxEqual<Precision>(Dist,kInfLength));
  assert(Dist < 0.);
  if (verbose) std::cout << "DistanceToIn for point inside and directing OUT: " << Dist << std::endl;

  Dist = b1.DistanceToOut(pointI, dirIO);
  assert(ApproxEqual<Precision>(Dist, 1));
  if (verbose) std::cout << "DistanceToOut for point inside and directing OUT: " << Dist << std::endl;

  // Point inside and directing in
  if (verbose) std::cout << "Testing point inside and directing IN" << std::endl;
  Vec_t dirII(3.618498437781364 / 8, -2.401810108299175 / 8, 6.718465394675017 / 8);
  Dist = b1.DistanceToIn(pointI, dirII);
  if (Dist >= kInfLength) Dist = kInfLength;
  // assert(ApproxEqual<Precision>(Dist,kInfLength));
  assert(Dist < 0.);

  Dist = b1.DistanceToOut(pointI, dirII);
  assert(ApproxEqual<Precision>(Dist, 17));

  // Testing Surface point generated by shape tester with Random direction
  Orb_t b3("Solid VecGeomOrb #3", 3);
  Vec_t tSrPoint(-1.704652541027918744, 0.92532982039202238411, -2.2886512267834184797);
  Vec_t tSrDir(0.55992127966252225324, -0.73207465430510332283, 0.38801399601708530529);
  Dist = b3.DistanceToIn(tSrPoint, tSrDir);
  assert(ApproxEqual<Precision>(Dist, 0));

  // Testing TAD Unqualified points by shape tester
  Orb_t b4("Solid VecGeomOrb #4", 3);

  // Trying T0  points by shape tester
  // Direction is random
  Orb_t b5("Solid VecGeomOrb #5", 10.);
  Vec_t t0UQualPoint(-89.066635313481228309, -125.66666807765292901, 105.36991856050774174);
  Vec_t t0UQualDir(0.45467415240829972545, 0.6924777322970485649, -0.56013034679843187735);
  double shiftDist = b5.DistanceToIn(t0UQualPoint, t0UQualDir);
  Vec_t t0UPoint(t0UQualPoint + shiftDist * t0UQualDir);
  Dist = b5.DistanceToIn(t0UPoint, t0UQualDir);
  assert(ApproxEqual<Precision>(Dist, 0));

  // Trying T0 Test of Shape tester.
  // Inside point chosen is 0.,0.,0.
  Vec_t t0KnownPoint(-89.066635313481228309, -125.66666807765292901, 105.36991856050774174);
  double vecLen = t0KnownPoint.Mag();
  Vec_t t0Dir((1 / vecLen) * (t0KnownPoint * (-1)));
  shiftDist = b5.DistanceToIn(t0KnownPoint, t0Dir);
  Vec_t t0Point(t0KnownPoint + shiftDist * t0Dir);
  Dist = b5.DistanceToIn(t0Point, t0Dir);
  assert(ApproxEqual<Precision>(Dist, 0));

  // Trying T0 test of Shape Tester
  Orb_t b6("Solid VecGeomOrb #6", 50.);
  Vec_t t0KnownPoint2(725.77514262389081523, -745.20572125476428482, 778.25732374796621116);
  vecLen = t0KnownPoint2.Mag();
  Vec_t t0Dir2((1 / vecLen) * (t0KnownPoint2 * (-1)));
  // Propagate to the solid first
  shiftDist = b6.GetUnplacedVolume()->ApproachSolid(t0KnownPoint2, 1 / t0Dir2);
  Vec_t t0Point2(t0KnownPoint2 + shiftDist * t0Dir2);
  double shiftDist2 = b6.DistanceToIn(t0Point2, t0Dir2);
  Vec_t t0Point3(t0Point2 + shiftDist2 * t0Dir2);
  Dist = b6.DistanceToIn(t0Point3, t0Dir2);
  assert(ApproxEqual<Precision>(Dist, 0));

  // Trying TS test of Shape Tester
  Vec_t tSKnownPoint(-46.418969559619192466, -17.038092391528927294, 7.4150301880992222081);
  Vec_t tSDir(0.84546676256780450842, -0.33558962646945966757, -0.41541010579811871173);
  Dist = b6.DistanceToIn(tSKnownPoint, tSDir);
  assert(ApproxEqual<Precision>(Dist, 0));

  return true;
}

int main(int argc, char *argv[])
{
  assert(TestOrb<vecgeom::SimpleOrb>());
  std::cout << "VecGeomOrb passed\n";

  return 0;
}
