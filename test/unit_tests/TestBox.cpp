//
// File:    TestBox.cpp
// Purpose: Unit tests for the box
//

//-- ensure asserts are compiled in
#undef NDEBUG

#include "base/Vector3D.h"
#include "volumes/Box.h"
#include "ApproxEqual.h"
#include <cmath>

#if defined(__GNUCC__) && !defined(__CLANG__)
#include <fenv.h>
#endif

using Vec_t = vecgeom::Vector3D<vecgeom::Precision>;

bool testvecgeom = false;

template <class Box_t>
bool Test_VECGEOM_431();

template <class Box_t, class Vec_t = vecgeom::Vector3D<vecgeom::Precision>>
bool TestBox()
{

  Vec_t pzero(0, 0, 0);
  Vec_t ponxside(20, 0, 0), ponyside(0, 30, 0), ponzside(0, 0, 40);
  Vec_t ponmxside(-20, 0, 0), ponmyside(0, -30, 0), ponmzside(0, 0, -40);
  Vec_t ponzsidey(0, 25, 40), ponmzsidey(0, 25, -40);

  Vec_t pbigx(100, 0, 0), pbigy(0, 100, 0), pbigz(0, 0, 100);
  Vec_t pbigmx(-100, 0, 0), pbigmy(0, -100, 0), pbigmz(0, 0, -100);

  Vec_t vx(1, 0, 0), vy(0, 1, 0), vz(0, 0, 1);
  Vec_t vmx(-1, 0, 0), vmy(0, -1, 0), vmz(0, 0, -1);
  Vec_t vxy(1 / std::sqrt(2.0), 1 / std::sqrt(2.0), 0);
  Vec_t vmxy(-1 / std::sqrt(2.0), 1 / std::sqrt(2.0), 0);
  Vec_t vmxmy(-1 / std::sqrt(2.0), -1 / std::sqrt(2.0), 0);
  Vec_t vxmy(1 / std::sqrt(2.0), -1 / std::sqrt(2.0), 0);
  Vec_t vxmz(1 / std::sqrt(2.0), 0, -1 / std::sqrt(2.0));

  double Dist;
  Vec_t norm;
  bool convex;

  Box_t b1("Test Box #1", 20, 30, 40);
  Box_t b2("Test Box #2", 10, 10, 10);
  Box_t box3("BABAR Box", 0.14999999999999999, 24.707000000000001, 22.699999999999999);

  // Check cubic volume

  assert(b2.Capacity() == 8000);
  assert(b1.Capacity() == 192000);

  // Check Surface area

  assert(b1.SurfaceArea() == 20800);
  assert(b2.SurfaceArea() == 6 * 20 * 20);

  // CalculateExtent

  Vec_t minExtent, maxExtent;
  b1.Extent(minExtent, maxExtent);
  assert(ApproxEqual(minExtent, Vec_t(-20, -30, -40)));
  assert(ApproxEqual(maxExtent, Vec_t(20, 30, 40)));
  b2.Extent(minExtent, maxExtent);
  assert(ApproxEqual(minExtent, Vec_t(-10, -10, -10)));
  assert(ApproxEqual(maxExtent, Vec_t(10, 10, 10)));

  // Check Surface Normal
  Vec_t normal;
  bool valid;
  // Normals on Surface
  valid = b1.Normal(ponxside, normal);
  assert(ApproxEqual(normal, Vec_t(1, 0, 0)));
  valid = b1.Normal(ponmxside, normal);
  assert(ApproxEqual(normal, Vec_t(-1, 0, 0)));
  valid = b1.Normal(ponyside, normal);
  assert(ApproxEqual(normal, Vec_t(0, 1, 0)));
  valid = b1.Normal(ponmyside, normal);
  assert(ApproxEqual(normal, Vec_t(0, -1, 0)));
  valid = b1.Normal(ponzside, normal);
  assert(ApproxEqual(normal, Vec_t(0, 0, 1)));
  valid = b1.Normal(ponmzside, normal);
  assert(ApproxEqual(normal, Vec_t(0, 0, -1)));
  valid = b1.Normal(ponzsidey, normal);
  assert(ApproxEqual(normal, Vec_t(0, 0, 1)));
  valid = b1.Normal(ponmzsidey, normal);
  assert(ApproxEqual(normal, Vec_t(0, 0, -1)));

  // Normals on Edges
  Vec_t edgeXY(20.0, 30., 0.0);
  Vec_t edgemXmY(-20.0, -30., 0.0);
  Vec_t edgeXmY(20.0, -30., 0.0);
  Vec_t edgemXY(-20.0, 30., 0.0);
  Vec_t edgeXZ(20.0, 0.0, 40.0);
  Vec_t edgemXmZ(-20.0, 0.0, -40.0);
  Vec_t edgeXmZ(20.0, 0.0, -40.0);
  Vec_t edgemXZ(-20.0, 0.0, 40.0);
  Vec_t edgeYZ(0.0, 30.0, 40.0);
  Vec_t edgemYmZ(0.0, -30.0, -40.0);
  Vec_t edgeYmZ(0.0, 30.0, -40.0);
  Vec_t edgemYZ(0.0, -30.0, 40.0);

  double invSqrt2 = 1.0 / std::sqrt(2.0);
  double invSqrt3 = 1.0 / std::sqrt(3.0);

  valid = b1.Normal(edgeXY, normal);
  assert(valid == true);
  assert(ApproxEqual(normal, Vec_t(invSqrt2, invSqrt2, 0.0)));
  valid = b1.Normal(edgemXmY, normal);
  assert(ApproxEqual(normal, Vec_t(-invSqrt2, -invSqrt2, 0.0)));
  valid = b1.Normal(edgeXmY, normal);
  assert(ApproxEqual(normal, Vec_t(invSqrt2, -invSqrt2, 0.0)));
  valid = b1.Normal(edgemXY, normal);
  assert(ApproxEqual(normal, Vec_t(-invSqrt2, invSqrt2, 0.0)));

  valid = b1.Normal(edgeXZ, normal);
  assert(ApproxEqual(normal, Vec_t(invSqrt2, 0.0, invSqrt2)));
  valid = b1.Normal(edgemXmZ, normal);
  assert(ApproxEqual(normal, Vec_t(-invSqrt2, 0.0, -invSqrt2)));
  valid = b1.Normal(edgeXmZ, normal);
  assert(ApproxEqual(normal, Vec_t(invSqrt2, 0.0, -invSqrt2)));
  valid = b1.Normal(edgemXZ, normal);
  assert(ApproxEqual(normal, Vec_t(-invSqrt2, 0.0, invSqrt2)));

  valid = b1.Normal(edgeYZ, normal);
  assert(ApproxEqual(normal, Vec_t(0.0, invSqrt2, invSqrt2)));
  valid = b1.Normal(edgemYmZ, normal);
  assert(ApproxEqual(normal, Vec_t(0.0, -invSqrt2, -invSqrt2)));
  valid = b1.Normal(edgeYmZ, normal);
  assert(ApproxEqual(normal, Vec_t(0.0, invSqrt2, -invSqrt2)));
  valid = b1.Normal(edgemYZ, normal);
  assert(ApproxEqual(normal, Vec_t(0.0, -invSqrt2, invSqrt2)));

  // Normals on corners
  Vec_t cornerXYZ(20.0, 30., 40.0);
  Vec_t cornermXYZ(-20.0, 30., 40.0);
  Vec_t cornerXmYZ(20.0, -30., 40.0);
  Vec_t cornermXmYZ(-20.0, -30., 40.0);
  Vec_t cornerXYmZ(20.0, 30., -40.0);
  Vec_t cornermXYmZ(-20.0, 30., -40.0);
  Vec_t cornerXmYmZ(20.0, -30., -40.0);
  Vec_t cornermXmYmZ(-20.0, -30., -40.0);

  valid = b1.Normal(cornerXYZ, normal);
  assert(ApproxEqual(normal, Vec_t(invSqrt3, invSqrt3, invSqrt3)));
  valid = b1.Normal(cornermXYZ, normal);
  assert(ApproxEqual(normal, Vec_t(-invSqrt3, invSqrt3, invSqrt3)));
  valid = b1.Normal(cornerXmYZ, normal);
  assert(ApproxEqual(normal, Vec_t(invSqrt3, -invSqrt3, invSqrt3)));
  valid = b1.Normal(cornermXmYZ, normal);
  assert(ApproxEqual(normal, Vec_t(-invSqrt3, -invSqrt3, invSqrt3)));
  valid = b1.Normal(cornerXYmZ, normal);
  assert(ApproxEqual(normal, Vec_t(invSqrt3, invSqrt3, -invSqrt3)));
  valid = b1.Normal(cornermXYmZ, normal);
  assert(ApproxEqual(normal, Vec_t(-invSqrt3, invSqrt3, -invSqrt3)));
  valid = b1.Normal(cornerXmYmZ, normal);
  assert(ApproxEqual(normal, Vec_t(invSqrt3, -invSqrt3, -invSqrt3)));
  valid = b1.Normal(cornermXmYmZ, normal);
  assert(ApproxEqual(normal, Vec_t(-invSqrt3, -invSqrt3, -invSqrt3)));

  // DistanceToOut(P,V) with asserts for norm and convex
  Dist = b1.DistanceToOut(pzero, vx, norm, convex);
  assert(ApproxEqual(Dist, 20) && ApproxEqual(norm, vx));
  if (!testvecgeom) assert(convex);
  Dist = b1.DistanceToOut(pzero, vmx, norm, convex);
  assert(ApproxEqual(Dist, 20) && ApproxEqual(norm, vmx));
  if (!testvecgeom) assert(convex);
  Dist = b1.DistanceToOut(pzero, vy, norm, convex);
  assert(ApproxEqual(Dist, 30) && ApproxEqual(norm, vy));
  if (!testvecgeom) assert(convex);
  Dist = b1.DistanceToOut(pzero, vmy, norm, convex);
  assert(ApproxEqual(Dist, 30) && ApproxEqual(norm, vmy));
  if (!testvecgeom) assert(convex);
  Dist = b1.DistanceToOut(pzero, vz, norm, convex);
  assert(ApproxEqual(Dist, 40) && ApproxEqual(norm, vz));
  if (!testvecgeom) assert(convex);
  Dist = b1.DistanceToOut(pzero, vmz, norm, convex);
  assert(ApproxEqual(Dist, 40) && ApproxEqual(norm, vmz));
  if (!testvecgeom) assert(convex);
  Dist = b1.DistanceToOut(pzero, vxy, norm, convex);
  assert(ApproxEqual(Dist, std::sqrt(800.)));
  if (!testvecgeom) assert(convex);

  // testing a few special directions ( with special zero values )
  // important since we sometimes operate with sign bits etc.
  Dist = b1.DistanceToOut(pzero, Vec_t(-0, -0, -1), norm, convex);
  assert(ApproxEqual(Dist, 40.));
  Dist = b1.DistanceToOut(pzero, Vec_t(0, -0, 1), norm, convex);
  assert(ApproxEqual(Dist, 40.));
  Dist = b1.DistanceToOut(pzero, Vec_t(0, 0, -1), norm, convex);
  assert(ApproxEqual(Dist, 40.));
  Dist = b1.DistanceToOut(pzero, Vec_t(0, 1, 0), norm, convex);
  assert(ApproxEqual(Dist, 30.));
  Dist = b1.DistanceToOut(pzero, Vec_t(0, 1, -0), norm, convex);
  assert(ApproxEqual(Dist, 30.));
  Dist = b1.DistanceToOut(pzero, Vec_t(-0, 1, 0), norm, convex);
  assert(ApproxEqual(Dist, 30.));
  Dist = b1.DistanceToOut(pzero, Vec_t(1, -0, 0), norm, convex);
  assert(ApproxEqual(Dist, 20.));
  Dist = b1.DistanceToOut(pzero, Vec_t(-1, 0, -0), norm, convex);
  assert(ApproxEqual(Dist, 20.));
  Dist = b1.DistanceToOut(pzero, Vec_t(-1, -0, -0), norm, convex);
  assert(ApproxEqual(Dist, 20.));

  Dist = b1.DistanceToOut(ponxside, vx, norm, convex);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(norm, vx));
  if (!testvecgeom) assert(convex);
  Dist = b1.DistanceToOut(ponxside, vmx, norm, convex);
  assert(ApproxEqual(Dist, 40) && ApproxEqual(norm, vmx));
  if (!testvecgeom) assert(convex);
  Dist = b1.DistanceToOut(ponmxside, vmx, norm, convex);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(norm, vmx));
  if (!testvecgeom) assert(convex);
  Dist = b1.DistanceToOut(ponyside, vy, norm, convex);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(norm, vy));
  if (!testvecgeom) assert(convex);
  Dist = b1.DistanceToOut(ponmyside, vmy, norm, convex);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(norm, vmy));
  if (!testvecgeom) assert(convex);
  Dist = b1.DistanceToOut(ponzside, vz, norm, convex);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(norm, vz));
  if (!testvecgeom) assert(convex);
  Dist = b1.DistanceToOut(ponmzside, vmz, norm, convex);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(norm, vmz));
  if (!testvecgeom) assert(convex);
  //#endif

  // Check Inside
  assert(b1.Inside(pzero) == vecgeom::EInside::kInside);
  // assert(b1.Inside(pzero)==VUSolid::eOutside);
  assert(b1.Inside(pbigz) == vecgeom::EInside::kOutside);
  assert(b1.Inside(ponxside) == vecgeom::EInside::kSurface);
  assert(b1.Inside(ponyside) == vecgeom::EInside::kSurface);
  assert(b1.Inside(ponzside) == vecgeom::EInside::kSurface);

  assert(b2.Inside(pzero) == vecgeom::EInside::kInside);
  assert(b2.Inside(pbigz) == vecgeom::EInside::kOutside);
  assert(b2.Inside(ponxside) == vecgeom::EInside::kOutside);
  assert(b2.Inside(ponyside) == vecgeom::EInside::kOutside);
  assert(b2.Inside(ponzside) == vecgeom::EInside::kOutside);
  assert(b2.Inside(Vec_t(10, 0, 0)) == vecgeom::EInside::kSurface);
  assert(b2.Inside(Vec_t(0, 10, 0)) == vecgeom::EInside::kSurface);
  assert(b2.Inside(Vec_t(0, 0, 10)) == vecgeom::EInside::kSurface);
  assert(b2.Inside(Vec_t(20, 20, 10)) == vecgeom::EInside::kOutside);
  assert(b2.Inside(Vec_t(100, 10, 30)) == vecgeom::EInside::kOutside);
  assert(b2.Inside(Vec_t(10, 20, 20)) == vecgeom::EInside::kOutside);

  // SafetyFromInside(P)
  Dist = b1.SafetyFromInside(pzero);
  assert(ApproxEqual(Dist, 20));
  Dist = b1.SafetyFromInside(vx);
  assert(ApproxEqual(Dist, 19));
  Dist = b1.SafetyFromInside(vy);
  assert(ApproxEqual(Dist, 20));
  Dist = b1.SafetyFromInside(vz);
  assert(ApproxEqual(Dist, 20));

  // Check DistanceToOut
  Dist = b1.DistanceToOut(pzero, vx, norm, convex);
  assert(ApproxEqual(Dist, 20));
  Dist = b1.DistanceToOut(pzero, vmx, norm, convex);
  assert(ApproxEqual(Dist, 20));
  Dist = b1.DistanceToOut(pzero, vy, norm, convex);
  assert(ApproxEqual(Dist, 30));
  Dist = b1.DistanceToOut(pzero, vmy, norm, convex);
  assert(ApproxEqual(Dist, 30));
  Dist = b1.DistanceToOut(pzero, vz, norm, convex);
  assert(ApproxEqual(Dist, 40));
  Dist = b1.DistanceToOut(pzero, vmz, norm, convex);
  assert(ApproxEqual(Dist, 40));
  Dist = b1.DistanceToOut(pzero, vxy, norm, convex);
  assert(ApproxEqual(Dist, std::sqrt(800.)));

  Dist = b1.DistanceToOut(ponxside, vx, norm, convex);
  assert(ApproxEqual(Dist, 0));
  Dist = b1.DistanceToOut(ponxside, vmx, norm, convex);
  assert(ApproxEqual(Dist, 40));
  Dist = b1.DistanceToOut(ponmxside, vmx, norm, convex);
  assert(ApproxEqual(Dist, 0));
  Dist = b1.DistanceToOut(ponyside, vy, norm, convex);
  assert(ApproxEqual(Dist, 0));
  Dist = b1.DistanceToOut(ponmyside, vmy, norm, convex);
  assert(ApproxEqual(Dist, 0));
  Dist = b1.DistanceToOut(ponzside, vz, norm, convex);
  assert(ApproxEqual(Dist, 0));
  Dist = b1.DistanceToOut(ponmzside, vmz, norm, convex);
  assert(ApproxEqual(Dist, 0));

  // SafetyFromOutside(P)
  Dist = b1.SafetyFromOutside(pbigx);
  assert(ApproxEqual(Dist, 80));
  Dist = b1.SafetyFromOutside(pbigmx);
  assert(ApproxEqual(Dist, 80));
  Dist = b1.SafetyFromOutside(pbigy);
  assert(ApproxEqual(Dist, 70));
  Dist = b1.SafetyFromOutside(pbigmy);
  assert(ApproxEqual(Dist, 70));
  Dist = b1.SafetyFromOutside(pbigz);
  assert(ApproxEqual(Dist, 60));
  Dist = b1.SafetyFromOutside(pbigmz);
  assert(ApproxEqual(Dist, 60));

  // DistanceToIn(P,V)
  Dist = b1.DistanceToIn(pbigx, vmx);
  assert(ApproxEqual(Dist, 80));
  Dist = b1.DistanceToIn(pbigmx, vx);
  assert(ApproxEqual(Dist, 80));
  Dist = b1.DistanceToIn(pbigy, vmy);
  assert(ApproxEqual(Dist, 70));
  Dist = b1.DistanceToIn(pbigmy, vy);
  assert(ApproxEqual(Dist, 70));
  Dist = b1.DistanceToIn(pbigz, vmz);
  assert(ApproxEqual(Dist, 60));
  Dist = b1.DistanceToIn(pbigmz, vz);
  assert(ApproxEqual(Dist, 60));

  // testing a few special directions ( with special zero values )
  // important since we sometimes operate with sign bits etc.
  Dist = b1.DistanceToIn(Vec_t(0, 0, -50), Vec_t(-0, -0, 1));
  assert(ApproxEqual(Dist, 10.));
  Dist = b1.DistanceToIn(Vec_t(0, 0, -50), Vec_t(0, -0, 1));
  assert(ApproxEqual(Dist, 10.));
  Dist = b1.DistanceToIn(Vec_t(0, 0, 50), Vec_t(0, 0, -1));
  assert(ApproxEqual(Dist, 10.));
  Dist = b1.DistanceToIn(Vec_t(0, -40, 0), Vec_t(0, 1, 0));
  assert(ApproxEqual(Dist, 10.));
  Dist = b1.DistanceToIn(Vec_t(0, -40, 0), Vec_t(0, 1, -0));
  assert(ApproxEqual(Dist, 10.));
  Dist = b1.DistanceToIn(Vec_t(0, -40, 0), Vec_t(-0, 1, 0));
  assert(ApproxEqual(Dist, 10.));
  Dist = b1.DistanceToIn(Vec_t(-30, 0, -0), Vec_t(1, -0, 0));
  assert(ApproxEqual(Dist, 10.));
  Dist = b1.DistanceToIn(Vec_t(30, 0., 0.), Vec_t(-1, 0, -0));
  assert(ApproxEqual(Dist, 10.));
  Dist = b1.DistanceToIn(Vec_t(-30, 0, 0), Vec_t(1, -0, -0));
  assert(ApproxEqual(Dist, 10.));

  Dist                                = b1.DistanceToIn(pbigx, vxy);
  if (Dist >= UUtils::kInfinity) Dist = UUtils::Infinity();
  assert(ApproxEqual(Dist, UUtils::Infinity()));
  Dist                                = b1.DistanceToIn(pbigmx, vxy);
  if (Dist >= UUtils::kInfinity) Dist = UUtils::Infinity();
  assert(ApproxEqual(Dist, UUtils::Infinity()));

  Vec_t pJohnXZ(9, 0, 12);
  Dist                                = b2.DistanceToIn(pJohnXZ, vxmz);
  if (Dist >= UUtils::kInfinity) Dist = UUtils::Infinity();
  assert(ApproxEqual(Dist, UUtils::Infinity()));

  Vec_t pJohnXY(12, 9, 0);
  Dist                                = b2.DistanceToIn(pJohnXY, vmxy);
  if (Dist >= UUtils::kInfinity) Dist = UUtils::Infinity();
  assert(ApproxEqual(Dist, UUtils::Infinity()));

  Dist = b2.DistanceToIn(pJohnXY, vmx);
  assert(ApproxEqual(Dist, 2));

  Vec_t pMyXY(32, -11, 0);
  Dist                                = b2.DistanceToIn(pMyXY, vmxy);
  if (Dist >= UUtils::kInfinity) Dist = UUtils::Infinity();
  assert(ApproxEqual(Dist, UUtils::Infinity()));

  Dist                                = b1.DistanceToIn(Vec_t(-25, -35, 0), vx);
  if (Dist >= UUtils::kInfinity) Dist = UUtils::Infinity();
  assert(ApproxEqual(Dist, UUtils::Infinity()));

  Dist                                = b1.DistanceToIn(Vec_t(-25, -35, 0), vy);
  if (Dist >= UUtils::kInfinity) Dist = UUtils::Infinity();
  assert(ApproxEqual(Dist, UUtils::Infinity()));

  Dist = b2.DistanceToIn(pJohnXY, vmx);
  assert(ApproxEqual(Dist, 2));

  Vec_t tempDir = Vec_t(-0.76165597579890043, 0.64364445891356026, -0.074515708658524193);
  Dist = box3.DistanceToIn(Vec_t(0.15000000000000185, -22.048743592955137, 2.4268539333219472), tempDir.Unit());

  assert(ApproxEqual(Dist, 0.0));

  /** testing tolerance of DistanceToIn **/
  Box_t b4("Box4", 5., 5., 5.);
  // a point very slightly inside should return 0
  tempDir = Vec_t(0.76315134679548990437, 0.53698876104646497964, -0.35950395323836459305);
  Dist = b4.DistanceToIn(Vec_t(-3.0087437277453119577, -4.9999999999999928946, 4.8935648380409944025), tempDir.Unit());
  assert(Dist <= 0.0);

  // a point on the surface pointing outside must return infinity length
  {
    auto d = b2.DistanceToIn(Vec_t(10, 0, 0), Vec_t(1, 0, 0));
    if (testvecgeom) {
      assert(d == vecgeom::InfinityLength<vecgeom::Precision>());
    } else {
      assert(d >= UUtils::kInfinity);
    }
  }
  {
    auto d = b2.DistanceToIn(Vec_t(10 - 0.9 * vecgeom::kHalfTolerance, 0, 0), Vec_t(1, 0, 0));
    assert(b2.Inside(Vec_t(10 - 0.9 * vecgeom::kHalfTolerance, 0, 0)) == vecgeom::EInside::kSurface);
    if (testvecgeom) {
      assert(d == vecgeom::InfinityLength<vecgeom::Precision>());
    } else {
      assert(d >= UUtils::kInfinity);
    }
  }
  /* **********************************************************
    */ /////////////////////////////////////////////////////

  bool ok = Test_VECGEOM_431<Box_t>();

  return ok;
}

template <class Box_t>
bool Test_VECGEOM_431()
{
  Vec_t vx(1, 0, 0), vy(0, 1, 0), vz(0, 0, 1);
  Vec_t norm;
  bool convex;
  double Dist;

  //=== add a couple test cases related to VECGEOM-431
  Box_t bx("Test Box #x", 1, 200, 200);
  // slightly outside of +x face and entering
  Vec_t testp = Vec_t(1, 0, 0) + 6.e-15 * vx;
  Vec_t testv = Vec_t(0, 1, 1) - 4.e-06 * vx;
  testv.Normalize();
  Dist = bx.DistanceToIn(testp, testv);
  assert(ApproxEqual(Dist, 0.0));
  Dist = bx.DistanceToOut(testp, testv, norm, convex);
  std::cout << "Line " << __LINE__ << ": point=" << testp << ", dir=" << testv << ", distOut=" << Dist
            << ", norm=" << norm << ", conv=" << convex << "\n";
  assert(ApproxEqual(Dist, 200. * sqrt(2.0)));
  assert(ApproxEqual(norm, (vy + vz).Normalized()));

  // and exiting
  Dist = bx.DistanceToOut(testp, vx, norm, convex);
  assert(ApproxEqual(Dist, 0.));
  assert(ApproxEqual(norm, vx));

  // slightly outside of -x face and entering
  testp = Vec_t(-1, 0, 0) - 6.e-15 * vx;
  testv = Vec_t(0, 1, 1) + 4.e-06 * vx;
  testv.Normalize();
  Dist = bx.DistanceToIn(testp, testv);
  assert(ApproxEqual(Dist, 0.0));
  Dist = bx.DistanceToOut(testp, testv, norm, convex);
  std::cout << "Line " << __LINE__ << ": point=" << testp << ", dir=" << testv << ", distOut=" << Dist
            << ", norm=" << norm << ", conv=" << convex << "\n";
  assert(ApproxEqual(Dist, 200. * sqrt(2.0)));
  assert(ApproxEqual(norm, (vy + vz).Normalized()));

  // and exiting
  Dist = bx.DistanceToOut(testp, -vx, norm, convex);
  assert(ApproxEqual(Dist, 0.));
  assert(ApproxEqual(norm, -vx));

  // slightly outside of +y face and entering
  Box_t by("Test Box #y", 200, 1, 200);
  testp = Vec_t(0, 1, 0) + 6.e-15 * vy;
  testv = Vec_t(1, 0, 1) - 4.e-6 * vy;
  testv.Normalize();
  Dist = by.DistanceToIn(testp, testv);
  assert(ApproxEqual(Dist, 0.0));
  Dist = by.DistanceToOut(testp, testv, norm, convex);
  std::cout << "Line " << __LINE__ << ": point=" << testp << ", dir=" << testv << ", distOut=" << Dist
            << ", norm=" << norm << ", conv=" << convex << "\n";
  assert(ApproxEqual(Dist, 200. * sqrt(2.0)));
  assert(ApproxEqual(norm, (vx + vz).Normalized()));

  // and exiting
  Dist = by.DistanceToOut(testp, vy, norm, convex);
  assert(ApproxEqual(Dist, 0.));
  assert(ApproxEqual(norm, vy));

  // slightly outside of -y face and entering
  testp = Vec_t(0, -1, 0) - 6.e-15 * vy;
  testv = Vec_t(1, 0, 1) + 4.e-6 * vy;
  testv.Normalize();
  Dist = by.DistanceToIn(testp, testv);
  assert(ApproxEqual(Dist, 0.0));
  Dist = by.DistanceToOut(testp, testv, norm, convex);
  assert(ApproxEqual(Dist, 200. * sqrt(2.0)));
  assert(ApproxEqual(norm, (vx + vz).Normalized()));

  // and exiting
  Dist = by.DistanceToOut(testp, -vy, norm, convex);
  assert(ApproxEqual(Dist, 0.));
  assert(ApproxEqual(norm, -vy));

  // slightly outside of +z face and entering
  Box_t bz("Test Box #z", 200, 200, 1);
  testp = Vec_t(0, 0, 1) + 6.e-15 * vz;
  testv = Vec_t(1, 1, 0) - 4.e-6 * vz;
  testv.Normalize();
  Dist = bz.DistanceToIn(testp, testv);
  assert(ApproxEqual(Dist, 0.0));
  Dist = bz.DistanceToOut(testp, testv, norm, convex);
  assert(ApproxEqual(Dist, 200 * sqrt(2.0)));
  assert(ApproxEqual(norm, (vx + vy).Normalized()));

  // and exiting
  Dist = bz.DistanceToOut(testp, vz, norm, convex);
  assert(ApproxEqual(Dist, 0.));
  assert(ApproxEqual(norm, vz));

  // slightly outside of -z face and entering
  testp = Vec_t(0, 0, -1) - 6.e-15 * vz;
  testv = Vec_t(1, 1, 0) + 4.e-6 * vz;
  testv.Normalize();
  Dist = bz.DistanceToIn(testp, testv);
  assert(ApproxEqual(Dist, 0.0));
  Dist = bz.DistanceToOut(testp, testv, norm, convex);
  assert(ApproxEqual(Dist, 200 * sqrt(2.0)));
  assert(ApproxEqual(norm, (vx + vy).Normalized()));

  // and exiting
  Dist = bz.DistanceToOut(testp, -vz, norm, convex);
  assert(ApproxEqual(Dist, 0.));
  assert(ApproxEqual(norm, -vz));

  //=== slightly inside of +x face
  testp = Vec_t(1, 0, 0) - 6.e-15 * vx;
  testv = Vec_t(0, 1, 1) + 4.e-6 * vx;
  testv.Normalize();
  Dist = bx.DistanceToIn(testp, testv);
  assert(ApproxEqual(Dist, UUtils::kInfinity));
  Dist = bx.DistanceToOut(testp, testv, norm, convex);
  assert(ApproxEqual(Dist, 0.0));

  // slightly inside of -x face
  testp = Vec_t(-1, 0, 0) + 6.e-15 * vx;
  testv = Vec_t(0, 1, 1) - 4.e-6 * vx;
  testv.Normalize();
  Dist = bx.DistanceToIn(testp, testv);
  assert(ApproxEqual(Dist, UUtils::kInfinity));
  Dist = bx.DistanceToOut(testp, testv, norm, convex);
  assert(ApproxEqual(Dist, 0.0));

  // slightly outside of +x face and exiting
  Box_t cube("cube", 10, 10, 10);
  testp = Vec_t(10, 0, 0) + 4.e-10 * vx;
  Dist  = cube.DistanceToOut(testp, vx, norm, convex);
  assert(ApproxEqual(Dist, 0.));
  assert(ApproxEqual(norm, vx.Normalized()));

  return true;
}

int main(int argc, char *argv[])
{
// enabling FPE exception
#if defined(__GNUCC__) && !defined(__CLANG__)
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
#endif

  testvecgeom = true;
  assert(TestBox<vecgeom::SimpleBox>());
  std::cout << "VecGeomBox passed\n";
  return 0;
}
