//
// File: TestTessellated.cpp
//
//    Ensure asserts are compiled in
//

//.. ensure asserts are compiled in
#undef NDEBUG

#include "base/Vector3D.h"
#include "volumes/Tessellated.h"
#include "ApproxEqual.h"
#ifdef VECGEOM_USOLIDS
#include "UTessellatedSolid.hh"
#include "UVector3.hh"
#endif
#include <cmath>

bool testvecgeom = true;

template <class Vec_t = vecgeom::Vector3D<vecgeom::Precision>>
vecgeom::SimpleTessellated *CreateTrdLikeTessellated(const char *name, double x1, double x2, double y1, double y2,
                                                     double z)
{
  // Create a tessellated solid from Trd parameters
  vecgeom::SimpleTessellated *stsl  = new vecgeom::SimpleTessellated(name);
  vecgeom::UnplacedTessellated *tsl = (vecgeom::UnplacedTessellated *)stsl->GetUnplacedVolume();
  // Top facet
  tsl->AddQuadrilateralFacet(Vec_t(-x2, y2, z), Vec_t(-x2, -y2, z), Vec_t(x2, -y2, z), Vec_t(x2, y2, z));
  // Bottom facet
  tsl->AddQuadrilateralFacet(Vec_t(-x1, y1, -z), Vec_t(x1, y1, -z), Vec_t(x1, -y1, -z), Vec_t(-x1, -y1, -z));
  // Front facet
  tsl->AddQuadrilateralFacet(Vec_t(-x2, -y2, z), Vec_t(-x1, -y1, -z), Vec_t(x1, -y1, -z), Vec_t(x2, -y2, z));
  // Right facet
  tsl->AddQuadrilateralFacet(Vec_t(x2, -y2, z), Vec_t(x1, -y1, -z), Vec_t(x1, y1, -z), Vec_t(x2, y2, z));
  // Behind facet
  tsl->AddQuadrilateralFacet(Vec_t(x2, y2, z), Vec_t(x1, y1, -z), Vec_t(-x1, y1, -z), Vec_t(-x2, y2, z));
  // Left facet
  tsl->AddQuadrilateralFacet(Vec_t(-x2, y2, z), Vec_t(-x1, y1, -z), Vec_t(-x1, -y1, -z), Vec_t(-x2, -y2, z));
  tsl->Close();
  return stsl;
}

template <typename Constants, class Tessellated_t, class Vec_t = vecgeom::Vector3D<vecgeom::Precision>>
bool TestTessellated()
{
  VUSolid::EnumInside inside;
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

  double Dist, vol, volCheck;
  Vec_t normal, norm;
  bool valid, convex;

  vecgeom::SimpleTessellated &tsl1 = *CreateTrdLikeTessellated<Vec_t>("Test Box #1", 20, 20, 30, 30, 40);
  vecgeom::SimpleTessellated &tsl2 = *CreateTrdLikeTessellated<Vec_t>("Test Trd", 10, 30, 20, 40, 40);
  vecgeom::SimpleTessellated &tsl3 =
      *CreateTrdLikeTessellated<Vec_t>("BABAR Trd", 0.14999999999999999, 0.14999999999999999, 24.707000000000001,
                                       24.707000000000001, 22.699999999999999);

  // check Cubic volume

  vol      = tsl1.Capacity();
  volCheck = 8 * 20 * 30 * 40;
  assert(ApproxEqual(vol, volCheck));

  // Check Surface area

  // std::cout<<"Trd Surface Area : " << tsl1.SurfaceArea()<<std::endl;
  assert(tsl1.SurfaceArea() == 20800);

  // Check Inside

  assert(tsl1.Inside(pzero) == vecgeom::EInside::kInside);
  assert(tsl1.Inside(pbigz) == vecgeom::EInside::kOutside);
  assert(tsl1.Inside(ponxside) == vecgeom::EInside::kSurface);
  assert(tsl1.Inside(ponyside) == vecgeom::EInside::kSurface);
  assert(tsl1.Inside(ponzside) == vecgeom::EInside::kSurface);

  inside = tsl1.Inside(Vec_t(20, 30, 40));
  //  std::cout << "tsl1.Inside((20,30,40)) = " << OutputInside(inside) << std::endl ;
  assert(inside == vecgeom::EInside::kSurface);

  inside = tsl1.Inside(Vec_t(-20, 30, 40));
  // std::cout << "tsl1.Inside((-20,30,40)) = " << OutputInside(inside) << std::endl ;
  assert(inside == vecgeom::EInside::kSurface);

  inside = tsl1.Inside(Vec_t(20, -30, 40));
  //  std::cout << "tsl1.Inside((20,-30,40)) = " << OutputInside(inside) << std::endl ;
  assert(inside == vecgeom::EInside::kSurface);

  inside = tsl1.Inside(Vec_t(20, 30, -40));
  // std::cout << "tsl1.Inside((20,30,-40)) = " << OutputInside(inside) << std::endl ;
  assert(inside == vecgeom::EInside::kSurface);

  inside = tsl1.Inside(Vec_t(20, 30, 0));
  // std::cout << "tsl1.Inside((20,30,0)) = " << OutputInside(inside) << std::endl ;
  assert(inside == vecgeom::EInside::kSurface);

  inside = tsl1.Inside(Vec_t(0, 30, 40));
  // std::cout << "tsl1.Inside((0,30,40)) = " << OutputInside(inside) << std::endl ;
  assert(inside == vecgeom::EInside::kSurface);

  inside = tsl1.Inside(Vec_t(20, 0, 40));
  // std::cout << "tsl1.Inside((20,0,40)) = " << OutputInside(inside) << std::endl ;
  assert(inside == vecgeom::EInside::kSurface);

  inside = tsl1.Inside(Vec_t(-20, -30, -40));
  // std::cout << "tsl1.Inside((-20,-30,-40)) = " << OutputInside(inside) << std::endl ;
  assert(inside == vecgeom::EInside::kSurface);

  assert(tsl2.Inside(pzero) == vecgeom::EInside::kInside);
  assert(tsl2.Inside(pbigz) == vecgeom::EInside::kOutside);
  assert(tsl2.Inside(ponxside) == vecgeom::EInside::kSurface);
  assert(tsl2.Inside(ponyside) == vecgeom::EInside::kSurface);
  assert(tsl2.Inside(ponzside) == vecgeom::EInside::kSurface);

  // Check Surface Normal

  (void)valid;
  valid = tsl1.Normal(ponxside, normal);
  assert(ApproxEqual(normal, Vec_t(1, 0, 0)));
  valid = tsl1.Normal(ponmxside, normal);
  assert(ApproxEqual(normal, Vec_t(-1, 0, 0)));
  valid = tsl1.Normal(ponyside, normal);
  assert(ApproxEqual(normal, Vec_t(0, 1, 0)));
  valid = tsl1.Normal(ponmyside, normal);
  assert(ApproxEqual(normal, Vec_t(0, -1, 0)));
  valid = tsl1.Normal(ponzside, normal);
  assert(ApproxEqual(normal, Vec_t(0, 0, 1)));
  valid = tsl1.Normal(ponmzside, normal);
  assert(ApproxEqual(normal, Vec_t(0, 0, -1)));
  valid = tsl1.Normal(ponzsidey, normal);
  assert(ApproxEqual(normal, Vec_t(0, 0, 1)));
  valid = tsl1.Normal(ponmzsidey, normal);
  assert(ApproxEqual(normal, Vec_t(0, 0, -1)));

  // Normals on Edges
  double cosa = 4 / std::sqrt(17.), sina = 1 / std::sqrt(17.);

  // SafetyFromInside(P)

  Dist = tsl1.SafetyFromInside(pzero);
  assert(ApproxEqual(Dist, 20));
  Dist = tsl1.SafetyFromInside(vx);
  assert(ApproxEqual(Dist, 19));
  Dist = tsl1.SafetyFromInside(vy);
  assert(ApproxEqual(Dist, 20));
  Dist = tsl1.SafetyFromInside(vz);
  assert(ApproxEqual(Dist, 20));

  Dist = tsl2.SafetyFromInside(pzero);
  assert(ApproxEqual(Dist, 20 * cosa));
  Dist = tsl2.SafetyFromInside(vx);
  assert(ApproxEqual(Dist, 19 * cosa));
  Dist = tsl2.SafetyFromInside(vy);
  assert(ApproxEqual(Dist, 20 * cosa));
  Dist = tsl2.SafetyFromInside(vz);
  assert(ApproxEqual(Dist, 20 * cosa + sina));

  // DistanceToOut(P,V)

  Dist = tsl1.DistanceToOut(pzero, vx, norm, convex);
  assert(ApproxEqual(Dist, 20) && ApproxEqual(norm, vx));
  if (!testvecgeom) assert(convex);
  Dist = tsl1.DistanceToOut(pzero, vmx, norm, convex);
  assert(ApproxEqual(Dist, 20) && ApproxEqual(norm, vmx));
  if (!testvecgeom) assert(convex);
  Dist = tsl1.DistanceToOut(pzero, vy, norm, convex);
  assert(ApproxEqual(Dist, 30) && ApproxEqual(norm, vy));
  if (!testvecgeom) assert(convex);
  Dist = tsl1.DistanceToOut(pzero, vmy, norm, convex);
  assert(ApproxEqual(Dist, 30) && ApproxEqual(norm, vmy));
  if (!testvecgeom) assert(convex);
  Dist = tsl1.DistanceToOut(pzero, vz, norm, convex);
  assert(ApproxEqual(Dist, 40) && ApproxEqual(norm, vz));
  if (!testvecgeom) assert(convex);
  Dist = tsl1.DistanceToOut(pzero, vmz, norm, convex);
  assert(ApproxEqual(Dist, 40) && ApproxEqual(norm, vmz));
  if (!testvecgeom) assert(convex);
  Dist = tsl1.DistanceToOut(pzero, vxy, norm, convex);
  assert(ApproxEqual(Dist, std::sqrt(800.)));
  if (!testvecgeom) assert(convex);

  Dist = tsl1.DistanceToOut(ponxside, vx, norm, convex);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(norm, vx));
  if (!testvecgeom) assert(convex);
  Dist = tsl1.DistanceToOut(ponmxside, vmx, norm, convex);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(norm, vmx));
  if (!testvecgeom) assert(convex);
  Dist = tsl1.DistanceToOut(ponyside, vy, norm, convex);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(norm, vy));
  if (!testvecgeom) assert(convex);
  Dist = tsl1.DistanceToOut(ponmyside, vmy, norm, convex);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(norm, vmy));
  if (!testvecgeom) assert(convex);
  Dist = tsl1.DistanceToOut(ponzside, vz, norm, convex);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(norm, vz));
  if (!testvecgeom) assert(convex);
  Dist = tsl1.DistanceToOut(ponmzside, vmz, norm, convex);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(norm, vmz));
  if (!testvecgeom) assert(convex);

  Dist = tsl2.DistanceToOut(pzero, vx, norm, convex);
  assert(ApproxEqual(Dist, 20) && ApproxEqual(norm, Vec_t(cosa, 0, -sina)));
  if (!testvecgeom) assert(convex);
  Dist = tsl2.DistanceToOut(pzero, vmx, norm, convex);
  assert(ApproxEqual(Dist, 20) && ApproxEqual(norm, Vec_t(-cosa, 0, -sina)));
  if (!testvecgeom) assert(convex);
  Dist = tsl2.DistanceToOut(pzero, vy, norm, convex);
  assert(ApproxEqual(Dist, 30) && ApproxEqual(norm, Vec_t(0, cosa, -sina)));
  if (!testvecgeom) assert(convex);
  Dist = tsl2.DistanceToOut(pzero, vmy, norm, convex);
  assert(ApproxEqual(Dist, 30) && ApproxEqual(norm, Vec_t(0, -cosa, -sina)));
  if (!testvecgeom) assert(convex);
  Dist = tsl2.DistanceToOut(pzero, vz, norm, convex);
  assert(ApproxEqual(Dist, 40) && ApproxEqual(norm, vz));
  if (!testvecgeom) assert(convex);
  Dist = tsl2.DistanceToOut(pzero, vmz, norm, convex);
  assert(ApproxEqual(Dist, 40) && ApproxEqual(norm, vmz));
  if (!testvecgeom) assert(convex);
  Dist = tsl2.DistanceToOut(pzero, vxy, norm, convex);
  assert(ApproxEqual(Dist, std::sqrt(800.)));
  if (!testvecgeom) assert(convex);

  Dist = tsl2.DistanceToOut(ponxside, vx, norm, convex);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(norm, Vec_t(cosa, 0, -sina)));
  if (!testvecgeom) assert(convex);
  Dist = tsl2.DistanceToOut(ponmxside, vmx, norm, convex);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(norm, Vec_t(-cosa, 0, -sina)));
  if (!testvecgeom) assert(convex);
  Dist = tsl2.DistanceToOut(ponyside, vy, norm, convex);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(norm, Vec_t(0, cosa, -sina)));
  if (!testvecgeom) assert(convex);
  Dist = tsl2.DistanceToOut(ponmyside, vmy, norm, convex);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(norm, Vec_t(0, -cosa, -sina)));
  if (!testvecgeom) assert(convex);
  Dist = tsl2.DistanceToOut(ponzside, vz, norm, convex);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(norm, vz));
  if (!testvecgeom) assert(convex);
  Dist = tsl2.DistanceToOut(ponmzside, vmz, norm, convex);
  assert(ApproxEqual(Dist, 0) && ApproxEqual(norm, vmz));
  if (!testvecgeom) assert(convex);

  // SafetyFromOutside(P)

  Dist = tsl1.SafetyFromOutside(pbigx);
  assert(ApproxEqual(Dist, 80));
  Dist = tsl1.SafetyFromOutside(pbigmx);
  assert(ApproxEqual(Dist, 80));
  Dist = tsl1.SafetyFromOutside(pbigy);
  assert(ApproxEqual(Dist, 70));
  Dist = tsl1.SafetyFromOutside(pbigmy);
  assert(ApproxEqual(Dist, 70));
  Dist = tsl1.SafetyFromOutside(pbigz);
  assert(ApproxEqual(Dist, 60));
  Dist = tsl1.SafetyFromOutside(pbigmz);
  assert(ApproxEqual(Dist, 60));

  Dist = tsl2.SafetyFromOutside(pbigx);
  assert(ApproxEqual(Dist, 80 * cosa));
  Dist = tsl2.SafetyFromOutside(pbigmx);
  assert(ApproxEqual(Dist, 80 * cosa));
  Dist = tsl2.SafetyFromOutside(pbigy);
  assert(ApproxEqual(Dist, 70 * cosa));
  Dist = tsl2.SafetyFromOutside(pbigmy);
  assert(ApproxEqual(Dist, 70 * cosa));
  Dist = tsl2.SafetyFromOutside(pbigz);
  assert(ApproxEqual(Dist, 60));
  Dist = tsl2.SafetyFromOutside(pbigmz);
  assert(ApproxEqual(Dist, 60));

  // DistanceToIn(P,V)

  Dist = tsl1.DistanceToIn(pbigx, vmx);
  assert(ApproxEqual(Dist, 80));
  Dist = tsl1.DistanceToIn(pbigmx, vx);
  assert(ApproxEqual(Dist, 80));
  Dist = tsl1.DistanceToIn(pbigy, vmy);
  assert(ApproxEqual(Dist, 70));
  Dist = tsl1.DistanceToIn(pbigmy, vy);
  assert(ApproxEqual(Dist, 70));
  Dist = tsl1.DistanceToIn(pbigz, vmz);
  assert(ApproxEqual(Dist, 60));
  Dist = tsl1.DistanceToIn(pbigmz, vz);
  assert(ApproxEqual(Dist, 60));
  Dist = tsl1.DistanceToIn(pbigx, vxy);
  assert(ApproxEqual(Dist, Constants::kInfLength));
  Dist = tsl1.DistanceToIn(pbigmx, vxy);
  assert(ApproxEqual(Dist, Constants::kInfLength));

  Dist = tsl2.DistanceToIn(pbigx, vmx);
  assert(ApproxEqual(Dist, 80));
  Dist = tsl2.DistanceToIn(pbigmx, vx);
  assert(ApproxEqual(Dist, 80));
  Dist = tsl2.DistanceToIn(pbigy, vmy);
  assert(ApproxEqual(Dist, 70));
  Dist = tsl2.DistanceToIn(pbigmy, vy);
  assert(ApproxEqual(Dist, 70));
  Dist = tsl2.DistanceToIn(pbigz, vmz);
  assert(ApproxEqual(Dist, 60));
  Dist = tsl2.DistanceToIn(pbigmz, vz);
  assert(ApproxEqual(Dist, 60));
  Dist = tsl2.DistanceToIn(pbigx, vxy);
  assert(ApproxEqual(Dist, Constants::kInfLength));
  Dist = tsl2.DistanceToIn(pbigmx, vxy);
  assert(ApproxEqual(Dist, Constants::kInfLength));

  Dist = tsl3.DistanceToIn(Vec_t(0.15000000000000185, -22.048743592955137, 2.4268539333219472),
                           Vec_t(-0.76165597579890043, 0.64364445891356026, -0.074515708658524193).Unit());

  //    std::cout<<"BABAR trd distance = "<<Dist<<std::ensl ;
  assert(ApproxEqual(Dist, 0.0));

  // return-value = 2.4415531753644804e-15

  // CalculateExtent
  Vec_t minExtent, maxExtent;
  tsl1.Extent(minExtent, maxExtent);
  // std::cout<<" min="<<minExtent<<" max="<<maxExtent<<std::endl;
  assert(ApproxEqual(minExtent, Vec_t(-20, -30, -40)));
  assert(ApproxEqual(maxExtent, Vec_t(20, 30, 40)));
  tsl2.Extent(minExtent, maxExtent);
  // std::cout<<" min="<<minExtent<<" max="<<maxExtent<<std::endl;
  assert(ApproxEqual(minExtent, Vec_t(-30, -40, -40)));
  assert(ApproxEqual(maxExtent, Vec_t(30, 40, 40)));

  return true;
}

#ifdef VECGEOM_USOLIDS
struct USOLIDSCONSTANTS {
  static constexpr double kInfLength = DBL_MAX; // UUSolids::kInfLength;
};
#endif
struct VECGEOMCONSTANTS {
  static constexpr double kInfLength = vecgeom::kInfLength;
};

int main(int argc, char *argv[])
{

  /*
    if (argc < 2) {
      std::cerr << "need to give argument: --usolids or --vecgeom\n";
      return 1;
    }

    if (!strcmp(argv[1], "--usolids")) {
  #ifndef VECGEOM_USOLIDS
      std::cerr << "VECGEOM_USOLIDS was not defined\n";
      return 2;
  #else
  #ifndef VECGEOM_REPLACE_USOLIDS
      TestTessellated<USOLIDSCONSTANTS, UTessellatedSolid>();
      std::cout << "USolids Tessellated passed\n";
  #else
      testvecgeom = true;
      TestTessellated<VECGEOMCONSTANTS, UTessellatedSolid>();
      std::cout << "USolids --> VecGeom Tessellated passed\n";
  #endif
  #endif
  */
  testvecgeom = true;
  TestTessellated<VECGEOMCONSTANTS, vecgeom::SimpleTessellated>();
  std::cout << "VecGeom Tessellated passed\n";
  /*
    } else {
      std::cerr << "need to give argument: --usolids or --vecgeom\n";
      return 1;
    }
  */
  return 0;
}