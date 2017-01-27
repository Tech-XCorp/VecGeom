/// \file UnplacedPolyhedron.cpp
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#include "base/Global.h"
#include "volumes/UnplacedPolyhedron.h"
#include "volumes/PlacedPolyhedron.h"
#include "volumes/SpecializedPolyhedron.h"
#include "volumes/utilities/GenerationUtilities.h"
#include "management/VolumeFactory.h"

#include <cmath>
#include <memory>

namespace vecgeom {

inline namespace VECGEOM_IMPL_NAMESPACE {

using namespace vecgeom::Polyhedron;

UnplacedPolyhedron::UnplacedPolyhedron(const int sideCount, const int zPlaneCount, Precision const zPlanes[],
                                       Precision const rMin[], Precision const rMax[])
    : UnplacedPolyhedron(0., kTwoPi, sideCount, zPlaneCount, zPlanes, rMin, rMax)
{
  DetectConvexity();
}

VECGEOM_CUDA_HEADER_BOTH
UnplacedPolyhedron::UnplacedPolyhedron(Precision phiStart, Precision phiDelta, const int sideCount,
                                       const int zPlaneCount, Precision const zPlanes[], Precision const rMin[],
                                       Precision const rMax[])
    : fPoly(phiStart, phiDelta, sideCount, zPlaneCount, zPlanes, rMin, rMax)
{
  DetectConvexity();
}

UnplacedPolyhedron::UnplacedPolyhedron(Precision phiStart, Precision phiDelta, const int sideCount,
                                       const int zPlaneCount,
                                       Precision const r[], // 2*zPlaneCount elements
                                       Precision const z[]  // ditto
                                       )
    : fPoly(phiStart, phiDelta, sideCount, zPlaneCount, r, z)
{
  DetectConvexity();
}

VECGEOM_CUDA_HEADER_BOTH
int UnplacedPolyhedron::GetNQuadrilaterals() const
{
  int count = 0;
  for (int i = 0; i < GetZSegmentCount(); ++i) {
    // outer
    count += GetZSegment(i).outer.size();
    // inner
    count += GetZSegment(i).inner.size();
    // phi
    count += GetZSegment(i).phi.size();
  }
  return count;
}

template <TranslationCode transCodeT, RotationCode rotCodeT>
VECGEOM_CUDA_HEADER_DEVICE
VPlacedVolume *UnplacedPolyhedron::Create(LogicalVolume const *const logical_volume,
                                          Transformation3D const *const transformation,
// const TranslationCode trans_code, const RotationCode rot_code,
#ifdef VECGEOM_NVCC
                                          const int id,
#endif
                                          VPlacedVolume *const placement)
{

  //#ifndef VECGEOM_NO_SPECIALIZATION

  UnplacedPolyhedron const *unplaced = static_cast<UnplacedPolyhedron const *>(logical_volume->GetUnplacedVolume());

  EInnerRadii innerRadii = unplaced->HasInnerRadii() ? EInnerRadii::kTrue : EInnerRadii::kFalse;

  EPhiCutout phiCutout = unplaced->HasPhiCutout()
                             ? (unplaced->HasLargePhiCutout() ? EPhiCutout::kLarge : EPhiCutout::kTrue)
                             : EPhiCutout::kFalse;

#ifndef VECGEOM_NVCC
#define POLYHEDRON_CREATE_SPECIALIZATION(INNER, PHI)                                              \
  return CreateSpecializedWithPlacement<SpecializedPolyhedron<transCodeT, rotCodeT, INNER, PHI>>( \
      logical_volume, transformation, placement)
#else
#define POLYHEDRON_CREATE_SPECIALIZATION(INNER, PHI)                                              \
  return CreateSpecializedWithPlacement<SpecializedPolyhedron<transCodeT, rotCodeT, INNER, PHI>>( \
      logical_volume, transformation, id, placement)
#endif

  if (innerRadii == EInnerRadii::kTrue) {
    if (phiCutout == EPhiCutout::kFalse) POLYHEDRON_CREATE_SPECIALIZATION(EInnerRadii::kTrue, EPhiCutout::kFalse);
    if (phiCutout == EPhiCutout::kTrue) POLYHEDRON_CREATE_SPECIALIZATION(EInnerRadii::kTrue, EPhiCutout::kTrue);
    if (phiCutout == EPhiCutout::kLarge) POLYHEDRON_CREATE_SPECIALIZATION(EInnerRadii::kTrue, EPhiCutout::kLarge);

  } else {
    if (phiCutout == EPhiCutout::kFalse) POLYHEDRON_CREATE_SPECIALIZATION(EInnerRadii::kFalse, EPhiCutout::kFalse);
    if (phiCutout == EPhiCutout::kTrue) POLYHEDRON_CREATE_SPECIALIZATION(EInnerRadii::kFalse, EPhiCutout::kTrue);
    if (phiCutout == EPhiCutout::kLarge) POLYHEDRON_CREATE_SPECIALIZATION(EInnerRadii::kFalse, EPhiCutout::kLarge);
  }

  // Return value in case of NO_SPECIALIZATION
  if (placement) {
    new (placement) SpecializedPolyhedron<transCodeT, rotCodeT, Polyhedron::EInnerRadii::kGeneric,
#ifdef VECGEOM_NVCC
                                          Polyhedron::EPhiCutout::kGeneric>(logical_volume, transformation, id);
#else
                                          Polyhedron::EPhiCutout::kGeneric>(logical_volume, transformation);
#endif
    return placement;
  }

  return new SpecializedPolyhedron<translation::kGeneric, rotation::kGeneric, Polyhedron::EInnerRadii::kGeneric,
#ifdef VECGEOM_NVCC
                                   Polyhedron::EPhiCutout::kGeneric>(logical_volume, transformation, id);
#else
                                   Polyhedron::EPhiCutout::kGeneric>(logical_volume, transformation);
#endif

#undef POLYHEDRON_CREATE_SPECIALIZATION
}

VECGEOM_CUDA_HEADER_DEVICE
VPlacedVolume *UnplacedPolyhedron::SpecializedVolume(LogicalVolume const *const volume,
                                                     Transformation3D const *const transformation,
                                                     const TranslationCode trans_code, const RotationCode rot_code,
#ifdef VECGEOM_NVCC
                                                     const int id,
#endif
                                                     VPlacedVolume *const placement) const
{

  return VolumeFactory::CreateByTransformation<UnplacedPolyhedron>(volume, transformation, trans_code, rot_code,
#ifdef VECGEOM_NVCC
                                                                   id,
#endif
                                                                   placement);
}

#ifndef VECGEOM_NVCC
void UnplacedPolyhedron::Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const
{
  aMin               = kInfLength;
  aMax               = -kInfLength;
  Precision phiStart = fPoly.fPhiStart;
  Precision phiDelta = fPoly.fPhiDelta;
  Precision sidePhi  = phiDelta / fPoly.fSideCount;
  // Specified radii are to the sides, not to the corners. Change these values,
  // as corners and not sides are used to compute the extent
  Precision conv = 1. / cos(0.5 * sidePhi);
  Vector3D<Precision> crt;
  // Loop all vertices and update min/max
  for (int iphi = 0; iphi <= fPoly.fSideCount; ++iphi) {
    Precision phi  = phiStart + iphi * sidePhi;
    Precision corx = conv * cos(phi);
    Precision cory = conv * sin(phi);
    for (int zPlaneCount = 0; zPlaneCount < fPoly.fZPlanes.size(); ++zPlaneCount) {
      // Do Rmin
      crt.Set(fPoly.fRMin[zPlaneCount] * corx, fPoly.fRMin[zPlaneCount] * cory, fPoly.fZPlanes[zPlaneCount]);
      for (int i = 0; i < 3; ++i) {
        aMin[i] = Min(aMin[i], crt[i]);
        aMax[i] = Max(aMax[i], crt[i]);
      }
      // Do Rmax
      crt.Set(fPoly.fRMax[zPlaneCount] * corx, fPoly.fRMax[zPlaneCount] * cory, fPoly.fZPlanes[zPlaneCount]);
      for (int i = 0; i < 3; ++i) {
        aMin[i] = Min(aMin[i], crt[i]);
        aMax[i] = Max(aMax[i], crt[i]);
      }
    }
  }
}

VECGEOM_CUDA_HEADER_BOTH
Precision UnplacedPolyhedron::DistanceSquarePointToSegment(Vector3D<Precision> &v1, Vector3D<Precision> &v2,
                                                           const Vector3D<Precision> &p) const
{

  Precision p1_p2_squareLength = (v1 - v2).Mag2();
  Precision dotProduct         = (p - v1).Dot(v1 - v2) / p1_p2_squareLength;
  if (dotProduct < 0) {
    return (p - v1).Mag2();
  } else if (dotProduct <= 1) {
    Precision p_p1_squareLength = (p - v1).Mag2();
    return p_p1_squareLength - dotProduct * dotProduct * p1_p2_squareLength;
  } else {
    return (p - v2).Mag2();
  }
}

VECGEOM_CUDA_HEADER_BOTH
bool UnplacedPolyhedron::InsideTriangle(Vector3D<Precision> &v1, Vector3D<Precision> &v2, Vector3D<Precision> &v3,
                                        const Vector3D<Precision> &p) const
{
  Precision epsilon_square = 0.00000001;
  Vector3D<Precision> vec1 = p - v1;
  Vector3D<Precision> vec2 = p - v2;
  Vector3D<Precision> vec3 = p - v3;

  bool sameSide1 = vec1.Dot(vec2) >= 0.;
  bool sameSide2 = vec1.Dot(vec3) >= 0.;
  bool sameSide3 = vec2.Dot(vec3) >= 0.;
  sameSide1      = sameSide1 && sameSide2 && sameSide3;

  if (sameSide1) return sameSide1;

  // If sameSide1 is false, point can be on the Surface or Outside
  // Use sqr of distance in order to check if point is on the Surface

  if (DistanceSquarePointToSegment(v1, v2, p) <= epsilon_square) return true;
  if (DistanceSquarePointToSegment(v1, v3, p) <= epsilon_square) return true;
  if (DistanceSquarePointToSegment(v2, v3, p) <= epsilon_square) return true;

  return false;
}

VECGEOM_CUDA_HEADER_BOTH
Precision UnplacedPolyhedron::GetTriangleArea(Vector3D<Precision> const &v1, Vector3D<Precision> const &v2,
                                              Vector3D<Precision> const &v3) const
{
  Vector3D<Precision> vec1 = v1 - v2;
  Vector3D<Precision> vec2 = v1 - v3;
  return 0.5 * (vec1.Cross(vec2)).Mag();
}

VECGEOM_CUDA_HEADER_BOTH
Vector3D<Precision> UnplacedPolyhedron::GetPointOnTriangle(Vector3D<Precision> const &v1, Vector3D<Precision> const &v2,
                                                           Vector3D<Precision> const &v3) const
{
  Precision alpha          = RNG::Instance().uniform(0.0, 1.0);
  Precision beta           = RNG::Instance().uniform(0.0, 1.0);
  Precision lambda1        = alpha * beta;
  Precision lambda0        = alpha - lambda1;
  Vector3D<Precision> vec1 = v2 - v1;
  Vector3D<Precision> vec2 = v3 - v1;
  return v1 + lambda0 * vec1 + lambda1 * vec2;
  /* A documented (slower) alternative that guarantees generating the points
     uniformly on the triangle surface:
     // http://www.cs.princeton.edu/~funk/tog02.pdf (section 4.2, formula 1)
     Precision sqr_r1 = sqrt(RNG::Instance().uniform(0.0, 1.0));
     Precision r2 = RNG::Instance().uniform(0.0, 1.0);
     return ( (1.-sqr_r1)*v1 + sqr_r1*(1.-r2)*v2 + sqr_r1*r2*v3 );
  */
}

VECGEOM_CUDA_HEADER_BOTH
Precision UnplacedPolyhedron::SurfaceArea() const
{
  if (fPoly.fSurfaceArea == 0.) {
    signed int j;
    Precision totArea = 0., area, aTop = 0., aBottom = 0.;

    // Below we generate the areas relevant to our solid
    // We are starting with ZSegments(lateral parts)

    for (j = 0; j < GetZSegmentCount(); ++j) {

      area = GetZSegment(j).outer.GetQuadrilateralArea(0) * GetSideCount();
      totArea += area;

      if (GetZSegment(j).hasInnerRadius) {
        area = GetZSegment(j).inner.GetQuadrilateralArea(0) * GetSideCount();
        totArea += area;
      }

      if (HasPhiCutout()) {
        area = GetZSegment(j).phi.GetQuadrilateralArea(0) * 2.0;
        totArea += area;
      }
    }

    // Must include top and bottom areas
    //

    Vector3D<Precision> point1 = GetZSegment(0).outer.GetCorners()[0][0];
    Vector3D<Precision> point2 = GetZSegment(0).outer.GetCorners()[1][0];
    Vector3D<Precision> point3, point4;
    if (GetZSegment(0).hasInnerRadius) {
      point3 = GetZSegment(0).inner.GetCorners()[0][0];
      point4 = GetZSegment(0).inner.GetCorners()[1][0];
      aTop   = GetSideCount() * (GetTriangleArea(point1, point2, point3) + GetTriangleArea(point3, point4, point2));

    } else {
      point3.Set(0.0, 0.0, GetZSegment(0).outer.GetCorners()[0][0].z());
      aTop = GetSideCount() * (GetTriangleArea(point1, point2, point3));
    }

    totArea += aTop;

    point1 = GetZSegment(GetZSegmentCount() - 1).outer.GetCorners()[2][0];
    point2 = GetZSegment(GetZSegmentCount() - 1).outer.GetCorners()[3][0];

    if (GetZSegment(GetZSegmentCount() - 1).hasInnerRadius) {
      point3  = GetZSegment(GetZSegmentCount() - 1).inner.GetCorners()[2][0];
      point4  = GetZSegment(GetZSegmentCount() - 1).inner.GetCorners()[3][0];
      aBottom = GetSideCount() * (GetTriangleArea(point1, point2, point3) + GetTriangleArea(point3, point4, point2));
    } else {
      point3.Set(0.0, 0.0, GetZSegment(GetZSegmentCount() - 1).outer.GetCorners()[2][0].z());
      aBottom = GetSideCount() * GetTriangleArea(point1, point2, point3);
    }

    totArea += aBottom;
    fPoly.fSurfaceArea = totArea;
  }
  return fPoly.fSurfaceArea;
}

Vector3D<Precision> UnplacedPolyhedron::GetPointOnSurface() const
{
  int j, Flag = 0;
  Precision chose, rnd, totArea = 0., Achose1, Achose2, area, aTop = 0., aBottom = 0.;

  Vector3D<Precision> p0, p1, p2, p3, pReturn;
  std::vector<Precision> aVector1;
  std::vector<Precision> aVector2;
  std::vector<Precision> aVector3;

  // Below we generate the areas relevant to our solid
  // We are starting with ZSegments(lateral parts)

  for (j = 0; j < GetZSegmentCount(); ++j) {
    area = GetZSegment(j).outer.GetQuadrilateralArea(0);
    totArea += area * GetSideCount();
    aVector1.push_back(area);
    if (GetZSegment(j).hasInnerRadius) {
      area = GetZSegment(j).inner.GetQuadrilateralArea(0);
      totArea += area * GetSideCount();
      aVector2.push_back(area);
    } else {
      aVector2.push_back(0.0);
    }

    if (HasPhiCutout()) {
      area = GetZSegment(j).phi.GetQuadrilateralArea(0);
      totArea += area * 2;
      aVector3.push_back(area);
    } else {
      aVector3.push_back(0.0);
    }
  }

  // Must include top and bottom areas
  //
  Vector3D<Precision> point1 = GetZSegment(0).outer.GetCorners()[0][0];
  Vector3D<Precision> point2 = GetZSegment(0).outer.GetCorners()[1][0];
  Vector3D<Precision> point3, point4;
  if (GetZSegment(0).hasInnerRadius) {
    point3 = GetZSegment(0).inner.GetCorners()[0][0];
    point4 = GetZSegment(0).inner.GetCorners()[1][0];
    aTop   = GetSideCount() * (GetTriangleArea(point1, point2, point3) + GetTriangleArea(point3, point4, point2));
  } else {
    point3.Set(0.0, 0.0, GetZSegment(0).outer.GetCorners()[0][0].z());
    aTop = GetSideCount() * (GetTriangleArea(point1, point2, point3));
  }

  totArea += aTop;
  point1 = GetZSegment(GetZSegmentCount() - 1).outer.GetCorners()[2][0];
  point2 = GetZSegment(GetZSegmentCount() - 1).outer.GetCorners()[3][0];

  if (GetZSegment(GetZSegmentCount() - 1).hasInnerRadius) {
    point3  = GetZSegment(GetZSegmentCount() - 1).inner.GetCorners()[2][0];
    point4  = GetZSegment(GetZSegmentCount() - 1).inner.GetCorners()[3][0];
    aBottom = GetSideCount() * (GetTriangleArea(point1, point2, point3) + GetTriangleArea(point3, point4, point2));
  } else {
    point3.Set(0.0, 0.0, GetZSegment(GetZSegmentCount() - 1).outer.GetCorners()[2][0].z());
    aBottom = GetSideCount() * GetTriangleArea(point1, point2, point3);
  }

  totArea += aBottom;

  // Chose area and Create Point on Surface

  Achose1 = 0.;
  Achose2 = GetSideCount() * (aVector1[0] + aVector2[0]) + 2 * aVector3[0];

  chose = RNG::Instance().uniform(0.0, totArea);
  // Point on Top or Bottom
  if ((chose >= 0.) && (chose < aTop + aBottom)) {

    chose = RNG::Instance().uniform(0.0, aTop + aBottom);
    Flag  = int(RNG::Instance().uniform(0.0, GetSideCount()));
    if ((chose >= 0.) && chose < aTop) {
      point1 = GetZSegment(GetZSegmentCount() - 1).outer.GetCorners()[2][Flag];
      point2 = GetZSegment(GetZSegmentCount() - 1).outer.GetCorners()[3][Flag];
      // Avoid generating points on degenerated triangles
      if (GetZSegment(GetZSegmentCount() - 1).hasInnerRadius) {
        point3 = GetZSegment(GetZSegmentCount() - 1).inner.GetCorners()[2][Flag];
        point4 = GetZSegment(GetZSegmentCount() - 1).inner.GetCorners()[3][Flag];
        if ((point4 - point3).Mag2() < kTolerance || RNG::Instance().uniform(0.0, 1.0) < 0.5)
          pReturn = GetPointOnTriangle(point3, point1, point2);
        else
          pReturn = GetPointOnTriangle(point4, point3, point2);
      } else {
        point3.Set(0.0, 0.0, GetZSegment(0).outer.GetCorners()[2][Flag].z());
        pReturn = GetPointOnTriangle(point3, point1, point2);
      }
      return pReturn;
    } else {
      point1 = GetZSegment(0).outer.GetCorners()[0][Flag];
      point2 = GetZSegment(0).outer.GetCorners()[1][Flag];

      if (GetZSegment(0).hasInnerRadius) {
        point3 = GetZSegment(0).inner.GetCorners()[0][Flag];
        point4 = GetZSegment(0).inner.GetCorners()[1][Flag];
        // Avoid generating points on degenerated triangles
        if ((point4 - point3).Mag2() < kTolerance || RNG::Instance().uniform(0.0, 1.0) < 0.5)
          pReturn = GetPointOnTriangle(point3, point1, point2);
        else
          pReturn = GetPointOnTriangle(point4, point3, point2);
      } else {
        point3.Set(0.0, 0.0, GetZSegment(0).outer.GetCorners()[0][Flag].z());
        pReturn = GetPointOnTriangle(point3, point1, point2);
      }
      return pReturn;
    }

  } else // Point on Lateral segment or Phi segment
  {

    for (j = 0; j < GetZSegmentCount(); j++) {
      if (((chose >= Achose1) && (chose < Achose2)) || (j == GetZSegmentCount() - 1)) {
        Flag = j;
        break;
      }
      Achose1 += GetSideCount() * (aVector1[j] + aVector2[j]) + 2. * aVector3[j];
      Achose2 = Achose1 + GetSideCount() * (aVector1[j + 1] + aVector2[j + 1]) + 2. * aVector3[j + 1];
    }
  }

  // At this point we have chosen a subsection
  // between to adjacent plane cuts

  j = Flag;

  rnd     = int(RNG::Instance().uniform(0.0, GetSideCount()));
  totArea = GetSideCount() * (aVector1[j] + aVector2[j]) + 2. * aVector3[j];
  chose   = RNG::Instance().uniform(0., totArea);
  Vector3D<Precision> RandVec;
  if ((chose >= 0.) && (chose <= GetSideCount() * aVector1[j])) {
    RandVec = (GetZSegment(j)).outer.GetPointOnFace(rnd);
    return RandVec;
  } else if ((chose >= 0.) && (chose <= GetSideCount() * (aVector1[j] + aVector2[j]))) {
    return (GetZSegment(j)).inner.GetPointOnFace(rnd);

  } else // Point on Phi segment
  {
    rnd     = int(RNG::Instance().uniform(0.0, 1.999999));
    RandVec = (GetZSegment(j)).phi.GetPointOnFace(rnd);
    return RandVec;
  }

  return Vector3D<Precision>(0, 0, 0); // error
}

// TODO: this functions seems to be neglecting the phi cut !!
Precision UnplacedPolyhedron::Capacity() const
{
  if (fPoly.fCapacity == 0.) {
    // Formula for section : V=h(f+F+sqrt(f*F))/3;
    // Fand f-areas of surfaces on +/-dz
    // h-heigh

    // a helper lambda for the volume calculation ( per quadrilaterial )
    auto VolumeHelperFunc = [&](Vector3D<Precision> const &a, Vector3D<Precision> const &b,
                                Vector3D<Precision> const &c, Vector3D<Precision> const &d) {
      Precision dz         = std::fabs(a.z() - c.z());
      Precision bottomArea = GetTriangleArea(a, b, Vector3D<Precision>(0., 0., a.z()));
      Precision topArea    = GetTriangleArea(c, d, Vector3D<Precision>(0., 0., c.z()));
      return dz * (bottomArea + topArea + std::sqrt(topArea * bottomArea));
    };

    for (int j = 0; j < GetZSegmentCount(); ++j) {
      // need to protect against empty segments because it could be
      // that the polyhedron makes a jump at this segment count
      if (GetZSegment(j).outer.size() > 0) {
        auto outercorners     = GetZSegment(j).outer.GetCorners();
        Vector3D<Precision> a = outercorners[0][0];
        Vector3D<Precision> b = outercorners[1][0];
        Vector3D<Precision> c = outercorners[2][0];
        Vector3D<Precision> d = outercorners[3][0];
        Precision volume      = VolumeHelperFunc(a, b, c, d); // outer volume

        if (GetZSegment(j).hasInnerRadius) {
          if (GetZSegment(j).inner.size() > 0) {
            auto innercorners = GetZSegment(j).inner.GetCorners();
            a                 = innercorners[0][0];
            b                 = innercorners[1][0];
            c                 = innercorners[2][0];
            d                 = innercorners[3][0];
            volume -= VolumeHelperFunc(a, b, c, d); // subtract inner volume
          }
        }
        fPoly.fCapacity += volume;
      }
    }
    fPoly.fCapacity *= GetSideCount() * (1. / 3.);
  }
  return fPoly.fCapacity;
}

VECGEOM_CUDA_HEADER_BOTH
bool UnplacedPolyhedron::Normal(Vector3D<Precision> const &point, Vector3D<Precision> &normal) const
{
  // Compute normal vector to closest surface
  return (
      PolyhedronImplementation<Polyhedron::EInnerRadii::kGeneric, Polyhedron::EPhiCutout::kGeneric>::ScalarNormalKernel(
          fPoly, point, normal));
}

#endif // !VECGEOM_NVCC

VECGEOM_CUDA_HEADER_BOTH
void UnplacedPolyhedron::Print() const
{
  printf("UnplacedPolyhedron {%i sides, phi %f to %f, %i segments}", fPoly.fSideCount, GetPhiStart() * kRadToDeg,
         GetPhiEnd() * kRadToDeg, fPoly.fZSegments.size());
  printf("}");
}

VECGEOM_CUDA_HEADER_BOTH
void UnplacedPolyhedron::PrintSegments() const
{
  printf("Printing %i polyhedron segments: ", fPoly.fZSegments.size());
  for (int i = 0, iMax = fPoly.fZSegments.size(); i < iMax; ++i) {
    printf("  Outer: ");
    fPoly.fZSegments[i].outer.Print();
    printf("\n");
    if (fPoly.fHasPhiCutout) {
      printf("  Phi: ");
      fPoly.fZSegments[i].phi.Print();
      printf("\n");
    }
    if (fPoly.fZSegments[i].hasInnerRadius) {
      printf("  Inner: ");
      fPoly.fZSegments[i].inner.Print();
      printf("\n");
    }
  }
}

void UnplacedPolyhedron::Print(std::ostream &os) const
{
  int oldprc = os.precision(16);
  int Nz     = fPoly.fZPlanes.size();
  os << "-----------------------------------------------------------\n"
     << "     *** Dump for solid - polyhedron ***\n"
     << "     ===================================================\n"
     << " Parameters:\n"
     << " Phi start= " << fPoly.fPhiStart * vecgeom::kRadToDeg
     << " deg, Phi delta= " << fPoly.fPhiDelta * vecgeom::kRadToDeg << " deg\n"
     << "     Number of segments along phi: " << fPoly.fSideCount << "\n"
     << "     N = number of Z-sections: " << Nz << "\n"
     << "     N+1 z-coordinates (in cm):\n";

  for (int i = 0; i < Nz; ++i) {
    os << "       at Z=" << fPoly.fZPlanes[i] << "cm:"
       << " Rmin=" << fPoly.fRMin[i] << "cm,"
       << " Rmax=" << fPoly.fRMax[i] << "cm\n";
  }
  os << "-----------------------------------------------------------\n";
  os.precision(oldprc);
}

VECGEOM_CUDA_HEADER_BOTH
void UnplacedPolyhedron::DetectConvexity()
{
  // Default safe convexity value
  fGlobalConvexity = false;

  if (fPoly.fConvexityPossible) {
    if (fPoly.fEqualRmax &&
        (fPoly.fPhiDelta <= kPi ||
         fPoly.fPhiDelta ==
             kTwoPi)) // In this case, Polycone become solid Cylinder, No need to check anything else, 100% convex
      fGlobalConvexity = true;
    else {
      if (fPoly.fPhiDelta <= kPi || fPoly.fPhiDelta == kTwoPi) {
        fGlobalConvexity = fPoly.fContinuousInSlope;
      }
    }
  }
}

#ifdef VECGEOM_CUDA_INTERFACE

DevicePtr<cuda::VUnplacedVolume> UnplacedPolyhedron::CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const gpuPtr) const
{

  // idea: reconstruct defining arrays: copy them to GPU; then construct the UnplacedPolycon object from scratch
  // on the GPU

  DevicePtr<Precision> zPlanesGpu;
  zPlanesGpu.Allocate(fPoly.fZPlanes.size());
  zPlanesGpu.ToDevice(&fPoly.fZPlanes[0], fPoly.fZPlanes.size());

  DevicePtr<Precision> rminGpu;
  rminGpu.Allocate(fPoly.fZPlanes.size());
  rminGpu.ToDevice(&fPoly.fRMin[0], fPoly.fZPlanes.size());

  DevicePtr<Precision> rmaxGpu;
  rmaxGpu.Allocate(fPoly.fZPlanes.size());
  rmaxGpu.ToDevice(&fPoly.fRMax[0], fPoly.fZPlanes.size());

  DevicePtr<cuda::VUnplacedVolume> gpupolyhedra = CopyToGpuImpl<UnplacedPolyhedron>(
      gpuPtr, fPoly.fPhiStart, fPoly.fPhiDelta, fPoly.fSideCount, fPoly.fZPlanes.size(), zPlanesGpu, rminGpu, rmaxGpu);

  zPlanesGpu.Deallocate();
  rminGpu.Deallocate();
  rmaxGpu.Deallocate();

  CudaAssertError();
  return gpupolyhedra;
}

DevicePtr<cuda::VUnplacedVolume> UnplacedPolyhedron::CopyToGpu() const
{
  return CopyToGpuImpl<UnplacedPolyhedron>();
}

#endif

} // End impl namespace

#ifdef VECGEOM_NVCC

namespace cxx {

template size_t DevicePtr<cuda::UnplacedPolyhedron>::SizeOf();
template void DevicePtr<cuda::UnplacedPolyhedron>::Construct(Precision phiStart, Precision phiDelta, int sideCount,
                                                             int zPlaneCount, DevicePtr<Precision> zPlanes,
                                                             DevicePtr<Precision> rMin,
                                                             DevicePtr<Precision> rMax) const;

} // End cxx namespace

#endif

} // End namespace vecgeom
