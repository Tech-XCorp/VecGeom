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
    : UnplacedPolyhedron(0, 360, sideCount, zPlaneCount, zPlanes, rMin, rMax)
{
}

VECGEOM_CUDA_HEADER_BOTH
bool UnplacedPolyhedron::CheckContinuityInSlope(const double rOuter[], const double zPlane[], const unsigned int nz)
{

  Precision prevSlope = kInfinity;
  for (unsigned int j = 0; j < nz - 1; ++j) {
    if (zPlane[j + 1] == zPlane[j]) {
      if (rOuter[j + 1] != rOuter[j]) return false;
    } else {
      Precision currentSlope = (rOuter[j + 1] - rOuter[j]) / (zPlane[j + 1] - zPlane[j]);
      if (currentSlope > prevSlope) return false;
      prevSlope = currentSlope;
    }
  }
  return true;
}

VECGEOM_CUDA_HEADER_BOTH
UnplacedPolyhedron::UnplacedPolyhedron(Precision phiStart, Precision phiDelta, const int sideCount,
                                       const int zPlaneCount, Precision const zPlanes[], Precision const rMin[],
                                       Precision const rMax[])
    : fSideCount(sideCount), fHasInnerRadii(false), fHasPhiCutout(phiDelta < 360), fHasLargePhiCutout(phiDelta < 180),
      fPhiStart(phiStart), fPhiDelta(phiDelta), fZSegments(zPlaneCount - 1), fZPlanes(zPlaneCount), fRMin(zPlaneCount),
      fRMax(zPlaneCount), fPhiSections(sideCount + 1), fBoundingTube(0, 1, 1, 0, kTwoPi), fSurfaceArea(0.),
      fCapacity(0.), fContinuousInSlope(true), fConvexityPossible(true), fEqualRmax(true)
{

  // initialize polyhedron internals
  Initialize(phiStart, phiDelta, sideCount, zPlaneCount, zPlanes, rMin, rMax);
  DetectConvexity();
}

VECGEOM_CUDA_HEADER_BOTH
void UnplacedPolyhedron::Initialize(Precision phiStart, Precision phiDelta, const int sideCount, const int zPlaneCount,
                                    Precision const zPlanes[], Precision const rMin[], Precision const rMax[])
{
  typedef Vector3D<Precision> Vec_t;

  // Sanity check of input parameters
  assert(zPlaneCount > 1);
  assert(fSideCount > 0);

  copy(zPlanes, zPlanes + zPlaneCount, &fZPlanes[0]);
  copy(rMin, rMin + zPlaneCount, &fRMin[0]);
  copy(rMax, rMax + zPlaneCount, &fRMax[0]);

  double startRmax = rMax[0];
  for (int i = 0; i < zPlaneCount; i++) {
    fConvexityPossible &= (rMin[i] == 0.);
    fEqualRmax &= (startRmax == rMax[i]);
  }
  fContinuousInSlope = CheckContinuityInSlope(rMax, zPlanes, zPlaneCount);

  // Initialize segments
  // sometimes there will be no quadrilaterals: for instance when
  // rmin jumps at some z and rmax remains continouus
  for (int i = 0; i < zPlaneCount - 1; ++i) {
    // Z-planes must be monotonically increasing
    assert(zPlanes[i] <= zPlanes[i + 1]);

    fZSegments[i].hasInnerRadius = rMin[i] > 0 || rMin[i + 1] > 0;

    int multiplier = (zPlanes[i] == zPlanes[i + 1] && rMax[i] == rMax[i + 1]) ? 0 : 1;

    // create quadrilaterals in a predefined place with placement new
    new (&fZSegments[i].outer) Quadrilaterals(sideCount * multiplier);

    // no phi segment here if degenerate z;
    if (fHasPhiCutout) {
      multiplier = (zPlanes[i] == zPlanes[i + 1]) ? 0 : 1;
      new (&fZSegments[i].phi) Quadrilaterals(2 * multiplier);
    }

    multiplier = (zPlanes[i] == zPlanes[i + 1] && rMin[i] == rMin[i + 1]) ? 0 : 1;

    if (fZSegments[i].hasInnerRadius) {
      new (&fZSegments[i].inner) Quadrilaterals(sideCount * multiplier);
      fHasInnerRadii = true;
    }
  }

  // Compute the cylindrical coordinate phi along which the corners are placed
  assert(phiDelta > 0);
  phiStart = NormalizeAngle<kScalar>(kDegToRad * phiStart);
  phiDelta *= kDegToRad;
  if (phiDelta > kTwoPi) phiDelta = kTwoPi;
  Precision sidePhi               = phiDelta / sideCount;
  vecgeom::unique_ptr<Precision[]> vertixPhi(new Precision[sideCount + 1]);
  for (int i = 0, iMax = sideCount + 1; i < iMax; ++i) {
    vertixPhi[i]                     = NormalizeAngle<kScalar>(phiStart + i * sidePhi);
    Vector3D<Precision> cornerVector = Vec_t::FromCylindrical(1., vertixPhi[i], 0).Normalized().FixZeroes();
    fPhiSections.set(i, cornerVector.Normalized().Cross(Vector3D<Precision>(0, 0, -1)));
  }
  if (!fHasPhiCutout) {
    // If there is no phi cutout, last phi is equal to the first
    vertixPhi[sideCount] = vertixPhi[0];
  }

  // Specified radii are to the sides, not to the corners. Change these values,
  // as corners and not sides are used to build the structure
  Precision cosHalfDeltaPhi = cos(0.5 * sidePhi);
  Precision innerRadius = kInfinity, outerRadius = -kInfinity;
  for (int i = 0; i < zPlaneCount; ++i) {
    // Use distance to side for minimizing inner radius of bounding tube
    if (rMin[i] < innerRadius) innerRadius = rMin[i];
    // rMin[i] /= cosHalfDeltaPhi;
    // rMax[i] /= cosHalfDeltaPhi;
    assert(rMin[i] >= 0 && rMax[i] > 0);
    // Use distance to corner for minimizing outer radius of bounding tube
    if (rMax[i] > outerRadius) outerRadius = rMax[i];
  }
  // need to convert from distance to planes to real radius in case of outerradius
  // the inner radius of the bounding tube is given by min(rMin[])
  outerRadius /= cosHalfDeltaPhi;

  // Create bounding tube with biggest outer radius and smallest inner radius
  Precision boundingTubeZ = 0.5 * (zPlanes[zPlaneCount - 1] - zPlanes[0] + kTolerance);
  // Make bounding tube phi range a bit larger to contain all points on phi boundaries
  const Precision kPhiTolerance = 100 * kTolerance;
  // The increase in the angle has to be large enough to contain most of
  // kSurface points. There will be some points close to the Z axis which will
  // not be contained. The value is empirical to satisfy ShapeTester
  Precision boundsPhiStart = !fHasPhiCutout ? 0 : phiStart - kPhiTolerance;
  Precision boundsPhiDelta = !fHasPhiCutout ? kTwoPi : phiDelta + 2 * kPhiTolerance;
  // correct inner and outer Radius with conversion factor
  // innerRadius /= cosHalfDeltaPhi;
  // outerRadius /= cosHalfDeltaPhi;

  fBoundingTube = UnplacedTube(innerRadius - kHalfTolerance, outerRadius + kHalfTolerance, boundingTubeZ,
                               boundsPhiStart, boundsPhiDelta);
  // The offset has to match the middle of the polyhedron
  fBoundingTubeOffset = 0.5 * (zPlanes[0] + zPlanes[zPlaneCount - 1]);

  // Ease indexing into twodimensional vertix array
  auto VertixIndex = [&sideCount](int plane, int corner) { return plane * (sideCount + 1) + corner; };

  // Precompute all vertices to ensure that there are no numerical cracks in the
  // surface.
  const int nVertices = zPlaneCount * (sideCount + 1);
  vecgeom::unique_ptr<Vec_t[]> outerVertices(new Vec_t[nVertices]);
  vecgeom::unique_ptr<Vec_t[]> innerVertices(new Vec_t[nVertices]);
  for (int i = 0; i < zPlaneCount; ++i) {
    for (int j = 0, jMax = sideCount + fHasPhiCutout; j < jMax; ++j) {
      int index            = VertixIndex(i, j);
      outerVertices[index] = Vec_t::FromCylindrical(rMax[i] / cosHalfDeltaPhi, vertixPhi[j], zPlanes[i]).FixZeroes();
      innerVertices[index] = Vec_t::FromCylindrical(rMin[i] / cosHalfDeltaPhi, vertixPhi[j], zPlanes[i]).FixZeroes();
    }
    // Non phi cutout case
    if (!fHasPhiCutout) {
      // Make last vertices identical to the first phi coordinate
      outerVertices[VertixIndex(i, sideCount)] = outerVertices[VertixIndex(i, 0)];
      innerVertices[VertixIndex(i, sideCount)] = innerVertices[VertixIndex(i, 0)];
    }
  }

  // Build segments by drawing quadrilaterals between vertices
  for (int iPlane = 0; iPlane < zPlaneCount - 1; ++iPlane) {

    auto WrongNormal = [](Vector3D<Precision> const &normal, Vector3D<Precision> const &corner) {
      return normal[0] * corner[0] + normal[1] * corner[1] < 0;
    };

    // Draw the regular quadrilaterals along phi
    for (int iSide = 0; iSide < fZSegments[iPlane].outer.size(); ++iSide) {
      fZSegments[iPlane].outer.Set(
          iSide, outerVertices[VertixIndex(iPlane, iSide)], outerVertices[VertixIndex(iPlane, iSide + 1)],
          outerVertices[VertixIndex(iPlane + 1, iSide + 1)], outerVertices[VertixIndex(iPlane + 1, iSide)]);
      // Normal has to point away from Z-axis
      if (WrongNormal(fZSegments[iPlane].outer.GetNormal(iSide), outerVertices[VertixIndex(iPlane, iSide)])) {
        fZSegments[iPlane].outer.FlipSign(iSide);
      }
    }
    for (int iSide = 0; iSide < fZSegments[iPlane].inner.size(); ++iSide) {
      if (fZSegments[iPlane].hasInnerRadius) {
        fZSegments[iPlane].inner.Set(
            iSide, innerVertices[VertixIndex(iPlane, iSide)], innerVertices[VertixIndex(iPlane, iSide + 1)],
            innerVertices[VertixIndex(iPlane + 1, iSide + 1)], innerVertices[VertixIndex(iPlane + 1, iSide)]);
        // Normal has to point away from Z-axis
        if (WrongNormal(fZSegments[iPlane].inner.GetNormal(iSide), innerVertices[VertixIndex(iPlane, iSide)])) {
          fZSegments[iPlane].inner.FlipSign(iSide);
        }
      }
    }

    if (fHasPhiCutout && fZSegments[iPlane].phi.size() == 2) {
      // If there's a phi cutout, draw two quadrilaterals connecting the four
      // corners (two inner, two outer) of the first and last phi coordinate,
      // respectively
      fZSegments[iPlane].phi.Set(0, innerVertices[VertixIndex(iPlane, 0)], innerVertices[VertixIndex(iPlane + 1, 0)],
                                 outerVertices[VertixIndex(iPlane + 1, 0)], outerVertices[VertixIndex(iPlane, 0)]);
      // Make sure normal points backwards along phi
      if (fZSegments[iPlane].phi.GetNormal(0).Dot(fPhiSections[0]) > 0) {
        fZSegments[iPlane].phi.FlipSign(0);
      }
      fZSegments[iPlane].phi.Set(
          1, outerVertices[VertixIndex(iPlane, sideCount)], outerVertices[VertixIndex(iPlane + 1, sideCount)],
          innerVertices[VertixIndex(iPlane + 1, sideCount)], innerVertices[VertixIndex(iPlane, sideCount)]);
      // Make sure normal points forwards along phi
      if (fZSegments[iPlane].phi.GetNormal(1).Dot(fPhiSections[fSideCount]) < 0) {
        fZSegments[iPlane].phi.FlipSign(1);
      }
    }

  } // End loop over segments
} // end constructor

UnplacedPolyhedron::UnplacedPolyhedron(Precision phiStart, Precision phiDelta, const int sideCount,
                                       const int zPlaneCount,
                                       Precision const r[], // 2*zPlaneCount elements
                                       Precision const z[]  // ditto
                                       )
    : fSideCount(sideCount), fHasInnerRadii(false), fHasPhiCutout(phiDelta < 360), fHasLargePhiCutout(phiDelta < 180),
      fPhiStart(phiStart), fPhiDelta(phiDelta), fZSegments(zPlaneCount - 1), fZPlanes(zPlaneCount), fRMin(zPlaneCount),
      fRMax(zPlaneCount), fPhiSections(sideCount + 1), fBoundingTube(0, 1, 1, 0, kTwoPi), fSurfaceArea(0.),
      fCapacity(0.)
{
  // data integrity checks
  for (int i = 0; i <= zPlaneCount; ++i) {
    assert(z[i] == z[2 * zPlaneCount - 1 - i] && "UnplPolyhedron ERROR: z[] array is not symmetrical, please fix.\n");
  }

  // reuse input array as argument, in ascending order
  int Nz                 = zPlaneCount;
  bool ascendingZ        = true;
  const Precision *zarg  = &z[0];
  const Precision *r1arg = r;
  if (z[0] > z[1]) {
    ascendingZ = false;
    zarg       = z + Nz; // second half of input z[] is ascending due to symmetry already verified
    r1arg      = r + Nz;
  }

  // reorganize remainder of r[] data in ascending-z order
  Precision r2arg[Nz];
  for (int i = 0; i < Nz; ++i)
    r2arg[i] = (ascendingZ ? r[2 * Nz - 1 - i] : r[Nz - 1 - i]);

  // identify which rXarg is rmax and rmin and ensure that Rmax > Rmin for all points provided
  const Precision *rmin = r1arg, *rmax = r2arg;
  if (r1arg[0] > r2arg[0]) {
    rmax = r1arg;
    rmin = r2arg;
  }

  // final data integrity cross-check
  for (int i = 0; i < Nz; ++i) {
    assert(rmax[i] > rmin[i] &&
           "UnplPolycone ERROR: r[] provided has problems of the Rmax < Rmin type, please check!\n");
  }

  // Delegate to full constructor
  Initialize(phiStart, phiDelta, sideCount, zPlaneCount, zarg, rmin, rmax);
}

// TODO: move this to HEADER; this is now stored as a member
VECGEOM_CUDA_HEADER_BOTH
Precision UnplacedPolyhedron::GetPhiStart() const
{
  return kRadToDeg * NormalizeAngle<kScalar>(fPhiSections[0].Cross(Vector3D<Precision>(0, 0, 1)).Phi());
}

// TODO: move this to HEADER; this is now stored as a member
VECGEOM_CUDA_HEADER_BOTH
Precision UnplacedPolyhedron::GetPhiEnd() const
{
  return !HasPhiCutout()
             ? 360
             : kRadToDeg *
                   NormalizeAngle<kScalar>(fPhiSections[GetSideCount()].Cross(Vector3D<Precision>(0, 0, 1)).Phi());
}

// TODO: move this to HEADER
VECGEOM_CUDA_HEADER_BOTH
Precision UnplacedPolyhedron::GetPhiDelta() const
{
  return fPhiDelta;
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

#if defined(VECGEOM_USOLIDS)
VECGEOM_CUDA_HEADER_BOTH
std::ostream &UnplacedPolyhedron::StreamInfo(std::ostream &os) const
{
  int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "     *** Dump for solid - polyhedron ***\n"
     << "     ===================================================\n"
     << " Solid type: " << GetEntityType() << "\n"
     << " Parameters:\n"
     << " Phi start=" << fPhiStart << "deg, Phi delta=" << fPhiDelta << "deg\n"
     << "     Number of segments along phi: " << fSideCount << "\n"
     << "     N = number of Z-sections: " << fZSegments.size() << "\n"
     << "     N+1 z-coordinates (in cm):\n";
  uint Nz = (uint)fZSegments.size();

  // check that we have the right number of planes
  assert((int)Nz == fZPlanes.size() - 1);

  for (uint i = 0; i < Nz; ++i) {
    os << "       at Z=" << fZPlanes[i] << "cm:"
       << " Rmin=" << fRMin[i] << "cm,"
       << " Rmax=" << fRMax[i] << "cm\n";
  }
  // print last plane info
  os << "       at Z=" << fZPlanes[Nz] << "cm:"
     << " Rmin=" << fRMin[Nz] << "cm,"
     << " Rmax=" << fRMax[Nz] << "cm\n";
  os << "-----------------------------------------------------------\n";
  os.precision(oldprc);
  return os;
}
#endif

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
  const UnplacedTube &bTube = GetBoundingTube();
  bTube.Extent(aMin, aMax);
  aMin.z() += fBoundingTubeOffset;
  aMax.z() += fBoundingTubeOffset;
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
  Precision fEpsilon_square = 0.00000001;
  Vector3D<Precision> vec1  = p - v1;
  Vector3D<Precision> vec2  = p - v2;
  Vector3D<Precision> vec3  = p - v3;

  bool sameSide1 = vec1.Dot(vec2) >= 0.;
  bool sameSide2 = vec1.Dot(vec3) >= 0.;
  bool sameSide3 = vec2.Dot(vec3) >= 0.;
  sameSide1      = sameSide1 && sameSide2 && sameSide3;

  if (sameSide1) return sameSide1;

  // If sameSide1 is false, point can be on the Surface or Outside
  // Use sqr of distance in order to check if point is on the Surface

  if (DistanceSquarePointToSegment(v1, v2, p) <= fEpsilon_square) return true;
  if (DistanceSquarePointToSegment(v1, v3, p) <= fEpsilon_square) return true;
  if (DistanceSquarePointToSegment(v2, v3, p) <= fEpsilon_square) return true;

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
  if (fSurfaceArea == 0.) {
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
    fSurfaceArea = totArea;
  }
  return fSurfaceArea;
}

VECGEOM_CUDA_HEADER_BOTH
Vector3D<Precision> UnplacedPolyhedron::GetPointOnSurface() const
{
  int j, numPlanes = GetZSegmentCount() + 1, Flag = 0;
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

    for (j = 0; j < numPlanes - 1; j++) {
      if (((chose >= Achose1) && (chose < Achose2)) || (j == numPlanes - 1)) {
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
  if (fCapacity == 0.) {
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
        fCapacity += volume;
      }
    }
    fCapacity *= GetSideCount() * (1. / 3.);
  }
  return fCapacity;
}

VECGEOM_CUDA_HEADER_BOTH
bool UnplacedPolyhedron::Normal(Vector3D<Precision> const &point, Vector3D<Precision> &normal) const
{
  // Compute normal vector to closest surface
  return (PolyhedronImplementation<translation::kGeneric, rotation::kGeneric, Polyhedron::EInnerRadii::kGeneric,
                                   Polyhedron::EPhiCutout::kGeneric>::ScalarNormalKernel(*this, point, normal));
}

#endif // !VECGEOM_NVCC

VECGEOM_CUDA_HEADER_BOTH
void UnplacedPolyhedron::Print() const
{
  printf("UnplacedPolyhedron {%i sides, phi %f to %f, %i segments}", fSideCount, GetPhiStart(), GetPhiEnd(),
         fZSegments.size());
  printf("}");
}

VECGEOM_CUDA_HEADER_BOTH
void UnplacedPolyhedron::PrintSegments() const
{
  printf("Printing %i polyhedron segments: ", fZSegments.size());
  for (int i = 0, iMax = fZSegments.size(); i < iMax; ++i) {
    printf("  Outer: ");
    fZSegments[i].outer.Print();
    printf("\n");
    if (fHasPhiCutout) {
      printf("  Phi: ");
      fZSegments[i].phi.Print();
      printf("\n");
    }
    if (fZSegments[i].hasInnerRadius) {
      printf("  Inner: ");
      fZSegments[i].inner.Print();
      printf("\n");
    }
  }
}

void UnplacedPolyhedron::Print(std::ostream &os) const
{
  os << "UnplacedPolyhedron {" << fSideCount << " sides, " << fZSegments.size() << " segments, "
     << ((fHasInnerRadii) ? "has inner radii" : "no inner radii") << "}";
}

VECGEOM_CUDA_HEADER_BOTH
void UnplacedPolyhedron::DetectConvexity()
{
  // Default safe convexity value
  fGlobalConvexity = false;

  if (fConvexityPossible) {
    if (fEqualRmax &&
        (fPhiDelta <= 180 ||
         fPhiDelta == 360)) // In this case, Polycone become solid Cylinder, No need to check anything else, 100% convex
      fGlobalConvexity = true;
    else {
      if (fPhiDelta <= 180 || fPhiDelta == 360) {
        fGlobalConvexity = fContinuousInSlope;
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
  zPlanesGpu.Allocate(fZPlanes.size());
  zPlanesGpu.ToDevice(&fZPlanes[0], fZPlanes.size());

  DevicePtr<Precision> rminGpu;
  rminGpu.Allocate(fZPlanes.size());
  rminGpu.ToDevice(&fRMin[0], fZPlanes.size());

  DevicePtr<Precision> rmaxGpu;
  rmaxGpu.Allocate(fZPlanes.size());
  rmaxGpu.ToDevice(&fRMax[0], fZPlanes.size());

  DevicePtr<cuda::VUnplacedVolume> gpupolyhedra = CopyToGpuImpl<UnplacedPolyhedron>(
      gpuPtr, fPhiStart, fPhiDelta, fSideCount, fZPlanes.size(), zPlanesGpu, rminGpu, rmaxGpu);

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
