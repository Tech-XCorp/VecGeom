//===-- kernel/TrapezoidImplementation.h ----------------------------*- C++ -*-===//
//===--------------------------------------------------------------------------===//
///
/// \file   kernel/TrapezoidImplementation.h
/// \author Guilherme Lima (lima@fnal.gov)
/// \brief This file implements the algorithms for the trapezoid
///
/// Implementation details: initially based on USolids algorithms and vectorized types.
///
//===--------------------------------------------------------------------------===//
///
/// 140520  G. Lima   Created from USolids' UTrap algorithms
/// 160722  G. Lima   Revision + moving to new backend structure

#ifndef VECGEOM_VOLUMES_KERNEL_TRAPEZOIDIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_TRAPEZOIDIMPLEMENTATION_H_

#include "base/Vector3D.h"
#include "volumes/TrapezoidStruct.h"
#include "volumes/kernel/GenericKernels.h"
#include <VecCore/VecCore>

#include <cstdio>

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(struct TrapezoidImplementation;);
VECGEOM_DEVICE_DECLARE_CONV(struct, TrapezoidImplementation);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedTrapezoid;
class UnplacedTrapezoid;

struct TrapezoidImplementation {

  using PlacedShape_t    = PlacedTrapezoid;
  using UnplacedStruct_t = TrapezoidStruct<double>;
  using UnplacedVolume_t = UnplacedTrapezoid;
#ifdef VECGEOM_PLANESHELL_DISABLE
  using TrapSidePlane = TrapezoidStruct<double>::TrapSidePlane;
#endif

  VECGEOM_CUDA_HEADER_BOTH
  static void PrintType()
  {
    // printf("SpecializedTrapezoid<%i, %i>", transCodeT, rotCodeT);
  }

  template <typename Stream>
  static void PrintType(Stream &st, int transCodeT = translation::kGeneric, int rotCodeT = rotation::kGeneric)
  {
    st << "SpecializedTrapezoid<" << transCodeT << "," << rotCodeT << ">";
  }

  template <typename Stream>
  static void PrintImplementationType(Stream &st)
  {
    (void)st;
    // st << "TrapezoidImplementation<" << transCodeT << "," << rotCodeT << ">";
  }

  template <typename Stream>
  static void PrintUnplacedType(Stream &st)
  {
    (void)st;
    // TODO: this is wrong
    st << "UnplacedTrapezoid";
  }

#ifdef VECGEOM_PLANESHELL_DISABLE
  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  static void EvaluateTrack(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point,
                            Vector3D<Real_v> const &dir, Real_v *pdist, Real_v *proj, Real_v *vdist)
  {
    TrapSidePlane const *fPlanes = unplaced.GetPlanes();
    // loop over side planes - find pdist,proj for each side plane
    // auto-vectorizable part of loop
    for (unsigned int i = 0; i < 4; ++i) {
      // Note: normal vector is pointing outside the volume (convention), therefore
      // pdist>0 if point is outside  and  pdist<0 means inside
      pdist[i] = fPlanes[i].fA * point.x() + fPlanes[i].fB * point.y() + fPlanes[i].fC * point.z() + fPlanes[i].fD;

      // proj is projection of dir over the normal vector of side plane, hence
      // proj > 0 if pointing ~same direction as normal and proj<0 if ~opposite to normal
      proj[i] = fPlanes[i].fA * dir.x() + fPlanes[i].fB * dir.y() + fPlanes[i].fC * dir.z();

      vdist[i] = -pdist[i] / NonZero(proj[i]);
    }
  }
#endif

  template <typename Real_v, typename Bool_v>
  VECGEOM_FORCE_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  static void Contains(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point, Bool_v &inside)
  {
    Bool_v unused(false), outside(false);
    GenericKernelForContainsAndInside<Real_v, Bool_v, false>(unplaced, point, unused, outside);
    inside = !outside;
  }

  // BIG QUESTION: DO WE WANT TO GIVE ALL 3 TEMPLATE PARAMETERS
  // -- OR -- DO WE WANT TO DEDUCE Bool_v, Index_t from Real_v???
  template <typename Real_v, typename Inside_t>
  VECGEOM_FORCE_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  static void Inside(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point, Inside_t &inside)
  {
    using Bool_v       = vecCore::Mask_v<Real_v>;
    using InsideBool_v = vecCore::Mask_v<Inside_t>;

    Bool_v completelyinside(false), completelyoutside(false);
    GenericKernelForContainsAndInside<Real_v, Bool_v, true>(unplaced, point, completelyinside, completelyoutside);
    inside = Inside_t(EInside::kSurface);
    vecCore::MaskedAssign(inside, (InsideBool_v)completelyoutside, Inside_t(EInside::kOutside));
    vecCore::MaskedAssign(inside, (InsideBool_v)completelyinside, Inside_t(EInside::kInside));
  }

  template <typename Real_v, typename Bool_v, bool ForInside>
  VECGEOM_FORCE_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  static void GenericKernelForContainsAndInside(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point,
                                                Bool_v &completelyInside, Bool_v &completelyOutside)
  {
    constexpr Precision trapSurfaceTolerance = kTolerance;
    // z-region
    completelyOutside = Abs(point[2]) > MakePlusTolerant<ForInside>(unplaced.fDz, trapSurfaceTolerance);
    if (vecCore::EarlyReturnAllowed() && vecCore::MaskFull(completelyOutside)) {
      return;
    }
    if (ForInside) {
      completelyInside = Abs(point[2]) < MakeMinusTolerant<ForInside>(unplaced.fDz, trapSurfaceTolerance);
    }

#ifndef VECGEOM_PLANESHELL_DISABLE
    unplaced.GetPlanes()->GenericKernelForContainsAndInside<Real_v, ForInside>(point, completelyInside,
                                                                               completelyOutside);
#else
    // here for PLANESHELL=OFF (disabled)
    TrapSidePlane const *fPlanes = unplaced.GetPlanes();
    Real_v dist[4];
    for (unsigned int i = 0; i < 4; ++i) {
      dist[i] = fPlanes[i].fA * point.x() + fPlanes[i].fB * point.y() + fPlanes[i].fC * point.z() + fPlanes[i].fD;
    }

    for (unsigned int i = 0; i < 4; ++i) {
      // is it outside of this side plane?
      completelyOutside = completelyOutside || dist[i] > MakePlusTolerant<ForInside>(0., trapSurfaceTolerance);
      if (ForInside) {
        completelyInside = completelyInside && dist[i] < MakeMinusTolerant<ForInside>(0., trapSurfaceTolerance);
      }
      if (vecCore::EarlyReturnAllowed() && vecCore::MaskFull(completelyOutside)) return;
    }
#endif

    return;
  }

  ////////////////////////////////////////////////////////////////////////////
  //
  // Calculate distance to shape from outside - return kInfLength if no
  // intersection.
  //
  // ALGORITHM: For each component (z-planes, side planes), calculate
  // pair of minimum (smin) and maximum (smax) intersection values for
  // which the particle is in the extent of the shape.  The point of
  // entrance (exit) is found by the largest smin (smallest smax).
  //
  //  If largest smin > smallest smax, the trajectory does not reach
  //  inside the shape.
  //
  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  static void DistanceToIn(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point, Vector3D<Real_v> const &dir,
                           Real_v const &stepMax, Real_v &distance)
  {
    (void)stepMax;
    using Bool_v = vecCore::Mask_v<Real_v>;
    distance     = kInfLength;

    //
    // Step 1: find range of distances along dir between Z-planes (smin, smax)
    //

    // step 1.a) input particle is moving away --> return infinity
    Real_v max = Sign(dir.z()) * unplaced.fDz - point.z(); // z-dist to farthest z-plane

    // done = done || (dir.z()>0.0 && max < MakePlusTolerant<true>(0.));  // check if moving away towards +z
    // done = done || (dir.z()<0.0 && max > MakeMinusTolerant<true>(0.)); // check if moving away towards -z
    Bool_v done(Sign(dir.z()) * max < MakePlusTolerant<true>(0.0)); // if outside + moving away towards +/-z

    // if all particles moving away, we're done
    if (vecCore::EarlyReturnAllowed() && vecCore::MaskFull(done)) return;

    Real_v invdir = Real_v(1.0) / NonZero(dir.z()); // convert distances from z to dir
    Real_v smax   = max * invdir;
    Real_v smin   = (-Sign(dir.z()) * unplaced.fDz - point.z()) * invdir;

#ifndef VECGEOM_PLANESHELL_DISABLE
    // If disttoplanes is such that smin < dist < smax, then distance=disttoplanes
    Real_v disttoplanes = unplaced.GetPlanes()->DistanceToIn(point, dir, smin, smax);
    vecCore::MaskedAssign(distance, !done, disttoplanes);

#else

    // here for VECGEOM_PLANESHELL_DISABLE

    // loop over side planes - find pdist,Comp for each side plane
    Real_v pdist[4], comp[4], vdist[4];
    // EvaluateTrack<Real_v>(unplaced, point, dir, pdist, comp, vdist);

    // auto-vectorizable part of loop
    TrapSidePlane const *fPlanes = unplaced.GetPlanes();
    for (unsigned int i = 0; i < 4; ++i) {
      // Note: normal vector is pointing outside the volume (convention), therefore
      // pdist>0 if point is outside  and  pdist<0 means inside
      pdist[i] = fPlanes[i].fA * point.x() + fPlanes[i].fB * point.y() + fPlanes[i].fC * point.z() + fPlanes[i].fD;

      // Comp is projection of dir over the normal vector of side plane, hence
      // Comp > 0 if pointing ~same direction as normal and Comp<0 if ~opposite to normal
      comp[i] = fPlanes[i].fA * dir.x() + fPlanes[i].fB * dir.y() + fPlanes[i].fC * dir.z();

      vdist[i] = -pdist[i] / NonZero(comp[i]);
    }

    // check special cases
    for (int i = 0; i < 4; ++i) {
      // points fully outside a plane and moving away or parallel to that plane
      done = done || pdist[i] > MakePlusTolerant<true>(0.0) && comp[i] >= 0.0;
      // points at a plane surface and exiting
      done = done || pdist[i] > MakeMinusTolerant<true>(0.0) && comp[i] > 0.0;
    }
    // if all particles moving away, we're done
    if (vecCore::EarlyReturnAllowed() && vecCore::MaskFull(done)) return;

    // this part does not auto-vectorize
    for (unsigned int i = 0; i < 4; ++i) {
      // if outside and moving away, return infinity
      Bool_v posPoint = pdist[i] > MakeMinusTolerant<true>(0.);
      Bool_v posDir   = comp[i] > 0;

      // check if trajectory will intercept plane within current range (smin,smax), otherwise track misses shape
      Bool_v interceptFromInside  = (!posPoint && posDir);
      Bool_v interceptFromOutside = (posPoint && !posDir);

      //.. If dist is such that smin < dist < smax, then adjust either smin or smax
      vecCore::MaskedAssign(smax, interceptFromInside && vdist[i] < smax, vdist[i]);
      vecCore::MaskedAssign(smin, interceptFromOutside && vdist[i] > smin, vdist[i]);
    }

    vecCore::MaskedAssign(distance, !done && smin <= smax, smin);
// vecCore::MaskedAssign(distance, !done && distance < MakeMinusTolerant<true>(0.0), Real_v(-1.0));
#endif
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  static void DistanceToOut(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point,
                            Vector3D<Real_v> const &dir, Real_v const &stepMax, Real_v &distance)
  {
    (void)stepMax;
    using Bool_v = vecCore::Mask_v<Real_v>;

    // step 0: if point is outside any plane --> return -1
    Bool_v outside = Abs(point.z()) > MakePlusTolerant<true>(unplaced.fDz);
    distance       = vecCore::Blend(outside, Real_v(-1.0), InfinityLength<Real_v>());
    Bool_v done(outside);
    if (vecCore::EarlyReturnAllowed() && vecCore::MaskFull(done)) return;

    //
    // Step 1: find range of distances along dir between Z-planes (smin, smax)
    //

    Real_v max   = Sign(dir.z()) * unplaced.fDz - point.z();
    Real_v distz = max / NonZero(dir.z());
    vecCore::MaskedAssign(distance, !done && dir.z() != Real_v(0.), distz);

//
// Step 2: find distances for intersections with side planes.
//

#ifndef VECGEOM_PLANESHELL_DISABLE
    Real_v disttoplanes = unplaced.GetPlanes()->DistanceToOut(point, dir);
    vecCore::MaskedAssign(distance, disttoplanes < distance, disttoplanes);

#else
    //=== Here for VECGEOM_PLANESHELL_DISABLE

    // loop over side planes - find pdist,Comp for each side plane
    Real_v pdist[4], comp[4], vdist[4];
    // EvaluateTrack<Real_v>(unplaced, point, dir, pdist, comp, vdist);

    TrapSidePlane const *fPlanes = unplaced.GetPlanes();
    for (unsigned int i = 0; i < 4; ++i) {
      // Note: normal vector is pointing outside the volume (convention), therefore
      // pdist>0 if point is outside  and  pdist<0 means inside
      pdist[i] = fPlanes[i].fA * point.x() + fPlanes[i].fB * point.y() + fPlanes[i].fC * point.z() + fPlanes[i].fD;

      // Comp is projection of dir over the normal vector of side plane, hence
      // Comp > 0 if pointing ~same direction as normal and Comp<0 if pointing ~opposite to normal
      comp[i] = fPlanes[i].fA * dir.x() + fPlanes[i].fB * dir.y() + fPlanes[i].fC * dir.z();

      vdist[i] = -pdist[i] / NonZero(comp[i]);
    }

    // early return if point is outside of plane
    for (unsigned int i = 0; i < 4; ++i) {
      done = done || (pdist[i] > MakePlusTolerant<true>(0.));
    }
    vecCore::MaskedAssign(distance, done, Real_v(-1.0));
    if (vecCore::EarlyReturnAllowed() && vecCore::MaskFull(done)) return;

    for (unsigned int i = 0; i < 4; ++i) {
      // if track is pointing towards plane and vdist<distance, then distance=vdist
      vecCore::MaskedAssign(distance, !done && comp[i] > 0.0 && vdist[i] < distance, vdist[i]);
    }
#endif
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  static void SafetyToIn(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point, Real_v &safety)
  {
    safety = Abs(point.z()) - unplaced.fDz;

#ifndef VECGEOM_PLANESHELL_DISABLE
    // Get safety over side planes
    unplaced.GetPlanes()->SafetyToIn(point, safety);
#else
    // Loop over side planes
    TrapSidePlane const *fPlanes = unplaced.GetPlanes();
    Real_v dist[4];
    for (int i = 0; i < 4; ++i) {
      dist[i] = fPlanes[i].fA * point.x() + fPlanes[i].fB * point.y() + fPlanes[i].fC * point.z() + fPlanes[i].fD;
    }

    // for (int i = 0; i < 4; ++i) {
    //   vecCore::MaskedAssign(safety, dist[i] > safety, Dist[i]);
    // }
    Real_v safmax = Max(Max(dist[0], dist[1]), Max(dist[2], dist[3]));
    vecCore::MaskedAssign(safety, safmax > safety, safmax);
#endif
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  static void SafetyToOut(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point, Real_v &safety)
  {
    // If point is outside (wrong-side) --> safety to negative value
    safety = unplaced.fDz - Abs(point.z());

    // If all test points are outside, we're done
    if (vecCore::EarlyReturnAllowed()) {
      if (vecCore::MaskFull(safety < kHalfTolerance)) return;
    }

#ifndef VECGEOM_PLANESHELL_DISABLE
    // Get safety over side planes
    unplaced.GetPlanes()->SafetyToOut(point, safety);
#else
    // Loop over side planes
    TrapSidePlane const *fPlanes = unplaced.GetPlanes();

    // auto-vectorizable loop
    Real_v dist[4];
    for (int i = 0; i < 4; ++i) {
      dist[i] = -(fPlanes[i].fA * point.x() + fPlanes[i].fB * point.y() + fPlanes[i].fC * point.z() + fPlanes[i].fD);
    }

    // unvectorizable loop
    // for (int i = 0; i < 4; ++i) {
    //   vecCore::MaskedAssign(safety, Dist[i] < safety, Dist[i]);
    // }

    Real_v safmin = Min(Min(dist[0], dist[1]), Min(dist[2], dist[3]));
    vecCore::MaskedAssign(safety, safmin < safety, safmin);
#endif
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  static Vector3D<Real_v> NormalKernel(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point,
                                       typename vecCore::Mask_v<Real_v> &valid)
  {
    using Bool_v  = vecCore::Mask_v<Real_v>;
    using Index_v = vecCore::Index<Real_v>;

    constexpr double delta = 1000. * kTolerance;
    Vector3D<Real_v> normal(0.);
    valid         = true;
    Real_v safety = InfinityLength<Real_v>();
    vecCore::Index<Real_v> iface;

#ifndef VECGEOM_PLANESHELL_DISABLE
    // Get normal from side planes -- PlaneShell case
    // normal = unplaced.GetPlanes()->NormalKernel(point, valid);
    iface = unplaced.GetPlanes()->ClosestFace(point, safety);

#else
    // Loop over side planes
    TrapSidePlane const *fPlanes = unplaced.GetPlanes();

    // vectorizable loop
    Real_v dist[4];
    for (int i = 0; i < 4; ++i) {
      dist[i] = -(fPlanes[i].fA * point.x() + fPlanes[i].fB * point.y() + fPlanes[i].fC * point.z() + fPlanes[i].fD);
    }

    // non-vectorizable part
    for (int i = 0; i < 4; ++i) {
      // valid becomes false if point is outside any plane by more than delta (dist < -delta here means outside)
      valid = valid && dist[i] >= -delta;

      const Bool_v closer = Abs(dist[i]) < safety;
      vecCore::MaskedAssign(safety, closer, Abs(dist[i]));
      vecCore::MaskedAssign(iface, closer, Index_v(i));
      // vecCore::MaskedAssign(normal, Abs(dist[i]) <= delta, normal + unplaced.normals[i]); // averaging
    }
#endif
    // check if normal is valid w.r.t. z-planes
    valid = valid && (unplaced.fDz - Abs(point[2])) >= -delta;

    if (vecCore::EarlyReturnAllowed() && vecCore::MaskEmpty(valid))
      return (normal + Vector3D<Real_v>(1.0)).Normalized();

    // then check z-planes
    // testing -fDZ (still keeping safety for next face)
    Real_v distTmp = Abs(point.z() + unplaced.fDz);
    Bool_v closer  = distTmp < safety;
    vecCore::MaskedAssign(safety, closer, distTmp);
    vecCore::MaskedAssign(iface, closer, Index_v(4));
    // vecCore::MaskedAssign(normal[2], distTmp < delta, normal[2] - Real_v(1.0));  // averaging

    // testing +fDZ (no need to update safety)
    distTmp = Abs(point.z() - unplaced.fDz);
    vecCore::MaskedAssign(iface, distTmp < safety, Index_v(5));
    // vecCore::MaskedAssign(normal[2], distTmp < delta, normal[2] + Real_v(1.0));  // averaging

    // vecCore::MaskedAssign(normal, normal.Mag2() < delta, unplaced.normals[iface]);
    normal = unplaced.normals[iface];

    // returned vector must be normalized
    // vecCore::MaskedAssign(normal, normal.Mag2() > 0., normal.Normalized());

    return normal;
  }
};

} // end of inline namespace
} // end of global namespace

#endif // VECGEOM_VOLUMES_KERNEL_TRAPEZOIDIMPLEMENTATION_H_
