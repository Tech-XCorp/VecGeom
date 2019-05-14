/*
 * ConeImplementation.h
 *
 *  Created on: May 14, 2014
 *      Author: swenzel
 */

/// History notes:
/// revision + moving to Vectorized Cone Kernels (Raman Sehgal)
/// May-June 2017: revision + moving to new Structure (Raman Sehgal)
/// 20180323 Guilherme Lima  Adapted to new UnplacedVolume factory

#ifndef VECGEOM_VOLUMES_KERNEL_CONEIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_CONEIMPLEMENTATION_H_

#include "base/Vector3D.h"
#include "volumes/kernel/GenericKernels.h"
#include "volumes/kernel/shapetypes/ConeTypes.h"
#include "volumes/ConeStruct.h"
#include <cstdio>
#include "volumes/ConeUtilities.h"
#define kConeTolerance 1e-7
#define kHalfConeTolerance 0.5 * kConeTolerance

namespace vecgeom {

VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE(struct, ConeImplementation, typename);

inline namespace VECGEOM_IMPL_NAMESPACE {

template <typename T>
class SPlacedCone;
template <typename T>
class SUnplacedCone;

template <typename coneTypeT>
struct ConeImplementation {

  using UnplacedStruct_t = ConeStruct<double>;
  using UnplacedVolume_t = SUnplacedCone<coneTypeT>;
  using PlacedShape_t    = SPlacedCone<UnplacedVolume_t>;

  VECCORE_ATT_HOST_DEVICE
  static void PrintType() {}

  template <typename Stream>
  static void PrintType(Stream &s, int transCodeT = translation::kGeneric, int rotCodeT = rotation::kGeneric)
  {
    s << "SpecializedCone<" << transCodeT << "," << rotCodeT << ">";
  }

  template <typename Stream>
  static void PrintImplementationType(Stream & /*s*/)
  {
  }

  template <typename Stream>
  static void PrintUnplacedType(Stream & /*s*/)
  {
  }

  /* A Function that will just check if the point is on the CONICAL (circle) edge
   * assuming that it is on either lowerZ or upperZ
   *
   * Beware : It will not do any checks on Z
   */
  template <typename Real_v, bool ForInnerSurface, bool ForLowerZ>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static typename vecCore::Mask_v<Real_v> IsOnRing(UnplacedStruct_t const &cone, Vector3D<Real_v> const &point)
  {
    using Bool_v = typename vecCore::Mask_v<Real_v>;

    Real_v rad2 = point.Perp2();
    Bool_v onRing(false);

    if (ForLowerZ) {
      if (ForInnerSurface) {
        onRing = (rad2 <= MakePlusTolerantSquare<true>(cone.fRmin1, cone.fRmin1 * cone.fRmin1)) &&
                 (rad2 >= MakeMinusTolerantSquare<true>(cone.fRmin1, cone.fRmin1 * cone.fRmin1));
      } else {
        onRing = (rad2 <= MakePlusTolerantSquare<true>(cone.fRmax1, cone.fRmax1 * cone.fRmax1)) &&
                 (rad2 >= MakeMinusTolerantSquare<true>(cone.fRmax1, cone.fRmax1 * cone.fRmax1));
      }
    } else {
      if (ForInnerSurface) {
        onRing = (rad2 <= MakePlusTolerantSquare<true>(cone.fRmin2, cone.fRmin2 * cone.fRmin2)) &&
                 (rad2 >= MakeMinusTolerantSquare<true>(cone.fRmin2, cone.fRmin2 * cone.fRmin2));
      } else {
        onRing = (rad2 <= MakePlusTolerantSquare<true>(cone.fRmax2, cone.fRmax2 * cone.fRmax2)) &&
                 (rad2 >= MakeMinusTolerantSquare<true>(cone.fRmax2, cone.fRmax2 * cone.fRmax2));
      }
    }

    return onRing;
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void Contains(UnplacedStruct_t const &cone, Vector3D<Real_v> const &point,
                       typename vecCore::Mask_v<Real_v> &inside)
  {
    typedef typename vecCore::Mask_v<Real_v> Bool_v;
    Bool_v unused;
    Bool_v outside;
    ConeHelpers<Real_v, coneTypeT>::template GenericKernelForContainsAndInside<false>(cone, point, unused, outside);
    inside = !outside;
  }

  template <typename Real_v, typename Inside_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void Inside(UnplacedStruct_t const &cone, Vector3D<Real_v> const &point, Inside_v &inside)
  {
    ConeHelpers<Real_v, coneTypeT>::template Inside<Inside_v>(cone, point, inside);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void DistanceToIn(UnplacedStruct_t const &cone, Vector3D<Real_v> const &point, Vector3D<Real_v> const &dir,
                           Real_v const & /*stepMax*/, Real_v &distance)
  {
    using namespace ConeUtilities;
    using namespace ConeTypes;
    typedef Real_v Float_t;
    typedef typename vecCore::Mask_v<Real_v> Bool_t;

    Bool_t done(false);

    //=== First, for points outside and moving away --> return infinity
    distance = kInfLength;

    // outside of Z range and going away?
    Float_t distz          = Abs(point.z()) - cone.fDz; // avoid a division for now
    Bool_t outZAndGoingOut = (distz > kHalfConeTolerance && (point.z() * dir.z()) >= 0) ||
                             (Abs(distz) < kHalfConeTolerance && (point.z() * dir.z()) > 0);
    done |= outZAndGoingOut;
    if (vecCore::MaskFull(done)) return;

    // outside of outer cone and going away?
    Float_t outerRadIrTol  = GetRadiusOfConeAtPoint<Real_v, false>(cone, point.z()) - cone.fOuterTolerance;
    Float_t outerRadIrTol2 = outerRadIrTol * outerRadIrTol;
    Float_t rsq            = point.Perp2(); // point.x()*point.x() + point.y()*point.y();

    done |= rsq > outerRadIrTol2 && (dir.Dot(GetNormal<Real_v, false>(cone, point)) >= 0.);
    if (vecCore::MaskFull(done)) return;

    //=== Next, check all dimensions of the cone, whether points are inside -->
    // return -1

    vecCore__MaskedAssignFunc(distance, !done, Float_t(-1.0));

    // For points inside z-range, return -1
    Bool_t inside = distz < -kHalfConeTolerance;

    inside &= rsq < outerRadIrTol2;

    if (checkRminTreatment<coneTypeT>(cone)) {
      Float_t innerRadIrTol  = GetRadiusOfConeAtPoint<Real_v, true>(cone, point.z()) + cone.fInnerTolerance;
      Float_t innerRadIrTol2 = innerRadIrTol * innerRadIrTol;
      inside &= rsq > innerRadIrTol2;
    }
    if (checkPhiTreatment<coneTypeT>(cone)) { // && !vecCore::MaskEmpty(inside)) {
      Bool_t insector;
      PointInCyclicalSector<Real_v, coneTypeT, false, false>(cone, point.x(), point.y(), insector);
      inside &= insector;
    }
    done |= inside;
    if (vecCore::MaskFull(done)) return;

    //=== Next step: check if z-plane is the right entry point (both r,phi
    // should be valid at z-plane crossing)
    vecCore__MaskedAssignFunc(distance, !done, Float_t(kInfLength));

    distz /= NonZero(Abs(dir.z()));

    Bool_t isOnZPlaneAndMovingInside(false);

    Bool_t isGoingUp          = dir.z() > 0.;
    isOnZPlaneAndMovingInside = ((isGoingUp && point.z() < 0. && Abs(distz) < kHalfTolerance) ||
                                 (!isGoingUp && point.z() > 0. && Abs(distz) < kHalfTolerance));
    vecCore__MaskedAssignFunc(distz, !done && isOnZPlaneAndMovingInside, Float_t(0.));

#ifdef EDGE_POINTS
    Bool_t newCond = (IsOnRing<Backend, false, true>(cone, point)) || (IsOnRing<Backend, true, true>(cone, point)) ||
                     (IsOnRing<Backend, false, false>(cone, point)) || (IsOnRing<Backend, true, false>(cone, point));

    vecCore::MaskedAssign(distz, newCond, kInfLength);
#endif

    Float_t hitx = point.x() + distz * dir.x();
    Float_t hity = point.y() + distz * dir.y();

    Float_t r2 = (hitx * hitx) + (hity * hity);

    Precision innerZTol         = cone.fTolIz;
    Bool_t isHittingTopPlane    = (point.z() >= innerZTol) && (r2 <= cone.fSqRmax2 + kTolerance);  // cone.fSqRmax2Tol
    Bool_t isHittingBottomPlane = (point.z() <= -innerZTol) && (r2 <= cone.fSqRmax1 + kTolerance); // cone.fSqRmax1Tol
    Bool_t okz                  = (isHittingTopPlane || isHittingBottomPlane);

    if (checkRminTreatment<coneTypeT>(cone)) {
      isHittingTopPlane &= (r2 >= cone.fSqRmin2 - kTolerance);    // cone.fSqRmin2Tol
      isHittingBottomPlane &= (r2 >= cone.fSqRmin1 - kTolerance); // cone.fSqRmin1Tol
      okz &= ((isHittingTopPlane || isHittingBottomPlane));
    }

    if (checkPhiTreatment<coneTypeT>(cone)) {
      Bool_t insector;
      PointInCyclicalSector<Real_v, coneTypeT, false>(cone, hitx, hity, insector);
      okz &= insector;
    }
    vecCore::MaskedAssign(distance, !done && okz, distz);
    done |= okz;
    if (vecCore::MaskFull(done)) return;

    Float_t dist_rOuter(kInfLength);
    Bool_t ok_outerCone =
        ConeHelpers<Real_v, coneTypeT>::template DetectIntersectionAndCalculateDistanceToConicalSurface<true, false>(
            cone, point, dir, dist_rOuter);
    ok_outerCone &= dist_rOuter < distance;
    vecCore::MaskedAssign(distance, !done && ok_outerCone, dist_rOuter);
    done |= ok_outerCone;
    if (vecCore::MaskFull(done)) return;

    Float_t dist_rInner(kInfLength);
    if (checkRminTreatment<coneTypeT>(cone)) {

      Bool_t ok_innerCone =
          ConeHelpers<Real_v, coneTypeT>::template DetectIntersectionAndCalculateDistanceToConicalSurface<true, true>(
              cone, point, dir, dist_rInner);
      ok_innerCone &= dist_rInner < distance;
      vecCore::MaskedAssign(distance, !done && ok_innerCone, dist_rInner);
    }

    if (checkPhiTreatment<coneTypeT>(cone)) {

      evolution::Wedge const &w = cone.fPhiWedge;

      Float_t dist_phi;
      Bool_t ok_phi;
      PhiPlaneTrajectoryIntersection<Real_v, coneTypeT, SectorType<coneTypeT>::value != kOnePi, true>(
          cone.fAlongPhi1x, cone.fAlongPhi1y, w.GetNormal1().x(), w.GetNormal1().y(), cone, point, dir, dist_phi,
          ok_phi);
      ok_phi &= dist_phi < distance;
      vecCore::MaskedAssign(distance, !done && ok_phi, dist_phi);
      done |= ok_phi;

      if (SectorType<coneTypeT>::value != kOnePi) {

        PhiPlaneTrajectoryIntersection<Real_v, coneTypeT, true, true>(cone.fAlongPhi2x, cone.fAlongPhi2y,
                                                                      w.GetNormal2().x(), w.GetNormal2().y(), cone,
                                                                      point, dir, dist_phi, ok_phi);

        vecCore::MaskedAssign(distance, ok_phi && dist_phi < distance, dist_phi);
      }
    }
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void DistanceToOut(UnplacedStruct_t const &cone, Vector3D<Real_v> const &point,
                            Vector3D<Real_v> const &direction, Real_v const & /*stepMax*/, Real_v &distance)
  {

    distance = kInfLength;
    using namespace ConeUtilities;
    using namespace ConeTypes;

    typedef typename vecCore::Mask_v<Real_v> Bool_t;

    Bool_t done(false);

    // Using this logic will improve performance of Scalar code
    Real_v distz         = Abs(point.z()) - cone.fDz;
    Real_v rsq           = point.Perp2();
    Real_v outerRadOrTol = ConeUtilities::GetRadiusOfConeAtPoint<Real_v, false>(cone, point.z()) + cone.fOuterTolerance;
    Real_v outerRadOrTol2 = outerRadOrTol * outerRadOrTol;

    //=== Next, check all dimensions of the cone, whether points are inside -->
    // return -1
    vecCore__MaskedAssignFunc(distance, !done, Real_v(-1.0));

    // For points inside z-range, return -1
    Bool_t outside = distz > kHalfConeTolerance || rsq > outerRadOrTol2;
    done |= outside;
    if (vecCore::MaskFull(done)) return;

    if (checkRminTreatment<coneTypeT>(cone) && !vecCore::MaskFull(outside)) {
      Real_v innerRadOrTol =
          ConeUtilities::GetRadiusOfConeAtPoint<Real_v, true>(cone, point.z()) - cone.fInnerTolerance;
      Real_v innerRadOrTol2 = innerRadOrTol * innerRadOrTol;
      outside |= rsq < innerRadOrTol2;
      done |= outside;
      if (vecCore::MaskFull(done)) return;
    }
    if (checkPhiTreatment<coneTypeT>(cone) && !vecCore::MaskEmpty(outside)) {
      Bool_t insector;
      ConeUtilities::PointInCyclicalSector<Real_v, coneTypeT, false, false>(cone, point.x(), point.y(), insector);
      outside |= !insector;
    }
    done |= outside;
    if (vecCore::MaskFull(done)) return;
    Bool_t isGoingUp   = direction.z() > 0.;
    Bool_t isGoingDown = direction.z() < 0.;
    Bool_t isOnZPlaneAndMovingOutside(false);
    isOnZPlaneAndMovingOutside = !outside && ((isGoingUp && point.z() > 0. && Abs(distz) < kHalfTolerance) ||
                                              (!isGoingUp && point.z() < 0. && Abs(distz) < kHalfTolerance));
    vecCore__MaskedAssignFunc(distance, !done && isOnZPlaneAndMovingOutside, Real_v(0.));
    done |= isOnZPlaneAndMovingOutside;
    if (vecCore::MaskFull(done)) return;

    //=== Next step: check if z-plane is the right entry point (both r,phi
    // should be valid at z-plane crossing)
    vecCore__MaskedAssignFunc(distance, !done, Real_v(kInfLength));

    Precision fDz  = cone.fDz;
    Real_v dirZInv = 1. / NonZero(direction.z());
    vecCore__MaskedAssignFunc(distance, isGoingUp, (fDz - point.z()) * dirZInv);
    vecCore__MaskedAssignFunc(distance, isGoingDown, (-fDz - point.z()) * dirZInv);

    Real_v dist_rOuter(kInfLength);
    Bool_t ok_outerCone =
        ConeHelpers<Real_v, coneTypeT>::template DetectIntersectionAndCalculateDistanceToConicalSurface<false, false>(
            cone, point, direction, dist_rOuter);

    vecCore::MaskedAssign(distance, !done && ok_outerCone && dist_rOuter < distance, dist_rOuter);

    Real_v dist_rInner(kInfLength);
    if (checkRminTreatment<coneTypeT>(cone)) {
      Bool_t ok_innerCone =
          ConeHelpers<Real_v, coneTypeT>::template DetectIntersectionAndCalculateDistanceToConicalSurface<false, true>(
              cone, point, direction, dist_rInner);
      vecCore::MaskedAssign(distance, !done && ok_innerCone && dist_rInner < distance, dist_rInner);
    }

    if (checkPhiTreatment<coneTypeT>(cone)) {

      Bool_t isOnStartPhi      = ConeUtilities::IsOnStartPhi<Real_v>(cone, point);
      Bool_t isOnEndPhi        = ConeUtilities::IsOnEndPhi<Real_v>(cone, point);
      Vector3D<Real_v> normal1 = cone.fPhiWedge.GetNormal1();
      Vector3D<Real_v> normal2 = cone.fPhiWedge.GetNormal2();
      Bool_t cond = (isOnStartPhi && direction.Dot(-normal1) > 0.) || (isOnEndPhi && direction.Dot(-normal2) > 0.);
      vecCore__MaskedAssignFunc(distance, !done && cond, Real_v(0.));
      done |= cond;
      if (vecCore::MaskFull(done)) return;

      Real_v dist_phi;
      Bool_t ok_phi;
      evolution::Wedge const &w = cone.fPhiWedge;
      PhiPlaneTrajectoryIntersection<Real_v, coneTypeT, SectorType<coneTypeT>::value != kOnePi, false>(
          cone.fAlongPhi1x, cone.fAlongPhi1y, w.GetNormal1().x(), w.GetNormal1().y(), cone, point, direction, dist_phi,
          ok_phi);
      ok_phi &= dist_phi < distance;
      vecCore::MaskedAssign(distance, !done && ok_phi, dist_phi);
      done |= ok_phi;

      if (SectorType<coneTypeT>::value != kOnePi) {
        ConeUtilities::PhiPlaneTrajectoryIntersection<Real_v, coneTypeT, true, false>(
            cone.fAlongPhi2x, cone.fAlongPhi2y, w.GetNormal2().x(), w.GetNormal2().y(), cone, point, direction,
            dist_phi, ok_phi);
        vecCore::MaskedAssign(distance, ok_phi && dist_phi < distance, dist_phi);
      }
    }
    vecCore__MaskedAssignFunc(distance, distance < 0. && Abs(distance) < kTolerance, Real_v(0.));
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void SafetyToIn(UnplacedStruct_t const &cone, Vector3D<Real_v> const &point, Real_v &safety)
  {
    using namespace ConeUtilities;
    using namespace ConeTypes;
    safety = -kInfLength;
    typedef typename vecCore::Mask_v<Real_v> Bool_t;
    typedef Real_v Float_t;

    Bool_t done(false);
    Precision fDz = cone.fDz;
    Float_t distz = Abs(point.z()) - fDz;

    // Next, check all dimensions of the cone, whether points are inside -->
    // return -1
    vecCore__MaskedAssignFunc(safety, !done, Float_t(-1.0));

    // For points inside z-range, return -1
    Bool_t inside = distz < -kHalfConeTolerance;

    // This logic to check if the point is inside is far better then
    // using GenericKernel and will improve performance.
    Float_t outerRadIrTol  = GetRadiusOfConeAtPoint<Real_v, false>(cone, point.z()) - kTolerance;
    Float_t outerRadIrTol2 = outerRadIrTol * outerRadIrTol;
    Float_t rsq            = point.Perp2();
    inside &= rsq < outerRadIrTol2;

    if (checkRminTreatment<coneTypeT>(cone)) {
      Float_t innerRadIrTol  = GetRadiusOfConeAtPoint<Real_v, true>(cone, point.z()) + kTolerance;
      Float_t innerRadIrTol2 = innerRadIrTol * innerRadIrTol;
      inside &= rsq > innerRadIrTol2;
    }
    if (checkPhiTreatment<coneTypeT>(cone) && !vecCore::MaskEmpty(inside)) {
      Bool_t insector;
      PointInCyclicalSector<Real_v, coneTypeT, false, false>(cone, point.x(), point.y(), insector);
      inside &= insector;
    }
    done |= inside;
    if (vecCore::MaskFull(done)) return;

    // Once it is checked that the point is inside or not, safety can be set to 0.
    // This will serve the case that the point is on the surface. So no need to check
    // that the point is really on surface.
    vecCore__MaskedAssignFunc(safety, !done, Float_t(0.));

    // Now if the point is neither inside nor on surface, then it should be outside
    // and the safety should be set to some finite value, which is done by below logic

    Float_t safeZ                = Abs(point.z()) - fDz;
    Float_t safeDistOuterSurface = -SafeDistanceToConicalSurface<Real_v, false>(cone, point);

    Float_t safeDistInnerSurface(-kInfLength);
    if (checkRminTreatment<coneTypeT>(cone)) {
      safeDistInnerSurface = -SafeDistanceToConicalSurface<Real_v, true>(cone, point);
    }

    vecCore__MaskedAssignFunc(safety, !done, Max(safeZ, Max(safeDistOuterSurface, safeDistInnerSurface)));

    if (checkPhiTreatment<coneTypeT>(cone)) {
      Float_t safetyPhi = cone.fPhiWedge.SafetyToIn<Real_v>(point);
      vecCore__MaskedAssignFunc(safety, !done, Max(safetyPhi, safety));
    }

    vecCore__MaskedAssignFunc(safety, vecCore::math::Abs(safety) < kTolerance, Float_t(0.));
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void SafetyToOut(UnplacedStruct_t const &cone, Vector3D<Real_v> const &point, Real_v &safety)
  {

    using namespace ConeUtilities;
    using namespace ConeTypes;
    safety = kInfLength;
    typedef typename vecCore::Mask_v<Real_v> Bool_t;
    typedef Real_v Float_t;

    Bool_t done(false);

    Float_t distz = Abs(point.z()) - cone.fDz;
    Float_t rsq   = point.Perp2();

    // This logic to check if the point is outside is far better then
    // using GenericKernel and will improve performance.
    Float_t outerRadOrTol  = GetRadiusOfConeAtPoint<Real_v, false>(cone, point.z()) + kTolerance;
    Float_t outerRadOrTol2 = outerRadOrTol * outerRadOrTol;

    //=== Next, check all dimensions of the cone, whether points are inside -->
    // return -1
    vecCore__MaskedAssignFunc(safety, !done, Float_t(-1.0));

    // For points outside z-range, return -1
    Bool_t outside = distz > kHalfConeTolerance;

    outside |= rsq > outerRadOrTol2;

    if (checkRminTreatment<coneTypeT>(cone)) {
      Float_t innerRadOrTol  = GetRadiusOfConeAtPoint<Real_v, true>(cone, point.z()) - kTolerance;
      Float_t innerRadOrTol2 = innerRadOrTol * innerRadOrTol;
      outside |= rsq < innerRadOrTol2;
    }

    if (checkPhiTreatment<coneTypeT>(cone) && !vecCore::MaskEmpty(outside)) {

      Bool_t insector;
      ConeUtilities::PointInCyclicalSector<Real_v, coneTypeT, false, false>(cone, point.x(), point.y(), insector);
      outside |= !insector;
    }
    done |= outside;
    if (vecCore::MaskFull(done)) return;

    // Once it is checked that the point is inside or not, safety can be set to 0.
    // This will serve the case that the point is on the surface. So no need to check
    // that the point is really on surface.
    vecCore__MaskedAssignFunc(safety, !done, Float_t(0.));

    // Now if the point is neither outside nor on surface, then it should be inside
    // and the safety should be set to some finite value, which is done by below logic

    Precision fDz = cone.fDz;
    Float_t safeZ = fDz - Abs(point.z());

    Float_t safeDistOuterSurface = SafeDistanceToConicalSurface<Real_v, false>(cone, point);
    Float_t safeDistInnerSurface(kInfLength);
    if (checkRminTreatment<coneTypeT>(cone)) {
      safeDistInnerSurface = SafeDistanceToConicalSurface<Real_v, true>(cone, point);
    } else {
    }

    vecCore__MaskedAssignFunc(safety, !done, Min(safeZ, Min(safeDistOuterSurface, safeDistInnerSurface)));

    if (checkPhiTreatment<coneTypeT>(cone)) {
      Float_t safetyPhi = cone.fPhiWedge.SafetyToOut<Real_v>(point);
      vecCore__MaskedAssignFunc(safety, !done, Min(safetyPhi, safety));
    }
    vecCore__MaskedAssignFunc(safety, vecCore::math::Abs(safety) < kTolerance, Float_t(0.));
  }

  template <typename Real_v, bool ForInnerSurface>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static Real_v SafeDistanceToConicalSurface(UnplacedStruct_t const &cone, Vector3D<Real_v> const &point)
  {

    typedef Real_v Float_t;
    Float_t rho = point.Perp();
    if (ForInnerSurface) {
      Float_t pRMin = cone.fTanRMin * point.z() + (cone.fRmin1 + cone.fRmin2) * 0.5; // cone.fRminAv;
      return (rho - pRMin) * cone.fInvSecRMin;
    } else {
    	Float_t pRMax(0.);
    	if (cone.fOriginalRmax1 == cone.fOriginalRmax2)
    	    pRMax = Real_v(cone.fOriginalRmax1);
    	else
    		pRMax = cone.fTanRMax * point.z() + (cone.fRmax1 + cone.fRmax2) * 0.5; // cone.fRmaxAv;
      return (pRMax - rho) * cone.fInvSecRMax;
    }
  }
};
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif
