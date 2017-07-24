/*
 * CutTubeStruct.h
 *
 *  Created on: 03.11.2016
 *      Author: mgheata
 */
#ifndef VECGEOM_CUTTUBESTRUCT_H_
#define VECGEOM_CUTTUBESTRUCT_H_

#include <VecCore/VecCore>
#include "base/Vector3D.h"

#include "TubeStruct.h"
#include "CutPlanes.h"

namespace vecgeom {

inline namespace VECGEOM_IMPL_NAMESPACE {

// a plain and lightweight struct to encapsulate data members of a cut tube
template <typename T = double>
struct CutTubeStruct {
  T fDz;                     //< Z half length
  TubeStruct<T> fTubeStruct; //< Tube parameters
  CutPlanes fCutPlanes;      //< Cut planes

  T fCosPhi1; //< Cosine of phi
  T fSinPhi1; //< Sine of phi
  T fCosPhi2; //< Cosine of phi+dphi
  T fSinPhi2; //< Sine of phi+dphi

  // constructors

  VECCORE_ATT_HOST_DEVICE
  CutTubeStruct() : fTubeStruct(0., 0., 0., 0., 0.), fCutPlanes() {}

  VECCORE_ATT_HOST_DEVICE
  CutTubeStruct(T const &rmin, T const &rmax, T const &z, T const &sphi, T const &dphi, Vector3D<T> const &bottomNormal,
                Vector3D<T> const &topNormal)
      : fDz(z), fTubeStruct(rmin, rmax, kInfLength, sphi, dphi), fCutPlanes()
  {
    fCutPlanes.Set(0, bottomNormal.Unit(), Vector3D<T>(0., 0., -z));
    fCutPlanes.Set(1, topNormal.Unit(), Vector3D<T>(0., 0., z));
    fCosPhi1 = vecCore::math::Cos(sphi);
    fSinPhi1 = vecCore::math::Sin(sphi);
    fCosPhi2 = vecCore::math::Cos(sphi + dphi);
    fSinPhi2 = vecCore::math::Sin(sphi + dphi);
  }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  TubeStruct<T> const &GetTubeStruct() const { return fTubeStruct; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  CutPlanes const &GetCutPlanes() const { return fCutPlanes; }
};
}
} // end global namespace

#endif
