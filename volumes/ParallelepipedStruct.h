/*
 * ParallelepipedStruct.h
 *
 *  Created on: 11.07.2016
 *      Author: mgheata
 */

#ifndef VECGEOM_VOLUMES_PARALLELEPIPEDSTRUCT_H_
#define VECGEOM_VOLUMES_PARALLELEPIPEDSTRUCT_H_
#include "base/Global.h"

namespace vecgeom {

inline namespace VECGEOM_IMPL_NAMESPACE {

// A plain struct without member functions to encapsulate just the parameters
// of a parallelepiped
template <typename T = double>
struct ParallelepipedStruct {
  Vector3D<T> fDimensions; /** Dimensions dx, dy, dx */
  T fAlpha;                /** Angle dx versus dy (degrees)*/
  T fTheta;                /** Theta angle of parallelepiped axis*/
  T fPhi;                  /** Phi angle of parallelepiped axis*/
  T fCtx;                  /** Cosine of xz angle */
  T fCty;                  /** Cosine of yz angle */
  Vector3D<T> fNormals[3]; /** Precomputed normals */

  // Precomputed values computed from parameters
  T fTanAlpha;       /** Tangent of alpha angle */
  T fTanThetaSinPhi; /** tan(theta)*sin(phi) */
  T fTanThetaCosPhi; /** tan(theta)*cos(phi) */

  VECGEOM_CUDA_HEADER_BOTH
  ParallelepipedStruct(Vector3D<T> const &dim, const T alpha, const T theta, const T phi)
      : fDimensions(dim), fAlpha(0), fTheta(0), fPhi(0), fCtx(0), fCty(0), fTanAlpha(0), fTanThetaSinPhi(0),
        fTanThetaCosPhi(0)
  {
    SetAlpha(alpha);
    SetThetaAndPhi(theta, phi);
  }

  VECGEOM_CUDA_HEADER_BOTH
  ParallelepipedStruct(const T x, const T y, const T z, const T alpha, const T theta, const T phi)
      : fDimensions(x, y, z), fAlpha(0), fTheta(0), fPhi(0), fCtx(0), fCty(0), fTanAlpha(0), fTanThetaSinPhi(0),
        fTanThetaCosPhi(0)
  {
    SetAlpha(alpha);
    SetThetaAndPhi(theta, phi);
  }

  VECGEOM_CUDA_HEADER_BOTH
  void SetAlpha(const T alpha)
  {
    fAlpha    = alpha;
    fTanAlpha = tan(kDegToRad * alpha);
    ComputeNormals();
  }

  VECGEOM_CUDA_HEADER_BOTH
  void SetTheta(const T theta) { SetThetaAndPhi(theta, fPhi); }

  VECGEOM_CUDA_HEADER_BOTH
  void SetPhi(const T phi) { SetThetaAndPhi(fTheta, phi); }

  VECGEOM_CUDA_HEADER_BOTH
  void SetThetaAndPhi(const T theta, const T phi)
  {
    fTheta          = theta;
    fPhi            = phi;
    fTanThetaCosPhi = tan(kDegToRad * fTheta) * cos(kDegToRad * fPhi);
    fTanThetaSinPhi = tan(kDegToRad * fTheta) * sin(kDegToRad * fPhi);
    ComputeNormals();
  }

  VECGEOM_CUDA_HEADER_BOTH
  void ComputeNormals()
  {
    Vector3D<T> v(sin(kDegToRad * fTheta) * cos(kDegToRad * fPhi), sin(kDegToRad * fTheta) * sin(kDegToRad * fPhi),
                  cos(kDegToRad * fTheta));
    Vector3D<T> vx(1., 0., 0.);
    Vector3D<T> vy(-sin(kDegToRad * fAlpha), -cos(kDegToRad * fAlpha), 0.);
    fNormals[0] = v.Cross(vy);
    fNormals[0].Normalize();
    fNormals[1] = v.Cross(vx);
    fNormals[1].Normalize();
    fNormals[2].Set(0., 0., 1.);
    fCtx = 1.0 / sqrt(1. + fTanAlpha * fTanAlpha + fTanThetaCosPhi * fTanThetaCosPhi);
    fCty = 1.0 / sqrt(1. + fTanThetaSinPhi * fTanThetaSinPhi);
  }
};
}
} // end

#endif