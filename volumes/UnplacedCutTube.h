/*
 * UnplacedCutTube.h
 *
 *  Created on: 03.11.2016
 *      Author: mgheata
 */

#ifndef VECGEOM_VOLUMES_UNPLACEDCUTTUBE_H_
#define VECGEOM_VOLUMES_UNPLACEDCUTTUBE_H_

#include "base/Global.h"
#include "base/AlignedBase.h"
#include "volumes/UnplacedVolume.h"
#include "volumes/CutTubeStruct.h"
#include "volumes/kernel/CutTubeImplementation.h"
#include "volumes/UnplacedVolumeImplHelper.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class UnplacedCutTube;);
VECGEOM_DEVICE_DECLARE_CONV(class, UnplacedCutTube);

inline namespace VECGEOM_IMPL_NAMESPACE {

class UnplacedCutTube : public SIMDUnplacedVolumeImplHelper<CutTubeImplementation>, public AlignedBase {

private:
  CutTubeStruct<double> fCutTube; //> The cut tube data holder

public:
  VECCORE_ATT_HOST_DEVICE
  UnplacedCutTube() : fCutTube() {}

  VECCORE_ATT_HOST_DEVICE
  UnplacedCutTube(Precision const &rmin, Precision const &rmax, Precision const &z, Precision const &sphi,
                  Precision const &dphi, Vector3D<Precision> const &bottomNormal, Vector3D<Precision> const &topNormal)
      : fCutTube(rmin, rmax, z, sphi, dphi, bottomNormal, topNormal)
  {
    // Constructor
    // Make sure that all points on top surface have z>0 and all points on
    // bottom surface have z<0. This avoids degenerations where top and bottom
    // intersect each other.
    if (bottomNormal.z() > 0 || bottomNormal.z() * bottomNormal.z() < rmax / (rmax + z) || topNormal.z() < 0 ||
        topNormal.z() * topNormal.z() < rmax / (rmax + z)) {
      Print();
#ifndef VECCORE_CUDA
      throw std::runtime_error("Illegal normal direction for cut planes");
#endif
    }
    DetectConvexity();
  }

  VECCORE_ATT_HOST_DEVICE
  UnplacedCutTube(Precision const &rmin, Precision const &rmax, Precision const &z, Precision const &sphi,
                  Precision const &dphi, Precision const &bx, Precision const &by, Precision const &bz,
                  Precision const &tx, Precision const &ty, Precision const &tz)
      : fCutTube(rmin, rmax, z, sphi, dphi, Vector3D<Precision>(bx, by, bz), Vector3D<Precision>(tx, ty, tz))
  {
    // Constructor
    if (bz > 0 || bz * bz < rmax / (rmax + z) || tz < 0 || tz * tz < rmax / (rmax + z)) {
      Print();
#ifndef VECCORE_CUDA
      throw std::runtime_error("Illegal normal direction for cut planes");
#endif
    }
    DetectConvexity();
  }

  VECCORE_ATT_HOST_DEVICE
  virtual ~UnplacedCutTube() = default;

  VECCORE_ATT_HOST_DEVICE
  UnplacedCutTube(UnplacedCutTube const &other) : fCutTube(other.fCutTube) {}

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  TubeStruct<double> const &GetTubeStruct() const { return fCutTube.fTubeStruct; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  CutTubeStruct<double> const &GetStruct() const { return fCutTube; }

  VECCORE_ATT_HOST_DEVICE
  void DetectConvexity();

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision rmin() const { return fCutTube.fTubeStruct.fRmin; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision rmax() const { return fCutTube.fTubeStruct.fRmax; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision z() const { return fCutTube.fDz; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision sphi() const { return fCutTube.fTubeStruct.fSphi; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision dphi() const { return fCutTube.fTubeStruct.fDphi; }

  VECGEOM_FORCE_INLINE
  void SetRMin(Precision const &_rmin) { fCutTube.fTubeStruct.SetRMin(_rmin); }

  VECGEOM_FORCE_INLINE
  void SetRMax(Precision const &_rmax) { fCutTube.fTubeStruct.SetRMax(_rmax); }

  VECGEOM_FORCE_INLINE
  void SetDz(Precision const &_z) { fCutTube.fTubeStruct.SetDz(_z); }

  VECGEOM_FORCE_INLINE
  void SetSPhi(Precision const &_sphi) { fCutTube.fTubeStruct.SetAndCheckSPhiAngle(_sphi); }

  VECGEOM_FORCE_INLINE
  void SetDPhi(Precision const &_dphi) { fCutTube.fTubeStruct.SetAndCheckDPhiAngle(_dphi); }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  evolution::Wedge const &GetWedge() const { return fCutTube.fTubeStruct.fPhiWedge; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Vector3D<Precision> BottomNormal() const { return fCutTube.fCutPlanes.GetNormal(0); }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Vector3D<Precision> TopNormal() const { return fCutTube.fCutPlanes.GetNormal(1); }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision ZlimitBottom(Precision rVal, Precision phiVal) const
  {
    return (-z() -
            (rVal / BottomNormal().z()) *
                (BottomNormal().x() * vecCore::math::Cos(phiVal) + BottomNormal().y() * vecCore::math::Sin(phiVal)));
  }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision ZlimitTop(Precision rVal, Precision phiVal) const
  {
    return (z() -
            (rVal / TopNormal().z()) *
                (TopNormal().x() * vecCore::math::Cos(phiVal) + TopNormal().y() * vecCore::math::Sin(phiVal)));
  }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetLateralArea(Precision rVal) const
  {
    return (2. * rVal * z() * dphi() -
            rVal * rVal * (((TopNormal().x() / TopNormal().z() - BottomNormal().x() / BottomNormal().z()) *
                            (fCutTube.fSinPhi2 - fCutTube.fSinPhi1)) -
                           ((TopNormal().y() / TopNormal().z() - BottomNormal().y() / BottomNormal().z()) *
                            (fCutTube.fCosPhi2 - fCutTube.fCosPhi1))));
  }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetTopArea() const
  {
    return 0.5 * (fCutTube.fTubeStruct.fRmax2 - fCutTube.fTubeStruct.fRmin2) * dphi() / TopNormal().z();
  }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetBottomArea() const
  {
    return -0.5 * (fCutTube.fTubeStruct.fRmax2 - fCutTube.fTubeStruct.fRmin2) * dphi() / BottomNormal().z();
  }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetLateralPhi1Area() const
  {
    if (dphi() == kTwoPi) return 0.;
    return (0.5 * (rmax() - rmin()) * ((ZlimitTop(rmin(), sphi()) - ZlimitBottom(rmin(), sphi())) +
                                       (ZlimitTop(rmax(), sphi()) - ZlimitBottom(rmax(), sphi()))));
  }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetLateralPhi2Area() const
  {
    if (dphi() == kTwoPi) return 0.;
    return (0.5 * (rmax() - rmin()) * ((ZlimitTop(rmin(), sphi() + dphi()) - ZlimitBottom(rmin(), sphi() + dphi())) +
                                       (ZlimitTop(rmax(), sphi() + dphi()) - ZlimitBottom(rmax(), sphi() + dphi()))));
  }

  VECCORE_ATT_HOST_DEVICE
  void Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const override;

  Vector3D<Precision> SamplePointOnSurface() const override;

  Precision volume() const;

  Precision Capacity() const { return volume(); }

  VECCORE_ATT_HOST_DEVICE
  Precision SurfaceArea() const
  {
    return (GetBottomArea() + GetTopArea() + GetLateralArea(rmin()) + GetLateralArea(rmax()) + GetLateralPhi1Area() +
            GetLateralPhi2Area());
  }

  VECCORE_ATT_HOST_DEVICE
  bool Normal(Vector3D<Precision> const &point, Vector3D<Precision> &normal) const override;

  VECCORE_ATT_HOST_DEVICE
  virtual void Print() const final;
  virtual void Print(std::ostream &os) const final;

  std::string GetEntityType() const { return "CutTube"; }

  virtual int MemorySize() const final { return sizeof(*this); }

#ifdef VECGEOM_CUDA_INTERFACE
  size_t DeviceSizeOf() const final { return DevicePtr<cuda::UnplacedCutTube>::SizeOf(); }
  DevicePtr<cuda::VUnplacedVolume> CopyToGpu() const final;
  DevicePtr<cuda::VUnplacedVolume> CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const gpu_ptr) const final;
#endif

  template <TranslationCode transCodeT, RotationCode rotCodeT>
  VECCORE_ATT_DEVICE
  static VPlacedVolume *Create(LogicalVolume const *const logical_volume, Transformation3D const *const transformation,
#ifdef VECCORE_CUDA
                               const int id,
#endif
                               VPlacedVolume *const placement = NULL);

private:
  VECCORE_ATT_DEVICE
  virtual VPlacedVolume *SpecializedVolume(LogicalVolume const *const volume,
                                           Transformation3D const *const transformation,
                                           const TranslationCode trans_code, const RotationCode rot_code,
#ifdef VECCORE_CUDA
                                           const int id,
#endif
                                           VPlacedVolume *const placement = NULL) const final;
};
} // end inline namespace
} // end vecgeom namespace

#endif // VECGEOM_VOLUMES_UNPLACEDCUTTUBE_H_
