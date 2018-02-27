/*
 * UnplacedCone.h
 *
 *  Created on: May 14, 2014
 *      Author: swenzel
 */

#ifndef VECGEOM_VOLUMES_UNPLACEDCONE_H_
#define VECGEOM_VOLUMES_UNPLACEDCONE_H_

#include "base/Global.h"
#include "base/AlignedBase.h"
#include "volumes/UnplacedVolume.h"
#include "volumes/ConeStruct.h"
#include "volumes/kernel/ConeImplementation.h"
#include "volumes/UnplacedVolumeImplHelper.h"
#include "volumes/kernel/shapetypes/ConeTypes.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class UnplacedCone;);
VECGEOM_DEVICE_DECLARE_CONV(class, UnplacedCone);

inline namespace VECGEOM_IMPL_NAMESPACE {

/**
 * Class representing an unplaced cone; Encapsulated parameters of a cone and
 * functions that do not depend on how the cone is placed in a reference frame
 *
 * The unplaced cone is represented by the following parameters
 *
 * Member Data:
 *
 * fCone.fDz half length in z direction;  ( the cone has height 2*fDz )
 * fCone.fRmin1  inside radius at  -fDz ( in internal coordinate system )
 * fCone.fRmin2  inside radius at  +fDz
 * fCone.fRmax1  outside radius at -fDz
 * fCone.fRmax2  outside radius at +fDz
 * fCone.fSPhi starting angle of the segment in radians
 * fCone.fDPhi delta angle of the segment in radians
 */
class UnplacedCone : public SIMDUnplacedVolumeImplHelper<ConeImplementation<ConeTypes::UniversalCone>>,
                     public AlignedBase {
private:
  ConeStruct<Precision> fCone;

public:
  VECCORE_ATT_HOST_DEVICE
  UnplacedCone(Precision rmin1, Precision rmax1, Precision rmin2, Precision rmax2, Precision dz, Precision phimin,
               Precision deltaphi)
      : fCone(rmin1, rmax1, rmin2, rmax2, dz, phimin, deltaphi)
  {

    DetectConvexity();
  }

  VECCORE_ATT_HOST_DEVICE
  ConeStruct<double> const &GetStruct() const { return fCone; }

  VECCORE_ATT_HOST_DEVICE
  Precision GetInvSecRMax() const { return fCone.fInvSecRMax; }

  VECCORE_ATT_HOST_DEVICE
  Precision GetInvSecRMin() const { return fCone.fInvSecRMin; }

  VECCORE_ATT_HOST_DEVICE
  Precision GetTolIz() const { return fCone.fTolIz; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetTolOz() const { return fCone.fTolOz; }

  VECCORE_ATT_HOST_DEVICE
  Precision GetConeTolerane() const { return fCone.fConeTolerance; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetSqRmin1() const { return fCone.fSqRmin1; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetSqRmin2() const { return fCone.fSqRmin2; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetSqRmax1() const { return fCone.fSqRmax1; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetSqRmax2() const { return fCone.fSqRmax2; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetTanRmax() const { return fCone.fTanRMax; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetTanRmin() const { return fCone.fTanRMin; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetSecRmax() const { return fCone.fSecRMax; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetSecRmin() const { return fCone.fSecRMin; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetZNormInner() const { return fCone.fZNormInner; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetZNormOuter() const { return fCone.fZNormOuter; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetInnerConeApex() const { return fCone.fInnerConeApex; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetTInner() const { return fCone.fTanInnerApexAngle; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetOuterConeApex() const { return fCone.fOuterConeApex; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetTOuter() const { return fCone.fTanOuterApexAngle; }

  VECCORE_ATT_HOST_DEVICE
  void DetectConvexity();
  VECCORE_ATT_HOST_DEVICE
  Precision GetRmin1() const { return fCone.fRmin1; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetRmax1() const { return fCone.fRmax1; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetRmin2() const { return fCone.fRmin2; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetRmax2() const { return fCone.fRmax2; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetDz() const { return fCone.fDz; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetSPhi() const { return fCone.fSPhi; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetDPhi() const { return fCone.fDPhi; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetInnerSlope() const { return fCone.fInnerSlope; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetOuterSlope() const { return fCone.fOuterSlope; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetInnerOffset() const { return fCone.fInnerOffset; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetOuterOffset() const { return fCone.fOuterOffset; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetAlongPhi1X() const { return fCone.fAlongPhi1x; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetAlongPhi1Y() const { return fCone.fAlongPhi1y; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetAlongPhi2X() const { return fCone.fAlongPhi2x; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetAlongPhi2Y() const { return fCone.fAlongPhi2y; }
  VECCORE_ATT_HOST_DEVICE
  evolution::Wedge const &GetWedge() const { return fCone.fPhiWedge; }

  VECCORE_ATT_HOST_DEVICE
  void SetAndCheckSPhiAngle(Precision sPhi);

  VECCORE_ATT_HOST_DEVICE
  void SetAndCheckDPhiAngle(Precision dPhi);

  void SetRmin1(Precision const &arg)
  {
    fCone.fRmin1 = arg;
    fCone.CalculateCached();
  }
  void SetRmax1(Precision const &arg)
  {
    fCone.fRmax1 = arg;
    fCone.CalculateCached();
  }
  void SetRmin2(Precision const &arg)
  {
    fCone.fRmin2 = arg;
    fCone.CalculateCached();
  }
  void SetRmax2(Precision const &arg)
  {
    fCone.fRmax2 = arg;
    fCone.CalculateCached();
  }
  void SetDz(Precision const &arg)
  {
    fCone.fDz = arg;
    fCone.CalculateCached();
  }
  void SetSPhi(Precision const &arg)
  {
    fCone.fSPhi = arg;
    fCone.SetAndCheckSPhiAngle(fCone.fSPhi);
    DetectConvexity();
  }
  void SetDPhi(Precision const &arg)
  {
    fCone.fDPhi = arg;
    fCone.SetAndCheckDPhiAngle(fCone.fDPhi);
    DetectConvexity();
  }

  VECCORE_ATT_HOST_DEVICE
  bool IsFullPhi() const { return fCone.fDPhi == kTwoPi; }

  virtual int MemorySize() const final { return sizeof(*this); }

  VECCORE_ATT_HOST_DEVICE
  virtual void Print() const final;
  virtual void Print(std::ostream &os) const final;

  std::string GetEntityType() const { return "Cone"; }
  std::ostream &StreamInfo(std::ostream &os) const;

  VECCORE_ATT_DEVICE
  virtual VPlacedVolume *SpecializedVolume(LogicalVolume const *const volume,
                                           Transformation3D const *const transformation,
                                           const TranslationCode trans_code, const RotationCode rot_code,
#ifdef VECCORE_CUDA
                                           const int id,
#endif
                                           VPlacedVolume *const placement = NULL) const final;

  template <TranslationCode transCodeT, RotationCode rotCodeT>
  VECCORE_ATT_DEVICE
  static VPlacedVolume *Create(LogicalVolume const *const logical_volume, Transformation3D const *const transformation,
#ifdef VECCORE_CUDA
                               const int id,
#endif
                               VPlacedVolume *const placement = NULL);

#ifdef VECGEOM_CUDA_INTERFACE
  virtual size_t DeviceSizeOf() const { return DevicePtr<cuda::UnplacedCone>::SizeOf(); }
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu() const;
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const gpu_ptr) const;
#endif

#ifndef VECCORE_CUDA
  Precision Capacity() const
  {
    return (fCone.fDz * fCone.fDPhi / 3.) *
           (fCone.fRmax1 * fCone.fRmax1 + fCone.fRmax2 * fCone.fRmax2 + fCone.fRmax1 * fCone.fRmax2 -
            fCone.fRmin1 * fCone.fRmin1 - fCone.fRmin2 * fCone.fRmin2 - fCone.fRmin1 * fCone.fRmin2);
  }

  Precision SurfaceArea() const
  {
    double mmin, mmax, dmin, dmax;
    mmin = (fCone.fRmin1 + fCone.fRmin2) * 0.5;
    mmax = (fCone.fRmax1 + fCone.fRmax2) * 0.5;
    dmin = (fCone.fRmin2 - fCone.fRmin1);
    dmax = (fCone.fRmax2 - fCone.fRmax1);

    return fCone.fDPhi * (mmin * std::sqrt(dmin * dmin + 4 * fCone.fDz * fCone.fDz) +
                          mmax * std::sqrt(dmax * dmax + 4 * fCone.fDz * fCone.fDz) +
                          0.5 * (fCone.fRmax1 * fCone.fRmax1 - fCone.fRmin1 * fCone.fRmin1 +
                                 fCone.fRmax2 * fCone.fRmax2 - fCone.fRmin2 * fCone.fRmin2));
  }

  void Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const;

  bool Normal(Vector3D<Precision> const &point, Vector3D<Precision> &normal) const;
  Vector3D<Precision> SamplePointOnSurface() const;

  // Helper funtion to detect edge points
  template <bool top>
  bool IsOnZPlane(Vector3D<Precision> const &point) const;
  template <bool start>
  bool IsOnPhiWedge(Vector3D<Precision> const &point) const;
  template <bool inner>
  bool IsOnConicalSurface(Vector3D<Precision> const &point) const;
  template <bool inner>
  Precision GetRadiusOfConeAtPoint(Precision const pointZ) const;

  bool IsOnEdge(Vector3D<Precision> &point) const;

#endif // !VECCORE_CUDA
};
}
} // End global namespace

#endif
