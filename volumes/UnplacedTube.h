#ifndef VECGEOM_VOLUMES_UNPLACEDTUBE_H_
#define VECGEOM_VOLUMES_UNPLACEDTUBE_H_

#include "base/Global.h"
#include "base/RNG.h"
#include "base/AlignedBase.h"
#include "base/Array.h"
#include "volumes/UnplacedVolume.h"
#include "volumes/TubeStruct.h"
#include "volumes/kernel/TubeImplementation.h"
#include "volumes/Wedge.h"
#include "volumes/UnplacedVolumeImplHelper.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class UnplacedTube;);
VECGEOM_DEVICE_DECLARE_CONV(class, UnplacedTube);
VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE(class, SUnplacedTube, typename);

inline namespace VECGEOM_IMPL_NAMESPACE {

// Introduce Intermediate class ( so that we can do typecasting )
class UnplacedTube : public VUnplacedVolume {
private:
  // tube defining parameters
  TubeStruct<Precision> fTube;

public:
  VECCORE_ATT_HOST_DEVICE
  UnplacedTube(Precision const &_rmin, Precision const &_rmax, Precision const &_z, Precision const &_sphi,
               Precision const &_dphi)
      : fTube(_rmin, _rmax, _z, _sphi, _dphi)
  {
    DetectConvexity();
  }

  VECCORE_ATT_HOST_DEVICE
  UnplacedTube(UnplacedTube const &other) = delete;

  VECCORE_ATT_HOST_DEVICE
  TubeStruct<double> const &GetStruct() const { return fTube; }

  VECCORE_ATT_HOST_DEVICE
  void DetectConvexity();

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision rmin() const { return fTube.fRmin; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision rmax() const { return fTube.fRmax; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision z() const { return fTube.fZ; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision sphi() const { return fTube.fSphi; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision dphi() const { return fTube.fDphi; }

  VECGEOM_FORCE_INLINE
  void SetRMin(Precision const &_rmin) { fTube.SetRMin(_rmin); }

  VECGEOM_FORCE_INLINE
  void SetRMax(Precision const &_rmax) { fTube.SetRMax(_rmax); }

  VECGEOM_FORCE_INLINE
  void SetDz(Precision const &_z) { fTube.SetDz(_z); }

  VECGEOM_FORCE_INLINE
  void SetSPhi(Precision const &_sphi) { fTube.SetAndCheckSPhiAngle(_sphi); /*weird name*/ }

  VECGEOM_FORCE_INLINE
  void SetDPhi(Precision const &_dphi) { fTube.SetAndCheckDPhiAngle(_dphi); /*weird name*/ }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  evolution::Wedge const &GetWedge() const { return fTube.fPhiWedge; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision volume() const { return fTube.fZ * (fTube.fRmax2 - fTube.fRmin2) * fTube.fDphi; }

  VECCORE_ATT_HOST_DEVICE
  void Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const override;

  Vector3D<Precision> GetPointOnSurface() const override;

  // VECCORE_ATT_HOST_DEVICE
  Precision Capacity() const { return volume(); }

  // VECCORE_ATT_HOST_DEVICE
  Precision SurfaceArea() const
  {
    return GetTopArea() + GetLateralPhiArea() + GetLateralROutArea() + GetLateralRInArea();
  }

  // VECCORE_ATT_HOST_DEVICE
  Precision GetTopArea() const
  { // Abhijit:: this is top and bottom circular area of tube
    return 2 * 0.5 * (fTube.fRmax2 - fTube.fRmin2) * fTube.fDphi;
  }

  // VECCORE_ATT_HOST_DEVICE
  Precision GetLateralPhiArea() const
  { // Abhijit:: this is vertical Phi_start and phi_end opening
    // factor of 2 since fZ is half length
    return (fTube.fDphi < kTwoPi) ? 4. * fTube.fZ * (fTube.fRmax - fTube.fRmin) : 0.;
  }

  // VECCORE_ATT_HOST_DEVICE
  Precision GetLateralRInArea() const
  { // Abhijit:: this is Inner surface of tube along Z
    // factor of 2 since fZ is half length
    return 2. * fTube.fZ * fTube.fRmin * fTube.fDphi;
  }

  // VECCORE_ATT_HOST_DEVICE
  Precision GetLateralROutArea() const
  { // Abhijit:: this is Outer surface of tube along Z
    // factor of 2 since fZ is half length
    return 2. * fTube.fZ * fTube.fRmax * fTube.fDphi;
  }

  //  This computes where the random point would be placed
  // 1::rTop, 2::rBot, 3::phiLeft, 4::phiRight, 5::zIn, 6::zOut
  // VECCORE_ATT_HOST_DEVICE
  int ChooseSurface() const;

  VECCORE_ATT_HOST_DEVICE
  bool Normal(Vector3D<Precision> const &point, Vector3D<Precision> &normal) const override;

  VECCORE_ATT_HOST_DEVICE
  virtual void Print() const override;

  std::string GetEntityType() const { return "Tube"; }

  template <TranslationCode transCodeT, RotationCode rotCodeT>
  VECCORE_ATT_DEVICE
  static VPlacedVolume *Create(LogicalVolume const *const logical_volume, Transformation3D const *const transformation,
#ifdef VECGEOM_NVCC
                               const int id,
#endif
                               VPlacedVolume *const placement = NULL);

#ifdef VECGEOM_CUDA_INTERFACE
  virtual size_t DeviceSizeOf() const override
  {
    return DevicePtr<cuda::SUnplacedTube<cuda::TubeTypes::UniversalTube>>::SizeOf();
  }
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu() const override;
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const gpu_ptr) const override;
#endif

private:
  virtual void Print(std::ostream &os) const override;

#ifndef VECGEOM_NVCC
  virtual VPlacedVolume *SpecializedVolume(LogicalVolume const *const volume,
                                           Transformation3D const *const transformation,
                                           const TranslationCode trans_code, const RotationCode rot_code,
                                           VPlacedVolume *const placement = NULL) const override;

#else
  __device__ virtual VPlacedVolume *SpecializedVolume(LogicalVolume const *const volume,
                                                      Transformation3D const *const transformation,
                                                      const TranslationCode trans_code, const RotationCode rot_code,
                                                      const int id,
                                                      VPlacedVolume *const placement = NULL) const override;

#endif
};

// this class finishes the implementation

template <typename TubeType = TubeTypes::UniversalTube>
class SUnplacedTube : public SIMDUnplacedVolumeImplHelper<TubeImplementation<TubeType>, UnplacedTube>,
                      public AlignedBase {
public:
  using BaseType_t = SIMDUnplacedVolumeImplHelper<TubeImplementation<TubeType>, UnplacedTube>;
  using BaseType_t::BaseType_t;
};

using GenericUnplacedTube = SUnplacedTube<TubeTypes::UniversalTube>;

} // end inline namespace
} // end vecgeom namespace

#endif // VECGEOM_VOLUMES_UNPLACEDTUBE_H_
