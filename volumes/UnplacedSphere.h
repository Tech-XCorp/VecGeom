/// \file UnplacedSphere.h
/// \author Raman Sehgal (raman.sehgal@cern.ch)

#ifndef VECGEOM_VOLUMES_UNPLACEDSPHERE_H_
#define VECGEOM_VOLUMES_UNPLACEDSPHERE_H_

#include "base/Global.h"
#include "base/AlignedBase.h"
#include "volumes/UnplacedVolume.h"
#include "volumes/SphereStruct.h"
#include "volumes/kernel/SphereImplementation.h"
#include "volumes/UnplacedVolumeImplHelper.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class UnplacedSphere;);
VECGEOM_DEVICE_DECLARE_CONV(class, UnplacedSphere);

inline namespace VECGEOM_IMPL_NAMESPACE {

class UnplacedSphere : public SIMDUnplacedVolumeImplHelper<SphereImplementation>, public AlignedBase {

private:
  SphereStruct<Precision> fSphere;

public:
  VECCORE_ATT_HOST_DEVICE
  SphereStruct<double> const &GetStruct() const { return fSphere; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  evolution::Wedge const &GetWedge() const { return fSphere.fPhiWedge; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  ThetaCone const &GetThetaCone() const { return fSphere.fThetaCone; }

  VECCORE_ATT_HOST_DEVICE
  UnplacedSphere(Precision pRmin, Precision pRmax, Precision pSPhi, Precision pDPhi, Precision pSTheta,
                 Precision pDTheta);

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetInsideRadius() const { return fSphere.fRmin; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetInnerRadius() const { return fSphere.fRmin; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetOuterRadius() const { return fSphere.fRmax; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetStartPhiAngle() const { return fSphere.fSPhi; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetDeltaPhiAngle() const { return fSphere.fDPhi; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetStartThetaAngle() const { return fSphere.fSTheta; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetDeltaThetaAngle() const { return fSphere.fDTheta; }

  // Functions to get Tolerance
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetFRminTolerance() const { return fSphere.fRminTolerance; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetMKTolerance() const { return fSphere.mkTolerance; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetAngTolerance() const { return fSphere.kAngTolerance; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  bool IsFullSphere() const { return fSphere.fFullSphere; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  bool IsFullPhiSphere() const { return fSphere.fFullPhiSphere; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  bool IsFullThetaSphere() const { return fSphere.fFullThetaSphere; }

  // All angle related functions
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetHDPhi() const { return fSphere.hDPhi; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetCPhi() const { return fSphere.cPhi; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetEPhi() const { return fSphere.ePhi; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetSinCPhi() const { return fSphere.sinCPhi; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetCosCPhi() const { return fSphere.cosCPhi; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetSinSPhi() const { return fSphere.sinSPhi; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetCosSPhi() const { return fSphere.cosSPhi; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetSinEPhi() const { return fSphere.sinEPhi; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetCosEPhi() const { return fSphere.cosEPhi; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetETheta() const { return fSphere.eTheta; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetSinSTheta() const { return fSphere.sinSTheta; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetCosSTheta() const { return fSphere.cosSTheta; }

  //****************************************************************
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetTanSTheta() const { return fSphere.tanSTheta; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetTanETheta() const { return fSphere.tanETheta; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetFabsTanSTheta() const { return fSphere.fabsTanSTheta; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetFabsTanETheta() const { return fSphere.fabsTanETheta; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetTanSTheta2() const { return fSphere.tanSTheta2; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetTanETheta2() const { return fSphere.tanETheta2; }
  //****************************************************************

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetSinETheta() const { return fSphere.sinETheta; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetCosETheta() const { return fSphere.cosETheta; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetCosHDPhiOT() const { return fSphere.cosHDPhiOT; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetCosHDPhiIT() const { return fSphere.cosHDPhiIT; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  void SetInsideRadius(Precision newRmin) { fSphere.SetInsideRadius(newRmin); }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  void SetInnerRadius(Precision newRmin) { SetInsideRadius(newRmin); }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  void SetOuterRadius(Precision newRmax) { fSphere.SetOuterRadius(newRmax); }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  void SetStartPhiAngle(Precision newSPhi, bool compute = true) { fSphere.SetStartPhiAngle(newSPhi, compute); }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  void SetDeltaPhiAngle(Precision newDPhi) { fSphere.SetDeltaPhiAngle(newDPhi); }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  void SetStartThetaAngle(Precision newSTheta) { fSphere.SetStartThetaAngle(newSTheta); }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  void SetDeltaThetaAngle(Precision newDTheta) { fSphere.SetDeltaThetaAngle(newDTheta); }

  // Old access functions
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetRmin() const { return fSphere.fRmin; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetRmax() const { return fSphere.fRmax; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetSPhi() const { return fSphere.fSPhi; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetDPhi() const { return fSphere.fDPhi; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetSTheta() const { return fSphere.fSTheta; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetDTheta() const { return fSphere.fDTheta; }

  VECCORE_ATT_HOST_DEVICE
  void CalcCapacity();

  VECCORE_ATT_HOST_DEVICE
  void CalcSurfaceArea();

  VECCORE_ATT_HOST_DEVICE
  void DetectConvexity();

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision Capacity() const { return fSphere.fCubicVolume; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision SurfaceArea() const { return fSphere.fSurfaceArea; }

#ifndef VECCORE_CUDA
  void Extent(Vector3D<Precision> &, Vector3D<Precision> &) const;
  bool Normal(Vector3D<Precision> const &point, Vector3D<Precision> &normal) const
  {

    bool valid;
    normal = SphereImplementation::Normal<Precision>(fSphere, point, valid);
    return valid;
  }

  Vector3D<Precision> SamplePointOnSurface() const;

  std::string GetEntityType() const;
#endif

  void GetParametersList(int aNumber, Precision *aArray) const;

#if defined(VECGEOM_USOLIDS)
  std::ostream &StreamInfo(std::ostream &os) const;
#endif

  VECCORE_ATT_HOST_DEVICE
  void ComputeBBox() const;

  // VECCORE_ATT_HOST_DEVICE
  // Precision sqr(Precision x) {return x*x;};

public:
  virtual int MemorySize() const final { return sizeof(*this); }

  VECCORE_ATT_HOST_DEVICE
  virtual void Print() const final;

  // VECCORE_ATT_HOST_DEVICE
  virtual void Print(std::ostream &os) const final;

#ifndef VECCORE_CUDA

  template <TranslationCode trans_code, RotationCode rot_code>
  static VPlacedVolume *Create(LogicalVolume const *const logical_volume, Transformation3D const *const transformation,
                               VPlacedVolume *const placement = NULL);

  static VPlacedVolume *CreateSpecializedVolume(LogicalVolume const *const volume,
                                                Transformation3D const *const transformation,
                                                const TranslationCode trans_code, const RotationCode rot_code,
                                                VPlacedVolume *const placement = NULL);

#else

  template <TranslationCode trans_code, RotationCode rot_code>
  __device__
  static VPlacedVolume *Create(LogicalVolume const *const logical_volume, Transformation3D const *const transformation,
                               const int id, VPlacedVolume *const placement = NULL);

  __device__ static VPlacedVolume *CreateSpecializedVolume(LogicalVolume const *const volume,
                                                           Transformation3D const *const transformation,
                                                           const TranslationCode trans_code,
                                                           const RotationCode rot_code, const int id,
                                                           VPlacedVolume *const placement = NULL);

#endif

#ifdef VECGEOM_CUDA_INTERFACE
  virtual size_t DeviceSizeOf() const { return DevicePtr<cuda::UnplacedSphere>::SizeOf(); }
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu() const;
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const gpu_ptr) const;
#endif

private:
#ifndef VECCORE_CUDA

  virtual VPlacedVolume *SpecializedVolume(LogicalVolume const *const volume,
                                           Transformation3D const *const transformation,
                                           const TranslationCode trans_code, const RotationCode rot_code,
                                           VPlacedVolume *const placement = NULL) const final
  {
    return CreateSpecializedVolume(volume, transformation, trans_code, rot_code, placement);
  }

#else

  __device__ virtual VPlacedVolume *SpecializedVolume(LogicalVolume const *const volume,
                                                      Transformation3D const *const transformation,
                                                      const TranslationCode trans_code, const RotationCode rot_code,
                                                      const int id, VPlacedVolume *const placement = NULL) const final
  {
    return CreateSpecializedVolume(volume, transformation, trans_code, rot_code, id, placement);
  }

#endif
};
}
} // End global namespace

#endif // VECGEOM_VOLUMES_UNPLACEDSPHERE_H_
