/// @file UnplacedTessellated.h
/// @author Mihaela Gheata (mihaela.gheata@cern.ch)

#ifndef VECGEOM_VOLUMES_UNPLACEDTESSELLATED_H_
#define VECGEOM_VOLUMES_UNPLACEDTESSELLATED_H_

#include "base/Global.h"
#include "base/AlignedBase.h"
#include "volumes/UnplacedVolume.h"
#include "TessellatedStruct.h"
#include "volumes/kernel/TessellatedImplementation.h"
#include "volumes/UnplacedVolumeImplHelper.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class UnplacedTessellated;);
VECGEOM_DEVICE_DECLARE_CONV(class, UnplacedTessellated);

inline namespace VECGEOM_IMPL_NAMESPACE {

class UnplacedTessellated : public LoopUnplacedVolumeImplHelper<TessellatedImplementation>, public AlignedBase {
protected:
  mutable TessellatedStruct<3, double> fTessellated; ///< Structure with Tessellated parameters

public:
  VECCORE_ATT_HOST_DEVICE
  UnplacedTessellated() : fTessellated() { fGlobalConvexity = false; }

  VECCORE_ATT_HOST_DEVICE
  TessellatedStruct<3, double> const &GetStruct() const { return fTessellated; }

  VECCORE_ATT_HOST_DEVICE
  bool AddTriangularFacet(Vector3D<double> const &vt0, Vector3D<double> const &vt1, Vector3D<double> const &vt2,
                          bool absolute = true)
  {
    return fTessellated.AddTriangularFacet(vt0, vt1, vt2, absolute);
  }

  VECCORE_ATT_HOST_DEVICE
  bool AddQuadrilateralFacet(Vector3D<double> const &vt0, Vector3D<double> const &vt1, Vector3D<double> const &vt2,
                             Vector3D<double> const &vt3, bool absolute = true)
  {
    return fTessellated.AddQuadrilateralFacet(vt0, vt1, vt2, vt3, absolute);
  }

  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  size_t GetNFacets() const { return fTessellated.fFacets.size(); }

  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  TriangleFacet<double> *GetFacet(int ifacet) const { return fTessellated.fFacets[ifacet]; }

  VECCORE_ATT_HOST_DEVICE
  void Close()
  {
    fTessellated.Close();
    PrecomputeInsideAndSafety();
  }

  VECCORE_ATT_HOST_DEVICE
  bool IsClosed() const { return fTessellated.fSolidClosed; }

  virtual int memory_size() const { return sizeof(*this); }

  VECCORE_ATT_HOST_DEVICE
  void Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const override { fTessellated.Extent(aMin, aMax); }

  // Computes capacity of the shape in [length^3]
  // VECCORE_ATT_HOST_DEVICE
  Precision Capacity() const override;

  // VECCORE_ATT_HOST_DEVICE
  Precision SurfaceArea() const override;

  VECCORE_ATT_HOST_DEVICE
  int ChooseSurface() const;

  Vector3D<Precision> SamplePointOnSurface() const override;

  VECCORE_ATT_HOST_DEVICE
  bool Normal(Vector3D<Precision> const &point, Vector3D<Precision> &normal) const override;

  VECCORE_ATT_HOST_DEVICE
  void PrecomputeInsideAndSafety();

  VECCORE_ATT_HOST_DEVICE
  virtual void Print() const override;

  std::string GetEntityType() const { return "Tessellated"; }

  template <TranslationCode transCodeT, RotationCode rotCodeT>
  VECCORE_ATT_DEVICE
  static VPlacedVolume *Create(LogicalVolume const *const logical_volume, Transformation3D const *const transformation,
#ifdef VECCORE_CUDA
                               const int id,
#endif
                               VPlacedVolume *const placement = NULL);

#ifdef VECGEOM_CUDA_INTERFACE
  virtual size_t DeviceSizeOf() const override { return DevicePtr<cuda::UnplacedTessellated>::SizeOf(); }
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu() const override;
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const gpu_ptr) const override;
#endif

  std::ostream &StreamInfo(std::ostream &os) const;

  virtual void Print(std::ostream &os) const override;

private:
  VECCORE_ATT_DEVICE
  virtual VPlacedVolume *SpecializedVolume(LogicalVolume const *const volume,
                                           Transformation3D const *const transformation,
                                           const TranslationCode trans_code, const RotationCode rot_code,
#ifdef VECCORE_CUDA
                                           const int id,
#endif
                                           VPlacedVolume *const placement = NULL) const override;
};
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_UNPLACEDTESSELLATED_H_
