#ifndef UNPLACEDBOOLEANVOLUME_H_
#define UNPLACEDBOOLEANVOLUME_H_

#include "base/Global.h"
#include "base/AlignedBase.h"
#include "base/Vector3D.h"
#include "volumes/UnplacedVolume.h"
#include "volumes/PlacedVolume.h"
#include "volumes/BooleanStruct.h"
#include "volumes/kernel/BooleanImplementation.h"
#include "volumes/UnplacedVolumeImplHelper.h"

namespace vecgeom {

VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_1v(class, UnplacedBooleanVolume, BooleanOperation, Arg1);

inline namespace VECGEOM_IMPL_NAMESPACE {

/**
 * A class representing a simple UNPLACED boolean volume A-B
 * It takes two template arguments:
 * 1.: the mother (or left) volume A in unplaced form
 * 2.: the (or right) volume B in placed form, acting on A with a boolean operation;
 * the placement is with respect to the left volume
 *
 *
 *
 * will be a boolean solid where two boxes are subtracted
 * and B is only translated (not rotated) with respect to A
 *
 */
template <BooleanOperation Op>
class UnplacedBooleanVolume : public LoopUnplacedVolumeImplHelper<BooleanImplementation<Op>>, public AlignedBase {

public:
  BooleanStruct fBoolean;
  using LoopUnplacedVolumeImplHelper<BooleanImplementation<Op>>::fGlobalConvexity;

  // the constructor
  VECCORE_ATT_HOST_DEVICE
  UnplacedBooleanVolume(BooleanOperation op, VPlacedVolume const *left, VPlacedVolume const *right)
      : fBoolean(op, left, right)
  {
    fGlobalConvexity = false;
#ifndef VECCORE_CUDA
    if (fBoolean.fLeftVolume->IsAssembly() || fBoolean.fRightVolume->IsAssembly()) {
      throw std::runtime_error("Trying to make boolean out of assembly which is not supported\n");
    }
#endif
  }

  virtual int MemorySize() const override { return sizeof(*this); }

#ifdef VECGEOM_CUDA_INTERFACE
  virtual size_t DeviceSizeOf() const override { return DevicePtr<cuda::UnplacedBooleanVolume<Op>>::SizeOf(); }
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu() const override;
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const gpu_ptr) const override;
#endif

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  BooleanOperation GetOp() const { return fBoolean.fOp; }

  VECCORE_ATT_HOST_DEVICE
  BooleanStruct const &GetStruct() const { return fBoolean; }

  VECCORE_ATT_HOST_DEVICE
  bool Normal(Vector3D<Precision> const &point, Vector3D<Precision> &normal) const override;

  Precision Capacity() const override
  {
    // TBDONE -- need some sampling
    return 0.;
  }
#if !defined(VECCORE_CUDA)
  VECGEOM_FORCE_INLINE
  Precision SurfaceArea() const
  {
    // TBDONE -- need some sampling
    return 0.;
  }
#endif // !VECCORE_CUDA

  VECCORE_ATT_HOST_DEVICE
  void Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const override;

  Vector3D<Precision> SamplePointOnSurface() const override
  {
    // TBDONE
    return Vector3D<Precision>();
  }

  std::string GetEntityType() const { return "BooleanVolume"; }

  VECCORE_ATT_HOST_DEVICE
  virtual void Print() const override{};

  virtual void Print(std::ostream & /*os*/) const override{};

  template <TranslationCode transCodeT, RotationCode rotCodeT>
  VECCORE_ATT_DEVICE
  static VPlacedVolume *Create(LogicalVolume const *const logical_volume, Transformation3D const *const transformation,
#ifdef VECCORE_CUDA
                               const int id,
#endif
                               VPlacedVolume *const placement = NULL);

  VPlacedVolume const *GetLeft() const { return fBoolean.fLeftVolume; }
  VPlacedVolume const *GetRight() const { return fBoolean.fRightVolume; }

private:
  VECCORE_ATT_DEVICE
  virtual VPlacedVolume *SpecializedVolume(LogicalVolume const *const volume,
                                           Transformation3D const *const transformation,
                                           const TranslationCode trans_code, const RotationCode rot_code,
#ifdef VECCORE_CUDA
                                           const int id,
#endif
                                           VPlacedVolume *const placement = NULL) const override;

  void SetLeft(VPlacedVolume const *pvol) { fBoolean.fLeftVolume = pvol; }
  void SetRight(VPlacedVolume const *pvol) { fBoolean.fRightVolume = pvol; }

  friend class GeoManager;
}; // End class
} // End impl namespace

} // End global namespace

#endif /* UNPLACEDBOOLEANVOLUME_H_ */
