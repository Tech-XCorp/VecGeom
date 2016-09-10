/// @file UnplacedScaledShape.h
/// @author Mihaela Gheata (mihaela.gheata@cern.ch)

#ifndef VECGEOM_VOLUMES_UNPLACEDSCALEDSHAPE_H_
#define VECGEOM_VOLUMES_UNPLACEDSCALEDSHAPE_H_

#include "base/Global.h"
#include "base/AlignedBase.h"
#include "base/Vector3D.h"
#include "base/Scale3D.h"
#include "ScaledShapeStruct.h"
#include "volumes/UnplacedVolume.h"
#include "volumes/PlacedVolume.h"
#include "volumes/LogicalVolume.h"
#include "volumes/kernel/ScaledShapeImplementation.h"
#include "volumes/UnplacedVolumeImplHelper.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class UnplacedScaledShape;);
VECGEOM_DEVICE_DECLARE_CONV(class, UnplacedScaledShape);

inline namespace VECGEOM_IMPL_NAMESPACE {

/**
 * The unplaced scaled shape class.
 */
class UnplacedScaledShape : public SIMDUnplacedVolumeImplHelper<ScaledShapeImplementation>, public AlignedBase {

public:
  ScaledShapeStruct<double> fScaled; /* The scaled shape structure */

public:
  /// Dummy ctor
  VECGEOM_CUDA_HEADER_BOTH
  UnplacedScaledShape() : fScaled() { fGlobalConvexity = fScaled.fPlaced->GetUnplacedVolume()->IsConvex(); }

  /// Constructor based on placed volume
  VECGEOM_CUDA_HEADER_BOTH
  UnplacedScaledShape(VPlacedVolume const *placed, Precision sx, Precision sy, Precision sz)
      : fScaled(placed, sx, sy, sz)
  {
    fGlobalConvexity = fScaled.fPlaced->GetUnplacedVolume()->IsConvex();
  }

#if defined(VECGEOM_NVCC)
  /// Constructor based on placed volume
  VECGEOM_CUDA_HEADER_BOTH
  UnplacedScaledShape(VPlacedVolume const *placed, Precision sx, Precision sy, Precision sz, bool globalConvexity)
      : fScaled(placed, sx, sy, sz)
  {
    /* assert(placed->GetTransformation()->IsIdentity());*/
    fGlobalConvexity = globalConvexity;
    /* We must have
         assert(globalConvexity == fPlaced->GetUnplacedVolume()->IsConvex())
       However due to the order we create the geometry on the GPU (i.e. all Unplaced *then* all
       Placed volume, we can not use this information (i.e. 'placed' points to uninitialized memory
       at time this constructor is callled.
    */
  }
#endif

/// Constructor based on unplaced volume
#if !defined(VECGEOM_NVCC)
  UnplacedScaledShape(VUnplacedVolume const *shape, Precision sx, Precision sy, Precision sz)
      : fScaled(nullptr, sx, sy, sz)
  {
    // We need to create a placement with identity transformation from the unplaced version
    // Hopefully we don't need to create a logical volume
    LogicalVolume *lvol = new LogicalVolume("", shape);
    fScaled.fPlaced     = lvol->Place();
    fGlobalConvexity    = fScaled.fPlaced->GetUnplacedVolume()->IsConvex();
  }
#endif

  /// Copy constructor
  //  VECGEOM_CUDA_HEADER_BOTH
  UnplacedScaledShape(UnplacedScaledShape const &other) : fScaled()
  {
    fScaled.fPlaced  = other.fScaled.fPlaced->GetLogicalVolume()->Place();
    fScaled.fScale   = other.fScaled.fScale;
    fGlobalConvexity = other.fGlobalConvexity;
  }

  /// Assignment operator
  //  VECGEOM_CUDA_HEADER_BOTH
  UnplacedScaledShape &operator=(UnplacedScaledShape const &other)
  {
    if (&other != this) {
      fScaled.fPlaced  = other.fScaled.fPlaced->GetLogicalVolume()->Place();
      fScaled.fScale   = other.fScaled.fScale;
      fGlobalConvexity = other.fGlobalConvexity;
    }
    return *this;
  }

  /// Destructor
  VECGEOM_CUDA_HEADER_BOTH
  virtual ~UnplacedScaledShape()
  {
    // The fPlaced was owned by the class, but now it gets deleted before
    // the destructor by GeoManager. The class will be restructured to use
    // the VUnplaceVolume new navigation interfaces after migration of all
    // shapes to VecCore, so this data member will dissapear.

    // delete fScaled.fPlaced;
  }

  /// Getter for the generic scaled shape structure
  VECGEOM_CUDA_HEADER_BOTH
  ScaledShapeStruct<double> const &GetStruct() const { return fScaled; }

  virtual int memory_size() const final { return (sizeof(*this) + fScaled.fPlaced->memory_size()); }

#ifdef VECGEOM_CUDA_INTERFACE
  virtual size_t DeviceSizeOf() const { return DevicePtr<cuda::UnplacedScaledShape>::SizeOf(); }
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu() const;
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const gpu_ptr) const;
#endif

  VECGEOM_FORCE_INLINE
  const VUnplacedVolume *UnscaledShape() const { return fScaled.fPlaced->GetLogicalVolume()->GetUnplacedVolume(); }

  VECGEOM_FORCE_INLINE
  VPlacedVolume const *GetPlaced() const { return fScaled.fPlaced; }

  VECGEOM_FORCE_INLINE
  Scale3D const &GetScale() const { return fScaled.fScale; }

#ifndef VECGEOM_NVCC
  VECGEOM_FORCE_INLINE
  Precision Volume() const
  {
    Precision capacity             = ((VPlacedVolume *)fScaled.fPlaced)->Capacity();
    const Vector3D<Precision> &scl = fScaled.fScale.Scale();
    capacity *= scl[0] * scl[1] * scl[2];
    return capacity;
  }

  VECGEOM_FORCE_INLINE
  Precision Capacity() { return Volume(); }

  VECGEOM_FORCE_INLINE
  Precision SurfaceArea() const
  {
    /// To do - not so easy as for the capacity...
    Precision area = ((VPlacedVolume *)fScaled.fPlaced)->SurfaceArea();
    return area;
  }
#endif // !VECGEOM_NVCC

  VECGEOM_CUDA_HEADER_BOTH
  bool Normal(Vector3D<Precision> const &point, Vector3D<Precision> &normal) const override;

  VECGEOM_CUDA_HEADER_BOTH
  void Extent(Vector3D<Precision> &min, Vector3D<Precision> &max) const override;

  Vector3D<Precision> GetPointOnSurface() const override;

  virtual std::string GetEntityType() const { return "ScaledShape"; }

  VECGEOM_CUDA_HEADER_BOTH
  virtual void Print() const final;

  virtual void Print(std::ostream &os) const final;

  template <TranslationCode trans_code, RotationCode rot_code>
  VECGEOM_CUDA_HEADER_DEVICE
  static VPlacedVolume *Create(LogicalVolume const *const logical_volume, Transformation3D const *const transformation,
#ifdef VECGEOM_NVCC
                               const int id,
#endif // !VECGEOM_NVCC
                               VPlacedVolume *const placement = NULL);

  VECGEOM_CUDA_HEADER_DEVICE
  static VPlacedVolume *CreateSpecializedVolume(LogicalVolume const *const volume,
                                                Transformation3D const *const transformation,
                                                const TranslationCode trans_code, const RotationCode rot_code,
#ifdef VECGEOM_NVCC
                                                const int id,
#endif // !VECGEOM_NVCC
                                                VPlacedVolume *const placement = NULL);

private:
  VECGEOM_CUDA_HEADER_DEVICE
  virtual VPlacedVolume *SpecializedVolume(LogicalVolume const *const volume,
                                           Transformation3D const *const transformation,
                                           const TranslationCode trans_code, const RotationCode rot_code,
#ifdef VECGEOM_NVCC
                                           const int id,
#endif
                                           VPlacedVolume *const placement = NULL) const final
  {
    return CreateSpecializedVolume(volume, transformation, trans_code, rot_code,
#ifdef VECGEOM_NVCC
                                   id,
#endif
                                   placement);
  }
  void SetPlaced(VPlacedVolume const *pvol) { fScaled.fPlaced = pvol; }
  void SetScale(Scale3D const &scale) { fScaled.fScale = scale; }

  friend class GeoManager;
};
}
} // End global namespace

#endif // VECGEOM_VOLUMES_UNPLACESCALEDSHAPE_H_
