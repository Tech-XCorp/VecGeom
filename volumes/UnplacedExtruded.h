/// @file UnplacedExtruded.h
/// @author Mihaela Gheata (mihaela.gheata@cern.ch)

#ifndef VECGEOM_VOLUMES_UNPLACEDEXTRUDED_H_
#define VECGEOM_VOLUMES_UNPLACEDEXTRUDED_H_

#include "base/Global.h"
#include "base/AlignedBase.h"
#include "volumes/UnplacedVolume.h"
#include "ExtrudedStruct.h"
#include "volumes/kernel/ExtrudedImplementation.h"
#include "volumes/UnplacedVolumeImplHelper.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class UnplacedExtruded;);
VECGEOM_DEVICE_DECLARE_CONV(class, UnplacedExtruded);

inline namespace VECGEOM_IMPL_NAMESPACE {

class UnplacedExtruded : public LoopUnplacedVolumeImplHelper<ExtrudedImplementation>, public AlignedBase {

  // template <typename U>
  // using vector_t = vecgeom::Vector<U>;
  template <typename U>
  using vector_t = std::vector<U>;

private:
  ExtrudedStruct fXtru; ///< Structure storing the data for the tessellated solid

public:
  /** @brief Dummy constructor */
  VECCORE_ATT_HOST_DEVICE
  UnplacedExtruded() : fXtru() {}

  /** @brief Constructor providing polygone vertices and sections */
  VECCORE_ATT_HOST_DEVICE
  UnplacedExtruded(int nvertices, XtruVertex2 const *vertices, int nsections, XtruSection const *sections)
      : fXtru(nvertices, vertices, nsections, sections)
  {
    fGlobalConvexity = (nsections == 2) && fXtru.IsConvexPolygon();
  }

  VECCORE_ATT_HOST_DEVICE
  ExtrudedStruct const &GetStruct() const { return fXtru; }

  /** @brief Initialize */
  VECCORE_ATT_HOST_DEVICE
  void Initialize(int nvertices, XtruVertex2 const *vertices, int nsections, XtruSection const *sections)
  {
    fXtru.Initialize(nvertices, vertices, nsections, sections);
    fGlobalConvexity = (nsections == 2) && fXtru.IsConvexPolygon();
  }

  /** @brief GetThe number of sections */
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  size_t GetNSections() const { return fXtru.GetNSections(); }

  /** @brief Get section i */
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  XtruSection GetSection(int i) const { return fXtru.GetSection(i); }

  /** @brief Get the number of vertices */
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  size_t GetNVertices() const { return fXtru.GetNVertices(); }

  /** @brief Get the polygone vertex i */
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  void GetVertex(int i, double &x, double &y) const { fXtru.GetVertex(i, x, y); }

  VECCORE_ATT_HOST_DEVICE
  void Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const override;

  // Computes capacity of the shape in [length^3]
  VECCORE_ATT_HOST_DEVICE
  Precision Capacity() const override;

  // VECCORE_ATT_HOST_DEVICE
  Precision SurfaceArea() const override;

  VECCORE_ATT_HOST_DEVICE
  int ChooseSurface() const;

  Vector3D<Precision> SamplePointOnSurface() const override;

  VECCORE_ATT_HOST_DEVICE
  bool Normal(Vector3D<Precision> const &point, Vector3D<Precision> &normal) const override;

  VECCORE_ATT_HOST_DEVICE
  virtual void Print() const final;

  virtual void Print(std::ostream &os) const final;

  virtual int memory_size() const final { return sizeof(*this); }
  std::string GetEntityType() const { return "Extruded"; }

  template <TranslationCode transCodeT, RotationCode rotCodeT>
  VECCORE_ATT_DEVICE
  static VPlacedVolume *Create(LogicalVolume const *const logical_volume, Transformation3D const *const transformation,
#ifdef VECCORE_CUDA
                               const int id,
#endif
                               VPlacedVolume *const placement = NULL);
#ifdef VECGEOM_CUDA_INTERFACE
  virtual size_t DeviceSizeOf() const override { return DevicePtr<cuda::UnplacedExtruded>::SizeOf(); }
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu() const override;
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const gpu_ptr) const override;
#endif

  std::ostream &StreamInfo(std::ostream &os) const;

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
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_UNPLACEDEXTRUDED_H_
