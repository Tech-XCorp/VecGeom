/// \file PlacedPolyhedron.h
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#ifndef VECGEOM_VOLUMES_PLACEDPOLYHEDRON_H_
#define VECGEOM_VOLUMES_PLACEDPOLYHEDRON_H_

#include "base/Global.h"
#include "backend/Backend.h"

#include "volumes/PlacedVolume.h"
#include "volumes/UnplacedPolyhedron.h"

namespace VECGEOM_NAMESPACE {

class PlacedPolyhedron : public VPlacedVolume {

public:

#ifndef VECGEOM_NVCC

  PlacedPolyhedron(char const *const label,
            LogicalVolume const *const logicalVolume,
            Transformation3D const *const transformation,
            PlacedBox const *const boundingBox)
      : VPlacedVolume(label, logicalVolume, transformation, boundingBox) {}

  PlacedPolyhedron(LogicalVolume const *const logicalVolume,
            Transformation3D const *const transformation,
            PlacedBox const *const boundingBox)
      : PlacedPolyhedron("", logicalVolume, transformation, boundingBox) {}

#else

  __device__
  PlacedPolyhedron(LogicalVolume const *const logicalVolume,
            Transformation3D const *const transformation,
            PlacedBox const *const boundingBox,
            const int id)
      : VPlacedVolume(logicalVolume, transformation, boundingBox, id) {}

#endif

  virtual ~PlacedPolyhedron() {}

  // Accessors

  VECGEOM_CUDA_HEADER_BOTH
  UnplacedPolyhedron const* GetUnplacedVolume() const {
    return static_cast<UnplacedPolyhedron const *>(
        logical_volume()->unplaced_volume());
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  int GetSideCount() const { return GetUnplacedVolume()->GetSideCount(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  int GetZSegmentCount() const {
    return GetUnplacedVolume()->GetZSegmentCount();
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  bool HasInnerRadii() const { return GetUnplacedVolume()->HasInnerRadii(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  bool HasPhiCutout() const { return GetUnplacedVolume()->HasPhiCutout(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  UnplacedPolyhedron::ZSegment const& GetZSegment(int index) const {
    return GetUnplacedVolume()->GetZSegment(index);
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Array<UnplacedPolyhedron::ZSegment> const& GetZSegments() const {
    return GetUnplacedVolume()->GetZSegments();
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Array<Precision> const& GetZPlanes() const {
    return GetUnplacedVolume()->GetZPlanes();
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Array<Precision> const& GetRMin() const {
    return GetUnplacedVolume()->GetRMin();
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Array<Precision> const& GetRMax() const {
    return GetUnplacedVolume()->GetRMax();
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Vector3D<Precision> GetPhiSection(int i) const {
    return GetUnplacedVolume()->GetPhiSection(i);
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  SOA3D<Precision> const& GetPhiSections() const {
    return GetUnplacedVolume()->GetPhiSections();
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision GetPhiStart() const {
    return GetUnplacedVolume()->GetPhiStart();
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision GetPhiEnd() const {
    return GetUnplacedVolume()->GetPhiEnd();
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision GetPhiDelta() const {
    return GetUnplacedVolume()->GetPhiDelta();
  }

  VECGEOM_CUDA_HEADER_BOTH
  int PhiSegmentIndex(Vector3D<Precision> const &point) const;

  // CUDA specific

  virtual int memory_size() const { return sizeof(*this); }

#ifdef VECGEOM_CUDA_INTERFACE
  virtual VPlacedVolume* CopyToGpu(
      LogicalVolume const *const logical_volume,
      Transformation3D const *const transformation,
      VPlacedVolume *const gpu_ptr) const;
  virtual VPlacedVolume* CopyToGpu(
      LogicalVolume const *const logical_volume,
      Transformation3D const *const transformation) const;
#endif

  // Comparison specific

#ifndef VECGEOM_NVCC
  virtual VPlacedVolume const* ConvertToUnspecialized() const;
#ifdef VECGEOM_ROOT
  virtual TGeoShape const* ConvertToRoot() const;
#endif
#ifdef VECGEOM_USOLIDS
  virtual ::VUSolid const* ConvertToUSolids() const;
#endif
#ifdef VECGEOM_GEANT4
  virtual G4VSolid const* ConvertToGeant4() const;
#endif
#endif // VECGEOM_NVCC

};

} // End global namespace

#endif // VECGEOM_VOLUMES_PLACEDPOLYHEDRON_H_
