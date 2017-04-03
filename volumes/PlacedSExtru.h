#ifndef VECGEOM_VOLUMES_PLACEDSEXTRU_H_
#define VECGEOM_VOLUMES_PLACEDSEXTRU_H_

#include "base/Global.h"
#include "volumes/PlacedVolume.h"
#include "volumes/PlacedVolImplHelper.h"
#include "volumes/UnplacedSExtruVolume.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class PlacedSExtru;);
VECGEOM_DEVICE_DECLARE_CONV(class, PlacedSExtru);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedSExtru : public PlacedVolumeImplHelper<UnplacedSExtruVolume, VPlacedVolume> {
  using Base = PlacedVolumeImplHelper<UnplacedSExtruVolume, VPlacedVolume>;

public:
#ifndef VECCORE_CUDA
  // constructor inheritance;
  using Base::Base;
  PlacedSExtru(char const *const label, LogicalVolume const *const logicalVolume,
               Transformation3D const *const transformation, vecgeom::PlacedBox const *const boundingBox)
      : Base(label, logicalVolume, transformation, boundingBox)
  {
  }

  PlacedSExtru(LogicalVolume const *const logicalVolume, Transformation3D const *const transformation,
               vecgeom::PlacedBox const *const boundingBox)
      : PlacedSExtru("", logicalVolume, transformation, boundingBox)
  {
  }
#else
  __device__ PlacedSExtru(LogicalVolume const *const logicalVolume, Transformation3D const *const transformation,
                          PlacedBox const *const boundingBox, const int id)
      : Base(logicalVolume, transformation, boundingBox, id)
  {
  }
#endif
  VECCORE_ATT_HOST_DEVICE
  virtual ~PlacedSExtru() {}

  VECCORE_ATT_HOST_DEVICE
  virtual void PrintType() const override;
  virtual void PrintType(std::ostream &os) const override;

  virtual Vector3D<Precision> GetPointOnSurface() const override;

// Comparison specific
#ifndef VECCORE_CUDA
  virtual VPlacedVolume const *ConvertToUnspecialized() const override;
#ifdef VECGEOM_ROOT
  virtual TGeoShape const *ConvertToRoot() const override;
#endif
#if defined(VECGEOM_USOLIDS) && !defined(VECGEOM_REPLACE_USOLIDS)
  virtual ::VUSolid const *ConvertToUSolids() const override;
#endif
#ifdef VECGEOM_GEANT4
  virtual G4VSolid const *ConvertToGeant4() const override;
#endif
#endif // VECCORE_CUDA
};

} // end inline namespace
} // End global namespace

#endif
