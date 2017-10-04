#ifndef VECGEOM_VOLUMES_PLACEDTORUS2_H_
#define VECGEOM_VOLUMES_PLACEDTORUS2_H_

#include "base/Global.h"
#include "volumes/PlacedVolume.h"
#include "volumes/UnplacedVolume.h"
#include "volumes/PlacedVolImplHelper.h"
#include "volumes/UnplacedTorus2.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class PlacedTorus2;);
VECGEOM_DEVICE_DECLARE_CONV(class, PlacedTorus2);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedTorus2 : public PlacedVolumeImplHelper<UnplacedTorus2, VPlacedVolume> {

  using Base = PlacedVolumeImplHelper<UnplacedTorus2, VPlacedVolume>;

public:
  // inherit base-class constructors
  using Base::Base;

#ifndef VECCORE_CUDA
  PlacedTorus2(char const *const label, LogicalVolume const *const logical_volume,
               Transformation3D const *const transformation, ::vecgeom::PlacedBox const *const boundingBox)
      : Base(label, logical_volume, transformation, boundingBox)
  {
  }

  PlacedTorus2(LogicalVolume const *const logical_volume, Transformation3D const *const transformation,
               ::vecgeom::PlacedBox const *const boundingBox)
      : PlacedTorus2("", logical_volume, transformation, boundingBox)
  {
  }
#else
  VECCORE_ATT_DEVICE
  PlacedTorus2(LogicalVolume const *const logical_volume, Transformation3D const *const transformation,
               PlacedBox const *const boundingBox, const int id)
      : Base(logical_volume, transformation, boundingBox, id)
  {
  }
#endif
  VECCORE_ATT_HOST_DEVICE
  virtual ~PlacedTorus2() {}

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
#endif
};

////  Maybe this is needed for TorusTypes specializations
// // a placed tube knowing abouts its volume/structural specialization
// template <typename UnplacedTorus_t>
// class SIMDPlacedTorus2 : public PlacedVolumeImplHelper<UnplacedTorus_t, PlacedTorus2> {
//   using Base = PlacedVolumeImplHelper<UnplacedTorus_t, PlacedTorus2>;

// public:
//   typedef UnplacedTorus2 UnplacedShape_t;
//   using Base::Base;
// };

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_PLACEDTORUS2_H_
