/*
 * @file   volumes/PlacedTrapezoid.h
 * @author Guilherme Lima (lima 'at' fnal 'dot' gov)
 *
 * 2014-05-01 - Created, based on the Parallelepiped draft
 * 2016-07-25 G.Lima     - Migrated to new helpers + VecCore
 */

#ifndef VECGEOM_VOLUMES_PLACEDTRAPEZOID_H_
#define VECGEOM_VOLUMES_PLACEDTRAPEZOID_H_

#include "base/Global.h"
#include "volumes/PlacedVolume.h"
#include "volumes/UnplacedVolume.h"
#include "volumes/kernel/TrapezoidImplementation.h"
#include "volumes/PlacedVolImplHelper.h"
#include "volumes/UnplacedTrapezoid.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class TGeoTrap;);
VECGEOM_DEVICE_FORWARD_DECLARE(class PlacedTrapezoid;);
VECGEOM_DEVICE_DECLARE_CONV(class, PlacedTrapezoid);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedTrapezoid : public PlacedVolumeImplHelper<UnplacedTrapezoid, VPlacedVolume> {

  using Base = PlacedVolumeImplHelper<UnplacedTrapezoid, VPlacedVolume>;

public:
#ifndef VECCORE_CUDA
  // constructor inheritance;
  using Base::Base;

  PlacedTrapezoid(char const *const label, LogicalVolume const *const logicalVolume,
                  Transformation3D const *const transformation, vecgeom::PlacedBox const *const boundingBox)
      : Base(label, logicalVolume, transformation, boundingBox)
  {
  }

  PlacedTrapezoid(LogicalVolume const *const logicalVolume, Transformation3D const *const transformation,
                  vecgeom::PlacedBox const *const boundingBox)
      : PlacedTrapezoid("", logicalVolume, transformation, boundingBox)
  {
  }

#else
  VECCORE_ATT_DEVICE PlacedTrapezoid(LogicalVolume const *const logicalVolume,
                                     Transformation3D const *const transformation, PlacedBox const *const boundingBox,
                                     const int id)
      : Base(logicalVolume, transformation, boundingBox, id)
  {
  }
#endif

  VECCORE_ATT_HOST_DEVICE
  virtual ~PlacedTrapezoid() {}

  VECCORE_ATT_HOST_DEVICE
  virtual void PrintType() const override;
  virtual void PrintType(std::ostream &os) const override;

// Comparison specific
#ifndef VECCORE_CUDA
  virtual VPlacedVolume const *ConvertToUnspecialized() const override;

#ifdef VECGEOM_ROOT
  virtual TGeoShape const *ConvertToRoot() const override;
#endif
#ifdef VECGEOM_GEANT4
  virtual G4VSolid const *ConvertToGeant4() const override;
#endif

#endif // VECCORE_CUDA
};

} // end inline namespace
} // End global namespace

#endif
