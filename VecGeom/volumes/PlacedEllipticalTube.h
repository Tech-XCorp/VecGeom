// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.
// For the full list of authors see CONTRIBUTORS.txt and `git log`.

/// Declaration of the Placed Elliptical Tube volume
/// @file volumes/PlacedEllipticalTube.h
/// @author Raman Sehgal, Evgueni Tcherniaev

#ifndef VECGEOM_VOLUMES_PLACEDELLIPTICALTUBE_H_
#define VECGEOM_VOLUMES_PLACEDELLIPTICALTUBE_H_

#include "VecGeom/base/Global.h"

#include "VecGeom/volumes/PlacedVolume.h"
#include "VecGeom/volumes/UnplacedVolume.h"
#include "VecGeom/volumes/kernel/EllipticalTubeImplementation.h"
#include "VecGeom/volumes/PlacedVolImplHelper.h"
#include "VecGeom/volumes/UnplacedEllipticalTube.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class PlacedEllipticalTube;);
VECGEOM_DEVICE_DECLARE_CONV(class, PlacedEllipticalTube);

inline namespace VECGEOM_IMPL_NAMESPACE {

/// Class for the positioned elliptical tube volume
class PlacedEllipticalTube : public PlacedVolumeImplHelper<UnplacedEllipticalTube, VPlacedVolume> {
  using Base = PlacedVolumeImplHelper<UnplacedEllipticalTube, VPlacedVolume>;

public:
#ifndef VECCORE_CUDA
  using Base::Base;

  /// Constructor
  /// @param label Name of logical volume
  /// @param logicalVolume The logical volume to be positioned
  /// @param transformation The positioning transformation
  PlacedEllipticalTube(char const *const label, LogicalVolume const *const logicalVolume,
                       Transformation3D const *const transformation)
      : Base(label, logicalVolume, transformation)
  {
  }

  /// Constructor
  /// @param logicalVolume The logical volume to be positioned
  /// @param transformation The positioning transformation.
  PlacedEllipticalTube(LogicalVolume const *const logicalVolume, Transformation3D const *const transformation)
      : PlacedEllipticalTube("", logicalVolume, transformation)
  {
  }
#else
  /// CUDA version of constructor
  VECCORE_ATT_DEVICE PlacedEllipticalTube(LogicalVolume const *const logicalVolume,
                                          Transformation3D const *const transformation, const int id, const int copy_no,
                                          const int child_id)
      : Base(logicalVolume, transformation, id, copy_no, child_id)
  {
  }
#endif
  /// Destructor
  VECCORE_ATT_HOST_DEVICE
  virtual ~PlacedEllipticalTube() {}

  /// Getter for x semi-axis
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetDx() const { return GetUnplacedVolume()->GetDx(); }

  /// Getter for y semi-axis
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetDy() const { return GetUnplacedVolume()->GetDy(); }

  /// Getter for the half length in z
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision GetDz() const { return GetUnplacedVolume()->GetDz(); }

  /// Setter for the elliptical tube parameters
  /// @param dx Length of x semi-axis
  /// @param dy Length of y semi-axis
  /// @param dz Half length in z
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  void SetParameters(Precision dx, Precision dy, Precision dz)
  {
    const_cast<UnplacedEllipticalTube *>(GetUnplacedVolume())->SetParameters(dx, dy, dz);
  };

  /// Set x semi-axis
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  void SetDx(Precision dx) { const_cast<UnplacedEllipticalTube *>(GetUnplacedVolume())->SetDx(dx); };

  /// Set y semi-axis
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  void SetDy(Precision dy) { const_cast<UnplacedEllipticalTube *>(GetUnplacedVolume())->SetDy(dy); };

  /// Set half length in z
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  void SetDz(Precision dz) { const_cast<UnplacedEllipticalTube *>(GetUnplacedVolume())->SetDz(dz); };

  VECCORE_ATT_HOST_DEVICE
  virtual void PrintType() const override;
  virtual void PrintType(std::ostream &os) const override;

// Comparison specific
#ifndef VECCORE_CUDA
#ifdef VECGEOM_ROOT
  virtual TGeoShape const *ConvertToRoot() const override;
#endif
#ifdef VECGEOM_GEANT4
  virtual G4VSolid const *ConvertToGeant4() const override;
#endif
#endif // VECCORE_CUDA
};

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif
