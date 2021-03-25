/*
 * PlacedPolycone.h
 *
 *  Created on: May 14, 2014
 *      Author: swenzel
 */
#ifndef VECGEOM_VOLUMES_PLACEDPOLYCONE_H_
#define VECGEOM_VOLUMES_PLACEDPOLYCONE_H_

#include "VecGeom/base/Global.h"
#ifndef VECCORE_CUDA
#include "VecGeom/base/RNG.h"
#include <cmath>
#endif
#include "VecGeom/volumes/PlacedVolume.h"
#include "VecGeom/volumes/UnplacedVolume.h"
#include "VecGeom/volumes/kernel/PolyconeImplementation.h"
#include "VecGeom/volumes/PlacedVolImplHelper.h"
#include "VecGeom/volumes/UnplacedPolycone.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class PlacedPolycone;);
VECGEOM_DEVICE_DECLARE_CONV(class, PlacedPolycone);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedPolycone : public PlacedVolumeImplHelper<UnplacedPolycone> {
  using Base = PlacedVolumeImplHelper<UnplacedPolycone>;

public:
  using Base::Base;
#ifndef VECCORE_CUDA
  PlacedPolycone(char const *const label, LogicalVolume const *const logicalVolume,
                 Transformation3D const *const transformation)
      : Base(label, logicalVolume, transformation)
  {
  }

  PlacedPolycone(LogicalVolume const *const logicalVolume, Transformation3D const *const transformation)
      : PlacedPolycone("", logicalVolume, transformation)
  {
  }
#else
  VECCORE_ATT_DEVICE PlacedPolycone(LogicalVolume const *const logicalVolume,
                                    Transformation3D const *const transformation, const int id, const int copy_no,
                                    const int child_id)
      : Base(logicalVolume, transformation, id, copy_no, child_id)
  {
  }
#endif

  VECCORE_ATT_HOST_DEVICE
  virtual ~PlacedPolycone() {}

  VECCORE_ATT_HOST_DEVICE
  virtual void PrintType() const override;
  virtual void PrintType(std::ostream &os) const override;

  VECCORE_ATT_HOST_DEVICE
  UnplacedPolycone const *GetUnplacedVolume() const
  {
    return static_cast<UnplacedPolycone const *>(GetLogicalVolume()->GetUnplacedVolume());
  }

  bool IsOpen() const { return (GetUnplacedVolume()->GetDeltaPhi() < kTwoPi); }
  Precision GetStartPhi() const { return GetUnplacedVolume()->GetStartPhi(); }
  Precision GetEndPhi() const { return GetUnplacedVolume()->GetEndPhi(); }
  int GetNumRZCorner() const { return 2 * (int)(GetUnplacedVolume()->GetNz()); } // in USolids nCorners = 2*nPlanes

  PolyconeHistorical *GetOriginalParameters() const { return GetUnplacedVolume()->GetOriginalParameters(); }

  void Reset() { const_cast<UnplacedPolycone *>(GetUnplacedVolume())->Reset(); }

#ifndef VECCORE_CUDA
#ifdef VECGEOM_ROOT
  virtual TGeoShape const *ConvertToRoot() const override;
#endif
#ifdef VECGEOM_GEANT4
  virtual G4VSolid const *ConvertToGeant4() const override;
#endif
#endif // VECGEOM_BENCHMARK

#if !defined(VECCORE_CUDA)
  virtual Precision SurfaceArea() const override { return GetUnplacedVolume()->SurfaceArea(); }
#endif

}; // end class

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_PLACEDPOLYCONE_H_
