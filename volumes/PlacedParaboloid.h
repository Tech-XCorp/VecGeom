//===-- volumes/PlacedParaboloid.h - Instruction class definition -------*- C++ -*-===//
///
/// \file volumes/PlacedParaboloid.h
/// \author Marilena Bandieramonte (marilena.bandieramonte@cern.ch)
/// \brief This file contains the declaration of the PlacedParaboloid class
//===----------------------------------------------------------------------===//
///
/// revision + moving to new backend structure : Raman Sehgal (raman.sehgal@cern.ch)

#ifndef VECGEOM_VOLUMES_PLACEDPARABOLOID_H_
#define VECGEOM_VOLUMES_PLACEDPARABOLOID_H_

#include "base/Global.h"
#include "backend/Backend.h"

#include "volumes/PlacedVolume.h"
#include "volumes/UnplacedVolume.h"
#include "volumes/kernel/ParaboloidImplementation.h"
#include "volumes/PlacedVolImplHelper.h"
#include "volumes/UnplacedParaboloid.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class PlacedParaboloid;);
VECGEOM_DEVICE_DECLARE_CONV(class, PlacedParaboloid);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedParaboloid : public PlacedVolumeImplHelper<UnplacedParaboloid, VPlacedVolume> {
  using Base = PlacedVolumeImplHelper<UnplacedParaboloid, VPlacedVolume>;

public:
#ifndef VECGEOM_NVCC
  // constructor inheritance;
  using Base::Base;
  PlacedParaboloid(char const *const label, LogicalVolume const *const logicalVolume,
                   Transformation3D const *const transformation, vecgeom::PlacedBox const *const boundingBox)
      : Base(label, logicalVolume, transformation, boundingBox)
  {
  }

  PlacedParaboloid(LogicalVolume const *const logicalVolume, Transformation3D const *const transformation,
                   vecgeom::PlacedBox const *const boundingBox)
      : PlacedParaboloid("", logicalVolume, transformation, boundingBox)
  {
  }
#else
  __device__ PlacedParaboloid(LogicalVolume const *const logicalVolume, Transformation3D const *const transformation,
                              PlacedBox const *const boundingBox, const int id)
      : Base(logicalVolume, transformation, boundingBox, id)
  {
  }
#endif
  VECGEOM_CUDA_HEADER_BOTH
  virtual ~PlacedParaboloid() {}

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Precision GetRlo() const { return GetUnplacedVolume()->GetRlo(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Precision GetRhi() const { return GetUnplacedVolume()->GetRhi(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Precision GetDz() const { return GetUnplacedVolume()->GetDz(); }

  VECGEOM_CUDA_HEADER_BOTH
  void SetRlo(Precision arg) { const_cast<UnplacedParaboloid *>(GetUnplacedVolume())->SetRlo(arg); }

  VECGEOM_CUDA_HEADER_BOTH
  void SetRhi(Precision arg) { const_cast<UnplacedParaboloid *>(GetUnplacedVolume())->SetRhi(arg); }

  VECGEOM_CUDA_HEADER_BOTH
  void SetDz(Precision arg) { const_cast<UnplacedParaboloid *>(GetUnplacedVolume())->SetDz(arg); }

  VECGEOM_CUDA_HEADER_BOTH
  virtual void PrintType() const override;
  virtual void PrintType(std::ostream &os) const override;

// Comparison specific
#ifndef VECGEOM_NVCC
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
#endif // VECGEOM_NVCC
};

} // end inline namespace
} // End global namespace

#endif
