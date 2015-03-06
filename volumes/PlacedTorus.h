/// \file PlacedTorus.h

#ifndef VECGEOM_VOLUMES_PLACEDTORUS_H_
#define VECGEOM_VOLUMES_PLACEDTORUS_H_

#include "base/Global.h"
#include "backend/Backend.h"
 
#include "volumes/PlacedVolume.h"
#include "volumes/UnplacedVolume.h"
#include "volumes/kernel/TorusImplementation.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE( class PlacedTorus; )
VECGEOM_DEVICE_DECLARE_CONV( PlacedTorus );

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedTorus : public VPlacedVolume {

public:

  typedef UnplacedTorus UnplacedShape_t;


#ifndef VECGEOM_NVCC

  PlacedTorus(char const *const label,
	      LogicalVolume const *const logical_volume,
	      Transformation3D const *const transformation,
	      PlacedBox const *const boundingBox)
      : VPlacedVolume(label, logical_volume, transformation, boundingBox) {}

  PlacedTorus(LogicalVolume const *const logical_volume,
	      Transformation3D const *const transformation,
	      PlacedBox const *const boundingBox)
      : PlacedTorus("", logical_volume, transformation, boundingBox) {}

#else

  __device__
  PlacedTorus(LogicalVolume const *const logical_volume,
	      Transformation3D const *const transformation,
	      PlacedBox const *const boundingBox, const int id)
      : VPlacedVolume(logical_volume, transformation, boundingBox, id) {}

#endif

  VECGEOM_CUDA_HEADER_BOTH
  virtual ~PlacedTorus() {}

  VECGEOM_CUDA_HEADER_BOTH
  UnplacedTorus const* GetUnplacedVolume() const {
    return static_cast<UnplacedTorus const *>(
        GetLogicalVolume()->unplaced_volume());
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision rmin() const { return GetUnplacedVolume()->rmin(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision rmax() const { return GetUnplacedVolume()->rmax(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision rtor() const { return GetUnplacedVolume()->rtor(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision sphi() const { return GetUnplacedVolume()->sphi(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision dphi() const { return GetUnplacedVolume()->dphi(); }


  VECGEOM_CUDA_HEADER_BOTH
  virtual Precision Capacity() { return GetUnplacedVolume()->volume(); }

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

} } // End global namespace

#endif // VECGEOM_VOLUMES_PLACEDTORUS_H_










