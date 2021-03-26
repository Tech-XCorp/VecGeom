/*
 * PlacedPolycone.cpp
 *
 *  Created on: Dec 9, 2014
 *      Author: swenzel
 */

#include "VecGeom/volumes/SpecializedPolycone.h"
#include <iostream>

#ifdef VECGEOM_ROOT
#include "TGeoPcon.h"
#endif

#ifdef VECGEOM_GEANT4
#include "G4Polycone.hh"
#endif

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

VECCORE_ATT_HOST_DEVICE
void PlacedPolycone::PrintType() const
{
  printf("PlacedPolycone");
}

void PlacedPolycone::PrintType(std::ostream &s) const
{
  s << "PlacedPolycone";
}

#ifndef VECCORE_CUDA
#ifdef VECGEOM_ROOT
TGeoShape const *PlacedPolycone::ConvertToRoot() const
{
  return GetUnplacedVolume()->ConvertToRoot(GetName());
}
#endif

#ifdef VECGEOM_GEANT4
G4VSolid const *PlacedPolycone::ConvertToGeant4() const
{
  return GetUnplacedVolume()->ConvertToGeant4(GetName());
}
#endif
#endif // ! VECCORE_CUDA
} // namespace VECGEOM_IMPL_NAMESPACE

#ifdef VECCORE_CUDA
VECGEOM_DEVICE_INST_PLACED_VOLUME_ALLSPEC_3(SpecializedPolycone, ConeTypes::UniversalCone)
#endif // VECCORE_CUDA

} // namespace vecgeom
