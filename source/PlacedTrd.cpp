// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.
// For the full list of authors see CONTRIBUTORS.txt and `git log`.

/// @file source/PlacedTrd.cpp
/// @author Georgios Bitzes

#include "VecGeom/volumes/Trd.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

VECCORE_ATT_HOST_DEVICE
void PlacedTrd::PrintType() const
{
  printf("PlacedTrd");
}

void PlacedTrd::PrintType(std::ostream &s) const
{
  s << "PlacedTrd";
}

#ifndef VECCORE_CUDA

#ifdef VECGEOM_ROOT
TGeoShape const *PlacedTrd::ConvertToRoot() const
{
  return GetUnplacedVolume()->ConvertToRoot(GetName());
}
#endif

#ifdef VECGEOM_GEANT4
G4VSolid const *PlacedTrd::ConvertToGeant4() const
{
  return GetUnplacedVolume()->ConvertToGeant4(GetName());
}
#endif

#endif // VECCORE_CUDA

} // namespace VECGEOM_IMPL_NAMESPACE

#ifdef VECCORE_CUDA

VECGEOM_DEVICE_INST_PLACED_VOLUME_ALLSPEC_3(SpecializedTrd, TrdTypes::UniversalTrd)

#endif

} // namespace vecgeom
