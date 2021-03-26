// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.
// For the full list of authors see CONTRIBUTORS.txt and `git log`.

/// @file source/PlacedParallelepiped.cpp
/// @author Johannes de Fine Licht

#include "VecGeom/volumes/PlacedParallelepiped.h"

#include "VecGeom/volumes/Parallelepiped.h"

#ifdef VECGEOM_ROOT
#include "TGeoPara.h"
#endif
#ifdef VECGEOM_GEANT4
#include "G4Para.hh"
#endif

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

#ifndef VECCORE_CUDA

#ifdef VECGEOM_ROOT
TGeoShape const *PlacedParallelepiped::ConvertToRoot() const
{
  return new TGeoPara(GetLabel().c_str(), GetX(), GetY(), GetZ(), GetAlpha() * kRadToDeg, GetTheta() * kRadToDeg,
                      GetPhi() * kRadToDeg);
}
#endif

#ifdef VECGEOM_GEANT4
G4VSolid const *PlacedParallelepiped::ConvertToGeant4() const
{
  return new G4Para(GetLabel(), GetX(), GetY(), GetZ(), GetAlpha(), GetTheta(), GetPhi());
}
#endif

#endif // VECCORE_CUDA

} // namespace VECGEOM_IMPL_NAMESPACE

#ifdef VECCORE_CUDA

VECGEOM_DEVICE_INST_PLACED_VOLUME_ALLSPEC(SpecializedParallelepiped)

#endif // VECCORE_CUDA

} // namespace vecgeom
