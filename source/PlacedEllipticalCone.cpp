// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.
// For the full list of authors see CONTRIBUTORS.txt and `git log`.

/// @file source/PlacedEllipticalCone.cpp
/// @author Raman Sehgal, Evgueni Tcherniaev

#include "VecGeom/volumes/PlacedEllipticalCone.h"
#include "VecGeom/volumes/SpecializedEllipticalCone.h"

#ifdef VECGEOM_ROOT
// Include header file for ROOT Elliptical Cone as done for Tube below
// #include "TGeoEltu.h"
#endif

#ifdef VECGEOM_GEANT4
#include "G4EllipticalCone.hh"
#endif

#include <stdio.h>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

VECCORE_ATT_HOST_DEVICE
void PlacedEllipticalCone::PrintType() const
{
  printf("PlacedEllipticalCone");
}

void PlacedEllipticalCone::PrintType(std::ostream &s) const
{
  s << "PlacedEllipticalCone";
}

#ifndef VECCORE_CUDA

#ifdef VECGEOM_ROOT
TGeoShape const *PlacedEllipticalCone::ConvertToRoot() const
{
  // Return ROOT Elliptical Cone
  return nullptr; // There is no suitable TGeo shape 2019.02.19
}
#endif

#ifdef VECGEOM_GEANT4
G4VSolid const *PlacedEllipticalCone::ConvertToGeant4() const
{
  return new G4EllipticalCone("", GetSemiAxisX(), GetSemiAxisY(), GetZMax(), GetZTopCut());
}
#endif

#endif // VECCORE_CUDA

} // namespace VECGEOM_IMPL_NAMESPACE

#ifdef VECCORE_CUDA

VECGEOM_DEVICE_INST_PLACED_VOLUME_ALLSPEC(SpecializedEllipticalCone)

#endif // VECCORE_CUDA

} // namespace vecgeom
