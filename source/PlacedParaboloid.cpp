/// \file PlacedParaboloid.cpp
/// \author Marilena Bandieramonte (marilena.bandieramonte@cern.ch)

#include "volumes/PlacedParaboloid.h"
#include "volumes/SpecializedParaboloid.h"
#ifdef VECGEOM_ROOT
#include "TGeoParaboloid.h"
#endif
#if defined(VECGEOM_USOLIDS) && !defined(VECGEOM_REPLACE_USOLIDS)
#include "UBox.hh"
#endif
#ifdef VECGEOM_GEANT4
#include "G4Paraboloid.hh"
#endif

#include <stdio.h>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

VECGEOM_CUDA_HEADER_BOTH
void PlacedParaboloid::PrintType() const
{
  printf("PlacedParaboloid");
}

void PlacedParaboloid::PrintType(std::ostream &s) const
{
  s << "PlacedParaboloid";
}

#ifndef VECGEOM_NVCC

VPlacedVolume const *PlacedParaboloid::ConvertToUnspecialized() const
{
  return new SimpleParaboloid(GetLabel().c_str(), logical_volume_, GetTransformation());
}

#ifdef VECGEOM_ROOT
TGeoShape const *PlacedParaboloid::ConvertToRoot() const
{
  return new TGeoParaboloid(GetLabel().c_str(), GetRlo(), GetRhi(), GetDz());
}
#endif

#if defined(VECGEOM_USOLIDS) && !defined(VECGEOM_REPLACE_USOLIDS)
::VUSolid const *PlacedParaboloid::ConvertToUSolids() const
{
  std::cerr << "**************************************************************\n";
  std::cerr << "WARNING: Paraboloid unsupported for USolids.; returning NULL\n";
  std::cerr << "**************************************************************\n";
  return nullptr;
}
#endif

#ifdef VECGEOM_GEANT4
G4VSolid const *PlacedParaboloid::ConvertToGeant4() const
{
  return new G4Paraboloid(GetLabel(), GetDz(), GetRlo(), GetRhi());
}
#endif

#endif // VECGEOM_NVCC

} // End impl namespace

#ifdef VECGEOM_NVCC

VECGEOM_DEVICE_INST_PLACED_VOLUME_ALLSPEC(SpecializedParaboloid)

#endif // VECGEOM_NVCC

} // End namespace vecgeom
