/*
 * PlacedPolycone.cpp
 *
 *  Created on: Dec 9, 2014
 *      Author: swenzel
 */

#include "volumes/SpecializedPolycone.h"
#include <iostream>

#ifdef VECGEOM_ROOT
#include "TGeoPcon.h"
#endif

#ifdef VECGEOM_GEANT4
#include "G4Polycone.hh"
#endif

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

#ifndef VECCORE_CUDA
VPlacedVolume const *PlacedPolycone::ConvertToUnspecialized() const
{
  return new SimplePolycone(GetLabel().c_str(), GetLogicalVolume(), GetTransformation());
}

#ifdef VECGEOM_ROOT
TGeoShape const *PlacedPolycone::ConvertToRoot() const
{
  UnplacedPolycone const *unplaced = GetUnplacedVolume();

  std::vector<double> rmin;
  std::vector<double> rmax;
  std::vector<double> z;
  unplaced->ReconstructSectionArrays(z, rmin, rmax);

  TGeoPcon *rootshape =
      new TGeoPcon(unplaced->GetStruct().fStartPhi * kRadToDeg, unplaced->GetStruct().fDeltaPhi * kRadToDeg, z.size());

  if (unplaced->GetStruct().fNz != z.size()) std::cout << "WARNING: Inconsistency in number of polycone sections\n";

  for (unsigned int i = 0; i < unplaced->GetStruct().fNz; ++i)
    rootshape->DefineSection(i, z[i], rmin[i], rmax[i]);

  return rootshape;
}
#endif

#ifdef VECGEOM_GEANT4
G4VSolid const *PlacedPolycone::ConvertToGeant4() const
{
  UnplacedPolycone const *unplaced = GetUnplacedVolume();

  std::vector<double> rmin;
  std::vector<double> rmax;
  std::vector<double> z;
  unplaced->ReconstructSectionArrays(z, rmin, rmax);

  G4Polycone *g4shape = new G4Polycone("", unplaced->GetStruct().fStartPhi, unplaced->GetStruct().fDeltaPhi, z.size(),
                                       &z[0], &rmin[0], &rmax[0]);

  return g4shape;
}
#endif
#endif // ! VECCORE_CUDA
}

#ifdef VECCORE_CUDA
VECGEOM_DEVICE_INST_PLACED_VOLUME_ALLSPEC(SpecializedPolycone)
#endif

} // end global namespace
