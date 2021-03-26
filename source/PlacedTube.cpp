/// \file PlacedTube.cpp
/// \author Georgios Bitzes (georgios.bitzes@cern.ch)

#include "VecGeom/volumes/PlacedTube.h"
#include "VecGeom/volumes/Tube.h"
#include "VecGeom/base/Vector3D.h"

#ifdef VECGEOM_ROOT
#include "TGeoTube.h"
#endif

#ifdef VECGEOM_GEANT4
#include "G4Tubs.hh"
#endif

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

VECCORE_ATT_HOST_DEVICE
void PlacedTube::PrintType() const
{
  printf("PlacedTube");
}

void PlacedTube::PrintType(std::ostream &s) const
{
  s << "PlacedTube";
}

#ifndef VECCORE_CUDA

#ifdef VECGEOM_ROOT
TGeoShape const *PlacedTube::ConvertToRoot() const
{
  UnplacedTube const *t = static_cast<UnplacedTube const *>(GetUnplacedVolume());
  if (t->dphi() >= 2 * M_PI) return new TGeoTube(GetLabel().c_str(), t->rmin(), t->rmax(), t->z());
  return new TGeoTubeSeg(GetLabel().c_str(), t->rmin(), t->rmax(), t->z(), t->sphi() * (180 / M_PI),
                         (t->sphi() + t->dphi()) * (180 / M_PI));
}
#endif

#ifdef VECGEOM_GEANT4
G4VSolid const *PlacedTube::ConvertToGeant4() const
{
  UnplacedTube const *t = static_cast<UnplacedTube const *>(GetUnplacedVolume());
  return new G4Tubs(GetLabel().c_str(), t->rmin(), t->rmax(), t->z(), t->sphi(), t->dphi());
}
#endif

#endif // VECCORE_CUDA

} // End impl namespace

#ifdef VECCORE_CUDA

VECGEOM_DEVICE_INST_PLACED_VOLUME(PlacedTube)

#endif

} // End global namespace
