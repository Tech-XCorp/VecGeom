/// @file PlacedSphere.cpp
/// @author Raman Sehgal (raman.sehgal@cern.ch)

#include "VecGeom/volumes/PlacedSphere.h"
#include "VecGeom/volumes/Sphere.h"
#include "VecGeom/volumes/SpecializedSphere.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

VECCORE_ATT_HOST_DEVICE
void PlacedSphere::PrintType() const
{
  printf("PlacedSphere");
}

void PlacedSphere::PrintType(std::ostream &s) const
{
  s << "PlacedSphere";
}

} // End impl namespace

#ifdef VECCORE_CUDA

VECGEOM_DEVICE_INST_PLACED_VOLUME_ALLSPEC(SpecializedSphere)

#endif // VECCORE_CUDA

} // End global namespace
