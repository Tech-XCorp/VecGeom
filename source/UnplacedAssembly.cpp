// LICENSING INFORMATION TBD

#include "volumes/UnplacedAssembly.h"
#include "volumes/PlacedAssembly.h"
#include "navigation/SimpleLevelLocator.h"
#include "management/ABBoxManager.h" // for Extent == bounding box calculation
#include "base/RNG.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

UnplacedAssembly::UnplacedAssembly() : fLogicalVolume(nullptr), fLowerCorner(-kInfLength), fUpperCorner(kInfLength)
{
  fIsAssembly = true;
}

UnplacedAssembly::~UnplacedAssembly()
{
}

void UnplacedAssembly::AddVolume(VPlacedVolume const *v)
{
  fLogicalVolume->PlaceDaughter(v);
}

VECCORE_ATT_HOST_DEVICE
void UnplacedAssembly::Print() const
{
  printf("UnplacedAssembly ");
}

void UnplacedAssembly::Print(std::ostream &os) const
{
  os << "UnplacedAssembly ";
}

//______________________________________________________________________________
void UnplacedAssembly::Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const
{
#ifndef VECGEOM_NVCC
  auto &abboxmgr = ABBoxManager::Instance();

  // Returns the full 3D cartesian extent of the solid.
  // Loop nodes and get their extent
  aMin.Set(kInfLength);
  aMax.Set(-kInfLength);
  for (VPlacedVolume const *pv : fLogicalVolume->GetDaughters()) {
    Vector3D<Precision> lower, upper;
    abboxmgr.ComputeABBox(pv, &lower, &upper);
    aMin.Set(std::min(lower.x(), aMin.x()), std::min(lower.y(), aMin.y()), std::min(lower.z(), aMin.z()));
    aMax.Set(std::max(upper.x(), aMax.x()), std::max(upper.y(), aMax.y()), std::max(upper.z(), aMax.z()));
  }
#endif
}

Vector3D<Precision> UnplacedAssembly::GetPointOnSurface() const
{
  // pick one of the constituents for now
  // should improve this to sample according to surface area
  const auto ndaughters = fLogicalVolume->GetDaughters().size();
  const size_t selected = RNG::Instance().uniform() * ndaughters;

  const auto selectedplaced = fLogicalVolume->GetDaughters()[selected];
  Vector3D<Precision> sp    = selectedplaced->GetPointOnSurface();
  // this is in the reference frame of the selected daughter
  // we need to return it in the reference of this assembly

  return selectedplaced->GetTransformation()->InverseTransform(sp);
}

Precision UnplacedAssembly::Capacity() const
{
  Precision capacity = 0.;
  // loop over nodes and sum the capacity of all their unplaced volumes
  for (VPlacedVolume const *pv : fLogicalVolume->GetDaughters()) {
    capacity += const_cast<VPlacedVolume *>(pv)->Capacity();
  }
  return capacity;
}

Precision UnplacedAssembly::SurfaceArea() const
{
  Precision area = 0.;
  // loop over nodes and sum the area of all their unplaced volumes
  // (this might be incorrect in case 2 constituents are touching)
  for (VPlacedVolume const *pv : fLogicalVolume->GetDaughters()) {
    area += const_cast<VPlacedVolume *>(pv)->SurfaceArea();
  }
  return area;
}

#ifndef VECGEOM_NVCC
VPlacedVolume *UnplacedAssembly::SpecializedVolume(LogicalVolume const *const volume,
                                                   Transformation3D const *const transformation,
                                                   const TranslationCode trans_code, const RotationCode rot_code,
                                                   VPlacedVolume *const placement) const
{
  if (placement) {
    return new (placement) PlacedAssembly("", volume, transformation);
  }
  return new PlacedAssembly("", volume, transformation);
}
#else
__device__ VPlacedVolume *UnplacedAssembly::SpecializedVolume(LogicalVolume const *const volume,
                                                              Transformation3D const *const transformation,
                                                              const TranslationCode trans_code,
                                                              const RotationCode rot_code, const int id,
                                                              VPlacedVolume *const placement) const
{
  if (placement) {
    return new (placement) PlacedAssembly("", volume, transformation, nullptr, id);
  }
  return new PlacedAssembly("", volume, transformation, nullptr, id);
}
#endif

#ifdef VECGEOM_CUDA_INTERFACE

DevicePtr<cuda::VUnplacedVolume> UnplacedAssembly::CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const in_gpu_ptr) const
{
  return CopyToGpuImpl<UnplacedAssembly>(in_gpu_ptr);
}

DevicePtr<cuda::VUnplacedVolume> UnplacedAssembly::CopyToGpu() const
{
  return CopyToGpuImpl<UnplacedAssembly>();
}

#endif // VECGEOM_CUDA_INTERFACE

} // End impl namespace

#ifdef VECGEOM_NVCC

namespace cxx {

template size_t DevicePtr<cuda::UnplacedAssembly>::SizeOf();
template void DevicePtr<cuda::UnplacedAssembly>::Construct() const;

} // End cxx namespace

#endif

} // End vecgeom namespace
