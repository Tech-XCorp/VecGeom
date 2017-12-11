/*
 * UnplacedBooleanVolume.cpp
 *
 *  Created on: 07.11.2014
 *      Author: swenzel
 */

#include "base/Global.h"
#include "volumes/UnplacedBooleanVolume.h"
#include "volumes/SpecializedBooleanVolume.h"
#include "management/VolumeFactory.h"
#include "volumes/utilities/GenerationUtilities.h"
#include "volumes/LogicalVolume.h"
#include "volumes/PlacedVolume.h"

#ifdef VECGEOM_CUDA_INTERFACE
#include "management/CudaManager.h"
#endif

namespace vecgeom {

inline namespace VECGEOM_IMPL_NAMESPACE {

template <>
template <TranslationCode transCodeT, RotationCode rotCodeT>
VECCORE_ATT_DEVICE
VPlacedVolume *UnplacedBooleanVolume<kSubtraction>::Create(LogicalVolume const *const logical_volume,
                                                           Transformation3D const *const transformation,
#ifdef VECCORE_CUDA
                                                           const int id,
#endif
                                                           VPlacedVolume *const placement)
{
  return CreateSpecializedWithPlacement<SpecializedBooleanVolume<kSubtraction, transCodeT, rotCodeT>>(
      logical_volume, transformation,
#ifdef VECCORE_CUDA
      id,
#endif
      placement); // TODO: add bounding box?
}

template <>
template <TranslationCode transCodeT, RotationCode rotCodeT>
VECCORE_ATT_DEVICE
VPlacedVolume *UnplacedBooleanVolume<kUnion>::Create(LogicalVolume const *const logical_volume,
                                                     Transformation3D const *const transformation,
#ifdef VECCORE_CUDA
                                                     const int id,
#endif
                                                     VPlacedVolume *const placement)
{
  return CreateSpecializedWithPlacement<SpecializedBooleanVolume<kUnion, transCodeT, rotCodeT>>(
      logical_volume, transformation,
#ifdef VECCORE_CUDA
      id,
#endif
      placement); // TODO: add bounding box?
}

template <>
template <TranslationCode transCodeT, RotationCode rotCodeT>
VECCORE_ATT_DEVICE
VPlacedVolume *UnplacedBooleanVolume<kIntersection>::Create(LogicalVolume const *const logical_volume,
                                                            Transformation3D const *const transformation,
#ifdef VECCORE_CUDA
                                                            const int id,
#endif
                                                            VPlacedVolume *const placement)
{
  return CreateSpecializedWithPlacement<SpecializedBooleanVolume<kIntersection, transCodeT, rotCodeT>>(
      logical_volume, transformation,
#ifdef VECCORE_CUDA
      id,
#endif
      placement); // TODO: add bounding box?
}

template <>
VECCORE_ATT_DEVICE
VPlacedVolume *UnplacedBooleanVolume<kSubtraction>::SpecializedVolume(LogicalVolume const *const volume,
                                                                      Transformation3D const *const transformation,
                                                                      const TranslationCode trans_code,
                                                                      const RotationCode rot_code,
#ifdef VECCORE_CUDA
                                                                      const int id,
#endif
                                                                      VPlacedVolume *const placement) const
{
#ifndef VECCORE_CUDA

  return VolumeFactory::CreateByTransformation<UnplacedBooleanVolume<kSubtraction>>(volume, transformation, trans_code,
                                                                                    rot_code,
#ifdef VECCORE_CUDA
                                                                                    id,
#endif
                                                                                    placement);

#else
  // Compiling the above code with nvcc 6.5 faile with the error:
  // nvcc error   : 'ptxas' died due to signal 11 (Invalid memory reference)
  // at least when optimized.
  return nullptr;
#endif
}

template <>
VECCORE_ATT_DEVICE
VPlacedVolume *UnplacedBooleanVolume<kUnion>::SpecializedVolume(LogicalVolume const *const volume,
                                                                Transformation3D const *const transformation,
                                                                const TranslationCode trans_code,
                                                                const RotationCode rot_code,
#ifdef VECCORE_CUDA
                                                                const int id,
#endif
                                                                VPlacedVolume *const placement) const
{
#ifndef VECCORE_CUDA

  return VolumeFactory::CreateByTransformation<UnplacedBooleanVolume<kUnion>>(volume, transformation, trans_code,
                                                                              rot_code,
#ifdef VECCORE_CUDA
                                                                              id,
#endif
                                                                              placement);

#else
  // Compiling the above code with nvcc 6.5 faile with the error:
  // nvcc error   : 'ptxas' died due to signal 11 (Invalid memory reference)
  // at least when optimized.
  return nullptr;
#endif
}

template <>
VECCORE_ATT_DEVICE
VPlacedVolume *UnplacedBooleanVolume<kIntersection>::SpecializedVolume(LogicalVolume const *const volume,
                                                                       Transformation3D const *const transformation,
                                                                       const TranslationCode trans_code,
                                                                       const RotationCode rot_code,
#ifdef VECCORE_CUDA
                                                                       const int id,
#endif
                                                                       VPlacedVolume *const placement) const
{
#ifndef VECCORE_CUDA

  return VolumeFactory::CreateByTransformation<UnplacedBooleanVolume<kIntersection>>(volume, transformation, trans_code,
                                                                                     rot_code,
#ifdef VECCORE_CUDA
                                                                                     id,
#endif
                                                                                     placement);

#else
  // Compiling the above code with nvcc 6.5 faile with the error:
  // nvcc error   : 'ptxas' died due to signal 11 (Invalid memory reference)
  // at least when optimized.
  return nullptr;
#endif
}

VECCORE_ATT_HOST_DEVICE
void TransformedExtent(VPlacedVolume const *pvol, Vector3D<Precision> &aMin, Vector3D<Precision> &aMax)
{
// CUDA does not have min and max in std:: namespace
#ifndef VECCORE_CUDA
  using std::min;
  using std::max;
#endif
  Vector3D<Precision> lower, upper;
  pvol->Extent(lower, upper);
  Vector3D<Precision> delta = upper - lower;
  Precision minx, miny, minz, maxx, maxy, maxz;
  minx         = kInfLength;
  miny         = kInfLength;
  minz         = kInfLength;
  maxx         = -kInfLength;
  maxy         = -kInfLength;
  maxz         = -kInfLength;
  auto *transf = pvol->GetTransformation();
  for (int x = 0; x <= 1; ++x)
    for (int y = 0; y <= 1; ++y)
      for (int z = 0; z <= 1; ++z) {
        Vector3D<Precision> corner;
        corner.x()                            = lower.x() + x * delta.x();
        corner.y()                            = lower.y() + y * delta.y();
        corner.z()                            = lower.z() + z * delta.z();
        Vector3D<Precision> transformedcorner = transf->InverseTransform(corner);
        minx                                  = min(minx, transformedcorner.x());
        miny                                  = min(miny, transformedcorner.y());
        minz                                  = min(minz, transformedcorner.z());
        maxx                                  = max(maxx, transformedcorner.x());
        maxy                                  = max(maxy, transformedcorner.y());
        maxz                                  = max(maxz, transformedcorner.z());
      }
  aMin.x() = minx;
  aMin.y() = miny;
  aMin.z() = minz;
  aMax.x() = maxx;
  aMax.y() = maxy;
  aMax.z() = maxz;
}

template <>
void UnplacedBooleanVolume<kSubtraction>::Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const
{
  Vector3D<Precision> minLeft, maxLeft, minRight, maxRight;
  // ATTENTION: Extent gives coordinates in the reference frame of the callee
  // therefore we have to calculate Extent in THIS frame using:
  TransformedExtent(fBoolean.fLeftVolume, minLeft, maxLeft);
  TransformedExtent(fBoolean.fRightVolume, minRight, maxRight);
  // rather than just
  // fLeftVolume->Extent(minLeft, maxLeft);
  // fRightVolume->Extent(minRight,maxRight);
  aMin = minLeft;
  aMax = maxLeft;
}

template <>
void UnplacedBooleanVolume<kUnion>::Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const
{
  Vector3D<Precision> minLeft, maxLeft, minRight, maxRight;
  // ATTENTION: Extent gives coordinates in the reference frame of the callee
  // therefore we have to calculate Extent in THIS frame using:
  TransformedExtent(fBoolean.fLeftVolume, minLeft, maxLeft);
  TransformedExtent(fBoolean.fRightVolume, minRight, maxRight);
  // rather than just
  // fLeftVolume->Extent(minLeft, maxLeft);
  // fRightVolume->Extent(minRight,maxRight);
  aMin = Vector3D<Precision>(Min(minLeft.x(), minRight.x()), Min(minLeft.y(), minRight.y()),
                             Min(minLeft.z(), minRight.z()));
  aMax = Vector3D<Precision>(Max(maxLeft.x(), maxRight.x()), Max(maxLeft.y(), maxRight.y()),
                             Max(maxLeft.z(), maxRight.z()));
}

template <>
void UnplacedBooleanVolume<kIntersection>::Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const
{
  Vector3D<Precision> minLeft, maxLeft, minRight, maxRight;
  // ATTENTION: Extent gives coordinates in the reference frame of the callee
  // therefore we have to calculate Extent in THIS frame using:
  TransformedExtent(fBoolean.fLeftVolume, minLeft, maxLeft);
  TransformedExtent(fBoolean.fRightVolume, minRight, maxRight);
  // rather than just
  // fLeftVolume->Extent(minLeft, maxLeft);
  // fRightVolume->Extent(minRight,maxRight);
  aMin = Vector3D<Precision>(Max(minLeft.x(), minRight.x()), Max(minLeft.y(), minRight.y()),
                             Max(minLeft.z(), minRight.z()));
  aMax = Vector3D<Precision>(Min(maxLeft.x(), maxRight.x()), Min(maxLeft.y(), maxRight.y()),
                             Min(maxLeft.z(), maxRight.z()));
}

template <>
VECCORE_ATT_HOST_DEVICE
bool UnplacedBooleanVolume<kSubtraction>::Normal(Vector3D<Precision> const &point, Vector3D<Precision> &normal) const
{
  // Compute normal vector to closest surface
  bool valid = false;
  BooleanImplementation<kSubtraction>::NormalKernel(GetStruct(), point, normal, valid);
  return valid;
}

template <>
VECCORE_ATT_HOST_DEVICE
bool UnplacedBooleanVolume<kUnion>::Normal(Vector3D<Precision> const &point, Vector3D<Precision> &normal) const
{
  // Compute normal vector to closest surface
  bool valid = false;
  BooleanImplementation<kUnion>::NormalKernel(GetStruct(), point, normal, valid);
  return valid;
}

template <>
VECCORE_ATT_HOST_DEVICE
bool UnplacedBooleanVolume<kIntersection>::Normal(Vector3D<Precision> const &point, Vector3D<Precision> &normal) const
{
  // Compute normal vector to closest surface
  bool valid = false;
  BooleanImplementation<kIntersection>::NormalKernel(GetStruct(), point, normal, valid);
  return valid;
}

#ifdef VECGEOM_CUDA_INTERFACE

// functions to copy data structures to GPU
template <>
DevicePtr<cuda::VUnplacedVolume> UnplacedBooleanVolume<kUnion>::CopyToGpu(
    DevicePtr<cuda::VUnplacedVolume> const in_gpu_ptr) const
{
  // here we have our recursion:
  // since UnplacedBooleanVolume has pointer members we need to copy/construct those members too
  // very brute force; because this might have been copied already
  // TODO: integrate this into CUDA MGR?

  // use CUDA Manager to lookup GPU pointer
  DevicePtr<cuda::VPlacedVolume> leftgpuptr  = CudaManager::Instance().LookupPlaced(GetLeft());
  DevicePtr<cuda::VPlacedVolume> rightgpuptr = CudaManager::Instance().LookupPlaced(GetRight());

  return CopyToGpuImpl<UnplacedBooleanVolume<kUnion>>(in_gpu_ptr, GetOp(), leftgpuptr, rightgpuptr);
}
template <>
DevicePtr<cuda::VUnplacedVolume> UnplacedBooleanVolume<kUnion>::CopyToGpu() const
{
  return CopyToGpuImpl<UnplacedBooleanVolume<kUnion>>();
}

// functions to copy data structures to GPU
template <>
DevicePtr<cuda::VUnplacedVolume> UnplacedBooleanVolume<kIntersection>::CopyToGpu(
    DevicePtr<cuda::VUnplacedVolume> const in_gpu_ptr) const
{
  // here we have our recursion:
  // since UnplacedBooleanVolume has pointer members we need to copy/construct those members too
  // very brute force; because this might have been copied already
  // TODO: integrate this into CUDA MGR?

  // use CUDA Manager to lookup GPU pointer
  DevicePtr<cuda::VPlacedVolume> leftgpuptr  = CudaManager::Instance().LookupPlaced(GetLeft());
  DevicePtr<cuda::VPlacedVolume> rightgpuptr = CudaManager::Instance().LookupPlaced(GetRight());

  return CopyToGpuImpl<UnplacedBooleanVolume<kIntersection>>(in_gpu_ptr, GetOp(), leftgpuptr, rightgpuptr);
}

template <>
DevicePtr<cuda::VUnplacedVolume> UnplacedBooleanVolume<kIntersection>::CopyToGpu() const
{
  return CopyToGpuImpl<UnplacedBooleanVolume<kIntersection>>();
}

template <>
// functions to copy data structures to GPU
DevicePtr<cuda::VUnplacedVolume> UnplacedBooleanVolume<kSubtraction>::CopyToGpu(
    DevicePtr<cuda::VUnplacedVolume> const in_gpu_ptr) const
{
  // here we have our recursion:
  // since UnplacedBooleanVolume has pointer members we need to copy/construct those members too
  // very brute force; because this might have been copied already
  // TODO: integrate this into CUDA MGR?

  // use CUDA Manager to lookup GPU pointer
  DevicePtr<cuda::VPlacedVolume> leftgpuptr  = CudaManager::Instance().LookupPlaced(GetLeft());
  DevicePtr<cuda::VPlacedVolume> rightgpuptr = CudaManager::Instance().LookupPlaced(GetRight());

  return CopyToGpuImpl<UnplacedBooleanVolume<kSubtraction>>(in_gpu_ptr, GetOp(), leftgpuptr, rightgpuptr);
}

template <>
DevicePtr<cuda::VUnplacedVolume> UnplacedBooleanVolume<kSubtraction>::CopyToGpu() const
{
  return CopyToGpuImpl<UnplacedBooleanVolume<kSubtraction>>();
}

#endif // VECGEOM_CUDA_INTERFACE

} // End impl namespace

#ifdef VECCORE_CUDA

namespace cxx {

template size_t DevicePtr<cuda::UnplacedBooleanVolume<kUnion>>::SizeOf();
template void DevicePtr<cuda::UnplacedBooleanVolume<kUnion>>::Construct(BooleanOperation op,
                                                                        DevicePtr<cuda::VPlacedVolume> left,
                                                                        DevicePtr<cuda::VPlacedVolume> right) const;
template size_t DevicePtr<cuda::UnplacedBooleanVolume<kIntersection>>::SizeOf();
template void DevicePtr<cuda::UnplacedBooleanVolume<kIntersection>>::Construct(
    BooleanOperation op, DevicePtr<cuda::VPlacedVolume> left, DevicePtr<cuda::VPlacedVolume> right) const;
template size_t DevicePtr<cuda::UnplacedBooleanVolume<kSubtraction>>::SizeOf();
template void DevicePtr<cuda::UnplacedBooleanVolume<kSubtraction>>::Construct(
    BooleanOperation op, DevicePtr<cuda::VPlacedVolume> left, DevicePtr<cuda::VPlacedVolume> right) const;

} // End cxx namespace

#endif

} // End namespace vecgeom
