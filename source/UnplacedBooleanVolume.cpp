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
#include "management/CudaManager.h"

namespace vecgeom {

inline namespace VECGEOM_IMPL_NAMESPACE {

template <TranslationCode transCodeT, RotationCode rotCodeT>
VECGEOM_CUDA_HEADER_DEVICE
VPlacedVolume* UnplacedBooleanVolume::Create(
    LogicalVolume const *const logical_volume,
    Transformation3D const *const transformation,
#ifdef VECGEOM_NVCC
    const int id,
#endif
    VPlacedVolume *const placement)
{
    // since this is a static function, we need to get instance of UnplacedBooleanVolume first of all from logical volume
   __attribute__((unused)) const UnplacedBooleanVolume &vol
            = static_cast<const UnplacedBooleanVolume&>( *(logical_volume->unplaced_volume()) );

   if( vol.GetOp() == kSubtraction ) {
      return CreateSpecializedWithPlacement<SpecializedBooleanVolume<kSubtraction, transCodeT, rotCodeT> >(
      logical_volume, transformation, 
#ifdef VECGEOM_NVCC
      id,
#endif 
      placement); // TODO: add bounding box?
   }
   else if ( vol.GetOp() == kUnion ) {
      return CreateSpecializedWithPlacement<SpecializedBooleanVolume<kUnion, transCodeT, rotCodeT> >(
              logical_volume, transformation,
        #ifdef VECGEOM_NVCC
              id,
        #endif
              placement); // TODO: add bounding box?
   }
   else if ( vol.GetOp() == kIntersection ){
      return CreateSpecializedWithPlacement<SpecializedBooleanVolume<kIntersection, transCodeT, rotCodeT> >(
              logical_volume, transformation,
        #ifdef VECGEOM_NVCC
              id,
        #endif
              placement); // TODO: add bounding box?
   }
   return nullptr;
}


VECGEOM_CUDA_HEADER_DEVICE
VPlacedVolume* UnplacedBooleanVolume::SpecializedVolume(
    LogicalVolume const *const volume,
    Transformation3D const *const transformation,
    const TranslationCode trans_code, const RotationCode rot_code,
#ifdef VECGEOM_NVCC
    const int id,
#endif
    VPlacedVolume *const placement) const {
   return VolumeFactory::CreateByTransformation<UnplacedBooleanVolume>(
    volume, transformation, trans_code, rot_code,
#ifdef VECGEOM_NVCC
    id,
#endif
    placement);
}




#ifdef VECGEOM_CUDA_INTERFACE

// functions to copy data structures to GPU
DevicePtr<cuda::VUnplacedVolume> UnplacedBooleanVolume::CopyToGpu(
   DevicePtr<cuda::VUnplacedVolume> const in_gpu_ptr) const
{
    // here we have our recursion:
    // since UnplacedBooleanVolume has pointer members we need to copy/construct those members too
    // very brute force; because this might have been copied already
    // TODO: integrate this into CUDA MGR?

    // use CUDA Manager to lookup GPU pointer
    DevicePtr<cuda::VPlacedVolume> leftgpuptr = CudaManager::Instance().LookupPlaced(fLeftVolume);
    DevicePtr<cuda::VPlacedVolume> rightgpuptr = CudaManager::Instance().LookupPlaced(fRightVolume);

    return CopyToGpuImpl<UnplacedBooleanVolume>(in_gpu_ptr, fOp, leftgpuptr, rightgpuptr);
}

DevicePtr<cuda::VUnplacedVolume> UnplacedBooleanVolume::CopyToGpu() const
{
   return CopyToGpuImpl<UnplacedBooleanVolume>();
}

#endif // VECGEOM_CUDA_INTERFACE

} // End impl namespace

#ifdef VECGEOM_NVCC

namespace cxx {

template size_t DevicePtr<cuda::UnplacedBooleanVolume>::SizeOf();
template void DevicePtr<cuda::UnplacedBooleanVolume>::Construct(
    BooleanOperation op,
    DevicePtr<cuda::VPlacedVolume> left,
    DevicePtr<cuda::VPlacedVolume> right) const;

} // End cxx namespace

#endif


} // End namespace vecgeom


