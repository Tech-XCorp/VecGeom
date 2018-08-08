/*
 * SpecializedPolycone.h
 *
 *  Created on: Dec 8, 2014
 *      Author: swenzel
 */

#ifndef VECGEOM_VOLUMES_SPECIALIZEDPOLYCONE_H_
#define VECGEOM_VOLUMES_SPECIALIZEDPOLYCONE_H_

#include "base/Global.h"

#include "volumes/kernel/PolyconeImplementation.h"
#include "volumes/PlacedPolycone.h"
#include "volumes/SpecializedPlacedVolImplHelper.h"
#include "volumes/UnplacedPolycone.h"

#include <stdio.h>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

template <TranslationCode transCodeT, RotationCode rotCodeT, typename polyconeTypeT>
using SpecializedPolycone = LoopSpecializedVolImplHelper<PolyconeImplementation<polyconeTypeT>, transCodeT, rotCodeT>;

using SimplePolycone = SpecializedPolycone<translation::kGeneric, rotation::kGeneric, ConeTypes::UniversalCone>;

template <typename Type>
template <TranslationCode transCodeT, RotationCode rotCodeT>
VECCORE_ATT_DEVICE
VPlacedVolume *SUnplacedPolycone<Type>::Create(LogicalVolume const *const logical_volume,
                                               Transformation3D const *const transformation,
#ifdef VECCORE_CUDA
                                               const int id,
#endif
                                               VPlacedVolume *const placement)
{
  (void)placement;
  return new SpecializedPolycone<transCodeT, rotCodeT, Type>(logical_volume, transformation
#ifdef VECCORE_CUDA
                                                             ,
                                                             id
#endif
  );
}

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif /* VECGEOM_VOLUMES_SPECIALIZEDPOLYCONE_H_ */
