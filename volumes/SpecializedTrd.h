/// @file SpecializedTrd.h
/// @author Georgios Bitzes (georgios.bitzes@cern.ch)

#ifndef VECGEOM_VOLUMES_SPECIALIZEDTRD_H_
#define VECGEOM_VOLUMES_SPECIALIZEDTRD_H_

#include "base/Global.h"

#include "volumes/kernel/TrdImplementation.h"
#include "volumes/PlacedTrd.h"
#include "volumes/ShapeImplementationHelper.h"
#include <stdio.h>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

template <TranslationCode transCodeT, RotationCode rotCodeT, typename trdTypeT>
using SpecializedTrd = ShapeImplementationHelper<TrdImplementation<transCodeT, rotCodeT, trdTypeT>>;

using SimpleTrd = SpecializedTrd<translation::kGeneric, rotation::kGeneric, TrdTypes::UniversalTrd>;
}
} // End global namespace

#endif // VECGEOM_VOLUMES_SPECIALIZEDTRD_H_
