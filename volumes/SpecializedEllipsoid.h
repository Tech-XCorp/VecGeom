// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.
// For the full list of authors see CONTRIBUTORS.txt and `git log`.

/// Declaration of the specialized ellipsoid volume
/// @file volumes/SpecializedEllipsoid.h
/// @author Evgueni Tcherniaev

#ifndef VECGEOM_VOLUMES_SPECIALIZEDELLIPSOID_H_
#define VECGEOM_VOLUMES_SPECIALIZEDELLIPSOID_H_

#include "base/Global.h"

#include "volumes/kernel/EllipsoidImplementation.h"
#include "volumes/PlacedEllipsoid.h"
#include "volumes/SpecializedPlacedVolImplHelper.h"
#include "volumes/UnplacedEllipsoid.h"

#include <stdio.h>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

template <TranslationCode transCodeT, RotationCode rotCodeT>
using SpecializedEllipsoid = SIMDSpecializedVolImplHelper<EllipsoidImplementation, transCodeT, rotCodeT>;

using SimpleEllipsoid = SpecializedEllipsoid<translation::kGeneric, rotation::kGeneric>;
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_SPECIALIZEDELLIPSOID_H_
