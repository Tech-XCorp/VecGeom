/*
 * BoxStruct.h
 *
 *  Created on: 09.10.2015
 *      Author: swenzel
 */

#ifndef VECGEOM_VOLUMES_BOXSTRUCT_H_
#define VECGEOM_VOLUMES_BOXSTRUCT_H_
#include "base/Global.h"

namespace vecgeom {

inline namespace VECGEOM_IMPL_NAMESPACE {

// we do something crazy : a plain struct without member functions to encapsulate just the parameters
// of a box
template <typename T = double>
struct BoxStruct {
  Vector3D<T> fDimensions; //<the HALF lengths of the box

  VECGEOM_CUDA_HEADER_BOTH
  BoxStruct(Vector3D<T> const &dim) : fDimensions(dim) {}

  VECGEOM_CUDA_HEADER_BOTH
  BoxStruct(const T dx, const T dy, const T dz) : fDimensions(dx, dy, dz) {}
};
}
} // end

#endif
