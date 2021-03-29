/// \file Transformation3D.h
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#ifndef VECGEOM_BASE_TRANSFORMATION3D_H_
#define VECGEOM_BASE_TRANSFORMATION3D_H_

#include "VecGeom/base/Global.h"

#include "VecGeom/base/Vector3D.h"
#include "VecGeom/backend/scalar/Backend.h"
#ifdef VECGEOM_ENABLE_CUDA
#include "VecGeom/backend/cuda/Interface.h"
#endif

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>

#ifdef VECGEOM_ROOT
class TGeoMatrix;
#endif

typedef int RotationCode;
typedef int TranslationCode;

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class Transformation3D;);

inline namespace VECGEOM_IMPL_NAMESPACE {

#ifndef VECCORE_CUDA
}
namespace cuda {
class Transformation3D;
}
inline namespace VECGEOM_IMPL_NAMESPACE {
// class vecgeom::cuda::Transformation3D;
#endif

class Transformation3D {

private:
  // TODO: it might be better to directly store this in terms of Vector3D<Precision> !!
  // and would allow for higher level abstraction
  Precision fTranslation[3];
  Precision fRotation[9];

public:
  VECCORE_ATT_HOST_DEVICE
  constexpr Transformation3D()
      : fTranslation{0., 0., 0.}, fRotation{1., 0., 0., 0., 1., 0., 0., 0., 1.}
  {
  }

  /**
   * Constructor for translation only.
   * @param tx Translation in x-coordinate.
   * @param ty Translation in y-coordinate.
   * @param tz Translation in z-coordinate.
   */
  VECCORE_ATT_HOST_DEVICE
  Transformation3D(const Precision tx, const Precision ty, const Precision tz)
      : fTranslation{tx, ty, tz}, fRotation{1., 0., 0., 0., 1., 0., 0., 0., 1.} 
  {
  }

  /**
   * @param tx Translation in x-coordinate.
   * @param ty Translation in y-coordinate.
   * @param tz Translation in z-coordinate.
   * @param phi Rotation angle about z-axis.
   * @param theta Rotation angle about new y-axis.
   * @param psi Rotation angle about new z-axis.
   */
  VECCORE_ATT_HOST_DEVICE
  Transformation3D(const Precision tx, const Precision ty, const Precision tz, const Precision phi,
                   const Precision theta, const Precision psi);

  /**
   * Constructor to manually set each entry. Used when converting from different
   * geometry.
   */
  VECCORE_ATT_HOST_DEVICE
  Transformation3D(const Precision tx, const Precision ty, const Precision tz, const Precision r0, const Precision r1,
                   const Precision r2, const Precision r3, const Precision r4, const Precision r5, const Precision r6,
                   const Precision r7, const Precision r8);

  /**
   * Constructor copying the translation and rotation from memory
   * geometry.
   */
  VECCORE_ATT_HOST_DEVICE
  Transformation3D(const Precision *trans, const Precision *rot, bool has_trans, bool has_rot)
  {
    this->Set(trans, rot, has_trans, has_rot);
  }

  /**
   * Constructor for a rotation based on a given direction
   * @param axis direction of the new z axis
   * @param inverse if true the origial axis will be rotated into (0,0,u)
                    if false a vector (0,0,u) will be rotated into the original axis
   */
  VECCORE_ATT_HOST_DEVICE
  Transformation3D(const Vector3D<Precision> &axis, bool inverse = true);

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Transformation3D(Transformation3D const &other);

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Transformation3D &operator=(Transformation3D const &rhs);

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  bool operator==(Transformation3D const &rhs) const;

  VECCORE_ATT_HOST_DEVICE
  ~Transformation3D() {}

  VECCORE_ATT_HOST_DEVICE
  void Clear()
  {
    fTranslation[0] = 0.;
    fTranslation[1] = 0.;
    fTranslation[2] = 0.;
    fRotation[0]    = 1.;
    fRotation[1]    = 0.;
    fRotation[2]    = 0.;
    fRotation[3]    = 0.;
    fRotation[4]    = 1.;
    fRotation[5]    = 0.;
    fRotation[6]    = 0.;
    fRotation[7]    = 0.;
    fRotation[8]    = 1.;
  }

  int MemorySize() const { return sizeof(*this); }

  VECCORE_ATT_HOST_DEVICE
  void FixZeroes()
  {
    for (unsigned int i = 0; i < 9; ++i) {
      if (std::abs(fRotation[i]) < vecgeom::kTolerance) fRotation[i] = 0.;
    }
    for (unsigned int i = 0; i < 3; ++i) {
      if (std::abs(fTranslation[i]) < vecgeom::kTolerance) fTranslation[i] = 0.;
    }
  }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Vector3D<Precision> Translation() const
  {
    return Vector3D<Precision>(fTranslation[0], fTranslation[1], fTranslation[2]);
  }

  /**
   * No safety against faulty indexing.
   * @param index Index of translation entry in the range [0-2].
   */
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision Translation(const int index) const { return fTranslation[index]; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision const *Rotation() const { return fRotation; }

  /**
   * No safety against faulty indexing.
   * \param index Index of rotation entry in the range [0-8].
   */
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Precision Rotation(const int index) const { return fRotation[index]; }

  VECCORE_ATT_HOST_DEVICE
  void Print() const;

  // print to a stream
  void Print(std::ostream &) const;

  // Mutators

  VECCORE_ATT_HOST_DEVICE
  void SetTranslation(const Precision tx, const Precision ty, const Precision tz);

  VECCORE_ATT_HOST_DEVICE
  void SetTranslation(Vector3D<Precision> const &vec);

  VECCORE_ATT_HOST_DEVICE
  void SetRotation(const Precision phi, const Precision theta, const Precision psi);

  VECCORE_ATT_HOST_DEVICE
  void SetRotation(Vector3D<Precision> const &vec);

  VECCORE_ATT_HOST_DEVICE
  void SetRotation(const Precision rot0, const Precision rot1, const Precision rot2, const Precision rot3,
                   const Precision rot4, const Precision rot5, const Precision rot6, const Precision rot7,
                   const Precision rot8);

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  void Set(const Precision *trans, const Precision *rot, bool has_trans, bool has_rot)
  {
    // Memory copy translation and rotation components.
    // Avoid memcpy that gives warning: non-null argument 2 expected
    // Avoid copy by de-referencing trans/rot which may not have alignof(Precision)
    auto trans_src  = reinterpret_cast<const char *>(trans);
    auto trans_dest = reinterpret_cast<char *>(fTranslation);
    auto rot_src    = reinterpret_cast<const char *>(rot);
    auto rot_dest   = reinterpret_cast<char *>(fRotation);
    for (size_t i = 0; i < 3 * sizeof(Precision); ++i)
      trans_dest[i] = trans_src[i];
    for (size_t i = 0; i < 9 * sizeof(Precision); ++i)
      rot_dest[i] = rot_src[i];
  }

private:
  template <typename InputType>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  void DoRotation(Vector3D<InputType> const &master, Vector3D<InputType> &local) const;

  template <typename InputType>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  void DoTranslation(Vector3D<InputType> const &master, Vector3D<InputType> &local) const;

  template <bool vectortransform, typename InputType>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  void InverseTransformKernel(Vector3D<InputType> const &local, Vector3D<InputType> &master) const;

public:
  // Transformation interface

  template <typename InputType>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  void Transform(Vector3D<InputType> const &master, Vector3D<InputType> &local) const;

  template <typename InputType>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  Vector3D<InputType> Transform(Vector3D<InputType> const &master) const;

  template <typename InputType>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  void TransformDirection(Vector3D<InputType> const &master, Vector3D<InputType> &local) const;

  template <typename InputType>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  Vector3D<InputType> TransformDirection(Vector3D<InputType> const &master) const;

  /** The inverse transformation ( aka LocalToMaster ) of an object transform like a point
   *  this does not need to currently template on placement since such a transformation is much less used
   */
  template <typename InputType>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  void InverseTransform(Vector3D<InputType> const &local, Vector3D<InputType> &master) const;

  template <typename InputType>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  Vector3D<InputType> InverseTransform(Vector3D<InputType> const &local) const;

  /** The inverse transformation of an object transforming like a vector */
  template <typename InputType>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  void InverseTransformDirection(Vector3D<InputType> const &master, Vector3D<InputType> &local) const;

  template <typename InputType>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  Vector3D<InputType> InverseTransformDirection(Vector3D<InputType> const &master) const;

  /** compose transformations - multiply transformations */
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  void MultiplyFromRight(Transformation3D const &rhs);

  /** compose transformations - multiply transformations */
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  void CopyFrom(Transformation3D const &rhs)
  {
    // not sure this compiles under CUDA
    copy(&rhs, &rhs + 1, this);
  }

  // stores the inverse of this matrix into inverse
  // taken from CLHEP implementation
  VECCORE_ATT_HOST_DEVICE
  void Inverse(Transformation3D &inverse) const
  {
    double xx_ = fRotation[0];
    double zz_ = fRotation[8];
    double yy_ = fRotation[4];
    double xy_ = fRotation[1];
    double xz_ = fRotation[2];
    double yx_ = fRotation[3];
    double yz_ = fRotation[5];
    double zx_ = fRotation[6];
    double zy_ = fRotation[7];
    double dx_ = fTranslation[0];
    double dy_ = fTranslation[1];
    double dz_ = fTranslation[2];

    double detxx = yy_ * zz_ - yz_ * zy_;
    double detxy = yx_ * zz_ - yz_ * zx_;
    double detxz = yx_ * zy_ - yy_ * zx_;
    double det   = xx_ * detxx - xy_ * detxy + xz_ * detxz;
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    if (det == 0) {
      std::cerr << "Transform3D::inverse error: zero determinant" << std::endl;
    }
#endif
    det = 1. / det;
    detxx *= det;
    detxy *= det;
    detxz *= det;
    double detyx            = (xy_ * zz_ - xz_ * zy_) * det;
    double detyy            = (xx_ * zz_ - xz_ * zx_) * det;
    double detyz            = (xx_ * zy_ - xy_ * zx_) * det;
    double detzx            = (xy_ * yz_ - xz_ * yy_) * det;
    double detzy            = (xx_ * yz_ - xz_ * yx_) * det;
    double detzz            = (xx_ * yy_ - xy_ * yx_) * det;
    inverse.fRotation[0]    = detxx;
    inverse.fRotation[1]    = -detyx;
    inverse.fRotation[2]    = detzx;
    inverse.fTranslation[0] = -detxx * dx_ + detyx * dy_ - detzx * dz_;
    inverse.fRotation[3] = -detxy, inverse.fRotation[4] = detyy, inverse.fRotation[5] = -detzy,
    inverse.fTranslation[1] = detxy * dx_ - detyy * dy_ + detzy * dz_;
    inverse.fRotation[6] = detxz, inverse.fRotation[7] = -detyz, inverse.fRotation[8] = detzz,
    inverse.fTranslation[2] = -detxz * dx_ + detyz * dy_ - detzz * dz_;
  }

  // Utility and CUDA

#ifdef VECGEOM_CUDA_INTERFACE
  size_t DeviceSizeOf() const { return DevicePtr<cuda::Transformation3D>::SizeOf(); }
  DevicePtr<cuda::Transformation3D> CopyToGpu() const;
  DevicePtr<cuda::Transformation3D> CopyToGpu(DevicePtr<cuda::Transformation3D> const gpu_ptr) const;
#endif

#ifdef VECGEOM_ROOT
  // function to convert this transformation to a TGeo transformation
  // mainly used for the benchmark comparisons with ROOT
  TGeoMatrix *ConvertToTGeoMatrix() const;
#endif

public:
  static const Transformation3D kIdentity;

}; // End class Transformation3D

VECCORE_ATT_HOST_DEVICE
Transformation3D::Transformation3D(Transformation3D const &other)
{
  *this = other;
}

VECCORE_ATT_HOST_DEVICE
Transformation3D &Transformation3D::operator=(Transformation3D const &rhs)
{
  copy(rhs.fTranslation, rhs.fTranslation + 3, fTranslation);
  copy(rhs.fRotation, rhs.fRotation + 9, fRotation);
  return *this;
}

VECCORE_ATT_HOST_DEVICE
bool Transformation3D::operator==(Transformation3D const &rhs) const
{
  return equal(fTranslation, fTranslation + 3, rhs.fTranslation) && equal(fRotation, fRotation + 9, rhs.fRotation);
}

/**
 * Rotates a vector to this transformation's frame of reference.
 * Templates on the RotationCode generated by GenerateTranslationCode() to
 * perform specialized rotation.
 * \sa GenerateTranslationCode()
 * \param master Vector in original frame of reference.
 * \param local Output vector rotated to the new frame of reference.
 */
template <typename InputType>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
void Transformation3D::DoRotation(Vector3D<InputType> const &master, Vector3D<InputType> &local) const
{
  local[0] = master[0] * fRotation[0];
  local[1] = master[0] * fRotation[1];
  local[2] = master[0] * fRotation[2];
  local[0] += master[1] * fRotation[3];
  local[1] += master[1] * fRotation[4];
  local[2] += master[1] * fRotation[5];
  local[0] += master[2] * fRotation[6];
  local[1] += master[2] * fRotation[7];
  local[2] += master[2] * fRotation[8];
}

template <typename InputType>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
void Transformation3D::DoTranslation(Vector3D<InputType> const &master, Vector3D<InputType> &local) const
{

  local[0] = master[0] - fTranslation[0];
  local[1] = master[1] - fTranslation[1];
  local[2] = master[2] - fTranslation[2];
}

/**
 * Transform a point to the local reference frame.
 * \param master Point to be transformed.
 * \param local Output destination. Should never be the same as the input
 *              vector!
 */
template <typename InputType>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
void Transformation3D::Transform(Vector3D<InputType> const &master, Vector3D<InputType> &local) const
{
  Vector3D<InputType> tmp;
  DoTranslation(master, tmp);
  DoRotation(tmp, local);
}

/**
 * Since transformation cannot be done in place, allows the transformed vector
 * to be constructed by Transform directly.
 * \param master Point to be transformed.
 * \return Newly constructed Vector3D with the transformed coordinates.
 */
template <typename InputType>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
Vector3D<InputType> Transformation3D::Transform(Vector3D<InputType> const &master) const
{
  Vector3D<InputType> local;
  Transform(master, local);
  return local;
}

template <bool transform_direction, typename InputType>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
void Transformation3D::InverseTransformKernel(Vector3D<InputType> const &local, Vector3D<InputType> &master) const
{

  // we are just doing the full stuff here ( LocalToMaster is less critical
  // than other way round )

  if (transform_direction) {
    master[0] = local[0] * fRotation[0];
    master[0] += local[1] * fRotation[1];
    master[0] += local[2] * fRotation[2];
    master[1] = local[0] * fRotation[3];
    master[1] += local[1] * fRotation[4];
    master[1] += local[2] * fRotation[5];
    master[2] = local[0] * fRotation[6];
    master[2] += local[1] * fRotation[7];
    master[2] += local[2] * fRotation[8];
  } else {
    master[0] = fTranslation[0];
    master[0] += local[0] * fRotation[0];
    master[0] += local[1] * fRotation[1];
    master[0] += local[2] * fRotation[2];
    master[1] = fTranslation[1];
    master[1] += local[0] * fRotation[3];
    master[1] += local[1] * fRotation[4];
    master[1] += local[2] * fRotation[5];
    master[2] = fTranslation[2];
    master[2] += local[0] * fRotation[6];
    master[2] += local[1] * fRotation[7];
    master[2] += local[2] * fRotation[8];
  }
}

template <typename InputType>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
void Transformation3D::InverseTransform(Vector3D<InputType> const &local, Vector3D<InputType> &master) const
{
  InverseTransformKernel<false, InputType>(local, master);
}

template <typename InputType>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
Vector3D<InputType> Transformation3D::InverseTransform(Vector3D<InputType> const &local) const
{
  Vector3D<InputType> tmp;
  InverseTransform(local, tmp);
  return tmp;
}

template <typename InputType>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
void Transformation3D::InverseTransformDirection(Vector3D<InputType> const &local, Vector3D<InputType> &master) const
{
  InverseTransformKernel<true, InputType>(local, master);
}

template <typename InputType>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
Vector3D<InputType> Transformation3D::InverseTransformDirection(Vector3D<InputType> const &local) const
{
  Vector3D<InputType> tmp;
  InverseTransformDirection(local, tmp);
  return tmp;
}

VECCORE_ATT_HOST_DEVICE
VECGEOM_FORCE_INLINE
void Transformation3D::MultiplyFromRight(Transformation3D const &rhs)
{
  // TODO: this code should directly operator on Vector3D and Matrix3D

  // ideal for fused multiply add
  fTranslation[0] += fRotation[0] * rhs.fTranslation[0];
  fTranslation[0] += fRotation[1] * rhs.fTranslation[1];
  fTranslation[0] += fRotation[2] * rhs.fTranslation[2];

  fTranslation[1] += fRotation[3] * rhs.fTranslation[0];
  fTranslation[1] += fRotation[4] * rhs.fTranslation[1];
  fTranslation[1] += fRotation[5] * rhs.fTranslation[2];

  fTranslation[2] += fRotation[6] * rhs.fTranslation[0];
  fTranslation[2] += fRotation[7] * rhs.fTranslation[1];
  fTranslation[2] += fRotation[8] * rhs.fTranslation[2];

  Precision tmpx = fRotation[0];
  Precision tmpy = fRotation[1];
  Precision tmpz = fRotation[2];

  // first row of matrix
  fRotation[0] = tmpx * rhs.fRotation[0];
  fRotation[1] = tmpx * rhs.fRotation[1];
  fRotation[2] = tmpx * rhs.fRotation[2];
  fRotation[0] += tmpy * rhs.fRotation[3];
  fRotation[1] += tmpy * rhs.fRotation[4];
  fRotation[2] += tmpy * rhs.fRotation[5];
  fRotation[0] += tmpz * rhs.fRotation[6];
  fRotation[1] += tmpz * rhs.fRotation[7];
  fRotation[2] += tmpz * rhs.fRotation[8];

  tmpx = fRotation[3];
  tmpy = fRotation[4];
  tmpz = fRotation[5];

  // second row of matrix
  fRotation[3] = tmpx * rhs.fRotation[0];
  fRotation[4] = tmpx * rhs.fRotation[1];
  fRotation[5] = tmpx * rhs.fRotation[2];
  fRotation[3] += tmpy * rhs.fRotation[3];
  fRotation[4] += tmpy * rhs.fRotation[4];
  fRotation[5] += tmpy * rhs.fRotation[5];
  fRotation[3] += tmpz * rhs.fRotation[6];
  fRotation[4] += tmpz * rhs.fRotation[7];
  fRotation[5] += tmpz * rhs.fRotation[8];

  tmpx = fRotation[6];
  tmpy = fRotation[7];
  tmpz = fRotation[8];

  // third row of matrix
  fRotation[6] = tmpx * rhs.fRotation[0];
  fRotation[7] = tmpx * rhs.fRotation[1];
  fRotation[8] = tmpx * rhs.fRotation[2];
  fRotation[6] += tmpy * rhs.fRotation[3];
  fRotation[7] += tmpy * rhs.fRotation[4];
  fRotation[8] += tmpy * rhs.fRotation[5];
  fRotation[6] += tmpz * rhs.fRotation[6];
  fRotation[7] += tmpz * rhs.fRotation[7];
  fRotation[8] += tmpz * rhs.fRotation[8];
}

/**
 * Only transforms by rotation, ignoring the translation part. This is useful
 * when transforming directions.
 * \param master Point to be transformed.
 * \param local Output destination of transformation.
 */
template <typename InputType>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
void Transformation3D::TransformDirection(Vector3D<InputType> const &master, Vector3D<InputType> &local) const
{
  DoRotation(master, local);
}

/**
 * Since transformation cannot be done in place, allows the transformed vector
 * to be constructed by TransformDirection directly.
 * \param master Point to be transformed.
 * \return Newly constructed Vector3D with the transformed coordinates.
 */
template <typename InputType>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
Vector3D<InputType> Transformation3D::TransformDirection(Vector3D<InputType> const &master) const
{

  Vector3D<InputType> local;
  TransformDirection(master, local);
  return local;
}

std::ostream &operator<<(std::ostream &os, Transformation3D const &trans);
}
} // namespace vecgeom

#endif // VECGEOM_BASE_TRANSFORMATION3D_H_
