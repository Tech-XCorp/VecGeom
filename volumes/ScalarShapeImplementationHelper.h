/// \file ScalarShapeImplementationHelper.h

#ifndef VECGEOM_VOLUMES_SCALARSHAPEIMPLEMENTATIONHELPER_H_
#define VECGEOM_VOLUMES_SCALARSHAPEIMPLEMENTATIONHELPER_H_

#include "base/Global.h"
#include "base/SOA3D.h"
#include "volumes/PlacedBox.h"

#include <algorithm>
#ifdef VECGEOM_DISTANCE_DEBUG
#include "volumes/utilities/ResultComparator.h"
#endif

#include <VecCore/VecCore>

namespace vecgeom {

VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE(class, ScalarShapeImplementationHelper, typename);

inline namespace VECGEOM_IMPL_NAMESPACE {

/**
 * A helper class implementing "repetetive" dispatching of high level interfaces to
 * actual implementations
 *
 * In contrast to the ordinary ShapeImplementationHelper,
 * the ScalarShapeImplementatioHelper
 * does not explicitely use vectorization; Hence the multi-particle interfaces
 * are dispatched to loops over scalar implementations
 *
 */
template <typename Specialization>
class ScalarShapeImplementationHelper : public Specialization::PlacedShape_t {

  using Real_v           = vecgeom::VectorBackend::Real_v;
  using PlacedShape_t    = typename Specialization::PlacedShape_t;
  using UnplacedShape_t  = typename Specialization::UnplacedShape_t;
  using Helper_t         = ScalarShapeImplementationHelper<Specialization>;
  using Implementation_t = Specialization;

public:
#ifndef VECGEOM_NVCC

  ScalarShapeImplementationHelper(char const *const label, LogicalVolume const *const logical_volume,
                                  Transformation3D const *const transformation, PlacedBox const *const boundingBox)
      : PlacedShape_t(label, logical_volume, transformation, boundingBox)
  {
  }

  ScalarShapeImplementationHelper(char const *const label, LogicalVolume const *const logical_volume,
                                  Transformation3D const *const transformation)
      : ScalarShapeImplementationHelper(label, logical_volume, transformation,
                                        details::UseIfSameType<PlacedShape_t, PlacedBox>::Get(this))
  {
  }

  ScalarShapeImplementationHelper(char const *const label, LogicalVolume *const logical_volume,
                                  Transformation3D const *const transformation, PlacedBox const *const boundingBox)
      : PlacedShape_t(label, logical_volume, transformation, boundingBox)
  {
  }

  ScalarShapeImplementationHelper(char const *const label, LogicalVolume *const logical_volume,
                                  Transformation3D const *const transformation)
      : ScalarShapeImplementationHelper(label, logical_volume, transformation,
                                        details::UseIfSameType<PlacedShape_t, PlacedBox>::Get(this))
  {
  }

  ScalarShapeImplementationHelper(LogicalVolume const *const logical_volume,
                                  Transformation3D const *const transformation, PlacedBox const *const boundingBox)
      : ScalarShapeImplementationHelper("", logical_volume, transformation, boundingBox)
  {
  }

  ScalarShapeImplementationHelper(LogicalVolume const *const logical_volume,
                                  Transformation3D const *const transformation)
      : ScalarShapeImplementationHelper("", logical_volume, transformation)
  {
  }

  template <typename... ArgTypes>
  ScalarShapeImplementationHelper(char const *const label, ArgTypes... params)
      : ScalarShapeImplementationHelper(label, new LogicalVolume(new UnplacedShape_t(params...)),
                                        &Transformation3D::kIdentity)
  {
  }

#else // Compiling for CUDA

  __device__ ScalarShapeImplementationHelper(LogicalVolume const *const logical_volume,
                                             Transformation3D const *const transformation,
                                             PlacedBox const *const boundingBox, const int id)
      : PlacedShape_t(logical_volume, transformation, boundingBox, id)
  {
  }

  __device__ ScalarShapeImplementationHelper(LogicalVolume const *const logical_volume,
                                             Transformation3D const *const transformation, const int id)
      : PlacedShape_t(logical_volume, transformation, details::UseIfSameType<PlacedShape_t, PlacedBox>::Get(this), id)
  {
  }
#endif
  using PlacedShape_t::SafetyToIn;
  using PlacedShape_t::SafetyToOut;
  using PlacedShape_t::DistanceToIn;
  using PlacedShape_t::DistanceToOut;

  virtual int memory_size() const override { return sizeof(*this); }

  VECGEOM_CUDA_HEADER_BOTH
  virtual void PrintType() const override { Specialization::PrintType(); }

  virtual void PrintType(std::ostream &os) const override { Specialization::PrintType(os); }
  virtual void PrintImplementationType(std::ostream &os) const override { Specialization::PrintImplementationType(os); }
  virtual void PrintUnplacedType(std::ostream &os) const override { Specialization::PrintUnplacedType(os); }

#ifdef VECGEOM_CUDA_INTERFACE

  virtual size_t DeviceSizeOf() const override { return DevicePtr<CudaType_t<Helper_t>>::SizeOf(); }

  DevicePtr<cuda::VPlacedVolume> CopyToGpu(DevicePtr<cuda::LogicalVolume> const logical_volume,
                                           DevicePtr<cuda::Transformation3D> const transform,
                                           DevicePtr<cuda::VPlacedVolume> const in_gpu_ptr) const override
  {
    DevicePtr<CudaType_t<Helper_t>> gpu_ptr(in_gpu_ptr);
    gpu_ptr.Construct(logical_volume, transform, DevicePtr<cuda::PlacedBox>(), this->id());
    CudaAssertError();
    // Need to go via the void* because the regular c++ compilation
    // does not actually see the declaration for the cuda version
    // (and thus can not determine the inheritance).
    return DevicePtr<cuda::VPlacedVolume>((void *)gpu_ptr);
  }

  DevicePtr<cuda::VPlacedVolume> CopyToGpu(DevicePtr<cuda::LogicalVolume> const logical_volume,
                                           DevicePtr<cuda::Transformation3D> const transform) const override
  {
    DevicePtr<CudaType_t<Helper_t>> gpu_ptr;
    gpu_ptr.Allocate();
    return CopyToGpu(logical_volume, transform, DevicePtr<cuda::VPlacedVolume>((void *)gpu_ptr));
  }

#endif // VECGEOM_CUDA_INTERFACE

  VECGEOM_CUDA_HEADER_BOTH
  virtual EnumInside Inside(Vector3D<Precision> const &point) const override
  {
    Inside_t output = EInside::kOutside;
    Specialization::template Inside<kScalar>(*this->GetUnplacedVolume(), *this->GetTransformation(), point, output);
    return (EnumInside)output;
  }

  VECGEOM_CUDA_HEADER_BOTH
  virtual bool Contains(Vector3D<Precision> const &point) const override
  {
    bool output = false;
    Vector3D<Precision> localPoint;
    Specialization::template Contains<kScalar>(*this->GetUnplacedVolume(), *this->GetTransformation(), point,
                                               localPoint, output);
    return output;
  }

  VECGEOM_CUDA_HEADER_BOTH
  virtual bool Contains(Vector3D<Precision> const &point, Vector3D<Precision> &localPoint) const override
  {
    bool output = false;
    Specialization::template Contains<kScalar>(*this->GetUnplacedVolume(), *this->GetTransformation(), point,
                                               localPoint, output);

#ifdef VECGEOM_DISTANCE_DEBUG
    DistanceComparator::CompareUnplacedContains(this, output, localPoint);
#endif

    return output;
  }

  VECGEOM_CUDA_HEADER_BOTH
  virtual bool UnplacedContains(Vector3D<Precision> const &point) const override
  {
    bool output = false;
    Specialization::template UnplacedContains<kScalar>(*this->GetUnplacedVolume(), point, output);

#ifdef VECGEOM_DISTANCE_DEBUG
    DistanceComparator::CompareUnplacedContains(this, output, point);
#endif

    return output;
  }

  VECGEOM_CUDA_HEADER_BOTH
  virtual Precision DistanceToIn(Vector3D<Precision> const &point, Vector3D<Precision> const &direction,
                                 const Precision stepMax = kInfinity) const override
  {
    Precision output = kInfinity;
    Specialization::template DistanceToIn<kScalar>(*this->GetUnplacedVolume(), *this->GetTransformation(), point,
                                                   direction, stepMax, output);

#ifdef VECGEOM_REPLACE_USOLIDS
    // apply USolids convention: convert negative values to zero
    MaskedAssign(output < kHalfTolerance, 0., &output);
#else
    // avoid distance values within tolerance
    MaskedAssign(Abs(output) < kHalfTolerance, 0., &output);
#endif

#ifdef VECGEOM_DISTANCE_DEBUG
    DistanceComparator::CompareDistanceToIn(this, output, point, direction, stepMax);
#endif
    return output;
  }

  VECGEOM_CUDA_HEADER_BOTH
  virtual Precision DistanceToOut(Vector3D<Precision> const &point, Vector3D<Precision> const &direction,
                                  const Precision stepMax = kInfinity) const override
  {
    Precision output = kInfinity;
    Specialization::template DistanceToOut<kScalar>(*this->GetUnplacedVolume(), point, direction, stepMax, output);

#ifdef VECGEOM_REPLACE_USOLIDS
    // apply USolids convention: convert negative values to zero
    MaskedAssign(output < kHalfTolerance, 0., &output);
#else
    // avoid distance values within tolerance
    MaskedAssign(Abs(output) < kHalfTolerance, 0., &output);
#endif

// detect -inf responses which are often an indication for a real bug
#ifndef VECGEOM_NVCC
    assert(!((output < 0.) && std::isinf(output)));
#endif

    return output;
  }

  VECGEOM_CUDA_HEADER_BOTH
  virtual Precision PlacedDistanceToOut(Vector3D<Precision> const &point, Vector3D<Precision> const &direction,
                                        const Precision stepMax = kInfinity) const override
  {
    Transformation3D const *t = this->GetTransformation();

    Precision output = kInfinity;
    Specialization::template DistanceToOut<kScalar>(
        *this->GetUnplacedVolume(), t->Transform<Specialization::transC, Specialization::rotC, Precision>(point),
        t->TransformDirection<Specialization::rotC, Precision>(direction), stepMax, output);

#ifdef VECGEOM_DISTANCE_DEBUG
    DistanceComparator::CompareDistanceToOut(this, output, point, direction, stepMax);
#endif

    return output;
  }

#ifdef VECGEOM_USOLIDS
  /*
   * WARNING: Trivial implementation for standard USolids interface
   * for DistanceToOut. The value for convex might be wrong
   */
  VECGEOM_CUDA_HEADER_BOTH
  virtual Precision DistanceToOut(Vector3D<Precision> const &point, Vector3D<Precision> const &direction,
                                  Vector3D<Precision> &normal, bool &convex, Precision step = kInfinity) const override
  {
    Precision d                  = DistanceToOut(point, direction, step);
    Vector3D<Precision> hitpoint = point + d * direction;
    PlacedShape_t::Normal(hitpoint, normal);

    // Lets the shape tell itself whether it is convex or not.
    // convex = PlacedShape_t::IsConvex;

    // Now Convexity is defined only for UnplacedVolume, not required for PlacedVolume
    convex = this->GetUnplacedVolume()->UnplacedShape_t::IsConvex();

    return d;
  }
#endif

  VECGEOM_CUDA_HEADER_BOTH
  virtual Precision SafetyToIn(Vector3D<Precision> const &point) const override
  {
    Precision output = kInfinity;
    Specialization::template SafetyToIn<kScalar>(*this->GetUnplacedVolume(), *this->GetTransformation(), point, output);

#ifdef VECGEOM_REPLACE_USOLIDS
    // apply USolids convention: convert negative values to zero
    MaskedAssign(output < kHalfTolerance, 0., &output);
#else
    // avoid distance values within tolerance
    MaskedAssign(Abs(output) < kHalfTolerance, 0., &output);
#endif

    return output;
  }

  VECGEOM_CUDA_HEADER_BOTH
  virtual Precision SafetyToOut(Vector3D<Precision> const &point) const override
  {
    Precision output = kInfinity;
    Specialization::template SafetyToOut<kScalar>(*this->GetUnplacedVolume(), point, output);

#ifdef VECGEOM_REPLACE_USOLIDS
    // apply USolids convention: convert negative values to zero
    MaskedAssign(output < kHalfTolerance, 0., &output);
#else
    // avoid distance values within tolerance
    MaskedAssign(Abs(output) < kHalfTolerance, 0., &output);
#endif

    return output;
  }

  template <class Container_t>
  void ContainsTemplate(Container_t const &points, bool *const output) const
  {
    for (int i = 0, i_max = points.size(); i < i_max; ++i) {
      Vector3D<Precision> localPoint;
      Specialization::template Contains<kScalar>(*this->GetUnplacedVolume(), *this->GetTransformation(), points[i],
                                                 localPoint, output[i]);
    }
  }

  template <class Container_t>
  void InsideTemplate(Container_t const &points, Inside_t *const output) const
  {
    for (int i = 0, i_max = points.size(); i < i_max; ++i) {
      Inside_t result = EInside::kOutside;
      Specialization::template Inside<kScalar>(*this->GetUnplacedVolume(), *this->GetTransformation(), points[i],
                                               result);
      output[i] = result;
    }
  }

  template <class Container_t>
  void DistanceToInTemplate(Container_t const &points, Container_t const &directions, Precision const *const stepMax,
                            Precision *const output) const
  {
    for (int i = 0, i_max = points.size(); i < i_max; ++i) {
      Specialization::template DistanceToIn<kScalar>(*this->GetUnplacedVolume(), *this->GetTransformation(), points[i],
                                                     directions[i], stepMax[i], output[i]);
    }
  }

  VECGEOM_FORCE_INLINE
  void DistanceToInMinimizeTemplate(SOA3D<Precision> const &points, SOA3D<Precision> const &directions, int daughterId,
                                    Precision *const currentDistance, int *const nextDaughterIdList) const
  {
    for (int i = 0, iMax = points.size(); i < iMax; ++i) {
      Precision stepMax = currentDistance[i];
      Precision result  = kInfinity;
      Specialization::template DistanceToIn<kScalar>(*this->GetUnplacedVolume(), *this->GetTransformation(), points[i],
                                                     directions[i], stepMax, result);
      if (result < currentDistance[i] && !IsInf(result)) {
        currentDistance[i]    = result;
        nextDaughterIdList[i] = daughterId;
      }
    }
  }

  template <class Container_t>
  void DistanceToOutTemplate(Container_t const &points, Container_t const &directions, Precision const *const stepMax,
                             Precision *const output) const
  {
    for (int i = 0, i_max = points.size(); i < i_max; ++i) {
      Specialization::template DistanceToOut<kScalar>(*this->GetUnplacedVolume(), points[i], directions[i], stepMax[i],
                                                      output[i]);
    }
  }

  VECGEOM_FORCE_INLINE
  void DistanceToOutTemplate(SOA3D<Precision> const &points, SOA3D<Precision> const &directions,
                             Precision const *const stepMax, Precision *const output, int *const nodeIndex) const
  {
    for (int i = 0, iMax = points.size(); i < iMax; ++i) {
      Specialization::template DistanceToOut<kScalar>(*this->GetUnplacedVolume(), points[i], directions[i], stepMax[i],
                                                      output[i]);
      if (output[i] < 0.) output[i] = vecgeom::kInfinity;
      nodeIndex[i]                  = (output[i] < stepMax[i]) ? -1 : -2;
    }
  }

  template <class Container_t>
  void SafetyToInTemplate(Container_t const &points, Precision *const output) const
  {
    for (int i = 0, i_max = points.size(); i < i_max; ++i) {
      Specialization::template SafetyToIn<kScalar>(*this->GetUnplacedVolume(), *this->GetTransformation(), points[i],
                                                   output[i]);
      // if(Abs(output[i]) < kHalfTolerance) {
      //   output[i] = 0.0;
      // }
    }
  }

  template <class Container_t>
  void SafetyToInMinimizeTemplate(Container_t const &points, Precision *const output) const
  {
    for (int i = 0, iMax = points.size(); i < iMax; ++i) {
      Precision result = 0;
      Specialization::template SafetyToIn<kScalar>(*this->GetUnplacedVolume(), *this->GetTransformation(), points[i],
                                                   result);
      output[i] = (result < output[i]) ? result : output[i];
    }
  }

  template <class Container_t>
  void SafetyToOutTemplate(Container_t const &points, Precision *const output) const
  {
    for (int i = 0, i_max = points.size(); i < i_max; ++i) {
      Specialization::template SafetyToOut<kScalar>(*this->GetUnplacedVolume(), points[i], output[i]);
      // if( Abs(output[i])<kHalfTolerance ) {
      //   output[i] = 0.0;
      // }
    }
  }

  template <class Container_t>
  void SafetyToOutMinimizeTemplate(Container_t const &points, Precision *const output) const
  {
    for (int i = 0, i_max = points.size(); i < i_max; ++i) {
      Precision result = 0;
      Specialization::template SafetyToOut<kScalar>(*this->GetUnplacedVolume(), points[i], result);
      output[i] = (result < output[i]) ? result : output[i];
    }
  }

  virtual void Contains(SOA3D<Precision> const &points, bool *const output) const override
  {
    ContainsTemplate(points, output);
  }

  virtual void Inside(SOA3D<Precision> const &points, Inside_t *const output) const override
  {
    InsideTemplate(points, output);
  }

  virtual void DistanceToIn(SOA3D<Precision> const &points, SOA3D<Precision> const &directions,
                            Precision const *const stepMax, Precision *const output) const override
  {
    DistanceToInTemplate(points, directions, stepMax, output);
  }

  virtual void DistanceToInMinimize(SOA3D<Precision> const &points, SOA3D<Precision> const &directions,
                                    int daughterindex, Precision *const output, int *const nextnodeids) const override
  {
    DistanceToInMinimizeTemplate(points, directions, daughterindex, output, nextnodeids);
  }

  virtual void DistanceToOut(SOA3D<Precision> const &points, SOA3D<Precision> const &directions,
                             Precision const *const stepMax, Precision *const output) const override
  {
    DistanceToOutTemplate(points, directions, stepMax, output);
  }

  virtual void DistanceToOut(SOA3D<Precision> const &points, SOA3D<Precision> const &directions,
                             Precision const *const stepMax, Precision *const output,
                             int *const nextNodeIndex) const override
  {
    DistanceToOutTemplate(points, directions, stepMax, output, nextNodeIndex);
  }

  virtual void SafetyToIn(SOA3D<Precision> const &points, Precision *const output) const override
  {
    SafetyToInTemplate(points, output);
  }

  // scalar fallback: dispatch a SIMD interface to a scalar kernel
  VECGEOM_FORCE_INLINE
  virtual Real_v SafetyToInVec(
      Vector3D<Real_v> const &position) const override
  {
    using vecCore::LaneAt;
    using vecCore::AssignLane;
    Real_v output = kInfinity;
    for (auto i = decltype(VECGEOM_BACKEND_PRECISION_TYPE_SIZE){0}; i < VECGEOM_BACKEND_PRECISION_TYPE_SIZE; ++i) {
      Precision tmp;
      Vector3D<Precision> pos(LaneAt(position.x(), i), LaneAt(position.y(), i), LaneAt(position.z(), i));
      Specialization::template SafetyToIn<kScalar>(*this->GetUnplacedVolume(), *this->GetTransformation(), pos, tmp);
      AssignLane(output, i, tmp);
    }
    return output;
  }

  // temporary until all shapes have migrated to VecCore backend
  VECGEOM_FORCE_INLINE
  virtual VECGEOM_BACKEND_PRECISION_TYPE SafetyToInVec(
      Vector3D<VECGEOM_BACKEND_PRECISION_TYPE> const &position) const override
  {
    using vecCore::LaneAt;
    using vecCore::AssignLane;
    VECGEOM_BACKEND_PRECISION_TYPE output = kInfinity;
    for (auto i = decltype(VECGEOM_BACKEND_PRECISION_TYPE_SIZE){0}; i < VECGEOM_BACKEND_PRECISION_TYPE_SIZE; ++i) {
      Precision tmp;
      Vector3D<Precision> pos(LaneAt(position.x(), i), LaneAt(position.y(), i), LaneAt(position.z(), i));
      Specialization::template SafetyToIn<kScalar>(*this->GetUnplacedVolume(), *this->GetTransformation(), pos, tmp);
      AssignLane(output, i, tmp);
    }
    return output;
  }

  VECGEOM_FORCE_INLINE
  virtual Real_v SafetyToOutVec(
      Vector3D<Real_v> const &position) const override
  {
    using vecCore::LaneAt;
    using vecCore::AssignLane;
    Real_v output = kInfinity;
    for (auto i = decltype(VECGEOM_BACKEND_PRECISION_TYPE_SIZE){0}; i < VECGEOM_BACKEND_PRECISION_TYPE_SIZE; ++i) {
      Precision tmp;
      Vector3D<Precision> pos(LaneAt(position.x(), i), LaneAt(position.y(), i), LaneAt(position.z(), i));
      Specialization::template SafetyToOut<kScalar>(*this->GetUnplacedVolume(), pos, tmp);
      AssignLane(output, i, tmp);
    }
    return output;
  }
  
  // temporary until all shapes have migrated to VecCore backend
  VECGEOM_FORCE_INLINE
  virtual VECGEOM_BACKEND_PRECISION_TYPE SafetyToOutVec(
      Vector3D<VECGEOM_BACKEND_PRECISION_TYPE> const &position) const override
  {
    using vecCore::LaneAt;
    using vecCore::AssignLane;
    VECGEOM_BACKEND_PRECISION_TYPE output = kInfinity;
    for (auto i = decltype(VECGEOM_BACKEND_PRECISION_TYPE_SIZE){0}; i < VECGEOM_BACKEND_PRECISION_TYPE_SIZE; ++i) {
      Precision tmp;
      Vector3D<Precision> pos(LaneAt(position.x(), i), LaneAt(position.y(), i), LaneAt(position.z(), i));
      Specialization::template SafetyToOut<kScalar>(*this->GetUnplacedVolume(), pos, tmp);
      AssignLane(output, i, tmp);
    }
    return output;
  }

  virtual Real_v DistanceToInVec(Vector3D<Real_v> const &position,
                                 Vector3D<Real_v> const &direction,
                                 const Real_v stepMax) const override
  {
    using vecCore::LaneAt;
    using vecCore::AssignLane;
    Real_v output = kInfinity;
    for (auto i = decltype(VECGEOM_BACKEND_PRECISION_TYPE_SIZE){0}; i < VECGEOM_BACKEND_PRECISION_TYPE_SIZE; ++i) {
      Precision tmp;
      Vector3D<Precision> pos(LaneAt(position.x(), i), LaneAt(position.y(), i), LaneAt(position.z(), i));
      Vector3D<Precision> dir(LaneAt(direction.x(), i), LaneAt(direction.y(), i), LaneAt(direction.z(), i));
      Specialization::template DistanceToIn<kScalar>(*this->GetUnplacedVolume(), *this->GetTransformation(), pos, dir,
                                                     LaneAt(stepMax, i), tmp);
      MaskedAssign(Abs(tmp) < kHalfTolerance, 0., &tmp);
      AssignLane(output, i, tmp);
    }
    return output;
  }
  

  // temporary until all shapes have migrated to VecCore backend
virtual VECGEOM_BACKEND_PRECISION_TYPE DistanceToInVec(Vector3D<VECGEOM_BACKEND_PRECISION_TYPE> const &position,
                                 Vector3D<VECGEOM_BACKEND_PRECISION_TYPE> const &direction,
                                 const VECGEOM_BACKEND_PRECISION_TYPE stepMax) const override
  {
    using vecCore::LaneAt;
    using vecCore::AssignLane;
    VECGEOM_BACKEND_PRECISION_TYPE output = kInfinity;
    for (auto i = decltype(VECGEOM_BACKEND_PRECISION_TYPE_SIZE){0}; i < VECGEOM_BACKEND_PRECISION_TYPE_SIZE; ++i) {
      Precision tmp;
      Vector3D<Precision> pos(LaneAt(position.x(), i), LaneAt(position.y(), i), LaneAt(position.z(), i));
      Vector3D<Precision> dir(LaneAt(direction.x(), i), LaneAt(direction.y(), i), LaneAt(direction.z(), i));
      Specialization::template DistanceToIn<kScalar>(*this->GetUnplacedVolume(), *this->GetTransformation(), pos, dir,
                                                     LaneAt(stepMax, i), tmp);
      MaskedAssign(Abs(tmp) < kHalfTolerance, 0., &tmp);
      AssignLane(output, i, tmp);
    }
    return output;
  }
  virtual Real_v DistanceToOutVec(Vector3D<Real_v> const &position,
                                  Vector3D<Real_v> const &direction,
                                  const Real_v stepMax) const override
  {
    using vecCore::LaneAt;
    using vecCore::AssignLane;
    Real_v output = kInfinity;
    for (auto i = decltype(VECGEOM_BACKEND_PRECISION_TYPE_SIZE){0}; i < VECGEOM_BACKEND_PRECISION_TYPE_SIZE; ++i) {
      Precision tmp;
      Vector3D<Precision> pos(LaneAt(position.x(), i), LaneAt(position.y(), i), LaneAt(position.z(), i));
      Vector3D<Precision> dir(LaneAt(direction.x(), i), LaneAt(direction.y(), i), LaneAt(direction.z(), i));
      Specialization::template DistanceToOut<kScalar>(*this->GetUnplacedVolume(), pos, dir, LaneAt(stepMax, i), tmp);
      MaskedAssign(Abs(tmp) < kHalfTolerance, 0., &tmp);
      AssignLane(output, i, tmp);
    }
    return output;
  }

  // temporary until all shapes have migrated to VecCore backend
  virtual VECGEOM_BACKEND_PRECISION_TYPE DistanceToOutVec(Vector3D<VECGEOM_BACKEND_PRECISION_TYPE> const &position,
                                  Vector3D<VECGEOM_BACKEND_PRECISION_TYPE> const &direction,
                                  const VECGEOM_BACKEND_PRECISION_TYPE stepMax) const override
  {
    using vecCore::LaneAt;
    using vecCore::AssignLane;
    VECGEOM_BACKEND_PRECISION_TYPE output = kInfinity;
    for (auto i = decltype(VECGEOM_BACKEND_PRECISION_TYPE_SIZE){0}; i < VECGEOM_BACKEND_PRECISION_TYPE_SIZE; ++i) {
      Precision tmp;
      Vector3D<Precision> pos(LaneAt(position.x(), i), LaneAt(position.y(), i), LaneAt(position.z(), i));
      Vector3D<Precision> dir(LaneAt(direction.x(), i), LaneAt(direction.y(), i), LaneAt(direction.z(), i));
      Specialization::template DistanceToOut<kScalar>(*this->GetUnplacedVolume(), pos, dir, LaneAt(stepMax, i), tmp);
      MaskedAssign(Abs(tmp) < kHalfTolerance, 0., &tmp);
      AssignLane(output, i, tmp);
    }
    return output;
  }
  virtual void SafetyToInMinimize(SOA3D<Precision> const &points, Precision *const safeties) const override
  {
    SafetyToInMinimizeTemplate(points, safeties);
  }

  virtual void SafetyToOut(SOA3D<Precision> const &points, Precision *const output) const override
  {
    SafetyToOutTemplate(points, output);
  }

  virtual void SafetyToOutMinimize(SOA3D<Precision> const &points, Precision *const safeties) const override
  {
    SafetyToOutMinimizeTemplate(points, safeties);
  }

}; // End class ScalarShapeImplementationHelper

} // End Impl namespace

} // End global namespace

#endif // VECGEOM_VOLUMES_SCALARSHAPEIMPLEMENTATIONHELPER_H_
