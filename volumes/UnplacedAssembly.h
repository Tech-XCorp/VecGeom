// LICENSING INFORMATION TBD

#ifndef VECGEOM_UNPLACEDASSEMBLY_H
#define VECGEOM_UNPLACEDASSEMBLY_H

#include "base/Global.h"
#include "base/Vector.h"
#include "navigation/VNavigator.h"
#include "navigation/NavigationState.h"
#include "base/Vector3D.h"
#include "volumes/UnplacedVolume.h"
#include "volumes/PlacedVolume.h"
#include "volumes/kernel/BoxImplementation.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class UnplacedAssembly;);
VECGEOM_DEVICE_DECLARE_CONV(class, UnplacedAssembly);

inline namespace VECGEOM_IMPL_NAMESPACE {

class NavigationState;

// An assembly volume offering navigation interfaces (Contains/Distances/...) in a loose/logical group of volumes
// An UnplacedAssembly is always strongly coupled to a logical volume because the later naturally manages the actual
// group of volumes that define the assembly

// The following construct marks a logical volume as an assembly:

// UnplacedAssembly *ass = new UnplacedAssembly
// LogicalVolume *lv = new LogicalVolume("assembly", ass); // this will implicitely couple ass to lv
// lv->PlacedDaughter(...)

class UnplacedAssembly : public VUnplacedVolume, public AlignedBase {

private:
  // back-reference to the logical volume
  // this get automacially set upon instantiation of a logical volume with an UnplacedAssembly
  LogicalVolume *fLogicalVolume;

  // caching the extent (bounding box)
  // these members are automatically updated whenever a new volume is added to the assembly
  Vector3D<double> fLowerCorner;
  Vector3D<double> fUpperCorner;

  void SetLogicalVolume(LogicalVolume *lv) { fLogicalVolume = lv; }
  void UpdateExtent() { UnplacedAssembly::Extent(fLowerCorner, fUpperCorner); }
  friend class LogicalVolume;

public:
  VECCORE_ATT_HOST_DEVICE
  UnplacedAssembly(); // the constructor

  VECCORE_ATT_HOST_DEVICE
  virtual ~UnplacedAssembly();

  LogicalVolume const *GetLogicalVolume() const { return fLogicalVolume; }

  // add content
  void AddVolume(VPlacedVolume const *);

  // get number of volumes
  size_t GetNVolumes() const { return fLogicalVolume->GetDaughters().size(); }

  // the extent function
  VECCORE_ATT_HOST_DEVICE
  void Extent(Vector3D<Precision> &, Vector3D<Precision> &) const override;

  // Getter to cached bounding box
  VECCORE_ATT_HOST_DEVICE
  Vector3D<Precision> GetLowerCorner() const { return fLowerCorner; }

  VECCORE_ATT_HOST_DEVICE
  Vector3D<Precision> GetUpperCorner() const { return fUpperCorner; }

  // the ordinary assembly function
  VECCORE_ATT_HOST_DEVICE
  bool Contains(Vector3D<Precision> const &point) const override
  {
    assert(fLogicalVolume);
    // check bound box first
    bool inBoundingBox;
    ABBoxImplementation::ABBoxContainsKernel(fLowerCorner, fUpperCorner, point, inBoundingBox);
    if (!inBoundingBox) return false;

    Vector3D<Precision> daughterlocalpoint;
    VPlacedVolume const *nextv;
    return fLogicalVolume->GetLevelLocator()->LevelLocate(fLogicalVolume, point, nextv, daughterlocalpoint);
  }

  // an extended contains function needed for navigation
  // if this function returns true it modifies the navigation state to point to the first non-assembly volume
  // the point is contained in
  // this function is not part of the generic UnplacedVolume interface but we could consider doing so
  VECCORE_ATT_HOST_DEVICE
  bool Contains(Vector3D<Precision> const &point, Vector3D<Precision> &daughterlocalpoint, NavigationState &state) const
  {
    assert(fLogicalVolume);
    // check bound box first
    bool inBoundingBox;
    ABBoxImplementation::ABBoxContainsKernel(fLowerCorner, fUpperCorner, point, inBoundingBox);
    if (!inBoundingBox) return false;

    return fLogicalVolume->GetLevelLocator()->LevelLocate(fLogicalVolume, point, state, daughterlocalpoint);
  }

  using VUnplacedVolume::DistanceToOut;
  // DistanceToOut does not make sense -- throw exeption
  VECCORE_ATT_HOST_DEVICE
  Precision DistanceToOut(Vector3D<Precision> const & /*p*/, Vector3D<Precision> const & /*d*/,
                          Precision /*step_max*/ = kInfLength) const override
  {
#ifndef VECCORE_CUDA
    throw std::runtime_error("Forbidden DistanceToOut in Assembly called");
#endif
    return -1.;
  }

  // DistanceToOut does not make sense -- throw exeption
  VECCORE_ATT_HOST_DEVICE
  Real_v DistanceToOutVec(Vector3D<Real_v> const & /*p*/, Vector3D<Real_v> const & /*d*/,
                          Real_v const & /*step_max*/) const override
  {
#ifndef VECCORE_CUDA
    throw std::runtime_error("Forbidden DistanceToOut in Assembly called");
#endif
    return Real_v(-1.);
  }

  using VUnplacedVolume::SafetyToOut;
  VECCORE_ATT_HOST_DEVICE
  Precision SafetyToOut(Vector3D<Precision> const &) const override
  {
#ifndef VECCORE_CUDA
    throw std::runtime_error("Forbidden SafetyToOut in Assembly called");
#endif
    return -1.;
  }

  // an explicit SIMD interface
  VECCORE_ATT_HOST_DEVICE
  Real_v SafetyToOutVec(Vector3D<Real_v> const &) const override
  {
#ifndef VECCORE_CUDA
    throw std::runtime_error("Forbidden SafetyToOut in Assembly called");
#endif
    return Real_v(-1.);
  }

  // ---------------- SafetyToIn functions -------------------------------------------------------
  VECCORE_ATT_HOST_DEVICE
  virtual Precision SafetyToIn(Vector3D<Precision> const &p) const override
  {
    return fLogicalVolume->GetSafetyEstimator()->ComputeSafetyToDaughtersForLocalPoint(p, fLogicalVolume);
  }

  // explicit SIMD interface
  VECCORE_ATT_HOST_DEVICE
  virtual Real_v SafetyToInVec(Vector3D<Real_v> const &) const override
  {
#ifndef VECCORE_CUDA
    throw std::runtime_error("SafetyToInVec in Assembly not yet implemented");
#endif
    return Real_v(-1.);
  }

  // ---------------- DistanceToIn functions -----------------------------------------------------

  VECCORE_ATT_HOST_DEVICE
  virtual Precision DistanceToIn(Vector3D<Precision> const &p, Vector3D<Precision> const &d,
                                 const Precision /*step_max*/ = kInfLength) const override
  {
    if (!BoxImplementation::Intersect(&fLowerCorner, p, d, 0, kInfLength)) return kInfLength;

    Precision step(kInfLength);
    VPlacedVolume const *pv;
    fLogicalVolume->GetNavigator()->CheckDaughterIntersections(fLogicalVolume, p, d, nullptr, nullptr, step, pv);
    return step;
  }

  VECCORE_ATT_HOST_DEVICE
  virtual Real_v DistanceToInVec(Vector3D<Real_v> const & /*p*/, Vector3D<Real_v> const & /*d*/,
                                 const Real_v & /*step_max*/ = Real_v(kInfLength)) const override
  {
#ifndef VECCORE_CUDA
    throw std::runtime_error("DistanceToInVec in Assembly not yet implemented");
#endif
    return Real_v(-1.);
  }

  Vector3D<Precision> GetPointOnSurface() const override;
  Precision Capacity() const;
  Precision SurfaceArea() const;

  // some dummy impl for virtual functions
  VECCORE_ATT_HOST_DEVICE
  virtual void Print() const override;
  virtual void Print(std::ostream &os) const override;
  virtual int MemorySize() const override { return sizeof(*this); }

#ifdef VECGEOM_CUDA_INTERFACE
  virtual size_t DeviceSizeOf() const override { return DevicePtr<cuda::UnplacedAssembly>::SizeOf(); }
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu() const override;
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const gpu_ptr) const override;
#endif

#ifndef VECCORE_CUDA
  virtual VPlacedVolume *SpecializedVolume(LogicalVolume const *const volume,
                                           Transformation3D const *const transformation,
                                           const TranslationCode trans_code, const RotationCode rot_code,
                                           VPlacedVolume *const placement = NULL) const override;
#else
  __device__ VPlacedVolume *SpecializedVolume(LogicalVolume const *const volume,
                                              Transformation3D const *const transformation,
                                              const TranslationCode trans_code, const RotationCode rot_code,
                                              const int id, VPlacedVolume *const placement) const override;
#endif
};
}
}

#endif // VECGEOM_UNPLACEDASSEMBLY_H
