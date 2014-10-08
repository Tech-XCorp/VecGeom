/// @file PlacedTrapezoid.cpp
/// @author Guilherme Lima (lima at fnal dot gov)

#include "volumes/PlacedTrapezoid.h"
#include "volumes/Trapezoid.h"
#include "volumes/PlacedBox.h"

#if defined(VECGEOM_BENCHMARK) && defined(VECGEOM_ROOT)
#include "TGeoArb8.h"
#endif
#ifdef VECGEOM_USOLIDS
#include "UTrap.hh"
#endif

#ifdef VECGEOM_GEANT4
#include "G4Trap.hh"
#endif

namespace VECGEOM_NAMESPACE {

PlacedTrapezoid::~PlacedTrapezoid() {}

#ifdef VECGEOM_BENCHMARK

VPlacedVolume const* PlacedTrapezoid::ConvertToUnspecialized() const {
  return new SimpleTrapezoid(GetLabel().c_str(), logical_volume(), transformation());
}

#ifdef VECGEOM_ROOT
TGeoShape const* PlacedTrapezoid::ConvertToRoot() const {
  return new TGeoTrap( GetLabel().c_str(), GetDz(), GetTheta()*kRadToDeg, GetPhi()*kRadToDeg,
                       GetDy1(), GetDx1(), GetDx2(), GetTanAlpha1(),
                       GetDy2(), GetDx3(), GetDx4(), GetTanAlpha2() );
}
#endif

#ifdef VECGEOM_USOLIDS
::VUSolid const* PlacedTrapezoid::ConvertToUSolids() const {
  return new ::UTrap(GetLabel().c_str(), GetDz(), GetTheta(), GetPhi(),
                     GetDy1(), GetDx1(), GetDx2(), GetAlpha1(),
                     GetDy2(), GetDx3(), GetDx4(), GetAlpha2());
}
#endif

#ifdef VECGEOM_USOLIDS
VUSolid* PlacedTrapezoid::Clone() const {
  return new ::UTrap(GetLabel().c_str(), GetDz(), GetTheta(), GetPhi(),
                     GetDy1(), GetDx1(), GetDx2(), GetAlpha1(),
                     GetDy2(), GetDx3(), GetDx4(), GetAlpha2());
}
#endif

#ifdef VECGEOM_GEANT4
G4VSolid const* PlacedTrapezoid::ConvertToGeant4() const {
  return new G4Trap(GetLabel().c_str(), GetDz(), GetTheta(), GetPhi(),
                     GetDy1(), GetDx1(), GetDx2(), GetAlpha1(),
                     GetDy2(), GetDx3(), GetDx4(), GetAlpha2());
}
#endif


#endif // VECGEOM_BENCHMARK

} // End global namespace

namespace vecgeom {

#ifdef VECGEOM_CUDA_INTERFACE

void PlacedTrapezoid_CopyToGpu(
    LogicalVolume const *const logical_volume,
    Transformation3D const *const transformation,
    const int id, VPlacedVolume *const gpu_ptr);

VPlacedVolume* PlacedTrapezoid::CopyToGpu(
    LogicalVolume const *const logical_volume,
    Transformation3D const *const transformation,
    VPlacedVolume *const gpu_ptr) const {
  PlacedTrapezoid_CopyToGpu(logical_volume, transformation, this->id(),
                             gpu_ptr);
  CudaAssertError();
  return gpu_ptr;
}

VPlacedVolume* PlacedTrapezoid::CopyToGpu(
    LogicalVolume const *const logical_volume,
    Transformation3D const *const transformation) const {
  VPlacedVolume *const gpu_ptr = vecgeom::AllocateOnGpu<PlacedTrapezoid>();
  return this->CopyToGpu(logical_volume, transformation, gpu_ptr);
}

#endif // VECGEOM_CUDA_INTERFACE

#ifdef VECGEOM_NVCC

class LogicalVolume;
class Transformation3D;
class VPlacedVolume;

__global__
void PlacedTrapezoid_ConstructOnGpu(
    LogicalVolume const *const logical_volume,
    Transformation3D const *const transformation,
    const int id, VPlacedVolume *const gpu_ptr) {
  new(gpu_ptr) vecgeom_cuda::SimpleTrapezoid(
    reinterpret_cast<vecgeom_cuda::LogicalVolume const*>(logical_volume),
    reinterpret_cast<vecgeom_cuda::Transformation3D const*>(transformation),
    NULL,
    id
  );
}

void PlacedTrapezoid_CopyToGpu(
    LogicalVolume const *const logical_volume,
    Transformation3D const *const transformation,
    const int id, VPlacedVolume *const gpu_ptr) {
  PlacedTrapezoid_ConstructOnGpu<<<1, 1>>>(logical_volume, transformation,
                                            id, gpu_ptr);
}

#endif // VECGEOM_NVCC

/*
void PlacedTrapezoid::ComputeBoundingBox() {
  Vector3D<Precision> aMin, aMax;
  GetUnplacedVolume()->Extent(aMin, aMax) ;

  // try a box with no rotation
  Vector3D<Precision> bbdims1 = 0.5*(aMax-aMin);
  Vector3D<Precision> center1 = 0.5*(aMax+aMin);
  UnplacedBox *box1 = new UnplacedBox(bbdims1);
  Precision vol1 = box1->volume();

  // try a box with a rotation by theta,phi
  Transformation3D* matrix2 =
    new Transformation3D(center1.x(), center1.y(), center1.z(),
                         this->GetTheta(), this->GetPhi(), 0);
  Vector3D<Precision> newMin, newMax;
  matrix2->Transform(aMin, newMin);
  matrix2->Transform(aMax, newMax);
  UnplacedBox *box2 = new UnplacedBox(0.5*(newMax-newMin));
  Precision vol2 = box2->volume();

  if(vol2>0.5*vol1) {
    // use box1
    bounding_box_ =
      new PlacedBox(new LogicalVolume(box1),
                    new Transformation3D(center1.x(), center1.y(), center1.z()),
                    SimpleBox(box1));
    delete box2, matrix2;
  }
  else {
    // use box2
    bounding_box_ = new PlacedBox(new LogicalVolume(box2), matrix2, 0);
    delete box1;
  }
*/
} // End namespace vecgeom
