/// \file PlacedBox.h
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#ifndef VECGEOM_VOLUMES_PLACEDBOX_H_
#define VECGEOM_VOLUMES_PLACEDBOX_H_

#include "base/Global.h"
#include "backend/Backend.h"
 
#include "volumes/PlacedVolume.h"
#include "volumes/UnplacedVolume.h"
#include "volumes/kernel/BoxImplementation.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE( class PlacedBox; )
VECGEOM_DEVICE_DECLARE_CONV( PlacedBox )

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedBox : public VPlacedVolume {

public:

#ifndef VECGEOM_NVCC

  PlacedBox(char const *const label,
            LogicalVolume const *const logicalVolume,
            Transformation3D const *const transformation,
            PlacedBox const *const boundingBox)
      : VPlacedVolume(label, logicalVolume, transformation, boundingBox) {}

  PlacedBox(LogicalVolume const *const logicalVolume,
            Transformation3D const *const transformation,
            PlacedBox const *const boundingBox)
      : PlacedBox("", logicalVolume, transformation, boundingBox) {}

#else

  __device__
  PlacedBox(LogicalVolume const *const logicalVolume,
            Transformation3D const *const transformation,
            PlacedBox const *const boundingBox,
            const int id)
      : VPlacedVolume(logicalVolume, transformation, boundingBox, id) {}

#endif
  VECGEOM_CUDA_HEADER_BOTH
  virtual ~PlacedBox() {}

  // Accessors

  VECGEOM_CUDA_HEADER_BOTH
  UnplacedBox const* GetUnplacedVolume() const {
    return static_cast<UnplacedBox const *>(
        GetLogicalVolume()->GetUnplacedVolume());
  }


  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Vector3D<Precision> const& dimensions() const {
    return GetUnplacedVolume()->dimensions();
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision x() const { return GetUnplacedVolume()->x(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision y() const { return GetUnplacedVolume()->y(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision z() const { return GetUnplacedVolume()->z(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision GetXHalfLength() const { return x(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  void SetXHalfLength(double dx) { const_cast<UnplacedBox*>(GetUnplacedVolume())->SetX(dx); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision GetYHalfLength() const { return y(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  void SetYHalfLength(double dy) { const_cast<UnplacedBox*>(GetUnplacedVolume())->SetY(dy); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision GetZHalfLength() const { return z(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  void SetZHalfLength(double dz) { const_cast<UnplacedBox*>(GetUnplacedVolume())->SetZ(dz); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  void Set(Precision xx, Precision yy, Precision zz,
           vecgeom::Transformation3D const*const transformation = &Transformation3D::kIdentity) {
    if(logical_volume_) delete logical_volume_;
    logical_volume_ = new LogicalVolume(new UnplacedBox(xx,yy,zz));
#ifdef VECGEOM_INPLACE_TRANSFORMATIONS
    fTransformation = *transformation;
#else
    if(fTransformation) delete fTransformation;
    fTransformation = new Transformation3D(*transformation);
#endif
  }

#if !defined(VECGEOM_NVCC)
  virtual Precision Capacity() override {
      return GetUnplacedVolume()->volume();
  }

  virtual
  void Extent(Vector3D<Precision> & aMin, Vector3D<Precision> & aMax) const override
  {
    GetUnplacedVolume()->Extent(aMin, aMax);
  }

  virtual
  bool Normal(Vector3D<Precision> const & point, Vector3D<Precision> & normal ) const override
  {
      bool valid;
      BoxImplementation<translation::kIdentity, rotation::kIdentity>::NormalKernel<kScalar>(
              *GetUnplacedVolume(),
              point,
              normal, valid);
      return valid;
  }

  virtual
  Vector3D<Precision> GetPointOnSurface() const override {
    return GetUnplacedVolume()->GetPointOnSurface();
  }

  virtual double SurfaceArea() override {
     return GetUnplacedVolume()->SurfaceArea();
  }

#if defined(VECGEOM_USOLIDS)
  virtual std::string GetEntityType() const override { return GetUnplacedVolume()->GetEntityType() ;}
#endif
#endif

  VECGEOM_CUDA_HEADER_BOTH
  virtual void PrintType() const override;

  // CUDA specific

  virtual int memory_size() const override { return sizeof(*this); }

  // Comparison specific

#if defined(VECGEOM_USOLIDS)
//  VECGEOM_CUDA_HEADER_BOTH
  std::ostream& StreamInfo(std::ostream &os) const override {
    return GetUnplacedVolume()->StreamInfo(os);
  }
#endif

#ifndef VECGEOM_NVCC
  virtual VPlacedVolume const* ConvertToUnspecialized() const override;
#ifdef VECGEOM_ROOT
  virtual TGeoShape const* ConvertToRoot() const override;
#endif
#if defined(VECGEOM_USOLIDS) && !defined(VECGEOM_REPLACE_USOLIDS)
  virtual ::VUSolid const* ConvertToUSolids() const override;
#endif
#ifdef VECGEOM_GEANT4
  virtual G4VSolid const* ConvertToGeant4() const override;
#endif
#endif // VECGEOM_NVCC

};

} } // End global namespace

#endif // VECGEOM_VOLUMES_PLACEDBOX_H_
