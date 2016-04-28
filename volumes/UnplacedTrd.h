/// @file UnplacedTrd.h
/// @author Georgios Bitzes (georgios.bitzes@cern.ch)

#ifndef VECGEOM_VOLUMES_UNPLACEDTRD_H_
#define VECGEOM_VOLUMES_UNPLACEDTRD_H_

#include "base/Global.h"
#include "base/AlignedBase.h"
#include "volumes/UnplacedVolume.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE( class UnplacedTrd; )
VECGEOM_DEVICE_DECLARE_CONV( UnplacedTrd )

inline namespace VECGEOM_IMPL_NAMESPACE {

class UnplacedTrd : public VUnplacedVolume, public AlignedBase {
private:
  // trd defining parameters
  Precision fDX1;   //Half-length along x at the surface positioned at -dz
  Precision fDX2;   //Half-length along x at the surface positioned at +dz
  Precision fDY1;   //Half-length along y at the surface positioned at -dz
  Precision fDY2;   //Half-length along y at the surface positioned at +dz
  Precision fDZ;    //Half-length along z axis

  // cached values
  Precision fX2minusX1;
  Precision fY2minusY1;
  Precision fHalfX1plusX2;
  Precision fHalfY1plusY2;
  Precision fCalfX, fCalfY;
  Precision fSecxz, fSecyz;
  Precision fToleranceX;    // Corrected tolerance for Inside checks on X
  Precision fToleranceY;    // Corrected tolerance for Inside checks on Y

  Precision fFx, fFy;

  VECGEOM_CUDA_HEADER_BOTH
  void calculateCached() {
    fX2minusX1 = fDX2 - fDX1;
    fY2minusY1 = fDY2 - fDY1;
    fHalfX1plusX2 = 0.5 * (fDX1 + fDX2);
    fHalfY1plusY2 = 0.5 * (fDY1 + fDY2);

    fFx = 0.5*(fDX1 - fDX2)/fDZ;
    fFy = 0.5*(fDY1 - fDY2)/fDZ;
    fSecxz = sqrt(1+fFx*fFx);
    fSecyz = sqrt(1+fFy*fFy);

    fCalfX = 1./Sqrt(1.0+fFx*fFx);
    fCalfY = 1./Sqrt(1.0+fFy*fFy);
    fToleranceX = kTolerance * Sqrt(fX2minusX1*fX2minusX1 + 4*fDZ*fDZ);
    fToleranceY = kTolerance * Sqrt(fX2minusX1*fX2minusX1 + 4*fDZ*fDZ);
  }

public:
  // special case Trd1 when dY1 == dY2
  VECGEOM_CUDA_HEADER_BOTH
  UnplacedTrd(const Precision x1, const Precision x2, const Precision y1, const Precision z)
    : fDX1(x1)
    , fDX2(x2)
    , fDY1(y1)
    , fDY2(y1)
    , fDZ(z)
    , fX2minusX1(0)
    , fY2minusY1(0)
    , fHalfX1plusX2(0)
    , fHalfY1plusY2(0)
    , fCalfX(0)
    , fCalfY(0)
    , fFx(0)
    , fFy(0)
  {
    calculateCached();
    fGlobalConvexity = true;
  }

  // general case
  VECGEOM_CUDA_HEADER_BOTH
  UnplacedTrd(const Precision x1, const Precision x2, const Precision y1, const Precision y2, const Precision z)
    : fDX1(x1)
    , fDX2(x2)
    , fDY1(y1)
    , fDY2(y2)
    , fDZ(z)
    , fX2minusX1(0)
    , fY2minusY1(0)
    , fHalfX1plusX2(0)
    , fHalfY1plusY2(0)
    , fCalfX(0)
    , fCalfY(0)
    , fFx(0)
    , fFy(0)
  {
    calculateCached();
    fGlobalConvexity = true;
  }

  VECGEOM_CUDA_HEADER_BOTH
  void SetAllParameters(Precision x1, Precision x2, Precision y1, Precision y2, Precision z) {
    fDX1 = x1;
    fDX2 = x2;
    fDY1 = y1;
    fDY2 = y2;
    fDZ  = z;
    calculateCached();
  }

  VECGEOM_CUDA_HEADER_BOTH
  void SetXHalfLength1(Precision arg) { fDX1 = arg; calculateCached(); }
  VECGEOM_CUDA_HEADER_BOTH
  void SetXHalfLength2(Precision arg) { fDX2 = arg; calculateCached(); }
  VECGEOM_CUDA_HEADER_BOTH
  void SetYHalfLength1(Precision arg) { fDY1 = arg; calculateCached(); }
  VECGEOM_CUDA_HEADER_BOTH
  void SetYHalfLength2(Precision arg) { fDY2 = arg; calculateCached(); }
  VECGEOM_CUDA_HEADER_BOTH
  void SetZHalfLength(Precision arg)  { fDZ  = arg; calculateCached(); }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision dx1() const { return fDX1; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision dx2() const { return fDX2; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision dy1() const { return fDY1; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision dy2() const { return fDY2; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision dz() const { return fDZ; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision x2minusx1() const { return fX2minusX1; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision y2minusy1() const { return fY2minusY1; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision halfx1plusx2() const { return fHalfX1plusX2; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision halfy1plusy2() const { return fHalfY1plusY2; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision fx() const { return fFx; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision fy() const { return fFy; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision calfx() const { return fCalfX; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision calfy() const { return fCalfY; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision ToleranceX() const { return fToleranceX; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision ToleranceY() const { return fToleranceY; }

  virtual int memory_size() const final { return sizeof(*this); }

  VECGEOM_CUDA_HEADER_BOTH
  void Extent(Vector3D<Precision> & aMin, Vector3D<Precision> & aMax ) const {
      aMin = Vector3D<Precision>(-Max(fDX1, fDX2), -Max(fDY1, fDY2), -fDZ);
      aMax = Vector3D<Precision>( Max(fDX1, fDX2), Max(fDY1, fDY2), fDZ);
  }


#ifndef VECGEOM_NVCC
  // Computes capacity of the shape in [length^3]
  Precision Capacity() const;

  Precision SurfaceArea() const;

  Precision GetPlusXArea() const { //  Area in +x direction 
      return 2*fDZ * (fDY1 + fDY2) * fSecyz;
  }

  Precision GetMinusXArea() const {  // Area in -x direction
      return GetPlusXArea();
  }

  Precision GetPlusYArea() const {    // Area in +y direction
      return 2*fDZ * (fDX1 + fDX2) * fSecxz;
  }

  Precision GetMinusYArea() const {  // Area in -y direction
      return GetPlusYArea();
  }  

  Precision GetPlusZArea() const {   // Area in +Z
      return 4 * fDX2 * fDY2;
  }
  
  Precision GetMinusZArea() const {  // Area in -Z
      return 4 * fDX1 * fDY1;
  }
 
  int ChooseSurface() const;

  Vector3D<Precision> GetPointOnSurface() const;

  bool Normal(Vector3D<Precision> const& point, Vector3D<Precision>& normal) const;

#endif

  VECGEOM_CUDA_HEADER_BOTH
  virtual void Print() const final;

  std::string GetEntityType() const { return "Trd";}

  template <TranslationCode transCodeT, RotationCode rotCodeT>
  VECGEOM_CUDA_HEADER_DEVICE
  static VPlacedVolume* Create(LogicalVolume const *const logical_volume,
                               Transformation3D const *const transformation,
#ifdef VECGEOM_NVCC
                               const int id,
#endif
                               VPlacedVolume *const placement = NULL);

#ifdef VECGEOM_CUDA_INTERFACE
  virtual size_t DeviceSizeOf() const { return DevicePtr<cuda::UnplacedTrd>::SizeOf(); }
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu() const;
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const gpu_ptr) const;
#endif

#if defined(VECGEOM_USOLIDS)
  std::ostream& StreamInfo(std::ostream &os) const;
#endif

private:

  virtual void Print(std::ostream &os) const final;

  VECGEOM_CUDA_HEADER_DEVICE
  virtual VPlacedVolume* SpecializedVolume(
      LogicalVolume const *const volume,
      Transformation3D const *const transformation,
      const TranslationCode trans_code, const RotationCode rot_code,
#ifdef VECGEOM_NVCC
      const int id,
#endif
      VPlacedVolume *const placement = NULL) const final;
};

} } // end global namespace

#endif // VECGEOM_VOLUMES_UNPLACEDTRD_H_
