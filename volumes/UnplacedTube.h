/// @file UnplacedTube.h
/// @author Georgios Bitzes (georgios.bitzes@cern.ch)

#ifndef VECGEOM_VOLUMES_UNPLACEDTUBE_H_
#define VECGEOM_VOLUMES_UNPLACEDTUBE_H_

#include "base/global.h"
#include "base/AlignedBase.h"
#include "volumes/unplaced_volume.h"

namespace VECGEOM_NAMESPACE {

class UnplacedTube : public VUnplacedVolume, AlignedBase {

private:
  // tube defining parameters
  Precision fRmin, fRmax, fZ, fSphi, fDphi;

  // cached values
  Precision fRmin2, fRmax2, fAlongPhi1x, fAlongPhi1y, fAlongPhi2x, fAlongPhi2y;
  Precision fTolIrmin2, fTolOrmin2, fTolIrmax2, fTolOrmax2, fTolIz, fTolOz;

  VECGEOM_CUDA_HEADER_BOTH
  static void GetAlongVectorToPhiSector(Precision phi, Precision &x, Precision &y) {
    x = std::cos(phi);
    y = std::sin(phi);
  }

  VECGEOM_CUDA_HEADER_BOTH
  void calculateCached() {
    fTolIz = fZ - kTolerance;
    fTolOz = fZ + kTolerance;

    fRmin2 = fRmin * fRmin;
    fRmax2 = fRmax * fRmax;

    fTolOrmin2 = (fRmin - kTolerance)*(fRmin - kTolerance);
    fTolIrmin2 = (fRmin + kTolerance)*(fRmin + kTolerance);
    
    fTolOrmax2 = (fRmax + kTolerance)*(fRmax + kTolerance);
    fTolIrmax2 = (fRmax - kTolerance)*(fRmax - kTolerance);

    GetAlongVectorToPhiSector(fSphi, fAlongPhi1x, fAlongPhi1y);
    GetAlongVectorToPhiSector(fSphi + fDphi, fAlongPhi2x, fAlongPhi2y);
  }

public:
  UnplacedTube(const Precision rmin, const Precision rmax, const Precision z,
               const Precision sphi, const Precision dphi) : fRmin(rmin), fRmax(rmax),
               fZ(z), fSphi(sphi), fDphi(dphi) {
    calculateCached();  
  }

  VECGEOM_CUDA_HEADER_BOTH
  UnplacedTube(UnplacedTube const &other) : fRmin(other.fRmin), fRmax(other.fRmax), fZ(other.fZ), fSphi(other.fSphi), fDphi(other.fDphi)  {
    calculateCached();
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision rmin() const { return fRmin; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision rmax() const { return fRmax; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision z() const { return fZ; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision sphi() const { return fSphi; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision dphi() const { return fDphi; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision rmin2() const { return fRmin2; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision rmax2() const { return fRmax2; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision alongPhi1x() const { return fAlongPhi1x; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision alongPhi1y() const { return fAlongPhi1y; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision alongPhi2x() const { return fAlongPhi2x; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision alongPhi2y() const { return fAlongPhi2y; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision tolIz() const { return fTolIz; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision tolOz() const { return fTolOz; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision tolOrmin2() const { return fTolOrmin2; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision tolIrmin2() const { return fTolIrmin2; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision tolOrmax2() const { return fTolOrmax2; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision tolIrmax2() const { return fTolIrmax2; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision volume() const {
    return fZ * (fRmax2 - fRmin2) * fDphi;
  }

  virtual int memory_size() const { return sizeof(*this); }

  VECGEOM_CUDA_HEADER_BOTH
  virtual void Print() const;

  template <TranslationCode transCodeT, RotationCode rotCodeT>
  VECGEOM_CUDA_HEADER_DEVICE
  static VPlacedVolume* Create(LogicalVolume const *const logical_volume,
                               Transformation3D const *const transformation,
#ifdef VECGEOM_NVCC
                               const int id,
#endif
                               VPlacedVolume *const placement = NULL);

#ifdef VECGEOM_CUDA_INTERFACE
  virtual VUnplacedVolume* CopyToGpu() const;
  virtual VUnplacedVolume* CopyToGpu(VUnplacedVolume *const gpu_ptr) const;
#endif

private:

  virtual void Print(std::ostream &os) const;

  VECGEOM_CUDA_HEADER_DEVICE
  virtual VPlacedVolume* SpecializedVolume(
      LogicalVolume const *const volume,
      Transformation3D const *const transformation,
      const TranslationCode trans_code, const RotationCode rot_code,
#ifdef VECGEOM_NVCC
      const int id,
#endif
      VPlacedVolume *const placement = NULL) const;

};

} // end global namespace

#endif // VECGEOM_VOLUMES_UNPLACEDTUBE_H_




