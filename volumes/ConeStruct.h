/*
 * ConeStruct.h
 *
 *  Created on: May 09, 2017
 *      Author: Raman Sehgal
 */
#ifndef VECGEOM_CONESTRUCT_H_
#define VECGEOM_CONESTRUCT_H_

#include "base/Global.h"
#include "volumes/Wedge_Evolution.h"

#define kConeTolerance 1e-7

namespace vecgeom {

inline namespace VECGEOM_IMPL_NAMESPACE {

// a plain and lightweight struct to encapsulate data members of a Cone
template <typename T = double>
struct ConeStruct {
  // Cone defining parameters
  T fRmin1;
  T fRmax1;
  T fRmin2;
  T fRmax2;
  T fDz;
  T fSPhi;
  T fDPhi;

  /* These new data members are introduced to store the original paramters of
   * Cone, which may change in the case where rmin is equal to rmax.
   * These are basically required by the Extent functions to do more accurate
   * bounding box calculations.
   */
  T _frmin1;
  T _frmin2;
  T _frmax1;
  T _frmax2;

  evolution::Wedge fPhiWedge;

  // vectors characterizing the normals of phi planes
  // makes task to detect phi sektors very efficient
  Vector3D<Precision> fNormalPhi1;
  Vector3D<Precision> fNormalPhi2;
  Precision fAlongPhi1x;
  Precision fAlongPhi1y;
  Precision fAlongPhi2x;
  Precision fAlongPhi2y;

  // Some Cached value, try to reduce them
  // Some precomputed values to avoid divisions etc
  Precision fInnerSlope; // "gradient" of inner surface in z direction
  Precision fOuterSlope; // "gradient" of outer surface in z direction
  Precision fInnerOffset;
  Precision fOuterOffset;
  Precision fInnerTolerance; // tolerance on radial direction for inner surface
  Precision fOuterTolerance; // tolerance on radial direction for outer surface
  // Values to be cached
  Precision fSqRmin1, fSqRmin2;
  Precision fSqRmax1, fSqRmax2;
  Precision fTolIz, fTolOz;
  Precision fInnerConeApex;
  Precision fTanInnerApexAngle;
  Precision fOuterConeApex;
  Precision fTanOuterApexAngle;

  Precision fSecRMin;
  Precision fSecRMax;
  Precision fInvSecRMin;
  Precision fInvSecRMax;
  Precision fTanRMin;
  Precision fTanRMax;
  Precision fZNormInner;
  Precision fZNormOuter;
  Precision fConeTolerance;

#ifndef VECCORE_CUDA
  Precision Capacity() const
  {
    return (fDz * fDPhi / 3.) *
           (fRmax1 * fRmax1 + fRmax2 * fRmax2 + fRmax1 * fRmax2 - fRmin1 * fRmin1 - fRmin2 * fRmin2 - fRmin1 * fRmin2);
  }
#endif

  VECCORE_ATT_HOST_DEVICE
  void CalculateCached()
  {
    if (fRmin1 == fRmax1) {
      fRmax1 += kConeTolerance;
    }

    if (fRmin2 == fRmax2) {
      fRmax2 += kConeTolerance;
    }

    fSqRmin1       = fRmin1 * fRmin1;
    fSqRmax1       = fRmax1 * fRmax1;
    fSqRmin2       = fRmin2 * fRmin2;
    fSqRmax2       = fRmax2 * fRmax2;
    fConeTolerance = 1e-7;

    fTanRMin    = (fRmin2 - fRmin1) * 0.5 / fDz;
    fSecRMin    = std::sqrt(1.0 + fTanRMin * fTanRMin);
    fInvSecRMin = 1. / NonZero(fSecRMin);
    fTanRMax    = (fRmax2 - fRmax1) * 0.5 / fDz;

    fSecRMax    = std::sqrt(1.0 + fTanRMax * fTanRMax);
    fInvSecRMax = 1. / NonZero(fSecRMax);

    // check this very carefully
    fInnerSlope     = -(fRmin1 - fRmin2) / (2. * fDz);
    fOuterSlope     = -(fRmax1 - fRmax2) / (2. * fDz);
    fInnerOffset    = fRmin2 - fInnerSlope * fDz;
    fOuterOffset    = fRmax2 - fOuterSlope * fDz;
    fInnerTolerance = kConeTolerance * fSecRMin;
    fOuterTolerance = kConeTolerance * fSecRMax;

    if (fRmin2 > fRmin1) {
      fInnerConeApex     = 2 * fDz * fRmin1 / (fRmin2 - fRmin1);
      fTanInnerApexAngle = fRmin2 / (2 * fDz + fInnerConeApex);
    } else { // Should we add a check if(fRmin1 > fRmin2)
      fInnerConeApex     = 2 * fDz * fRmin2 / NonZero(fRmin1 - fRmin2);
      fTanInnerApexAngle = fRmin1 / (2 * fDz + fInnerConeApex);
    }

    if (fRmin1 == 0. || fRmin2 == 0.) fInnerConeApex = 0.;

    if (fRmin1 == 0.) fTanInnerApexAngle = fRmin2 / (2 * fDz);
    if (fRmin2 == 0.) fTanInnerApexAngle = fRmin1 / (2 * fDz);

    if (fRmax2 > fRmax1) {
      fOuterConeApex     = 2 * fDz * fRmax1 / (fRmax2 - fRmax1);
      fTanOuterApexAngle = fRmax2 / (2 * fDz + fOuterConeApex);
    } else { // Should we add a check if(fRmax1 > fRmax2)
      fOuterConeApex     = 2 * fDz * fRmax2 / NonZero(fRmax1 - fRmax2);
      fTanOuterApexAngle = fRmax1 / (2 * fDz + fOuterConeApex);
    }

    if (fRmax1 == 0. || fRmax2 == 0.) fOuterConeApex = 0.;

    if (fRmax1 == 0.) fTanOuterApexAngle = fRmax2 / (2 * fDz);
    if (fRmax2 == 0.) fTanOuterApexAngle = fRmax1 / (2 * fDz);

    fZNormInner = fTanRMin / NonZero(fSecRMin);
    fZNormOuter = -fTanRMax / NonZero(fSecRMax);

    fTolIz = fDz - kHalfTolerance;
    fTolOz = fDz + kHalfTolerance;

    // DetectConvexity();
  }

  VECCORE_ATT_HOST_DEVICE
  void Print() const
  {
    printf("UnplacedCone {rmin1 %.2f, rmax1 %.2f, rmin2 %.2f, "
           "rmax2 %.2f, dz %.2f, phistart %.2f, deltaphi %.2f}",
           fRmin1, fRmax1, fRmin2, fRmax2, fDz, fSPhi, fDPhi);
  }

  void Print(std::ostream &os) const { os << "UnplacedCone; please implement Print to outstream\n"; }

  VECCORE_ATT_HOST_DEVICE
  bool IsFullPhi() const { return fDPhi == kTwoPi; }

  bool Normal(Vector3D<Precision> const &p, Vector3D<Precision> &norm) const
  {
    int noSurfaces = 0;
    Precision rho, pPhi;
    Precision distZ, distRMin, distRMax;
    Precision distSPhi = kInfLength, distEPhi = kInfLength;
    Precision pRMin, widRMin;
    Precision pRMax, widRMax;

    // const double kHalfTolerance = 0.5 * kTolerance;

    Vector3D<Precision> sumnorm(0., 0., 0.), nZ = Vector3D<Precision>(0., 0., 1.);
    Vector3D<Precision> nR, nr(0., 0., 0.), nPs, nPe;
    norm = sumnorm;

    // do not use an extra fabs here -- negative/positive distZ tells us when point is outside or inside
    distZ = std::fabs(p.z()) - fDz;
    rho   = std::sqrt(p.x() * p.x() + p.y() * p.y());

    pRMin    = rho - p.z() * fTanRMin;
    widRMin  = fRmin2 - fDz * fTanRMin;
    distRMin = (pRMin - widRMin) / fSecRMin;

    pRMax    = rho - p.z() * fTanRMax;
    widRMax  = fRmax2 - fDz * fTanRMax;
    distRMax = (pRMax - widRMax) / fSecRMax;

    bool inside = distZ < kTolerance && distRMax < kTolerance;
    if (fRmin1 || fRmin2) inside &= distRMin > -kTolerance;

    distZ    = std::fabs(distZ);
    distRMax = std::fabs(distRMax);
    distRMin = std::fabs(distRMin);

    // keep track of nearest normal, needed in case point is not on a surface
    double distNearest              = distZ;
    Vector3D<Precision> normNearest = nZ;
    if (p.z() < 0.) normNearest.Set(0, 0, -1.);

    if (!IsFullPhi()) {
      if (rho) { // Protected against (0,0,z)
        pPhi = std::atan2(p.y(), p.x());

        if (pPhi < fSPhi - kHalfTolerance)
          pPhi += 2 * kPi;
        else if (pPhi > fSPhi + fDPhi + kHalfTolerance)
          pPhi -= 2 * kPi;

        distSPhi = rho * (pPhi - fSPhi);
        distEPhi = rho * (pPhi - fSPhi - fDPhi);
        inside   = inside && (distSPhi > -kTolerance) && (distEPhi < kTolerance);
        distSPhi = std::abs(distSPhi);
        distEPhi = std::abs(distEPhi);
      }

      else if (!(fRmin1) || !(fRmin2)) {
        distSPhi = 0.;
        distEPhi = 0.;
      }
      nPs = Vector3D<Precision>(std::sin(fSPhi), -std::cos(fSPhi), 0);
      nPe = Vector3D<Precision>(-std::sin(fSPhi + fDPhi), std::cos(fSPhi + fDPhi), 0);
    }

    if (rho > kHalfTolerance) {
      nR = Vector3D<Precision>(p.x() / rho / fSecRMax, p.y() / rho / fSecRMax, -fTanRMax / fSecRMax);
      if (fRmin1 || fRmin2) {
        nr = Vector3D<Precision>(-p.x() / rho / fSecRMin, -p.y() / rho / fSecRMin, fTanRMin / fSecRMin);
      }
    }

    if (inside && distZ <= kHalfTolerance) {
      noSurfaces++;
      if (p.z() >= 0.)
        sumnorm += nZ;
      else
        sumnorm.Set(0, 0, -1.);
    }

    if (inside && distRMax <= kHalfTolerance) {
      noSurfaces++;
      sumnorm += nR;
    } else if (noSurfaces == 0 && distRMax < distNearest) {
      distNearest = distRMax;
      normNearest = nR;
    }

    if (fRmin1 || fRmin2) {
      if (inside && distRMin <= kHalfTolerance) {
        noSurfaces++;
        sumnorm += nr;
      } else if (noSurfaces == 0 && distRMin < distNearest) {
        distNearest = distRMin;
        normNearest = nr;
      }
    }

    if (!IsFullPhi()) {
      if (inside && distSPhi <= kHalfTolerance) {
        noSurfaces++;
        sumnorm += nPs;
      } else if (noSurfaces == 0 && distSPhi < distNearest) {
        distNearest = distSPhi;
        normNearest = nPs;
      }
      if (inside && distEPhi <= kHalfTolerance) {
        noSurfaces++;
        sumnorm += nPe;
      } else if (noSurfaces == 0 && distEPhi < distNearest) {
        // No more check on distNearest, no need to assign to it.
        // distNearest = distEPhi;
        normNearest = nPe;
      }
    }
    // Final checks
    if (noSurfaces == 0)
      norm = normNearest;
    else if (noSurfaces == 1)
      norm = sumnorm;
    else
      norm = sumnorm.Unit();
    return noSurfaces != 0;
  }

  VECCORE_ATT_HOST_DEVICE
  void SetAndCheckSPhiAngle(Precision sPhi)
  {
    // Ensure fSphi in 0-2PI or -2PI-0 range if shape crosses 0
    if (sPhi < 0) {
      fSPhi = kTwoPi - std::fmod(std::fabs(sPhi), kTwoPi);
    } else {
      fSPhi = std::fmod(sPhi, kTwoPi);
    }
    if (fSPhi + fDPhi > kTwoPi) {
      fSPhi -= kTwoPi;
    }

    // Update Wedge
    fPhiWedge.SetStartPhi(fSPhi);
    // Update cached values.
    GetAlongVectorToPhiSector(fSPhi, fAlongPhi1x, fAlongPhi1y);
    GetAlongVectorToPhiSector(fSPhi + fDPhi, fAlongPhi2x, fAlongPhi2y);
  }

  VECCORE_ATT_HOST_DEVICE
  void SetAndCheckDPhiAngle(Precision dPhi)
  {
    if (dPhi >= kTwoPi - 0.5 * kAngTolerance) {
      fDPhi = kTwoPi;
      fSPhi = 0;
    } else {
      if (dPhi > 0) {
        fDPhi = dPhi;
      } else {
        //        std::ostringstream message;
        //        message << "Invalid dphi.\n"
        //                << "Negative or zero delta-Phi (" << dPhi << ")\n";
        //        std::cout<<"UnplacedTube::CheckDPhiAngle(): Fatal error: "<< message.str().c_str() <<"\n";
      }
    }
    // Update Wedge
    fPhiWedge.SetDeltaPhi(fDPhi);
    // Update cached values.
    GetAlongVectorToPhiSector(fSPhi, fAlongPhi1x, fAlongPhi1y);
    GetAlongVectorToPhiSector(fSPhi + fDPhi, fAlongPhi2x, fAlongPhi2y);
  }

  VECCORE_ATT_HOST_DEVICE
  static void GetAlongVectorToPhiSector(Precision phi, Precision &x, Precision &y)
  {
    x = std::cos(phi);
    y = std::sin(phi);
  }

  void SetRmin1(Precision const &arg)
  {
    fRmin1 = arg;
    CalculateCached();
  }
  void SetRmax1(Precision const &arg)
  {
    fRmax1 = arg;
    CalculateCached();
  }
  void SetRmin2(Precision const &arg)
  {
    fRmin2 = arg;
    CalculateCached();
  }
  void SetRmax2(Precision const &arg)
  {
    fRmax2 = arg;
    CalculateCached();
  }
  void SetDz(Precision const &arg)
  {
    fDz = arg;
    CalculateCached();
  }
  void SetSPhi(Precision const &arg)
  {
    fSPhi = arg;
    SetAndCheckSPhiAngle(fSPhi);
    // DetectConvexity();
  }
  void SetDPhi(Precision const &arg)
  {
    fDPhi = arg;
    SetAndCheckDPhiAngle(fDPhi);
    // DetectConvexity();
  }

  VECCORE_ATT_HOST_DEVICE
  Precision GetTolIz() const { return fTolIz; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetTolOz() const { return fTolOz; }

  VECCORE_ATT_HOST_DEVICE
  Precision GetConeTolerane() const { return fConeTolerance; }

  VECCORE_ATT_HOST_DEVICE
  evolution::Wedge const &GetWedge() const { return fPhiWedge; }

  // constructors
  VECCORE_ATT_HOST_DEVICE
  ConeStruct(T const &_rmin1, T const &_rmax1, T const &_rmin2, T const &_rmax2, T const &_z, T const &_sphi,
             T const &_dphi)
      : fRmin1(_rmin1 < 0.0 ? 0.0 : _rmin1), fRmax1(_rmax1), fRmin2(_rmin2 < 0.0 ? 0.0 : _rmin2), fRmax2(_rmax2),
        fDz(_z), fSPhi(_sphi), fDPhi(_dphi), _frmin1(_rmin1), _frmin2(_rmin2), _frmax1(_rmax1), _frmax2(_rmax2),
        fPhiWedge(_dphi, _sphi)
  {

    SetAndCheckDPhiAngle(_dphi);
    SetAndCheckSPhiAngle(_sphi);
    CalculateCached();

    // DetectConvexity();
  }
};
}
} // end global namespace

#endif
