/*
 * UnplacedCone.cpp
 *
 *  Created on: Jun 18, 2014
 *      Author: swenzel
 */

#include "volumes/UnplacedCone.h"
#include "volumes/SpecializedCone.h"
#include "volumes/utilities/VolumeUtilities.h"
#include "volumes/utilities/GenerationUtilities.h"
#ifndef VECGEOM_NVCC
#include "base/RNG.h"
#endif
#include "management/VolumeFactory.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

void UnplacedCone::Print() const
{
  printf("UnplacedCone {rmin1 %.2f, rmax1 %.2f, rmin2 %.2f, "
         "rmax2 %.2f, dz %.2f, phistart %.2f, deltaphi %.2f}",
         fRmin1, fRmax1, fRmin2, fRmax2, fDz, fSPhi, fDPhi);
}

void UnplacedCone::Print(std::ostream &os) const
{
  os << "UnplacedCone; please implement Print to outstream\n";
}

VECGEOM_CUDA_HEADER_BOTH
void UnplacedCone::DetectConvexity()
{

  // Default safe convexity value
  fGlobalConvexity = false;

  // Logic to calculate the convexity
  if (fRmin1 == 0. && fRmin2 == 0.) { // Implies Solid cone
    if (fDPhi <= kPi || fDPhi == kTwoPi) fGlobalConvexity = true;
  }
}

#if (0)
// Simplest Extent definition, that does not take PHI into consideration
void UnplacedCone::void Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const
{
  Precision max = fRmax1 > fRmax2 ? fRmax1 : fRmax2;
  aMin          = Vector3D<Precision>(-max, -max, -fDz);
  aMax          = Vector3D<Precision>(max, max, fDz);
}
#endif

#if (1)
// Improved Extent definition, that takes PHI also into consideration
void UnplacedCone::Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const
{
  // most general case

  Precision max = fRmax1 > fRmax2 ? fRmax1 : fRmax2;
  Precision min = fRmin1 > fRmin2 ? fRmin2 : fRmin1;

  aMin = Vector3D<Precision>(-max, -max, -fDz);
  aMax = Vector3D<Precision>(max, max, fDz);

  /* Below logic borrowed from Tube.
  **
  ** But it would be great, if its possible to directly call Extent of Tube.
  ** because in that case we can avoid code replication.
  */

  if (fDPhi == kTwoPi) return;

  // check how many of phi=90, 180, 270, 360deg are outside this tube
  auto Rin       = 0.5 * (max + min);
  bool phi0out   = !GetWedge().Contains<kScalar>(Vector3D<Precision>(Rin, 0, 0));
  bool phi90out  = !GetWedge().Contains<kScalar>(Vector3D<Precision>(0, Rin, 0));
  bool phi180out = !GetWedge().Contains<kScalar>(Vector3D<Precision>(-Rin, 0, 0));
  bool phi270out = !GetWedge().Contains<kScalar>(Vector3D<Precision>(0, -Rin, 0));

  // if none of those 4 phis is outside, largest box still required
  if (!(phi0out || phi90out || phi180out || phi270out)) return;

  // some extent(s) of box will be reduced
  // --> think of 4 points A,B,C,D such that A,B are at Rmin, C,D at Rmax
  //     and A,C at startPhi (fSphi), B,D at endPhi (fSphi+fDphi)
  auto Cx = max * cos(fSPhi);
  auto Dx = max * cos(fSPhi + fDPhi);
  auto Cy = max * sin(fSPhi);
  auto Dy = max * sin(fSPhi + fDPhi);

  // then rewrite box sides whenever each one of those phis are not contained in the tube section
  if (phi0out) aMax.x()   = Max(Cx, Dx);
  if (phi90out) aMax.y()  = Max(Cy, Dy);
  if (phi180out) aMin.x() = Min(Cx, Dx);
  if (phi270out) aMin.y() = Min(Cy, Dy);

  if (fDPhi >= kPi) return;

  auto Ax = min * cos(fSPhi);
  auto Bx = min * cos(fSPhi + fDPhi);
  auto Ay = min * sin(fSPhi);
  auto By = min * sin(fSPhi + fDPhi);

  Precision temp;
  temp     = Max(Ax, Bx);
  aMax.x() = temp > aMax.x() ? temp : aMax.x();

  temp     = Max(Ay, By);
  aMax.y() = temp > aMax.y() ? temp : aMax.y();

  temp     = Min(Ax, Bx);
  aMin.x() = temp < aMin.x() ? temp : aMin.x();

  temp     = Min(Ay, By);
  aMin.y() = temp < aMin.y() ? temp : aMin.y();

  return;
}
#endif

#if !defined(VECGEOM_NVCC)
bool UnplacedCone::Normal(Vector3D<Precision> const &p, Vector3D<Precision> &norm) const
{
  int noSurfaces = 0;
  Precision rho, pPhi;
  Precision distZ, distRMin, distRMax;
  Precision distSPhi = kInfLength, distEPhi = kInfLength;
  Precision pRMin, widRMin;
  Precision pRMax, widRMax;

  const double kHalfTolerance = 0.5 * kTolerance;

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
      distNearest = distEPhi;
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

Vector3D<Precision> UnplacedCone::GetPointOnSurface() const
{
  // implementation taken from UCons; not verified
  //
  double Aone, Atwo, Athree, Afour, Afive, slin, slout, phi;
  double zRand, cosu, sinu, rRand1, rRand2, chose, rone, rtwo, qone, qtwo;
  rone = (fRmax1 - fRmax2) / (2. * fDz);
  rtwo = (fRmin1 - fRmin2) / (2. * fDz);
  qone = 0.;
  qtwo = 0.;
  if (fRmax1 != fRmax2) {
    qone = fDz * (fRmax1 + fRmax2) / (fRmax1 - fRmax2);
  }
  if (fRmin1 != fRmin2) {
    qtwo = fDz * (fRmin1 + fRmin2) / (fRmin1 - fRmin2);
  }
  slin   = Sqrt((fRmin1 - fRmin2) * (fRmin1 - fRmin2) + 4. * fDz * fDz);
  slout  = Sqrt((fRmax1 - fRmax2) * (fRmax1 - fRmax2) + 4. * fDz * fDz);
  Aone   = 0.5 * fDPhi * (fRmax2 + fRmax1) * slout;
  Atwo   = 0.5 * fDPhi * (fRmin2 + fRmin1) * slin;
  Athree = 0.5 * fDPhi * (fRmax1 * fRmax1 - fRmin1 * fRmin1);
  Afour  = 0.5 * fDPhi * (fRmax2 * fRmax2 - fRmin2 * fRmin2);
  Afive  = fDz * (fRmax1 - fRmin1 + fRmax2 - fRmin2);

  phi    = RNG::Instance().uniform(fSPhi, fSPhi + fDPhi);
  cosu   = std::cos(phi);
  sinu   = std::sin(phi);
  rRand1 = volumeUtilities::GetRadiusInRing(fRmin1, fRmin2);
  rRand2 = volumeUtilities::GetRadiusInRing(fRmax1, fRmax2);

  if ((fSPhi == 0.) && IsFullPhi()) {
    Afive = 0.;
  }
  chose = RNG::Instance().uniform(0., Aone + Atwo + Athree + Afour + 2. * Afive);

  if ((chose >= 0.) && (chose < Aone)) {
    if (fRmin1 != fRmin2) {
      zRand = RNG::Instance().uniform(-1. * fDz, fDz);
      return Vector3D<Precision>(rtwo * cosu * (qtwo - zRand), rtwo * sinu * (qtwo - zRand), zRand);
    } else {
      return Vector3D<Precision>(fRmin1 * cosu, fRmin2 * sinu, RNG::Instance().uniform(-1. * fDz, fDz));
    }
  } else if ((chose >= Aone) && (chose <= Aone + Atwo)) {
    if (fRmax1 != fRmax2) {
      zRand = RNG::Instance().uniform(-1. * fDz, fDz);
      return Vector3D<Precision>(rone * cosu * (qone - zRand), rone * sinu * (qone - zRand), zRand);
    } else {
      return Vector3D<Precision>(fRmax1 * cosu, fRmax2 * sinu, RNG::Instance().uniform(-1. * fDz, fDz));
    }
  } else if ((chose >= Aone + Atwo) && (chose < Aone + Atwo + Athree)) {
    return Vector3D<Precision>(rRand1 * cosu, rRand1 * sinu, -1 * fDz);
  } else if ((chose >= Aone + Atwo + Athree) && (chose < Aone + Atwo + Athree + Afour)) {
    return Vector3D<Precision>(rRand2 * cosu, rRand2 * sinu, fDz);
  } else if ((chose >= Aone + Atwo + Athree + Afour) && (chose < Aone + Atwo + Athree + Afour + Afive)) {
    zRand  = RNG::Instance().uniform(-1. * fDz, fDz);
    rRand1 = RNG::Instance().uniform(fRmin2 - ((zRand - fDz) / (2. * fDz)) * (fRmin1 - fRmin2),
                                     fRmax2 - ((zRand - fDz) / (2. * fDz)) * (fRmax1 - fRmax2));
    return Vector3D<Precision>(rRand1 * std::cos(fSPhi), rRand1 * std::sin(fSPhi), zRand);
  } else {
    zRand  = RNG::Instance().uniform(-1. * fDz, fDz);
    rRand1 = RNG::Instance().uniform(fRmin2 - ((zRand - fDz) / (2. * fDz)) * (fRmin1 - fRmin2),
                                     fRmax2 - ((zRand - fDz) / (2. * fDz)) * (fRmax1 - fRmax2));
    return Vector3D<Precision>(rRand1 * std::cos(fSPhi + fDPhi), rRand1 * std::sin(fSPhi + fDPhi), zRand);
  }
}
#endif // VECGEOM_NVCC

template <TranslationCode transCodeT, RotationCode rotCodeT>
VECGEOM_CUDA_HEADER_DEVICE
VPlacedVolume *UnplacedCone::Create(LogicalVolume const *const logical_volume,
                                    Transformation3D const *const transformation,
#ifdef VECGEOM_NVCC
                                    const int id,
#endif
                                    VPlacedVolume *const placement)
{

  using namespace ConeTypes;
  __attribute__((unused)) const UnplacedCone &cone =
      static_cast<const UnplacedCone &>(*(logical_volume->GetUnplacedVolume()));

#ifdef VECGEOM_NVCC
#define RETURN_SPECIALIZATION(coneTypeT)                                                   \
  return CreateSpecializedWithPlacement<SpecializedCone<transCodeT, rotCodeT, coneTypeT>>( \
      logical_volume, transformation, id, placement)
#else
#define RETURN_SPECIALIZATION(coneTypeT)                                                                  \
  return CreateSpecializedWithPlacement<SpecializedCone<transCodeT, rotCodeT, coneTypeT>>(logical_volume, \
                                                                                          transformation, placement)
#endif

#ifdef GENERATE_CONE_SPECIALIZATIONS
  if (cone.GetRmin1() <= 0 && cone.GetRmin2() <= 0) {
    if (cone.GetDPhi() >= 2 * M_PI) RETURN_SPECIALIZATION(NonHollowCone);
    if (cone.GetDPhi() == M_PI) RETURN_SPECIALIZATION(NonHollowConeWithPiSector); // == M_PI ???

    if (cone.GetDPhi() < M_PI) RETURN_SPECIALIZATION(NonHollowConeWithSmallerThanPiSector);
    if (cone.GetDPhi() > M_PI) RETURN_SPECIALIZATION(NonHollowConeWithBiggerThanPiSector);
  } else if (cone.GetRmin1() > 0 || cone.GetRmin2() > 0) {
    if (cone.GetDPhi() >= 2 * M_PI) RETURN_SPECIALIZATION(HollowCone);
    if (cone.GetDPhi() == M_PI) RETURN_SPECIALIZATION(HollowConeWithPiSector); // == M_PI ???
    if (cone.GetDPhi() < M_PI) RETURN_SPECIALIZATION(HollowConeWithSmallerThanPiSector);
    if (cone.GetDPhi() > M_PI) RETURN_SPECIALIZATION(HollowConeWithBiggerThanPiSector);
  }
#endif

  RETURN_SPECIALIZATION(UniversalCone);

#undef RETURN_SPECIALIZATION
}

#if defined(VECGEOM_USOLIDS)
VECGEOM_CUDA_HEADER_BOTH
std::ostream &UnplacedCone::StreamInfo(std::ostream &os) const
{
  int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "     *** Dump for solid - " << GetEntityType() << " ***\n"
     << "     ===================================================\n"
     << " Solid type: Cone\n"
     << " Parameters: \n"
     << "     Cone Radii Rmin1, Rmax1: " << fRmin1 << "mm, " << fRmax1 << "mm\n"
     << "                Rmin2, Rmax2: " << fRmin2 << "mm, " << fRmax2 << "mm\n"
     << "     Half-length Z = " << fDz << "mm\n";
  if (fDPhi < kTwoPi) {
    os << "     Wedge starting angles: fSPhi=" << fSPhi * kRadToDeg << "deg, "
       << ", fDphi=" << fDPhi * kRadToDeg << "deg\n";
  }
  os << "-----------------------------------------------------------\n";
  os.precision(oldprc);
  return os;
}
#endif

// this is repetitive code:

VECGEOM_CUDA_HEADER_DEVICE
VPlacedVolume *UnplacedCone::SpecializedVolume(LogicalVolume const *const volume,
                                               Transformation3D const *const transformation,
                                               const TranslationCode trans_code, const RotationCode rot_code,
#ifdef VECGEOM_NVCC
                                               const int id,
#endif
                                               VPlacedVolume *const placement) const
{

  return VolumeFactory::CreateByTransformation<UnplacedCone>(volume, transformation, trans_code, rot_code,
#ifdef VECGEOM_NVCC
                                                             id,
#endif
                                                             placement);
}

#ifdef VECGEOM_CUDA_INTERFACE

DevicePtr<cuda::VUnplacedVolume> UnplacedCone::CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const in_gpu_ptr) const
{
  return CopyToGpuImpl<UnplacedCone>(in_gpu_ptr, GetRmin1(), GetRmax1(), GetRmin2(), GetRmax2(), GetDz(), GetSPhi(),
                                     GetDPhi());
}

DevicePtr<cuda::VUnplacedVolume> UnplacedCone::CopyToGpu() const
{
  return CopyToGpuImpl<UnplacedCone>();
}

#endif // VECGEOM_CUDA_INTERFACE

} // End impl namespace

#ifdef VECGEOM_NVCC

namespace cxx {

template size_t DevicePtr<cuda::UnplacedCone>::SizeOf();
template void DevicePtr<cuda::UnplacedCone>::Construct(const Precision rmin1, const Precision rmax1,
                                                       const Precision rmin2, const Precision rmax2, const Precision z,
                                                       const Precision sphi, const Precision dphi) const;

} // End cxx namespace

#endif

} // End namespace vecgeom
