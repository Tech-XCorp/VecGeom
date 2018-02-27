/// \file UnplacedSphere.cpp
/// \author Raman Sehgal (raman.sehgal@cern.ch)

#include "volumes/SphereUtilities.h"
#include "volumes/UnplacedSphere.h"
#include "volumes/SpecializedSphere.h"
#include "volumes/utilities/VolumeUtilities.h"
#include "volumes/utilities/GenerationUtilities.h"
#ifndef VECCORE_CUDA
#include "base/RNG.h"
#endif
#include "management/VolumeFactory.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

VECCORE_ATT_HOST_DEVICE
UnplacedSphere::UnplacedSphere(Precision pRmin, Precision pRmax, Precision pSPhi, Precision pDPhi, Precision pSTheta,
                               Precision pDTheta)
    : fSphere(pRmin, pRmax, pSPhi, pDPhi, pSTheta, pDTheta)
{

  DetectConvexity();
}

VECCORE_ATT_HOST_DEVICE
void UnplacedSphere::DetectConvexity()
{
  // Default Convexity set to false
  fGlobalConvexity = false;
  // Logic to calculate the convexity
  if (fSphere.fRmin == 0.) {
    if (((fSphere.fDPhi == kTwoPi) && (fSphere.fSTheta == 0.) && (fSphere.eTheta == kPi)) ||
        ((fSphere.fDPhi <= kPi) && (fSphere.fSTheta == 0) && (fSphere.eTheta == kPi)) ||
        ((fSphere.fDPhi == kTwoPi) && (fSphere.fSTheta == 0) && (fSphere.eTheta <= kPi / 2)) ||
        ((fSphere.fDPhi == kTwoPi) && (fSphere.fSTheta >= kPi / 2) && (fSphere.eTheta == kPi)))
      fGlobalConvexity = true;
  }
}

#if (0)
// Simplest Implementation of Extent
void UnplacedSphere::Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const
{
  // Returns the full 3D cartesian extent of the solid.
  aMin.Set(-fRmax);
  aMax.Set(fRmax);
}
#endif

#ifndef VECCORE_CUDA
#if (1)
// Sophisticated Implementation taking into account the PHI and THETA cut also.
void UnplacedSphere::Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const
{
  // most general case
  aMin.Set(-fSphere.fRmax);
  aMax.Set(fSphere.fRmax);
  Precision eTheta = fSphere.fSTheta + fSphere.fDTheta;

  if (fSphere.fFullSphere) return;

  Precision st1 = 0.;
  Precision st2 = 0.;

  if (!fSphere.fFullThetaSphere) {
    // Simplified logic suggested by Evgueni Tcherniaev
    aMax.z() = ((fSphere.fSTheta <= kPi / 2.) ? fSphere.fRmax : fSphere.fRmin) * Cos(fSphere.fSTheta);
    aMin.z() = ((eTheta <= kPi / 2.) ? fSphere.fRmin : fSphere.fRmax) * Cos(eTheta);

    st1 = Sin(fSphere.fSTheta);
    st2 = Sin(fSphere.eTheta);

    if (fSphere.fSTheta <= kPi / 2.) {
      if ((eTheta) < kPi / 2.) {
        aMax.x() = fSphere.fRmax * st2;
        aMin.x() = -fSphere.fRmax * st2;
        aMax.y() = fSphere.fRmax * st2;
        aMin.y() = -fSphere.fRmax * st2;
      }
    } else {
      aMax.x() = fSphere.fRmax * st1;
      aMin.x() = -fSphere.fRmax * st1;
      aMax.y() = fSphere.fRmax * st1;
      aMin.y() = -fSphere.fRmax * st1;
    }
  }
  if (!fSphere.fFullPhiSphere) {
    Precision Rmax = fSphere.fRmax;
    if (fSphere.fSTheta > kPi / 2.) Rmax *= st1;
    if (eTheta < kPi / 2.) Rmax *= st2;
    Precision Rmin = fSphere.fRmin * Min(st1, st2);
    // These newly calculated Rmin and Rmax will be used by PHI section

    // Borrowed PHI logic from Tube
    // check how many of phi=90, 180, 270, 360deg are outside this tube

    auto Rin       = 0.5 * (Rmax + Rmin);
    bool phi0out   = !fSphere.fPhiWedge.Contains(Vector3D<Precision>(Rin, 0, 0));
    bool phi90out  = !fSphere.fPhiWedge.Contains(Vector3D<Precision>(0, Rin, 0));
    bool phi180out = !fSphere.fPhiWedge.Contains(Vector3D<Precision>(-Rin, 0, 0));
    bool phi270out = !fSphere.fPhiWedge.Contains(Vector3D<Precision>(0, -Rin, 0));

    // if none of those 4 phis is outside, largest box still required
    if (!(phi0out || phi90out || phi180out || phi270out)) return;

    // some extent(s) of box will be reduced
    // --> think of 4 points A,B,C,D such that A,B are at Rmin, C,D at Rmax
    //     and A,C at startPhi (fSphi), B,D at endPhi (fSphi+fDphi)
    auto Cx = Rmax * cos(fSphere.fSPhi);
    auto Dx = Rmax * cos(fSphere.fSPhi + fSphere.fDPhi);
    auto Cy = Rmax * sin(fSphere.fSPhi);
    auto Dy = Rmax * sin(fSphere.fSPhi + fSphere.fDPhi);

    // then rewrite box sides whenever each one of those phis are not contained in the tube section
    if (phi0out) aMax.x()   = Max(Cx, Dx);
    if (phi90out) aMax.y()  = Max(Cy, Dy);
    if (phi180out) aMin.x() = Min(Cx, Dx);
    if (phi270out) aMin.y() = Min(Cy, Dy);

    if (fSphere.fDPhi >= kPi) return;

    auto Ax = Rmin * cos(fSphere.fSPhi);
    auto Bx = Rmin * cos(fSphere.fSPhi + fSphere.fDPhi);
    auto Ay = Rmin * sin(fSphere.fSPhi);
    auto By = Rmin * sin(fSphere.fSPhi + fSphere.fDPhi);

    Precision temp;
    temp     = Max(Ax, Bx);
    aMax.x() = temp > aMax.x() ? temp : aMax.x();

    temp     = Max(Ay, By);
    aMax.y() = temp > aMax.y() ? temp : aMax.y();

    temp     = Min(Ax, Bx);
    aMin.x() = temp < aMin.x() ? temp : aMin.x();

    temp     = Min(Ay, By);
    aMin.y() = temp < aMin.y() ? temp : aMin.y();
  }
  return;
}
#endif

#endif // !VECCORE_CUDA

void UnplacedSphere::GetParametersList(int, double *aArray) const
{
  aArray[0] = GetInnerRadius();
  aArray[1] = GetOuterRadius();
  aArray[2] = GetStartPhiAngle();
  aArray[3] = GetDeltaPhiAngle();
  aArray[4] = GetStartThetaAngle();
  aArray[5] = GetDeltaThetaAngle();
}

#ifndef VECCORE_CUDA
Vector3D<Precision> UnplacedSphere::SamplePointOnSurface() const
{

  Precision zRand, aOne, aTwo, aThr, aFou, aFiv, chose, phi, sinphi, cosphi;
  Precision height1, height2, slant1, slant2, costheta, sintheta, rRand;

  height1 = (fSphere.fRmax - fSphere.fRmin) * fSphere.cosSTheta;
  height2 = (fSphere.fRmax - fSphere.fRmin) * fSphere.cosETheta;
  slant1  = std::sqrt(sqr((fSphere.fRmax - fSphere.fRmin) * fSphere.sinSTheta) + height1 * height1);
  slant2  = std::sqrt(sqr((fSphere.fRmax - fSphere.fRmin) * fSphere.sinETheta) + height2 * height2);
  rRand   = GetRadiusInRing(fSphere.fRmin, fSphere.fRmax);

  aOne = fSphere.fRmax * fSphere.fRmax * fSphere.fDPhi * (fSphere.cosSTheta - fSphere.cosETheta);
  aTwo = fSphere.fRmin * fSphere.fRmin * fSphere.fDPhi * (fSphere.cosSTheta - fSphere.cosETheta);
  aThr = fSphere.fDPhi * ((fSphere.fRmax + fSphere.fRmin) * fSphere.sinSTheta) * slant1;
  aFou = fSphere.fDPhi * ((fSphere.fRmax + fSphere.fRmin) * fSphere.sinETheta) * slant2;
  aFiv = 0.5 * fSphere.fDTheta * (fSphere.fRmax * fSphere.fRmax - fSphere.fRmin * fSphere.fRmin);

  phi      = RNG::Instance().uniform(fSphere.fSPhi, fSphere.ePhi);
  cosphi   = std::cos(phi);
  sinphi   = std::sin(phi);
  costheta = RNG::Instance().uniform(fSphere.cosETheta, fSphere.cosSTheta);
  sintheta = std::sqrt(1. - sqr(costheta));

  if (fSphere.fFullPhiSphere) {
    aFiv = 0;
  }
  if (fSphere.fSTheta == 0) {
    aThr = 0;
  }
  if (fSphere.eTheta == kPi) {
    aFou = 0;
  }
  if (fSphere.fSTheta == kPi / 2) {
    aThr = kPi * (fSphere.fRmax * fSphere.fRmax - fSphere.fRmin * fSphere.fRmin);
  }
  if (fSphere.eTheta == kPi / 2) {
    aFou = kPi * (fSphere.fRmax * fSphere.fRmax - fSphere.fRmin * fSphere.fRmin);
  }

  chose = RNG::Instance().uniform(0., aOne + aTwo + aThr + aFou + 2. * aFiv);
  if ((chose >= 0.) && (chose < aOne)) {
    return Vector3D<Precision>(fSphere.fRmax * sintheta * cosphi, fSphere.fRmax * sintheta * sinphi,
                               fSphere.fRmax * costheta);
  } else if ((chose >= aOne) && (chose < aOne + aTwo)) {
    return Vector3D<Precision>(fSphere.fRmin * sintheta * cosphi, fSphere.fRmin * sintheta * sinphi,
                               fSphere.fRmin * costheta);
  } else if ((chose >= aOne + aTwo) && (chose < aOne + aTwo + aThr)) {
    if (fSphere.fSTheta != kPi / 2) {
      zRand = RNG::Instance().uniform(fSphere.fRmin * fSphere.cosSTheta, fSphere.fRmax * fSphere.cosSTheta);
      return Vector3D<Precision>(fSphere.tanSTheta * zRand * cosphi, fSphere.tanSTheta * zRand * sinphi, zRand);
    } else {
      return Vector3D<Precision>(rRand * cosphi, rRand * sinphi, 0.);
    }
  } else if ((chose >= aOne + aTwo + aThr) && (chose < aOne + aTwo + aThr + aFou)) {
    if (fSphere.eTheta != kPi / 2) {
      zRand = RNG::Instance().uniform(fSphere.fRmin * fSphere.cosETheta, fSphere.fRmax * fSphere.cosETheta);
      return Vector3D<Precision>(fSphere.tanETheta * zRand * cosphi, fSphere.tanETheta * zRand * sinphi, zRand);
    } else {
      return Vector3D<Precision>(rRand * cosphi, rRand * sinphi, 0.);
    }
  } else if ((chose >= aOne + aTwo + aThr + aFou) && (chose < aOne + aTwo + aThr + aFou + aFiv)) {
    return Vector3D<Precision>(rRand * sintheta * fSphere.cosSPhi, rRand * sintheta * fSphere.sinSPhi,
                               rRand * costheta);
  } else {
    return Vector3D<Precision>(rRand * sintheta * fSphere.cosEPhi, rRand * sintheta * fSphere.sinEPhi,
                               rRand * costheta);
  }
}

std::string UnplacedSphere::GetEntityType() const
{
  return "Sphere\n";
}

#endif // !VECCORE_CUDA

VECCORE_ATT_HOST_DEVICE
void UnplacedSphere::ComputeBBox() const
{
}

/*UnplacedSphere *UnplacedSphere::Clone() const
{
  return new UnplacedSphere(fSphere.fRmin, fSphere.fRmax, fSphere.fSPhi, fSphere.fDPhi, fSphere.fSTheta,
fSphere.fDTheta);
}*/

std::ostream &UnplacedSphere::StreamInfo(std::ostream &os) const
// Definition taken from USphere
{

  int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     //  << "     *** Dump for solid - " << GetName() << " ***\n"
     //  << "     ===================================================\n"

     << " Solid type: VecGeomSphere\n"
     << " Parameters: \n"

     << "       outer radius: " << fSphere.fRmax << " mm \n"
     << "               Inner radius: " << fSphere.fRmin << "mm\n"
     << "               Start Phi Angle: " << fSphere.fSPhi << "\n"
     << "               Delta Phi Angle: " << fSphere.fDPhi << "\n"
     << "               Start Theta Angle: " << fSphere.fSTheta << "\n"
     << "               Delta Theta Angle: " << fSphere.fDTheta << "\n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

void UnplacedSphere::Print() const
{
  printf("UnplacedSphere {%.2f , %.2f , %.2f , %.2f , %.2f , %.2f}", GetInnerRadius(), GetOuterRadius(),
         GetStartPhiAngle(), GetDeltaPhiAngle(), GetStartThetaAngle(), GetDeltaThetaAngle());
}

void UnplacedSphere::Print(std::ostream &os) const
{
  os << "UnplacedSphere { " << GetInnerRadius() << " " << GetOuterRadius() << " " << GetStartPhiAngle() << " "
     << GetDeltaPhiAngle() << " " << GetStartThetaAngle() << " " << GetDeltaThetaAngle() << " }";
}

#ifndef VECCORE_CUDA

template <TranslationCode trans_code, RotationCode rot_code>
VPlacedVolume *UnplacedSphere::Create(LogicalVolume const *const logical_volume,
                                      Transformation3D const *const transformation, VPlacedVolume *const placement)
{
  if (placement) {
    new (placement) SpecializedSphere<trans_code, rot_code>(logical_volume, transformation);
    return placement;
  }
  return new SpecializedSphere<trans_code, rot_code>(logical_volume, transformation);
}

VPlacedVolume *UnplacedSphere::CreateSpecializedVolume(LogicalVolume const *const volume,
                                                       Transformation3D const *const transformation,
                                                       const TranslationCode trans_code, const RotationCode rot_code,
                                                       VPlacedVolume *const placement)
{
  return VolumeFactory::CreateByTransformation<UnplacedSphere>(volume, transformation, trans_code, rot_code, placement);
}

#else

template <TranslationCode trans_code, RotationCode rot_code>
VECCORE_ATT_DEVICE
VPlacedVolume *UnplacedSphere::Create(LogicalVolume const *const logical_volume,
                                      Transformation3D const *const transformation, const int id,
                                      VPlacedVolume *const placement)
{
  if (placement) {
    new (placement) SpecializedSphere<trans_code, rot_code>(logical_volume, transformation, NULL, id);
    return placement;
  }
  return new SpecializedSphere<trans_code, rot_code>(logical_volume, transformation, NULL, id);
}

VECCORE_ATT_DEVICE VPlacedVolume *UnplacedSphere::CreateSpecializedVolume(LogicalVolume const *const volume,
                                                                          Transformation3D const *const transformation,
                                                                          const TranslationCode trans_code,
                                                                          const RotationCode rot_code, const int id,
                                                                          VPlacedVolume *const placement)
{
  return VolumeFactory::CreateByTransformation<UnplacedSphere>(volume, transformation, trans_code, rot_code, id,
                                                               placement);
}

#endif

#ifdef VECGEOM_CUDA_INTERFACE

DevicePtr<cuda::VUnplacedVolume> UnplacedSphere::CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const in_gpu_ptr) const
{
  return CopyToGpuImpl<UnplacedSphere>(in_gpu_ptr, GetInnerRadius(), GetOuterRadius(), GetStartPhiAngle(),
                                       GetDeltaPhiAngle(), GetStartThetaAngle(), GetDeltaThetaAngle());
}

DevicePtr<cuda::VUnplacedVolume> UnplacedSphere::CopyToGpu() const
{
  return CopyToGpuImpl<UnplacedSphere>();
}

#endif // VECGEOM_CUDA_INTERFACE

} // End impl namespace

#ifdef VECCORE_CUDA

namespace cxx {

template size_t DevicePtr<cuda::UnplacedSphere>::SizeOf();
template void DevicePtr<cuda::UnplacedSphere>::Construct(const Precision rmin, const Precision rmax,
                                                         const Precision sphi, const Precision dphi,
                                                         const Precision stheta, const Precision dtheta) const;

} // End cxx namespace

#endif

} // End global namespace
