/// \file UnplacedTube.cpp
/// \author Georgios Bitzes (georgios.bitzes@cern.ch)

#include "volumes/UnplacedTube.h"
#include "volumes/SpecializedTube.h"
#ifndef VECCORE_CUDA
#include "base/RNG.h"
#include <cmath>
#include <iostream>
#endif

#include "volumes/utilities/GenerationUtilities.h"
#include "management/VolumeFactory.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

void UnplacedTube::Print() const
{
  printf("UnplacedTube {%.2f, %.2f, %.2f, %.2f, %.2f}", rmin(), rmax(), z(), sphi(), dphi());
}

void UnplacedTube::Print(std::ostream &os) const
{
  os << "UnplacedTube {" << rmin() << ", " << rmax() << ", " << z() << ", " << sphi() << ", " << dphi() << "}\n";
}

// template <TranslationCode transCodeT, RotationCode rotCodeT>
// VECCORE_ATT_DEVICE
// VPlacedVolume *UnplacedTube::Create(LogicalVolume const *const logical_volume,
//                                    Transformation3D const *const transformation,
//#ifdef VECCORE_CUDA
//                                    const int id,
//#endif
//                                    VPlacedVolume *const placement)
//{

//  using namespace TubeTypes;
//  __attribute__((unused)) const UnplacedTube &tube =
//      static_cast<const UnplacedTube &>(*(logical_volume->GetUnplacedVolume()));

////#ifdef VECCORE_CUDA
////#define RETURN_SPECIALIZATION(tubeTypeT)				"\"
////  return CreateSpecializedWithPlacement<SpecializedTube<transCodeT, rotCodeT, tubeTypeT>>( "\"
//      logical_volume, transformation, id, placement)
//#else
////#define RETURN_SPECIALIZATION(tubeTypeT)				"\"
////  return CreateSpecializedWithPlacement<SpecializedTube<transCodeT, rotCodeT, tubeTypeT>>(logical_volume, "\"
//                                                                                          transformation, placement)
//#endif

//#ifdef GENERATE_TUBE_SPECIALIZATIONS
//  if (tube.rmin() <= 0) {
//    if (tube.dphi() >= 2 * M_PI) RETURN_SPECIALIZATION(NonHollowTube);
//    if (tube.dphi() == M_PI) RETURN_SPECIALIZATION(NonHollowTubeWithPiSector); // == M_PI ???

//    if (tube.dphi() < M_PI) RETURN_SPECIALIZATION(NonHollowTubeWithSmallerThanPiSector);
//    if (tube.dphi() > M_PI) RETURN_SPECIALIZATION(NonHollowTubeWithBiggerThanPiSector);
//  } else if (tube.rmin() > 0) {
//    if (tube.dphi() >= 2 * M_PI) RETURN_SPECIALIZATION(HollowTube);
//    if (tube.dphi() == M_PI) RETURN_SPECIALIZATION(HollowTubeWithPiSector); // == M_PI ???

//    if (tube.dphi() < M_PI) RETURN_SPECIALIZATION(HollowTubeWithSmallerThanPiSector);
//    if (tube.dphi() > M_PI) RETURN_SPECIALIZATION(HollowTubeWithBiggerThanPiSector);
//  }
//#endif

//  RETURN_SPECIALIZATION(UniversalTube);

//#undef RETURN_SPECIALIZATION
//}

#ifndef VECCORE_CUDA
VPlacedVolume *UnplacedTube::SpecializedVolume(LogicalVolume const *const volume,
                                               Transformation3D const *const transformation,
                                               const TranslationCode trans_code, const RotationCode rot_code,
                                               VPlacedVolume *const placement) const
{
  return VolumeFactory::CreateByTransformation<UnplacedTube>(volume, transformation, trans_code, rot_code, placement);
}
#else
VECCORE_ATT_DEVICE VPlacedVolume *UnplacedTube::SpecializedVolume(LogicalVolume const *const volume,
                                                                  Transformation3D const *const transformation,
                                                                  const TranslationCode trans_code,
                                                                  const RotationCode rot_code, const int id,
                                                                  VPlacedVolume *const placement) const
{
  return VolumeFactory::CreateByTransformation<UnplacedTube>(volume, transformation, trans_code, rot_code, id,
                                                             placement);
}
#endif

int UnplacedTube::ChooseSurface() const
{
  int choice = 0; // 0 = rTop, 1 = rBot, 2 = phiLeft, 3 = phiRight, 4 = zIn, 5 = zOut
  Precision S[6], Stotal = 0.0;

  S[0] = S[1] = 0.5 * GetTopArea();        // 50% divide into top and bottom
  S[2] = S[3] = 0.5 * GetLateralPhiArea(); // 50% divide into left and right
  S[4]        = GetLateralRInArea();       // inner tube surface area
  S[5]        = GetLateralROutArea();      // outer tube surface area

  for (int i = 0; i < 6; ++i)
    Stotal += S[i];

  /* random value to choose surface to place the point */
  Precision rand = RNG::Instance().uniform() * Stotal;

  while (rand > S[choice])
    rand -= S[choice], choice++;

  assert(choice < 6);

  return choice;
}

Vector3D<Precision> UnplacedTube::SamplePointOnSurface() const
{
  int surface      = ChooseSurface();
  Precision rVal   = RNG::Instance().uniform(rmin(), rmax());
  Precision phiVal = RNG::Instance().uniform(sphi(), sphi() + dphi());
  Precision zVal   = RNG::Instance().uniform() * 2.0 * z() - z();

  switch (surface) {
  case 0:
    zVal = z();
    break;
  case 1:
    zVal = -z();
    break;
  case 2:
    phiVal = sphi();
    break;
  case 3:
    phiVal = sphi() + dphi();
    break;
  case 4:
    rVal = rmin();
    break;
  case 5:
    rVal = rmax();
    break;
  }

  Precision xVal = rVal * cos(phiVal);
  Precision yVal = rVal * sin(phiVal);

  return Vector3D<Precision>(xVal, yVal, zVal);
}

bool UnplacedTube::Normal(Vector3D<Precision> const &point, Vector3D<Precision> &norm) const
{
  bool valid = true;
  TubeImplementation<TubeTypes::UniversalTube>::NormalKernel<double, bool>(fTube, point, norm, valid);
  return valid;
}

/*
  VECCORE_ATT_HOST_DEVICE
  Precision UnplacedTube::SurfaceArea () const {
    Precision area = fDphi * (rmin() + rmax()) * (2 * fZ + rmax() - rmin());
    if (fDphi<kTwoPi) {
      area += 4 * fZ * (rmax() - rmin());
    }
    return area;
  }

  */

VECCORE_ATT_HOST_DEVICE
void UnplacedTube::DetectConvexity()
{

  // Default safe convexity value
  fGlobalConvexity = false;

  // Logic to calculate the convexity
  if (rmin() == 0.) {
    if (dphi() <= kPi || dphi() == kTwoPi) fGlobalConvexity = true;
  }
}

void UnplacedTube::Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const
{
  // most general case
  aMin = Vector3D<Precision>(-rmax(), -rmax(), -z());
  aMax = Vector3D<Precision>(rmax(), rmax(), z());

  if (dphi() == kTwoPi) return;

  // check how many of phi=90, 180, 270, 360deg are outside this tube
  auto Rin       = 0.5 * (rmax() + rmin());
  bool phi0out   = !GetWedge().Contains(Vector3D<Precision>(Rin, 0, 0));
  bool phi90out  = !GetWedge().Contains(Vector3D<Precision>(0, Rin, 0));
  bool phi180out = !GetWedge().Contains(Vector3D<Precision>(-Rin, 0, 0));
  bool phi270out = !GetWedge().Contains(Vector3D<Precision>(0, -Rin, 0));

  // if none of those 4 phis is outside, largest box still required
  if (!(phi0out || phi90out || phi180out || phi270out)) return;

  // some extent(s) of box will be reduced
  // --> think of 4 points A,B,C,D such that A,B are at Rmin, C,D at Rmax
  //     and A,C at startPhi (fSphi), B,D at endPhi (fSphi+fDphi)
  auto Cx = rmax() * cos(sphi());
  auto Dx = rmax() * cos(sphi() + dphi());
  auto Cy = rmax() * sin(sphi());
  auto Dy = rmax() * sin(sphi() + dphi());

  // then rewrite box sides whenever each one of those phis are not contained in the tube section
  if (phi0out) aMax.x()   = Max(Cx, Dx);
  if (phi90out) aMax.y()  = Max(Cy, Dy);
  if (phi180out) aMin.x() = Min(Cx, Dx);
  if (phi270out) aMin.y() = Min(Cy, Dy);

  if (dphi() >= kPi) return;

  auto Ax = rmin() * cos(sphi());
  auto Bx = rmin() * cos(sphi() + dphi());
  auto Ay = rmin() * sin(sphi());
  auto By = rmin() * sin(sphi() + dphi());

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

#ifdef VECGEOM_CUDA_INTERFACE

DevicePtr<cuda::VUnplacedVolume> UnplacedTube::CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const in_gpu_ptr) const
{
  return CopyToGpuImpl<SUnplacedTube<TubeTypes::UniversalTube>>(in_gpu_ptr, rmin(), rmax(), z(), sphi(), dphi());
}

DevicePtr<cuda::VUnplacedVolume> UnplacedTube::CopyToGpu() const
{
  return CopyToGpuImpl<SUnplacedTube<TubeTypes::UniversalTube>>();
}

#endif // VECGEOM_CUDA_INTERFACE

} // End impl namespace

#ifdef VECCORE_CUDA

namespace cxx {

template size_t DevicePtr<cuda::SUnplacedTube<cuda::TubeTypes::UniversalTube>>::SizeOf();
template void DevicePtr<cuda::SUnplacedTube<cuda::TubeTypes::UniversalTube>>::Construct(
    const Precision rmin, const Precision rmax, const Precision z, const Precision sphi, const Precision dphi) const;

} // End cxx namespace

#endif

} // End global namespace
