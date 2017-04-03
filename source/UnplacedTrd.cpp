/// @file UnplacedTrd.cpp
/// @author Georgios Bitzes (georgios.bitzes@cern.ch)

#include "volumes/UnplacedTrd.h"
#include "volumes/SpecializedTrd.h"
#include "volumes/utilities/GenerationUtilities.h"
#include "base/RNG.h"

#include "management/VolumeFactory.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

void UnplacedTrd::Print() const
{
  printf("UnplacedTrd {%.2f, %.2f, %.2f, %.2f, %.2f}", dx1(), dx2(), dy1(), dy2(), dz());
}

void UnplacedTrd::Print(std::ostream &os) const
{
  os << "UnplacedTrd {" << dx1() << ", " << dx2() << ", " << dy1() << ", " << dy2() << ", " << dz();
}

#ifndef VECGEOM_NVCC
Precision UnplacedTrd::Capacity() const
{
  return 2 * (fTrd.fDX1 + fTrd.fDX2) * (fTrd.fDY1 + fTrd.fDY2) * fTrd.fDZ +
         (2. / 3.) * (fTrd.fDX1 - fTrd.fDX2) * (fTrd.fDY1 - fTrd.fDY2) * fTrd.fDZ;
}

Precision UnplacedTrd::SurfaceArea() const
{
  Precision dz = 2 * fTrd.fDZ;
  bool xvert   = (fTrd.fDX1 == fTrd.fDX2) ? true : false;
  Precision SA = 0.0;

  // Sum of area for planes Perp. to +X and -X
  Precision ht = (xvert) ? dz : Sqrt((fTrd.fDX1 - fTrd.fDX2) * (fTrd.fDX1 - fTrd.fDX2) + dz * dz);
  SA += 2.0 * 2.0 * 0.5 * (fTrd.fDY1 + fTrd.fDY2) * ht;

  // Sum of area for planes Perp. to +Y and -Y
  SA += 2.0 * 2.0 * 0.5 * (fTrd.fDX1 + fTrd.fDX2) * ht; // if xvert then topology forces to become yvert for closing

  // Sum of area for top and bottom planes +Z and -Z
  SA += 4. * (fTrd.fDX1 * fTrd.fDY1) + 4. * (fTrd.fDX2 * fTrd.fDY2);

  return SA;
}

int UnplacedTrd::ChooseSurface() const
{
  int choice = 0; // 0 = zm, 1 = zp, 2 = ym, 3 = yp, 4 = xm, 5 = xp
  Precision S[6], Stotal = 0.0;

  S[0] = S[1] = GetPlusXArea();
  S[2] = S[3] = GetPlusYArea();
  S[4]        = GetMinusZArea();
  S[5]        = GetPlusZArea();

  for (int i = 0; i < 6; ++i)
    Stotal += S[i];

  // random value to choose surface to place the point
  Precision rand = RNG::Instance().uniform() * Stotal;

  while (rand > S[choice])
    rand -= S[choice], choice++;

  assert(choice < 6);
  return choice;
}

Vector3D<Precision> UnplacedTrd::GetPointOnSurface() const
{
  int surface = ChooseSurface();

  Precision invHeight = 0.5 / fTrd.fDZ;
  Precision slopeX    = (fTrd.fDX2 - fTrd.fDX1) * invHeight;
  Precision slopeY    = (fTrd.fDY2 - fTrd.fDY1) * invHeight;
  Precision midX      = 0.5 * (fTrd.fDX1 + fTrd.fDX2);
  Precision midY      = 0.5 * (fTrd.fDY1 + fTrd.fDY2);

  // Note that randoms are drawn in range [-1,1], then we just need scales, e.g. xval *= xmax for x = [-xmax, xmax] and
  // so on
  Precision xval = RNG::Instance().uniform(-1.0, 1.0);
  Precision yval = RNG::Instance().uniform(-1.0, 1.0);
  Precision zval = RNG::Instance().uniform(-1.0, 1.0);
  Precision vmax;
  switch (surface) {
  case 0:              // -Z face
    xval *= fTrd.fDX1; // = ran[-fDX1, fDX1]
    yval *= fTrd.fDY1;
    zval = -fTrd.fDZ;
    break;
  case 1:              // +Z face
    xval *= fTrd.fDX2; // = ran[-fDX2, fDX2]
    yval *= fTrd.fDY2;
    zval = fTrd.fDZ;
    break;
  case 2: // +X face
    xval = midX + slopeX * zval;
    vmax = midY + slopeY * zval;
    yval *= vmax; // = rand[-vmax,vmax]
    break;
  case 3: // -X face
    xval = -(midX + slopeX * zval);
    vmax = midY + slopeY * zval;
    yval *= vmax; // = rand[-vmax,vmax]
    break;
  case 4: // +Y face
    vmax = midX + slopeX * zval;
    xval *= vmax; // = rand[-vmax,vmax]
    yval = midY + slopeY * zval;
    break;
  case 5: // -Y face
    vmax = midX + slopeX * zval;
    xval *= vmax; // = rand[-vmax,vmax]
    yval = -(midY + slopeY * zval);
    break;
  }

  return Vector3D<Precision>(xval, yval, zval);
}

bool UnplacedTrd::Normal(Vector3D<Precision> const &point, Vector3D<Precision> &norm) const
{

  int noSurfaces = 0;
  Vector3D<Precision> sumnorm(0., 0., 0.), vecnorm(0., 0., 0.);
  Precision distz;

  distz = std::fabs(std::fabs(point[2]) - fTrd.fDZ);

  Precision xnorm = 1.0 / sqrt(4 * fTrd.fDZ * fTrd.fDZ + (fTrd.fDX2 - fTrd.fDX1) * (fTrd.fDX2 - fTrd.fDX1));
  Precision ynorm = 1.0 / sqrt(4 * fTrd.fDZ * fTrd.fDZ + (fTrd.fDY2 - fTrd.fDY1) * (fTrd.fDY2 - fTrd.fDY1));

  Precision distmx =
      -2.0 * fTrd.fDZ * point[0] - (fTrd.fDX2 - fTrd.fDX1) * point[2] - fTrd.fDZ * (fTrd.fDX1 + fTrd.fDX2);
  distmx *= xnorm;

  Precision distpx =
      2.0 * fTrd.fDZ * point[0] - (fTrd.fDX2 - fTrd.fDX1) * point[2] - fTrd.fDZ * (fTrd.fDX1 + fTrd.fDX2);
  distpx *= xnorm;

  Precision distmy =
      -2.0 * fTrd.fDZ * point[1] - (fTrd.fDY2 - fTrd.fDY1) * point[2] - fTrd.fDZ * (fTrd.fDY1 + fTrd.fDY2);
  distmy *= ynorm;

  Precision distpy =
      2.0 * fTrd.fDZ * point[1] - (fTrd.fDY2 - fTrd.fDY1) * point[2] - fTrd.fDZ * (fTrd.fDY1 + fTrd.fDY2);
  distpy *= ynorm;

  if (fabs(distmx) <= kHalfTolerance) {
    noSurfaces++;
    sumnorm += Vector3D<Precision>(-2.0 * fTrd.fDZ, 0.0, -(fTrd.fDX2 - fTrd.fDX1)) * xnorm;
  }
  if (fabs(distpx) <= kHalfTolerance) {
    noSurfaces++;
    sumnorm += Vector3D<Precision>(2.0 * fTrd.fDZ, 0.0, -(fTrd.fDX2 - fTrd.fDX1)) * xnorm;
  }
  if (fabs(distpy) <= kHalfTolerance) {
    noSurfaces++;
    sumnorm += Vector3D<Precision>(0.0, 2.0 * fTrd.fDZ, -(fTrd.fDY2 - fTrd.fDY1)) * ynorm;
  }
  if (fabs(distmy) <= kHalfTolerance) {
    noSurfaces++;
    sumnorm += Vector3D<Precision>(0.0, -2.0 * fTrd.fDZ, -(fTrd.fDY2 - fTrd.fDY1)) * ynorm;
  }

  if (std::fabs(distz) <= kHalfTolerance) {
    noSurfaces++;
    if (point[2] >= 0.)
      sumnorm += Vector3D<Precision>(0., 0., 1.);
    else
      sumnorm += Vector3D<Precision>(0., 0., -1.);
  }
  if (noSurfaces == 0) {
#ifdef UDEBUG
    UUtils::Exception("UnplacedTrapezoid::SurfaceNormal(point)", "GeomSolids1002", Warning, 1,
                      "Point is not on surface.");
#endif
    // vecnorm = ApproxSurfaceNormal( Vector3D<Precision>(point[0],point[1],point[2]) );
    vecnorm =
        Vector3D<Precision>(0., 0., 1.); // any plane will do it, since false is returned, so save the CPU cycles...
  } else if (noSurfaces == 1)
    vecnorm = sumnorm;
  else
    vecnorm = sumnorm.Unit();

  norm[0] = vecnorm[0];
  norm[1] = vecnorm[1];
  norm[2] = vecnorm[2];

  return noSurfaces != 0;
}

#endif

template <TranslationCode transCodeT, RotationCode rotCodeT>
VECCORE_ATT_DEVICE
VPlacedVolume *UnplacedTrd::Create(LogicalVolume const *const logical_volume,
                                   Transformation3D const *const transformation,
#ifdef VECGEOM_NVCC
                                   const int id,
#endif
                                   VPlacedVolume *const placement)
{

  using namespace TrdTypes;

#ifndef VECGEOM_NO_SPECIALIZATION

  __attribute__((unused)) const UnplacedTrd &trd =
      static_cast<const UnplacedTrd &>(*(logical_volume->GetUnplacedVolume()));

#define GENERATE_TRD_SPECIALIZATIONS
#ifdef GENERATE_TRD_SPECIALIZATIONS
  if (trd.dy1() == trd.dy2()) {
    //          std::cout << "trd1" << std::endl;
    return CreateSpecializedWithPlacement<SpecializedTrd<transCodeT, rotCodeT, TrdTypes::Trd1>>(logical_volume,
                                                                                                transformation
#ifdef VECGEOM_NVCC
                                                                                                ,
                                                                                                id
#endif
                                                                                                ,
                                                                                                placement);
  } else {
    //          std::cout << "trd2" << std::endl;
    return CreateSpecializedWithPlacement<SpecializedTrd<transCodeT, rotCodeT, TrdTypes::Trd2>>(logical_volume,
                                                                                                transformation
#ifdef VECGEOM_NVCC
                                                                                                ,
                                                                                                id
#endif
                                                                                                ,
                                                                                                placement);
  }
#endif

#endif // VECGEOM_NO_SPECIALIZATION

  //    std::cout << "universal trd" << std::endl;
  return CreateSpecializedWithPlacement<SpecializedTrd<transCodeT, rotCodeT, TrdTypes::UniversalTrd>>(logical_volume,
                                                                                                      transformation
#ifdef VECGEOM_NVCC
                                                                                                      ,
                                                                                                      id
#endif
                                                                                                      ,
                                                                                                      placement);
}

VECCORE_ATT_DEVICE
VPlacedVolume *UnplacedTrd::SpecializedVolume(LogicalVolume const *const volume,
                                              Transformation3D const *const transformation,
                                              const TranslationCode trans_code, const RotationCode rot_code,
#ifdef VECGEOM_NVCC
                                              const int id,
#endif
                                              VPlacedVolume *const placement) const
{

  return VolumeFactory::CreateByTransformation<UnplacedTrd>(volume, transformation, trans_code, rot_code,
#ifdef VECGEOM_NVCC
                                                            id,
#endif
                                                            placement);
}

#if defined(VECGEOM_USOLIDS)
VECCORE_ATT_HOST_DEVICE
std::ostream &UnplacedTrd::StreamInfo(std::ostream &os) const
{
  int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "     *** Dump for solid - " << GetEntityType() << " ***\n"
     << "     ===================================================\n"
     << " Solid type: Trd\n"
     << " Parameters: \n"
     << "     half lengths X1,X2: " << fTrd.fDX1 << "mm, " << fTrd.fDX2 << "mm \n"
     << "     half lengths Y1,Y2: " << fTrd.fDY1 << "mm, " << fTrd.fDY2 << "mm \n"
     << "     half length Z: " << fTrd.fDZ << "mm \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);
  return os;
}
#endif

#ifdef VECGEOM_CUDA_INTERFACE

DevicePtr<cuda::VUnplacedVolume> UnplacedTrd::CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const in_gpu_ptr) const
{
  return CopyToGpuImpl<UnplacedTrd>(in_gpu_ptr, dx1(), dx2(), dy1(), dy2(), dz());
}

DevicePtr<cuda::VUnplacedVolume> UnplacedTrd::CopyToGpu() const
{
  return CopyToGpuImpl<UnplacedTrd>();
}

#endif // VECGEOM_CUDA_INTERFACE

} // End impl namespace

#ifdef VECGEOM_NVCC

namespace cxx {

template size_t DevicePtr<cuda::UnplacedTrd>::SizeOf();
template void DevicePtr<cuda::UnplacedTrd>::Construct(const Precision dx1, const Precision dx2, const Precision dy1,
                                                      const Precision dy2, const Precision z) const;

} // End cxx namespace

#endif

} // End global namespace
