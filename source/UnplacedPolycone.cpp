/*
 * UnplacedPolycone.cpp
 *
 *  Created on: Dec 8, 2014
 *      Author: swenzel
 */

#include "volumes/UnplacedPolycone.h"
#include "volumes/UnplacedCone.h"
#include "volumes/SpecializedPolycone.h"
#include "volumes/kernel/ConeImplementation.h"
#include "management/VolumeFactory.h"
#ifndef VECGEOM_NVCC
#include "base/RNG.h"
#endif
#include <iostream>
#include <cstdio>
#include <vector>
#include "base/Vector.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

//
// Constructor (GEANT3 style parameters)
//
void UnplacedPolycone::Init(double phiStart, double phiTotal, unsigned int numZPlanes, const double zPlane[],
                            const double rInner[], const double rOuter[])
{
  // Conversion for angles
  if (phiTotal <= 0. || phiTotal > kTwoPi - kTolerance) {
    // phiIsOpen=false;
    fStartPhi = 0;
    fDeltaPhi = kTwoPi;
  } else {
    //
    // Convert phi into our convention
    //
    fStartPhi = phiStart;
    while (fStartPhi < 0)
      fStartPhi += kTwoPi;
  }

  // Calculate RMax of Polycone in order to determine convexity of sections
  //
  double RMaxextent = rOuter[0];

  Vector<Precision> newROuter, newZPlane, newRInner;
  fContinuityOverAll &= CheckContinuity(rOuter, rInner, zPlane, newROuter, newRInner, newZPlane);
  fConvexityPossible &= (newRInner[0] == 0.);

  Precision startRmax = newROuter[0];
  for (unsigned int j = 1; j < newROuter.size(); j++) {
    fEqualRmax &= (startRmax == newROuter[j]);
    startRmax = newROuter[j];
    fConvexityPossible &= (newRInner[j] == 0.);
  }

  for (unsigned int j = 1; j < numZPlanes; j++) {

    if (rOuter[j] > RMaxextent) RMaxextent = rOuter[j];

    if (rInner[j] > rOuter[j]) {
#ifndef VECGEOM_NVCC
      std::cerr << "Cannot create Polycone with rInner > rOuter for the same Z"
                << "\n"
                << "        rInner > rOuter for the same Z !\n"
                << "        rMin[" << j << "] = " << rInner[j] << " -- rMax[" << j << "] = " << rOuter[j];
#endif
    }
  }

  double prevZ = zPlane[0], prevRmax = 0, prevRmin = 0;
  int dirZ                        = 1;
  if (zPlane[1] < zPlane[0]) dirZ = -1;

  for (unsigned int i = 0; i < numZPlanes; ++i) {
    if ((i < numZPlanes - 1) && (zPlane[i] == zPlane[i + 1])) {
      if ((rInner[i] > rOuter[i + 1]) || (rInner[i + 1] > rOuter[i])) {
#ifndef VECGEOM_NVCC
        std::cerr << "Cannot create a Polycone with no contiguous segments." << std::endl
                  << "                Segments are not contiguous !" << std::endl
                  << "                rMin[" << i << "] = " << rInner[i] << " -- rMax[" << i + 1
                  << "] = " << rOuter[i + 1] << std::endl
                  << "                rMin[" << i + 1 << "] = " << rInner[i + 1] << " -- rMax[" << i
                  << "] = " << rOuter[i];
#endif
      }
    }

    double rMin = rInner[i];

    double rMax = rOuter[i];
    double z    = zPlane[i];

    // i has to be at least one to complete a section
    if (i > 0) {
      if (((z > prevZ) && (dirZ > 0)) || ((z < prevZ) && (dirZ < 0))) {
        if (dirZ * (z - prevZ) < 0) {
#ifndef VECGEOM_NVCC
          std::cerr << "Cannot create a Polycone with different Z directions.Use GenericPolycone." << std::endl
                    << "              ZPlane is changing direction  !" << std::endl
                    << "  zPlane[0] = " << zPlane[0] << " -- zPlane[1] = " << zPlane[1] << std::endl
                    << "  zPlane[" << i - 1 << "] = " << zPlane[i - 1] << " -- rPlane[" << i << "] = " << zPlane[i];
#endif
        }

        UnplacedCone *solid;

        double dz = (z - prevZ) / 2;

        solid = new UnplacedCone(prevRmin, prevRmax, rMin, rMax, dz, phiStart, phiTotal);

        fZs.push_back(z);
        int zi       = fZs.size() - 1;
        double shift = fZs[zi - 1] + 0.5 * (fZs[zi] - fZs[zi - 1]);

        PolyconeSection section;
        section.fShift = shift;
        section.fSolid = solid;

        section.fConvex = !((rMax < prevRmax) || (rMax < RMaxextent) || (prevRmax < RMaxextent));

        fSections.push_back(section);
      }
    } else { // for i == 0 just push back first z plane
      fZs.push_back(z);
    }

    prevZ    = z;
    prevRmin = rMin;
    prevRmax = rMax;
  }
  DetectConvexity();
}

// Alternative constructor, required for integration with Geant4.
// Input must be such that r[i],z[i] describe the outer,inner or inner,outer envelope of the polycone, after
// connecting all adjacent points, and closing the polygon by connecting last -> first point.
// Hence z[] array must be symmetrical: z[0..Nz] = z[2Nz, 2Nz-1, ..., Nz+1], where 2*Nz = numRz.
VECGEOM_CUDA_HEADER_BOTH
UnplacedPolycone::UnplacedPolycone(Precision phiStart, // initial phi starting angle
                                   Precision phiTotal, // total phi angle
                                   int numRZ,          // number corners in r,z space (must be an even number)
                                   Precision const *r, // r coordinate of these corners
                                   Precision const *z) // z coordinate of these corners
    : fStartPhi(phiStart),
      fDeltaPhi(phiStart + phiTotal),
      fNz(numRZ / 2),
      fSections(),
      fZs(numRZ / 2)
{
  // data integrity checks
  int Nz = numRZ / 2;
  assert(numRZ % 2 == 0 && "UnplPolycone ERROR: r[],z[] arrays provided contain odd number of points, please fix.\n");
  for (int i = 0; i < numRZ / 2; ++i) {
    assert(z[i] == z[numRZ - 1 - i] && "UnplPolycone ERROR: z[] array is not symmetrical, please fix.\n");
  }

  // reuse input array as argument, in ascending order
  bool ascendingZ        = true;
  const Precision *zarg  = z;
  const Precision *r1arg = r;
  if (z[0] > z[1]) {
    ascendingZ = false;
    zarg       = z + Nz; // second half of input z[] is ascending due to symmetry already verified
    r1arg      = r + Nz;
  }

  // reorganize remainder of r[] data in ascending-z order
  Precision *r2arg = new Precision[Nz];
  for (int i = 0; i < Nz; ++i)
    r2arg[i] = (ascendingZ ? r[2 * Nz - 1 - i] : r[Nz - 1 - i]);

  // identify which rXarg is rmax and rmin and ensure that Rmax > Rmin for all points provided
  const Precision *rmin = r1arg, *rmax = r2arg;
  if (r1arg[0] > r2arg[0]) {
    rmax = r1arg;
    rmin = r2arg;
  }
  delete[] r2arg;

  // final data integrity cross-check
  for (int i = 0; i < Nz; ++i) {
    assert(rmax[i] > rmin[i] &&
           "UnplPolycone ERROR: r[] provided has problems of the Rmax < Rmin type, please check!\n");
  }

  // init internal members
  Init(phiStart, phiTotal, Nz, zarg, rmin, rmax);
}

template <TranslationCode transCodeT, RotationCode rotCodeT>
VECGEOM_CUDA_HEADER_DEVICE
VPlacedVolume *UnplacedPolycone::Create(LogicalVolume const *const logical_volume,
                                        Transformation3D const *const transformation,
#ifdef VECGEOM_NVCC
                                        const int id,
#endif
                                        VPlacedVolume *const placement)
{

  if (placement) {
    new (placement) SpecializedPolycone<transCodeT, rotCodeT>(logical_volume, transformation
#ifdef VECGEOM_NVCC
                                                              ,
                                                              NULL, id
#endif
                                                              );
    return placement;
  }
  return new SpecializedPolycone<transCodeT, rotCodeT>(logical_volume, transformation
#ifdef VECGEOM_NVCC
                                                       ,
                                                       NULL, id
#endif
                                                       );
}

void UnplacedPolycone::Print() const
{
  printf("UnplacedPolycone {%.2f, %.2f, %d}\n", fStartPhi, fDeltaPhi, fNz);
  printf("have %zu size Z\n", fZs.size());
  printf("------- z planes follow ---------\n");
  for (size_t p = 0; p < fZs.size(); ++p) {
    printf(" plane %zu at z pos %lf\n", p, fZs[p]);
  }

  printf("have %zu size fSections\n", fSections.size());
  printf("------ sections follow ----------\n");
  for (int s = 0; s < GetNSections(); ++s) {
    printf("## section %d, shift %lf\n", s, fSections[s].fShift);
    fSections[s].fSolid->Print();
    printf("\n");
  }
}

void UnplacedPolycone::Print(std::ostream &os) const
{
  os << "UnplacedPolycone output to string not implemented -- calling Print() instead:\n";
  Print();
}

#if defined(VECGEOM_USOLIDS)
VECGEOM_CUDA_HEADER_BOTH
std::ostream &UnplacedPolycone::StreamInfo(std::ostream &os) const
{
  int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "     *** Dump for solid - " << GetEntityType() << " ***\n"
     << "     ===================================================\n"
     << " Solid type: Polycone\n"
     << " Parameters: \n"
     << "     N = number of Z-sections: " << fSections.size() << ", # Z-coords=" << fZs.size() << "\n"
     << "     z-coordinates:\n";

  uint nz = fZs.size();
  for (uint j = 0; j < (nz - 1) / 5 + 1; ++j) {
    os << "       [ ";
    for (uint i = 0; i < 5; ++i) {
      uint ind = 5 * j + i;
      if (ind < fNz) os << fZs[ind] << "; ";
    }
    os << " ]\n";
  }
  if (fDeltaPhi < kTwoPi) {
    os << "     Wedge starting angles: fSphi=" << fStartPhi * kRadToDeg << "deg, "
       << ", fDphi=" << fDeltaPhi * kRadToDeg << "deg\n";
  }

  size_t nsections = fSections.size();
  os << "\n    # cone sections: " << nsections << "\n";
  for (size_t i = 0; i < nsections; ++i) {
    UnplacedCone *subcone = fSections[i].fSolid;
    os << "     cone #" << i << " Rmin1=" << subcone->GetRmin1() << " Rmax1=" << subcone->GetRmax1()
       << " Rmin2=" << subcone->GetRmin2() << " Rmax2=" << subcone->GetRmax2() << " HalfZ=" << subcone->GetDz()
       << " from z=" << fZs[i] << " to z=" << fZs[i + 1] << "mm\n";
  }
  os << "-----------------------------------------------------------\n";
  os.precision(oldprc);
  return os;
}
#endif

VECGEOM_CUDA_HEADER_DEVICE
VPlacedVolume *UnplacedPolycone::SpecializedVolume(LogicalVolume const *const volume,
                                                   Transformation3D const *const transformation,
                                                   const TranslationCode trans_code, const RotationCode rot_code,
#ifdef VECGEOM_NVCC
                                                   const int id,
#endif
                                                   VPlacedVolume *const placement) const
{

  // TODO: for the Polycone this might be overkill
  return VolumeFactory::CreateByTransformation<UnplacedPolycone>(volume, transformation, trans_code, rot_code,
#ifdef VECGEOM_NVCC
                                                                 id,
#endif
                                                                 placement);
}

#ifdef VECGEOM_CUDA_INTERFACE

DevicePtr<cuda::VUnplacedVolume> UnplacedPolycone::CopyToGpu() const
{
  return CopyToGpuImpl<UnplacedPolycone>();
}

DevicePtr<cuda::VUnplacedVolume> UnplacedPolycone::CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const gpu_ptr) const
{

  // idea: reconstruct defining arrays: copy them to GPU; then construct the UnplacedPolycon object from scratch
  // on the GPU
  std::vector<Precision> rmin, z, rmax;
  ReconstructSectionArrays(z, rmin, rmax);

  // somehow this does not work:
  //        Precision *z_gpu_ptr = AllocateOnGpu<Precision>( (z.size() + rmin.size() + rmax.size())*sizeof(Precision) );
  //        Precision *rmin_gpu_ptr = z_gpu_ptr + sizeof(Precision)*z.size();
  //        Precision *rmax_gpu_ptr = rmin_gpu_ptr + sizeof(Precision)*rmin.size();

  Precision *z_gpu_ptr    = AllocateOnGpu<Precision>(z.size() * sizeof(Precision));
  Precision *rmin_gpu_ptr = AllocateOnGpu<Precision>(rmin.size() * sizeof(Precision));
  Precision *rmax_gpu_ptr = AllocateOnGpu<Precision>(rmax.size() * sizeof(Precision));

  vecgeom::CopyToGpu(&z[0], z_gpu_ptr, sizeof(Precision) * z.size());
  vecgeom::CopyToGpu(&rmin[0], rmin_gpu_ptr, sizeof(Precision) * rmin.size());
  vecgeom::CopyToGpu(&rmax[0], rmax_gpu_ptr, sizeof(Precision) * rmax.size());

  vecgeom::CopyToGpu(&z[0], z_gpu_ptr, sizeof(Precision) * z.size());
  vecgeom::CopyToGpu(&rmin[0], rmin_gpu_ptr, sizeof(Precision) * rmin.size());
  vecgeom::CopyToGpu(&rmax[0], rmax_gpu_ptr, sizeof(Precision) * rmax.size());

  int s = z.size();

  // attention here z.size() might be different than fNz due to compactification during Reconstruction
  DevicePtr<cuda::VUnplacedVolume> gpupolycon =
      CopyToGpuImpl<UnplacedPolycone>(gpu_ptr, fStartPhi, fDeltaPhi, s, z_gpu_ptr, rmin_gpu_ptr, rmax_gpu_ptr);

  // remove temporary space from GPU
  FreeFromGpu(z_gpu_ptr);
  FreeFromGpu(rmin_gpu_ptr);
  FreeFromGpu(rmax_gpu_ptr);

  return gpupolycon;
}

#endif // VECGEOM_CUDA_INTERFACE

#ifndef VECGEOM_NVCC
/////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface
//
// GetPointOnCone
//
// Auxiliary method for Get Point On Surface
//

Vector3D<Precision> UnplacedPolycone::GetPointOnCone(Precision fRmin1, Precision fRmax1, Precision fRmin2,
                                                     Precision fRmax2, Precision zOne, Precision zTwo,
                                                     Precision &totArea) const
{
  // declare working variables
  //
  Precision Aone, Atwo, Afive, phi, zRand, fDPhi, cosu, sinu;
  Precision rRand1, rmin, rmax, chose, rone, rtwo, qone, qtwo;
  Precision fDz = (zTwo - zOne) / 2., afDz = std::fabs(fDz);
  Vector3D<Precision> point, offset        = Vector3D<Precision>(0., 0., 0.5 * (zTwo + zOne));
  fDPhi = GetDeltaPhi();
  rone  = (fRmax1 - fRmax2) / (2. * fDz);
  rtwo  = (fRmin1 - fRmin2) / (2. * fDz);
  if (fRmax1 == fRmax2) {
    qone = 0.;
  } else {
    qone = fDz * (fRmax1 + fRmax2) / (fRmax1 - fRmax2);
  }
  if (fRmin1 == fRmin2) {
    qtwo = 0.;
  } else {
    qtwo = fDz * (fRmin1 + fRmin2) / (fRmin1 - fRmin2);
  }
  Aone    = 0.5 * fDPhi * (fRmax2 + fRmax1) * ((fRmin1 - fRmin2) * (fRmin1 - fRmin2) + (zTwo - zOne) * (zTwo - zOne));
  Atwo    = 0.5 * fDPhi * (fRmin2 + fRmin1) * ((fRmax1 - fRmax2) * (fRmax1 - fRmax2) + (zTwo - zOne) * (zTwo - zOne));
  Afive   = fDz * (fRmax1 - fRmin1 + fRmax2 - fRmin2);
  totArea = Aone + Atwo + 2. * Afive;

  phi  = RNG::Instance().uniform(GetStartPhi(), GetEndPhi());
  cosu = std::cos(phi);
  sinu = std::sin(phi);

  if (GetDeltaPhi() >= kTwoPi) {
    Afive = 0;
  }
  chose = RNG::Instance().uniform(0., Aone + Atwo + 2. * Afive);
  if ((chose >= 0) && (chose < Aone)) {
    if (fRmax1 != fRmax2) {
      zRand = RNG::Instance().uniform(-1. * afDz, afDz);
      point = Vector3D<Precision>(rone * cosu * (qone - zRand), rone * sinu * (qone - zRand), zRand);
    } else {
      point = Vector3D<Precision>(fRmax1 * cosu, fRmax1 * sinu, RNG::Instance().uniform(-1. * afDz, afDz));
    }
  } else if (chose >= Aone && chose < Aone + Atwo) {
    if (fRmin1 != fRmin2) {
      zRand = RNG::Instance().uniform(-1. * afDz, afDz);
      point = Vector3D<Precision>(rtwo * cosu * (qtwo - zRand), rtwo * sinu * (qtwo - zRand), zRand);

    } else {
      point = Vector3D<Precision>(fRmin1 * cosu, fRmin1 * sinu, RNG::Instance().uniform(-1. * afDz, afDz));
    }
  } else if ((chose >= Aone + Atwo + Afive) && (chose < Aone + Atwo + 2. * Afive)) {
    zRand  = RNG::Instance().uniform(-afDz, afDz);
    rmin   = fRmin2 - ((zRand - fDz) / (2. * fDz)) * (fRmin1 - fRmin2);
    rmax   = fRmax2 - ((zRand - fDz) / (2. * fDz)) * (fRmax1 - fRmax2);
    rRand1 = std::sqrt(RNG::Instance().uniform(0., 1.) * (rmax * rmax - rmin * rmin) + rmin * rmin);
    point  = Vector3D<Precision>(rRand1 * std::cos(GetStartPhi()), rRand1 * std::sin(GetStartPhi()), zRand);
  } else {
    zRand  = RNG::Instance().uniform(-1. * afDz, afDz);
    rmin   = fRmin2 - ((zRand - fDz) / (2. * fDz)) * (fRmin1 - fRmin2);
    rmax   = fRmax2 - ((zRand - fDz) / (2. * fDz)) * (fRmax1 - fRmax2);
    rRand1 = std::sqrt(RNG::Instance().uniform(0., 1.) * (rmax * rmax - rmin * rmin) + rmin * rmin);
    point  = Vector3D<Precision>(rRand1 * std::cos(GetEndPhi()), rRand1 * std::sin(GetEndPhi()), zRand);
  }

  return point + offset;
}

//
// GetPointOnTubs
//
// Auxiliary method for GetPoint On Surface
//
Vector3D<Precision> UnplacedPolycone::GetPointOnTubs(Precision fRMin, Precision fRMax, Precision zOne, Precision zTwo,
                                                     Precision &totArea) const
{
  Precision xRand, yRand, zRand, phi, cosphi, sinphi, chose, aOne, aTwo, aFou, rRand, fDz, fSPhi, fDPhi;
  fDz   = std::fabs(0.5 * (zTwo - zOne));
  fSPhi = GetStartPhi();
  fDPhi = GetDeltaPhi();

  aOne    = 2. * fDz * fDPhi * fRMax;
  aTwo    = 2. * fDz * fDPhi * fRMin;
  aFou    = 2. * fDz * (fRMax - fRMin);
  totArea = aOne + aTwo + 2. * aFou;
  phi     = RNG::Instance().uniform(GetStartPhi(), GetEndPhi());
  cosphi  = std::cos(phi);
  sinphi  = std::sin(phi);
  rRand   = fRMin + (fRMax - fRMin) * std::sqrt(RNG::Instance().uniform(0., 1.));

  if (GetDeltaPhi() >= 2 * kPi) aFou = 0;

  chose = RNG::Instance().uniform(0., aOne + aTwo + 2. * aFou);
  if ((chose >= 0) && (chose < aOne)) {
    xRand = fRMax * cosphi;
    yRand = fRMax * sinphi;
    zRand = RNG::Instance().uniform(-1. * fDz, fDz);
    return Vector3D<Precision>(xRand, yRand, zRand + 0.5 * (zTwo + zOne));
  } else if ((chose >= aOne) && (chose < aOne + aTwo)) {
    xRand = fRMin * cosphi;
    yRand = fRMin * sinphi;
    zRand = RNG::Instance().uniform(-1. * fDz, fDz);
    return Vector3D<Precision>(xRand, yRand, zRand + 0.5 * (zTwo + zOne));
  } else if ((chose >= aOne + aTwo) && (chose < aOne + aTwo + aFou)) {
    xRand = rRand * std::cos(fSPhi + fDPhi);
    yRand = rRand * std::sin(fSPhi + fDPhi);
    zRand = RNG::Instance().uniform(-1. * fDz, fDz);
    return Vector3D<Precision>(xRand, yRand, zRand + 0.5 * (zTwo + zOne));
  }

  // else

  xRand = rRand * std::cos(fSPhi + fDPhi);
  yRand = rRand * std::sin(fSPhi + fDPhi);
  zRand = RNG::Instance().uniform(-1. * fDz, fDz);
  return Vector3D<Precision>(xRand, yRand, zRand + 0.5 * (zTwo + zOne));
}

//
// GetPointOnRing
//
// Auxiliary method for GetPoint On Surface
//
Vector3D<Precision> UnplacedPolycone::GetPointOnRing(Precision fRMin1, Precision fRMax1, Precision fRMin2,
                                                     Precision fRMax2, Precision zOne) const
{
  Precision xRand, yRand, phi, cosphi, sinphi, rRand1, rRand2, A1, Atot, rCh;
  phi    = RNG::Instance().uniform(GetStartPhi(), GetEndPhi());
  cosphi = std::cos(phi);
  sinphi = std::sin(phi);

  if (fRMin1 == fRMin2) {
    rRand1 = fRMin1;
    A1     = 0.;
  } else {
    rRand1 = RNG::Instance().uniform(fRMin1, fRMin2);
    A1     = std::fabs(fRMin2 * fRMin2 - fRMin1 * fRMin1);
  }
  if (fRMax1 == fRMax2) {
    rRand2 = fRMax1;
    Atot   = A1;
  } else {
    rRand2 = RNG::Instance().uniform(fRMax1, fRMax2);
    Atot   = A1 + std::fabs(fRMax2 * fRMax2 - fRMax1 * fRMax1);
  }
  rCh = RNG::Instance().uniform(0., Atot);

  if (rCh > A1) {
    rRand1 = rRand2;
  }

  xRand = rRand1 * cosphi;
  yRand = rRand1 * sinphi;

  return Vector3D<Precision>(xRand, yRand, zOne);
}

//
// GetPointOnCut
//
// Auxiliary method for Get Point On Surface
//
Vector3D<Precision> UnplacedPolycone::GetPointOnCut(Precision fRMin1, Precision fRMax1, Precision fRMin2,
                                                    Precision fRMax2, Precision zOne, Precision zTwo,
                                                    Precision &totArea) const
{
  if (zOne == zTwo) {
    return GetPointOnRing(fRMin1, fRMax1, fRMin2, fRMax2, zOne);
  }
  if ((fRMin1 == fRMin2) && (fRMax1 == fRMax2)) {
    return GetPointOnTubs(fRMin1, fRMax1, zOne, zTwo, totArea);
  }
  return GetPointOnCone(fRMin1, fRMax1, fRMin2, fRMax2, zOne, zTwo, totArea);
}

//
// GetPointOnSurface
//
Vector3D<Precision> UnplacedPolycone::GetPointOnSurface() const
{
  Precision Area = 0, totArea = 0, Achose1 = 0, Achose2 = 0, phi, cosphi, sinphi, rRand;
  int i         = 0;
  int numPlanes = GetNSections();

  phi    = RNG::Instance().uniform(GetStartPhi(), GetEndPhi());
  cosphi = std::cos(phi);
  sinphi = std::sin(phi);
  std::vector<Precision> areas;
  PolyconeSection const &sec0 = GetSection(0);
  areas.push_back(
      kPi * (sec0.fSolid->GetRmax1() * sec0.fSolid->GetRmax1() - sec0.fSolid->GetRmin1() * sec0.fSolid->GetRmin1()));
  rRand = sec0.fSolid->GetRmin1() +
          ((sec0.fSolid->GetRmax1() - sec0.fSolid->GetRmin1()) * std::sqrt(RNG::Instance().uniform(0., 1.)));

  areas.push_back(
      kPi * (sec0.fSolid->GetRmax1() * sec0.fSolid->GetRmax1() - sec0.fSolid->GetRmin1() * sec0.fSolid->GetRmin1()));

  for (i = 0; i < numPlanes; i++) {
    PolyconeSection const &sec = GetSection(i);
    Area                       = (sec.fSolid->GetRmin1() + sec.fSolid->GetRmin2()) *
           std::sqrt((sec.fSolid->GetRmin1() - sec.fSolid->GetRmin2()) *
                         (sec.fSolid->GetRmin1() - sec.fSolid->GetRmin2()) +
                     4. * sec.fSolid->GetDz() * sec.fSolid->GetDz());

    Area += (sec.fSolid->GetRmax1() + sec.fSolid->GetRmax2()) *
            std::sqrt((sec.fSolid->GetRmax1() - sec.fSolid->GetRmax2()) *
                          (sec.fSolid->GetRmax1() - sec.fSolid->GetRmax2()) +
                      4. * sec.fSolid->GetDz() * sec.fSolid->GetDz());

    Area *= 0.5 * GetDeltaPhi();

    if (GetDeltaPhi() < kTwoPi) {
      Area += std::fabs(2 * sec.fSolid->GetDz()) *
              (sec.fSolid->GetRmax1() + sec.fSolid->GetRmax2() - sec.fSolid->GetRmin1() - sec.fSolid->GetRmin2());
    }

    areas.push_back(Area);
    totArea += Area;
  }
  PolyconeSection const &secn = GetSection(numPlanes - 1);
  areas.push_back(
      kPi * (secn.fSolid->GetRmax2() * secn.fSolid->GetRmax2() - secn.fSolid->GetRmin2() * secn.fSolid->GetRmin2()));

  totArea += (areas[0] + areas[numPlanes + 1]);
  Precision chose = RNG::Instance().uniform(0., totArea);

  if ((chose >= 0.) && (chose < areas[0])) {
    return Vector3D<Precision>(rRand * cosphi, rRand * sinphi, fZs[0]);
  }

  for (i = 0; i < numPlanes; i++) {
    Achose1 += areas[i];
    Achose2 = (Achose1 + areas[i + 1]);
    if (chose >= Achose1 && chose < Achose2) {
      PolyconeSection const &sec = GetSection(i);
      return GetPointOnCut(sec.fSolid->GetRmin1(), sec.fSolid->GetRmax1(), sec.fSolid->GetRmin2(),
                           sec.fSolid->GetRmax2(), fZs[i], fZs[i + 1], Area);
    }
  }

  rRand = secn.fSolid->GetRmin2() +
          ((secn.fSolid->GetRmax2() - secn.fSolid->GetRmin2()) * std::sqrt(RNG::Instance().uniform(0., 1.)));

  return Vector3D<Precision>(rRand * cosphi, rRand * sinphi, fZs[numPlanes]);
}

bool UnplacedPolycone::Normal(Vector3D<Precision> const &point, Vector3D<Precision> &norm) const
{
  bool valid = true;
  int index  = GetSectionIndex(point.z() - kTolerance);

  if (index < 0) {
    valid                 = false;
    if (index == -1) norm = Vector3D<Precision>(0., 0., -1.);
    if (index == -2) norm = Vector3D<Precision>(0., 0., 1.);
    return valid;
  }
  PolyconeSection const &sec = GetSection(index);
  valid                      = sec.fSolid->Normal(point - Vector3D<Precision>(0, 0, sec.fShift), norm);

  // if point is within tolerance of a Z-plane between 2 sections, get normal from other section too
  if (size_t(index + 1) < fSections.size() && std::abs(point.z() - fZs[index + 1]) < kTolerance) {
    PolyconeSection const &sec2 = GetSection(index + 1);
    bool valid2                 = false;
    Vector3D<Precision> norm2;
    valid2 = sec2.fSolid->Normal(point - Vector3D<Precision>(0, 0, sec2.fShift), norm2);

    if (!valid && valid2) {
      norm  = norm2;
      valid = valid2;
    }

    // if both valid && valid2 true, norm and norm2 should be added...
    if (valid && valid2) {

      // discover exiting direction by moving point a bit (it would be good to have track direction here)
      // if(sec.fSolid->Contains(point + kTolerance*10*norm - Vector3D<Precision>(0, 0, sec.fShift))){

      //}
      bool c2;
      using CI = ConeImplementation<translation::kIdentity, rotation::kIdentity, ConeTypes::UniversalCone>;
      CI::UnplacedContains<kScalar>(*sec2.fSolid,
                                    point + kTolerance * 10 * norm2 - Vector3D<Precision>(0, 0, sec2.fShift), c2);
      if (c2) {
        norm = norm2;
      } else {
        norm = norm + norm2;
        // but we might be in the interior of the polycone, and norm2=(0,0,-1) and norm1=(0,0,1) --> norm=(0,0,0)
        // quick fix:  set a normal pointing to the input point, but setting its z=0 (radial)
        if (norm.Mag2() < kTolerance) norm = Vector3D<Precision>(point.x(), point.y(), 0);
      }
    }
  }
  if (valid) norm /= norm.Mag();
  return valid;
}

Precision UnplacedPolycone::SurfaceArea() const
{
  Precision Area = 0, totArea = 0;
  int i                  = 0;
  int numPlanes          = GetNSections();
  Precision fSurfaceArea = 0;

  Vector<Precision> areas; // (numPlanes+1);

  PolyconeSection const &sec0 = GetSection(0);
  areas.push_back(
      kPi * (sec0.fSolid->GetRmax1() * sec0.fSolid->GetRmax1() - sec0.fSolid->GetRmin1() * sec0.fSolid->GetRmin1()));
  for (i = 0; i < numPlanes; i++) {
    PolyconeSection const &sec = GetSection(i);
    Area                       = (sec.fSolid->GetRmin1() + sec.fSolid->GetRmin2()) *
           std::sqrt((sec.fSolid->GetRmin1() - sec.fSolid->GetRmin2()) *
                         (sec.fSolid->GetRmin1() - sec.fSolid->GetRmin2()) +
                     4. * sec.fSolid->GetDz() * sec.fSolid->GetDz());

    Area += (sec.fSolid->GetRmax1() + sec.fSolid->GetRmax2()) *
            std::sqrt((sec.fSolid->GetRmax1() - sec.fSolid->GetRmax2()) *
                          (sec.fSolid->GetRmax1() - sec.fSolid->GetRmax2()) +
                      4. * sec.fSolid->GetDz() * sec.fSolid->GetDz());

    Area *= 0.5 * GetDeltaPhi();

    if (GetDeltaPhi() < kTwoPi) {
      Area += std::fabs(2 * sec.fSolid->GetDz()) *
              (sec.fSolid->GetRmax1() + sec.fSolid->GetRmax2() - sec.fSolid->GetRmin1() - sec.fSolid->GetRmin2());
    }
    areas.push_back(Area);
    totArea += Area;
  }
  PolyconeSection const &secn = GetSection(numPlanes - 1);
  areas.push_back(
      kPi * (secn.fSolid->GetRmax2() * secn.fSolid->GetRmax2() - secn.fSolid->GetRmin2() * secn.fSolid->GetRmin2()));

  totArea += (areas[0] + areas[numPlanes + 1]);
  fSurfaceArea = totArea;

  return fSurfaceArea;
}

void UnplacedPolycone::Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const
{

  int i          = 0;
  Precision maxR = 0;

  for (i = 0; i < GetNSections(); i++) {
    PolyconeSection const &sec              = GetSection(i);
    if (maxR < sec.fSolid->GetRmax1()) maxR = sec.fSolid->GetRmax1();
    if (maxR < sec.fSolid->GetRmax2()) maxR = sec.fSolid->GetRmax2();
  }

  aMin.x() = -maxR;
  aMin.y() = -maxR;
  aMin.z() = fZs[0];
  aMax.x() = maxR;
  aMax.y() = maxR;
  aMax.z() = fZs[GetNSections()];
}
#endif // !VECGEOM_NVCC

bool UnplacedPolycone::CheckContinuityInRmax(const Vector<Precision> &rOuter)
{

  bool continuous  = true;
  unsigned int len = rOuter.size();
  if (len > 2) {
    for (unsigned int j = 1; j < len;) {
      if (j != (len - 1)) continuous &= (rOuter[j] == rOuter[j + 1]);
      j = j + 2;
    }
  }
  return continuous;
}

/*Improvising the convexity detection algo.
* Much Better implementation, currently passes all
* the existing test and newly added test also.
*
* Algo: Instead of putting complex checks on slope, especially the
* conditions where z[i]==z[i+1], convert the polycone to reduced
* polycone, so that it is not having those section, which are
* basically not required for convexity calculation.
*
* Once reduced polycone is calculated, pass it to do the check
* for Slope continuity and Rmax continuity
*
* return the output of AND operation of SlopeContinuity and RmaxContinuity
*
* Much cleaner implementation, and easy to debug.
*/

bool UnplacedPolycone::CheckContinuity(const double rOuter[], const double rInner[], const double zPlane[],
                                       Vector<Precision> &newROuter, Vector<Precision> &newRInner,
                                       Vector<Precision> &newZPlane)
{

  Vector<Precision> rOut, rIn;
  Vector<Precision> zPl;
  rOut.push_back(rOuter[0]);
  rIn.push_back(rInner[0]);
  zPl.push_back(zPlane[0]);
  for (unsigned int j = 1; j < fNz; j++) {

    if (j == fNz - 1) {
      rOut.push_back(rOuter[j]);
      rIn.push_back(rInner[j]);
      zPl.push_back(zPlane[j]);
    } else {
      if ((zPlane[j] != zPlane[j + 1]) || (rOuter[j] != rOuter[j + 1])) {
        rOut.push_back(rOuter[j]);
        rOut.push_back(rOuter[j]);

        zPl.push_back(zPlane[j]);
        zPl.push_back(zPlane[j]);

        rIn.push_back(rInner[j]);
        rIn.push_back(rInner[j]);

      } else {
        rOut.push_back(rOuter[j]);
        zPl.push_back(zPlane[j]);
        rIn.push_back(rInner[j]);
      }
    }
  }

  if (rOut.size() % 2 != 0) {
    // fNz is odd, the adding of the last item did not happen in the loop.
    rOut.push_back(rOut[rOut.size() - 1]);
    rIn.push_back(rIn[rIn.size() - 1]);
    zPl.push_back(zPl[zPl.size() - 1]);
  }

  /* Creating a new temporary Reduced polycone with desired data elements,
  *  which makes sure that denominator will never be zero (hence avoiding FPE(division by zero)),
  *  while calculating slope.
  *
  *  This will be the minimum polycone,i.e. no extra section which
  *  affect its shape
  */

  for (size_t j = 0; j < rOut.size();) {

    if (zPl[j] != zPl[j + 1]) {
      newZPlane.push_back(zPl[j]);
      newZPlane.push_back(zPl[j + 1]);
      newROuter.push_back(rOut[j]);
      newROuter.push_back(rOut[j + 1]);
      newRInner.push_back(rIn[j]);
      newRInner.push_back(rIn[j + 1]);
    }

    j = j + 2;
  }
  // Minimum polycone construction over

  // Checking Slope continuity and Rmax Continuity
  bool contRmax  = CheckContinuityInRmax(newROuter);
  bool contSlope = CheckContinuityInSlope(newROuter, newZPlane);

  // If both are true then the polycone can be convex
  // but still final convexity depends on Inner Radius also.
  return (contRmax && contSlope);
}

/* Cleaner CheckContinuityInSlope.
 * Because of new design, it will not get the case of FPE exception
 * (division by zero)
 */
bool UnplacedPolycone::CheckContinuityInSlope(const Vector<Precision> &rOuter, const Vector<Precision> &zPlane)
{

  bool continuous      = true;
  Precision startSlope = kInfLength;

  // Doing the actual slope calculation here, and checking continuity,
  for (size_t j = 0; j < rOuter.size(); j = j + 2) {
    Precision currentSlope = (rOuter[j + 1] - rOuter[j]) / (zPlane[j + 1] - zPlane[j]);
    continuous &= (currentSlope <= startSlope);
    startSlope = currentSlope;
  }
  return continuous;
}

VECGEOM_CUDA_HEADER_BOTH
void UnplacedPolycone::DetectConvexity()
{
  // Default safe convexity value
  fGlobalConvexity = false;

  if (fConvexityPossible) {
    if (fEqualRmax && (fDeltaPhi <= kPi || fDeltaPhi == kTwoPi))
      // In this case, Polycone become solid Cylinder, No need to check anything else, 100% convex
      fGlobalConvexity = true;
    else {
      if (fDeltaPhi <= kPi || fDeltaPhi == kTwoPi) {
        fGlobalConvexity = fContinuityOverAll;
      }
    }
  }

  // return convexity;
}
} // End impl namespace

#ifdef VECGEOM_NVCC

namespace cxx {

template size_t DevicePtr<cuda::UnplacedPolycone>::SizeOf();
template void DevicePtr<cuda::UnplacedPolycone>::Construct(Precision, Precision, int, Precision *, Precision *,
                                                           Precision *) const;

} // End cxx namespace

#endif

} // end global namespace
