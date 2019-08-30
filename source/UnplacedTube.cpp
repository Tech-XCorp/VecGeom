/// \file UnplacedTube.cpp
/// \author Georgios Bitzes (georgios.bitzes@cern.ch)

#include "volumes/UnplacedTube.h"
#include "volumes/SpecializedTube.h"
#include "base/RNG.h"
#ifndef VECCORE_CUDA
#include <cmath>
#include <iostream>
#endif

#include "volumes/utilities/GenerationUtilities.h"
#include "management/VolumeFactory.h"
#include "volumes/UnplacedEllipticalTube.h"

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

#ifndef VECCORE_CUDA
SolidMesh *UnplacedTube::CreateMesh3D(Transformation3D const &trans, size_t nSegments) const
{
  typedef Vector3D<double> Vec_t;

  SolidMesh *sm = new SolidMesh();

  sm->ResetMesh(4*(nSegments + 1), 4*nSegments + 2);

  // fill vertex array
  Vec_t *vertices = new Vec_t[4*(nSegments + 1)];
  double angle, step, rcos, rsin;



  size_t idx  = 0;
  size_t idx1 = (nSegments + 1);
  size_t idx2 = 2 * (nSegments + 1);
  size_t idx3 = 3 * (nSegments + 1);

  double phi      = sphi();
  double phi_step = dphi() / nSegments;

  double x, y;
  for (size_t i = 0; i <= nSegments; i++, phi += phi_step) {
    x               = rmax() * std::cos(phi);
    y               = rmax() * std::sin(phi);
    vertices[idx++] = Vec_t(x, y, z()); // top outer
    vertices[idx1++] = Vec_t(x, y, -z()); // bottom outer
    x                = rmin() * std::cos(phi);
    y                = rmin() * std::sin(phi);
    vertices[idx2++] = Vec_t(x, y, z()); // top inner
    vertices[idx3++] = Vec_t(x, y, -z()); // bottom inner
  }



  sm->SetVertices(vertices, 4*(nSegments + 1));
  delete[] vertices;
  sm->TransformVertices(trans);


  for (size_t i = 0, j = nSegments + 1; i < nSegments; i++, j++) {
    sm->AddPolygon(4, {i, j, j + 1, i + 1}, true); // OUTER
  }

  for (size_t i = 0, j = 2 * (nSegments + 1), k = j + nSegments + 1; i < nSegments; i++, j++, k++) {
    sm->AddPolygon(4, {j, j + 1, k + 1, k}, true); // inner
  }

  for (size_t i = 0, j = (nSegments + 1), k = j + 2 * (nSegments + 1); i < nSegments; i++, j++, k++) {
    sm->AddPolygon(4, {j, k, k + 1, j + 1}, true); // lower
  }

  for (size_t i = 0, j = 0, k = j + 2 * (nSegments + 1); i < nSegments; i++, j++, k++) {
    sm->AddPolygon(4, {j, j + 1, k + 1, k}, true); // Upper
  }

  if (dphi() != kTwoPi) {
    sm->AddPolygon(4, {0, 2 * (nSegments + 1), 3 * (nSegments + 1), nSegments + 1}, true);
    sm->AddPolygon(
        4, {nSegments, nSegments + nSegments + 1, nSegments + 3 * (nSegments + 1), nSegments + 2 * (nSegments + 1)},
        true);
  }

  // case 1: use elliptical tube with circle base


  /*
  if (rmin() == 0. && dphi() == 2 * kPi) {
    VUnplacedVolume *tube = new UnplacedEllipticalTube(rmax(), rmax(), z());
    sm                    = tube->CreateMesh3D(trans, nFaces);
    delete tube;
  }
  // case 2: partial cylinder (not hollow)
  else if (rmin() == 0. && dphi() != 2 * kPi) {

    size_t nMeshVertices = (nFaces - 3) * 2 + 2;
    sm                   = new SolidMesh();
    sm->ResetMesh(nMeshVertices, nFaces);

    // fill vertex array
    Vec_t *vertices = new Vec_t[nMeshVertices];
    double angle, step, rcos, rsin;
    step  = dphi() / (nMeshVertices - 4);
    angle = sphi();
    for (size_t i = 0; i < nMeshVertices - 2; i += 2, angle += 2 * step) {
      rcos            = rmax() * std::cos(angle);
      rsin            = rmax() * std::sin(angle);
      vertices[i]     = Vec_t(rcos, rsin, -z()); // even indices hold lower vertices
      vertices[i + 1] = Vec_t(rcos, rsin, z());  // odd indices hold upper vertices
    }
    vertices[nMeshVertices - 2] = Vec_t(0, 0, -z());
    vertices[nMeshVertices - 1] = Vec_t(0, 0, z());
    sm->SetVertices(vertices, nMeshVertices);
    delete[] vertices;
    sm->TransformVertices(trans);

    // fill polygon array

    Utils3D::vector_t<size_t> indices;
    indices.reserve(nFaces - 2);

    // lower surface
    for (size_t i = nFaces - 3; i > 0; i--) {
      indices.push_back(2 * (i - 1));
    }
    indices.push_back(nMeshVertices - 2);
    sm->AddPolygon(nFaces - 2, indices, true);

    // upper surface
    indices.clear();
    for (size_t i = 0; i < nMeshVertices; i += 2) {
      indices.push_back(i + 1);
    }
    indices.push_back(nMeshVertices - 1);
    sm->AddPolygon(nFaces - 2, indices, true);

    // surfaces due to cylinder not being full

    sm->AddPolygon(4, {nMeshVertices - 1, nMeshVertices - 3, nMeshVertices - 4, nMeshVertices - 2}, true);
    sm->AddPolygon(4, {0, 1, nMeshVertices - 1, nMeshVertices - 2}, true);

    // outer surface

    for (size_t i = 0; i < nFaces - 4; i++) {
      sm->AddPolygon(4, {2 * i, 2 * i + 2, 2 * i + 3, 2 * i + 1}, true);
    }
    sm->InitPolygons();

  } // case 3, 4: rmin > 0
  else {

    size_t nMeshPerRegion = std::ceil(nFaces / 4.);
    size_t nMeshVertices  = nMeshPerRegion * 4 + 4;
    sm                    = new SolidMesh();
    sm->ResetMesh(nMeshVertices, 4 * nMeshPerRegion + 2);

    // fill vertex array
    Vec_t *vertices = new Vec_t[nMeshVertices];
    double angle, step, rmincos, rminsin, rmaxcos, rmaxsin;
    step  = dphi() / (nMeshVertices - 4);
    angle = sphi();
    for (size_t i = 0; i < nMeshVertices; i += 4, angle += 4 * step) {
      rmincos         = rmin() * std::cos(angle);
      rminsin         = rmin() * std::sin(angle);
      rmaxcos         = rmax() * std::cos(angle);
      rmaxsin         = rmax() * std::sin(angle);
      vertices[i]     = Vec_t(rmincos, rminsin, -z()); // rmin low
      vertices[i + 1] = Vec_t(rmincos, rminsin, z());  // rmin high
      vertices[i + 2] = Vec_t(rmaxcos, rmaxsin, -z()); // rmax low
      vertices[i + 3] = Vec_t(rmaxcos, rmaxsin, z());  // rmax high
    }
    sm->SetVertices(vertices, nMeshVertices);
    delete[] vertices;
    sm->TransformVertices(trans);

    // fill polygons

    // lower surface
    for (size_t i = 0; i < nMeshPerRegion; i++) {
      sm->AddPolygon(4, {4 * i, 4 * (i + 1), 4 * (i + 1) + 2, 4 * i + 2}, true);
    }

    // upper surface
    for (size_t i = 0; i < nMeshPerRegion; i++) {
      sm->AddPolygon(4, {4 * i + 3, 4 * (i + 1) + 3, 4 * (i + 1) + 1, 4 * i + 1}, true);
    }

    //  inner surface
    for (size_t i = 0; i < nMeshPerRegion; i++) {
      sm->AddPolygon(4, {4 * (i + 1), 4 * i, 4 * i + 1, 4 * (i + 1) + 1}, true);
    }

    // outer surface
    for (size_t i = 0; i < nMeshPerRegion; i++) {
      sm->AddPolygon(4, {4 * (i + 1) + 3, 4 * i + 3, 4 * i + 2, 4 * (i + 1) + 2}, true);
    }

    // surfaces due to cylinder not being full
    if (dphi() != 2 * kPi) {
      sm->AddPolygon(4, {2, 3, 1, 0}, true);
      sm->AddPolygon(4, {nMeshVertices - 4, nMeshVertices - 3, nMeshVertices - 1, nMeshVertices - 2}, true);
    }
    sm->InitPolygons();
  }
  */

  return sm;
}
#endif

template <>
UnplacedTube *Maker<UnplacedTube>::MakeInstance(const Precision &rmin, const Precision &rmax, const Precision &z,
                                                const Precision &sphi, const Precision &dphi)
{
#ifndef VECGEOM_NO_SPECIALIZATION
  if (rmin <= 0) {
    if (dphi >= 2 * M_PI) return new SUnplacedTube<TubeTypes::NonHollowTube>(rmin, rmax, z, sphi, dphi);
    if (dphi == M_PI) return new SUnplacedTube<TubeTypes::NonHollowTubeWithPiSector>(rmin, rmax, z, sphi, dphi);
    if (dphi < M_PI)
      return new SUnplacedTube<TubeTypes::NonHollowTubeWithSmallerThanPiSector>(rmin, rmax, z, sphi, dphi);
    if (dphi > M_PI)
      return new SUnplacedTube<TubeTypes::NonHollowTubeWithBiggerThanPiSector>(rmin, rmax, z, sphi, dphi);
  } else if (rmin > 0) {
    if (dphi >= 2 * M_PI) return new SUnplacedTube<TubeTypes::HollowTube>(rmin, rmax, z, sphi, dphi);
    if (dphi == M_PI) return new SUnplacedTube<TubeTypes::HollowTubeWithPiSector>(rmin, rmax, z, sphi, dphi);
    if (dphi < M_PI) return new SUnplacedTube<TubeTypes::HollowTubeWithSmallerThanPiSector>(rmin, rmax, z, sphi, dphi);
    if (dphi > M_PI) return new SUnplacedTube<TubeTypes::HollowTubeWithBiggerThanPiSector>(rmin, rmax, z, sphi, dphi);
  }
  // just here to trigger symbol creation (because it might be used explicitly elsewhere)
  return new SUnplacedTube<TubeTypes::UniversalTube>(rmin, rmax, z, sphi, dphi);
#else
  // if nothing matches return the most general case
  // in principle this should never happen
  return new SUnplacedTube<TubeTypes::UniversalTube>(rmin, rmax, z, sphi, dphi);
#endif
}

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
  if (phi0out) aMax.x() = Max(Cx, Dx);
  if (phi90out) aMax.y() = Max(Cy, Dy);
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

} // namespace VECGEOM_IMPL_NAMESPACE

#ifdef VECCORE_CUDA

namespace cxx {

template size_t DevicePtr<cuda::SUnplacedTube<cuda::TubeTypes::UniversalTube>>::SizeOf();
template void DevicePtr<cuda::SUnplacedTube<cuda::TubeTypes::UniversalTube>>::Construct(
    const Precision rmin, const Precision rmax, const Precision z, const Precision sphi, const Precision dphi) const;

} // namespace cxx

#endif

} // namespace vecgeom
