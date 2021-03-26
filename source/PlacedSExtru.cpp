#include "VecGeom/volumes/PlacedSExtru.h"
#include "VecGeom/base/SOA3D.h"

#ifdef VECGEOM_ROOT
#include "TGeoXtru.h"
#endif

#ifdef VECGEOM_GEANT4
#include "G4ExtrudedSolid.hh"
#include "G4TwoVector.hh"
#include <vector>
#endif

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

VECCORE_ATT_HOST_DEVICE
void PlacedSExtru::PrintType() const
{
  printf("PlacedSExtru");
}

void PlacedSExtru::PrintType(std::ostream &s) const
{
  s << "PlacedSExtru";
}

// Comparison specific
#ifndef VECCORE_CUDA
#ifdef VECGEOM_ROOT
TGeoShape const *PlacedSExtru::ConvertToRoot() const
{
  TGeoXtru *s = new TGeoXtru(2);

  // get vertices and make array to construct the ROOT volume
  auto &vertices = GetUnplacedStruct()->GetPolygon().GetVertices();
  const auto N   = GetUnplacedStruct()->GetPolygon().GetNVertices();
  double *vertdx = new double[N];
  double *vertdy = new double[N];
  for (unsigned int i = 0; i < N; ++i) {
    vertdx[i] = vertices.x(i);
    vertdy[i] = vertices.y(i);
  }
  s->DefinePolygon(N, vertdx, vertdy);
  s->DefineSection(0, GetUnplacedStruct()->GetLowerZ(), 0., 0., 1.);
  s->DefineSection(1, GetUnplacedStruct()->GetUpperZ(), 0., 0., 1.);
  delete[] vertdx;
  delete[] vertdy;
  return s;
};
#endif
#ifdef VECGEOM_GEANT4
G4VSolid const *PlacedSExtru::ConvertToGeant4() const
{
  std::vector<G4TwoVector> G4vertices;
  using ZSection = G4ExtrudedSolid::ZSection;
  std::vector<ZSection> zsections;
  auto &vertices = GetUnplacedStruct()->GetPolygon().GetVertices();
  for (size_t i = 0; i < vertices.size(); ++i) {
    G4vertices.push_back(G4TwoVector(vertices.x(i), vertices.y(i)));
  }
  zsections.push_back(ZSection(GetUnplacedStruct()->GetLowerZ(), 0., 1.));
  zsections.push_back(ZSection(GetUnplacedStruct()->GetUpperZ(), 0., 1.));
  G4String s("g4extru");
  return new G4ExtrudedSolid(s, G4vertices, zsections);
};
#endif
#endif // VECCORE_CUDA

} // namespace VECGEOM_IMPL_NAMESPACE

#ifdef VECCORE_CUDA

VECGEOM_DEVICE_INST_PLACED_VOLUME(PlacedSExtru)

#endif // VECCORE_CUDA

} // End namespace vecgeom
