/// \file Tile.cpp
/// \author Mihaela Gheata (mihaela.gheata@cern.ch)

#include <ostream>
#include "volumes/Tile.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

#ifndef VECCORE_CUDA
std::ostream &operator<<(std::ostream &os, TriangleFacet<double> const &facet)
{
#ifndef VECCORE_ENABLE_UMESIMD
  os << " triangle facet:\n";
  os << "    vertices: {" << facet.fVertices[0] << ", " << facet.fVertices[1] << ", " << facet.fVertices[2] << "}\n";
  os << "    indices:  {" << facet.fIndices << "}\n";
  os << "    normal: {" << facet.fNormal << "}\n";
  os << "    distance: " << facet.fDistance << "}";
#endif
  return os;
}

std::ostream &operator<<(std::ostream &os, QuadrilateralFacet<double> const &facet)
{
#ifndef VECCORE_ENABLE_UMESIMD
  os << " quadrilateral facet:\n";
  os << "    vertices: {" << facet.fVertices[0] << ", " << facet.fVertices[1] << ", " << facet.fVertices[2] << ", "
     << facet.fVertices[3] << "}\n";
  os << "    indices:  {" << facet.fIndices[0] << ", " << facet.fIndices[1] << ", " << facet.fIndices[2] << ", "
     << facet.fIndices[3] << "}\n";
  os << "    normal: {" << facet.fNormal << "}\n";
  os << "    distance: " << facet.fDistance << "}";
#endif
  return os;
}
#endif
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom