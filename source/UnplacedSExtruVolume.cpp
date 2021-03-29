#include "VecGeom/volumes/UnplacedSExtruVolume.h"
#include "VecGeom/volumes/PlacedSExtru.h"
#include "VecGeom/base/RNG.h"
#include <stdio.h>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

VECCORE_ATT_HOST_DEVICE
void UnplacedSExtruVolume::Print() const
{
  printf("UnplacedSExtruVolume\n");
}

void UnplacedSExtruVolume::Print(std::ostream &os) const
{
  os << "UnplacedSExtruVolume";
}

#ifndef VECCORE_CUDA
SolidMesh *UnplacedSExtruVolume::CreateMesh3D(Transformation3D const &trans, size_t nSegments) const
{
  const PlanarPolygon &polygon = fPolyShell.GetPolygon();
  const Precision lowerZ       = fPolyShell.GetLowerZ();
  const Precision upperZ       = fPolyShell.GetUpperZ();

  const Precision *verticesX = polygon.GetVertices().x();
  const Precision *verticesY = polygon.GetVertices().y();
  size_t nVertices           = polygon.GetVertices().size();

  size_t nMeshVertices = nVertices * 2;
  size_t nMeshPolygons = 2 + nVertices;

  SolidMesh *sm = new SolidMesh();
  sm->ResetMesh(nMeshVertices, nMeshPolygons);

  typedef Vector3D<Precision> Vec_t;
  Vec_t *const vertices = new Vec_t[nMeshVertices];
  for (size_t i = 0; i < nVertices; i++) {
    vertices[2 * i]     = Vec_t(verticesX[i], verticesY[i], lowerZ); // even lower vertices
    vertices[2 * i + 1] = Vec_t(verticesX[i], verticesY[i], upperZ); // odd upper vertices
  }

  sm->SetVertices(vertices, nMeshVertices);
  delete[] vertices;

  sm->TransformVertices(trans);

  std::vector<size_t> indices;
  indices.reserve(nVertices);

  // lower surface
  for (size_t i = 0; i < 2 * nVertices; i += 2) {
    indices.push_back(i);
  }
  sm->AddPolygon(nVertices, indices, polygon.IsConvex());

  indices.clear();
  // upper surface
  for (size_t i = 2 * nVertices; i > 0; i -= 2) {
    indices.push_back(i - 1);
  }
  sm->AddPolygon(nVertices, indices, polygon.IsConvex());

  // lateral
  for (size_t i = 0, j = 0; i < (nVertices - 1); i++, j += 2) {
    sm->AddPolygon(4, {j, j + 1, j + 3, j + 2}, true);
  }
  sm->AddPolygon(4, {0, nMeshVertices - 2, nMeshVertices - 1, 1}, true);

  // sm->InitSExtruVolume(nMeshVertices, nMeshPolygons, polygon.IsConvex());

  return sm;
}
#endif

VECCORE_ATT_DEVICE
VPlacedVolume *UnplacedSExtruVolume::PlaceVolume(LogicalVolume const *const logical_volume, Transformation3D const *const transformation,
#ifdef VECCORE_CUDA
                                                 const int id, const int copy_no, const int child_id,
#endif
                                                 VPlacedVolume *const placement) const
{
  if (placement) {
    return new (placement) PlacedSExtru(logical_volume, transformation
#ifdef VECCORE_CUDA
                                        , id, copy_no, child_id
#endif
    );
    return placement;
  }
  return new PlacedSExtru(logical_volume, transformation
#ifdef VECCORE_CUDA
                          , id, copy_no, child_id
#endif
  );
}

#ifdef VECGEOM_CUDA_INTERFACE

DevicePtr<cuda::VUnplacedVolume> UnplacedSExtruVolume::CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const gpu_ptr) const
{
  auto &vertices         = fPolyShell.fPolygon.GetVertices();
  Precision const *x_cpu = vertices.x();
  Precision const *y_cpu = vertices.y();
  const auto Nvert       = vertices.size();

  // copying the arrays needed for the constructor
  Precision *x_gpu_ptr = AllocateOnGpu<Precision>(Nvert * sizeof(Precision));
  Precision *y_gpu_ptr = AllocateOnGpu<Precision>(Nvert * sizeof(Precision));
  vecgeom::CopyToGpu(x_cpu, x_gpu_ptr, sizeof(Precision) * Nvert);
  vecgeom::CopyToGpu(y_cpu, y_gpu_ptr, sizeof(Precision) * Nvert);

  DevicePtr<cuda::VUnplacedVolume> gpusextru = CopyToGpuImpl<UnplacedSExtruVolume>(
      gpu_ptr, (int)Nvert, x_gpu_ptr, y_gpu_ptr, fPolyShell.fLowerZ, fPolyShell.fUpperZ);

  // remove temporary space from GPU
  FreeFromGpu(x_gpu_ptr);
  FreeFromGpu(y_gpu_ptr);

  return gpusextru;
}

DevicePtr<cuda::VUnplacedVolume> UnplacedSExtruVolume::CopyToGpu() const
{
  return CopyToGpuImpl<UnplacedSExtruVolume>();
}

#endif // VECGEOM_CUDA_INTERFACE

} // namespace VECGEOM_IMPL_NAMESPACE

#ifdef VECCORE_CUDA

namespace cxx {

template size_t DevicePtr<cuda::UnplacedSExtruVolume>::SizeOf();
template void DevicePtr<cuda::UnplacedSExtruVolume>::Construct(int, Precision *, Precision *, Precision,
                                                               Precision) const;

} // namespace cxx

#endif

} // namespace vecgeom
