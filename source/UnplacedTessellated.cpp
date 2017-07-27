/// @file UnplacedTessellated.cpp
/// @author Mihaela Gheata (mihaela.gheata@cern.ch)

#include "volumes/UnplacedTessellated.h"
#include "volumes/SpecializedTessellated.h"
#include "volumes/utilities/GenerationUtilities.h"
#include "base/RNG.h"

#include "management/VolumeFactory.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

void UnplacedTessellated::Print() const
{
  printf("UnplacedTessellated {%lu facets}", fTessellated.fFacets.size());
}

void UnplacedTessellated::Print(std::ostream &os) const
{
  os << "UnplacedTessellated {" << fTessellated.fFacets.size() << " facets " << std::endl;
}

Precision UnplacedTessellated::Capacity() const
{
  if (fTessellated.fCubicVolume != 0.) return fTessellated.fCubicVolume;

  // For explanation of the following algorithm see:
  // https://en.wikipedia.org/wiki/Polyhedron#Volume
  // http://wwwf.imperial.ac.uk/~rn/centroid.pdf

  int size = fTessellated.fFacets.size();
  for (int i = 0; i < size; ++i) {
    TriangleFacet<double> &facet = *fTessellated.fFacets[i];
    double area                  = facet.fSurfaceArea;
    fTessellated.fCubicVolume += area * (facet.fVertices[0].Dot(facet.fNormal));
  }
  fTessellated.fCubicVolume /= 3.;
  return fTessellated.fCubicVolume;
}

Precision UnplacedTessellated::SurfaceArea() const
{
  if (fTessellated.fSurfaceArea != 0.) return fTessellated.fSurfaceArea;

  int size = fTessellated.fFacets.size();
  for (int i = 0; i < size; ++i) {
    TriangleFacet<double> *facet = fTessellated.fFacets[i];
    fTessellated.fSurfaceArea += facet->fSurfaceArea;
  }
  return fTessellated.fSurfaceArea;
}

int UnplacedTessellated::ChooseSurface() const
{
  int choice       = 0; // 0 = zm, 1 = zp, 2 = ym, 3 = yp, 4 = xm, 5 = xp
  Precision Stotal = SurfaceArea();

  // random value to choose surface to place the point
  Precision rand = RNG::Instance().uniform() * Stotal;

  while (rand > fTessellated.fFacets[choice]->fSurfaceArea)
    rand -= fTessellated.fFacets[choice]->fSurfaceArea, choice++;

  return choice;
}

Vector3D<Precision> UnplacedTessellated::SamplePointOnSurface() const
{
  int surface    = ChooseSurface();
  double alpha   = RNG::Instance().uniform(0., 1.);
  double beta    = RNG::Instance().uniform(0., 1.);
  double lambda1 = alpha * beta;
  double lambda0 = alpha - lambda1;
  auto facet     = fTessellated.fFacets[surface];
  return (facet->fVertices[0] + lambda0 * (facet->fVertices[1] - facet->fVertices[0]) +
          lambda1 * (facet->fVertices[2] - facet->fVertices[1]));
}

bool UnplacedTessellated::Normal(Vector3D<Precision> const &point, Vector3D<Precision> &norm) const
{
  // Redirect to normal implementation
  bool valid = false;
  norm       = TessellatedImplementation::NormalKernel<Precision>(fTessellated, point, valid);
  return valid;
}

#ifdef VECCORE_CUDA
template <TranslationCode transCodeT, RotationCode rotCodeT>
VECCORE_ATT_DEVICE
VPlacedVolume *UnplacedTessellated::Create(LogicalVolume const *const logical_volume,
                                           Transformation3D const *const transformation, const int id,
                                           VPlacedVolume *const placement)
{
  if (placement) {
    new (placement) SpecializedTessellated<transCodeT, rotCodeT>(logical_volume, transformation, id);
    return placement;
  }
  return new SpecializedTessellated<transCodeT, rotCodeT>(logical_volume, transformation, id);
}
#else
template <TranslationCode transCodeT, RotationCode rotCodeT>
VPlacedVolume *UnplacedTessellated::Create(LogicalVolume const *const logical_volume,
                                           Transformation3D const *const transformation, VPlacedVolume *const placement)
{
  if (placement) {
    new (placement) SpecializedTessellated<transCodeT, rotCodeT>(logical_volume, transformation);
    return placement;
  }
  return new SpecializedTessellated<transCodeT, rotCodeT>(logical_volume, transformation);
}
#endif

VECCORE_ATT_DEVICE
VPlacedVolume *UnplacedTessellated::SpecializedVolume(LogicalVolume const *const volume,
                                                      Transformation3D const *const transformation,
                                                      const TranslationCode trans_code, const RotationCode rot_code,
#ifdef VECCORE_CUDA
                                                      const int id,
#endif
                                                      VPlacedVolume *const placement) const
{

  return VolumeFactory::CreateByTransformation<UnplacedTessellated>(volume, transformation, trans_code, rot_code,
#ifdef VECCORE_CUDA
                                                                    id,
#endif
                                                                    placement);
}

#if defined(VECGEOM_USOLIDS)
std::ostream &UnplacedTessellated::StreamInfo(std::ostream &os) const
{
  int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "     *** Dump for solid - " << GetEntityType() << " ***\n"
     << "     ===================================================\n"
     << " Solid type: Trd\n"
     << " Parameters: \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);
  return os;
}
#endif

#ifdef VECGEOM_CUDA_INTERFACE

DevicePtr<cuda::VUnplacedVolume> UnplacedTessellated::CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const in_gpu_ptr) const
{
  return CopyToGpuImpl<UnplacedTessellated>(in_gpu_ptr);
}

DevicePtr<cuda::VUnplacedVolume> UnplacedTessellated::CopyToGpu() const
{
  return CopyToGpuImpl<UnplacedTessellated>();
}

#endif // VECGEOM_CUDA_INTERFACE

} // End impl namespace

#ifdef VECCORE_CUDA

namespace cxx {

template size_t DevicePtr<cuda::UnplacedTessellated>::SizeOf();
template void DevicePtr<cuda::UnplacedTessellated>::Construct() const;

} // End cxx namespace

#endif

} // End global namespace
