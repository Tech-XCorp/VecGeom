/// \file PlacedVolume.cpp
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#include "volumes/PlacedVolume.h"
#include "base/Vector3D.h"
#include "base/RNG.h"
#include "management/GeoManager.h"

#include <stdio.h>
#include <cassert>

#ifdef VECGEOM_USOLIDS
#include "volumes/USolidsInterfaceHelper.h"
#endif

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

int VPlacedVolume::g_id_count = 0;

#ifndef VECGEOM_NVCC
VPlacedVolume::VPlacedVolume(char const *const label,
                             LogicalVolume const *const logical_volume,
                             Transformation3D const *const transformation,
                             PlacedBox const *const bounding_box)
  :
#ifdef VECGEOM_USOLIDS
    USolidsInterfaceHelper(label),
#endif
     id_(), label_(NULL), logical_volume_(logical_volume), transformation_(transformation),
      bounding_box_(bounding_box) {
  id_ = g_id_count++;
  GeoManager::Instance().RegisterPlacedVolume(this);
  label_ = new std::string(label);
}

VPlacedVolume::VPlacedVolume(VPlacedVolume const & other) : id_(), label_(NULL), logical_volume_(), transformation_(),
    bounding_box_() {
  assert( 0 && "COPY CONSTRUCTOR FOR PlacedVolumes NOT IMPLEMENTED");
}

VPlacedVolume * VPlacedVolume::operator=( VPlacedVolume const & other )
{
  printf("ASSIGNMENT OPERATOR FOR VPlacedVolumes NOT IMPLEMENTED");
  return NULL;
}
#endif

VECGEOM_CUDA_HEADER_BOTH
VPlacedVolume::~VPlacedVolume() {
#ifndef VECGEOM_NVCC_DEVICE
  delete label_;
#endif
}

VECGEOM_CUDA_HEADER_BOTH
void VPlacedVolume::Print(const int indent) const {
  for (int i = 0; i < indent; ++i) printf("  ");
  PrintType();
  printf(" [%i]", id_);
#ifndef VECGEOM_NVCC
  if (label_->size()) {
    printf(" \"%s\"", label_->c_str());
  }
#endif
  printf(": \n");
  for (int i = 0; i <= indent; ++i) printf("  ");
  transformation_->Print();
  printf("\n");
  logical_volume_->Print(indent+1);
}

VECGEOM_CUDA_HEADER_BOTH
void VPlacedVolume::PrintContent(const int indent) const {
  Print(indent);
  if (daughters().size() > 0) {
    printf(":");
    for (VPlacedVolume const **vol = daughters().begin(),
         **volEnd = daughters().end(); vol != volEnd; ++vol) {
      printf("\n");
      (*vol)->PrintContent(indent+3);
    }
  }
}

VECGEOM_CUDA_HEADER_HOST
std::ostream& operator<<(std::ostream& os, VPlacedVolume const &vol) {
  os << "(" << (*vol.unplaced_volume()) << ", " << (*vol.GetTransformation())
     << ")";
  return os;
}

// implement a default function for surface area
// based on the method of G4
Precision VPlacedVolume::SurfaceArea() {
    std::cout << "WARNING : Sampling SurfaceArea called \n";
 int nStat = 100000;
 double ell = -1.;
 Vector3D<Precision> p;
 Vector3D<Precision> minCorner;
 Vector3D<Precision> maxCorner;
 Vector3D<Precision> delta;

 // min max extents of pSolid along X,Y,Z
 this->Extent(minCorner,maxCorner);

 // limits
 delta = maxCorner - minCorner;

 if(ell<=0.)          // Automatic definition of skin thickness
 {
   Precision minval = delta.x();
   if(delta.y() < delta.x()) { minval= delta.y(); }
   if(delta.z() < minval) { minval= delta.z(); }
   ell=.01*minval;
 }

 Precision dd=2*ell;
 minCorner.x()-=ell;
 minCorner.y()-=ell;
 minCorner.z()-=ell;
 delta.x()+=dd;
 delta.y()+=dd;
 delta.z()+=dd;

 int inside=0;
 for(int i = 0; i < nStat; ++i )
 {
   p = minCorner + Vector3D<Precision>( delta.x()*RNG::Instance(). uniform(),
           delta.y()*RNG::Instance(). uniform(),
           delta.z()*RNG::Instance(). uniform() );
   if( this->UnplacedContains(p) ) {
     if( this->SafetyToOut(p)<ell) { inside++; }
   }
   else{
     if( this->SafetyToIn(p)<ell) { inside++; }
   }
}
 // @@ The conformal correction can be upgraded
 return delta.x()*delta.y()*delta.z()*inside/dd/nStat;
}

} // End impl namespace

#ifdef VECGEOM_NVCC

namespace cxx {

template size_t DevicePtr<cuda::VPlacedVolume const*>::SizeOf();
template size_t DevicePtr<char>::SizeOf();
template size_t DevicePtr<Precision>::SizeOf();
// template void DevicePtr<cuda::PlacedBox>::Construct(
//    DevicePtr<cuda::LogicalVolume> const logical_volume,
//    DevicePtr<cuda::Transformation3D> const transform,
//    const int id) const;

} // End cxx namespace

#endif // VECGEOM_NVCC

} // End global namespace
