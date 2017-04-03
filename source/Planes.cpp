#include "volumes/Planes.h"
#include "backend/Backend.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

VECCORE_ATT_HOST_DEVICE
Planes::Planes(int size) : fNormals(size), fDistances(size)
{
}

#ifdef VECCORE_CUDA
__device__ Planes::Planes() : fNormals(), fDistances()
{
}
__device__ Planes::Planes(Precision *a, Precision *b, Precision *c, Precision *d, int size)
    : fNormals(a, b, c, size), fDistances(d, size)
{
}
#endif

VECCORE_ATT_HOST_DEVICE
Planes::~Planes()
{
}

VECCORE_ATT_HOST_DEVICE
Planes &Planes::operator=(Planes const &rhs)
{
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  fNormals   = rhs.fNormals;
  fDistances = rhs.fDistances;
#else
  fNormals = SOA3D<Precision>(const_cast<Precision *>(rhs.fNormals.x()), const_cast<Precision *>(rhs.fNormals.y()),
                              const_cast<Precision *>(rhs.fNormals.z()), rhs.fNormals.size());
  fDistances = Array<Precision>(const_cast<Precision *>(&rhs.fDistances[0]), rhs.fDistances.size());
#endif
  return *this;
}

VECCORE_ATT_HOST_DEVICE
void Planes::Set(int index, Vector3D<Precision> const &normal, Vector3D<Precision> const &x0)
{
  Vector3D<Precision> fixedNormal(normal);
  fixedNormal.FixZeroes();
  Precision inverseLength = 1. / fixedNormal.Mag();
  fNormals.set(index, inverseLength * fixedNormal);
  fDistances[index] = inverseLength * -fixedNormal.Dot(x0);
}

VECCORE_ATT_HOST_DEVICE
void Planes::Set(int index, Vector3D<Precision> const &normal, Precision distance)
{
  fNormals.set(index, normal);
  fDistances[index] = distance;
}

void Planes::FlipSign(int index)
{
  fNormals.set(index, -fNormals[index]);
  fDistances[index] = -fDistances[index];
}

std::ostream &operator<<(std::ostream &os, Planes const &planes)
{
  for (int i = 0, iMax = planes.size(); i < iMax; ++i) {
    os << "{" << planes.GetNormal(i) << ", " << planes.GetDistance(i) << "}\n";
  }
  return os;
}

} // End inline implementation namespace

} // End global namespace
