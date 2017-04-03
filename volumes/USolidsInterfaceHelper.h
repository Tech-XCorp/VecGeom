/// \file USolidsInterfaceHelper.h
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#ifndef VECGEOM_VOLUMES_USOLIDSINTERFACEHELPER_H_
#define VECGEOM_VOLUMES_USOLIDSINTERFACEHELPER_H_

#include "base/Global.h"

#ifndef VECGEOM_USOLIDS

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {
struct USolidsInterfaceHelper {
  VECCORE_ATT_HOST_DEVICE
  virtual ~USolidsInterfaceHelper() {}
};
}
}

#else // Compiling with USolids compatibility

#include "base/Vector3D.h"
#include "VUSolid.hh"
#include "UVector3.hh"
#include <string>

// class UVector3;

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

/// \brief USolids compatibility signatures.
///
/// These do not necessarily provide all the return values promised by the
/// interface, so use volumes in this way with caution.
class USolidsInterfaceHelper : public VUSolid {

public:
  //   VUSolid(const std::string &name);
  USolidsInterfaceHelper(const std::string &name) : VUSolid(name) {}
  USolidsInterfaceHelper() : VUSolid() {}

  VECCORE_ATT_HOST_DEVICE
  virtual Precision DistanceToOut(Vector3D<Precision> const &position, Vector3D<Precision> const &direction,
                                  Precision stepMax = kInfLength) const = 0;

  VECCORE_ATT_HOST_DEVICE
  virtual ~USolidsInterfaceHelper() {}

  VECCORE_ATT_HOST_DEVICE
  virtual Precision SafetyToOut(Vector3D<Precision> const &position) const = 0;

  VECCORE_ATT_HOST_DEVICE
  virtual Precision SafetyToIn(Vector3D<Precision> const &position) const = 0;

  // these function names are specific to USolids but can be reimplemented in terms of
  // other interfaces:

  virtual double SafetyFromOutside(Vector3D<double> const &point, bool accurate = false) const
  {
    (void)accurate; // fix 'unused variable' warning
    return SafetyToIn(point);
  }

  virtual double SafetyFromInside(Vector3D<double> const &point, bool accurate = false) const
  {
    (void)accurate; // fix 'unused variable' warning
    return SafetyToOut(point);
  }

  // the following functions are somewhat USolid specific
  // a default (asserting) implementation is implemented in USolidsInterfaceHelper.cpp
  virtual double DistanceToOut(Vector3D<double> const &point, Vector3D<double> const &direction,
                               Vector3D<double> &normal, bool &convex, double stepMax = kInfLength) const;

  virtual std::string GetEntityType() const;

  virtual void GetParametersList(int number, double *array) const;

  virtual VUSolid *Clone() const;

#if defined(VECGEOM_USOLIDS)
  virtual std::ostream &StreamInfo(std::ostream &os) const;
#endif

  virtual UVector3 GetPointOnSurface() const;

  virtual void ComputeBBox(UBBox *aBox, bool aStore = false);
};
}
} // End global namespace

#endif // USolids defined

#endif // VECGEOM_VOLUMES_USOLIDSINTERFACEHELPER_H_
