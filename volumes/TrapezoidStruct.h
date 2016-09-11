/// @file: TrapezoidStruct.h
/// @author Guilherme Lima (lima@fnal.gov)
//
//  2016-07-22 Guilherme Lima  Created
//

#ifndef VECGEOM_VOLUMES_TRAPEZOIDSTRUCT_H_
#define VECGEOM_VOLUMES_TRAPEZOIDSTRUCT_H_
#include "base/Global.h"
#include "base/PlaneShell.h"
#include "VecCore/VecMath.h"

namespace vecgeom {

// using namespace vecCore::math;

inline namespace VECGEOM_IMPL_NAMESPACE {

/*
 * A Trapezoid struct to encapsulate the parameters and some other cached values
 * related to Trapezoid that are required in Implementation
 */
template <typename T = double>
struct TrapezoidStruct {

  struct TrapSidePlane {
    Precision fA, fB, fC, fD;
    // Plane equation: Ax+By+Cz+D=0, where
    // normal unit vector nvec=(A,B,C)  and offset=D is the distance from origin to plane

    VECGEOM_CUDA_HEADER_BOTH
    TrapSidePlane() : fA(0.0), fB(0.0), fC(0.0), fD(0.0) {}

    VECGEOM_CUDA_HEADER_BOTH
    TrapSidePlane(Precision a, Precision b, Precision c, Precision d) : fA(a), fB(b), fC(c), fD(d) {}
  };

#ifndef VECGEOM_PLANESHELL_DISABLE
  typedef PlaneShell<4, Precision> Planes;
#endif

  T fDz;
  T fTheta;
  T fPhi;
  T fDy1;
  T fDx1;
  T fDx2;
  T fTanAlpha1;
  T fDy2;
  T fDx3;
  T fDx4;
  T fTanAlpha2;

  // Values computed from parameters, to be cached
  T fTthetaCphi;
  T fTthetaSphi;

#ifndef VECGEOM_PLANESHELL_DISABLE
  Planes fPlanes;
#else
  TrapSidePlane fPlanes[4];
#endif

  T sideAreas[6]; // including z-planes

public:
  /// \brief Constructors
  /// @{
  // \brief General constructor.  All other constructors should delegate to it
  VECGEOM_CUDA_HEADER_BOTH
  TrapezoidStruct(const T pDz, const T pTheta, const T pPhi, const T pDy1, const T pDx1, const T pDx2,
                  const T pTanAlpha1, const T pDy2, const T pDx3, const T pDx4, const T pTanAlpha2)
      : fDz(pDz), fTheta(pTheta), fPhi(pPhi), fDy1(pDy1), fDx1(pDx1), fDx2(pDx2), fTanAlpha1(pTanAlpha1), fDy2(pDy2),
        fDx3(pDx3), fDx4(pDx4), fTanAlpha2(pTanAlpha2)
  {
    CalculateCached();
  }

  /// \brief Constructor for a "default" trapezoid whose parameters are to be set later
  VECGEOM_CUDA_HEADER_BOTH
  TrapezoidStruct() : TrapezoidStruct(0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.) {}

  /// \brief Copy constructor
  VECGEOM_CUDA_HEADER_BOTH
  TrapezoidStruct(TrapezoidStruct const &other)
      : TrapezoidStruct(other.fDz, other.fTheta, other.fPhi, other.fDy1, other.fDx1, other.fDx2, other.fTanAlpha1,
                        other.fDy2, other.fDx3, other.fDx4, other.fTanAlpha2)
  {
  }
  /// @}

  /// \brief Assignment operator
  VECGEOM_CUDA_HEADER_BOTH
  TrapezoidStruct &operator=(TrapezoidStruct const &other);

  /// \brief Destructor
  VECGEOM_CUDA_HEADER_BOTH
  virtual ~TrapezoidStruct() {}

  VECGEOM_CUDA_HEADER_BOTH
  void CalculateCached()
  {
    fTthetaCphi = vecCore::math::Tan(fTheta) * vecCore::math::Cos(fPhi);
    fTthetaSphi = vecCore::math::Tan(fTheta) * vecCore::math::Sin(fPhi);
    // ComputeBoundingBox();
  }

public:
#ifndef VECGEOM_PLANESHELL_DISABLE

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  Planes const *GetPlanes() const { return &fPlanes; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  TrapSidePlane GetPlane(unsigned int i) const
  {
    return TrapSidePlane(fPlanes.fA[i], fPlanes.fB[i], fPlanes.fC[i], fPlanes.fD[i]);
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  void SetPlane(unsigned int i, Precision a, Precision b, Precision c, Precision d)
  {
    fPlanes.fA[i] = a;
    fPlanes.fB[i] = b;
    fPlanes.fC[i] = c;
    fPlanes.fD[i] = d;
  }

#else

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  TrapSidePlane const *GetPlanes() const { return fPlanes; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  TrapSidePlane GetPlane(unsigned int i) const { return fPlanes[i]; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  void SetPlane(unsigned int i, Precision a, Precision b, Precision c, Precision d)
  {
    fPlanes[i].fA = a;
    fPlanes[i].fB = b;
    fPlanes[i].fC = c;
    fPlanes[i].fD = d;
  }
#endif
};

} // inline NS
} // vecgeom NS

#endif
