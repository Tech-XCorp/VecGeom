#ifndef VECGEOM_SURFACE_MODEL_H_
#define VECGEOM_SURFACE_MODEL_H_

#include <VecGeom/surfaces/Equations.h>
#include <VecGeom/navigation/NavStateIndex.h>

namespace vgbrep {
///< Supported surface types
enum SurfaceType {kPlanar, kCylindrical, kConical, kSpherical, kTorus, kGenSecondOrder};

///< Supported frame types
enum FrameType {kRangeZ, kRangeCyl, kRangeSph, kWindow, kTriangle};

template <typename Real_t>
using Vector3D = vecgeom::Vector3D<Real_t>;

template <typename Real_t>
struct SurfData;

///< 
/* 
   Unplaced half-space surface type. The actual surface data pointed by the surface id
   is stored in a separate SurfData structure
   - All unplaced planes are (xOy), having the oriented normal on positive z
   - Unplaced cylinders, cones and tori have the z axis as axis of symmetry. Normal pointing outwards
   - Unplaced spheres have the origin as center

   To get the distance from a global point/direction to the surface one needs to first convert
   those to the local frame
 */
struct UnplacedSurface {
  SurfaceType        type {kPlanar}; ///< surface type
  int                id {-1};        ///< surface id

  UnplacedSurface() = default;
  UnplacedSurface(SurfaceType stype, int sid) { type = stype; id = sid; }
  
  /// Find positive distance to next intersection from local point
  template <typename Real_t>
  Real_t Intersect(Vector3D<Real_t> const &point, Vector3D<Real_t> const &dir,
                   SurfData<Real_t> const &storage)
  {
    QuadraticCoef<Real_t> coef;
    QuadraticRoots<Real_t> roots;
    Real_t dist = -1;

    switch (type) {
      case kPlanar:
        // Just need to propagate to (xOy) plane
        dist = -point[2] / dir[2];
        return dist;
      case kCylindrical:
        // Intersect with the cylindrical surface having Z as axis of symmetry
        CylinderEq<Real_t>(point, dir, storage.GetCylData(id), coef);
        break;
      case kConical:
        // Intersect with the conical surface having Z as axis of symmetry
        ConeEq<Real_t>(point, dir, storage.GetConeData(id), coef);
        break;
      case kSpherical:
        // Intersect with the sphere having the center in the origin
        SphereEq<Real_t>(point, dir, storage.GetSphereData(id), coef);
        break;
      case kTorus:
      case kGenSecondOrder:
        std::cout << "kTorus, kGenSecondOrder unhandled\n";
        return -1;
    };

    QuadraticSolver(coef, roots);
    if (roots.numroots > 0)
      return roots.xmin;
    return -1;
  }
};

/* An frame delimiting the real surface on an infinite half-space */
struct Frame {
  FrameType type {kWindow};  ///< frame type
  int       id {-1};         ///< frame data id

  Frame() = default;
  Frame(FrameType mtype, int mid) : type(mtype), id(mid) {}
    
  template <typename Real_t>
  bool Inside(Vector3D<Real_t> const &local, SurfData<Real_t> const &storage)
  {
    auto data = storage.GetFrameData(*this);
    Real_t rsq {0};
    switch (type) {
      case kRangeZ:
        return (local[2] > vecgeom::MakeMinusTolerant<Real_t>(data[0]) 
             && local[2] < vecgeom::MakePlusTolerant<Real_t>(data[1]));
      case kRangeCyl:
        rsq = local[0] * local[0] + local[1] * local[1];
        return (rsq > vecgeom::MakeMinusTolerantSquare<Real_t>(data[0])
             && rsq < vecgeom::MakePlusTolerantSquare<Real_t>(data[1]));
      case kRangeSph:
        rsq = local.Mag2();
        return (rsq > vecgeom::MakeMinusTolerantSquare<Real_t>(data[0])
             && rsq < vecgeom::MakePlusTolerantSquare<Real_t>(data[1]));
      case kWindow:
        return (Abs(local[0]) < vecgeom::MakePlusTolerant<Real_t>(data[0])
             && Abs(local[1]) < vecgeom::MakePlusTolerant<Real_t>(data[1]));
      case kTriangle:
        std::cout << "Triangle masks not handled yet\n";
    };
    return false;
  }
};

// This holds the transformation of the surface
// with respect to the frame of the ancestor volume onto which this surface is flattened.
//
// Example of a 4-level hierarchy flattened on 2 levels:
// A, B, C, D, E, F = logical volumes
// t1, t2, t3, t4, t5 = local volume transformations
//
//            A                     A + C(t2) + E(t2*t4)   (one scene given by volume A)
//      (t1) / \ (t2)                      |
//          B   C                     (t1) |
//    (t3) /     \ (t4)   -->              |
//        D       E                 B + D(t3) + F(t3*t5)   (another scene given by volume B)
//   (t5) |
//        F
//
// In the above example, say F has a surface S positioned with a local ransformation (ts).
// The global corresponding to S will have a total transformation (t1 * t3 * t5 * ts).
// However, it the tree is flattened on two scenes as above (one level is A and the other level is B)
// then the local transformation of S will be just (t3 * t5 * ts). Its global transformation 
// will be just (t1). In this approach, the transformation for a global surface is always obtained
// by multiplying the full scene transformation with the local surface transformation.
//
// The advantage of this approach is that it gives full flexibility for chosing the flattened
// volumes, and a given local surface can be referenced by multiple portals (less memory)

/* A placed surface on a scene having a frame and a navigation state associated to a touchable */
struct FramedSurface {
  UnplacedSurface  fSurface;     ///< Surface identifier
  Frame            fFrame;       ///< Frame
  int              fTrans {-1};  ///< Transformation of the surface in the compacted sub-hierarchy top volume frame
  NavIndex_t       fState {0};   ///< sub-path navigation state id in the parent scene

  FramedSurface() = default;
  FramedSurface(UnplacedSurface const &unplaced, Frame const &frame, int trans, NavIndex_t index = 0)
    : fSurface(unplaced), fFrame(frame), fTrans(trans), fState(index) {}

  /// Sorting by decreasing state depth and increasing state index 
  bool operator<(FramedSurface const &other)
  {
    using vecgeom::NavStateIndex;
    auto level1 = NavStateIndex::GetLevelImpl(fState);
    auto level2 = NavStateIndex::GetLevelImpl(other.fState);
    if (level1 > level2)
      return true;
    else if (level1 < level2)
      return false;
    if (fState < other.fState)
      return true;
    return false;
  }

  /// Transform point and direction to the local frame
  template <typename Real_t>
  void Transform(Vector3D<Real_t> const &point, Vector3D<Real_t> const &dir,
                 Vector3D<Real_t> &localpoint, Vector3D<Real_t> &localdir,
                 SurfData<Real_t>const &storage)
  {
    auto &localRef = storage.LocalT(fTrans);
    localpoint = localRef.Transform(point);
    localdir   = localRef.TransformDirection(dir);
  }

  ///< This finds the distance to intersecting the half-space, without checking the mask
  // Tte point and direction are in the reference frame of the scene
  template <typename Real_t>
  Real_t Intersect(Vector3D<Real_t> const &point, Vector3D<Real_t> const &dir,
                   SurfData<Real_t>const &storage)
  {
    Vector3D<Real_t> localpoint, localdir;
    Transform(point, dir, localpoint, localdir);
    return fSurface.Intersect<Real_t>(localpoint, localdir, storage);
  }
 
  ///< Check if the propagated point on surface is within the frame
  template <typename Real_t>
  bool InsideFrame(Vector3D<Real_t> const &point, Vector3D<Real_t> const &dir,
                     SurfData<Real_t> const &storage)
  {
    Vector3D<Real_t> localpoint, localdir;
    Transform(point, dir, localpoint, localdir);
    return fFrame.Inside(point, storage);
  }
  
};

///< A list of candidate surfaces
struct Candidates {
  int  fNcand {0};                ///< Number of candidate surfaces
  int *fCandidates {nullptr};     ///< [fNcand] Array of candidates

  Candidates() = default;
};

///< A side represents all common placed surfaces
struct Side {
  int  fNsurf {0};                ///< Number of placed surfaces on this side
  int *fSurfaces {nullptr};       ///< [fNsurf] Array of placed surfaces on this side

  Side() = default;

  // Add existing placed surface to this side
  int AddSurface(int isurf)
  {
    // Re-allocate policy to keep memory footprint low
    // Sides have to be relocated for GPU in contiguous memory
    int *surfaces = new int[fNsurf + 1];
    for (auto i = 0; i < fNsurf; ++i)
      surfaces[i] = fSurfaces[i];
    surfaces[fNsurf++] = isurf;
    delete [] fSurfaces;
    fSurfaces = surfaces;
    return fNsurf - 1;
  }

  size_t size() const { return sizeof(Side) + fNsurf * sizeof(int); }

  void CopyTo(char *buffer)
  {
    // to be implemented
  }

};

struct CommonSurface {
  SurfaceType fType {kPlanar}; ///< Type of surface
  int              fTrans {-1};  ///< Transformation of the first left frame

  NavIndex_t fDefaultState {0};   ///< The default state for this surface (deepest mother)
  // The portal property can be embedded in the general surface type
  
  Side             fLeftSide;  ///< Left-side portal side id (behind normal)
  Side             fRightSide; ///< Right-side portal side id (alongside normal)

  CommonSurface() = default;
  
  CommonSurface(SurfaceType type, int global_surf) : fType(type)
  {
    // Add by default the first surface to the left side
    fLeftSide.AddSurface(global_surf);
  };  
};

class BVH;

// A level of the geometry setup with a coordinate system and multiple volumes
// Currently called 'Universe' in Orange (name taken from MCNP, other codes)
// A detector or setup will be composed of one or two levels of Scene
struct Scene 
{
  CommonSurface  *fSurfaces {nullptr};    // Doors to other scenes
  BVH             *fNavigator {nullptr};
};
// How we decompose scene in hierarchical Scenes
// 

// If we introduce masks that are generated by 'boolean solids', we will likely need:
struct CompositeMask { /// ???   They must share the Surface Id !! (JA 2022/02/02 )
  int              fNumMasks;
  int             *fSimpleMasks;
  char            *fBooleanMaskString;  // "18 17 + 32 - 99 97 U *" // "Tube17mask Box333mask -" 
  // 
};

template <typename Real_t>
struct RangeMask {
  Real_t range[2];
};

template <typename Real_t, typename Real_s = Real_t>
struct CylData {
  Real_t           radius{0};     ///< Cylinder radius

  CylData() = default;
  CylData(Real_s rad) : radius(rad) {}
};


template <typename Real_t, typename Real_s = Real_t>
using SphData = CylData<Real_t, Real_s>;

template <typename Real_t, typename Real_s = Real_t>
struct ConeData {
  Real_t           radius{0};       ///< Cone radus at Z = 0
  Real_t           slope{0};     ///< Cone slope  --> for cyl extension this would be 0

  ConeData() = default;
  ConeData(Real_s rad, Real_s slope) : radius(rad), slope(slope) {}
};

// We need replicable surfaces, for memory reasons. All surfaces generated by a solid are positioned
// using identical local transformations, while the solid itself will use a different global transformation
// depending on the placement. So it makes sense to use tuples of global + local transformation. LocalSurface
// will use the local component, while portals will use the full tuple.
struct CombiTrans {
  int fLocalT{-1};      ///< index of the local surface transformation
  int fTouchableT{-1};  ///< index of touchable (containing the surface) transformation
};

/*
 * The main surface storage utility, providing access by index to:
 *    * global and local transformations applied to surfaces
 *    * surface data per surface type
 *    * mask data per mask type
 *    * imprint data (list of masks)
 */
template <typename Real_t>
struct SurfData {

  using Transformation = vecgeom::Transformation3D;
  using CylData_t = CylData<Real_t>;
  using ConeData_t = ConeData<Real_t>;
  using SphData_t = SphData<Real_t>;
  using RangeMask_t = RangeMask<Real_t>;

  int fNglobalTrans {0};
  int fNglobalSurf {0};
  int fNcommonSurf {0};
  int fNcylsph {0};
  int fNcone {0};
  int fNrange {0};

  /// Transformations. A portal transformation is a tuple global + local
  Transformation *fGlobalTrans {nullptr};   ///< Touchable global transformations

  /// Cylindrical surface data (radius)
  CylData_t      *fCylSphData {nullptr};       ///< Cyl and sphere data
  ConeData_t     *fConeData {nullptr};         ///< Cone data

  /// Mask data
  RangeMask_t    *fRangeData {nullptr};

  FramedSurface  *fGlobalSurfaces {nullptr}; ///< global surfaces
  CommonSurface  *fCommonSurfaces {nullptr}; ///< common surfaces
  Candidates     *fCandidates;               ///< candidate surfaces per navigation state
  int            *fSides {nullptr};          ///< side surface indices
  int            *fCandList {nullptr};       ///< global list of candidate indices

  
  SurfData() = default;  
  /// Transformation getters
  //Transformation const &LocalT(int id) const { return fLocalTrans[id]; }
  //Transformation const &TouchableT(int id) const { return fTouchableTrans[id]; }
  //void                  GlobalT(CombiTrans const &ct, Transformation &global)
  //{
  //  global = TouchableT(ct.fTouchableT);
  //  global.MultiplyFromRight(LocalT(ct.fLocalT));
  //}

  /// Surface data accessors
  CylData_t const        &GetCylData(int id) const { return fCylSphData[id]; }
  SphData_t const        &GetSphData(int id) const { return fCylSphData[id]; }
  ConeData_t const       &GetConeData(int id) const { return fConeData[id]; }
  RangeMask_t const      &GetRangeMask(int id) const { return fRangeData[id]; }

};

} // namespace vgbrep

#endif

