/// \file RootGeoManager.h
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#ifndef VECGEOM_MANAGEMENT_ROOTMANAGER_H_
#define VECGEOM_MANAGEMENT_ROOTMANAGER_H_

#include "base/Global.h"

#include "base/TypeMap.h"

class TGeoNode;
class TGeoVolume;
#include "TGeoNode.h"
class TGeoMatrix;

namespace VECGEOM_NAMESPACE {

class LogicalVolume;
class Transformation3D;
class UnplacedBox;
class VPlacedVolume;
class VUnplacedVolume;

/// \brief Manager to handle interaction with ROOT geometry.
/// \details Allows integration with ROOT geometries for compatability reasons.
///          Is not necessary for VecGeom library, and will only have source
///          compiled if the VECGEOM_ROOT flag is set by the compiler, activated
///          with -DROOT=ON in CMake.
class RootGeoManager {

private:

  /** Remember pointer to generated world from imported ROOT geometry. */
  VPlacedVolume const* fWorld;

  BidirectionalTypeMap<VPlacedVolume const*, TGeoNode const*> fPlacedVolumeMap;
  BidirectionalTypeMap<VUnplacedVolume const*, TGeoShape const*> fUnplacedVolumeMap;
  BidirectionalTypeMap<LogicalVolume const*, TGeoVolume const*> fLogicalVolumeMap;
  BidirectionalTypeMap<Transformation3D const*, TGeoMatrix const*> fTransformationMap;

  int fVerbose;

public:

  /// Access singleton instance.
  static RootGeoManager& Instance() {
    static RootGeoManager instance;
    return instance;
  }

  /**
   * @return Most recently generated world from ROOT geometry. Will return NULL
   *         if no ROOT geometry has been imported.
   * @sa LoadRootGeometry()
   */
  VPlacedVolume const* world() const { return fWorld; }

  int GetVerboseLevel() const { return fVerbose; }

  TGeoNode const * tgeonode( VPlacedVolume const * p ) const
  {
      if( p==NULL ) return NULL;
      return (fPlacedVolumeMap[const_cast<VPlacedVolume*>(p)]);
  }

  /// Get placed volume that corresponds to a TGeoNode
  VPlacedVolume const * GetPlacedVolume( TGeoNode const * n ) const {
      if( n==NULL ) return NULL;
      return (fPlacedVolumeMap[n]);
  }

  char const *  GetName( VPlacedVolume const * p  ) const { return tgeonode(p)->GetName(); }
  void PrintNodeTable() const;

  void set_verbose(const int verbose) { fVerbose = verbose; }

  /**
   * Queries the global ROOT GeoManager for the top volume and recursively
   * imports and converts to VecGeom geometry.
   * Will register the imported ROOT geometry as the new world of the VecGeom
   * GeoManager singleton.
   *
   *
   * requires already initialized gGeoManager object in ROOT
   */
  void LoadRootGeometry();

  /**
   * same but can take a root file as input
   * initializes gGeoManager object in ROOT at the same time
   */
  void LoadRootGeometry(std::string);

  /**
   * Exports a VecGeom geometry to a ROOT geometry
   */
  void ExportToROOTGeometry( VPlacedVolume const *, std::string );

  /**
   * @brief Deletes all VecGeom geometry generated by this class.
   */
  void Clear();

  /// Converts a TGeoNode to a VPlacedVolume, recursively converting daughters.
  /// Will take care not to convert anything twice by checking the birectional
  /// map between ROOT and VecGeom geometry.
  VPlacedVolume* Convert(TGeoNode const *const node);

  VUnplacedVolume* Convert(TGeoShape const *const shape);

  LogicalVolume* Convert(TGeoVolume const *const volume);

  Transformation3D* Convert(TGeoMatrix const *const trans);

  // inverse process
  TGeoNode* Convert(VPlacedVolume const *const node);

  TGeoVolume* Convert(VPlacedVolume const *const, LogicalVolume const *const);

  TGeoMatrix* Convert(Transformation3D const *const trans);



private:

 RootGeoManager() : fWorld(NULL), fPlacedVolumeMap(), fUnplacedVolumeMap(),
    fLogicalVolumeMap(), fTransformationMap(), fVerbose(0) {}
  RootGeoManager(RootGeoManager const&);
  RootGeoManager& operator=(RootGeoManager const&);

};

} // End global namespace

#endif // VECGEOM_MANAGEMENT_ROOTMANAGER_H_
