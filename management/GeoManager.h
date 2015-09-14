/// \file GeoManager.h
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#ifndef VECGEOM_MANAGEMENT_GEOMANAGER_H_
#define VECGEOM_MANAGEMENT_GEOMANAGER_H_

#ifdef OFFLOAD_MODE
#pragma offload_attribute(push, target(mic))
#endif

#include "base/Global.h"

#include "volumes/PlacedVolume.h"
#include "volumes/LogicalVolume.h"

#include <map>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

// probably don't need apply to be virtual
template<typename Container>
class GeoVisitor
{
protected:
   Container & c_;
public:
   GeoVisitor( Container & c ) : c_(c) {};

   virtual void apply( VPlacedVolume *, int level=0 ) = 0;
   virtual ~GeoVisitor(){}
};


class NavigationState;
template<typename Container>
class GeoVisitorWithAccessToPath
{
protected:
   Container & c_;
public:
   GeoVisitorWithAccessToPath( Container & c ) : c_(c) {};

   virtual void apply( NavigationState * state, int level=0 ) = 0;
   virtual ~GeoVisitorWithAccessToPath( ) { }
};


template<typename Container>
class SimpleLogicalVolumeVisitor : public GeoVisitor<Container>
{
public:
   SimpleLogicalVolumeVisitor( Container & c ) : GeoVisitor<Container>(c) {}
   virtual void apply( VPlacedVolume * vol, int /*level*/ ){
      LogicalVolume const *lvol = vol->GetLogicalVolume();
      if( std::find( this->c_.begin(), this->c_.end(), lvol ) == this->c_.end() )
      {
         this->c_.push_back( const_cast<LogicalVolume *>(lvol) );
      }
   }
   virtual ~SimpleLogicalVolumeVisitor(){}
};

template<typename Container>
class SimplePlacedVolumeVisitor : public GeoVisitor<Container>
{
public:
   SimplePlacedVolumeVisitor( Container & c) : GeoVisitor<Container>(c) {}
   virtual void apply( VPlacedVolume * vol, int /* level */ ){
      this->c_.push_back( vol );
   }
   virtual ~SimplePlacedVolumeVisitor(){}
};

class GetMaxDepthVisitor
{
private:
   int maxdepth_;
public:
   GetMaxDepthVisitor() : maxdepth_(0) {}
   void apply( VPlacedVolume * /* vol */, int level )
   {
      maxdepth_ = (level>maxdepth_) ? level : maxdepth_;
   }
   int getMaxDepth( ) const {return maxdepth_;}
};

class GetTotalNodeCountVisitor
{
private:
    int fTotalNodeCount;
public:
    GetTotalNodeCountVisitor() : fTotalNodeCount(0) {}
    void apply( VPlacedVolume *, int /* level */ )
    {
        fTotalNodeCount++;
    }
    int GetTotalNodeCount() const {return fTotalNodeCount;}
};




/**
 * @brief Knows about the current world volume.
 */
class GeoManager {

private:

  int fVolumeCount;
  int fTotalNodeCount; // total number of nodes in the geometry tree
  VPlacedVolume const *fWorld;

  // consider making these things rvalues
  std::map<unsigned int, VPlacedVolume *> fPlacedVolumesMap;
  std::map<unsigned int, LogicalVolume *> fLogicalVolumesMap;
  std::map<VPlacedVolume const *, unsigned int> fVolumeToIndexMap;
  int fMaxDepth;
  bool fIsClosed;

  // traverses the geometry tree of placed volumes and applies injected Visitor
  template<typename Visitor>
  void visitAllPlacedVolumes(VPlacedVolume const *, Visitor * visitor, int level=1 ) const;

  // traverse the geometry tree keeping track of the state context ( volume path or navigation state )
  // apply the injected Visitor
  template<typename Visitor>
  void visitAllPlacedVolumesWithContext(VPlacedVolume const *, Visitor * visitor, NavigationState * state, int level=1 ) const;


public:

  static VPlacedVolume *gCompactPlacedVolBuffer;

  static GeoManager& Instance() {
    static GeoManager instance;
    return instance;
  }

  /**
   * mark the current detector geometry as finished and initialize
   * important cached variables such as the maximum tree depth etc.
   */
  void CloseGeometry();

  // compactify memory space
  // an internal method which should be called by ClosedGeometry
  // it analyses the geometry and puts objects in contiguous buffers
  // it also fixes resulting pointer inconsistencies
  void CompactifyMemory();

  void SetWorld(VPlacedVolume const *const w) { fWorld = w; }

  VPlacedVolume const* GetWorld() const { return fWorld; }

  // initialize geometry from a precompiled shared library
  // ( such as obtained from the CppExporter )
  // this function sets the world and closes the geometry
  void LoadGeometryFromSharedLib( std::string );

  VPlacedVolume const *Convert( unsigned int index ) {
      return fPlacedVolumesMap[index];
  }
  unsigned int Convert( VPlacedVolume const * pvol ) {
      return fVolumeToIndexMap[pvol];
  }

  /**
   *  give back container containing all logical volumes in detector
   *  Container is supposed to be any Container that can store pointers to
   */
  template<typename Container>
  void getAllLogicalVolumes( Container & c ) const;

  /**
   *  give back container containing all logical volumes in detector
   */
  template<typename Container>
  void getAllPlacedVolumes( Container & c ) const;

  /**
   *  give back container containing all logical volumes in detector
   *  Container has to be an (stl::)container keeping pointer to NavigationStates
   */
  template<typename Container>
  void getAllPathForLogicalVolume( LogicalVolume const * lvol, Container & c ) const;

  /**
   *  return max depth of volume hierarchy
   */
  int getMaxDepth() const {
      Assert( fMaxDepth > 0, "geometry not closed" );
      return fMaxDepth;
  }

  void RegisterPlacedVolume(VPlacedVolume *const placed_volume);

  void RegisterLogicalVolume(LogicalVolume *const logical_volume);

  void DeregisterPlacedVolume(const int id);

  void DeregisterLogicalVolume(const int id);

  /**
   * \return Volume with passed id, or NULL is the id wasn't found.
   */
  VPlacedVolume* FindPlacedVolume(const int id);

  /**
   * \return First occurrence of volume with passed label. If multiple volumes
   *         are found, their id will be printed to standard output.
   */
  VPlacedVolume* FindPlacedVolume(char const *const label);


  /**
   * \return Volume with passed id, or NULL is the id wasn't found.
   */
  LogicalVolume* FindLogicalVolume(const int id);


  /**
   * \return First occurrence of volume with passed label. If multiple volumes
   *         are found, their id will be printed to standard output.
   */
  LogicalVolume* FindLogicalVolume(char const *const label);

  /**
   * Clear/reset the manager
   *
   */
   void Clear();


   int GetPlacedVolumesCount() const {return fPlacedVolumesMap.size();}
   int GetLogicalVolumesCount() const {return fLogicalVolumesMap.size();}
   int GetTotalNodeCount() const {return fTotalNodeCount;}

protected:

private:
 GeoManager() : fVolumeCount(0), fTotalNodeCount(0), fWorld(NULL), fPlacedVolumesMap(),
     fLogicalVolumesMap(), fVolumeToIndexMap(), fMaxDepth(-1), fIsClosed(false)
 {}

  GeoManager(GeoManager const&);
  GeoManager& operator=(GeoManager const&);
};

template<typename Visitor>
void
GeoManager::visitAllPlacedVolumes( VPlacedVolume const * currentvolume, Visitor * visitor, int level ) const
{
   if( currentvolume != NULL )
   {
      visitor->apply( const_cast<VPlacedVolume *>(currentvolume), level );
      int size = currentvolume->GetDaughters().size();
      for( int i=0; i<size; ++i )
      {
         visitAllPlacedVolumes( currentvolume->GetDaughters().operator[](i), visitor, level+1 );
      }
   }
}

template<typename Container>
void GeoManager::getAllLogicalVolumes( Container & c ) const
{
   c.clear();
   // walk all the volume hierarchy and insert
   // logical volume if not already in the container
   SimpleLogicalVolumeVisitor<Container> lv(c);
   visitAllPlacedVolumes( GetWorld(), &lv );
}


template<typename Container>
void GeoManager::getAllPlacedVolumes( Container & c ) const
{
   c.clear();
   // walk all the volume hierarchy and insert
   // placed volumes if not already in the container
   SimplePlacedVolumeVisitor<Container> pv(c);
   visitAllPlacedVolumes( GetWorld(), &pv );
}






} } // End global namespace

#ifdef OFFLOAD_MODE
#pragma offload_attribute(pop)
#endif

#endif // VECGEOM_MANAGEMENT_GEOMANAGER_H_
