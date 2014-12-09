/// \file NavigationState.h
/// \author Sandro Wenzel (sandro.wenzel@cern.ch)
/// \date 12.03.2014

#ifndef VECGEOM_NAVIGATION_NAVIGATIONSTATE_H_
#define VECGEOM_NAVIGATION_NAVIGATIONSTATE_H_

#include "backend/Backend.h"
#include "base/Transformation3D.h"
#include "volumes/PlacedVolume.h"

#ifdef VECGEOM_ROOT
#include "management/RootGeoManager.h"
#endif

#include <iostream>
#include <string>

class TGeoBranchArray;


namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

/**
 * a class describing a current geometry state
 * likely there will be such an object for each
 * particle/track currently treated
 */
class NavigationState
{
private:
   int fMaxlevel;
   int fCurrentLevel;
   // add other navigation state here, stuff like:
   bool fOnBoundary; // flag indicating whether track is on boundary of the "Top()" placed volume
   mutable Transformation3D global_matrix_;

   // pointer data follows; has to be last
   VPlacedVolume const * * fPath; //
   VPlacedVolume const * fBuffer[1]; // the real buffer to store the path

   // constructors and assignment operators are private
   // states have to be constructed using MakeInstance() function
   VECGEOM_INLINE
   VECGEOM_CUDA_HEADER_BOTH
   NavigationState( int );

   VECGEOM_INLINE
   VECGEOM_CUDA_HEADER_BOTH
   NavigationState( NavigationState const & rhs );

   // some private management methods
   VECGEOM_INLINE
   VECGEOM_CUDA_HEADER_BOTH
   void InitInternalStorage();

private:

   // The data start should point to the address of the first data member,
   // after the virtual table
   // the purpose is probably for the Copy function
   const void*  DataStart() const {return (const void*)&fMaxlevel;}
   const void*  ObjectStart() const {return (const void*)this;}
   void*  DataStart() {return (void*)&fMaxlevel;}
   void*  ObjectStart() {return (void*)this;}


     // The actual size of the data for an instance, excluding the virtual table
   size_t      DataSize() const {
      return SizeOf() + (size_t)ObjectStart() - (size_t)DataStart();}


public:

   VECGEOM_CUDA_HEADER_BOTH
   // produces a compact navigation state object of a certain depth
   // the caller can give a memory address where the object will
   // be placed
   // the caller has to make sure that the size of the external memory
   // is >= sizeof(NavigationState) + sizeof(VPlacedVolume*)*maxlevel
   //
   // TODO: It would probably be clearer to the user to have 2 functions
   // MakeInstance and MakeInstanceAt
   static NavigationState* MakeInstance(int maxlevel, void *addr=0);

   // The equivalent of the copy constructor
   // MakeCopy + MakeCopyAt?
   static NavigationState* MakeCopy(NavigationState const & other, void *addr=0);


   // static function to release an instance
   static void ReleaseInstance( NavigationState *instance ){}

   // returns the size in bytes of a NavigationState object with internal
   // path depth maxlevel
   VECGEOM_CUDA_HEADER_BOTH
   static int SizeNeeded( int maxlevel ){
     return sizeof(NavigationState) + sizeof(VPlacedVolume*)*(maxlevel+1);
   }

   // same as above ( to be equivalent with TGeoBranchArray in ROOT )
   VECGEOM_CUDA_HEADER_BOTH
   static int SizeOf( int maxlevel ){
     return sizeof(NavigationState) + sizeof(VPlacedVolume*)*(maxlevel+1);
   }

   //
   VECGEOM_CUDA_HEADER_BOTH
   int GetObjectSize() const {
        return NavigationState::SizeNeeded(fMaxlevel);
   }

   VECGEOM_CUDA_HEADER_BOTH
   int SizeOf() const {
       return NavigationState::SizeOf(fMaxlevel);
   }

   VECGEOM_INLINE
   VECGEOM_CUDA_HEADER_BOTH
   NavigationState & operator=( NavigationState const & rhs );

   // another copy method (inverse direction )
   VECGEOM_CUDA_HEADER_BOTH
   void CopyTo( NavigationState * other ) const {
       // Raw memcpy of the content to another existing state.
       //
       // in case NavigationState was a virtual class: change to
       // std::memcpy(other->DataStart(), DataStart(), DataSize());
       std::memcpy(other, this, this->SizeOf());
       other->fPath = &(other->fBuffer[0]);
   }


#ifdef VECGEOM_ROOT
   TGeoBranchArray * ToTGeoBranchArray() const;
   NavigationState & operator=( TGeoBranchArray const & rhs );
#endif

   VECGEOM_INLINE
   VECGEOM_CUDA_HEADER_BOTH
   ~NavigationState( );


   // what else: operator new etc...

   VECGEOM_INLINE
   VECGEOM_CUDA_HEADER_BOTH
   int GetMaxLevel() const {return fMaxlevel;}

   VECGEOM_INLINE
   VECGEOM_CUDA_HEADER_BOTH
   int GetCurrentLevel() const {return fCurrentLevel;}

   // better to use pop and push
   VECGEOM_INLINE
   VECGEOM_CUDA_HEADER_BOTH
   void
   Push(VPlacedVolume const *);

   VECGEOM_INLINE
   VECGEOM_CUDA_HEADER_BOTH
   VPlacedVolume const *
   Top() const;

   VECGEOM_INLINE
   VECGEOM_CUDA_HEADER_BOTH
   VPlacedVolume const *
   At(int level) const {return fPath[level];}

   VECGEOM_INLINE
   VECGEOM_CUDA_HEADER_BOTH
   Transformation3D const &
   TopMatrix() const;

   VECGEOM_INLINE
   VECGEOM_CUDA_HEADER_BOTH
   Vector3D<Precision>
   GlobalToLocal(Vector3D<Precision> const &);

   VECGEOM_INLINE
   VECGEOM_CUDA_HEADER_BOTH
   void Pop();

   VECGEOM_INLINE
   VECGEOM_CUDA_HEADER_BOTH
   int Distance( NavigationState const & ) const;
//   int Distance(NavigationState const &) const;

   // clear all information
   VECGEOM_INLINE
   VECGEOM_CUDA_HEADER_BOTH
   void Clear();

   VECGEOM_INLINE
   void Print() const;

   VECGEOM_INLINE
   VECGEOM_CUDA_HEADER_BOTH
   bool HasSamePathAsOther( NavigationState const & other ) const
   {
        if( other.fCurrentLevel != fCurrentLevel ) return false;
        for( int i= fCurrentLevel-1; i>=0; --i ){
            if( fPath[i] != other.fPath[i] ) return false;
        }
        return true;
   }

#ifdef VECGEOM_ROOT
   VECGEOM_INLINE
   void printVolumePath() const;

   /**
    * returns the number of FILLED LEVELS such that
    * state.GetNode( state.GetLevel() ) == state.Top()
    */
   VECGEOM_INLINE
   int GetLevel() const {return fCurrentLevel-1;}

   TGeoNode const * GetNode(int level) const {return
		   RootGeoManager::Instance().tgeonode( fPath[level] );}
#endif

   /**
     function returning whether the point (current navigation state) is outside the detector setup
   */
   VECGEOM_INLINE
   VECGEOM_CUDA_HEADER_BOTH
   bool IsOutside() const { return !(fCurrentLevel>0); }


   VECGEOM_INLINE
   VECGEOM_CUDA_HEADER_BOTH
   bool IsOnBoundary() const { return fOnBoundary; }

   VECGEOM_INLINE
   VECGEOM_CUDA_HEADER_BOTH
   void SetBoundaryState( bool b ) { fOnBoundary = b; }

#ifdef VECGEOM_ROOT
   /**
    * function return the ROOT TGeoNode object which is equivalent to calling Top()
    * function included for convenience; to make porting Geant-V easier; we should eventually get rid of this function
    */
   VECGEOM_INLINE
   TGeoNode const * GetCurrentNode() const
   {
      return RootGeoManager::Instance().tgeonode(this->Top());
   }
#endif

   //void GetGlobalMatrixFromPath( Transformation3D *const m ) const;
   //Transformation3D const * GetGlobalMatrixFromPath() const;
};

inline
NavigationState * NavigationState::MakeInstance(int maxlevel , void * addr)
{
  // Make an instance of the class which allocates the node array. To be
  // released using ReleaseInstance. If addr is non-zero, the user promised that
  // addr contains at least that many bytes:  size_t needed = SizeOf(maxlevel);
  NavigationState* ba = 0;
  if (!addr) {
     size_t needed = NavigationState::SizeNeeded(maxlevel);
     char *ptr = new char[ needed ];
     if (!ptr) return 0;
     new (ptr) NavigationState(maxlevel);
     ba = reinterpret_cast<NavigationState*>(ptr);
     // ba->SetBit(kBASelfAlloc, kTRUE);
  }
  else {
     new (addr) NavigationState(maxlevel);
     ba = reinterpret_cast<NavigationState*>(addr);
     // ba->SetBit(kBASelfAlloc, kFALSE);
  }
  return ba;
}

inline
NavigationState * NavigationState::MakeCopy( NavigationState const & rhs , void * addr)
{
  // this might be doing too much work
  NavigationState* copy = NavigationState::MakeInstance( rhs.fMaxlevel, addr );
  // copy info
  rhs.CopyTo(copy);
  return copy;
}


NavigationState & NavigationState::operator=( NavigationState const & rhs )
{
   fCurrentLevel=rhs.fCurrentLevel;
   fMaxlevel = rhs.fMaxlevel;
   fOnBoundary = rhs.fOnBoundary;
   // what about the matrix????
   std::memcpy(fPath, rhs.fPath, sizeof(*fPath)*fCurrentLevel);
   return *this;
}

/*
NavigationState::NavigationState( NavigationState const & rhs ) :
        fMaxlevel(rhs.fMaxlevel),
        fCurrentLevel(rhs.fCurrentLevel),
        fOnBoundary(rhs.fOnBoundary),
        global_matrix_() ,
        fPath(&fBuffer[0])
{
   InitInternalStorage();
   std::memcpy(fPath, rhs.fPath, sizeof(*fPath)*rhs.fCurrentLevel );
}
*/

// private implementation of standard constructor
 NavigationState::NavigationState( int maxlevel ) :
         fMaxlevel(maxlevel),
         fCurrentLevel(0),
         fOnBoundary(false),
         global_matrix_(),
         fPath(&fBuffer[0]) /* point to immediately following buffer */
{
   // clear the buffer
   std::memset(fBuffer, 0, fMaxlevel*sizeof(VPlacedVolume*));
}

  VECGEOM_CUDA_HEADER_BOTH
NavigationState::~NavigationState()
{
   delete[] fPath;
}


void
NavigationState::Pop()
{
   if(fCurrentLevel > 0){
       fPath[--fCurrentLevel]=0;
   }
}

void
NavigationState::Clear()
{
   fCurrentLevel=0;
   fOnBoundary=false;
}

void
NavigationState::Push( VPlacedVolume const * v )
{
#ifdef DEBUG
  assert( fCurrentLevel < fMaxlevel );
#endif
   fPath[fCurrentLevel++]=v;
}

VPlacedVolume const *
NavigationState::Top() const
{
   return (fCurrentLevel > 0 )? fPath[fCurrentLevel-1] : 0;
}

VECGEOM_INLINE
VECGEOM_CUDA_HEADER_BOTH
Transformation3D const &
NavigationState::TopMatrix() const
{
// this could be actually cached in case the path does not change ( particle stays inside a volume )
   global_matrix_.CopyFrom( *(fPath[0]->transformation()) );
   for(int i=1;i<fCurrentLevel;++i)
   {
      global_matrix_.MultiplyFromRight( *(fPath[i]->transformation()) );
   }
   return global_matrix_;
}

/**
 * function that transforms a global point to local point in reference frame of deepest volume in current navigation state
 * ( equivalent to using a global matrix )
 */
VECGEOM_INLINE
VECGEOM_CUDA_HEADER_BOTH
Vector3D<Precision>
NavigationState::GlobalToLocal(Vector3D<Precision> const & globalpoint)
{
   Vector3D<Precision> tmp=globalpoint;
   Vector3D<Precision> current;
   for(int level=0;level<fCurrentLevel;++level)
   {
      Transformation3D const *m = fPath[level]->transformation();
      current = m->Transform( tmp );
      tmp = current;
   }
   return tmp;
}

VECGEOM_INLINE
void NavigationState::Print() const
{
   std::cerr << "NavState: Level(cur/max)=" << fCurrentLevel <<'/'<< fMaxlevel
             <<" onBoundary="<< fOnBoundary
             <<" topVol="<< Top() <<" this="<< this
             << std::endl;
   // std::cerr << "maxlevel " << fMaxlevel << std::endl;
   // std::cerr << "currentlevel " << fCurrentLevel << std::endl;
   // std::cerr << "onboundary " << fOnBoundary << std::endl;
   // std::cerr << "deepest volume " << Top() << std::endl;
}


#ifdef VECGEOM_ROOT
VECGEOM_INLINE
/**
 * prints the path of the track as a verbose string ( like TGeoBranchArray in ROOT )
 * (uses internal root representation for the moment)
 */
void NavigationState::printVolumePath() const
{
   for(int i=0; i < fCurrentLevel; ++i)
   {
    std::cout << "/" << RootGeoManager::Instance().tgeonode( fPath[i] )->GetName();
   }
   std::cout << "\n";
}
#endif

/**
 * calculates if other navigation state takes a different branch in geometry path or is on same branch
 * ( two states are on same branch if one can connect the states just by going upwards or downwards ( or do nothing ))
 */
VECGEOM_INLINE
VECGEOM_CUDA_HEADER_BOTH
int NavigationState::Distance( NavigationState const & other ) const
{
   int lastcommonlevel=0;
   int maxlevel = Max( GetCurrentLevel() , other.GetCurrentLevel() );

   //  algorithm: start on top and go down until paths split
   for(int i=0; i < maxlevel; i++)
   {
      VPlacedVolume const *v1 = this->fPath[i];
      VPlacedVolume const *v2 = other.fPath[i];
      if( v1 == v2 )
      {
         lastcommonlevel = i;
      }
      else
      {
         break;
      }
   }
   return (GetCurrentLevel()-lastcommonlevel) + ( other.GetCurrentLevel() - lastcommonlevel ) - 2;
}


} } // End global namespace


#endif // VECGEOM_NAVIGATION_NAVIGATIONSTATE_H_
