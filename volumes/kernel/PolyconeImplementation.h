/*
 * PolyconeImplementation.h
 *
 *  Created on: Dec 8, 2014
 *      Author: swenzel
 */


#ifndef VECGEOM_VOLUMES_KERNEL_POLYCONEIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_POLYCONEIMPLEMENTATION_H_

#include "base/Global.h"
#include "base/Transformation3D.h"
#include "volumes/UnplacedPolycone.h"
#include "volumes/kernel/ConeImplementation.h"

#include <cassert>
#include <cstdio>

namespace vecgeom {

  //VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_2v(PolyconeImplementation,
  //     TranslationCode,transCodeT, RotationCode,rotCodeT)

  VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_2v(PolyconeImplementation, TranslationCode, translation::kGeneric, RotationCode, rotation::kGeneric)

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedPolycone;

template <TranslationCode transCodeT, RotationCode rotCodeT>
struct PolyconeImplementation {

  static const int transC = transCodeT;
  static const int rotC   = rotCodeT;

  using PlacedShape_t = PlacedPolycone;
  using UnplacedShape_t = UnplacedPolycone;

  // here put all the implementations
  VECGEOM_CUDA_HEADER_BOTH
  static void PrintType() {
      printf("SpecializedPolycone<%i, %i>", transCodeT, rotCodeT);
  }

  /////GenericKernel Contains/Inside implementation for a section of the polycone
  template <typename Backend, bool ForInside>
  VECGEOM_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  static void GenericKernelForASection(UnplacedPolycone const &unplaced, int isect,
                                       Vector3D<typename Backend::precision_v> const &polyconePoint,
                                       typename Backend::bool_v &secFullyInside,
                                       typename Backend::bool_v &secFullyOutside) {

    if(isect<0) {
      secFullyInside = false;
      secFullyOutside = true;
      return;
    }

    PolyconeSection const & sec = unplaced.GetSection(isect);
    Vector3D<Precision> secLocalp = polyconePoint - Vector3D<Precision>(0,0,sec.fShift);
#ifdef POLYCONEDEBUG
    std::cout<<" isect="<< isect <<"/"<< unplaced.GetNSections()
             <<" secLocalP="<< secLocalp <<", secShift="<< sec.fShift
             <<" sec.fSolid="<< sec.fSolid
             << std::endl;
    if( sec.fSolid ) sec.fSolid->StreamInfo(std::cout);
#endif
    ConeImplementation< translation::kIdentity, rotation::kIdentity,
                        ConeTypes::UniversalCone>::GenericKernelForContainsAndInside<Backend,ForInside>(
                          *sec.fSolid, secLocalp, secFullyInside, secFullyOutside );
  }

  /////GenericKernel Contains/Inside implementation
  template <typename Backend, bool ForInside>
  VECGEOM_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  static void GenericKernelForContainsAndInside(UnplacedPolycone const &unplaced,
                                                Vector3D<typename Backend::precision_v> const &localPoint,
                                                typename Backend::bool_v &completelyInside,
                                                typename Backend::bool_v &completelyOutside)
    {
      // ??? Can I assume that this code will never be used in vector mode?  If yes, some of this code is not needed (done, IsFull(), etc)
      // static_assert(!std::is_same<Backend::bool_v,bool>::value, "Must use scalar interface!\n");

#ifdef POLYCONEDEBUG
      std::cout<<"=== Polycone::GenKern4ContainsAndInside(): localPoint = "<< localPoint
               <<" and kTolerance = "<< kTolerance <<"\n";
#endif
      typedef typename Backend::precision_v Float_t;
      typedef typename Backend::bool_v Bool_t;
      Bool_t done = Bool_t(false);

      // z-region
      int Nz = unplaced.GetNz();
      completelyInside = Bool_t(false);
      completelyOutside = localPoint[2] < MakeMinusTolerant<ForInside>( unplaced.GetZAtPlane(0) );
      completelyOutside |= localPoint[2] > MakePlusTolerant<ForInside>( unplaced.GetZAtPlane(Nz-1) );
      done |= completelyOutside;

#ifdef POLYCONEDEBUG
      std::cout<<" spot 2: checked z-region - compOut="<< completelyOutside
               <<", zAtPlane[0]="<< unplaced.GetZAtPlane(0)
               <<", zAtPlane[Nz-1]="<< unplaced.GetZAtPlane(Nz-1)
               <<", done="<< done <<"\n";
#endif

      // ??? Can I assume that this code will never be used in vector mode?  If yes, some of this code is not needed (done, IsFull(), etc)
      if ( Backend::early_returns && IsFull(completelyOutside) ) return;
#ifdef POLYCONEDEBUG
      if ( Backend::early_returns && IsFull(completelyOutside) ) {
        std::cout<<"Polycone::GenKern4ConAndInside() (spot 2a): returning compIn="<< completelyInside <<", compOut="<< completelyOutside <<"\n";
      }
#endif

      if(ForInside) {
        // z-check only so far
        completelyInside = ( localPoint[2] > MakePlusTolerant<ForInside>( unplaced.GetZAtPlane(0) ) &&
                             localPoint[2] < MakeMinusTolerant<ForInside>( unplaced.GetZAtPlane(Nz-1) ) );
#ifdef POLYCONEDEBUG
        std::cout<<" spot 2b: ForInside only: compIn="<< completelyInside <<", compOut="<< completelyOutside <<"\n";
#endif
      }

      // find section
      int isec = unplaced.GetSectionIndex( localPoint.z() );
      Bool_t secIn=false, secOut=false;
      // test if point is inside this section (note that surface may be very tricky!)
#ifdef POLYCONEDEBUG
      std::cout<<"Polycone::GenKern4ContAndInside(): spot 3: ";
#endif

      // is it inside of isec?
      GenericKernelForASection<Backend,ForInside>(unplaced, isec, localPoint, secIn, secOut );
      // if fully inside that section, we may be done
      Bool_t update = (!done && secIn);
      MaskedAssign(update, secIn, &completelyInside);
      done |= update;

#ifdef POLYCONEDEBUG
      if ( Backend::early_returns && IsFull(done) ) {
        std::cout<<"Polycone::GenKern4ConAndInside() (spot 2a): returning compIn="<< completelyInside <<", compOut="<< completelyOutside <<"\n";
      }
#endif
      if( Backend::early_returns && IsFull(done) ) return;

      // once here, need to check if near z-plane
      //========================================
      //-- Check a possible scenario that can produce incorrect results:
      //  I have seen cases where isec returns outside, but point is right at the Z-plane
      //  between isec & isec+1 -- and isec+1 returns Surface.  We must check it here first!
      Float_t zplus = unplaced.GetZAtPlane(isec);
      Bool_t nearPlusZ = localPoint.z() > zplus-kTolerance;
      Bool_t nextSectionPossible = nearPlusZ && isec+1<unplaced.GetNSections();

      // away from z-plane and outside isec
      update = (!done && !nextSectionPossible && secOut);
      completelyOutside |= update;
      done |= update;

      // away from z-plane and on isec surface
      update = (!done && !nextSectionPossible && !secOut && !secIn);
      MaskedAssign( update, false, &completelyInside );
      done |= update;

#ifdef POLYCONEDEBUG
      if ( Backend::early_returns && IsFull(done) ) {
        std::cout<<"Polycone::GenKern4ConAndInside() (spot 2b): returning compIn="<< completelyInside <<", compOut="<< completelyOutside <<"\n";
      }
#endif
      if( Backend::early_returns && IsFull(done) ) return;

      // ok, we are indeed near z-plane - need to check next surface
      assert( IsFull(nextSectionPossible) );

      Bool_t secIn2=false, secOut2=false;
      // very uncommon case -- OK to pay the penalty of a branch.  Any other options?!
#ifdef POLYCONEDEBUG
      std::cout<<"Polycone::GenKern4ContAndInside(): spot pre-3 (RARE): ";
#endif
      GenericKernelForASection<Backend,ForInside>(unplaced, isec+1, localPoint, secIn2, secOut2 );

      // near z-plane: outside in both
      update = (!done && nextSectionPossible && secOut && secOut2);
      completelyOutside |= update;
      done |= update;

      if(ForInside) {
        // near z-plane: check if inside radii on *both* sections
        Float_t localx = localPoint.x();
        Float_t localy = localPoint.y();
        Float_t rho2 = localx*localx + localy*localy;

        // for isec...
        Float_t rmaxtol = unplaced.GetRmaxAtPlane(isec) - kTolerance;
        Bool_t insideRmaxA = (rho2 < rmaxtol*rmaxtol);
        Float_t rmin = unplaced.GetRminAtPlane(isec);
        Bool_t insideRminA = (rmin == 0);
        // take tolerance into account
        rmin += kTolerance;
        insideRminA |= (rho2 > rmin*rmin);

        // for isec+1...
        rmaxtol = unplaced.GetRmaxAtPlane(isec+1) - kTolerance;
        Bool_t insideRmaxB = (rho2 < rmaxtol*rmaxtol);
        rmin = unplaced.GetRminAtPlane(isec+1);
        Bool_t insideRminB = (rmin == 0);
        // take tolerance into account
        rmin += kTolerance;
        insideRminB |= (rho2 > rmin*rmin);

        // near z-plane AND inside radii on *both* sections
        update = (!done && insideRminA && insideRmaxA && insideRminB && insideRmaxB );
        completelyInside |= update;
        done |= update;

        // .... hopefully last case left...  surface -- 
        MaskedAssign(!done, false, &completelyInside);
        MaskedAssign(!done, false, &completelyOutside);
      } // end if(ForInside)

#ifdef POLYCONEDEBUG
      if ( Backend::early_returns && IsFull(done) )  {
          std::cout<<"Section::GenKern4ConAndInside() spot 4: returning from section - secIn="<< secIn
                   <<", secOut="<< secOut
                   <<", done="<< done
                   <<", compIn="<< completelyInside
                   <<", compOut="<< completelyOutside
                   <<"\n";
      }
#endif

#ifdef POLYCONEDEBUG
      std::cout<<" Polycone::GenKern4ConAndInside(): spot FINAL - near surf: compIn="<< completelyInside <<", compOut="<< completelyOutside <<"\n";
#endif
    } // end of GenericKernelForContainsAndInside()

    template <class Backend>
    VECGEOM_CUDA_HEADER_BOTH
    VECGEOM_INLINE
    static void ContainsKernel(
        UnplacedPolycone const &polycone,
        Vector3D<typename Backend::precision_v> const &point,
        typename Backend::bool_v &contains)
    {
        // add z check
        if( point.z() < polycone.fZs[0] || point.z() > polycone.fZs[polycone.GetNSections()] )
        {
            contains = Backend::kFalse;
            return;
        }

        // now we have to find a section
        PolyconeSection const & sec = polycone.GetSection(point.z());

        Vector3D<Precision> localp;
        ConeImplementation< translation::kIdentity, rotation::kIdentity, ConeTypes::UniversalCone>::Contains<Backend>(
                *sec.fSolid,
                Transformation3D(),
                point - Vector3D<Precision>(0,0,sec.fShift),
                localp,
                contains
        );
        return;
    }

    template <class Backend>
    VECGEOM_CUDA_HEADER_BOTH
    VECGEOM_INLINE
    static void UnplacedContains( UnplacedPolycone const &polycone,
        Vector3D<typename Backend::precision_v> const &point,
        typename Backend::bool_v &contains) {
      // TODO: do this generically WITH a generic contains/inside kernel
      // forget about sector for the moment
      ContainsKernel<Backend>(polycone, point, contains);

      typedef typename Backend::bool_v Bool_t;
      Bool_t unused, outside;
      GenericKernelForContainsAndInside<Backend, false>(polycone, point, unused, outside);
      Bool_t contains2 = !outside;
#if POLYCONEDEBUG
      if(contains != contains2) {
        polycone.StreamInfo(std::cout);
        std::cout<<" Comparing Contains() fails!!!  old="<< contains <<", new="<< contains2 << std::endl;
      }
#endif
      assert( contains == contains2 );
    }

    template <typename Backend>
    VECGEOM_INLINE
    VECGEOM_CUDA_HEADER_BOTH
    static void Contains(
        UnplacedPolycone const &unplaced,
        Transformation3D const &transformation,
        Vector3D<typename Backend::precision_v> const &point,
        Vector3D<typename Backend::precision_v> &localPoint,
        typename Backend::bool_v &contains){

        localPoint = transformation.Transform<transCodeT, rotCodeT>(point);
        UnplacedContains<Backend>(unplaced, localPoint, contains);
    }

    template <class Backend>
    VECGEOM_CUDA_HEADER_BOTH
    VECGEOM_INLINE
    static void Inside(UnplacedPolycone const &unplaced,
                       Transformation3D const &transformation,
                       Vector3D<typename Backend::precision_v> const &masterPoint,
                       typename Backend::inside_v &inside) {

      typedef typename Backend::bool_v Bool_t;

      // convert from master to local coordinates
      Vector3D<typename Backend::precision_v> localPoint =
        transformation.Transform<transCodeT, rotCodeT>(masterPoint);

      std::cout<<" Polycone::Inside(): localz="<< localPoint.z() <<"\n";

      Bool_t fullyInside, fullyOutside;
      GenericKernelForContainsAndInside<Backend,true>( unplaced, localPoint, fullyInside, fullyOutside);

      inside = EInside::kSurface;
      MaskedAssign(fullyInside,  EInside::kInside,  &inside);
      MaskedAssign(fullyOutside, EInside::kOutside, &inside);
  }


      template <class Backend>
      VECGEOM_CUDA_HEADER_BOTH
      VECGEOM_INLINE
      static void DistanceToIn(
          UnplacedPolycone const &polycone,
          Transformation3D const &transformation,
          Vector3D<typename Backend::precision_v> const &point,
          Vector3D<typename Backend::precision_v> const &direction,
          typename Backend::precision_v const &stepMax,
          typename Backend::precision_v &distance) {

         Vector3D<typename Backend::precision_v> p = transformation.Transform<transCodeT,rotCodeT>(point);
         Vector3D<typename Backend::precision_v> v = transformation.TransformDirection<rotCodeT>(direction);
        // TODO
        // TODO: add bounding box check maybe??

         distance=kInfinity;
        int increment = (v.z() > 0) ? 1 : -1;
        if (std::fabs(v.z()) < kTolerance) increment = 0;
        int index = polycone.GetSectionIndex(p.z());
        if(index == -1) index = 0;
        if(index == -2) index = polycone.GetNSections()-1;

        do{
        // now we have to find a section

         PolyconeSection const & sec = polycone.GetSection(index);
        
         ConeImplementation< translation::kIdentity, rotation::kIdentity, ConeTypes::UniversalCone>::DistanceToIn<Backend>(
                *sec.fSolid,
                Transformation3D(),
                p - Vector3D<Precision>(0,0,sec.fShift),
                v,
                stepMax,
                distance);
     // std::cout<<"dist="<<distance<<" index="<<index<<" p.z()="<<p.z()-sec.fShift<<std::endl;

        if (distance < kInfinity || !increment)
        break;
        index += increment;
      
       }
       while (index >= 0 && index < polycone.GetNSections());
        
       return;
    }

    template <class Backend>
    VECGEOM_CUDA_HEADER_BOTH
    VECGEOM_INLINE
    static void DistanceToOut(
        UnplacedPolycone const &polycone,
        Vector3D<typename Backend::precision_v> const &point,
        Vector3D<typename Backend::precision_v> const &dir,
        typename Backend::precision_v const &stepMax,
        typename Backend::precision_v &distance) {

    Vector3D<typename Backend::precision_v>  pn(point);

    // specialization for N==1 ??? It should be a cone in the first place
    if (polycone.GetNSections()==1) {
     const PolyconeSection& section = polycone.GetSection(0);
       
     ConeImplementation< translation::kIdentity, rotation::kIdentity,
        ConeTypes::UniversalCone>::DistanceToOut<Backend>(
                *section.fSolid,
                point - Vector3D<Precision>(0,0,section.fShift),dir,stepMax,distance);
     return;
    }

    int indexLow = polycone.GetSectionIndex(point.z()-kTolerance);
    int indexHigh = polycone.GetSectionIndex(point.z()+kTolerance);
    int index = 0;
 
    if ( indexLow != indexHigh && (indexLow >= 0 )) {
      //we are close to Surface, section has to be identified
      const PolyconeSection& section = polycone.GetSection(indexLow);
      
      bool inside;
      Vector3D<Precision> localp;
      ConeImplementation< translation::kIdentity, rotation::kIdentity, ConeTypes::UniversalCone>::Contains<Backend>(
                *section.fSolid,
                Transformation3D(),
                point - Vector3D<Precision>(0,0,section.fShift),
                localp,
                inside);
      if(!inside){index=indexHigh;}
      else{index=indexLow;}
    
    }
    else{
        index=indexLow;
        if(index<0)index=polycone.GetSectionIndex(point.z());
    }
    if(index < 0 ){distance = 0.; return; }

    Precision totalDistance = 0.;
    Precision dist;
    int increment = (dir.z() > 0) ? 1 : -1;
    if (std::fabs(dir.z()) < kTolerance) increment = 0;


    // What is the relevance of istep?
    int istep = 0;

    do
    {
        const PolyconeSection& section = polycone.GetSection(index);
    
        if ( (totalDistance!=0) || (istep < 2)) {
            pn = point + totalDistance*dir; // point must be shifted, so it could eventually get into another solid
            pn.z() -= section.fShift;
            typename Backend::int_v inside;
            ConeImplementation< translation::kIdentity, rotation::kIdentity,
                ConeTypes::UniversalCone>::Inside<Backend>(
                    *section.fSolid,
                    Transformation3D(),
                    pn,
                    inside);
            if (inside == EInside::kOutside) {
                break;
            }
        }
        else pn.z() -= section.fShift;

        istep++;
    
        ConeImplementation< translation::kIdentity, rotation::kIdentity,
            ConeTypes::UniversalCone>::DistanceToOut<Backend>( *section.fSolid,
                pn, dir, stepMax, dist );

        //Section Surface case
        if(std::fabs(dist) < 0.5*kTolerance) {
            int index1 = index;
            if(( index > 0) && ( index < polycone.GetNSections()-1 )){
                index1 += increment;}
            else{
                if((index == 0) && ( increment > 0 ))
                    index1 += increment;
                if((index == polycone.GetNSections()-1) && (increment<0))
                    index1 += increment;
            }

            Vector3D<Precision> pte = point+(totalDistance+dist)*dir;
            const PolyconeSection& section1 = polycone.GetSection(index1);
            bool inside1;
            pte.z() -= section1.fShift;
            Vector3D<Precision> localp;
            ConeImplementation< translation::kIdentity, rotation::kIdentity,
                    ConeTypes::UniversalCone>::Contains<Backend>( *section1.fSolid,
                Transformation3D(), pte, localp, inside1 );
            if (!inside1) {
                break;
            }
        } // end if surface case

    /* this was on the master:
        Vector3D<Precision> pte = point+(totalDistance+dist)*dir;
        const PolyconeSection& section1 = polycone.GetSection(index1);
        bool inside1;
        pte.z() -= section1.fShift;
        Vector3D<Precision> localp;
        ConeImplementation< translation::kIdentity, rotation::kIdentity, ConeTypes::UniversalCone>::Contains<Backend>(
                *section1.fSolid,
                Transformation3D(),
                pte,
                localp,
                inside1);
        if ( (!inside1) || (increment == 0) )
        {
         break;
        }
    }
   
    totalDistance += dist;
    index += increment;
    */
   
        totalDistance += dist;
        index += increment;
  }
  while ( increment!=0 && index >= 0 && index < polycone.GetNSections());
 
  distance=totalDistance;
 
  return;
}

    

    template <class Backend>
    VECGEOM_CUDA_HEADER_BOTH
    VECGEOM_INLINE
    static void SafetyToIn(UnplacedPolycone const &polycone,
                           Transformation3D const &transformation,
                           Vector3D<typename Backend::precision_v> const &point,
                           typename Backend::precision_v &safety) {

    Vector3D<typename Backend::precision_v> p = transformation.Transform<transCodeT,rotCodeT>(point);
   
    int index = polycone.GetSectionIndex(p.z());
    bool needZ = false;
    if(index < 0)
      {needZ = true;
       if(index == -1) index = 0;
       if(index == -2) index = polycone.GetNSections()-1;
      }
    Precision minSafety=0 ;//= SafetyFromOutsideSection(index, p);
    PolyconeSection const & sec = polycone.GetSection(index);
    // safety to current segment
    if(needZ)
      {
      safety = ConeImplementation< translation::kIdentity, rotation::kIdentity, ConeTypes::UniversalCone>::SafetyToInUSOLIDS<Backend,false>(
                *sec.fSolid,
                Transformation3D(),
                p - Vector3D<Precision>(0,0,sec.fShift));
      return;
      }
    else
       safety = ConeImplementation< translation::kIdentity, rotation::kIdentity, ConeTypes::UniversalCone>::SafetyToInUSOLIDS<Backend,true>(
                *sec.fSolid,
                Transformation3D(),
                p - Vector3D<Precision>(0,0,sec.fShift));
      
    if (safety < kTolerance) return ;
    minSafety=safety;
    //std::cout<<" entering STI minSaf="<<minSafety<<" index="<<index<<std::endl;
    Precision zbase = polycone.fZs[index + 1];
    // going right
    for (int i = index + 1; i < polycone.GetNSections(); ++i)
    {
     Precision dz = polycone.fZs[i] - zbase;
     if (dz >= minSafety) break;
     
     PolyconeSection const & sec = polycone.GetSection(i);
     safety = ConeImplementation< translation::kIdentity, rotation::kIdentity, ConeTypes::UniversalCone>::SafetyToInUSOLIDS<Backend,false>(
                *sec.fSolid,
                Transformation3D(),
                p - Vector3D<Precision>(0,0,sec.fShift));
     //std::cout<<"Din="<<safety<<" i="<<i<<"dz="<<dz<<std::endl;
     if (safety < minSafety) minSafety = safety;
    }

    // going left if this is possible
    if (index > 0) {
      zbase = polycone.fZs[index - 1];
      for (int i = index - 1; i >= 0; --i) {
        Precision dz = zbase - polycone.fZs[i];
        if (dz >= minSafety)
          break;
        PolyconeSection const &sec = polycone.GetSection(i);
        safety =
            ConeImplementation<translation::kIdentity, rotation::kIdentity,
                               ConeTypes::UniversalCone>::SafetyToInUSOLIDS<Backend, false>(*sec.fSolid,
                                                                                            Transformation3D(),
                                                                                            p - Vector3D<Precision>(
                                                                                                    0, 0, sec.fShift));
        // std::cout<<"Din-1="<<safety<<" i="<<i<<"dz="<<dz<<std::endl;
        if (safety < minSafety)
          minSafety = safety;
      }
    }
    safety = minSafety;
    return ;
    /*
     if (!aAccurate)
    return enclosingCylinder->SafetyFromOutside(p);

  int index = GetSection(p.z);
  double minSafety = SafetyFromOutsideSection(index, p);
  if (minSafety < 1e-6) return minSafety;

  double zbase = fZs[index + 1];
  for (int i = index + 1; i <= fMaxSection; ++i)
  {
    double dz = fZs[i] - zbase;
    if (dz >= minSafety) break;
    double safety = SafetyFromOutsideSection(i, p);
    if (safety < minSafety) minSafety = safety;
  }

  zbase = fZs[index - 1];
  for (int i = index - 1; i >= 0; --i)
  {
    double dz = zbase - fZs[i];
    if (dz >= minSafety) break;
    double safety = SafetyFromOutsideSection(i, p);
    if (safety < minSafety) minSafety = safety;
  }
  return minSafety;
    */
    }

    template <class Backend>
    VECGEOM_CUDA_HEADER_BOTH
    VECGEOM_INLINE
    static void SafetyToOut(UnplacedPolycone const &polycone,
                            Vector3D<typename Backend::precision_v> const &point,
                            typename Backend::precision_v &safety) {

    
 int index = polycone.GetSectionIndex(point.z());
   
    if (index < 0 ){ safety=0;return ;}
  
  
  
  PolyconeSection const & sec = polycone.GetSection(index);
  Vector3D<typename Backend::precision_v> p = point - Vector3D<Precision>(0,0,sec.fShift);
  safety = ConeImplementation< translation::kIdentity, rotation::kIdentity, ConeTypes::UniversalCone>::SafetyToOutUSOLIDS<Backend,false>(
       *sec.fSolid,p);
  Precision minSafety =safety;
  //std::cout<<"Sout0="<<safety<<" index="<<index<<std::endl;
  if (minSafety == kInfinity) {safety = 0.;return ;}
  if (minSafety < kTolerance) {safety = 0.; return ;}

  Precision zbase = polycone.fZs[index + 1];
  for (int i = index + 1; i < polycone.GetNSections(); ++i)
  {
    Precision dz = polycone.fZs[i] - zbase;
    if (dz >= minSafety) break;
    PolyconeSection const & sec = polycone.GetSection(i);
    p = point - Vector3D<Precision>(0,0,sec.fShift);
    safety = ConeImplementation< translation::kIdentity, rotation::kIdentity, ConeTypes::UniversalCone>::SafetyToInUSOLIDS<Backend,false>(
        *sec.fSolid,  Transformation3D(),p);
    if(safety < minSafety)minSafety =safety;
    // std::cout<<"Sout+1="<<safety<<" i="<<i<<"dz="<<dz<<std::endl;
  }

  if (index > 0)
  {
    zbase = polycone.fZs[index - 1];
    for (int i = index - 1; i >= 0; --i)
    {
    Precision dz = zbase - polycone.fZs[i];
    if (dz >= minSafety) break;
    PolyconeSection const & sec = polycone.GetSection(i);
    p = point - Vector3D<Precision>(0,0,sec.fShift);
    safety = ConeImplementation< translation::kIdentity, rotation::kIdentity, ConeTypes::UniversalCone>::SafetyToInUSOLIDS<Backend,false>(
         *sec.fSolid,  Transformation3D(),p);
    if(safety < minSafety)minSafety =safety;
    //std::cout<<"Sout-1="<<safety<<" i="<<i<<"dz="<<dz<<std::endl;
    }
  }

  safety=minSafety;
  return;
}



}; // end PolyconeImplementation

}} // end namespace


#endif /* POLYCONEIMPLEMENTATION_H_ */
