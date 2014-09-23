
/// @file SphereImplementation.h
/// @author Raman Sehgal (raman.sehgal@cern.ch)

#ifndef VECGEOM_VOLUMES_KERNEL_SPHEREIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_SPHEREIMPLEMENTATION_H_
#include "math.h"
#include "base/Global.h"

#include "base/Transformation3D.h"
#include "volumes/kernel/BoxImplementation.h"
#include "volumes/UnplacedSphere.h"
#include "base/Vector3D.h"
#include "volumes/kernel/GenericKernels.h"
#include <iomanip>
#include <Vc/Vc>
//#include "volumes/SphereUtilities.h"
namespace VECGEOM_NAMESPACE { 

template <TranslationCode transCodeT, RotationCode rotCodeT>
struct SphereImplementation {

template <class Backend>
VECGEOM_CUDA_HEADER_BOTH
VECGEOM_INLINE    
static typename Backend::precision_v fabs(typename Backend::precision_v &v)
  {
      typedef typename Backend::precision_v Double_t;
      Double_t mone(-1.);
      Double_t ret(0);
      MaskedAssign( (v<0), mone*v , &ret );
      return ret;
  }

    
template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void UnplacedContains(
      UnplacedSphere const &unplaced,
      Vector3D<typename Backend::precision_v> const &point,
      typename Backend::bool_v &inside);//{}
    
template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void Contains(
      UnplacedSphere const &unplaced,
      Transformation3D const &transformation,
      Vector3D<typename Backend::precision_v> const &point,
      Vector3D<typename Backend::precision_v> &localPoint,
      typename Backend::bool_v &inside);//{}

template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void Inside(UnplacedSphere const &unplaced,
                     Transformation3D const &transformation,
                     Vector3D<typename Backend::precision_v> const &point,
                     typename Backend::inside_v &inside);//{}

template <typename Backend, bool ForInside>
VECGEOM_CUDA_HEADER_BOTH
VECGEOM_INLINE
static void GenericKernelForContainsAndInside(
    UnplacedSphere const &unplaced,
    Vector3D<typename Backend::precision_v> const &localPoint,
    typename Backend::bool_v &completelyinside,
    typename Backend::bool_v &completelyoutside);// {}

template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void DistanceToIn(
      UnplacedSphere const &unplaced,
      Transformation3D const &transformation,
      Vector3D<typename Backend::precision_v> const &point,
      Vector3D<typename Backend::precision_v> const &direction,
      typename Backend::precision_v const &stepMax,
      typename Backend::precision_v &distance){}

template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void DistanceToOut(
      UnplacedSphere const &unplaced,
      Vector3D<typename Backend::precision_v> const &point,
      Vector3D<typename Backend::precision_v> const &direction,
      typename Backend::precision_v const &stepMax,
      typename Backend::precision_v &distance){}

template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void SafetyToIn(UnplacedSphere const &unplaced,
                         Transformation3D const &transformation,
                         Vector3D<typename Backend::precision_v> const &point,
                         typename Backend::precision_v &safety);//{}

 template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void SafetyToOut(UnplacedSphere const &unplaced,
                          Vector3D<typename Backend::precision_v> const &point,
                          typename Backend::precision_v &safety);//{}
 

template <typename Backend>
VECGEOM_CUDA_HEADER_BOTH
VECGEOM_INLINE
static void ContainsKernel(
    UnplacedSphere const &unplaced,
    Vector3D<typename Backend::precision_v> const &localPoint,
    typename Backend::bool_v &inside);//{}
  

template <class Backend>
VECGEOM_CUDA_HEADER_BOTH
VECGEOM_INLINE
static void InsideKernel(
    UnplacedSphere const &unplaced,
    Vector3D<typename Backend::precision_v> const &point,
    typename Backend::inside_v &inside);// {}

 
template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void DistanceToInKernel(
      UnplacedSphere const &unplaced,
      Vector3D<typename Backend::precision_v> const &point,
      Vector3D<typename Backend::precision_v> const &direction,
      typename Backend::precision_v const &stepMax,
      typename Backend::precision_v &distance){}


template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void DistanceToOutKernel(
      UnplacedSphere const &unplaced,
      Vector3D<typename Backend::precision_v> const &point,
      Vector3D<typename Backend::precision_v> const &direction,
      typename Backend::precision_v const &stepMax,
      typename Backend::precision_v &distance){}


template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void SafetyToInKernel(UnplacedSphere const &unplaced,
                         Vector3D<typename Backend::precision_v> const &point,
                         typename Backend::precision_v &safety);//{}
  


  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void SafetyToOutKernel(UnplacedSphere const &unplaced,
                          Vector3D<typename Backend::precision_v> const &point,
                          typename Backend::precision_v &safety);//{}
  
  
  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void Normal(
       UnplacedSphere const &unplaced,
       Vector3D<typename Backend::precision_v> const &point,
       Vector3D<typename Backend::precision_v> &normal,
       typename Backend::bool_v &valid );//{}
  
  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void NormalKernel(
       UnplacedSphere const &unplaced,
       Vector3D<typename Backend::precision_v> const &point,
       Vector3D<typename Backend::precision_v> &normal,
       typename Backend::bool_v &valid );//{}
  
  
  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void ApproxSurfaceNormalKernel(
       UnplacedSphere const &unplaced,
       Vector3D<typename Backend::precision_v> const &point,
       Vector3D<typename Backend::precision_v> &normal);
  
};

template <TranslationCode transCodeT, RotationCode rotCodeT>
template <typename Backend>
VECGEOM_CUDA_HEADER_BOTH
void SphereImplementation<transCodeT, rotCodeT>::Normal(
       UnplacedSphere const &unplaced,
       Vector3D<typename Backend::precision_v> const &point,
       Vector3D<typename Backend::precision_v> &normal,
       typename Backend::bool_v &valid ){

    NormalKernel<Backend>(unplaced, point, normal, valid);
}

template <TranslationCode transCodeT, RotationCode rotCodeT>
template <typename Backend>
VECGEOM_CUDA_HEADER_BOTH
void SphereImplementation<transCodeT, rotCodeT>::ApproxSurfaceNormalKernel(
       UnplacedSphere const &unplaced,
       Vector3D<typename Backend::precision_v> const &point,
       Vector3D<typename Backend::precision_v> &norm){

    //std::cout<<"Entered AppoxNormalKernel"<<std::endl;
        
    typedef typename Backend::precision_v Double_t;
    typedef typename Backend::bool_v      Bool_t;
    
    Double_t kNRMin(0.), kNRMax(1.), kNSPhi(2.), kNEPhi(3.), kNSTheta(4.), kNETheta(5.);
    Double_t side(10.);
 
    Double_t rho, rho2, radius;
    Double_t distRMin, distRMax, distSPhi, distEPhi, distSTheta, distETheta, distMin;
    Double_t zero(0.),mone(-1.);
    Double_t temp=zero;
    
    Vector3D<Double_t> localPoint;
    localPoint = point;
    
    Double_t radius2=localPoint.Mag2();
    radius = Sqrt(radius2);
    rho2 = localPoint.x() * localPoint.x() + localPoint.y() * localPoint.y();
    rho = Sqrt(rho2);
    
    Double_t fRmax = unplaced.GetOuterRadius();
    Double_t fRmin = unplaced.GetInnerRadius();
    Double_t fSPhi = unplaced.GetStartPhiAngle();
    Double_t fDPhi = unplaced.GetDeltaPhiAngle();
    Double_t ePhi = fSPhi + fDPhi;
    Double_t fSTheta = unplaced.GetStartThetaAngle();
    Double_t fDTheta = unplaced.GetDeltaThetaAngle();
    Double_t eTheta = fSTheta + fDTheta;
    Double_t pPhi = localPoint.Phi();
    Double_t pTheta(std::atan2(rho,localPoint.z()));
    Double_t sinSPhi = std::sin(fSPhi);
    Double_t cosSPhi = std::cos(fSPhi);
    Double_t sinEPhi = std::sin(ePhi);
    Double_t cosEPhi = std::cos(ePhi);
    Double_t sinSTheta = std::sin(fSTheta);
    Double_t cosSTheta = std::cos(fSTheta);
    Double_t sinETheta = std::sin(eTheta);
    Double_t cosETheta = std::cos(eTheta);
        
    //
    // Distance to r shells
    //
    temp = radius - fRmax;
    MaskedAssign((temp<zero),mone*temp,&temp);
    distRMax = temp;
    
    temp = radius - fRmin;
    MaskedAssign((temp<zero),mone*temp,&temp);
    distRMin = temp;
    Double_t prevDistMin(zero);
    prevDistMin = distMin;
    MaskedAssign( ( (fRmin > zero) && (distRMin < distRMax) ) , distRMin, &distMin );
    MaskedAssign( ( (fRmin > zero) && (distRMin < distRMax) ) , kNRMin, &side ); //ENorm issue : Resolved hopefully
    
    prevDistMin = distMin;
    MaskedAssign( ( (fRmin > zero) && !(distRMin < distRMax) ) , distRMax, &distMin );
    MaskedAssign( ( (fRmin > zero) && !(distRMin < distRMax) ) , kNRMax, &side );//ENorm issue : Resolved hopefully
    
    MaskedAssign( !(fRmin > zero), distRMax, &distMin );
    MaskedAssign( !(fRmin > zero), kNRMax, &side );
    
    //
    // Distance to phi planes
    //
    // Protected against (0,0,z)
    
    MaskedAssign( (pPhi < zero) ,(pPhi+(2*kPi)) , &pPhi);
    
    temp = pPhi - (fSPhi + 2 * kPi);
    MaskedAssign((temp<zero),(mone*temp),&temp);
    MaskedAssign( (!unplaced.IsFullPhiSphere() && (rho>zero) && (fSPhi < zero) ) ,(temp*rho)  ,&distSPhi);
    
    temp = pPhi - fSPhi;
    MaskedAssign((temp<zero),(mone*temp),&temp);
    MaskedAssign( (!unplaced.IsFullPhiSphere() && (rho>zero) && !(fSPhi < zero) ) ,(temp*rho)  ,&distSPhi);
    
    temp=pPhi - fSPhi - fDPhi;
    MaskedAssign((temp<zero),(mone*temp),&temp);
    MaskedAssign((!unplaced.IsFullPhiSphere() && (rho>zero)), temp*rho, &distEPhi); //distEPhi = temp * rho;
    
    
    prevDistMin = distMin;
    //std::cout<<"--------------------------------------------------------------------------------"<<std::endl;
    MaskedAssign( ( !unplaced.IsFullPhiSphere() && (rho>zero) && (distSPhi < distEPhi) && (distSPhi < distMin)),distSPhi ,&distMin );
    MaskedAssign( ( !unplaced.IsFullPhiSphere() && (rho>zero) && (distSPhi < distEPhi) && (distSPhi < prevDistMin)),kNSPhi ,&side ); //CULPRIT
    
    prevDistMin = distMin;
    MaskedAssign( ( !unplaced.IsFullPhiSphere() && (rho>zero) && !(distSPhi < distEPhi) && (distEPhi < distMin)),distEPhi ,&distMin );
    MaskedAssign( ( !unplaced.IsFullPhiSphere() && (rho>zero) && !(distSPhi < distEPhi) && (distEPhi < prevDistMin)),kNEPhi ,&side );
    
    //
    // Distance to theta planes
    //
    temp = pTheta - fSTheta;
    MaskedAssign((temp<zero),(mone*temp),&temp);
    MaskedAssign( ( !unplaced.IsFullThetaSphere() && (radius > zero)) , temp*radius ,&distSTheta);
    
    temp = pTheta - fSTheta - fDTheta;
    MaskedAssign((temp<zero),(mone*temp),&temp);
    MaskedAssign( ( !unplaced.IsFullThetaSphere() && (radius > zero)) , temp*radius ,&distETheta);
    
    prevDistMin = distMin;
    MaskedAssign( ( !unplaced.IsFullThetaSphere() && (radius > zero) && (distSTheta < distETheta) && (distSTheta < distMin)), distSTheta, &distMin); 
    MaskedAssign( ( !unplaced.IsFullThetaSphere() && (radius > zero) && (distSTheta < distETheta) && (distSTheta < prevDistMin)), kNSTheta, &side); 
    //std::cout<<" (distSTheta < distMin) : "<<distSTheta<<" , "<<distETheta<<" :: "<<distMin<<" :: side :"<<side<<std::endl;
    
    prevDistMin = distMin;
    MaskedAssign( ( !unplaced.IsFullThetaSphere() && (radius > zero) && !(distSTheta < distETheta) && (distETheta < distMin)), distETheta, &distMin); 
    MaskedAssign( ( !unplaced.IsFullThetaSphere() && (radius > zero) && !(distSTheta < distETheta) && (distETheta < prevDistMin)), kNETheta, &side);
    // std::cout<<" !(distSTheta < distMin) : "<<distSTheta<<" , "<<distETheta<<" :: "<<distMin<<" :: side :"<<side<<std::endl;
    
    //Switching
    //std::cout<<"Side : "<<side<<std::endl;
    Bool_t done(false);
    done |= (side == kNRMin);
    MaskedAssign( (side == kNRMin), Vector3D<Double_t>(-localPoint.x() / radius, -localPoint.y() / radius, -localPoint.z() / radius),&norm);
   
    if( IsFull(done) )return ;//{  std::cout<<"----- 1 ------"<<std::endl; return;}
    
    done |= (side == kNRMax);
    MaskedAssign( (side == kNRMax),Vector3D<Double_t>(localPoint.x() / radius, localPoint.y() / radius, localPoint.z() / radius),&norm);
    //std::cout<<"----- 2 ------"<<std::endl;
    if( IsFull(done) )return ;//{  std::cout<<"----- 2 ------"<<std::endl; return;}
    
    done |= (side == kNSPhi);
    MaskedAssign( (side == kNSPhi),Vector3D<Double_t>(sinSPhi, -cosSPhi, 0),&norm);
    //std::cout<<"----- 3 ------"<<std::endl;
    if( IsFull(done) )return ;//{  std::cout<<"----- 3 ------"<<std::endl; return;}
    
    done |= (side == kNEPhi);
    MaskedAssign( (side == kNEPhi),Vector3D<Double_t>(-sinEPhi, cosEPhi, 0),&norm);
    //std::cout<<"----- 4 ------"<<std::endl;
    if( IsFull(done) )return ;//{  std::cout<<"----- 4 ------"<<std::endl; return;}
    
    done |= (side == kNSTheta);
    MaskedAssign( (side == kNSTheta),Vector3D<Double_t>(-cosSTheta * std::cos(pPhi), -cosSTheta * std::sin(pPhi),sinSTheta),&norm);
    //std::cout<<"----- 5 ------"<<std::endl;
    if( IsFull(done) )return ;//{  std::cout<<"----- 5 ------"<<std::endl; return;}
    
    done |= (side == kNETheta);
    MaskedAssign( (side == kNETheta),Vector3D<Double_t>(cosETheta * std::cos(pPhi), cosETheta * std::sin(pPhi), sinETheta),&norm);
    //std::cout<<"----- 6 ------"<<std::endl;
    if( IsFull(done) )return ;//{  std::cout<<"----- 6 ------"<<std::endl; return;}
    
    /*
    switch (side)
  {
    case kNRMin:      // Inner radius
      norm = Vector3D<Double_t>(-localPoint.x() / radius, -localPoint.y() / radius, -localPoint.z() / radius);
      break;
    case kNRMax:      // Outer radius
      norm = Vector3D<Double_t>(localPoint.x() / radius, localPoint.y() / radius, localPoint.z() / radius);
      break;
    case kNSPhi:
      norm = Vector3D<Double_t>(sinSPhi, -cosSPhi, 0);
      break;
    case kNEPhi:
      norm = Vector3D<Double_t>(-sinEPhi, cosEPhi, 0);
      break;
    case kNSTheta:
      norm = Vector3D<Double_t>(-cosSTheta * std::cos(pPhi), -cosSTheta * std::sin(pPhi),sinSTheta);
      break;
    case kNETheta:
      norm = Vector3D<Double_t>(cosETheta * std::cos(pPhi), cosETheta * std::sin(pPhi), sinETheta);
      break;
    default:          // Should never reach this case ...

        std::cout<<"Undefined side for valid surface normal to solid."<<std::endl;
      break;
  }
  */
    
}

template <TranslationCode transCodeT, RotationCode rotCodeT>
template <typename Backend>
VECGEOM_CUDA_HEADER_BOTH
void SphereImplementation<transCodeT, rotCodeT>::NormalKernel(
       UnplacedSphere const &unplaced,
       Vector3D<typename Backend::precision_v> const &point,
       Vector3D<typename Backend::precision_v> &normal,
       typename Backend::bool_v &valid ){

    typedef typename Backend::precision_v Double_t;
    typedef typename Backend::bool_v      Bool_t;

    Vector3D<Double_t> localPoint;
    localPoint = point;
    
    Double_t radius2=localPoint.Mag2();
    Double_t radius=Sqrt(radius2);
    Double_t rho2 = localPoint.x() * localPoint.x() + localPoint.y() * localPoint.y();
    Double_t rho = Sqrt(rho2);
    Double_t fRmax = unplaced.GetOuterRadius();
    Double_t fRmin = unplaced.GetInnerRadius();
    Double_t fSPhi = unplaced.GetStartPhiAngle();
    Double_t fDPhi = unplaced.GetDeltaPhiAngle();
    Double_t ePhi = fSPhi + fDPhi;
    Double_t fSTheta = unplaced.GetStartThetaAngle();
    Double_t fDTheta = unplaced.GetDeltaThetaAngle();
    Double_t eTheta = fSTheta + fDTheta;
    Double_t pPhi = localPoint.Phi();
    Double_t pTheta(std::atan2(rho,localPoint.z()));
    Double_t sinSPhi = std::sin(fSPhi);
    Double_t cosSPhi = std::cos(fSPhi);
    Double_t sinEPhi = std::sin(ePhi);
    Double_t cosEPhi = std::cos(ePhi);
    Double_t sinSTheta = std::sin(fSTheta);
    Double_t cosSTheta = std::cos(fSTheta);
    Double_t sinETheta = std::sin(eTheta);
    Double_t cosETheta = std::cos(eTheta);
    
    Precision fRminTolerance = unplaced.GetFRminTolerance();
    Precision kAngTolerance = unplaced.GetAngTolerance();
    Precision halfAngTolerance = (0.5 * kAngTolerance);
    
    Double_t distSPhi(kInfinity),distSTheta(kInfinity);
    Double_t distEPhi(kInfinity),distETheta(kInfinity);
    Double_t distRMax(kInfinity); 
    Double_t distRMin(kInfinity);
    
    Vector3D<Double_t> nR, nPs, nPe, nTs, nTe, nZ(0., 0., 1.);
    Vector3D<Double_t> norm, sumnorm(0., 0., 0.);
    
    Bool_t fFullPhiSphere = unplaced.IsFullPhiSphere();
    //std::cout<<"Full Phi Sphere Check : "<<fFullPhiSphere<<std::endl;
    Bool_t fFullThetaSphere = unplaced.IsFullThetaSphere();
    //std::cout<<"Full Theta Sphere Check : "<<fFullThetaSphere<<std::endl;

    Double_t zero(0.);
    Double_t mone(-1.);
    
    Double_t temp=0;
    temp=radius - fRmax;
    MaskedAssign((temp<zero),mone*temp,&temp);
    //distRMax = fabs<Backend>(temp);
    distRMax = temp;
    
    distRMin = 0;
    temp = radius - fRmin;
    MaskedAssign((temp<zero),mone*temp,&temp);
    //MaskedAssign( (fRmin > 0) , fabs<Backend>(temp), &distRMin);
    MaskedAssign( (fRmin > 0) , temp, &distRMin);
    
    MaskedAssign( ((rho>zero) && !(unplaced.IsFullSphere()) && (pPhi < fSPhi - halfAngTolerance)) , (pPhi+(2*kPi)), &pPhi);
    MaskedAssign( ((rho>zero) && !(unplaced.IsFullSphere()) && (pPhi > ePhi + halfAngTolerance)) , (pPhi-(2*kPi)), &pPhi);
    
    //Phi Stuff
    temp = pPhi - fSPhi;
    MaskedAssign((temp<zero),mone*temp,&temp);
    //MaskedAssign( (!unplaced.IsFullPhiSphere() && (rho>zero) ), (fabs<Backend>(temp)),&distSPhi);
    MaskedAssign( (!unplaced.IsFullPhiSphere() && (rho>zero) ), temp,&distSPhi);
    temp = pPhi - ePhi;
    MaskedAssign((temp<zero),mone*temp,&temp);
    //MaskedAssign( (!unplaced.IsFullPhiSphere() && (rho>zero) ), (fabs<Backend>(temp)),&distEPhi);
    MaskedAssign( (!unplaced.IsFullPhiSphere() && (rho>zero) ), temp,&distEPhi);
    
    MaskedAssign( (!unplaced.IsFullPhiSphere() && (!fRmin) ), zero ,&distSPhi);
    MaskedAssign( (!unplaced.IsFullPhiSphere() && (!fRmin) ), zero ,&distEPhi);
    MaskedAssign( (!unplaced.IsFullPhiSphere()), Vector3D<Double_t>(sinSPhi, -cosSPhi, 0),&nPs );
    MaskedAssign( (!unplaced.IsFullPhiSphere()), Vector3D<Double_t>(-sinEPhi, cosEPhi, 0),&nPe );
    
    //Theta Stuff
    temp = pTheta - fSTheta;
    MaskedAssign((temp<zero),mone*temp,&temp);
    //MaskedAssign( ( !unplaced.IsFullThetaSphere() && (rho>zero) ) , fabs<Backend>(temp) , &distSTheta ) ;
    MaskedAssign( ( !unplaced.IsFullThetaSphere() && (rho>zero) ) , temp , &distSTheta ) ;
    
    temp = pTheta - eTheta;
    MaskedAssign((temp<zero),mone*temp,&temp);
    //MaskedAssign( ( !unplaced.IsFullThetaSphere() && (rho>zero) ) , fabs<Backend>(temp) , &distETheta ) ;
    MaskedAssign( ( !unplaced.IsFullThetaSphere() && (rho>zero) ) , temp , &distETheta ) ;
    
    MaskedAssign( ( !unplaced.IsFullThetaSphere() && (rho>zero) ) , Vector3D<Double_t>(-cosSTheta * localPoint.x() / rho , -cosSTheta * localPoint.y() / rho, sinSTheta ) , &nTs ) ;
    MaskedAssign( ( !unplaced.IsFullThetaSphere() && (rho) ) , Vector3D<Double_t>(cosETheta * localPoint.x() / rho , cosETheta * localPoint.y() / rho , -sinETheta) , &nTe ) ;
    
    MaskedAssign( ( !unplaced.IsFullThetaSphere() && (!fRmin) && (fSTheta)) , zero , &distSTheta );
    MaskedAssign( ( !unplaced.IsFullThetaSphere() && (!fRmin) && (fSTheta)) , Vector3D<Double_t>(0., 0., -1.) , &nTs );
    MaskedAssign( ( !unplaced.IsFullThetaSphere() && (!fRmin) && (eTheta < kPi)) , zero , &distETheta );
    MaskedAssign( ( !unplaced.IsFullThetaSphere() && (!fRmin) && (eTheta < kPi)) , Vector3D<Double_t>(0., 0., 1.) , &nTe );
    
    MaskedAssign( (radius), Vector3D<Double_t>(localPoint.x() / radius, localPoint.y() / radius, localPoint.z() / radius) ,&nR);
    
    
    Double_t noSurfaces(0);
    Double_t halfCarTolerance(0.5 * 1e-9);
    MaskedAssign((distRMax <= halfCarTolerance) , noSurfaces+1 ,&noSurfaces);
    MaskedAssign((distRMax <= halfCarTolerance) , (sumnorm+nR) ,&sumnorm);
    
    MaskedAssign((fRmin && (distRMin <= halfCarTolerance)) , noSurfaces+1 ,&noSurfaces);
    MaskedAssign((fRmin && (distRMin <= halfCarTolerance)) , (sumnorm-nR) ,&sumnorm);
    
    MaskedAssign( (!unplaced.IsFullPhiSphere() && (distSPhi <= halfAngTolerance)) , noSurfaces+1 ,&noSurfaces);
    MaskedAssign( (!unplaced.IsFullPhiSphere() && (distSPhi <= halfAngTolerance)) , (sumnorm+nPs) ,&sumnorm);
    MaskedAssign( (!unplaced.IsFullPhiSphere() && (distEPhi <= halfAngTolerance)) , noSurfaces+1 ,&noSurfaces);
    MaskedAssign( (!unplaced.IsFullPhiSphere() && (distEPhi <= halfAngTolerance)) , (sumnorm+nPe) ,&sumnorm);
    
    MaskedAssign( (!unplaced.IsFullThetaSphere() && (distSTheta <= halfAngTolerance) && (fSTheta > zero)), noSurfaces+1 ,&noSurfaces);
    MaskedAssign( (!unplaced.IsFullThetaSphere() && (distSTheta <= halfAngTolerance) && (fSTheta > zero) && ((radius <= halfCarTolerance) && fFullPhiSphere)), (sumnorm+nZ) ,&sumnorm);
    MaskedAssign( (!unplaced.IsFullThetaSphere() && (distSTheta <= halfAngTolerance) && (fSTheta > zero) && !((radius <= halfCarTolerance) && fFullPhiSphere)), (sumnorm+nTs) ,&sumnorm);
    MaskedAssign( (!unplaced.IsFullThetaSphere() && (distETheta <= halfAngTolerance) && (eTheta < kPi)), noSurfaces+1 ,&noSurfaces);
    MaskedAssign( (!unplaced.IsFullThetaSphere() && (distETheta <= halfAngTolerance) && (eTheta < kPi) && ((radius <= halfCarTolerance) && fFullPhiSphere)), (sumnorm-nZ) ,&sumnorm);
    MaskedAssign( (!unplaced.IsFullThetaSphere() && (distETheta <= halfAngTolerance) && (eTheta < kPi) && !((radius <= halfCarTolerance) && fFullPhiSphere)), (sumnorm+nTe) ,&sumnorm);
    MaskedAssign( (!unplaced.IsFullThetaSphere() && (distETheta <= halfAngTolerance) && (eTheta < kPi) && (sumnorm.z() == zero)), (sumnorm+nZ) ,&sumnorm);
    
    //Now considering case of ApproxSurfaceNormal
    if(noSurfaces == 0)
        ApproxSurfaceNormalKernel<Backend>(unplaced,point,norm);
    
    MaskedAssign((noSurfaces == 1),sumnorm,&norm);
    MaskedAssign((!(noSurfaces == 1) && (noSurfaces !=0 )),(sumnorm*1./sumnorm.Mag()),&norm);
    MaskedAssign(true,norm,&normal);
    
    
    valid = (noSurfaces>zero);
   
}

template <TranslationCode transCodeT, RotationCode rotCodeT>
template <typename Backend>
VECGEOM_CUDA_HEADER_BOTH
void SphereImplementation<transCodeT, rotCodeT>::Contains(UnplacedSphere const &unplaced,
      Transformation3D const &transformation,
      Vector3D<typename Backend::precision_v> const &point,
      Vector3D<typename Backend::precision_v> &localPoint,
      typename Backend::bool_v &inside){

    localPoint = transformation.Transform<transCodeT, rotCodeT>(point);
    UnplacedContains<Backend>(unplaced, localPoint, inside);
}


template <TranslationCode transCodeT, RotationCode rotCodeT>
template <typename Backend>
VECGEOM_CUDA_HEADER_BOTH
void SphereImplementation<transCodeT, rotCodeT>::UnplacedContains(
      UnplacedSphere const &unplaced,
      Vector3D<typename Backend::precision_v> const &point,
      typename Backend::bool_v &inside){

      ContainsKernel<Backend>(unplaced, point, inside);
}    


template <TranslationCode transCodeT, RotationCode rotCodeT>
template <typename Backend>
VECGEOM_CUDA_HEADER_BOTH
void SphereImplementation<transCodeT, rotCodeT>::ContainsKernel(UnplacedSphere const &unplaced,
    Vector3D<typename Backend::precision_v> const &localPoint,
    typename Backend::bool_v &inside) {

  typedef typename Backend::bool_v Bool_t;
  Bool_t unused;
  Bool_t outside;
  GenericKernelForContainsAndInside<Backend, false>(unplaced,
    localPoint, unused, outside);
  inside=!outside;
}  

template <TranslationCode transCodeT, RotationCode rotCodeT>
template <typename Backend, bool ForInside>
VECGEOM_CUDA_HEADER_BOTH
void SphereImplementation<transCodeT, rotCodeT>::GenericKernelForContainsAndInside(
UnplacedSphere const &unplaced,
    Vector3D<typename Backend::precision_v> const &localPoint,
    typename Backend::bool_v &completelyinside,
    typename Backend::bool_v &completelyoutside) {

    typedef typename Backend::precision_v Double_t;
    typedef typename Backend::bool_v      Bool_t;	

    Precision fRmin = unplaced.GetInnerRadius();
    Precision fRminTolerance = unplaced.GetFRminTolerance();
    Precision kAngTolerance = unplaced.GetAngTolerance();
    Precision halfAngTolerance = (0.5 * kAngTolerance);
    Precision fRmax = unplaced.GetOuterRadius();
        
    Double_t rad2 = localPoint.Mag2();
    Double_t tolRMin = fRmin + (0.5 * fRminTolerance); //rMinPlus
    Double_t tolRMin2 = tolRMin * tolRMin;
    Double_t tolRMax = fRmax - (0.5 * fRminTolerance); //rMaxMinus
    Double_t tolRMax2 = tolRMax * tolRMax;
     
    // Check radial surfaces
    //Radial check for GenericKernel Start
    completelyinside = (rad2 <= tolRMax*tolRMax) && (rad2 >= tolRMin*tolRMin);
    
    tolRMin = fRmin - (0.5 * fRminTolerance); //rMinMinus
    tolRMax = fRmax + (0.5 * fRminTolerance); //rMaxPlus
    
    completelyoutside = (rad2 <= tolRMin*tolRMin) || (rad2 >= tolRMax*tolRMax);
    if( IsFull(completelyoutside) )return;
    //return; //Later on remove it and should be only at the end when checks for PHI and THETA finishes
    
    //Radial Check for GenericKernel Over
    
    
    
    Double_t pPhi = localPoint.Phi();
    Double_t fSPhi(unplaced.GetStartPhiAngle());
    Double_t fDPhi(unplaced.GetDeltaPhiAngle());
    Double_t ePhi = fSPhi+fDPhi;
    
    //*******************************
    //Very important but needs to understand
    MaskedAssign((pPhi<(fSPhi - halfAngTolerance)),pPhi+(2.*kPi),&pPhi);
    MaskedAssign((pPhi>(ePhi + halfAngTolerance)),pPhi-(2.*kPi),&pPhi);
    
    //*******************************
    
    Double_t tolAngMin = fSPhi + halfAngTolerance;
    Double_t tolAngMax = ePhi - halfAngTolerance;
    
    // Phi boundaries  : Do not check if it has no phi boundary!
    if(!unplaced.IsFullPhiSphere()) //Later May be done using MaskedAssign
    {
    completelyinside &= (pPhi <= tolAngMax) && (pPhi >= tolAngMin);
    
    tolAngMin = fSPhi - halfAngTolerance;
    tolAngMax = ePhi + halfAngTolerance;
    
    //std::cout<<std::setprecision(20)<<tolAngMin<<" : "<<tolAngMax<<" : "<<pPhi<<std::endl;
    completelyoutside |= (pPhi < tolAngMin) || (pPhi > tolAngMax);
    if( IsFull(completelyoutside) )return;
    }
    //Phi Check for GenericKernel Over
         
    // Theta bondaries
    //Double_t pTheta = localPoint.Theta();
    Double_t pTheta = ATan2(Sqrt(localPoint.x()*localPoint.x() + localPoint.y()*localPoint.y()), localPoint.z()); //This needs to be implemented in Vector3D.h as Theta() function
    Double_t fSTheta(unplaced.GetStartThetaAngle());
    Double_t fDTheta(unplaced.GetDeltaThetaAngle());
    Double_t eTheta = fSTheta + fDTheta;
    
    tolAngMin = fSTheta + halfAngTolerance;
    tolAngMax = eTheta - halfAngTolerance;
    
    if(!unplaced.IsFullThetaSphere())
    {
        completelyinside &= (pTheta <= tolAngMax) && (pTheta >= tolAngMin);
       
        tolAngMin = fSTheta - halfAngTolerance;
        tolAngMax = eTheta + halfAngTolerance;
    
        //std::cout<<std::setprecision(20)<<tolAngMin<<" : "<<tolAngMax<<" : "<<pPhi<<std::endl;
        completelyoutside |= (pTheta < tolAngMin) || (pTheta > tolAngMax);
        if( IsFull(completelyoutside) )return;
        
    }
    
    return;
}


template <TranslationCode transCodeT, RotationCode rotCodeT>
template <typename Backend>
VECGEOM_CUDA_HEADER_BOTH
void SphereImplementation<transCodeT, rotCodeT>::Inside(UnplacedSphere const &unplaced,
                     Transformation3D const &transformation,
                     Vector3D<typename Backend::precision_v> const &point,
                     typename Backend::inside_v &inside){

    InsideKernel<Backend>(unplaced,
                        transformation.Transform<transCodeT, rotCodeT>(point),
                        inside);

}

template <TranslationCode transCodeT, RotationCode rotCodeT>
template <class Backend>
VECGEOM_CUDA_HEADER_BOTH
void SphereImplementation<transCodeT, rotCodeT>::InsideKernel(UnplacedSphere const &unplaced,
    Vector3D<typename Backend::precision_v> const &point,
    typename Backend::inside_v &inside) {
  
  typedef typename Backend::bool_v      Bool_t;
  Bool_t completelyinside, completelyoutside;
  GenericKernelForContainsAndInside<Backend,true>(
      unplaced, point, completelyinside, completelyoutside);
  inside=EInside::kSurface;
  MaskedAssign(completelyoutside, EInside::kOutside, &inside);
  MaskedAssign(completelyinside, EInside::kInside, &inside);
}



template <TranslationCode transCodeT, RotationCode rotCodeT>
template <class Backend>
VECGEOM_CUDA_HEADER_BOTH
VECGEOM_INLINE
void SphereImplementation<transCodeT, rotCodeT>::SafetyToIn(UnplacedSphere const &unplaced,
                         Transformation3D const &transformation,
                         Vector3D<typename Backend::precision_v> const &point,
                         typename Backend::precision_v &safety){
    
    SafetyToInKernel<Backend>(
    unplaced,
    transformation.Transform<transCodeT, rotCodeT>(point),
    safety
  );
}


template <TranslationCode transCodeT, RotationCode rotCodeT>
template <class Backend>
VECGEOM_CUDA_HEADER_BOTH
void SphereImplementation<transCodeT, rotCodeT>::SafetyToInKernel(UnplacedSphere const &unplaced,
                         Vector3D<typename Backend::precision_v> const &point,
                         typename Backend::precision_v &safety){

    typedef typename Backend::precision_v Double_t;
    typedef typename Backend::bool_v      Bool_t;

    Double_t safe=Backend::kZero;
    Double_t zero=Backend::kZero; 

    Vector3D<Double_t> localPoint;
    localPoint=point;

    //General Precalcs
    Double_t rad2    = localPoint.Mag2();
    Double_t rad = Sqrt(rad2);
    Double_t rho2 = localPoint.x() * localPoint.x() + localPoint.y() * localPoint.y();
    Double_t rho = std::sqrt(rho2);
    
    //Distance to r shells
    Double_t fRmin = unplaced.GetInnerRadius();
    Double_t safeRMin = unplaced.GetInnerRadius() - rad;
    Double_t safeRMax = rad - unplaced.GetOuterRadius();
    MaskedAssign(((fRmin > 0)&& (safeRMin > safeRMax)),safeRMin,&safety);
    MaskedAssign(((fRmin > 0)&& (!(safeRMin > safeRMax))),safeRMax,&safety);
    MaskedAssign((!(fRmin > 0)),safeRMax,&safety);
    //Distance to r shells over
    
    //Some Precalc
    Double_t fSPhi = unplaced.GetStartPhiAngle();
    Double_t fDPhi = unplaced.GetDeltaPhiAngle();
    Double_t hDPhi = (0.5 * fDPhi);
    Double_t cPhi  = fSPhi + hDPhi;
    Double_t ePhi  = fSPhi + fDPhi;
    Double_t sinCPhi   = std::sin(cPhi);
    Double_t cosCPhi   = std::cos(cPhi);
    Double_t sinSPhi = std::sin(fSPhi);
    Double_t cosSPhi = std::cos(fSPhi);
    Double_t sinEPhi = std::sin(ePhi);
    Double_t cosEPhi = std::cos(ePhi);
    Double_t safePhi = zero;
    
    Double_t mone(-1.);
    
    Double_t cosPsi = (localPoint.x() * cosCPhi + localPoint.y() * sinCPhi) / rho; 
    //
    // Distance to phi extent
    //
    if(!unplaced.IsFullPhiSphere())
    {
        Double_t test1=((localPoint.x() * sinSPhi - localPoint.y() * cosSPhi));
        MaskedAssign((test1<0),mone*test1,&test1); //Facing the issue with abs function of Vc, Its actually giving the absolute value of floor function
        Double_t test2=((localPoint.x() * sinEPhi - localPoint.y() * cosEPhi));
        MaskedAssign((test2<0),mone*test2,&test2);
        MaskedAssign(((cosPsi < std::cos(hDPhi)) && ((localPoint.y() * cosCPhi - localPoint.x() * sinCPhi) <= 0)),test1,&safePhi);
        MaskedAssign(((cosPsi < std::cos(hDPhi)) && !((localPoint.y() * cosCPhi - localPoint.x() * sinCPhi) <= 0)),test2,&safePhi);
        MaskedAssign(((cosPsi < std::cos(hDPhi)) && (safePhi > safety)),safePhi,&safety);
                
    }
    //
    // Distance to Theta extent
    //
    Double_t KPI(kPi);
    Double_t rds = localPoint.Mag();
    Double_t piby2(kPi/2);
    Double_t pTheta = piby2 - asin(localPoint.z() / rds);
    
    MaskedAssign((pTheta<zero),pTheta+KPI,&pTheta);
    
    Double_t fSTheta = unplaced.GetStartThetaAngle();
    Double_t fDTheta = unplaced.GetDeltaThetaAngle();
    Double_t eTheta  = fSTheta + fDTheta;
    Double_t dTheta1 = fSTheta - pTheta;
    Double_t dTheta2 = pTheta-eTheta;
    Double_t sinDTheta1(std::sin(dTheta1));
    Double_t sinDTheta2(std::sin(dTheta2));
    Double_t safeTheta = zero;
    
    if(!unplaced.IsFullThetaSphere())
    {
        MaskedAssign(((dTheta1 > dTheta2) && (dTheta1 >= zero)),(rds * sinDTheta1),&safeTheta);
        MaskedAssign((!(dTheta1 > dTheta2) && (dTheta2 >= zero)),(rds * sinDTheta2),&safeTheta);
        MaskedAssign(((dTheta1 > dTheta2) && (dTheta1 >= zero) && (safety <= safeTheta)),safeTheta,&safety);
        MaskedAssign((!(dTheta1 > dTheta2) && (dTheta2 >= zero) && (safety <= safeTheta)),safeTheta,&safety);
    }
    
    //Last line
    MaskedAssign( (safety < zero) , zero, &safety);
    
}



template <TranslationCode transCodeT, RotationCode rotCodeT>
template <class Backend>
VECGEOM_CUDA_HEADER_BOTH
VECGEOM_INLINE
void SphereImplementation<transCodeT, rotCodeT>::SafetyToOut(UnplacedSphere const &unplaced,
                          Vector3D<typename Backend::precision_v> const &point,
                          typename Backend::precision_v &safety){
    SafetyToOutKernel<Backend>(
    unplaced,
    point,
    safety
  );
}

template <TranslationCode transCodeT, RotationCode rotCodeT>
template <class Backend>
VECGEOM_CUDA_HEADER_BOTH
void SphereImplementation<transCodeT, rotCodeT>::SafetyToOutKernel(UnplacedSphere const &unplaced,
                          Vector3D<typename Backend::precision_v> const &point,
                          typename Backend::precision_v &safety){
    //std::cout<<"Safety to OUT Kernel call"<<std::endl;
    typedef typename Backend::precision_v Double_t;
    typedef typename Backend::bool_v      Bool_t;

    Double_t safe=Backend::kZero;
    Double_t zero=Backend::kZero; 

    Vector3D<Double_t> localPoint;
    localPoint = point;
    
    //General Precalcs
    Double_t rad2    = localPoint.Mag2();
    Double_t rad = Sqrt(rad2);
    Double_t rho2 = localPoint.x() * localPoint.x() + localPoint.y() * localPoint.y();
    Double_t rho = std::sqrt(rho2);
    
    //Distance to r shells
    Double_t fRmin = unplaced.GetInnerRadius();
    Double_t fRmax = unplaced.GetOuterRadius();
    Double_t safeRMin = rad - fRmin;
    Double_t safeRMax = fRmax - rad ;
        
    //
    // Distance to r shells
    //
    MaskedAssign( ( (fRmin > zero) && (safeRMin < safeRMax) ),safeRMin,&safety);
    MaskedAssign( ( (fRmin > zero) && !(safeRMin < safeRMax) ),safeRMax,&safety);
    MaskedAssign( ( !(fRmin > zero) ),safeRMax,&safety);
    
    //
    // Distance to phi extent
    //
    
    //Some Precalc
    Double_t fSPhi = unplaced.GetStartPhiAngle();
    Double_t fDPhi = unplaced.GetDeltaPhiAngle();
    Double_t hDPhi = (0.5 * fDPhi);
    Double_t cPhi  = fSPhi + hDPhi;
    Double_t ePhi  = fSPhi + fDPhi;
    Double_t sinCPhi   = std::sin(cPhi);
    Double_t cosCPhi   = std::cos(cPhi);
    Double_t sinSPhi = std::sin(fSPhi);
    Double_t cosSPhi = std::cos(fSPhi);
    Double_t sinEPhi = std::sin(ePhi);
    Double_t cosEPhi = std::cos(ePhi);
    Double_t safePhi = zero;
    
    Double_t mone(-1.);
    
    if(!unplaced.IsFullPhiSphere())
    {
        Double_t test1 = (localPoint.y() * cosCPhi - localPoint.x() * sinCPhi);
        //Double_t test2 = (localPoint.x() * sinEPhi - localPoint.y() * cosEPhi);
        MaskedAssign( (test1<=zero) ,(mone*(localPoint.x() * sinSPhi - localPoint.y() * cosSPhi)), &safePhi);
        MaskedAssign( (!(test1<=zero)) ,(localPoint.x() * sinEPhi - localPoint.y() * cosEPhi), &safePhi);
        MaskedAssign( (safePhi < safety) ,safePhi , &safety);
        
    }
    
    //
    // Distance to Theta extent
    //
    Double_t KPI(kPi);
    Double_t rds = localPoint.Mag();
    Double_t piby2(kPi/2);
    Double_t pTheta = piby2 - asin(localPoint.z() / rds);
    
    MaskedAssign((pTheta<zero),pTheta+KPI,&pTheta);
    
    Double_t fSTheta = unplaced.GetStartThetaAngle();
    Double_t fDTheta = unplaced.GetDeltaThetaAngle();
    Double_t eTheta  = fSTheta + fDTheta;
    Double_t dTheta1 =  pTheta - fSTheta;
    Double_t dTheta2 = eTheta - pTheta;
    Double_t sinDTheta1(std::sin(dTheta1));
    Double_t sinDTheta2(std::sin(dTheta2));
    Double_t safeTheta = zero;
    
    if(!unplaced.IsFullThetaSphere())
    {
        MaskedAssign(( (dTheta1 < dTheta2) ),(rds * sinDTheta1),&safeTheta);
        MaskedAssign((!(dTheta1 < dTheta2)),(rds * sinDTheta2),&safeTheta);
        MaskedAssign( (safety > safeTheta ),safeTheta,&safety);
    }
    
    
    //Last line
    MaskedAssign( (safety < zero) , zero, &safety);
}


 

/*
template <TranslationCode transCodeT, RotationCode rotCodeT>
template <class Backend>
VECGEOM_CUDA_HEADER_BOTH
void SphereImplementation<transCodeT, rotCodeT>::DistanceToIn(UnplacedSphere const &unplaced,
      Transformation3D const &transformation,
      Vector3D<typename Backend::precision_v> const &point,
      Vector3D<typename Backend::precision_v> const &direction,
      typename Backend::precision_v const &stepMax,
      typename Backend::precision_v &distance){

    DistanceToInKernel<Backend>(
            unplaced,
            transformation.Transform<transCodeT, rotCodeT>(point),
            transformation.TransformDirection<rotCodeT>(direction),
            stepMax,
            distance);
}

template <TranslationCode transCodeT, RotationCode rotCodeT>
template <class Backend>
VECGEOM_CUDA_HEADER_BOTH
void SphereImplementation<transCodeT, rotCodeT>::DistanceToOut(UnplacedSphere const &unplaced,
      Vector3D<typename Backend::precision_v> const &point,
      Vector3D<typename Backend::precision_v> const &direction,
      typename Backend::precision_v const &stepMax,
      typename Backend::precision_v &distance){

    DistanceToOutKernel<Backend>(
    unplaced,
    point,
    direction,
    stepMax,
    distance
  );
}


template <TranslationCode transCodeT, RotationCode rotCodeT>
template <class Backend>
VECGEOM_CUDA_HEADER_BOTH
void SphereImplementation<transCodeT, rotCodeT>::DistanceToInKernel(
      UnplacedSphere const &unplaced,
      Vector3D<typename Backend::precision_v> const &point,
      Vector3D<typename Backend::precision_v> const &direction,
      typename Backend::precision_v const &stepMax,
      typename Backend::precision_v &distance){
    
    typedef typename Backend::precision_v Double_t;
    typedef typename Backend::bool_v      Bool_t;

    Vector3D<Double_t> localPoint;
    localPoint = point;

    Vector3D<Double_t> localDir;
    localDir = direction;

    //General Precalcs
    Double_t rad2 = localPoint.Mag2();
    Double_t rad = Sqrt(rad2);
    Double_t pDotV3d = localPoint.Dot(localDir);

    Double_t radius2 = unplaced.GetRadius() * unplaced.GetRadius();
    Double_t c = rad2 - radius2;
    Double_t d2 = pDotV3d * pDotV3d - c;

    Double_t pos_dot_dir_x = localPoint.x()*localDir.x();
    Double_t pos_dot_dir_y = localPoint.y()*localDir.y();
    Double_t pos_dot_dir_z = localPoint.z()*localDir.z();

    Bool_t done(false);
    distance = kInfinity;
    Double_t zero=Backend::kZero;

    //Is the point Inside
    Bool_t isInside = ((rad < unplaced.GetfRTolI()));
    done |= isInside;
    MaskedAssign( isInside, kInfinity, &distance );
    if( IsFull(done) )return;

    //On the Surface and Moving In
    Bool_t isInsideOuterToleranceAndMovingIn=((rad >= unplaced.GetfRTolI()) && (rad <= unplaced.GetfRTolO()) && (pDotV3d < 0));
    done |= isInsideOuterToleranceAndMovingIn;
    MaskedAssign(isInsideOuterToleranceAndMovingIn,zero,&distance);
    if( IsFull(done) )return;

    //On the Surface and Moving Out
    Bool_t isInsideOuterToleranceAndMovingOut=((rad >= unplaced.GetfRTolI()) && (rad <= unplaced.GetfRTolO()) && (pDotV3d >= 0));
    done |= isInsideOuterToleranceAndMovingOut;
    MaskedAssign(isInsideOuterToleranceAndMovingIn,kInfinity,&distance);
    if( IsFull(done) )return;
    
    //Outside the Surface and Moving In
    Bool_t isOutsideOuterToleranceAndMovingIn=( (rad > unplaced.GetfRTolO()) && (pDotV3d < 0));
    done |= isOutsideOuterToleranceAndMovingIn;
    MaskedAssign(isOutsideOuterToleranceAndMovingIn,(-pDotV3d - Sqrt(d2)),&distance);
    if( IsFull(done) )return;
    
    //Outside the Surface and Moving Out
    Bool_t isOutsideOuterToleranceAndMovingOut=( (rad > unplaced.GetfRTolO()) && (pDotV3d >= 0));
    done |= isOutsideOuterToleranceAndMovingOut;
    MaskedAssign(isOutsideOuterToleranceAndMovingOut,kInfinity,&distance);
    if( IsFull(done) )return;
}

template <TranslationCode transCodeT, RotationCode rotCodeT>
template <class Backend>
VECGEOM_CUDA_HEADER_BOTH
void SphereImplementation<transCodeT, rotCodeT>::DistanceToOutKernel(UnplacedSphere const &unplaced,
      Vector3D<typename Backend::precision_v> const &point,
      Vector3D<typename Backend::precision_v> const &direction,
      typename Backend::precision_v const &stepMax,
      typename Backend::precision_v &distance){

   typedef typename Backend::precision_v Double_t;
    typedef typename Backend::bool_v      Bool_t;

    distance = kInfinity;  
    Double_t zero=Backend::kZero;
    distance = zero;
    
    Vector3D<Double_t> localPoint;
    localPoint = point;
    
    Vector3D<Double_t> localDir;
    localDir =  direction;
    
    //General Precalcs
    Double_t rad2    = localPoint.Mag2();
    Double_t rad = Sqrt(rad2);
    Double_t pDotV3d = localPoint.Dot(localDir);

    Double_t radius2 = unplaced.GetRadius() * unplaced.GetRadius();
    Double_t c = rad2 - radius2;
    Double_t d2 = pDotV3d * pDotV3d - c;
  
    Bool_t done(false);
    distance = kInfinity;
    distance = zero;
    
    //checking if the point is outside
    Double_t tolRMax = unplaced.GetfRTolO();
    Double_t tolRMax2 = tolRMax * tolRMax;
    Bool_t isOutside = ( rad2 > tolRMax2);
    done|= isOutside;
    if( IsFull(done) )return;

    Bool_t isInsideAndWithinOuterTolerance = ((rad <= tolRMax) && (c < (kTolerance * unplaced.GetRadius())));
    Bool_t isInsideAndOnTolerantSurface = ((c > (-2*kTolerance*unplaced.GetRadius())) && ( (pDotV3d >= 0) || (d2 < 0) ));

    Bool_t onSurface=(isInsideAndWithinOuterTolerance && isInsideAndOnTolerantSurface );
    MaskedAssign(onSurface , zero, &distance);
    done|=onSurface;
    if( IsFull(done) )return;

    Bool_t notOnSurface=(isInsideAndWithinOuterTolerance && !isInsideAndOnTolerantSurface );
    MaskedAssign(notOnSurface , (-pDotV3d + Sqrt(d2)), &distance);
    done|=notOnSurface;
    if( IsFull(done) )return;
    
}



*/


} // End global namespace

#endif // VECGEOM_VOLUMES_KERNEL_SPHEREIMPLEMENTATION_H_
