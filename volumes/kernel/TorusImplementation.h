/// @file TorusImplementation.h

#ifndef VECGEOM_VOLUMES_KERNEL_TORUSIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_TORUSIMPLEMENTATION_H_


#include "base/Global.h"
#include "base/Transformation3D.h"
#include "volumes/kernel/GenericKernels.h"
#include "volumes/kernel/TubeImplementation.h"
#include "volumes/UnplacedTorus.h"
#include <math.h>


#ifndef VECGEOM_NVCC
#include <iomanip>
#include <Vc/Vc>
#endif

namespace VECGEOM_NAMESPACE {


#ifndef VECGEOM_NVCC
using namespace Vc;

inline
double_v Vccbrt( double_v x )
{
  Vc::double_v xorig=x;
  x=Vc::abs(x);
  Vc::double_v tmp= Vc::exp(0.33333333333333333333*Vc::log(x) );
  return tmp.copySign(xorig);
}
#endif


template<class T>
class Complex
{
 private:
  T fR, fI;

 public:

 VECGEOM_CUDA_HEADER_BOTH
 Complex() : fR(T(0.)), fI(T(0.)) {}

 VECGEOM_CUDA_HEADER_BOTH
 Complex(T x, T y) : fR(x), fI(y) {}
 VECGEOM_CUDA_HEADER_BOTH
  ~Complex(){}

  // copy constructor
 VECGEOM_CUDA_HEADER_BOTH
  Complex( const Complex & other )
    {
      fR=other.fR;
      fI=other.fI;
    }

 VECGEOM_CUDA_HEADER_BOTH
 Complex& operator=( const Complex<T>& x)
   {
     fR=x.fR;
     fI=x.fI;
     return *this;
   }

 template <typename S>
 VECGEOM_CUDA_HEADER_BOTH
 Complex& operator=( const S& x)
   {
     fR=x;
     fI=0;
     return *this;
   }

 VECGEOM_CUDA_HEADER_BOTH
 T real() const {return fR;}
 VECGEOM_CUDA_HEADER_BOTH
 T imag() const {return fI;}

 // problem: how can we vary the math functions here
 // for Vc:: argument dependent lookup might help
 VECGEOM_CUDA_HEADER_BOTH
 T carg() const {return ATan2(fI,fR);}
 VECGEOM_CUDA_HEADER_BOTH
 T cabs() const {return Sqrt(fR*fR+fI*fI);}

 VECGEOM_CUDA_HEADER_BOTH
 Complex conj() const { return Complex(fR,-fI);}

 VECGEOM_CUDA_HEADER_BOTH
 Complex operator+=( Complex const & x ) const
 {
   return Complex( fR + x.real(), fI + x.imag() );
 }

 VECGEOM_CUDA_HEADER_BOTH
 Complex operator-() const
 {
   return Complex(-fR,-fI);
 }

}; // end class complex

template <typename T>
inline
 VECGEOM_CUDA_HEADER_BOTH
Complex<T> operator+( Complex<T> const & x, Complex<T> const & y )
{
  return Complex<T>( x.real() + y.real(), x.imag() + y.imag() );
}

template <typename T>
inline
 VECGEOM_CUDA_HEADER_BOTH
Complex<T> operator+( Complex<T> const & x, T const & y )
{
  return Complex<T>( x.real()+y , x.imag()  );
}

template <typename T>
inline
 VECGEOM_CUDA_HEADER_BOTH
Complex<T> operator+( T const & y, Complex<T> const & x )
{
  return Complex<T>( x.real()+y , x.imag()  );
}


template <typename T>
inline
 VECGEOM_CUDA_HEADER_BOTH
Complex<T> operator-( Complex<T> const & x, Complex<T> const & y )
{
  return Complex<T>( x.real() - y.real(), x.imag() - y.imag() );
}

template <typename T>
inline
 VECGEOM_CUDA_HEADER_BOTH
Complex<T> operator-( Complex<T> const & x, T const & y )
{
  return Complex<T>( x.real() - y, x.imag() );
}

template <typename T>
inline
 VECGEOM_CUDA_HEADER_BOTH
Complex<T> operator-(  T const & y,  Complex<T> const & x)
{
  return Complex<T>( y - x.real() , -x.imag() );
}


template <typename T>
inline
 VECGEOM_CUDA_HEADER_BOTH
Complex<T> operator*( Complex<T> const & x, Complex<T> const & y )
{
  return Complex<T>( x.real()*y.real() - x.imag()*y.imag(), x.real()*y.imag() + x.imag()*y.real() );
}

template <typename T, typename Other>
inline
 VECGEOM_CUDA_HEADER_BOTH
Complex<T> operator*( Complex<T> const & x, Other const & y )
{
  return Complex<T>( x.real()*y , x.imag()*y );
}

template <typename T, typename Other>
inline
 VECGEOM_CUDA_HEADER_BOTH
Complex<T> operator*( Other const & y, Complex<T> const & x )
{
  return Complex<T>( x.real()*y , x.imag()*y );
}


// division by T
template <typename T>
inline
 VECGEOM_CUDA_HEADER_BOTH
Complex<T> operator/( Complex<T> const & x, T const & y )
{
  // could use some fast math here
  // division by zero is not treated
  T invy = 1./y;
  return Complex<T>( x.real()*invy, x.imag()*invy );
}


template <typename T>
inline
 VECGEOM_CUDA_HEADER_BOTH
Complex<T> operator/( T const & x, Complex<T> const & y )
{
  // multiply by conjugate
  Complex<T> tmp = y*y.conj();
  return Complex<T>( x*y.conj()) / tmp.real();
}

template <typename T>
inline
 VECGEOM_CUDA_HEADER_BOTH
Complex<T> operator/( Complex<T> const & x, Complex<T> const & y )
{
  // multiply by conjugate
  Complex<T> tmp = y*y.conj();
  return Complex<T>( x*y.conj()) / tmp.real();
}


// standalone function for Sqrt
template <typename T>
inline
 VECGEOM_CUDA_HEADER_BOTH
Complex<T> csqrt( const Complex<T>& x )
{
  T r = x.real();
  T i = x.imag();
  T l = Sqrt( r*r + i*i );
  return Complex<T>( Sqrt( 0.5*(l+r) ), Sqrt( 0.5*(l-r) ));
}


// standalone function for sqrt
template <typename T>
inline
 VECGEOM_CUDA_HEADER_BOTH
Complex<T> csqrtrealargument( const T & x )
{
  T imcoef   = (x>=0.) ? 0. : 1.;
  T realcoef = (x>=0.) ? 1. : 0.;
  T l = Sqrt( fabs(x) );
  return Complex<T>( realcoef * l , imcoef * l );
}

#ifndef VECGEOM_NVCC
// template specialization for Vc
typedef Vc::double_v VCT;
template <>
inline
Complex<VCT> csqrtrealargument( const VCT & x )
{
  VCT::Mask real = (x>=0.);
  VCT l = Vc::sqrt( Vc::abs(x) );
  VCT realpart(Vc::Zero);
  VCT impart(Vc::Zero);
  realpart(real) = l;
  impart(!real) = l;
  return Complex<VCT>( realpart , impart );
}
#endif

// we need standalone function for cubic root
template <typename T>
VECGEOM_CUDA_HEADER_BOTH
inline
Complex<T> cbrt( const Complex<T>& x )
{
  // use the sin/cosine formula??
  T r;

  T sinnewangle, cosnewangle;
  if( x.imag() != 0.0 )
    {
      r = x.cabs();
      T angle = x.carg();

      T newangle = angle/3.;
      sinnewangle=sin(newangle);
      cosnewangle=cos(newangle);
      //sincos(newangle, &sinnewangle, &cosnewangle);
    }
  else
    {
      r = x.real();
      sinnewangle=0.;
      cosnewangle=1.;
    }
  // use ordinary cubic root function ( is it in Vc ?? )
  T rcbrt = Pow(r,1./3);//cbrt( r );
  return Complex<T>(  rcbrt*cosnewangle, rcbrt*sinnewangle );
}

#ifndef VECGEOM_NVCC
// template specialization for Vc
// we need standalone function for cubic root
template <>
inline
Complex<Vc::double_v> cbrt( const Complex<Vc::double_v>& x )
{
  // use the sin/cosine formula??
  Vc::double_v r;
  Vc::double_v sinnewangle, cosnewangle;
  Vc::double_m imaginary = x.imag() != 0.0;
  if( ! imaginary.isEmpty() )
    {
      r = x.cabs();
      Vc::double_v angle = x.carg();
      angle( ! imaginary ) = 0.;
      r( ! imaginary ) = r.copySign( x.real() );
      Vc::double_v newangle = 0.33333333333333333*angle;
      sincos(newangle, &sinnewangle, &cosnewangle);
    }
  else
    {
      r = x.real();
      sinnewangle=Vc::Zero;
      cosnewangle=Vc::One;
    }
  // use cubic root function defined above
  Vc::double_v rcbrt = Vccbrt( r );
  return Complex<Vc::double_v>(  rcbrt*cosnewangle, rcbrt*sinnewangle );
}
#endif


#ifndef VECGEOM_CUDA
template <typename T, typename stream>
stream& operator<<(stream& s, const Complex<T> & x)
{
  s << "(" << x.real() << " , " << x.imag() << ")";
  return s;
}
#endif

static const double inv256 = 1./256;
static const double inv16 = 1./16;
static const double inv3 = 1./3;
static const double inv128 = 1./128;
static const double inv27 = 1./27;
static const double inv6 = 1./6;
static const double inv108 = 1./108;
static const double inv12 = 1./12;


template <typename CT>
VECGEOM_CUDA_HEADER_BOTH
void solveQuartic(double a, double b, double c, double d, double e, CT * roots)
{
  // Uses Ferrari's Method; this will vectorize trivially ( at least when we treat the real and imaginary parts separately )
  double aa = a*a, aaa=aa*a, bb=b*b, bbb=bb*b;
  double alpha = -3.0*bb/(8.0*aa)   + c/a, alpha2 = alpha * alpha;
  double beta  =    bbb/(8.0*aaa) + b*c/(-2.0*aa) + d/a;
  double gamma = -3.0*bbb*b/(256.0*aaa*a) + c*bb/(16.0*aaa) + b*d/(-4.0*aa) + e/a;


  CT P = -alpha2/12.0 - gamma;
  CT Q = -alpha2*alpha/108.0 + alpha*gamma/3.0 - beta*beta/8.0;
  CT R = Q*0.5 + Sqrt(Q*Q*0.25 + P*P*P/27.0);
  CT U = Pow(R,1./3);
  CT y = -5.0*alpha/6.0 - U;
  if(U != 0.0) y += P/(3.0*U);
  CT W = Sqrt(alpha + y + y);

  CT aRoot;
   double firstPart = b/(-4.0*a);
  CT secondPart = -3.0*alpha - 2.0*y;
  CT thirdPart = 2.0*beta/(W);

  aRoot = firstPart + 0.5 * (-W - Sqrt(secondPart + thirdPart));
  roots[0] = aRoot;

  aRoot = firstPart + 0.5 * (-W + Sqrt(secondPart + thirdPart));
  roots[1] = aRoot;

  aRoot = firstPart + 0.5 * (W - Sqrt(secondPart - thirdPart));
  roots[2] = aRoot;

  aRoot = firstPart + 0.5 * (W + Sqrt(secondPart - thirdPart));
  roots[3] = aRoot;

}

inline
double foo(double a, double b)
{
    return -3.0*inv256*a + b*a;
}
// finding all real roots
// this is Ferrari's method which is not really numerically stable but very elegant
// CT == complextype
//template <typename CT>
typedef Complex<double> CT;
VECGEOM_CUDA_HEADER_BOTH
inline
void solveQuartic2(double a, double b, double c, double d, double e, CT * roots)
{

  // Uses Ferrari's Method; this will vectorize trivially ( at least when we treat the real and imaginary parts separately )
  //  double inva=1./a;
  //  double invaa = inva*inva;
  //  double invaaa=invaa*inva;
  double bb=b*b;
  double bbb=bb*b;
  double alpha = -3.0*0.125*bb   + c;
  double beta  =    0.125*bbb - 0.5*b*c + d;
  double alpha2 = alpha * alpha;
  double gamma = -3.0*bbb*b*inv256 + c*bb*inv16 - 0.25*b*d + e;

  /* std::cerr << alpha << "\n"; */
  /* std::cerr << alpha2 << "\n"; */
  /* std::cerr << beta << "\n"; */
  /* std::cerr << gamma << "\n"; */

  double P = -alpha2*inv12 - gamma;
  double Q = -alpha2*alpha*inv108 + alpha*gamma*inv3 - 0.125*beta*beta;
  //   std::cerr << "P " << P << "\n"; 
  //   std::cerr << "Q " << Q << "\n"; 

  double tmp = 0.25*Q*Q + P*P*P*inv27;
  CT R = Q*0.5 + csqrtrealargument(tmp);
  //    std::cerr << "R " << R << "\n";
  CT U = cbrt(R);
  //    std::cerr << "U " << U << "\n";
  //    std::cerr << "U*U*U " << U*U*U << "\n";

  CT y = -5.0*alpha*inv6 - U;
  //    std::cerr << "y " << y << "\n";
  y = y + P/(3.*U);
  //    std::cerr << "y " << y << "\n";
  CT W = csqrt((alpha + y) + y);
  //    std::cerr << "W " << W << "\n";

  CT aRoot;
  double firstPart = -0.25*b;
  CT secondPart = -3.0*alpha - 2.0*y;
  CT thirdPart = (2.0*beta)/W;

  aRoot = firstPart + 0.5 * (-W - csqrt(secondPart + thirdPart));
  roots[0] = aRoot;
  aRoot = firstPart + 0.5 * (-W + csqrt(secondPart + thirdPart));
  roots[1] = aRoot;
  aRoot = firstPart + 0.5 * (W - csqrt(secondPart - thirdPart));
  roots[2] = aRoot;
  aRoot = firstPart + 0.5 * (W + csqrt(secondPart - thirdPart));
  roots[3] = aRoot;
}

// finding all real roots
// this is Ferrari's method which is not really numerically stable but very elegant
// CT == complextype
//typedef Vc::double_v VCT2;
#ifndef VECGEOM_NVCC
typedef Complex<VCT> CVCT;
inline
void solveQuartic2(VCT a, VCT b, VCT c, VCT d, VCT e, CVCT * roots)
{
  //  VCT inva=1./a;
  //  VCT invaa = inva*inva;
  //  VCT invaaa=invaa*inva;
  VCT bb=b*b;
  VCT bbb=bb*b;
  VCT alpha = -3.0*0.125*bb + c, alpha2 = alpha * alpha;
  VCT beta  = 0.125*bbb - 0.5*b*c+ d;
  VCT gamma = -3.0*bbb*b*inv256 + c*bb*inv16 - 0.25*b*d + e;

  VCT P = -alpha2*inv12 - gamma;
  VCT Q = -alpha2*alpha*inv108 + alpha*gamma*inv3 - 0.125*beta*beta;
  //   std::cerr << "P " << P << "\n"; 
  //   std::cerr << "Q " << Q << "\n"; 


  VCT tmp = 0.25*Q*Q + P*P*P*inv27;
  CVCT R = Q*0.5 + csqrtrealargument(tmp);
  CVCT U = cbrt(R);
  //    std::cerr << "R " << R << "\n";
  //    std::cerr << "U " << U << "\n";
  //    std::cerr << "U*U*U " << U*U*U << "\n";
  CVCT y = -5.0*alpha*inv6 - U;
  y = y + P/(3.*U );
  //    std::cerr << "y " << y << "\n";
  CVCT W = csqrt((alpha + y) + y);
  //    std::cerr << "W " << W << "\n";

  CVCT aRoot;

  VCT firstPart = -0.25*b;
  CVCT secondPart = -3.0*alpha - 2.0*y;
  CVCT thirdPart = (2.0*beta)/(W);

  aRoot = firstPart + 0.5 * (-W - csqrt(secondPart + thirdPart));
  roots[0] = aRoot;

  aRoot = firstPart + 0.5 * (-W + csqrt(secondPart + thirdPart));
  roots[1] = aRoot;

  aRoot = firstPart + 0.5 * (W - csqrt(secondPart - thirdPart));
  roots[2] = aRoot;

  aRoot = firstPart + 0.5 * (W + csqrt(secondPart - thirdPart));
  roots[3] = aRoot;
}
#endif

template <TranslationCode transCodeT, RotationCode rotCodeT>
struct TorusImplementation {
  /*
  // First Implementation of Contains/Inside without GenericKernel
  // Tested with Benchmarker
  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void UnplacedContains( UnplacedTorus const &torus,
      Vector3D<typename Backend::precision_v> const &point,
      typename Backend::bool_v &inside) {

   typedef typename Backend::precision_v Float_t;

   // TODO: do this generically WITH a generic contains/inside kernel
   // forget about sector for the moment

   Float_t rxy = Sqrt(point[0]*point[0] + point[1]*point[1]);
   Float_t radsq = ( rxy - torus.rtor() ) * (rxy - torus.rtor() ) + point[2]*point[2];
   inside = radsq > torus.rmin2() && radsq < torus.rmax2();
   
}
*/
  /*
  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void Contains(
      UnplacedTorus const &torus,
      Transformation3D const &transformation,
      Vector3D<typename Backend::precision_v> const &point,
      Vector3D<typename Backend::precision_v> &localPoint,
      typename Backend::bool_v &inside) {

    localPoint = transformation.Transform<transCodeT, rotCodeT>(point);
    UnplacedContains<Backend>(torus, localPoint, inside);
  }
  */
  /*
  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void Inside(UnplacedTorus const &torus,
                     Transformation3D const &transformation,
                     Vector3D<typename Backend::precision_v> const &point,
                     typename Backend::inside_v &inside) {
    
 
   // TODO: No Phi sector for the moment
   typedef typename Backend::precision_v Float_t;
   typedef typename Backend::bool_v      Bool_t;
 
   Float_t rxy = Sqrt(point[0]*point[0] + point[1]*point[1]);
   Float_t radsq = ( rxy - torus.rtor() ) * (rxy - torus.rtor() ) + point[2]*point[2];
   
   // Geant4 code with Surface included
 
   Float_t r2, tolRMin = 0.0, tolRMax, fRminTolerance = kTolerance, fRmaxTolerance = kTolerance ;
   
   MaskedAssign( torus.rmin(),torus.rmin() + fRminTolerance , &tolRMin );
   tolRMax = torus.rmax() - fRmaxTolerance;

   Bool_t ifOutside = (radsq < tolRMin*tolRMin) || (radsq > tolRMax*tolRMax);
   MaskedAssign ( ifOutside, EInside::kOutside, &inside); //exit early return?
   //std::cout<<"ifOutside="<<ifOutside<<" rad2="<<radsq<<" p="<<point<<std::endl;
   if (Backend::early_returns && IsFull(ifOutside) )
   {
        inside = EInside::kOutside;
        return;
   }
 
   tolRMin = torus.rmin() - fRminTolerance ;
   Bool_t ifNegativeTolerance = (tolRMin < 0);// =RMin=0
  
   tolRMax = torus.rmax() + fRmaxTolerance ;
    
   Bool_t ifSurface = ( (!ifOutside) && (!ifNegativeTolerance) && (radsq>= tolRMin*tolRMin) ) && (radsq >= tolRMax*tolRMax) ;
   //  std::cout<<"ifSurface="<<ifSurface<<" rad2="<<radsq<<" p="<<point<<std::endl;
   MaskedAssign ( ifSurface, EInside::kSurface, &inside); //exit early return?
   if (Backend::early_returns && IsFull(ifSurface) )
   {
       inside = EInside::kSurface;
       return;
   }
   Bool_t ifInside = (!ifSurface) && (!ifOutside);
   // std::cout<<"conditionINside ="<< ifInside<<" rad2="<<radsq<<" p="<<point<<std::endl;
   MaskedAssign( ifInside, EInside::kInside, &inside);
   
  
   
  }

  */

  /////GenericKernel Contains/Inside implementation
  template <typename Backend, bool ForInside>
  VECGEOM_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  static void GenericKernelForContainsAndInside(UnplacedTorus const &torus,
          Vector3D<typename Backend::precision_v> const &point,
          typename Backend::bool_v &completelyinside,
          typename Backend::bool_v &completelyoutside)

   {
   // using vecgeom::GenericKernels;
   // here we are explicitely unrolling the loop since  a for statement will likely be a penality
   // check if second call to Abs is compiled away
   // and it can anyway not be vectorized
    /* rmax */
    typedef typename Backend::precision_v Float_t;
    Float_t rxy = Sqrt(point[0]*point[0] + point[1]*point[1]);
    Float_t radsq = ( rxy - torus.rtor() ) * (rxy -  torus.rtor() ) + point[2]*point[2];

    // completelyoutside = radsq > MakePlusTolerant<ForInside>( torus.rmax2() );//rmax
    completelyoutside = radsq > MakePlusTolerantSquare<ForInside>( torus.rmax(),torus.rmax2() );//rmax
    //std::cout<<"Kernelrmax  point="<<point<<" radsq="<<radsq<<" Out?"<<completelyoutside<<std::endl;
    if (ForInside)
    {
      // completelyinside = radsq < MakeMinusTolerant<ForInside>( torus.rmax2() );
      completelyinside = radsq < MakeMinusTolerantSquare<ForInside>( torus.rmax(),torus.rmax2() );
    }
    if (Backend::early_returns) {
      if ( IsFull(completelyoutside) ) {
        return;
      }
    }
    /* rmin */
    //completelyoutside |= radsq < MakePlusTolerant<ForInside>( torus.rmin2() );//rmin
   completelyoutside |= radsq < MakePlusTolerantSquare<ForInside>( torus.rmin(),torus.rmin2() );//rmin
   // std::cout<<"Kernelrmin  point="<<point<<" radsq="<<radsq<<" Out?"<<completelyoutside<<std::endl;
    if (ForInside)
    {
      //completelyinside &= radsq > MakeMinusTolerant<ForInside>( torus.rmin2() );
      completelyinside &= radsq > MakeMinusTolerantSquare<ForInside>( torus.rmin(),torus.rmin2() );
    }

    // NOT YET NEEDED WHEN NOT PHI TREATMENT
    //        if (Backend::early_returns) {
    //      if ( IsFull(completelyoutside) ) {
    //        return;
    //      }
    //    }
    /* phi TOO DO*/

  }

  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void ContainsKernel(
      UnplacedTorus const &torus,
      Vector3D<typename Backend::precision_v> const &point,
      typename Backend::bool_v &inside)
  {
   typedef typename Backend::bool_v Bool_t;
    Bool_t unused;
    Bool_t outside;
    GenericKernelForContainsAndInside<Backend, false>(torus,
    point, unused, outside);
    inside = !outside;
 }
  //template <TranslationCode transCodeT, RotationCode rotCodeT>
  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  static void InsideKernel(
     UnplacedTorus const &torus,
     Vector3D<typename Backend::precision_v> const &point,
     typename Backend::inside_v &inside) {

  typedef typename Backend::bool_v      Bool_t;
  inside = EInside::kOutside;
  //Check Bounding tube first
   Bool_t inBounds;
      TubeImplementation<translation::kIdentity,
      rotation::kIdentity, TubeTypes::HollowTube>::
      UnplacedContains<Backend>(
      torus.GetBoundingTube(), point, inBounds);
      if( !inBounds ){
        if (Backend::early_returns) {
                    return;  
        }
      }
  //
  Bool_t completelyinside, completelyoutside;
  GenericKernelForContainsAndInside<Backend,true>(torus, 
  point, completelyinside, completelyoutside);
  inside = EInside::kSurface;
  MaskedAssign(completelyoutside, EInside::kOutside, &inside);
  MaskedAssign(completelyinside, EInside::kInside, &inside);
}
  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void UnplacedContains( UnplacedTorus const &torus,
      Vector3D<typename Backend::precision_v> const &point,
      typename Backend::bool_v &inside) {

 
    // TODO: do this generically WITH a generic contains/inside kernel
    // forget about sector for the moment
  
     ContainsKernel<Backend>(torus, point, inside);
}
  template <typename Backend>
  VECGEOM_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  static void Contains(
      UnplacedTorus const &unplaced,
      Transformation3D const &transformation,
      Vector3D<typename Backend::precision_v> const &point,
      Vector3D<typename Backend::precision_v> &localPoint,
      typename Backend::bool_v &inside){
   localPoint = transformation.Transform<transCodeT, rotCodeT>(point);
   UnplacedContains<Backend>(unplaced, localPoint, inside);

  }
  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void Inside(UnplacedTorus const &torus,
                     Transformation3D const &transformation,
                     Vector3D<typename Backend::precision_v> const &point,
                     typename Backend::inside_v &inside) {
  
    InsideKernel<Backend>(torus, transformation.Transform<transCodeT, rotCodeT>(point),  inside);

}
  

  /////End GenericKernel Contains/Inside implementation


    template <class T>
    VECGEOM_CUDA_HEADER_BOTH
    static
    //__attribute__((noinline))
    T CheckZero(T b, T c, T d, T e, T x)
    {
        T x2=x*x;
        return x2*x2 + b*x2*x + c*x2 + d*x + e;
    }

    template <class T>
    VECGEOM_CUDA_HEADER_BOTH
    static
    T NewtonIter(T b, T c, T d, T e, T x, T fold)
    {
        T x2 = x*x;
        T fprime = 4*x2*x + 3*b*x2 + 2*c*x + d;
        return x - fold/fprime;
    }


   template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static
  typename Backend::precision_v ToBoundary(
          UnplacedTorus const& torus,
          Vector3D<typename Backend::precision_v> const &point,
          Vector3D<typename Backend::precision_v> const &dir,
          Precision radius )
  {
      // Returns distance to the surface or the torus (fR,r) from a point, along
      // a direction. Point is close enough to the boundary so that the distance
      // to the torus is decreasing while moving along the given direction.

     // Compute coefficients of the quartic polynomial
     typedef typename Backend::precision_v Float_t;
     typedef typename Backend::bool_v Bool_t;

     Float_t tol = kTolerance;

     // actually a scalar product
     Float_t r0sq  = point[0]*point[0]+point[1]*point[1]+point[2]*point[2];

     // actually a scalar product
     Float_t rdotn = point[0]*dir[0]+point[1]*dir[1]+point[2]*dir[2];

     // can be precomputed
     Float_t rsumsq = torus.rtor2() + radius*radius;
     Float_t a = 4.*rdotn;
     Float_t b = 2.*(r0sq+2.*rdotn*rdotn-rsumsq+2.*torus.rtor2()*dir[2]*dir[2]);
     Float_t c = 4.*(r0sq*rdotn-rsumsq*rdotn+2.*torus.rtor2()*point[2]*dir[2]);
     Float_t d = r0sq*r0sq-2.*r0sq*rsumsq+4.*torus.rtor2()*point[2]*point[2]+(torus.rtor2()-radius*radius)*(torus.rtor2()-radius*radius);

     //std::cerr << "#a " << a << "\n"; 
     // std::cerr << "#b " << b << "\n"; 
     //  std::cerr << "#c " << c << "\n"; 
     //  std::cerr << "#d " << d << "\n"; 

     //   std::cerr << "#torus par " << torus.rtor2()<<" "<<radius<<" "<<r0sq  << "\n"; 
     

     // the 4 complex roots
     Complex<Float_t> roots[4];
     // get the roots
     solveQuartic2(Float_t(1.),a,b,c,d,roots);

     //    std::cerr << "#ROOTS " << roots[0] << "\n"; 
     //    std::cerr << "#ROOTS " << roots[1] << "\n"; 
     //    std::cerr << "#ROOTS " << roots[2] << "\n"; 
     //    std::cerr << "#ROOTS " << roots[3] << "\n"; 

     Float_t validdistance = kInfinity;
     Bool_t havevalidsolution = Abs(roots[0].imag()) < 1E-10 && roots[0].real() > 0.;
     MaskedAssign( havevalidsolution, roots[0].real(), &validdistance );

     havevalidsolution = Abs(roots[1].imag()) < 1E-10 && roots[1].real() > 0.;
     MaskedAssign( havevalidsolution, Min(roots[1].real(), validdistance), &validdistance );

     havevalidsolution = Abs(roots[2].imag()) < 1E-10 && roots[2].real() > 0.;
     MaskedAssign( havevalidsolution, Min(roots[2].real(), validdistance), &validdistance );

     havevalidsolution = Abs(roots[3].imag()) < 1E-10 && roots[3].real() > 0.;
     MaskedAssign( havevalidsolution, Min(roots[3].real(), validdistance), &validdistance );

     // TODO: only do this in case there is any finite real solution
     havevalidsolution = (validdistance < kInfinity);
     if( havevalidsolution){
      validdistance = NewtonIter(a,b,c,d,validdistance,CheckZero(a,b,c,d,validdistance));
      validdistance = NewtonIter(a,b,c,d,validdistance,CheckZero(a,b,c,d,validdistance));
     
     }
     //std::cerr << std::setprecision(20);
     //std::cerr << "#DISTANCE " << validdistance;
     //Float_t fold = CheckZero(a, b, c, d, validdistance);
     //std::cerr << point << "\n";
     //std::cerr << " ZERO CHECK " << fold << " dist " << validdistance << "\n";
     //Float_t newdist = NewtonIter(a, b, c, d, validdistance, fold);
     //std::cerr << " NEWDIST " << newdist;
     //std::cerr << " " << CheckZero( a,b,c,d, newdist);
     //std::cerr << "\n";
     //std::cerr << std::setprecision(5);
     
     return validdistance;
  }
  
    template <class Backend>
    VECGEOM_CUDA_HEADER_BOTH
    VECGEOM_INLINE
    static void DistanceToIn(
        UnplacedTorus const &torus,
        Transformation3D const &transformation,
        Vector3D<typename Backend::precision_v> const &point,
        Vector3D<typename Backend::precision_v> const &direction,
        typename Backend::precision_v const &stepMax,
        typename Backend::precision_v &distance) {
      
      // TODO
      //No Phi Section for the moment
      typedef typename Backend::precision_v Float_t;
      typedef typename Backend::bool_v Bool_t;
      
      Vector3D<Float_t> localPoint = transformation.Transform<transCodeT, rotCodeT>(point);
      Vector3D<Float_t> localDirection = transformation.Transform<transCodeT, rotCodeT>(direction);

      ////////First naive implementation
      distance = kInfinity;

      //Check Bounding Cylinder first
      Bool_t inBounds;
      Bool_t missTorus;
      Bool_t done;
      Float_t tubeDistance = kInfinity;

      // call the tube functionality -- first of all we check whether we are inside
      // bounding volume
      TubeImplementation<translation::kIdentity,
        rotation::kIdentity, TubeTypes::HollowTube>::UnplacedContains<Backend>(
            torus.GetBoundingTube(), localPoint, inBounds);


      // only need to do this check if all particles (in vector) are outside ( otherwise useless )
      TubeImplementation<translation::kIdentity,
            rotation::kIdentity, TubeTypes::HollowTube>::DistanceToIn<Backend>(
                torus.GetBoundingTube(), transformation, localPoint,
                localDirection, stepMax, tubeDistance);

       //std::cout<<"tubeDistance="<<tubeDistance<<std::endl;
       done = (!inBounds && tubeDistance == kInfinity);

       if (Backend::early_returns) {
           if ( IsFull(done) ) {
            return;
       }
       }

    // move close the points which where outside bounding tube
    Vector3D<Float_t> forwardedPoint = localPoint + localDirection*tubeDistance;

    MaskedAssign( !inBounds, forwardedPoint[0], &localPoint[0]);
    MaskedAssign( !inBounds, forwardedPoint[1], &localPoint[1]);
    MaskedAssign( !inBounds, forwardedPoint[2], &localPoint[2]);

    // move points closer to torus ( for better numerical stability )
    //localPoint = localPoint+localDirection*tubeDistance;

    Float_t dout = ToBoundary<Backend>(torus,localPoint,localDirection,torus.rmax());

    Float_t din(kInfinity);
    bool hasrmin = (torus.rmin() > 0);
    if( hasrmin ){
      din = ToBoundary<Backend>(torus,localPoint,localDirection,torus.rmin());
    }

    MaskedAssign( !inBounds && !done, Min(dout+tubeDistance, din+tubeDistance), &distance);
    MaskedAssign(  inBounds && !done, Min(dout, din), &distance);

    bool hasphi = (torus.dphi() < kTwoPi );
    if( hasphi )
    {
     // TODO
    }
    /*bool checkSolution=1;
    if ( checkSolution && (distance < kInfinity ) )
      {
           Vector3D<typename Backend::precision_v> ptemp=point+direction*distance;
            Float_t rxy = Sqrt(ptemp[0]*ptemp[0] + ptemp[1]*ptemp[1]);
            Float_t radsq = ( rxy - torus.rtor() ) * (rxy - torus.rtor() ) + ptemp[2]*ptemp[2];
            Float_t dif=Abs(radsq-torus.rmax2());
	    if( dif > 1E-6)std::cerr<<"radsq="<<radsq<<" rmax2="<<torus.rmax2()<<" dif="<<radsq-torus.rmax2()<<std::endl;

            //Float_t douut=ToBoundary<Backend>(torus,ptemp,direction,torus.rmax());;
            // std::cerr << point << "\n";   
            //std::cerr<<"VC ptempZ="<<ptemp[2]<<" Dout="<<douut<<std::endl;
        
           }
    */
  }

  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void DistanceToOut(
      UnplacedTorus const &torus,
      Vector3D<typename Backend::precision_v> const &point,
      Vector3D<typename Backend::precision_v> const &dir,
      typename Backend::precision_v const &stepMax,
      typename Backend::precision_v &distance) {

    typedef typename Backend::precision_v Float_t;
    distance = kInfinity;

    // TODO
    // Compute distance from inside point to surface of the torus.
    bool hasphi = (torus.dphi()<kTwoPi);
    bool hasrmin = (torus.rmin()>0);

    Float_t dout = ToBoundary<Backend>(torus,point,dir,torus.rmax());
    Float_t din(kInfinity);
    if( hasrmin )
    {
       din = ToBoundary<Backend>(torus,point,dir,torus.rmin());
    }
    distance = Min(dout,din);

    if( hasphi )
    {
    // TODO
    }
  }

  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void SafetyToIn(UnplacedTorus const &torus,
                         Transformation3D const &transformation,
                         Vector3D<typename Backend::precision_v> const &point,
                         typename Backend::precision_v &safety) {

    typedef typename Backend::precision_v Float_t;

    Vector3D<Float_t> localPoint = transformation.Transform<transCodeT, rotCodeT>(point);
   
    // implementation taken from TGeoTorus
    
    Float_t rxy = Sqrt( localPoint[0]*localPoint[0] + localPoint[1]*localPoint[1] );
    Float_t rad = Sqrt( (rxy - torus.rtor())*(rxy-torus.rtor()) + localPoint[2]*localPoint[2] );
    safety = rad - torus.rmax();
    if( torus.rmin() )
      {safety = Max( torus.rmin()- rad, rad - torus.rmax() );}
  
    
    //std::cerr << "#SAF IN "<<" rxy="<<rxy <<" rad="<< rad<<" saf="<<safety<<std::endl;
    // TODO: extend implementation for phi sector case

  }

  template <class Backend>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static void SafetyToOut(UnplacedTorus const &torus,
                          Vector3D<typename Backend::precision_v> const &point,
                          typename Backend::precision_v &safety) {

    typedef typename Backend::precision_v Float_t;
    // implementation taken from TGeoTorus
    
    Float_t rxy = Sqrt( point[0]*point[0] + point[1]*point[1] );
    Float_t rad = Sqrt( (rxy - torus.rtor())*(rxy-torus.rtor()) + point[2]*point[2] );
    safety= torus.rmax() - rad;
    if( torus.rmin() )
       { safety = Min( rad - torus.rmin(), torus.rmax() - rad );}
    
    // TODO: extend implementation for phi sector case
    
}

}; // end struct

} // end namespace


#endif // VECGEOM_VOLUMES_KERNEL_TORUSIMPLEMENTATION_H_
