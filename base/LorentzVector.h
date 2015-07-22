/// \file vectorlorentz.h
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch) + fca

#ifndef VECGEOM_BASE_LORENTZVECTOR_H_
#define VECGEOM_BASE_LORENTZVECTOR_H_

#include "base/Global.h"
#include "base/Vector3D.h"

#include "backend/Backend.h"
#ifndef VECGEOM_NVCC
  #if (defined(VECGEOM_VC) || defined(VECGEOM_VC_ACCELERATION))
    #include <Vc/Vc>
  #endif
#endif
#include "base/AlignedBase.h"
#include "base/Vector3D.h"
#include "base/RNG.h"

#include <cstdlib>
#include <ostream>
#include <string>
#include <iostream>

using std::vector;

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE( template <typename T> class LorentzVector; )

inline namespace VECGEOM_IMPL_NAMESPACE {


/**
 * @brief Lorentz dimensional vector class supporting most arithmetic operations.
 * @details If vector acceleration is enabled, the scalar template instantiation
 *          will use vector instructions for operations when possible.
 */
template <typename T>
class LorentzVector : public AlignedBase {

  typedef LorentzVector<T> VecType;

private:

  T fVec[4];

public:

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  LorentzVector(const T a, const T b, const T c, const T d): 
  fVec{a,b,c,d} {}

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  LorentzVector() : fVec{0,0,0,0} {}

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  LorentzVector(const T a) : fVec{0,0,0,0} {}

  template <typename U>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  LorentzVector(LorentzVector<U> const &other) : 
  fVec{other[0], other[1], other[2], other[3]} {}

  template <typename U>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  LorentzVector(Vector3D<U> const &other, T t) :
  fVec{other[0], other[1], other[2], t} {}

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  LorentzVector& operator=(LorentzVector const &other) {
    if(this != &other) {
       fVec[0] = other[0];
       fVec[1] = other[1];
       fVec[2] = other[2];
       fVec[3] = other[3];
    }
    return *this;
  }

  operator Vector3D<T>&() {return reinterpret_cast<Vector3D<T>&>(*this);}

  /**
   * Constructs a vector from an std::string of the same format as output by the
   * "<<"-operator for outstreams.
   * @param str String formatted as "(%d, %d, %d %d)".
   */
  VECGEOM_CUDA_HEADER_HOST
  LorentzVector(std::string const &str) {
    int begin = str.find("(")+1, end = str.find(",")-1;
    fVec[0] = std::atof(str.substr(begin, end-begin+1).c_str());
    begin = end + 2;
    end = str.find(",", begin)-1;
    fVec[1] = std::atof(str.substr(begin, end-begin+1).c_str());
    begin = end + 2;
    end = str.find(",", begin)-1;
    fVec[2] = std::atof(str.substr(begin, end-begin+1).c_str());
    begin = end + 2;
    end = str.find(")", begin)-1;
    fVec[3] = std::atof(str.substr(begin, end-begin+1).c_str());
  }

  /**
   * Contains no check for correct indexing to avoid impairing performance.
   * @param index Index of content in the range [0-2].
   */
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T& operator[](const int index) {
    return fVec[index];
  }

  /**
   * Contains no check for correct indexing to avoid impairing performance.
   * @param index Index of content in the range [0-2].
   */
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T const& operator[](const int index) const {
    return fVec[index];
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T& x() { return fVec[0]; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T const& x() const { return fVec[0]; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T& y() { return fVec[1]; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T const& y() const { return fVec[1]; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T& z() { return fVec[2]; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T const& z() const { return fVec[2]; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T& t() { return fVec[3]; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T const& t() const { return fVec[3]; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  void Set(T const &a, T const &b, T const &c, T const &d) {
    fVec[0] = a;
    fVec[1] = b;
    fVec[2] = c;
    fVec[3] = d;
  }

  template <typename U, typename V>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
     void Set(const Vector3D<U> &v, V t) {
    fVec[0] = v[0];
    fVec[1] = v[1];
    fVec[2] = v[2];
    fVec[3] = t;
  }

  template <typename U>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  void Set(const Vector3D<U> &v) {
    fVec[0] = v[0];
    fVec[1] = v[1];
    fVec[2] = v[2];
  }

  VECGEOM_CUDA_HEADER_BOTH
  void Set(const T a) {
    Set(a, a, a, a);
  }

  /// \return the length squared perpendicular to z direction
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T Perp2() const {
    return fVec[0]*fVec[0]+fVec[1]*fVec[1];
  }

  /// \return the length perpendicular to z direction
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T Perp() const {
    return Sqrt(Perp2());
  }

  template <typename U>
  ///The dot product of two LorentzVector<T> objects
  /// \return T (where T is float, double, or various SIMD vector types)
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static
  T Dot(LorentzVector<T> const &left, LorentzVector<U> const &right) {
     return - left[0]*right[0] - left[1]*right[1] - left[2]*right[2] + left[3]*right[3];
  }

  template <typename U>
  /// The dot product of two LorentzVector<T> objects
  /// \return T (where T is float, double, or various SIMD vector types)
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T Dot(LorentzVector<U> const &right) const {
    return Dot(*this, right);
  }

  // For UVector3 compatibility. Is equal to normal multiplication.
  // TODO: check if there are implicit dot products in USolids...
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  VecType MultiplyByComponents(VecType const &other) const {
    return *this * other;
  }

  /// \return Squared magnitude of the vector.
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T Mag2() const {
    return Dot(*this, *this);
  }

  /// \return Magnitude of the vector.
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T Mag() const {
    return Sqrt(Mag2());
  }

  /// \return Azimuthal angle between -pi and pi.
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T Phi() const {
    //T output = 0;
    //vecgeom::MaskedAssign(vec[0] != 0. || vec[1] != 0.,
    //                      ATan2(vec[1], vec[0]), &output);
    //return output;
    return ATan2(fVec[1], fVec[0]);
  }

  /// \return Polar angle between 0 and pi.
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T Theta() const {
    return ACos(fVec[2]/SpaceVector<T>().Mag());
  }

  /// Maps each vector entry to a function that manipulates the entry type.
  /// \param f A function of type "T f(const T&)" to map over entries.
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  void Map(T (*f)(const T&)) {
    fVec[0] = f(fVec[0]);
    fVec[1] = f(fVec[1]);
    fVec[2] = f(fVec[2]);
    fVec[3] = f(fVec[3]);
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  LorentzVector<T> Abs() const {
    return LorentzVector<T>(vecgeom::Abs(fVec[0]),
                          vecgeom::Abs(fVec[1]),
                          vecgeom::Abs(fVec[2]),
                          vecgeom::Abs(fVec[3]));
  }

  template <typename U>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  void MaskedAssign(LorentzVector<U> const &condition,
                    LorentzVector<T> const &value) {
    fVec[0] = (condition[0]) ? value[0] : fVec[0];
    fVec[1] = (condition[1]) ? value[1] : fVec[1];
    fVec[2] = (condition[2]) ? value[2] : fVec[2];
    fVec[3] = (condition[3]) ? value[3] : fVec[3];
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
     static VecType FromCylindrical(T r, T phi, T z, T t) {
     return VecType(r*cos(phi), r*sin(phi), z, t);
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  VecType& FixZeroes() {
    for (int i = 0; i < 4; ++i) {
      vecgeom::MaskedAssign(vecgeom::Abs(fVec[i]) < kTolerance, 0., &fVec[i]);
    }
    return *this;
  }

  // Inplace binary operators

  #define LORENTZVECTOR_TEMPLATE_INPLACE_BINARY_OP(OPERATOR) \
  VECGEOM_CUDA_HEADER_BOTH \
  VECGEOM_INLINE \
  VecType& operator OPERATOR(const VecType &other) { \
    fVec[0] OPERATOR other.fVec[0]; \
    fVec[1] OPERATOR other.fVec[1]; \
    fVec[2] OPERATOR other.fVec[2]; \
    fVec[3] OPERATOR other.fVec[3]; \
    return *this; \
  } \
  template <typename V> \
  VECGEOM_CUDA_HEADER_BOTH \
  VECGEOM_INLINE \
  VecType& operator OPERATOR(const LorentzVector<V> &other) { \
    fVec[0] OPERATOR other[0]; \
    fVec[1] OPERATOR other[1]; \
    fVec[2] OPERATOR other[2]; \
    fVec[3] OPERATOR other[3]; \
    return *this; \
  } \
  VECGEOM_CUDA_HEADER_BOTH \
  VECGEOM_INLINE \
  VecType& operator OPERATOR(const T &scalar) { \
    fVec[0] OPERATOR scalar; \
    fVec[1] OPERATOR scalar; \
    fVec[2] OPERATOR scalar; \
    fVec[3] OPERATOR scalar; \
    return *this; \
  }
  LORENTZVECTOR_TEMPLATE_INPLACE_BINARY_OP(+=)
  LORENTZVECTOR_TEMPLATE_INPLACE_BINARY_OP(-=)
  LORENTZVECTOR_TEMPLATE_INPLACE_BINARY_OP(*=)
  LORENTZVECTOR_TEMPLATE_INPLACE_BINARY_OP(/=)
  LORENTZVECTOR_TEMPLATE_INPLACE_BINARY_OP(*)

  #undef LORENTZVECTOR_TEMPLATE_INPLACE_BINARY_OP

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  operator bool() const {
    return fVec[0] && fVec[1] && fVec[2] && fVec[3];
  }

  template <typename U> 
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Vector3D<U> SpaceVector() const {
     return Vector3D<U>(fVec[0], fVec[1], fVec[2]);
  }

  template <typename U> 
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  VecType Boost(const Vector3D<U> &beta) const {
    double beta2 = beta.Mag2();
    if(beta2<=1e-16) return LorentzVector(*this);
    double gamma = Sqrt(1./(1-beta2));
    double bdotv = beta.Dot(SpaceVector());
    return LorentzVector(SpaceVector() + 
			 ((gamma-1)/beta2*bdotv-gamma*fVec[3]) * beta,
			 gamma*(fVec[3]-bdotv));
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T Beta2() const {
     return (SpaceVector<T>()/fVec[4]).Mag2();}

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T Beta() const {
     return Sqrt(Beta2());}

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T Gamma2() const {
     return 1/(1-Beta2());}

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T Gamma() const {
     return Sqrt(Gamma2());}

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T Rapidity() const {
     return 0.5*Log((fVec[3]+fVec[2])/(fVec[3]-fVec[2]));}

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  T PseudoRapidity() const {
     return -Log(Tan(0.5*Theta()));}

  /**
   * This method generates a two body decay 
   @param masses [input] a vector of lenght two with the masses of the products
   @param t1     [out]   LorentzVector of the first decay product
   @param t2     [out]   LorentzVector of the second decay product
   @param func   [in]    function defining the angular distribution of the particles if not uniform
  */
  template<typename U, typename V, typename W>
  VECGEOM_CUDA_HEADER_BOTH
  void TwoBodyDecay(const U masses[], LorentzVector<V> &t1, LorentzVector<W> &t2, void (*func)(double &mthe, double &mphi)=0) const {
   // Make two body decay
   double theta=0;
   double phi=0;
   if(func) func(theta,phi); 
      else {
	 phi = RNG::Instance().uniform() * kTwoPi;
	 theta = RNG::Instance().uniform() * kPi;
      }
   double mass2 = Mag2();
   if(mass2 > -1e-6 && mass2 < 0) mass2=0;
   double dmass = masses[0]+masses[1];
   if(mass2*(1.+1.e-7) < dmass*dmass || (mass2==0 && dmass==0)) {
      if(dmass !=0)
	 std::cout << "WARNING::LorentzVector::TwoBodyDecay: daughters' mass square " 
		   << dmass*dmass << " larger than parent " << mass2 << ", no decay " << std::endl;
      t1.Set(0,0,0,masses[0]);
      t2.Set(0,0,0,masses[1]);
      return;
   }
   const double mass = sqrt(mass2);
	     
   double en1 = 0.5*(mass*mass + masses[0]*masses[0] - masses[1]*masses[1])/mass;
   double en2 = 0.5*(mass*mass - masses[0]*masses[0] + masses[1]*masses[1])/mass;
   if((en1-masses[0]<-1e-10) || (en2-masses[1]<-1e-10)) {
      std::cout << "WARNING::LorentzVector::TwoBodyDecay en<0  en1-masses[0] " << en1-masses[0] 
	   << " en2-masses[1] " << en2-masses[1] << std::endl;
   }
   double p1=0;
   if(en1<masses[0]) 
      en1 = masses[0];
  else
      p1 = sqrt((en1+masses[0])*(en1-masses[0]));
   double p2=0;
   if(en2<masses[1]) 
      en2 = masses[1];
   else
      p2 = sqrt((en2+masses[1])*(en2-masses[1]));
   double p = 0.5*(p1+p2);
   // debug
   if(fabs(p1-p2)/mass > 1e-4) std::cout << "WARNING::LorentzVector::TwoBodyDecay: momentum imbalance err " 
					 << fabs(p1-p2)/mass << " p1-p2 " << p1-p2 << " p1 " << p1 << " p2 " << p2 << " mass " << mass << std::endl;
   if(fabs(mass-en1-en2)/(mass+en1+en2)/3 > 1e-7) std::cout << "WARNING::LorentzVector::TwoBodyDecay: mass imbalance mass - en1 - en2 " 
									<< mass - en1 - en2 << std::endl;
   //debug

   // Set the two vectors
   const Vector3D<T> space(p1*sin(theta)*cos(phi),p1*sin(theta)*sin(phi),p1*cos(theta));
   t1.Set( space, en1);
   t2.Set(-space, en2);

   if(func) {
   // This part only if there is a matrix element, otherwise isotropic is isotropic
      const Transformation3D tr(SpaceVector<T>(),false);
      tr.Transform(space,static_cast<Vector3D<T>& >(t1));
      tr.Transform(-space,static_cast<Vector3D<T>& >(t2));
   }
   
   Vector3D<T> beta(-fVec[0]/fVec[3],-fVec[1]/fVec[3],-fVec[2]/fVec[3]);
   t1 = t1.Boost(beta);
   t2 = t2.Boost(beta);
   // dedbug
   static const LorentzVector<T> metric(-1,-1,-1,1);
   const LorentzVector<T> vdiff = (*this) - t1 -t2;
   const double mdiff = sqrt(vdiff.Dot(vdiff*metric)/(Dot((*this)*metric)+t1.Dot(t1*metric)+t2.Dot(t2*metric)))/3;
   if(mdiff > 1e-7) std::cout << "WARNING::LorentzVector::TwoBodyDecay error " << mdiff << " in momentum energy balance this=" 
			      << (*this) << " difference=" << (*this) - t1 - t2 << std::endl;
   // debug
}

  /**
   * This method generate recursively the phase space via two body decays
   * The method is robust, however it provides a biased distribution
   * this will be fixed in due time, for the moment it is good enough
   @param[in] masses the masses of the products
   @param[out] daughters a vector of LorentzVectors with the decay products
  */
 template<typename U>
 VECGEOM_CUDA_HEADER_BOTH
 void PhaseSpace(const vector<U> &masses, vector<LorentzVector<T> > &daughters) {
   //Relativistic phase space
   LorentzVector<T> pa(*this);
   LorentzVector<T> d1;
   LorentzVector<T> d2;

   //
   int ndec = masses.size();
   daughters.clear();
   double mass2=0;
   for(int j=0; j<ndec; ++j) mass2+=masses[j];
   if(mass2>Mag()) {
      std::cout << "LorentzVector::PhaseSpace: cannot decay: parent mass " << Mag() << " sum of daughter mass " << mass2 << std::endl;
      return;
   }
   for(int j=0; j<ndec-1;++j) {
      mass2=masses[j+1];
      if(j<ndec-2) {
	 for(int i=j+2; i<ndec; ++i) mass2+=masses[i];
	 double dm = 0;
	 double pam2 = pa.Mag2();
	 if(pam2<0) pam2=0;
	 double freen = sqrt(pam2)-masses[j]-mass2;
	 if(mass2 > 0) dm = freen*std::pow(RNG::Instance().uniform(),masses[j]/(mass2-0.5*(masses[ndec-1]+masses[ndec-2])));
	 while(dm < 1e-7*freen) dm = freen*RNG::Instance().uniform();
	 mass2+=dm;
      }
      double vmass[2]={masses[j],mass2};
      pa.TwoBodyDecay(vmass,d1,d2);
      pa = d2;
      daughters.push_back(d1);
   }
   daughters.push_back(d2);
}

};

template <typename T>
std::ostream& operator<<(std::ostream& os, LorentzVector<T> const &vec) {
  os << "(" << vec[0] << ", " << vec[1] << ", " << vec[2] << ", " << vec[3] << ")";
  return os;
}

#ifdef VECGEOM_VC_ACCELERATION

/// This is a template specialization of class LorentzVector<double> or
/// LorentzVector<float> that can provide internal vectorization of common vector
/// operations.
template <>
class LorentzVector<Precision> : public AlignedBase {

  typedef LorentzVector<Precision> VecType;
  typedef LorentzVector<bool> BoolType;
  typedef Vc::Vector<Precision> Base_t;

private:

  Vc::Memory<Vc::Vector<Precision>, 4> fMem;

public:

  Precision* AsArray() {
    return &fMem[0];
  }

 LorentzVector(const Precision a, const Precision b, const Precision c, const Precision d) : fMem() {
    fMem[0] = a;
    fMem[1] = b;
    fMem[2] = c;
    fMem[3] = d;
  }

  // Performance issue in Vc with: mem = a;
  VECGEOM_INLINE
  LorentzVector(const Precision a) : LorentzVector(a, a, a, a) {}

  VECGEOM_INLINE
  LorentzVector() : LorentzVector(0, 0, 0, 0) {}

  template <typename U>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
     LorentzVector(Vector3D<U> const &other, Precision t) {
    fMem[0] = other[0];
    fMem[1] = other[1];
    fMem[2] = other[2];
    fMem[3] = t;
  }

  VECGEOM_INLINE
    LorentzVector(LorentzVector const &other) : fMem() {
    //for( int i=0; i < 1 + 3/Base_t::Size; i++ )
//      {
         //Base_t v1 = other.mem.vector(i);
         //this->mem.vector(i)=v1;
       //}
     fMem[0]=other.fMem[0];
     fMem[1]=other.fMem[1];
     fMem[2]=other.fMem[2];
     fMem[3]=other.fMem[3];
  }

  VECGEOM_INLINE
  LorentzVector & operator=( LorentzVector const & rhs )
   {
      if(this != &rhs) {
	 //for( int i=0; i < 1 + 3/Base_t::Size; i++ )
	 //{
         //Base_t v1 = rhs.mem.vector(i);
         //this->mem.vector(i)=v1;
	 //}
	 // the following line must not be used: this is a bug in Vc
	 // this->mem = rhs.mem;
	 fMem[0]=rhs.fMem[0];
	 fMem[1]=rhs.fMem[1];
	 fMem[2]=rhs.fMem[2];
	 fMem[3]=rhs.fMem[3];
      }
      return *this;
   }

  LorentzVector(std::string const &str) : fMem() {
    int begin = str.find("(")+1, end = str.find(",")-1;
    fMem[0] = std::atof(str.substr(begin, end-begin+1).c_str());
    begin = end + 2;
    end = str.find(",", begin)-1;
    fMem[1] = std::atof(str.substr(begin, end-begin+1).c_str());
    begin = end + 2;
    end = str.find(",", begin)-1;
    fMem[2] = std::atof(str.substr(begin, end-begin+1).c_str());
    begin = end + 2;
    end = str.find(")", begin)-1;
    fMem[3] = std::atof(str.substr(begin, end-begin+1).c_str());
  }

  operator Vector3D<Precision>&() {return reinterpret_cast<Vector3D<Precision>&>(*this);}

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision& operator[](const int index) {
    return fMem[index];
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  const Precision& operator[](const int index) const {
    return fMem[index];
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision& x() { return fMem[0]; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  const Precision& x() const { return fMem[0]; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision& y() { return fMem[1]; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  const Precision& y() const { return fMem[1]; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision& z() { return fMem[2]; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  const Precision& z() const { return fMem[2]; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision& t() { return fMem[3]; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  const Precision& t() const { return fMem[4]; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
    void Set(const Precision in_x, const Precision in_y, const Precision in_z, const Precision in_t) {
    fMem[0] = in_x;
    fMem[1] = in_y;
    fMem[2] = in_z;
    fMem[3] = in_t;
  }

  template <typename U, typename V>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
     void Set(const Vector3D<U> &v, V t) {
    fMem[0] = v[0];
    fMem[1] = v[1];
    fMem[2] = v[2];
    fMem[3] = t;
  }

  template <typename U>
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  void Set(const Vector3D<U> &v) {
    fMem[0] = v[0];
    fMem[1] = v[1];
    fMem[2] = v[2];
  }

  VECGEOM_CUDA_HEADER_BOTH
  void Set(const Precision in_x) {
     Set(in_x, in_x, in_x, in_x);
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision Mag2() const {
      return Dot(*this,*this);
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision Mag() const {
    return Sqrt(Mag2());
  }

  // TODO: study if we gain from internal vectorization here.
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision Perp2() const {
    return fMem[0]*fMem[0] + fMem[1]*fMem[1];
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision Perp() const {
    return Sqrt(Perp2());
  }

  VECGEOM_INLINE
  void Map(Precision (*f)(const Precision&)) {
    fMem[0] = f(fMem[0]);
    fMem[1] = f(fMem[1]);
    fMem[2] = f(fMem[2]);
    fMem[3] = f(fMem[3]);
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  void MaskedAssign(LorentzVector<bool> const &condition,
                    LorentzVector<Precision> const &value) {
    fMem[0] = (condition[0]) ? value[0] : fMem[0];
    fMem[1] = (condition[1]) ? value[1] : fMem[1];
    fMem[2] = (condition[2]) ? value[2] : fMem[2];
    fMem[3] = (condition[3]) ? value[3] : fMem[3];
  }

  /// \return Azimuthal angle between -pi and pi.
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision Phi() const {
    return (fMem[0] != 0. || fMem[1] != 0.) ? ATan2(fMem[1], fMem[0]) : 0.;
  }

  /// \return Polar angle between 0 and pi.
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision Theta() const {
     Precision theta = ATan2(Sqrt(fMem[0]*fMem[0]+fMem[1]*fMem[1]),fMem[2]);
     if(theta < 0) theta += kPi;
     return theta;
  }

  /// \return The dot product of two LorentzVector<Precision> objects.
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  static
  Precision Dot(LorentzVector<Precision> const &left,
                LorentzVector<Precision> const &right) {
    // TODO: This function should be internally vectorized (if proven to be
    //       beneficial)

    // To avoid to initialize the padding component, we can not use mem.vector's
    // multiplication and addition since it would accumulate also the (random) padding component
    return - left.fMem[0]*right.fMem[0] - left.fMem[1]*right.fMem[1] - left.fMem[2]*right.fMem[2] + left.fMem[3]*right.fMem[3];
  }

  /// \return The dot product with another LorentzVector<Precision> object.
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision Dot(LorentzVector<Precision> const &right) const {
    return Dot(*this, right);
  }

  // For UVector3 compatibility. Is equal to normal multiplication.
  // TODO: check if there are implicit dot products in USolids...
  VECGEOM_INLINE
  VecType MultiplyByComponents(VecType const &other) const;

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  LorentzVector<Precision>& FixZeroes() {
    for (int i = 0; i < 4; ++i) {
      if (std::abs(fMem.scalar(i)) < kTolerance) fMem.scalar(i) = 0;
    }
    return *this;
  }

  // Inplace binary operators

  #define LORENTZVECTOR_ACCELERATED_INPLACE_BINARY_OP(OPERATOR) \
  VECGEOM_CUDA_HEADER_BOTH \
  VecType& operator OPERATOR(const VecType &other) { \
    for (unsigned i = 0; i < 1 + 4/Vc::Vector<Precision>::Size; ++i) { \
      this->fMem.vector(i) OPERATOR other.fMem.vector(i); \
    } \
    return *this; \
  } \
  VECGEOM_CUDA_HEADER_BOTH \
  VecType& operator OPERATOR(const Precision &scalar) { \
    for (unsigned i = 0; i < 1 + 4/Vc::Vector<Precision>::Size; ++i) { \
      this->fMem.vector(i) OPERATOR scalar; \
    } \
    return *this; \
  }
  LORENTZVECTOR_ACCELERATED_INPLACE_BINARY_OP(+=)
  LORENTZVECTOR_ACCELERATED_INPLACE_BINARY_OP(-=)
  LORENTZVECTOR_ACCELERATED_INPLACE_BINARY_OP(*=)
  LORENTZVECTOR_ACCELERATED_INPLACE_BINARY_OP(/=)
  #undef LORENTZVECTOR_ACCELERATED_INPLACE_BINARY_OP

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
     static VecType FromCylindrical(Precision r, Precision phi, Precision z, Precision t) {
     return VecType(r*cos(phi), r*sin(phi), z, t);
  }

  template <typename U> 
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Vector3D<U> SpaceVector() const {
     return Vector3D<U>(fMem[0], fMem[1], fMem[2]);
  }

  template <typename U> 
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  VecType Boost(const Vector3D<U> &beta) const {
    Precision beta2 = beta.Mag2();
    if(beta2 <= 1e-16) return LorentzVector(*this);
    Precision gamma = Sqrt(1./(1-beta2));
    Precision bdotv = beta.Dot(SpaceVector<Precision>());
    return LorentzVector(SpaceVector<Precision>() + 
			 ((gamma-1)/beta2*bdotv-gamma*fMem[3]) * beta,
			 gamma*(fMem[3]-bdotv));
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision Beta2() const {
     return (SpaceVector<Precision>()/fMem[4]).Mag2();}

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision Beta() const {
     return Sqrt(Beta2());}

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision Gamma2() const {
     return 1/(1-Beta2());}

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision Gamma() const {
     return Sqrt(Gamma2());}

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision Rapidity() const {
     return 0.5*Log((fMem[3]+fMem[2])/(fMem[3]-fMem[2]));}

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  Precision PseudoRapidity() const {
     return -Log(tan(0.5*Theta()));}

  VECGEOM_CUDA_HEADER_BOTH
  VecType operator*(const VecType &other) const {
     return LorentzVector(fMem[0]*other.fMem[0], fMem[1]*other.fMem[1],
			  fMem[2]*other.fMem[2], fMem[3]*other.fMem[3]);
  }

  template <typename U, typename V, typename W> 
  VECGEOM_CUDA_HEADER_BOTH
  void TwoBodyDecay(const U masses[], LorentzVector<V> &t1, LorentzVector<W> &t2, void (*func)(double &mthe, double &mphi)=0) const {
   // Make two body decay
   double theta=0;
   double phi=0;
   if(func) func(theta,phi); 
      else {
	 phi = RNG::Instance().uniform() * kTwoPi;
	 theta = RNG::Instance().uniform() * kPi;
      }
   double mass2 = Mag2();
   if(mass2 > -1e-6 && mass2 < 0) mass2=0;
   double dmass = masses[0]+masses[1];
   if(mass2*(1.+1.e-7) < dmass*dmass || (mass2==0 && dmass==0)) {
      if(dmass !=0)
	 std::cout << "WARNING::LorentzVector::TwoBodyDecay: daughters' mass square " 
		   << dmass*dmass << " larger than parent " << mass2 << ", no decay " << std::endl;
      t1.Set(0,0,0,masses[0]);
      t2.Set(0,0,0,masses[1]);
      return;
   }
   const double mass = sqrt(mass2);
	     
   double en1 = 0.5*(mass*mass + masses[0]*masses[0] - masses[1]*masses[1])/mass;
   double en2 = 0.5*(mass*mass - masses[0]*masses[0] + masses[1]*masses[1])/mass;
   if((en1-masses[0]<-1e-10) || (en2-masses[1]<-1e-10)) {
      std::cout << "WARNING::LorentzVector::TwoBodyDecay en<0  en1-masses[0] " << en1-masses[0] 
	   << " en2-masses[1] " << en2-masses[1] << std::endl;
   }
   double p1=0;
   if(en1<masses[0]) 
      en1 = masses[0];
  else
      p1 = sqrt((en1+masses[0])*(en1-masses[0]));
   double p2=0;
   if(en2<masses[1]) 
      en2 = masses[1];
   else
      p2 = sqrt((en2+masses[1])*(en2-masses[1]));
   double p = 0.5*(p1+p2);
   // debug
   if(fabs(p1-p2)/mass > 1e-4) std::cout << "WARNING::LorentzVector::TwoBodyDecay: momentum imbalance err " 
					 << fabs(p1-p2)/mass << " p1-p2 " << p1-p2 << " p1 " << p1 << " p2 " << p2 << " mass " << mass << std::endl;
   if(fabs(mass-en1-en2)/(mass+en1+en2)/3 > 1e-7) std::cout << "WARNING::LorentzVector::TwoBodyDecay: mass imbalance mass - en1 - en2 " 
									<< mass - en1 - en2 << std::endl;
   // debug

   // Set the two output vectors
   const Vector3D<Precision> space(p*sin(theta)*cos(phi),p*sin(theta)*sin(phi),p*cos(theta));
   t1.Set( space, en1);
   t2.Set(-space, en2);

   if(func) {
   // This part only if there is a matrix element, otherwise isotropic is isotropic
      const Transformation3D tr(SpaceVector<Precision>(),false);
      tr.Transform(space,static_cast<Vector3D<Precision>& >(t1));
      tr.Transform(-space,static_cast<Vector3D<Precision>& >(t2));
   }
   
   Vector3D<Precision> beta(-fMem[0]/fMem[3],-fMem[1]/fMem[3],-fMem[2]/fMem[3]);
   t1 = t1.Boost(beta);
   t2 = t2.Boost(beta);

   // dedbug
   static const LorentzVector<Precision> metric(-1,-1,-1,1);
   const LorentzVector<Precision> vdiff = (*this) - t1 -t2;
   const double mdiff = sqrt(vdiff.Dot(vdiff*metric)/(Dot((*this)*metric)+t1.Dot(t1*metric)+t2.Dot(t2*metric)))/3;
   if(mdiff > 1e-7) std::cout << "WARNING::LorentzVector::TwoBodyDecay error " << mdiff << " in momentum energy balance this=" 
			      << (*this) << " difference=" << (*this) - t1 - t2 << std::endl;
   // debug
}

  /**
   * This method generate recursively the phase space via two body decays
   * The method is robust, however it provides a biased distribution
   * this will be fixed in due time, for the moment it is good enough
   @param[in] masses the masses of the products
   @param[out] daughters a vector of LorentzVectors with the decay products
  */
 template<typename U>
 VECGEOM_CUDA_HEADER_BOTH
 void PhaseSpace(const vector<U> &masses, vector<LorentzVector<Precision> > &daughters) {
   //Relativistic phase space
   LorentzVector<Precision> pa(*this);
   LorentzVector<Precision> d1;
   LorentzVector<Precision> d2;

   //
   int ndec = masses.size();
   daughters.clear();
   double mass2=0;
   for(int j=0; j<ndec; ++j) mass2+=masses[j];
   if(mass2>Mag()) {
      std::cout << "LorentzVector::PhaseSpace: cannot decay: parent mass " << Mag() << " sum of daughter mass " << mass2 << std::endl;
      return;
   }
   for(int j=0; j<ndec-1;++j) {
      mass2=masses[j+1];
      if(j<ndec-2) {
	 for(int i=j+2; i<ndec; ++i) mass2+=masses[i];
	 double dm = 0;
	 double pam2 = pa.Mag2();
	 if(pam2<0) pam2=0;
	 double freen = sqrt(pam2)-masses[j]-mass2;
	 if(mass2 > 0) dm = freen*std::pow(RNG::Instance().uniform(),masses[j]/(mass2-0.5*(masses[ndec-1]+masses[ndec-2])));
	 while(dm < 1e-7*freen) dm = freen*RNG::Instance().uniform();
	 mass2+=dm;
      }
      double vmass[2]={masses[j],mass2};
      pa.TwoBodyDecay(vmass,d1,d2);
      pa = d2;
      daughters.push_back(d1);
   }
   daughters.push_back(d2);
}

};
#else
//#pragma message "using normal LorentzVector.h"
#endif // VECGEOM_VC_ACCELERATION


#define LORENTZVECTOR_BINARY_OP(OPERATOR, INPLACE) \
template <typename T, typename V> \
VECGEOM_INLINE \
VECGEOM_CUDA_HEADER_BOTH \
LorentzVector<T> operator OPERATOR(const LorentzVector<T> &lhs, \
                                 const LorentzVector<V> &rhs) { \
  LorentzVector<T> result(lhs); \
  result INPLACE rhs; \
  return result; \
} \
template <typename T, typename V> \
VECGEOM_INLINE \
VECGEOM_CUDA_HEADER_BOTH \
LorentzVector<T> operator OPERATOR(LorentzVector<T> const &lhs, \
                                 const V rhs) { \
  LorentzVector<T> result(lhs); \
  result INPLACE rhs; \
  return result; \
} \
template <typename T, typename V> \
VECGEOM_INLINE \
VECGEOM_CUDA_HEADER_BOTH \
LorentzVector<T> operator OPERATOR(const V lhs, \
                                 LorentzVector<T> const &rhs) { \
  LorentzVector<T> result(lhs); \
  result INPLACE rhs; \
  return result; \
}
LORENTZVECTOR_BINARY_OP(+, +=)
LORENTZVECTOR_BINARY_OP(-, -=)
LORENTZVECTOR_BINARY_OP(*, *=)
LORENTZVECTOR_BINARY_OP(/, /=)
#undef LORENTZVECTOR_BINARY_OP

VECGEOM_INLINE
VECGEOM_CUDA_HEADER_BOTH
bool operator==(
    LorentzVector<Precision> const &lhs,
    LorentzVector<Precision> const &rhs) {
  return Abs(lhs[0] - rhs[0]) < kTolerance &&
         Abs(lhs[1] - rhs[1]) < kTolerance &&
         Abs(lhs[2] - rhs[2]) < kTolerance &&
         Abs(lhs[3] - rhs[3]) < kTolerance;
}

VECGEOM_INLINE
VECGEOM_CUDA_HEADER_BOTH
LorentzVector<bool> operator!=(
    LorentzVector<Precision> const &lhs,
    LorentzVector<Precision> const &rhs) {
  return !(lhs == rhs);
}

template <typename T>
VECGEOM_CUDA_HEADER_BOTH
VECGEOM_INLINE
LorentzVector<T> operator-(LorentzVector<T> const &vec) {
   return LorentzVector<T>(-vec[0], -vec[1], -vec[2], -vec[3]);
}

VECGEOM_CUDA_HEADER_BOTH
VECGEOM_INLINE
LorentzVector<bool> operator!(LorentzVector<bool> const &vec) {
   return LorentzVector<bool>(!vec[0], !vec[1], !vec[2], !vec[3]);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#define LORENTZVECTOR_SCALAR_BOOLEAN_LOGICAL_OP(OPERATOR) \
VECGEOM_CUDA_HEADER_BOTH \
VECGEOM_INLINE \
LorentzVector<bool> operator OPERATOR(LorentzVector<bool> const &lhs, \
                                 LorentzVector<bool> const &rhs) { \
  return LorentzVector<bool>(lhs[0] OPERATOR rhs[0], \
                        lhs[1] OPERATOR rhs[1], \
			lhs[2] OPERATOR rhs[2],	\
                        lhs[3] OPERATOR rhs[3]); \
}
LORENTZVECTOR_SCALAR_BOOLEAN_LOGICAL_OP(&&)
LORENTZVECTOR_SCALAR_BOOLEAN_LOGICAL_OP(||)
#undef LORENTZVECTOR_SCALAR_BOOLEAN_LOGICAL_OP
#pragma GCC diagnostic pop

#ifdef VECGEOM_VC

VECGEOM_INLINE
LorentzVector<VcBool> operator!(LorentzVector<VcBool> const &vec) {
   return LorentzVector<VcBool>(!vec[0], !vec[1], !vec[2], !vec[3]);
}

VECGEOM_INLINE
VcBool operator==(
    LorentzVector<VcPrecision> const &lhs,
    LorentzVector<VcPrecision> const &rhs) {
  return Abs(lhs[0] - rhs[0]) < kTolerance &&
         Abs(lhs[1] - rhs[1]) < kTolerance &&
         Abs(lhs[2] - rhs[2]) < kTolerance &&
         Abs(lhs[3] - rhs[3]) < kTolerance;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#define LORENTZVECTOR_VC_BOOLEAN_LOGICAL_OP(OPERATOR) \
VECGEOM_INLINE \
LorentzVector<VcBool> operator OPERATOR( \
    LorentzVector<VcBool> const &lhs, \
    LorentzVector<VcBool> const &rhs) { \
  return LorentzVector<VcBool>(lhs[0] OPERATOR rhs[0], \
                          lhs[1] OPERATOR rhs[1], \
                          lhs[2] OPERATOR rhs[2], \
                          lhs[3] OPERATOR rhs[3]); \
}
LORENTZVECTOR_VC_BOOLEAN_LOGICAL_OP(&&)
LORENTZVECTOR_VC_BOOLEAN_LOGICAL_OP(||)
#undef LORENTZVECTOR_VC_BOOLEAN_LOGICAL_OP
#pragma GCC diagnostic pop

#endif // VECGEOM_VC


#ifdef VECGEOM_VC_ACCELERATION

LorentzVector<Precision> LorentzVector<Precision>::MultiplyByComponents(
    VecType const &other) const {
  return (*this) * other;
}

#endif

static const LorentzVector<double> LorentzMetric(-1,-1,-1,1);


} // End inline namespace

} // End global namespace

#endif // VECGEOM_BASE_LORENTZVECTOR_H_
