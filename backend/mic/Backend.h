/// \file mic/backend.h

#ifndef VECGEOM_BACKEND_MICBACKEND_H_
#define VECGEOM_BACKEND_MICBACKEND_H_

#include "base/Global.h"
#include "backend/scalar/Backend.h"

#if defined(DEBUG) || defined(NDEBUG) || defined(_DEBUG)
#include <iostream>
#endif
#include <mic/micvec.h>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

class MicMask;
class MicIntegerVector;
class MicDoubleVector;

struct kMic {
  typedef MicIntegerVector int_v;
  typedef MicDoubleVector precision_v;
  typedef MicMask bool_v;
  typedef MicIntegerVector inside_v;
  constexpr static bool early_returns = false;
  const static precision_v kOne;
  const static precision_v kZero;
  const static bool_v kTrue;
  const static bool_v kFalse;
  // alternative typedefs ( might supercede above typedefs )
  typedef MicIntegerVector Int_t;
  typedef MicDoubleVector Double_t;
  typedef MicMask Bool_t;
  typedef MicDoubleVector Index_t;
};

#ifdef kVectorSize
#undef kVectorSize
#endif
constexpr int kVectorSize = 8;
#define VECGEOM_BACKEND_TYPE         kMic
#define VECGEOM_BACKEND_PRECISION    MicPrecision
#define VECGEOM_BACKEND_BOOL         MicBool
#define VECGEOM_BACKEND_INSIDE       kMic::inside_v

typedef kMic::int_v       MicInt;
typedef kMic::precision_v MicPrecision;
typedef kMic::bool_v      MicBool;
typedef kMic::inside_v    MicInside;

// Class, operators and functions for Mask

class MicMask : public VecMask16 {

public:
  MicMask() : VecMask16() {}
  MicMask(__mmask mm) : VecMask16(mm) {}
  MicMask(int mm) : VecMask16(mm) {}
  MicMask(VecMask16 mm) : VecMask16(mm) {}

  VECGEOM_INLINE
  bool operator[](size_t index) const { return static_cast<bool>(m & (1 << index)); }

#if defined(_IOSTREAM_) || defined(_CPP_IOSTREAM) || defined(_GLIBCXX_IOSTREAM)
  friend std::ostream& operator <<(std::ostream &os, const MicMask &a) {
    int size = 8;
    unsigned int num = a;
    int i;
    os << "{";
    for(i=0;i<size;i++) {
      os << (num&0x01);
      num = num >> 1;
    }
    os << "}";
    return os;
  }
#endif
};

VECGEOM_INLINE
MicMask operator !(MicMask const &val) {
  MicMask r(val);
  return ~r;
}

VECGEOM_INLINE
MicMask operator &&(MicMask const &val1, MicMask const &val2) {
  return val1 & val2;
}

VECGEOM_INLINE
MicMask operator &&(bool const &val1, MicMask const &val2) {
  return MicMask(val1) & val2;
}

VECGEOM_INLINE
MicMask operator &&(MicMask const &val1, bool const &val2) {
  return val1 & MicMask(val2);
}

VECGEOM_INLINE
bool IsFull(MicMask const &cond) {
  return  _mm512_kortestc(cond, cond);
}

// Class, operators and functions for Integer

class MicIntegerVector : public Is32vec16 {

public:
  MicIntegerVector() : Is32vec16() {}
  MicIntegerVector(__m512i mm) : Is32vec16(mm) {}
  MicIntegerVector(const int i) : Is32vec16(i) {}
  MicIntegerVector(const int i[16]) : Is32vec16((int*)i) {}
  MicIntegerVector(Is32vec16 mm) : MicIntegerVector((__m512i)mm) {}

};

// Class, operators and functions for Double

class MicDoubleVector : public F64vec8 {

public:
  MicDoubleVector() : F64vec8() {}
  MicDoubleVector(__m512d m) : F64vec8(m) {}
  MicDoubleVector(const double d) : F64vec8(d) {}
  MicDoubleVector(const double d[8]) : F64vec8((double*)d) {}
  MicDoubleVector(F64vec8 m) : MicDoubleVector((__m512d)m) {}

  VECGEOM_INLINE
  MicDoubleVector operator =(Precision const &val) { return MicDoubleVector(vec = _mm512_set1_pd(val)); }

  VECGEOM_INLINE
  MicDoubleVector operator =(MicDoubleVector const &val) { return MicDoubleVector(vec = val.vec); }

  VECGEOM_INLINE
  MicDoubleVector operator -() const { return MicPrecision(0.0) - vec; }

  VECGEOM_INLINE
  MicDoubleVector operator -(Precision const &val) const { return (*this) - MicDoubleVector(val); }

};

// Operators and Functions for MicDoubleVector/MicPrecision

VECGEOM_INLINE
MicMask operator ==(MicPrecision const &val1, MicPrecision const &val2) {
  return cmpeq(val1,val2);
}

VECGEOM_INLINE
MicMask operator ==(MicPrecision const &val1, Precision const &val2) {
  return cmpeq(val1,MicPrecision(val2));
}

VECGEOM_INLINE
MicMask operator !=(MicPrecision const &val1, MicPrecision const &val2) {
  return val1 != val2;
}

VECGEOM_INLINE
MicMask operator !=(MicPrecision const &val1, Precision const &val2) {
  return val1 != MicPrecision(val2);
}

VECGEOM_INLINE
MicMask operator >(MicPrecision const &val1, MicPrecision const &val2) {
  return cmpnle(val1, val2);
}

VECGEOM_INLINE
MicMask operator >(MicPrecision const &val1, Precision const &val2) {
  return val1 > MicPrecision(val2);
}

VECGEOM_INLINE
MicMask operator >=(MicPrecision const &val1, MicPrecision const &val2) {
  return cmpnlt(val1, val2);
}

VECGEOM_INLINE
MicMask operator >=(MicPrecision const &val1, Precision const &val2) {
  return val1 >= MicPrecision(val2);
}

VECGEOM_INLINE
MicMask operator <(MicPrecision const &val1, MicPrecision const &val2) {
  return cmplt(val1, val2);
}

VECGEOM_INLINE
MicMask operator <(MicPrecision const &val1, double const &val2) {
  return val1 < MicPrecision(val2);
}

VECGEOM_INLINE
MicMask operator <=(MicPrecision const &val1, MicPrecision const &val2) {
  return cmple(val1, val2);
}

VECGEOM_INLINE
MicMask operator <=(MicPrecision const &val1, Precision const &val2) {
  return val1 <= MicPrecision(val2);
}

VECGEOM_INLINE
MicMask operator <=(Precision const &val1, MicPrecision const &val2) {
  return MicPrecision(val1) <= val2;
}

VECGEOM_INLINE
MicPrecision operator -(double const &val1, MicPrecision const &val2) {
  return MicPrecision(MicPrecision((double)val1) - val2);
}

VECGEOM_INLINE
MicPrecision operator +(MicPrecision const &val1, Precision const &val2) {
  return val1 + MicPrecision(val2);
}

VECGEOM_INLINE
MicPrecision operator +(Precision const &val1, MicPrecision const &val2) {
  return MicPrecision(val1) + val2;
}

VECGEOM_INLINE
MicPrecision operator *(Precision const &val1, MicPrecision const &val2) {
  return MicPrecision(val1) * val2;
}

VECGEOM_INLINE
MicPrecision operator /(double const &val1, MicPrecision const &val2) {
  return MicPrecision((double)val1) / val2;
}

VECGEOM_INLINE
MicPrecision Abs(F64vec8 const &val) {
  MicPrecision _v(-0.0);
  return val & (val ^_v);
}

VECGEOM_INLINE
MicPrecision Abs(MicPrecision const &val) {
  MicPrecision _v(-0.0);
  return val & (val ^ _v);
}

VECGEOM_INLINE
MicPrecision Sqrt(MicPrecision const &val) {
  return sqrt(val);
}

VECGEOM_INLINE
MicPrecision Sqrt(F64vec8 const &val) {
  return MicPrecision(sqrt(val));
}

VECGEOM_INLINE
MicPrecision Max(MicPrecision const &val1, MicPrecision const &val2) {
  return max(val1, val2);
}

VECGEOM_INLINE
MicPrecision Min(MicPrecision const &val1, MicPrecision const &val2) {
  return min(val1, val2);
}

VECGEOM_INLINE
MicPrecision ATan2(MicPrecision const &y, MicPrecision const &x) {
  return atan2(x, y);
}

// MaskedAssign

VECGEOM_INLINE
void MaskedAssign(MicBool const &cond,
                  MicInside const &thenval,
                  MicInside * const output) {
  (*output) = _mm512_mask_or_epi32(*output, cond, thenval, thenval);
}

VECGEOM_INLINE
void MaskedAssign(MicBool const &cond,
                  Inside_t const &thenval,
                  MicInside * const output) {
  MaskedAssign(cond, MicInside(thenval), output);
}

VECGEOM_INLINE
void MaskedAssign(MicBool const &cond,
                  MicBool const &thenval,
                  MicBool * const output) {
  MicBool not_cond(~cond);
  MicBool p1(not_cond & (*output));
  MicBool p2(cond & thenval);
  (*output) = p1 | p2;
}

VECGEOM_INLINE
void MaskedAssign(MicBool const &cond,
                  MicPrecision const &thenval,
                  MicPrecision * const output) {
  (*output) = _mm512_castsi512_pd(_mm512_mask_or_epi64(_mm512_castpd_si512(*output), cond, _mm512_castpd_si512(thenval), _mm512_castpd_si512(thenval)));
}

VECGEOM_INLINE
void MaskedAssign(MicBool const &cond,
                  MicPrecision const &thenval,
                  double *const output) {
  MicPrecision _m(output);
  _m = _mm512_castsi512_pd(_mm512_mask_or_epi64(_mm512_castpd_si512(_m), cond, _mm512_castpd_si512(thenval), _mm512_castpd_si512(thenval)));
  _mm512_store_pd(output,_m);
}

} // End inline namespace

} // End global namespace


#endif // VECGEOM_BACKEND_VCBACKEND_H_
