//
// File:    test/core/TestEarlyReturns.cpp
// Purpose: Tests for EarlyReturnAllowed() and EarlyReturnMaxLength() constexpr's
//

//-- ensure asserts are compiled in
#undef NDEBUG

#include <iostream>
//#include "Vc/Vc"
#include "VecCore/VecCore"
#include "backend/Backend.h"

using namespace vecCore;

template <class Backend>
VECCORE_FORCE_NOINLINE
void TestType(const char *type)
{
  using Real_v = typename Backend::Float_v;
  std::cout << "Testing type " << type << ", VectorSize=" << VectorSize<Real_v>() << "\n";

  constexpr size_t vsize = VectorSize<Real_v>();
  Real_v flagReal_v;

  // returns false for vsize-1
  if (!EarlyReturnMaxLength(flagReal_v, vsize - 1))
    assert(true);
  else
    assert(false);

  // returns true for vsize
  if (EarlyReturnMaxLength(flagReal_v, vsize))
    assert(true);
  else
    assert(false);
}

int main()
{
#ifdef VECGEOM_NVCC
  // always false for CUDA
  std::cout << "Testing early returns for CUDA=ON...\n";
  assert(!EarlyReturnAllowed());
  double flagDouble;
  assert(!EarlyReturnMaxLength(, 1);
#else
  std::cout << "Testing early returns for CUDA=OFF...\n";

  // always true for non-CUDA
  assert(EarlyReturnAllowed());

  // always returns true for floats or floats
  std::cout << "Testing early returns for float...\n";
  float flagFloat;
  if (EarlyReturnMaxLength(flagFloat, 1))
    assert(true);
  else
    assert(false);
  if (EarlyReturnMaxLength(flagFloat, 2))
    assert(true);
  else
    assert(false);

  std::cout << "Testing early returns for double...\n";
  double flagDouble;
  if (EarlyReturnMaxLength(flagDouble, 1))
    assert(true);
  else
    assert(false);
  if (EarlyReturnMaxLength(flagDouble, 2))
    assert(true);
  else
    assert(false);

#ifdef VECGEOM_ENABLE_VC
  // tests for Float_v
  using Float_v = backend::VcVector::Float_v;
  Float_v flagFloat_v;
  std::cout << "Testing early returns for Float_v of size " << 8 * sizeof(Float_v)
            << " bits..., VectorSize=" << VectorSize<Float_v>() << "\n";
  constexpr size_t vsize1 = VectorSize<Float_v>();
  if (!EarlyReturnMaxLength(flagFloat_v, vsize1 - 1)
    assert(true);
  else
    assert(false);
  if (EarlyReturnMaxLength(flagFloat_v, vsize1)
    assert(true);
  else
    assert(false);

  // tests for Double_v
  using Double_v = backend::VcVector::Double_v;
  Double_v flagDouble_v;
  std::cout << "Testing early returns for Double_v of size " << 8 * sizeof(Double_v)
            << " bits..., VectorSize=" << VectorSize<Double_v>() << "\n";
  constexpr size_t vsize2 = VectorSize<Double_v>();
  if (!EarlyReturnMaxLength(flagDouble_v, vsize2 - 1)
    assert(true);
  else
    assert(false);
  if (EarlyReturnMaxLength(flagDouble_v, vsize2)
    assert(true);
  else
    assert(false);
#endif

  TestType<backend::Scalar>("Scalar");

#ifdef VECCORE_ENABLE_VC
  TestType<backend::VcScalar>("VcScalar");
  TestType<backend::VcVector>("VcVector");
  TestType<backend::VcSimdArray<8>>("VcSimdArray<8>");
  TestType<backend::VcSimdArray<16>>("VcSimdArray<16>");
  TestType<backend::VcSimdArray<32>>("VcSimdArray<32>");
#endif

#ifdef VECCORE_ENABLE_UMESIMD
  TestType<backend::UMESimd>("UME::SIMD");
  TestType<backend::UMESimdArray<8>>("UME::SIMD<8>");
  TestType<backend::UMESimdArray<16>>("UME::SIMD<16>");
  TestType<backend::UMESimdArray<32>>("UME::SIMD<32>");
#endif

#endif

  std::cout << "TestEarlyReturns: passed\n";
  return 0;
}
