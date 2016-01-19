#ifndef VECGEOM_CUDA_H
#define VECGEOM_CUDA_H

#if (defined(__CUDACC__) || defined(__NVCC__))
  // Compiling with nvcc
  #define VECGEOM_NVCC
  #ifdef __CUDA_ARCH__
    // Compiling device code
    #define VECGEOM_NVCC_DEVICE
  #endif
  #define HAVECUDANAMESPACE
  #define VECGEOM_IMPL_NAMESPACE cuda
  #define VECGEOM_NAMESPACE ::vecgeom
  #define VECGEOM_CUDA_HEADER_HOST __host__
  #define VECGEOM_CUDA_HEADER_DEVICE __device__
  #define VECGEOM_CUDA_HEADER_BOTH __host__ __device__
  #define VECGEOM_CUDA_HEADER_GLOBAL __global__
  #define VECGEOM_ALIGNED __align__((64))
  #define VECGEOM_HOST_FORWARD_DECLARE(X) namespace cxx { X }
  #define VECGEOM_DEVICE_FORWARD_DECLARE(X)
  #define VECGEOM_DEVICE_DECLARE_CONV(X)
  #define VECGEOM_DEVICE_DECLARE_NS_CONV(NS,X,Def)
  #define VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE(X,ArgType)
  #define VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_2v(X,ArgType1,Def1,ArgType2,Def2)
  #define VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_2v_1t(X,ArgType1,Def1,ArgType2,Def2,ArgType3)
  #define VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_3v(X,ArgType1,Def1,ArgType2,Def2,ArgType3,Def3)
  #define VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_4v(X,ArgType1,Def1,ArgType2,Def2,ArgType3,Def3,ArgType4,Def4)
  #undef VECGEOM_VC
  #undef VECGEOM_CILK
  #undef VECGEOM_ROOT
  #undef VECGEOM_USOLIDS
  #undef VECGEOM_GEANT4
  #undef VECGEOM_BENCHMARK
#else
  // Not compiling with NVCC
#define HAVENORMALNAMESPACE
  #define VECGEOM_IMPL_NAMESPACE cxx
  #define VECGEOM_NAMESPACE ::vecgeom
  #define VECGEOM_CUDA_HEADER_HOST
  #define VECGEOM_CUDA_HEADER_DEVICE
  #define VECGEOM_CUDA_HEADER_BOTH
  #define VECGEOM_CUDA_HEADER_GLOBAL
  #ifdef VECGEOM_CUDA
    // CUDA is enabled, but currently compiling regular C++ code.
    // This enables methods that interface between C++ and CUDA environments
    #define VECGEOM_CUDA_INTERFACE
  #endif
  namespace vecgeom {
     template <typename DataType> struct kCudaType;
     template <typename DataType> using CudaType_t = typename kCudaType<DataType>::type_t;
     template <> struct kCudaType<float> { using type_t = float; };
     template <> struct kCudaType<double> { using type_t = double; };
     template <> struct kCudaType<int> { using type_t = int; };
  }
  #define VECGEOM_HOST_FORWARD_DECLARE(X)
  #define VECGEOM_DEVICE_FORWARD_DECLARE(X)  namespace cuda { X }

  #define VECGEOM_DEVICE_DECLARE_CONV(X) \
     namespace cuda { class X; } \
     inline namespace cxx  { class X; } \
     template <> struct kCudaType<cxx::X> { using type_t = cuda::X; };

  #define VECGEOM_DEVICE_DECLARESTRUCT_CONV(X) \
     namespace cuda { struct X; } \
     inline namespace cxx  { struct X; } \
     template <> struct kCudaType<cxx::X> { using type_t = cuda::X; };

  #define VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE(X,ArgType) \
     namespace cuda { template <ArgType Arg> class X; } \
     inline namespace cxx  { template <ArgType Arg> class X; } \
     template <ArgType Arg> struct kCudaType<cxx::X<Arg> > \
     { using type_t = cuda::X<CudaType_t<Arg> >; };

#ifdef VECGEOM_CUDA_VOLUME_SPECIALIZATION

  #define VECGEOM_DEVICE_DECLARE_NS_CONV(NS,X,Def)     \
     namespace cuda { namespace NS { struct X; } } \
     inline namespace cxx { namespace NS { struct X; } } \
     template <> struct kCudaType<cxx::NS::X> { using type_t = cuda::NS::X; };

  #define VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_2v(X,ArgType1,Def1,ArgType2,Def2) \
     namespace cuda { template <ArgType1 Arg1,ArgType2 Arg2> struct X; } \
     inline namespace cxx  { template <ArgType1 Arg1,ArgType2 Arg2> struct X; } \
     template <ArgType1 Arg1,ArgType2 Arg2> struct kCudaType<cxx::X<Arg1,Arg2> > \
     { using type_t = cuda::X<Arg1,Arg2 >; };
  #define VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_2v_1t(X,ArgType1,Def1,ArgType2,Def2,ArgType3) \
     namespace cuda { template <ArgType1 Arg1,ArgType2 Arg2,ArgType3 Arg3> struct X; } \
     inline namespace cxx  { template <ArgType1 Arg1,ArgType2 Arg2,ArgType3 Arg3> struct X; } \
     template <ArgType1 Arg1,ArgType2 Arg2,ArgType3 Arg3> struct kCudaType<cxx::X<Arg1,Arg2,Arg3> > \
     { using type_t = cuda::X<Arg1, Arg2, CudaType_t<Arg3> >; };
  #define VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_3v(X,ArgType1,Def1,ArgType2,Def2,ArgType3,Def3) \
     namespace cuda { template <ArgType1 Arg1,ArgType2 Arg2,ArgType3 Arg3> struct X; } \
     inline namespace cxx  { template <ArgType1 Arg1,ArgType2 Arg2,ArgType3 Arg3> struct X; } \
     template <ArgType1 Arg1,ArgType2 Arg2,ArgType3 Arg3> struct kCudaType<cxx::X<Arg1,Arg2,Arg3> > \
     { using type_t = cuda::X<Arg1,Arg2,Arg3 >; };

#define VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_4v(X,ArgType1,Def1,ArgType2,Def2,ArgType3,Def3,ArgType4,Def4) \
    namespace cuda { template <ArgType1 Arg1,ArgType2 Arg2,ArgType3 Arg3,ArgType4 Arg4> struct X; } \
    inline namespace cxx  { template <ArgType1 Arg1,ArgType2 Arg2,ArgType3 Arg3,ArgType4 Arg4> struct X; } \
    template <ArgType1 Arg1,ArgType2 Arg2,ArgType3 Arg3,ArgType4 Arg4> struct kCudaType<cxx::X<Arg1,Arg2,Arg3,Arg4> > \
    { using type_t = cuda::X<Arg1,Arg2,Arg3,Arg4 >; };

#else // VECGEOM_CUDA_VOLUME_SPECIALIZATION

  #define VECGEOM_DEVICE_DECLARE_NS_CONV(NS,X,Def)     \
     namespace cuda { namespace NS { struct Def; } } \
     inline namespace cxx { namespace NS { struct X; } } \
     template <> struct kCudaType<cxx::NS::X> { using type_t = cuda::NS::Def; };

  #define VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_2v(X,ArgType1,Def1,ArgType2,Def2) \
     namespace cuda { template <ArgType1 Arg1,ArgType2 Arg2> struct X; } \
     inline namespace cxx  { template <ArgType1 Arg1,ArgType2 Arg2> struct X; } \
     template <ArgType1 Arg1,ArgType2 Arg2> struct kCudaType<cxx::X<Arg1,Arg2> > \
     { using type_t = cuda::X<Def1, Def2 >; };
  #define VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_2v_1t(X,ArgType1,Def1,ArgType2,Def2,ArgType3) \
     namespace cuda { template <ArgType1 Arg1,ArgType2 Arg2,ArgType3 Arg3> struct X; } \
     inline namespace cxx  { template <ArgType1 Arg1,ArgType2 Arg2,ArgType3 Arg3> struct X; } \
     template <ArgType1 Arg1,ArgType2 Arg2,ArgType3 Arg3> struct kCudaType<cxx::X<Arg1,Arg2,Arg3> > \
     { using type_t = cuda::X<Def2, Def2, CudaType_t<Arg3> >; };
  #define VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_3v(X,ArgType1,Def1,ArgType2,Def2,ArgType3,Def3) \
     namespace cuda { template <ArgType1 Arg1,ArgType2 Arg2,ArgType3 Arg3> struct X; } \
     inline namespace cxx  { template <ArgType1 Arg1,ArgType2 Arg2,ArgType3 Arg3> struct X; } \
     template <ArgType1 Arg1,ArgType2 Arg2,ArgType3 Arg3> struct kCudaType<cxx::X<Arg1,Arg2,Arg3> > \
     { using type_t = cuda::X<Def1,Def2,Def3 >; };
#define VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_4v(X,ArgType1,Def1,ArgType2,Def2,ArgType3,Def3,ArgType4,Def4) \
     namespace cuda { template <ArgType1 Arg1,ArgType2 Arg2,ArgType3 Arg3,ArgType4 Arg4> struct X; } \
     inline namespace cxx  { template <ArgType1 Arg1,ArgType2 Arg2,ArgType3 Arg3,ArgType4 Arg4> struct X; } \
     template <ArgType1 Arg1,ArgType2 Arg2,ArgType3 Arg3,ArgType4 Arg4> struct kCudaType<cxx::X<Arg1,Arg2,Arg3,Arg4> > \
     { using type_t = cuda::X<Def1,Def2,Def3,Def4 >; };

#endif // VECGEOM_CUDA_VOLUME_SPECIALIZATION

/* Instead of multiple macro, when we have auto expansion of Template pack we could use:
template <typename... Arguments>
struct kCudaType<cxx::BoxImplementation<Arguments...>  >
   { using type_t = typename cuda::BoxImplementation<CudaType_t<Arguments...> >; };
*/
#endif

#endif
