find_package(CUDA REQUIRED)

set(VECCORE_ENABLE_CUDA True)

if(CMAKE_CXX_STANDARD STREQUAL 11)
  list(APPEND CUDA_NVCC_FLAGS -std=c++11)
else()
  message(FATAL_ERROR "CUDA compilation supports only ISO C++ 2011 standard")
endif()

set(CUDA_NVCC_FLAGS_DEBUG          "-O2 -g -G")
set(CUDA_NVCC_FLAGS_RELEASE        "-O3 -use_fast_math")
set(CUDA_NVCC_FLAGS_RELWITHDEBINFO "-O3 -g -G -use_fast_math")

list(APPEND CUDA_NVCC_FLAGS "-gencode arch=compute_30,code=sm_30")

message(STATUS "Compiling with CUDA enabled")
