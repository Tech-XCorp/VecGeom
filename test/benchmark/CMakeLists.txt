
#define some benchmark executables

set(TEST_EXECUTABLES_BENCHMARK
  ${CMAKE_SOURCE_DIR}/test/benchmark/BoxBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/SExtruBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/ConcaveSExtruBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/ParaboloidBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/ParaboloidScriptBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/ParallelepipedBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/PolyhedronBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/TubeBenchmark.cpp
#  ${CMAKE_SOURCE_DIR}/test/benchmark/BoxMinusHollowTubeBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/BoxMinusHollowTubeBenchmark_virtual.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/BoxUnionHollowTubeBenchmark_virtual.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/BoxIntersectHollowTubeBenchmark_virtual.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/BoxUnionBoxBenchmark_virtual.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/3LevelSubtractionBenchmark_virtual.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/TorusBenchmark2.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/TrapezoidBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/TrapezoidBenchmarkScript.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/OrbBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/SphereBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/HypeBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/TrdBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/ConeBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/PolyconeBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/ScaledBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/BoxScaledBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/GenTrapBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/GenTrapBenchmark_planar.cpp
  PARENT_SCOPE
)

#define which ones should be run under ctest

set(THISCTESTS
# TESTS MARKED # FAIL AND NEED ATTENTION
  ${CMAKE_SOURCE_DIR}/test/benchmark/BoxBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/SExtruBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/ConcaveSExtruBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/ParaboloidBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/ParaboloidScriptBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/ParallelepipedBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/PolyhedronBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/TubeBenchmark.cpp
# ${CMAKE_SOURCE_DIR}/test/benchmark/BoxMinusHollowTubeBenchmark_virtual.cpp
# ${CMAKE_SOURCE_DIR}/test/benchmark/BoxUnionHollowTubeBenchmark_virtual.cpp
# ${CMAKE_SOURCE_DIR}/test/benchmark/BoxIntersectHollowTubeBenchmark_virtual.cpp
# ${CMAKE_SOURCE_DIR}/test/benchmark/BoxUnionBoxBenchmark_virtual.cpp
# ${CMAKE_SOURCE_DIR}/test/benchmark/3LevelSubtractionBenchmark_virtual.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/TorusBenchmark2.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/TrapezoidBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/TrapezoidBenchmarkScript.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/OrbBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/SphereBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/HypeBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/TrdBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/ConeBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/GenTrapBenchmark.cpp
  ${CMAKE_SOURCE_DIR}/test/benchmark/GenTrapBenchmark_planar.cpp
)

# these tests fail because of a problem in USOLIDS:
IF(NOT USOLIDS)
  set(THISCTESTS ${THISCTESTS}
      ${CMAKE_SOURCE_DIR}/test/benchmark/PolyconeBenchmark.cpp
  )
ENDIF()

# these fail on CUDA
IF(NOT CUDA)
  set(THISCTESTS ${THISCTESTS}
      ${CMAKE_SOURCE_DIR}/test/benchmark/ScaledBenchmark.cpp
      ${CMAKE_SOURCE_DIR}/test/benchmark/BoxScaledBenchmark.cpp
  )
ENDIF()

# set the variable use in the parent cmake
set(CTESTS_BENCHMARK ${THISCTESTS} PARENT_SCOPE)