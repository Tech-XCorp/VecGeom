/// @file Benchmarker.h
/// @author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#ifndef VECGEOM_BENCHMARKING_BENCHMARKER_H_
#define VECGEOM_BENCHMARKING_BENCHMARKER_H_


#include "base/global.h"

#include "base/soa3d.h"
#include "benchmarking/BenchmarkResult.h"
#include "management/volume_pointers.h"
#include "volumes/placed_volume.h"

#include <list>

namespace vecgeom {

/// \brief Benchmarks geometry methods of arbitrary volumes for different
///        backends and compares to any included external libraries.
///
/// In order to run a benchmark, a world volume must be provided to the
/// benchmarker. This volume must have an available bounding box, and can
/// contain any number of daughter volumes. When the benchmark is run, points
/// will be sampled with a bias in regard to these daughter volumes. Deeper
/// hierarchies are not considered. For any level of verbosity above zero, the
/// benchmarker will output results to standard output. However, result structs
/// are available and can be retrieved, containing all information related to a
/// specific run.
/// \sa BenchmarkResult
class Benchmarker {

private:

  VPlacedVolume const *fWorld;
  unsigned fPointCount;
  unsigned fPoolMultiplier;
  unsigned fRepetitions;
  std::list<VolumePointers> fVolumes;
  std::list<BenchmarkResult> fResults;
  int fVerbosity;
  double fToInBias, fInsideBias;
  SOA3D<Precision> *fPointPool;
  SOA3D<Precision> *fDirectionPool;
  Precision *fStepMax;

public:

  /// \brief Runs all geometry benchmarks.
  void RunBenchmark();

  /// \brief Runs a benchmark of the Inside method.
  ///
  /// The fraction of sampled points that will be located inside of daughter
  /// volume is specified by calling SetInsideBias().
  /// \sa SetInsideBias(const double)
  void RunInsideBenchmark();

  /// \brief Runs a benchmark of the DistanceToIn and SafetyToIn methods.
  ///
  /// The fraction of sampled points that should be hitting a daughter volume is
  /// specified by calling SetToInBias().
  /// \sa SetToInBias(const double)
  void RunToInBenchmark();

  /// \brief Runs a benchmark of the DistanceToOut and SafetyToOut methods.
  void RunToOutBenchmark();

  /// \param world Mother volume containing daughters that will be benchmarked.
  ///              The mother volume must have an available bounding box, as it
  ///              is used in the sampling process.
  Benchmarker(VPlacedVolume const *const world);

  ~Benchmarker();

  /// \return Amount of points and directions sampled for each benchmark
  ///         iteration.
  unsigned GetPointCount() const { return fPointCount; }

  /// \return Bias with which directions for DistanceToIn and SafetyToIn are
  ///         sampled.
  double GetToInBias() const { return fToInBias; }

  /// \return Bias with which the points for Inside are sampled.
  double GetInsideBias() const { return fInsideBias; }

  /// \return Multiplier of point and direction pool to use for simulating
  ///         random memory access.
  unsigned GetPoolMultiplier() const { return fPoolMultiplier; }

  /// \return Level of verbosity to standard output.
  int GetVerbosity() const { return fVerbosity; }

  /// \return Amount of iterations the benchmark is run for.
  unsigned GetRepetitions() const { return fRepetitions; }

  /// \return World whose daughters are benchmarked.
  VPlacedVolume const* GetWorld() const { return fWorld; }

  /// \param pointCount Amount of points to benchmark in each iteration.
  void SetPointCount(const unsigned pointCount) { fPointCount = pointCount; }

  /// \param toInBias Fraction of sampled particles that should be facing a
  /// daughter volume.
  void SetToInBias(const double toInBias) { fToInBias = toInBias; }

  /// \param insideBias Fraction of sampled particles that should be contained
  ///                   in a daughter volume.
  void SetInsideBias(const double insideBias) { fInsideBias = insideBias; }

  /// \param Multiplier for the pool of sampled points and directions.
  ///
  /// Can be increased to simulate more random access of memory, but will
  /// disable comparison of output.     
  void SetPoolMultiplier(const unsigned poolMultiplier);

  /// \param verbosity Regulates the amount of information printed to standard
  ///                  output.
  ///
  /// If set to zero nothing is printed, but results are stored and can be
  /// retrieved using the GetResults() or PopResults() methods.
  /// \sa GetResults()
  /// \sa PopResults()
  void SetVerbosity(const int verbosity) { fVerbosity = verbosity; }

  /// \param Amount of iterations to run the benchmarks.
  void SetRepetitions(const unsigned repetitions) {
    fRepetitions = repetitions;
  }

  /// \param World volume containing daughters to be benchmarked.
  void SetWorld(VPlacedVolume const *const world);

  /// \return List of results of previously performed benchmarks.
  std::list<BenchmarkResult> const& GetResults() const { return fResults; }

  /// \return List of results of previously performed benchmarks. Clears the
  ///         internal history.
  std::list<BenchmarkResult> PopResults();

private:
    
  void GenerateVolumePointers(VPlacedVolume const *const vol);

  BenchmarkResult GenerateBenchmarkResult(const double elapsed,
                                          const EBenchmarkedMethod method,
                                          const EBenchmarkedLibrary library,
                                          const double bias) const;

  void RunInsideSpecialized(bool *const inside);
  void RunToInSpecialized(Precision *const distances,
                          Precision *const safeties);
  void RunToOutSpecialized(Precision *const distances,
                           Precision *const safeties);

  void RunInsideVectorized(bool *const inside);
  void RunToInVectorized(Precision *const distances, Precision *const safeties);
  void RunToOutVectorized(Precision *const distances,
                          Precision *const safeties);

  void RunInsideUnspecialized(bool *const inside);
  void RunToInUnspecialized(Precision *const distances,
                            Precision *const safeties);
  void RunToOutUnspecialized(Precision *const distances,
                             Precision *const safeties);

#ifdef VECGEOM_USOLIDS
  void RunInsideUSolids(bool *const inside);
  void RunToInUSolids(double *const distances, Precision *const safeties);
  void RunToOutUSolids(double *const distances, Precision *const safeties);
#endif
#ifdef VECGEOM_ROOT
  void RunInsideRoot(bool *const inside);
  void RunToInRoot(double *const distances, Precision *const safeties);
  void RunToOutRoot(double *const distances, Precision *const safeties);
#endif
#ifdef VECGEOM_CUDA
  void RunInsideCuda(
    Precision *const pos_x, Precision *const pos_y,
    Precision *const pos_z, bool *const inside);
  void RunToInCuda(
    Precision *const pos_x, Precision *const pos_y,
    Precision *const pos_z, Precision *const dir_x,
    Precision *const dir_y, Precision *const dir_z,
    Precision *const distances, Precision *const safeties);
  void RunToOutCuda(
    Precision *const pos_x, Precision *const pos_y,
    Precision *const pos_z, Precision *const dir_x,
    Precision *const dir_y, Precision *const dir_z,
    Precision *const distances, Precision *const safeties);
#endif

  template <typename Type>
  Type* AllocateAligned() const;

  template <typename Type>
  static void FreeAligned(Type *const distance);

  void CompareDistances(
    Precision const *const specialized,
    Precision const *const vectorized,
    Precision const *const unspecialized,
#ifdef VECGEOM_ROOT
    Precision const *const root,
#endif
#ifdef VECGEOM_USOLIDS
    Precision const *const usolids,
#endif
#ifdef VECGEOM_CUDA
    Precision const *const cuda,
#endif
    char const *const method) const;
  
};

} // End namespace vecgeom

#endif // VECGEOM_BENCHMARKING_BENCHMARKER_H_