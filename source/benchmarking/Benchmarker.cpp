/// \file Benchmarker.cpp
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#include "benchmarking/Benchmarker.h"

#include "base/SOA3D.h"
#include "base/Stopwatch.h"
#include "base/Transformation3D.h"
#include "volumes/LogicalVolume.h"
#include "volumes/PlacedBox.h"
#include "volumes/utilities/VolumeUtilities.h"

#ifdef VECGEOM_USOLIDS
#include "VUSolid.hh"
#include "UUtils.hh"
#include "UVector3.hh"
#endif

#ifdef VECGEOM_ROOT
#include "TGeoShape.h"
#endif

#include <cassert>
#include <random>
#include <sstream>

namespace vecgeom {

Benchmarker::Benchmarker() : Benchmarker(NULL) {}

Benchmarker::Benchmarker(
    VPlacedVolume const *const world)
    : fPointCount(128), fPoolMultiplier(8), fRepetitions(1<<13),
      fVerbosity(1), fToInBias(0.8), fInsideBias(0.5), fPointPool(NULL),
      fDirectionPool(NULL), fStepMax(NULL) {
  SetWorld(world);
}

Benchmarker::~Benchmarker() {
  if (fPointPool) delete fPointPool;
  if (fDirectionPool) delete fDirectionPool;
  FreeAligned(fStepMax);
}

void Benchmarker::SetWorld(VPlacedVolume const *const world) {
  fVolumes.clear();
  fWorld = world;
  if (!world) return;
  GenerateVolumePointers(fWorld);
  if (fVerbosity > 2) {
    printf("Found %lu volumes in world volume to be used for benchmarking.\n",
           fVolumes.size());
  }
}

void Benchmarker::SetPoolMultiplier(
    const unsigned poolMultiplier) {
  assert(poolMultiplier >= 1
         && "Pool multiplier for benchmarker must be >= 1.");
  fPoolMultiplier = poolMultiplier;
}

std::list<BenchmarkResult> Benchmarker::PopResults() {
  std::list<BenchmarkResult> results = fResults;
  fResults.clear();
  return results;
}

void Benchmarker::GenerateVolumePointers(VPlacedVolume const *const vol) {
  for (auto i = vol->daughters().begin(), iEnd = vol->daughters().end();
       i != iEnd; ++i) {
    fVolumes.push_back(*i);
    GenerateVolumePointers(*i);
  }
}

BenchmarkResult Benchmarker::GenerateBenchmarkResult(
    Precision elapsed, const EBenchmarkedMethod method,
    const EBenchmarkedLibrary library, const double bias) const {
  const BenchmarkResult benchmark = {
    .elapsed = elapsed,
    .method = method,
    .library = library,
    .repetitions = fRepetitions,
    .volumes = static_cast<unsigned>(fVolumes.size()),
    .points = fPointCount,
    .bias = bias
  };
  return benchmark;
}

void Benchmarker::CompareDistances(
    SOA3D<Precision> *points,
    SOA3D<Precision> *directions,
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
    char const *const method) const {

  static char const *const outputLabels =
      "Specialized / Vectorized / Unspecialized"
#ifdef VECGEOM_ROOT
      " / ROOT"
#endif
#ifdef VECGEOM_USOLIDS
      " / USolids"
#endif
#ifdef VECGEOM_CUDA
      " / CUDA"
#endif
      ;

  if (fPoolMultiplier == 1 && fVerbosity > 0) {

    printf("Comparing %s results...\n", method);
    if (fVerbosity > 2) printf("%s\n", outputLabels);

    // Compare results
    int mismatches = 0;
    for (unsigned i = 0; i < fPointCount; ++i) {
      bool mismatch = false;
      std::stringstream mismatchOutput;
      if (fVerbosity > 2) {
        mismatchOutput << specialized[i] << " / "
                       << vectorized[i] << " / "
                       << unspecialized[i];
      }
      if (std::fabs(specialized[i] - vectorized[i]) > kTolerance
          && !(specialized[i] == kInfinity && vectorized[i] == kInfinity)) {
        mismatch = true;
      }
      if (std::fabs(specialized[i] - unspecialized[i]) > kTolerance
          && !(specialized[i] == kInfinity && unspecialized[i] == kInfinity)) {
        mismatch = true;
      }
#ifdef VECGEOM_ROOT
      if (std::fabs(specialized[i] - root[i]) > kTolerance
          && !(specialized[i] == kInfinity && root[i] == 1e30)) {
        mismatch = true;
      }
      if (fVerbosity > 2) mismatchOutput << " / " << root[i];
#endif
#ifdef VECGEOM_USOLIDS
      if (std::fabs(specialized[i] - usolids[i]) > kTolerance
          && !(specialized[i] == kInfinity
               && usolids[i] == UUtils::kInfinity)) {
        mismatch = true;
      }
      if (fVerbosity > 2) mismatchOutput << " / " << usolids[i];
#endif
#ifdef VECGEOM_CUDA
      if (std::fabs(specialized[i] - cuda[i]) > kTolerance
          && !(specialized[i] == kInfinity && cuda[i] == kInfinity)) {
        mismatch = true;
      }
      if (fVerbosity > 2) mismatchOutput << " / " << cuda[i];
#endif
      mismatches += mismatch;

      if ((mismatch && fVerbosity > 2) || fVerbosity > 4) {
        printf("Point (%f, %f, %f)", points->x(i), points->y(i),
               points->z(i));
        if (directions != NULL) {
          printf(", Direction (%f, %f, %f)", directions->x(i), directions->y(i),
                 directions->z(i));
        }
        printf(": ");
      }

      if ((mismatch && fVerbosity > 2) || fVerbosity > 3) {
        printf("%s\n", mismatchOutput.str().c_str());
      }
    }
    printf("%i / %i mismatches detected.\n", mismatches, fPointCount);

  }

}

void Benchmarker::RunBenchmark() {
  Assert(fWorld, "No world specified to benchmark.\n");
  RunInsideBenchmark();
  RunToInBenchmark();
  RunToOutBenchmark();
}

void Benchmarker::RunInsideBenchmark() {

  assert(fWorld);

  if (fVerbosity > 0) {
    printf("Running Contains and Inside benchmark for %i points for "
           "%i repetitions.\n", fPointCount, fRepetitions);
  }
#ifndef VECGEOM_SCALAR
  if (fVerbosity > 1) {
    printf("Vector instruction size is %i doubles.\n", kVectorSize);
  }
#endif

  if (fPointPool) delete fPointPool;
  fPointPool = new SOA3D<Precision>(fPointCount*fPoolMultiplier);

  if (fVerbosity > 1) printf("Generating points with bias %f... ", fInsideBias);

  volumeUtilities::FillContainedPoints(*fWorld, fInsideBias, *fPointPool, true);

  if (fVerbosity > 1) printf("Done.\n");

  std::stringstream outputLabelsContains, outputLabelsInside;
  outputLabelsContains << "Specialized - Vectorized - Unspecialized";
  outputLabelsInside   << "Specialized - Vectorized - Unspecialized";

  // Allocate memory
  bool *const containsSpecialized = AllocateAligned<bool>();
  bool *const containsVectorized = AllocateAligned<bool>();
  bool *const containsUnspecialized = AllocateAligned<bool>();
  Inside_t *const insideSpecialized = AllocateAligned<Inside_t>();
  Inside_t *const insideVectorized = AllocateAligned<Inside_t>();
  Inside_t *const insideUnspecialized = AllocateAligned<Inside_t>();
#ifdef VECGEOM_ROOT
  bool *const containsRoot = AllocateAligned<bool>();
  outputLabelsContains << " - ROOT";
#endif
#ifdef VECGEOM_USOLIDS
  ::VUSolid::EnumInside *const insideUSolids =
      AllocateAligned< ::VUSolid::EnumInside>();
  outputLabelsInside << " - USolids";
#endif
#ifdef VECGEOM_CUDA
  bool *const containsCuda = AllocateAligned<bool>();
  Inside_t *const insideCuda = AllocateAligned<Inside_t>();
  outputLabelsContains << " - CUDA";
  outputLabelsInside << " - CUDA";
#endif

  // Run all benchmarks
  RunInsideSpecialized(containsSpecialized, insideSpecialized);
  RunInsideVectorized(containsVectorized, insideVectorized);
  RunInsideUnspecialized(containsUnspecialized, insideUnspecialized);
#ifdef VECGEOM_USOLIDS
  RunInsideUSolids(insideUSolids);
#endif
#ifdef VECGEOM_ROOT
  RunInsideRoot(containsRoot);
#endif
#ifdef VECGEOM_CUDA
  RunInsideCuda(fPointPool->x(), fPointPool->y(), fPointPool->z(),
                containsCuda, insideCuda);
#endif

  if (fPoolMultiplier == 1 && fVerbosity > 0) {

    // Compare Contains results

    printf("Comparing Contains results:\n");
    if (fVerbosity > 2) printf("%s\n", outputLabelsContains.str().c_str());

    // Compare results
    int mismatches = 0;
    for (unsigned i = 0; i < fPointCount; ++i) {
      bool mismatch = false;
      std::stringstream mismatchOutput;
      if (fVerbosity > 2) {
        mismatchOutput << containsSpecialized[i] << " / "
                       << containsVectorized[i] << " / "
                       << containsUnspecialized[i];
      }
      if (containsSpecialized[i] != containsVectorized[i]) mismatch = true;
      if (containsSpecialized[i] != containsUnspecialized[i]) mismatch = true;
#ifdef VECGEOM_ROOT
      if (containsSpecialized[i] != containsRoot[i]) mismatch = true;
      if (fVerbosity > 2) mismatchOutput << " / " << containsRoot[i];
#endif
#ifdef VECGEOM_CUDA
      if (insideSpecialized[i] != containsCuda[i]) mismatch = true;
      if (fVerbosity > 2) mismatchOutput << " / " << containsCuda[i];
#endif
      mismatches += mismatch;
      if ((mismatch && fVerbosity > 2) || fVerbosity > 3) {
        printf("%s\n", mismatchOutput.str().c_str());
      }
    }
    if (fVerbosity > 2 && mismatches > 100) {
      printf("%s\n", outputLabelsContains.str().c_str());
    }
    printf("%i / %i mismatches detected.\n", mismatches, fPointCount);


    // Compare Inside results

    printf("Comparing Inside results:\n");
    if (fVerbosity > 2) printf("%s\n", outputLabelsInside.str().c_str());

    mismatches = 0;
    for (unsigned i = 0; i < fPointCount; ++i) {
      bool mismatch = false;
      std::stringstream mismatchOutput;
      if (fVerbosity > 2) {
        mismatchOutput << insideSpecialized[i] << " / "
                       << insideVectorized[i] << " / "
                       << insideUnspecialized[i];
      }
      if (insideSpecialized[i] != insideVectorized[i]) mismatch = true;
      if (insideSpecialized[i] != insideUnspecialized[i]) mismatch = true;
#ifdef VECGEOM_USOLIDS
      if (insideSpecialized[i] != insideUSolids[i]) mismatch = true;
      if (fVerbosity > 2) mismatchOutput << " / " << insideUSolids[i];
#endif
#ifdef VECGEOM_CUDA
      if (insideSpecialized[i] != insideCuda[i]) mismatch = true;
      if (fVerbosity > 2) mismatchOutput << " / " << insideCuda[i];
#endif
      mismatches += mismatch;
      if ((mismatch && fVerbosity > 2) || fVerbosity > 3) {
        printf("%s\n", mismatchOutput.str().c_str());
      }
    }
    if (fVerbosity > 2 && mismatches > 100) {
      printf("%s\n", outputLabelsInside.str().c_str());
    }
    printf("%i / %i mismatches detected.\n", mismatches, fPointCount);

  }

  // Clean up memory
  FreeAligned(containsSpecialized);
  FreeAligned(containsVectorized);
  FreeAligned(containsUnspecialized);
  FreeAligned(insideSpecialized);
  FreeAligned(insideVectorized);
  FreeAligned(insideUnspecialized);
#ifdef VECGEOM_USOLIDS
  FreeAligned(insideUSolids);
#endif
#ifdef VECGEOM_ROOT
  FreeAligned(containsRoot);
#endif
#ifdef VECGEOM_CUDA
  FreeAligned(containsCuda);
  FreeAligned(insideCuda);
#endif

}

void Benchmarker::RunToInBenchmark() {

  assert(fWorld);

  if (fVerbosity > 0) {
    printf("Running DistanceToIn and SafetyToIn benchmark for %i points for "
           "%i repetitions.\n", fPointCount, fRepetitions);
  }
  if (fVerbosity > 1) {
#ifndef VECGEOM_SCALAR
    printf("Vector instruction size is %i doubles.\n", kVectorSize);
#endif
  }

  // Allocate memory
  if (fPointPool) delete fPointPool;
  if (fDirectionPool) delete fDirectionPool;
  if (fStepMax) delete fStepMax;
  fPointPool = new SOA3D<Precision>(fPointCount*fPoolMultiplier);
  fDirectionPool = new SOA3D<Precision>(fPointCount*fPoolMultiplier);
  fStepMax = AllocateAligned<Precision>();
  for (unsigned i = 0; i < fPointCount; ++i) fStepMax[i] = kInfinity;
  
  if (fVerbosity > 1) printf("Generating points with bias %f...", fToInBias);

  // Generate points not contained in any daughters and set the fraction hitting
  // a daughter to the specified bias.
  volumeUtilities::FillUncontainedPoints(*fWorld, *fPointPool);
  volumeUtilities::FillBiasedDirections(*fWorld, *fPointPool, fToInBias,
                                        *fDirectionPool);

  if (fVerbosity > 1) printf(" Done.\n");

  fPointPool->resize(fPointCount*fPoolMultiplier);
  fDirectionPool->resize(fPointCount*fPoolMultiplier);

  std::stringstream outputLabels;
  outputLabels << "Specialized - Vectorized - Unspecialized";

  // Allocate output memory
  Precision *const distancesSpecialized = AllocateAligned<Precision>();
  Precision *const safetiesSpecialized = AllocateAligned<Precision>();
  Precision *const distancesVectorized = AllocateAligned<Precision>();
  Precision *const safetiesVectorized = AllocateAligned<Precision>();
  Precision *const distancesUnspecialized = AllocateAligned<Precision>();
  Precision *const safetiesUnspecialized = AllocateAligned<Precision>();
#ifdef VECGEOM_USOLIDS
  Precision *const distancesUSolids = AllocateAligned<Precision>();
  Precision *const safetiesUSolids = AllocateAligned<Precision>();
  outputLabels << " - USolids";
#endif
#ifdef VECGEOM_ROOT
  Precision *const distancesRoot = AllocateAligned<Precision>();
  Precision *const safetiesRoot = AllocateAligned<Precision>();
  outputLabels << " - ROOT";
#endif
#ifdef VECGEOM_CUDA
  Precision *const distancesCuda = AllocateAligned<Precision>();
  Precision *const safetiesCuda = AllocateAligned<Precision>();
  outputLabels << " - CUDA";
#endif

  // Run all benchmarks
  RunToInSpecialized(distancesSpecialized, safetiesSpecialized);
  RunToInVectorized(distancesVectorized, safetiesVectorized);
  RunToInUnspecialized(distancesUnspecialized, safetiesUnspecialized);
#ifdef VECGEOM_USOLIDS
  RunToInUSolids(distancesUSolids, safetiesUSolids);
#endif
#ifdef VECGEOM_ROOT
  RunToInRoot(distancesRoot, safetiesRoot);
#endif
#ifdef VECGEOM_CUDA
  RunToInCuda(fPointPool->x(),     fPointPool->y(),     fPointPool->z(),
              fDirectionPool->x(), fDirectionPool->y(), fDirectionPool->z(),
              distancesCuda, safetiesCuda);
#endif

  CompareDistances(
    fPointPool,
    fDirectionPool,
    distancesSpecialized,
    distancesVectorized,
    distancesUnspecialized,
#ifdef VECGEOM_ROOT
    distancesRoot,
#endif
#ifdef VECGEOM_USOLIDS
    distancesUSolids,
#endif
#ifdef VECGEOM_CUDA
    distancesCuda,
#endif
    "DistanceToIn");

  // Clean up memory
  FreeAligned(distancesSpecialized);
  FreeAligned(distancesUnspecialized);
  FreeAligned(distancesVectorized);
#ifdef VECGEOM_USOLIDS
  FreeAligned(distancesUSolids);
#endif
#ifdef VECGEOM_ROOT
  FreeAligned(distancesRoot);
#endif
#ifdef VECGEOM_CUDA
  FreeAligned(distancesCuda);
#endif

  CompareDistances(
    fPointPool,
    NULL,
    safetiesSpecialized,
    safetiesVectorized,
    safetiesUnspecialized,
#ifdef VECGEOM_ROOT
    safetiesRoot,
#endif
#ifdef VECGEOM_USOLIDS
    safetiesUSolids,
#endif
#ifdef VECGEOM_CUDA
    safetiesCuda,
#endif
    "SafetyToIn");

  FreeAligned(safetiesSpecialized);
  FreeAligned(safetiesUnspecialized);
  FreeAligned(safetiesVectorized);
#ifdef VECGEOM_USOLIDS
  FreeAligned(safetiesUSolids);
#endif
#ifdef VECGEOM_ROOT
  FreeAligned(safetiesRoot);
#endif
#ifdef VECGEOM_CUDA
  FreeAligned(safetiesCuda);
#endif

}

void Benchmarker::RunToOutBenchmark() {

  assert(fWorld);
  
  if (fVerbosity > 0) {
    printf("Running DistanceToOut and SafetyToOut benchmark for %i points for "
           "%i repetitions.\n", fPointCount, fRepetitions);
  }
  if (fVerbosity > 1) {
#ifndef VECGEOM_SCALAR
    printf("Vector instruction size is %i doubles.\n", kVectorSize);
#endif
  }

  // Allocate memory
  if (fPointPool) delete fPointPool;
  if (fDirectionPool) delete fDirectionPool;
  if (fStepMax) delete fStepMax;
  fPointPool = new SOA3D<Precision>(fPointCount*fPoolMultiplier);
  fDirectionPool = new SOA3D<Precision>(fPointCount*fPoolMultiplier);
  fStepMax = AllocateAligned<Precision>();
  for (unsigned i = 0; i < fPointCount; ++i) fStepMax[i] = kInfinity;
  
  if (fVerbosity > 1) printf("Generating points...");

  // Generate points not contained in any daughters and set the fraction hitting
  // a daughter to the specified bias.
  volumeUtilities::FillContainedPoints(*fWorld, *fPointPool, false);
  volumeUtilities::FillRandomDirections(*fDirectionPool);

  if (fVerbosity > 1) printf(" Done.\n");

  fPointPool->resize(fPointCount*fPoolMultiplier);
  fDirectionPool->resize(fPointCount*fPoolMultiplier);

  std::stringstream outputLabels;
  outputLabels << "Specialized - Vectorized - Unspecialized";

  // Allocate output memory
  Precision *const distancesSpecialized = AllocateAligned<Precision>();
  Precision *const safetiesSpecialized = AllocateAligned<Precision>();
  Precision *const distancesVectorized = AllocateAligned<Precision>();
  Precision *const safetiesVectorized = AllocateAligned<Precision>();
  Precision *const distancesUnspecialized = AllocateAligned<Precision>();
  Precision *const safetiesUnspecialized = AllocateAligned<Precision>();
#ifdef VECGEOM_USOLIDS
  Precision *const distancesUSolids = AllocateAligned<Precision>();
  Precision *const safetiesUSolids = AllocateAligned<Precision>();
  outputLabels << " - USolids";
#endif
#ifdef VECGEOM_ROOT
  Precision *const distancesRoot = AllocateAligned<Precision>();
  Precision *const safetiesRoot = AllocateAligned<Precision>();
  outputLabels << " - ROOT";
#endif
#ifdef VECGEOM_CUDA
  Precision *const distancesCuda = AllocateAligned<Precision>();
  Precision *const safetiesCuda = AllocateAligned<Precision>();
  outputLabels << " - CUDA";
#endif

  // Run all benchmarks
  RunToOutSpecialized(distancesSpecialized, safetiesSpecialized);
  RunToOutVectorized(distancesVectorized, safetiesVectorized);
  RunToOutUnspecialized(distancesUnspecialized, safetiesUnspecialized);
#ifdef VECGEOM_USOLIDS
  RunToOutUSolids(distancesUSolids, safetiesUSolids);
#endif
#ifdef VECGEOM_ROOT
  RunToOutRoot(distancesRoot, safetiesRoot);
#endif
#ifdef VECGEOM_CUDA
  RunToOutCuda(fPointPool->x(),     fPointPool->y(),     fPointPool->z(),
               fDirectionPool->x(), fDirectionPool->y(), fDirectionPool->z(),
               distancesCuda, safetiesCuda);
#endif

  CompareDistances(
    fPointPool,
    fDirectionPool,
    distancesSpecialized,
    distancesVectorized,
    distancesUnspecialized,
#ifdef VECGEOM_ROOT
    distancesRoot,
#endif
#ifdef VECGEOM_USOLIDS
    distancesUSolids,
#endif
#ifdef VECGEOM_CUDA
    distancesCuda,
#endif
    "DistanceToOut");

  // Clean up memory
  FreeAligned(distancesSpecialized);
  FreeAligned(distancesUnspecialized);
  FreeAligned(distancesVectorized);
#ifdef VECGEOM_USOLIDS
  FreeAligned(distancesUSolids);
#endif
#ifdef VECGEOM_ROOT
  FreeAligned(distancesRoot);
#endif
#ifdef VECGEOM_CUDA
  FreeAligned(distancesCuda);
#endif

  CompareDistances(
    fPointPool,
    NULL,
    safetiesSpecialized,
    safetiesVectorized,
    safetiesUnspecialized,
#ifdef VECGEOM_ROOT
    safetiesRoot,
#endif
#ifdef VECGEOM_USOLIDS
    safetiesUSolids,
#endif
#ifdef VECGEOM_CUDA
    safetiesCuda,
#endif
    "SafetyToOut");

  FreeAligned(safetiesSpecialized);
  FreeAligned(safetiesUnspecialized);
  FreeAligned(safetiesVectorized);
#ifdef VECGEOM_USOLIDS
  FreeAligned(safetiesUSolids);
#endif
#ifdef VECGEOM_ROOT
  FreeAligned(safetiesRoot);
#endif
#ifdef VECGEOM_CUDA
  FreeAligned(safetiesCuda);
#endif

}

void Benchmarker::RunInsideSpecialized(bool *contains, Inside_t *inside) {
  if (fVerbosity > 0) printf("Specialized   - ");
  Stopwatch timer;
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      for (unsigned i = 0; i < fPointCount; ++i) {
        contains[i] = v->specialized()->Contains((*fPointPool)[index + i]);
      }
    }
  }
  Precision elapsedContains = timer.Stop();
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      for (unsigned i = 0; i < fPointCount; ++i) {
        inside[i] = v->specialized()->Inside((*fPointPool)[index + i]);
      }
    }
  }
  Precision elapsedInside = timer.Stop();
  if (fVerbosity > 0) {
    printf("Inside: %.6fs (%.6fs), Contains: %.6fs (%.6fs), "
           "Inside/Contains: %.2f\n",
           elapsedInside, elapsedInside/fVolumes.size(),
           elapsedContains, elapsedContains/fVolumes.size(),
           elapsedInside/elapsedContains);
  }
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedContains, kBenchmarkContains, kBenchmarkSpecialized, fInsideBias
    )
  );
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedInside, kBenchmarkInside, kBenchmarkSpecialized, fInsideBias
    )
  );
}

void Benchmarker::RunToInSpecialized(
    Precision *distances, Precision *safeties) {
  if (fVerbosity > 0) printf("Specialized   - ");
  Stopwatch timer;
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      for (unsigned i = 0; i < fPointCount; ++i) {
        int p = index + i;
        distances[i] = v->specialized()->DistanceToIn(
          (*fPointPool)[p], (*fDirectionPool)[p]
        );
      }
    }
  }
  Precision elapsedDistance = timer.Stop();
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      for (unsigned i = 0; i < fPointCount; ++i) {
        int p = index + i;
        safeties[i] = v->specialized()->SafetyToIn(
          (*fPointPool)[p]
        );
      }
    }
  }
  Precision elapsedSafety = timer.Stop();
  if (fVerbosity > 0) {
    printf("DistanceToIn: %.6fs (%.6fs), SafetyToIn: %.6fs (%.6fs), "
           "DistanceToIn/SafetyToIn: %.2f\n",
           elapsedDistance, elapsedDistance/fVolumes.size(),
           elapsedSafety, elapsedSafety/fVolumes.size(),
           elapsedDistance/elapsedSafety);
  }
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedDistance, kBenchmarkDistanceToIn, kBenchmarkSpecialized, fToInBias
    )
  );
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedSafety, kBenchmarkSafetyToIn, kBenchmarkSpecialized, fToInBias
    )
  );
}

void Benchmarker::RunToOutSpecialized(
    Precision *distances, Precision *safeties) {
  if (fVerbosity > 0) printf("Specialized   - ");
  Stopwatch timer;
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      for (unsigned i = 0; i < fPointCount; ++i) {
        int p = index + i;
        distances[i] = v->specialized()->DistanceToOut(
          (*fPointPool)[p], (*fDirectionPool)[p]
        );
      }
    }
  }
  Precision elapsedDistance = timer.Stop();
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      for (unsigned i = 0; i < fPointCount; ++i) {
        int p = index + i;
        safeties[i] = v->specialized()->SafetyToOut(
          (*fPointPool)[p]
        );
      }
    }
  }
  Precision elapsedSafety = timer.Stop();
  if (fVerbosity > 0) {
    printf("DistanceToOut: %.6fs (%.6fs), SafetyToOut: %.6fs (%.6fs), "
           "DistanceToOut/SafetyToOut: %.2f\n",
           elapsedDistance, elapsedDistance/fVolumes.size(),
           elapsedSafety, elapsedSafety/fVolumes.size(),
           elapsedDistance/elapsedSafety);
  }
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedDistance, kBenchmarkDistanceToOut, kBenchmarkSpecialized, 1
    )
  );
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedSafety, kBenchmarkSafetyToOut, kBenchmarkSpecialized, fToInBias
    )
  );
}

void Benchmarker::RunInsideVectorized(bool *contains, Inside_t *inside) {
  if (fVerbosity > 0) {
    printf("Vectorized    - ");
  }
  Stopwatch timer;
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    SOA3D<Precision> points(fPointPool->x()+index, fPointPool->y()+index,
                            fPointPool->z()+index, fPointCount);
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      v->specialized()->Contains(points, contains);
    }
  }
  Precision elapsedContains = timer.Stop();
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    SOA3D<Precision> points(fPointPool->x()+index, fPointPool->y()+index,
                            fPointPool->z()+index, fPointCount);
    for (auto v = fVolumes.begin(), v_end = fVolumes.end(); v != v_end; ++v) {
      v->specialized()->Inside(points, inside);
    }
  }
  Precision elapsedInside = timer.Stop();
  if (fVerbosity > 0) {
    printf("Inside: %.6fs (%.6fs), Contains: %.6fs (%.6fs), "
           "Inside/Contains: %.2f\n",
           elapsedInside, elapsedInside/fVolumes.size(),
           elapsedContains, elapsedContains/fVolumes.size(),
           elapsedInside/elapsedContains);
  }
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedContains, kBenchmarkContains, kBenchmarkVectorized, fInsideBias
    )
  );
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedInside, kBenchmarkInside, kBenchmarkVectorized, fInsideBias
    )
  );
}

void Benchmarker::RunToInVectorized(
    Precision *const distances, Precision *const safeties) {
  if (fVerbosity > 0) {
    printf("Vectorized    - ");
  }
  Stopwatch timer;
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    SOA3D<Precision> points(fPointPool->x()+index, fPointPool->y()+index,
                            fPointPool->z()+index, fPointCount);
    SOA3D<Precision> directions(fDirectionPool->x() + index,
                                fDirectionPool->y() + index,
                                fDirectionPool->z() + index, fPointCount);
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      v->specialized()->DistanceToIn(points, directions, fStepMax, distances);
    }
  }
  Precision elapsedDistance = timer.Stop();
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    SOA3D<Precision> points(fPointPool->x()+index, fPointPool->y()+index,
                            fPointPool->z()+index, fPointCount);
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      v->specialized()->SafetyToIn(points, safeties);
    }
  }
  Precision elapsedSafety = timer.Stop();
  if (fVerbosity > 0) {
    printf("DistanceToIn: %.6fs (%.6fs), SafetyToIn: %.6fs (%.6fs), "
           "DistanceToIn/SafetyToIn: %.2f\n",
           elapsedDistance, elapsedDistance/fVolumes.size(),
           elapsedSafety, elapsedSafety/fVolumes.size(),
           elapsedDistance/elapsedSafety);
  }
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedDistance, kBenchmarkDistanceToIn, kBenchmarkVectorized, fToInBias
    )
  );
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedSafety, kBenchmarkSafetyToIn, kBenchmarkVectorized, fToInBias
    )
  );
}

void Benchmarker::RunToOutVectorized(
    Precision *const distances, Precision *const safeties) {
  if (fVerbosity > 0) {
    printf("Vectorized    - ");
  }
  Stopwatch timer;
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    SOA3D<Precision> points(fPointPool->x()+index, fPointPool->y()+index,
                            fPointPool->z()+index, fPointCount);
    SOA3D<Precision> directions(fDirectionPool->x()+index,
                                fDirectionPool->y()+index,
                                fDirectionPool->z()+index, fPointCount);
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      v->specialized()->DistanceToOut(points, directions, fStepMax, distances);
    }
  }
  Precision elapsedDistance = timer.Stop();
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    SOA3D<Precision> points(fPointPool->x()+index, fPointPool->y()+index,
                            fPointPool->z()+index, fPointCount);
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      v->specialized()->SafetyToOut(points, safeties);
    }
  }
  Precision elapsedSafety = timer.Stop();
  if (fVerbosity > 0) {
    printf("DistanceToOut: %.6fs (%.6fs), SafetyToOut: %.6fs (%.6fs), "
           "DistanceToOut/SafetyToOut: %.2f\n",
           elapsedDistance, elapsedDistance/fVolumes.size(),
           elapsedSafety, elapsedSafety/fVolumes.size(),
           elapsedDistance/elapsedSafety);
  }
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedDistance, kBenchmarkDistanceToOut, kBenchmarkVectorized, 1
    )
  );
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedSafety, kBenchmarkSafetyToOut, kBenchmarkVectorized, 1
    )
  );
}

void Benchmarker::RunInsideUnspecialized(bool *contains, Inside_t *inside) {
  if (fVerbosity > 0) printf("Unspecialized - ");
  Stopwatch timer;
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      for (unsigned i = 0; i < fPointCount; ++i) {
        int p = index + i;
        contains[i] = v->unspecialized()->Contains((*fPointPool)[p]);
      }
    }
  }
  Precision elapsedContains = timer.Stop();
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      for (unsigned i = 0; i < fPointCount; ++i) {
        int p = index + i;
        inside[i] = v->unspecialized()->Inside((*fPointPool)[p]);
      }
    }
  }
  Precision elapsedInside = timer.Stop();
  if (fVerbosity > 0) {
    printf("Inside: %.6fs (%.6fs), Contains: %.6fs (%.6fs), "
           "Inside/Contains: %.2f\n",
           elapsedInside, elapsedInside/fVolumes.size(),
           elapsedContains, elapsedContains/fVolumes.size(),
           elapsedInside/elapsedContains);
  }
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedContains, kBenchmarkContains, kBenchmarkUnspecialized, fInsideBias
    )
  );
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedInside, kBenchmarkInside, kBenchmarkUnspecialized, fInsideBias
    )
  );
}

void Benchmarker::RunToInUnspecialized(
    Precision *const distances, Precision *const safeties) {
  if (fVerbosity > 0) printf("Unspecialized - ");
  Stopwatch timer;
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      for (unsigned i = 0; i < fPointCount; ++i) {
        int p = index + i;
        distances[i] = v->unspecialized()->DistanceToIn(
          (*fPointPool)[p], (*fDirectionPool)[p]
        );
      }
    }
  }
  Precision elapsedDistance = timer.Stop();
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      for (unsigned i = 0; i < fPointCount; ++i) {
        int p = index + i;
        safeties[i] = v->unspecialized()->SafetyToIn((*fPointPool)[p]);
      }
    }
  }
  Precision elapsedSafety = timer.Stop();
  if (fVerbosity > 0) {
    printf("DistanceToIn: %.6fs (%.6fs), SafetyToIn: %.6fs (%.6fs), "
           "DistanceToIn/SafetyToIn: %.2f\n",
           elapsedDistance, elapsedDistance/fVolumes.size(),
           elapsedSafety, elapsedSafety/fVolumes.size(),
           elapsedDistance/elapsedSafety);
  }
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedDistance, kBenchmarkDistanceToIn, kBenchmarkUnspecialized,
      fToInBias
    )
  );
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedSafety, kBenchmarkSafetyToIn, kBenchmarkUnspecialized, fToInBias
    )
  );
}

void Benchmarker::RunToOutUnspecialized(
    Precision *const distances, Precision *const safeties) {
  if (fVerbosity > 0) printf("Unspecialized - ");
  Stopwatch timer;
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      for (unsigned i = 0; i < fPointCount; ++i) {
        int p = index + i;
        distances[i] = v->unspecialized()->DistanceToOut(
          (*fPointPool)[p], (*fDirectionPool)[p]
        );
      }
    }
  }
  Precision elapsedDistance = timer.Stop();
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      for (unsigned i = 0; i < fPointCount; ++i) {
        int p = index + i;
        safeties[i] = v->unspecialized()->SafetyToOut((*fPointPool)[p]);
      }
    }
  }
  Precision elapsedSafety = timer.Stop();
  if (fVerbosity > 0) {
    printf("DistanceToOut: %.6fs (%.6fs), SafetyToOut: %.6fs (%.6fs), "
           "DistanceToOut/SafetyToOut: %.2f\n",
           elapsedDistance, elapsedDistance/fVolumes.size(),
           elapsedSafety, elapsedSafety/fVolumes.size(),
           elapsedDistance/elapsedSafety);
  }
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedDistance, kBenchmarkDistanceToOut, kBenchmarkUnspecialized, 1
    )
  );
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedSafety, kBenchmarkSafetyToOut, kBenchmarkUnspecialized, 1
    )
  );
}

#ifdef VECGEOM_USOLIDS
void Benchmarker::RunInsideUSolids(::VUSolid::EnumInside *const inside) {
  if (fVerbosity > 0) printf("USolids       - ");
  Stopwatch timer;
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      Transformation3D const *transformation =
          v->unspecialized()->transformation();
      for (unsigned i = 0; i < fPointCount; ++i) {
        int p = index + i;
        const Vector3D<Precision> point =
            transformation->Transform((*fPointPool)[p]);
        inside[i] = v->usolids()->Inside(
          UVector3(point[0], point[1], point[2])
        );
      }
    }
  }
  Precision elapsed = timer.Stop();
  if (fVerbosity > 0) {
    printf("Inside: %.6fs (%.6fs), Contains: -.------s (-.------s), "
           "Inside/Contains: -.--\n",
           elapsed, elapsed/fVolumes.size());
  }
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsed, kBenchmarkInside, kBenchmarkUSolids, fInsideBias
    )
  );
}

void Benchmarker::RunToInUSolids(
    Precision *const distances, Precision *const safeties) {
  if (fVerbosity > 0) printf("USolids       - ");
  Stopwatch timer;
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      Transformation3D const *transformation =
          v->unspecialized()->transformation();
      for (unsigned i = 0; i < fPointCount; ++i) {
        int p = index + i;
        const Vector3D<Precision> point =
            transformation->Transform((*fPointPool)[p]);
        const Vector3D<Precision> dir =
            transformation->TransformDirection((*fDirectionPool)[p]);
        distances[i] = v->usolids()->DistanceToIn(
          UVector3(point[0], point[1], point[2]),
          UVector3(dir[0], dir[1], dir[2])
        );
      }
    }
  }
  Precision elapsedDistance = timer.Stop();
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      Transformation3D const *transformation =
          v->unspecialized()->transformation();
      for (unsigned i = 0; i < fPointCount; ++i) {
        int p = index + i;
        const Vector3D<Precision> point =
            transformation->Transform((*fPointPool)[p]);
        safeties[i] = v->usolids()->SafetyFromOutside(
          UVector3(point[0], point[1], point[2])
        );
      }
    }
  }
  Precision elapsedSafety = timer.Stop();
  if (fVerbosity > 0) {
    printf("DistanceToIn: %.6fs (%.6fs), SafetyToIn: %.6fs (%.6fs), "
           "DistanceToIn/SafetyToIn: %.2f\n",
           elapsedDistance, elapsedDistance/fVolumes.size(),
           elapsedSafety, elapsedSafety/fVolumes.size(),
           elapsedDistance/elapsedSafety);
  }
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedDistance, kBenchmarkDistanceToIn, kBenchmarkUSolids, fToInBias
    )
  );
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedSafety, kBenchmarkSafetyToIn, kBenchmarkUSolids, fToInBias
    )
  );
}

void Benchmarker::RunToOutUSolids(
    Precision *const distances, Precision *const safeties) {
  if (fVerbosity > 0) printf("USolids       - ");
  Stopwatch timer;
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      for (unsigned i = 0; i < fPointCount; ++i) {
        int p = index + i;
        const Vector3D<Precision> point = (*fPointPool)[p];
        const Vector3D<Precision> dir = (*fDirectionPool)[p];
        UVector3 normal;
        bool convex;
        distances[i] = v->usolids()->DistanceToOut(
          UVector3(point[0], point[1], point[2]),
          UVector3(dir[0], dir[1], dir[2]),
          normal,
          convex
        );
      }
    }
  }
  Precision elapsedDistance = timer.Stop();
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      for (unsigned i = 0; i < fPointCount; ++i) {
        int p = index + i;
        const Vector3D<Precision> point = (*fPointPool)[p];
        safeties[i] = v->usolids()->SafetyFromInside(
          UVector3(point[0], point[1], point[2])
        );
      }
    }
  }
  Precision elapsedSafety = timer.Stop();
  if (fVerbosity > 0) {
    printf("DistanceToOut: %.6fs (%.6fs), SafetyToOut: %.6fs (%.6fs), "
           "DistanceToOut/SafetyToOut: %.2f\n",
           elapsedDistance, elapsedDistance/fVolumes.size(),
           elapsedSafety, elapsedSafety/fVolumes.size(),
           elapsedDistance/elapsedSafety);
  }
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedDistance, kBenchmarkDistanceToOut, kBenchmarkUSolids, 1
    )
  );
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedSafety, kBenchmarkSafetyToOut, kBenchmarkUSolids, 1
    )
  );
}
#endif

#ifdef VECGEOM_ROOT
void Benchmarker::RunInsideRoot(bool *inside) {
  if (fVerbosity > 0) printf("ROOT          - ");
  Stopwatch timer;
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      Transformation3D const *transformation =
          v->unspecialized()->transformation();
      for (unsigned i = 0; i < fPointCount; ++i) {
        int p = index + i;
        Vector3D<Precision> point =
            transformation->Transform((*fPointPool)[p]);
        inside[i] = v->root()->Contains(&point[0]);
      }
    }
  }
  Precision elapsed = timer.Stop();
  if (fVerbosity > 0) {
    printf("Inside: -.------s (-.------s), Contains: %.6fs (%.6fs), "
           "Inside/Contains: -.--\n",
           elapsed, elapsed/fVolumes.size());
  }
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsed, kBenchmarkContains, kBenchmarkRoot, fInsideBias
    )
  );
}

void Benchmarker::RunToInRoot(
    Precision *const distances, Precision *const safeties) {
  if (fVerbosity > 0) printf("ROOT          - ");
  Stopwatch timer;
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      Transformation3D const *transformation =
          v->unspecialized()->transformation();
      for (unsigned i = 0; i < fPointCount; ++i) {
        int p = index + i;
        Vector3D<Precision> point =
            transformation->Transform((*fPointPool)[p]);
        Vector3D<Precision> dir =
            transformation->TransformDirection((*fDirectionPool)[p]);
        distances[i] = v->root()->DistFromOutside(&point[0], &dir[0]);
      }
    }
  }
  Precision elapsedDistance = timer.Stop();
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      Transformation3D const *transformation =
          v->unspecialized()->transformation();
      for (unsigned i = 0; i < fPointCount; ++i) {
        int p = index + i;
        Vector3D<Precision> point =
            transformation->Transform((*fPointPool)[p]);
        safeties[i] = v->root()->Safety(&point[0], false);
      }
    }
  }
  Precision elapsedSafety = timer.Stop();
  if (fVerbosity > 0) {
    printf("DistanceToIn: %.6fs (%.6fs), SafetyToIn: %.6fs (%.6fs), "
           "DistanceToIn/SafetyToIn: %.2f\n",
           elapsedDistance, elapsedDistance/fVolumes.size(),
           elapsedSafety, elapsedSafety/fVolumes.size(),
           elapsedDistance/elapsedSafety);
  }
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedDistance, kBenchmarkDistanceToIn, kBenchmarkRoot, fToInBias
    )
  );
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedSafety, kBenchmarkSafetyToIn, kBenchmarkRoot, fToInBias
    )
  );
}
void Benchmarker::RunToOutRoot(
    Precision *const distances, Precision *const safeties) {
  if (fVerbosity > 0) printf("ROOT          - ");
  Stopwatch timer;
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      for (unsigned i = 0; i < fPointCount; ++i) {
        int p = index + i;
        Vector3D<Precision> point = (*fPointPool)[p];
        Vector3D<Precision> dir = (*fDirectionPool)[p];
        distances[i] = v->root()->DistFromInside(&point[0], &dir[0]);
      }
    }
  }
  Precision elapsedDistance = timer.Stop();
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    int index = (rand() % fPoolMultiplier) * fPointCount;
    for (auto v = fVolumes.begin(), vEnd = fVolumes.end(); v != vEnd; ++v) {
      for (unsigned i = 0; i < fPointCount; ++i) {
        int p = index + i;
        Vector3D<Precision> point = (*fPointPool)[p];
        safeties[i] = v->root()->Safety(&point[0], true);
      }
    }
  }
  Precision elapsedSafety = timer.Stop();
  if (fVerbosity > 0) {
    printf("DistanceToOut: %.6fs (%.6fs), SafetyToOut: %.6fs (%.6fs), "
           "DistanceToOut/SafetyToOut: %.2f\n",
           elapsedDistance, elapsedDistance/fVolumes.size(),
           elapsedSafety, elapsedSafety/fVolumes.size(),
           elapsedDistance/elapsedSafety);
  }
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedDistance, kBenchmarkDistanceToOut, kBenchmarkRoot, 1
    )
  );
  fResults.push_back(
    GenerateBenchmarkResult(
      elapsedSafety, kBenchmarkSafetyToOut, kBenchmarkRoot, 
1    )
  );
}
#endif

template <typename Type>
Type* Benchmarker::AllocateAligned() const {
  return (Type*) _mm_malloc(fPointCount*sizeof(Type), kAlignmentBoundary);
}

template <typename Type>
void Benchmarker::FreeAligned(Type *const distance) {
  if (distance) _mm_free(distance);
}

} // End namespace vecgeom
