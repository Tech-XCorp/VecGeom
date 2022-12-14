/*
 * RaytraceBenchmark.cpp
 *
 *  Created on: May 8, 2020
 *      Author: andrei.gheata@cern.ch
 */

#include <iomanip>
#include <VecGeom/management/GeoManager.h>
#include <VecGeom/management/NavIndexTable.h>
#include <VecGeom/navigation/NavigationState.h>
#include <VecGeom/base/Stopwatch.h>
#include "ArgParser.h"

#ifdef VECGEOM_GDML
#include "Frontend.h"
#endif

using namespace vecgeom;

int TestNavIndexCPU(vecgeom::cxx::VPlacedVolume const *const world, int maxdepth);

#ifdef VECGEOM_ENABLE_CUDA
int TestNavIndexGPU(vecgeom::cxx::VPlacedVolume const *const world, int maxdepth);
#endif

namespace visitorcxx {

class GlobalToLocalVisitor {
private:
  int fError = 0; ///< error code

public:
  GlobalToLocalVisitor() {}

  VECCORE_ATT_HOST_DEVICE
  int GetError() const { return fError; }

  VECCORE_ATT_HOST_DEVICE
  void apply(NavStatePath *state, NavIndex_t nav_index)
  {
    unsigned char level            = state->GetLevel();
    int dind                       = 0;
    NavIndex_t nav_ind             = 1;
    VPlacedVolume const *pdaughter = nullptr;
    for (int i = 1; i < level + 1; ++i) {
      pdaughter = state->At(i);
      dind      = pdaughter->GetChildId();
      if (dind < 0) {
        fError = 1;
        return;
      }
      nav_ind = NavStateIndex::PushImpl(nav_ind, pdaughter);
    }

    // Check if navigation index matches input
    if (nav_ind != nav_index) {
      fError = 2;
      return;
    }

    // Check if the physical volume is correct
    if (NavStateIndex::TopImpl(nav_ind) != state->Top()) {
      fError = 3;
      return;
    }

    // Check if the current level is valid
    if (level != NavStateIndex::GetLevelImpl(nav_ind)) {
      fError = 4;
      return;
    }

    // Check if mother navigation index is consistent
    if (level > 0 && nav_ind != NavStateIndex::PushImpl(NavStateIndex::PopImpl(nav_ind), pdaughter)) {
      fError = 5;
      return;
    }

    // Check if the number of daughters is correct
    if (NavStateIndex::GetNdaughtersImpl(nav_ind) != state->Top()->GetDaughters().size()) {
      fError = 6;
      return;
    }

    Transformation3D trans, trans_nav_ind;
    state->TopMatrix(trans);
    NavStateIndex::TopMatrixImpl(nav_ind, trans_nav_ind);
    if (!trans.operator==(trans_nav_ind)) {
      fError = 7;
      return;
    }

    // success
    fError = 0;
  }
};

/// Traverses the geometry tree keeping track of the state context (volume path or navigation state)
/// and applies the injected Visitor
VECCORE_ATT_HOST_DEVICE
template <typename Visitor>
int visitAllPlacedVolumesPassNavIndex(VPlacedVolume const *currentvolume, Visitor *visitor, NavStatePath *state,
                                      NavIndex_t nav_ind)
{
  const char *errcodes[] = {"incompatible daughter pointer",
                            "navigation index mismatch",
                            "top placed volume pointer mismatch",
                            "level mismatch",
                            "navigation index inconsistency for Push/Pop",
                            "number of daughters mismatch",
                            "transformation matrix mismatch"};
  if (currentvolume != NULL) {
    state->Push(currentvolume);
    visitor->apply(state, nav_ind);
    auto ierr = visitor->GetError();
    if (ierr) {
      printf("=== EEE === TestNavIndex: %s\n", errcodes[ierr - 1]);
      return ierr;
    }
    for (auto daughter : currentvolume->GetDaughters()) {
      auto nav_ind_d = NavStateIndex::PushImpl(nav_ind, daughter);
      ierr           = visitAllPlacedVolumesPassNavIndex(daughter, visitor, state, nav_ind_d);
      if (ierr) return ierr;
    }
    state->Pop();
  }
  return 0;
}

} // namespace visitorcxx

int TestNavIndexCPU(vecgeom::cxx::VPlacedVolume const *const world, int maxdepth)
{
  // Check performance
  using namespace visitorcxx;

  Stopwatch timer;
  NavStatePath *state = NavStatePath::MakeInstance(maxdepth);
  state->Clear();
  GlobalToLocalVisitor visitor;

  NavIndex_t nav_ind_top = 1; // The navigation index corresponding to the world

  timer.Start();
  auto ierr      = visitAllPlacedVolumesPassNavIndex(world, &visitor, state, nav_ind_top);
  auto tvalidate = timer.Stop();

  if (!ierr) std::cout << "=== Info navigation table validation on CPU took: " << tvalidate << " sec.\n";

  NavStatePath::ReleaseInstance(state);
  return ierr;
}

int main(int argc, char *argv[])
{
  OPTION_STRING(gdml_name, "default.gdml");
  OPTION_INT(max_depth, 0);
  OPTION_INT(on_gpu, 0);
#ifndef VECGEOM_GDML
  (void)max_depth;
  std::cout << "### VecGeom must be compiled with GDML support to run this.\n";
  return 1;
#endif

  Stopwatch timer;
  // Try to open the input file
#ifdef VECGEOM_GDML
  GeoManager::Instance().SetTransformationCacheDepth(max_depth);
  auto load = vgdml::Frontend::Load(gdml_name.c_str(), false);
  if (!load) return 2;
#endif

  auto world = GeoManager::Instance().GetWorld();
  if (!world) return 3;
  int maxdepth = GeoManager::Instance().getMaxDepth();

  auto ierr = 0;
  timer.Start();
  if (on_gpu) {
#ifdef VECGEOM_ENABLE_CUDA
    ierr = TestNavIndexGPU(GeoManager::Instance().GetWorld(), maxdepth);
#else
    std::cout << "=== Cannot run the test on GPU since VecGeom CUDA support not compiled.\n";
    return 1;
#endif
  } else {
    ierr = TestNavIndexCPU(GeoManager::Instance().GetWorld(), maxdepth);
  }
  auto validation_time = timer.Stop();
  if (ierr)
    std::cout << "TestNavIndex FAILED\n";
  else {
    std::cout << "Navigation index table validation took " << validation_time << " seconds\n";
    std::cout << "TestNavIndex PASSED\n";
  }

  return ierr;
}
