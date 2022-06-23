#ifndef VECGEOM_SURFACE_NAVIGATOR_H_
#define VECGEOM_SURFACE_NAVIGATOR_H_

#include <VecGeom/surfaces/Model.h>
#include <VecGeom/navigation/NavStateIndex.h>

namespace vgbrep {
namespace protonav {

///< Returns the distance to the next surface starting from a global point located in in_state, optionally
///< skipping the surface exit_surf. Computes the new state after crossing.
template <typename Real_t>
Real_t ComputeStepAndHit(vecgeom::Vector3D<Real_t> const &point, vecgeom::Vector3D<Real_t> const &direction,
                         vecgeom::NavStateIndex const &in_state, vecgeom::NavStateIndex &out_state,
                         SurfData<Real_t> const &surfdata, int &exit_surf)
{
  // Get the list of candidate surfaces for in_state
  out_state            = in_state;
  auto skip_surf       = exit_surf;
  exit_surf            = 0;
  Real_t distance      = vecgeom::InfinityLength<Real_t>();
  int isurfcross       = 0;
  NavIndex_t in_navind = in_state.GetNavIndex();
  auto const &cand     = surfdata.fCandidates[in_state.GetId()];
  // simple loop on all candidates. This should be in future optimized to give the reduced list that may be crossed
  for (auto icand = 0; icand < cand.fNcand; ++icand) {
    bool can_hit = true;
    int isurf    = std::abs(cand[icand]);
    if (std::abs(isurf) == std::abs(skip_surf)) continue;
    auto const &surf = surfdata.fCommonSurfaces[isurf];
    // Convert point and direction to surface frame
    vecgeom::Transformation3D const &trans = surfdata.fGlobalTrans[surf.fTrans];
    Vector3D<Real_t> local                 = trans.Transform(point);
    Vector3D<Real_t> localdir              = trans.TransformDirection(direction);
    // Compute distance to surface
    auto dist = surfdata.GetUnplaced(isurf).Intersect(local, localdir, surfdata);
    if (dist < vecgeom::kTolerance || dist > distance) continue;
    Vector3D<Real_t> onsurf = local + dist * localdir;

    // should index sides rather than left/right...
    auto const &exit_side  = (cand[icand] > 0) ? surf.fLeftSide : surf.fRightSide;
    auto const &entry_side = (cand[icand] > 0) ? surf.fRightSide : surf.fLeftSide;

    // Check if the current state is exited on this surface. This is true if
    // the in_state does not match the default state for the surface.
    bool exiting = surf.fDefaultState != in_navind;

    if (exiting) {
      // This is an exiting surface for in_state
      for (auto ind = 0; ind < exit_side.fNsurf; ++ind) {
        auto const &framedsurf = surfdata.fFramedSurf[exit_side.fSurfaces[ind]];
        // we need to check the current framed surface if its state is either the same or
        // is an ancestor of the in_state
        bool same_state = framedsurf.fState == in_navind;
        bool is_sibling = same_state || vecgeom::NavStateIndex::IsDescendentImpl(framedsurf.fState, in_navind);
        if (is_sibling) {
          // convert to local frame
          Vector3D<Real_t> local_frame = onsurf;
          if (framedsurf.fTrans) local_frame = surfdata.fGlobalTrans[framedsurf.fTrans].Transform(onsurf);
          bool inframe = framedsurf.fFrame.Inside(local_frame, surfdata);
          // When exiting the frame with same state, the exit point must be in frame
          // When exiting a daughter volume frame, the exit point must NOT be in the daughter frame
          if (inframe != same_state) {
            can_hit = false;
            break; // next candidate
          }
        }
      }
      if (!can_hit) continue;
      // the current state is correctly exited, so there is a transition on this surface
      distance   = dist;
      isurfcross = cand[icand];
      // backup exited state
      out_state.SetLastExited();
      // the default next navigation index is the one of the common state for the surface
      out_state.SetNavIndex(surf.fDefaultState);
    }

    // Now check if something is being entered

    // the onsurf local point can be used as input for a side search optimization structure.
    // for now just loop candidates in order. Since candidates are sorted by depth, the first
    // frame entry is the good one.
    for (auto ind = 0; ind < entry_side.fNsurf; ++ind) {
      auto const &framedsurf = surfdata.fFramedSurf[entry_side.fSurfaces[ind]];
      // convert to local frame
      Vector3D<Real_t> local_frame = onsurf;
      if (framedsurf.fTrans) local_frame = surfdata.fGlobalTrans[framedsurf.fTrans].Transform(onsurf);
      bool inframe = framedsurf.fFrame.Inside(local_frame, surfdata);
      if (inframe) {
        // the first hit frame is the good one
        distance   = dist;
        isurfcross = cand[icand];
        out_state.SetLastExited();
        out_state.SetNavIndex(framedsurf.fState);
        break;
      }
    }
  } // end loop on candidates

  assert(isurfcross != 0); // something MUST be hit (since no step limitation at this point)
  exit_surf = isurfcross;

  return distance;
}

} // namespace protonav
} // namespace vgbrep
#endif