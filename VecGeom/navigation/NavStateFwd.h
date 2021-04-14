/// \file NavStateFwd.h
/// \author Andrei Gheata (andrei.gheata@cern.ch)
/// \date 12.03.2014

#ifndef VECGEOM_NAVIGATION_NAVSTATEFWD_H_
#define VECGEOM_NAVIGATION_NAVSTATEFWD_H_

#ifdef VECGEOM_USE_NAVINDEX
#define NavigationStateImpl NavStateIndex
#else
#define NavigationStateImpl NavStatePath
#endif
namespace vecgeom {

VECGEOM_HOST_FORWARD_DECLARE(class NavigationStateImpl;);
VECGEOM_HOST_FORWARD_DECLARE(using NavigationState = NavigationStateImpl;);

VECGEOM_DEVICE_FORWARD_DECLARE(class NavigationStateImpl;);
VECGEOM_DEVICE_FORWARD_DECLARE(using NavigationState = NavigationStateImpl;);

inline namespace VECGEOM_IMPL_NAMESPACE {

class NavStateIndex;
class NavStatePath;   // Needed even when not the navigation state, at least for GeoVisitor and NavIndexTable.
using NavigationState = NavigationStateImpl;

} // namespace VECGEOM_IMPL_NAMESPACE

} // namespace vecgeom

#undef NavigationStateImpl

#endif // VECGEOM_NAVIGATION_NAVSTATEFWD_H_
