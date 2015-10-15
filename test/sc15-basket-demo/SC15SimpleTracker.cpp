/// author: Andrei Gheata (andrei.gheata@cern.ch) Oct 14, 2015
#include <sstream>
#include <string>
#include "volumes/Tube.h"
#include "volumes/LogicalVolume.h"
#include "management/GeoManager.h"
#include "SC15Navigators.h"
#include "navigation/SimpleNavigator.h"
#include "base/Vector3D.h"
#include <iostream>
using namespace vecgeom;

#define SETINNERNAV(depth) \
  if(layer==depth) vol->SetUserExtensionPtr( (void*) InnerMostTubeNavigator<depth+1>::Instance() );

#define SETLAYERNAV(depth) \
  if(layer==depth) vol->SetUserExtensionPtr( (void*) LayerNavigator<depth+1>::Instance() );

#define SETWORLDNAV(depth) \
  if(layer==0) vol->SetUserExtensionPtr( (void*) WorldNavigator<true>::Instance() );

static VNavigator const* GetNavigator(LogicalVolume const*lvol){
  return (VNavigator const*)lvol->GetUserExtensionPtr();
}

void AssignNavigatorToVolume( LogicalVolume *vol, int layer, int maxlayers ){
  assert( maxlayers <= 10 ); // we are only listing 10 template specializations up to depth 10 here
  if( layer == maxlayers - 1 ) // it is the last layer
  {
    SETINNERNAV(0);
    SETINNERNAV(1);
    SETINNERNAV(2);
    SETINNERNAV(3);
    SETINNERNAV(4);
    SETINNERNAV(5);
    SETINNERNAV(6);
    SETINNERNAV(7);
    SETINNERNAV(8);
    SETINNERNAV(9);
    SETINNERNAV(10);
  }
  if( layer < maxlayers - 1 ){
    SETLAYERNAV(0);
    SETLAYERNAV(1);
    SETLAYERNAV(2);
    SETLAYERNAV(3);
    SETLAYERNAV(4);
    SETLAYERNAV(5);
    SETLAYERNAV(6);
    SETLAYERNAV(7);
    SETLAYERNAV(8);
    SETLAYERNAV(9);
    SETLAYERNAV(10);
  }
}


VPlacedVolume *CreateSimpleTracker(int nlayers) {
  std::cout << "Creating SimpleTracker geometry having " << nlayers << " layers"
            << std::endl;
  // World size
  const double world_size = 50.;

  // Top volume
  UnplacedBox *uTop = new UnplacedBox(world_size + 2, world_size + 2, world_size + 2);
  LogicalVolume *top = new LogicalVolume("world", uTop);
  top->SetUserExtensionPtr( (void *) WorldNavigator<true>::Instance() );
  
// Cylindrical layers
  double rmax;
  std::string layerBase = "layer_";
  double deltaR = world_size/nlayers;

  LogicalVolume *mother = top;
  for (auto layer=0; layer<nlayers; ++layer) {
    rmax = world_size - layer * deltaR;
    std::ostringstream layerName;
    layerName << layerBase << layer;
    UnplacedTube *uLayer = new UnplacedTube(0, rmax, world_size, 0, kTwoPi);
    LogicalVolume *layerVol = new LogicalVolume(layerName.str().c_str(), uLayer);
    
    AssignNavigatorToVolume(layerVol, layer, nlayers);

    // Place in mother
    mother->PlaceDaughter(layerName.str().c_str(), layerVol, &Transformation3D::kIdentity);
    // change mother to current lvol (to make a real hierarchy)
    mother = layerVol;
  }  

  VPlacedVolume *world = top->Place();
  GeoManager::Instance().SetWorld(world);
  GeoManager::Instance().CloseGeometry();
  return world;
}

// a test function to verify correct functioning of the navigators
// for a simple test case
void TestScalarNavigation() {
  // setup point and direction in world
  Vector3D<Precision> p(-51.,0,0);
  Vector3D<Precision> dir(1.,0,0);

  // init navstates
  NavigationState * curnavstate = NavigationState::MakeInstance(GeoManager::Instance().getMaxDepth());
  NavigationState * newnavstate = NavigationState::MakeInstance(GeoManager::Instance().getMaxDepth());
  
  SimpleNavigator nav;
  nav.LocatePoint( GeoManager::Instance().GetWorld(), p, *curnavstate, true );

  while( ! curnavstate->IsOutside() ) {
    //
    std::cout << "tracking in " << curnavstate->Top()->GetLogicalVolume()->GetName() << "\n";
      // get navigator object and move point
    VNavigator const *specialnav = GetNavigator(curnavstate->Top()->GetLogicalVolume());
    double step = specialnav->ComputeStepAndPropagatedState(p, dir, kInfinity, *curnavstate, *newnavstate);
    std::cout << "step " << step << "\n";
    p = p + dir * (step + 1E-6);

    // pointer swap is enough
    auto *tmp = curnavstate;
    curnavstate = newnavstate;
    newnavstate = tmp;
  }
}

void TestVectorNavigation() {
  // to be filled in

}

// reproducing the pixel-by-pixel XRayBenchmark
// target: show speed gain from specialized navigators
// in comparision to current XRay-Benchmark
void XRayBenchmark() {
  // Sofia?

}

void BasketBasedXRayBenchmark() {
  // to be filled in by Andrei

}


int main(int argc, char* argv[]) {
  int nlayers = 10;
  if (argc > 1) {
    nlayers = std::atoi(argv[1]);
  }
  CreateSimpleTracker(nlayers);

  // loop over all logical volumes and print navigator
  std::vector<LogicalVolume *> lvols;
  GeoManager::Instance().GetAllLogicalVolumes(lvols);

  for( auto v : lvols ){
    auto nav = (VNavigator*) v->GetUserExtensionPtr();
    std::cerr << v->GetName() << " has navigator " << nav->GetName() << "\n";
  }


  // test tracking
  TestScalarNavigation();
}
