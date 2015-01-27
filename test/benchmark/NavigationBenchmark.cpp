/*
 * testGPUNavigation.cpp
 *
 *  Created on: Oct 25, 2014
 *      Author: swenzel, lima
 */

#include "benchmarking/NavigationBenchmarker.h"
#include "ArgParser.h"
#include "volumes/utilities/VolumeUtilities.h"

#ifdef VECGEOM_ROOT
#include "management/RootGeoManager.h"
#include "utilities/Visualizer.h"
#endif

#include "management/GeoManager.h"
#include "volumes/Box.h"
#include "volumes/Orb.h"
#include "volumes/Trapezoid.h"

using namespace VECGEOM_NAMESPACE;

VPlacedVolume* SetupGeometry() {

  UnplacedBox *worldUnplaced = new UnplacedBox(10, 10, 10);
  UnplacedTrapezoid *trapUnplaced = new UnplacedTrapezoid(4,0,0,4,4,4,0,4,4,4,0);
  UnplacedBox *boxUnplaced = new UnplacedBox(2.5, 2.5, 2.5);
  UnplacedOrb *orbUnplaced = new UnplacedOrb(2.4);

  LogicalVolume *world = new LogicalVolume("world",worldUnplaced);
  LogicalVolume *trap  = new LogicalVolume("trap",trapUnplaced);
  LogicalVolume *box   = new LogicalVolume("box",boxUnplaced);
  LogicalVolume *orb   = new LogicalVolume("orb",orbUnplaced);

  Transformation3D *ident = new Transformation3D( 0,  0,  0,  0,  0,  0);
  box->PlaceDaughter("orb1",orb, ident);
  trap->PlaceDaughter("box1",box, ident);

  Transformation3D *placement1 = new Transformation3D( 5,  5,  5,  0,  0,  0);
  Transformation3D *placement2 = new Transformation3D(-5,  5,  5,  0,  0,  0); //45,  0,  0);
  Transformation3D *placement3 = new Transformation3D( 5, -5,  5,  0,  0,  0); // 0, 45,  0);
  Transformation3D *placement4 = new Transformation3D( 5,  5, -5,  0,  0,  0); // 0,  0, 45);
  Transformation3D *placement5 = new Transformation3D(-5, -5,  5,  0,  0,  0); //45, 45,  0);
  Transformation3D *placement6 = new Transformation3D(-5,  5, -5,  0,  0,  0); //45,  0, 45);
  Transformation3D *placement7 = new Transformation3D( 5, -5, -5,  0,  0,  0); // 0, 45, 45);
  Transformation3D *placement8 = new Transformation3D(-5, -5, -5,  0,  0,  0); //45, 45, 45);

  world->PlaceDaughter("trap1",trap, placement1);
  world->PlaceDaughter("trap2",trap, placement2);
  world->PlaceDaughter("trap3",trap, placement3);
  world->PlaceDaughter("trap4",trap, placement4);
  world->PlaceDaughter("trap5",trap, placement5);
  world->PlaceDaughter("trap6",trap, placement6);
  world->PlaceDaughter("trap7",trap, placement7);
  world->PlaceDaughter("trap8",trap, placement8);

  VPlacedVolume  * w = world->Place();
  GeoManager::Instance().SetWorld(w);
  GeoManager::Instance().CloseGeometry();
  return w;
}


int main(int argc, char* argv[])
{
  OPTION_INT(npoints, 10000);
  OPTION_INT(nreps, 3);
#ifdef VECGEOM_ROOT
  OPTION_BOOL(vis, false);
#endif

  VPlacedVolume *world = SetupGeometry();

#ifdef VECGEOM_ROOT
  if(vis) {  // note that visualization block returns, excluding the rest of benchmark
    Visualizer visualizer;
    visualizer.AddVolume(*world);

    Vector<Daughter> const* daughters = world->logical_volume()->daughtersp();
    for(int i=0; i<daughters->size(); ++i) {
      VPlacedVolume const* daughter = (*daughters)[i];
      visualizer.AddVolume(*daughter, *(daughter->transformation()));
    }

    visualizer.Show();
    return 0;
  }
#endif

  testVectorSafety(world);

#ifdef VECGEOM_ROOT
  // Exporting to ROOT file
  RootGeoManager::Instance().ExportToROOTGeometry( world, "navBench.root" );
  RootGeoManager::Instance().Clear();

  // Now try to read back in.  This is needed to make comparisons to VecGeom easily,
  // since it builds VecGeom geometry based on the ROOT geometry and its TGeoNodes.
  RootGeoManager::Instance().set_verbose(0);
  RootGeoManager::Instance().LoadRootGeometry("navBench.root");
#endif

  std::cout<<"\n*** Validating VecGeom navigation..."<< std::endl;
  Precision fraction = 0.7;


  //*** These are for Ex03.root
  //fraction=0.8;  const VPlacedVolume* startVolume = GeoManager::Instance().GetWorld();  // OK
  //fraction=0.0;  const VPlacedVolume* startVolume = GeoManager::Instance().FindPlacedVolume("Lead_1");  // ok with a fraction = 0.0
  //fraction=0.0;  const VPlacedVolume* startVolume = GeoManager::Instance().FindPlacedVolume("Calorimeter_13"); // can't find uncontained points
  //fraction=0.0;  const VPlacedVolume* startVolume = GeoManager::Instance().FindPlacedVolume("Layer_5");   // can't find uncontained points
  //fraction=0.0; const VPlacedVolume* startVolume = GeoManager::Instance().FindPlacedVolume("liquidArgon_2"); // can't find uncontained points

  //*** Here for navBench.root
  //fraction=0.8;  const VPlacedVolume* startVolume = GeoManager::Instance().GetWorld();  // OK
  fraction=0.8; const VPlacedVolume* startVolume = GeoManager::Instance().FindPlacedVolume("trap_0"); // OK
  //fraction=0.5; const VPlacedVolume* startVolume = GeoManager::Instance().FindPlacedVolume("box_0");  // OK
  //fraction=0.0; const VPlacedVolume* startVolume = GeoManager::Instance().FindPlacedVolume("orb_0");  // OK

  printf("startVolume = %s\n", startVolume->GetLabel().c_str());

  int np = Min( npoints, 1000 );  // no more than 1000 points used for validation
  SOA3D<Precision> points(np);
  SOA3D<Precision> dirs(np);
  SOA3D<Precision> locpts(np);

  vecgeom::volumeUtilities::FillGlobalPointsAndDirectionsForLogicalVolume( startVolume->logical_volume(), locpts, points, dirs, fraction, np);
  //vecgeom::volumeUtilities::FillRandomPoints( *startVolume, points );
  //vecgeom::volumeUtilities::FillRandomDirections( dirs );

  bool ok = validateVecGeomNavigation(np, points, dirs);

  // Must be validated before being benchmarked
  if(!ok) {
    std::cout<<"VecGeom validation failed."<< std::endl;
    return 1;
  }

  std::cout<<"VecGeom validation passed."<< std::endl;

  // on mic.fnal.gov CPUs, loop execution takes ~70sec for npoints=10M
  while(npoints<=10000000) {
    std::cout<<"\n*** Running navigation benchmarks with npoints="<<npoints<<" and nreps="<< nreps <<".\n";
    runNavigationBenchmarks(startVolume, npoints, nreps);
    npoints*=10;
  }


/*
// GPU part
  int nDevice;
  cudaGetDeviceCount(&nDevice);

  if(nDevice > 0) {
    cudaDeviceReset();
  }
  else {
    std::cout << "No Cuda Capable Device ... " << std::endl;
    return 0;
  }
*/

  return 0;
}
