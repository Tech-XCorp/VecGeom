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
  UnplacedOrb *orbUnplaced = new UnplacedOrb(2.8);

  LogicalVolume *world = new LogicalVolume("world",worldUnplaced);
  LogicalVolume *trap  = new LogicalVolume("trap",trapUnplaced);
  LogicalVolume *box   = new LogicalVolume("box",boxUnplaced);
  LogicalVolume *orb   = new LogicalVolume("orb",orbUnplaced);

  Transformation3D *ident = new Transformation3D( 0,  0,  0,  0,  0,  0);
  orb->PlaceDaughter("orb1",box, ident);
  trap->PlaceDaughter("box1",orb, ident);

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
  OPTION_STRING(geometry, "navBench.root");
  OPTION_STRING(testVolume, "world");
  OPTION_DOUBLE(fraction, 0.8f);
#ifdef VECGEOM_ROOT
  OPTION_BOOL(vis, false);
#endif

  // default values used above are always printed.  If help true, stop now, so user will know which options
  // are available, and what the default values are.
  OPTION_BOOL(help, false);
  if(help) return 0;

  const VPlacedVolume *world = NULL;
  if(geometry.compare("navBench.root")==0) {
    world = SetupGeometry();

#ifdef VECGEOM_ROOT
    // Exporting to ROOT file
    RootGeoManager::Instance().ExportToROOTGeometry( world, "navBench.root" );
    RootGeoManager::Instance().Clear();
#endif
  }

  // Now try to read back in.  This is needed to make comparisons to VecGeom easily,
  // since it builds VecGeom geometry based on the ROOT geometry and its TGeoNodes.
#ifdef VECGEOM_ROOT
  RootGeoManager::Instance().set_verbose(0);
  RootGeoManager::Instance().LoadRootGeometry(geometry.c_str());
#endif

  // Visualization
#ifdef VECGEOM_ROOT
  if(vis) {  // note that visualization block returns, excluding the rest of benchmark
    Visualizer visualizer;
    const VPlacedVolume* world = GeoManager::Instance().GetWorld();
    world = GeoManager::Instance().FindPlacedVolume(testVolume.c_str());
    visualizer.AddVolume( *world );

    Vector<Daughter> const* daughters = world->logical_volume()->daughtersp();
    for(int i=0; i<daughters->size(); ++i) {
      VPlacedVolume const* daughter = (*daughters)[i];
      Transformation3D const& trf1 = *(daughter->transformation());
      visualizer.AddVolume(*daughter, trf1);

      // Vector<Daughter> const* daughters2 = daughter->logical_volume()->daughtersp();
      // for(int ii=0; ii<daughters2->size(); ++ii) {
      //   VPlacedVolume const* daughter2 = (*daughters2)[ii];
      //   Transformation3D const& trf2 = *(daughter2->transformation());
      //   Transformation3D comb = trf1;
      //   comb.MultiplyFromRight(trf2);
      //   visualizer.AddVolume(*daughter2, comb);
      // }
    }

    visualizer.Show();
    return 0;
  }
#endif

  //testVectorSafety(world);

  std::cout<<"\n*** Validating VecGeom navigation..."<< std::endl;

  const VPlacedVolume* startVolume = GeoManager::Instance().GetWorld();
  if( testVolume.compare("world")!=0 ) {
    startVolume = GeoManager::Instance().FindPlacedVolume(testVolume.c_str());
  }

  std::cout<<"NavigationBenchmark: testVolume=<"<< testVolume
           <<">, startVolume="<< (startVolume ? startVolume->GetLabel() : NULL)
           <<" - "<< *startVolume <<"\n";

  int np = Min( npoints, 1000 );  // no more than 1000 points used for validation
  SOA3D<Precision> points(np);
  SOA3D<Precision> dirs(np);
  SOA3D<Precision> locpts(np);

  vecgeom::volumeUtilities::FillGlobalPointsAndDirectionsForLogicalVolume( startVolume->logical_volume(), locpts, points, dirs, fraction, np);

  bool ok = validateVecGeomNavigation(np, points, dirs);

  // Must be validated before being benchmarked
  if(!ok) {
    std::cout<<"VecGeom validation failed."<< std::endl;
    return 1;
  }

  std::cout<<"VecGeom validation passed."<< std::endl;

  // on mic.fnal.gov CPUs, loop execution takes ~70sec for npoints=10M
  while(npoints<=10000) {
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
