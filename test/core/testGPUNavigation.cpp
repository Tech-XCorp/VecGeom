/*
 * testGPUNavigation.cpp
 *
 *  Created on: Oct 25, 2014
 *      Author: swenzel, lima
 */

#include "volumes/utilities/VolumeUtilities.h"
#include "volumes/Box.h"
#include "volumes/PlacedVolume.h"
#include "base/Transformation3D.h"
#include "base/SOA3D.h"
#include "management/RootGeoManager.h"
#include "navigation/NavigationState.h"
#include "navigation/SimpleNavigator.h"
#include "management/GeoManager.h"
#include "base/Global.h"
#include "test/benchmark/ArgParser.h"
#include "base/Stopwatch.h"

#include "navigationgpu.h"
#include "TGeoNavigator.h"
#include "TGeoManager.h"

#include <cuda.h>
#include <cuda_runtime.h>

using namespace VECGEOM_NAMESPACE;

VPlacedVolume* SetupBoxGeometry() {
  UnplacedBox *worldUnplaced = new UnplacedBox(10, 10, 10);
  UnplacedBox *boxUnplaced = new UnplacedBox(0.5, 0.5, 0.5);
  Transformation3D *placement1 = new Transformation3D( 2,  2,  2,  0,  0,  0);
  Transformation3D *placement2 = new Transformation3D(-2,  2,  2, 45,  0,  0);
  Transformation3D *placement3 = new Transformation3D( 2, -2,  2,  0, 45,  0);
  Transformation3D *placement4 = new Transformation3D( 2,  2, -2,  0,  0, 45);
  Transformation3D *placement5 = new Transformation3D(-2, -2,  2, 45, 45,  0);
  Transformation3D *placement6 = new Transformation3D(-2,  2, -2, 45,  0, 45);
  Transformation3D *placement7 = new Transformation3D( 2, -2, -2,  0, 45, 45);
  Transformation3D *placement8 = new Transformation3D(-2, -2, -2, 45, 45, 45);
  LogicalVolume *world = new LogicalVolume("world",worldUnplaced);
  LogicalVolume *box = new LogicalVolume("box",boxUnplaced);
  world->PlaceDaughter(box, placement1);
  world->PlaceDaughter(box, placement2);
  world->PlaceDaughter(box, placement3);
  world->PlaceDaughter(box, placement4);
  world->PlaceDaughter(box, placement5);
  world->PlaceDaughter(box, placement6);
  world->PlaceDaughter(box, placement7);
  world->PlaceDaughter(box, placement8);
  VPlacedVolume  * w = world->Place();
  GeoManager::Instance().SetWorld(w);
  GeoManager::Instance().CloseGeometry();
  return w;
}

// function to test safety
void testVectorSafety( VPlacedVolume* world ){
   SOA3D<Precision> points(1024);
   SOA3D<Precision> workspace(1024);
   Precision * safeties = (Precision *) _mm_malloc(sizeof(Precision)*1024,32);
   vecgeom::volumeUtilities::FillUncontainedPoints( *world, points );

   // now setup all the navigation states
   NavigationState ** states = new NavigationState*[1024];
   vecgeom::SimpleNavigator nav;
   for (int i=0;i<1024;++i){
       states[i]=new NavigationState( GeoManager::Instance().getMaxDepth() );
       nav.LocatePoint( world, points[i], *states[i], true);
   }

    // calculate safeties with vector interface
    nav.GetSafeties(points, states, workspace, safeties );

    // verify against serial interface
    for (int i=0;i<1024;++i){
        vecgeom::Assert( safeties[i] == nav.GetSafety( points[i], *states[i] ), ""
                " Problem in VectorSafety (in SimpleNavigator)" );
    }
    std::cout << "Safety test passed\n";
   _mm_free(safeties);
}



// function to test vector navigator
void testVectorNavigator( VPlacedVolume* world, int np ){
   SOA3D<Precision> points(np);
   SOA3D<Precision> dirs(np);
   SOA3D<Precision> workspace1(np);
   SOA3D<Precision> workspace2(np);

   Precision * steps = (Precision *) _mm_malloc(sizeof(Precision)*np,32);
   Precision * pSteps = (Precision *) _mm_malloc(sizeof(Precision)*np,32);
   Precision * GPUSteps = (Precision *) _mm_malloc(sizeof(Precision)*np,32);
   Precision * safeties = (Precision *) _mm_malloc(sizeof(Precision)*np,32);

   int * intworkspace = (int *) _mm_malloc(sizeof(int)*np,32);

   vecgeom::volumeUtilities::FillUncontainedPoints( *world, points );
   vecgeom::volumeUtilities::FillRandomDirections( dirs );

   vecgeom::SimpleNavigator nav;
   Stopwatch timer;

   // now setup all the navigation states
   NavigationState ** states = new NavigationState*[np];
   //NavigationState ** newstates = new NavigationState*[np];

   for (int i=0;i<np;++i){
     // pSteps[i] = kInfinity;
     pSteps[i] = (i%2)? 1 : kInfinity;
     states[i] = new NavigationState( GeoManager::Instance().getMaxDepth() );
   //    newstates[i] = new NavigationState( GeoManager::Instance().getMaxDepth() );
   }

   for (int i=0;i<np;++i) {
     nav.LocatePoint( world, points[i], *states[i], true);
   }

   timer.Start();
   for(int i=0; i<np; ++i) {
     steps[i] = (Precision)i;
   }

   // calculate steps with vector interface
   // nav.FindNextBoundaryAndStep( points, dirs, workspace1, workspace2,
   //         states, newstates, pSteps, safeties, steps, intworkspace );
   Precision elapsedCPU = timer.Stop();
   printf("CPU elapsed time: %f ms\n", 1000.*elapsedCPU);

   std::cout<<" testGPUNav: gGeoManager = "<< ::gGeoManager << std::endl;
   TGeoNavigator * rootnav = ::gGeoManager->GetCurrentNavigator();

   // verify against serial interface
   for (int i=0;i<np;++i) {
     Precision s;
     NavigationState cmp( GeoManager::Instance().getMaxDepth() );
     cmp.Clear();
     nav.FindNextBoundaryAndStep( points[i], dirs[i], *states[i],
               cmp, pSteps[i], s );
     // vecgeom::Assert( steps[i] == s ,
     //         " Problem in VectorNavigation (steps) (in SimpleNavigator)" );
     // vecgeom::Assert( cmp.Top() == newstates[i]->Top() ,
     //                " Problem in VectorNavigation (states) (in SimpleNavigator)" );
     // vecgeom::Assert( cmp.IsOnBoundary() == newstates[i]->IsOnBoundary(),
     //                " Problem in VectorNavigation (boundary) (in SimpleNavigator)" );
     // vecgeom::Assert( safeties[i] == nav.GetSafety( points[i], *states[i] ),
     //         " Problem with safety " );
   // }

//=== ROOT Navigation

   // for (int i=0;i<np;++i) {
     Vector3D<Precision> const& pos = points[i];
     Vector3D<Precision> const& dir = dirs[i];
     rootnav->ResetState();

     //TGeoNode * node = 
     rootnav->FindNode( pos.x(), pos.y(), pos.z() );

     rootnav->SetCurrentPoint( pos.x(), pos.y(), pos.z() );
     rootnav->SetCurrentDirection( dir.x(), dir.y(), dir.z() );
//     path->UpdateNavigator( rootnav );
     rootnav->FindNextBoundaryAndStep( 1E30 );

     // if( cmp.Top() != NULL ) {
       // if( rootnav->GetCurrentNode()
       //     != RootGeoManager::Instance().tgeonode( cmp.Top() ) ) {
             std::cerr << "ERROR ON ITERATION " << i;
             std::cerr <<" pos = " << pos;
             std::cerr <<" dir = " << dir << "\n";
             std::cerr << "ROOT GOES HERE: " << rootnav->GetCurrentNode()->GetName();
             std::cerr << " with step: "<< rootnav->GetStep() << "\n";
             if(cmp.Top()==NULL) std::cerr<<"VECGEOM GOES TO <NULL> !!!\n";
             else {
               std::cerr << "VECGEOM GOES TO <" << RootGeoManager::Instance().GetName( cmp.Top() ) << ">\n";
               nav.InspectEnvironmentForPointAndDirection( pos, dir, *states[i] );
             }
        // }
      // }
   }

//=== Comparing with 

#ifdef VECGEOM_CUDA
   printf("Start GPU\n");
   RunNavigationCuda(world, np,
                     points.x(),  points.y(), points.z(),
                     dirs.x(), dirs.y(), dirs.z(), pSteps, GPUSteps );
#endif
   unsigned mismatches=0;
   for (int i=0;i<np;++i) {
     if( abs(steps[i]-GPUSteps[i]) > kTolerance ) {
       ++mismatches;
       std::cerr << "i " << i << " steps " << steps[i] << " CUDA steps " << GPUSteps[i]
                 << " - diff="<< (steps[i]-GPUSteps[i]) <<"\n";
     //  vecgeom::Assert( steps[i] == GPUSteps[i], " Problem in CUDA Navigator " );
     }
   }

   if(mismatches>0) std::cout << "Navigation test failed with "<< mismatches <<" mismatches\n";
   else std::cout<<"Navigation test passed.\n";
   _mm_free(steps);
   _mm_free(intworkspace);
   _mm_free(pSteps);
   _mm_free(safeties);
}


int main(int argc, char* argv[])
{
  OPTION_INT(npoints, 100);
  VPlacedVolume *world = SetupBoxGeometry();

  testVectorSafety(world);

  // exporting to ROOT file
  RootGeoManager::Instance().ExportToROOTGeometry( world, "geom1.root" );

  assert( ::gGeoManager->GetNNodes() == ntotalnodes1 );
  assert( ::gGeoManager->GetListOfVolumes()->GetEntries() == mlv1 );

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

  testVectorNavigator(world, npoints);
}
