/// @file ParaboloidBenchmark.cpp
/// @author Marilena Bandieramonte (marilena.bandieramonte@cern.ch)

#include "volumes/logical_volume.h"
#include "volumes/box.h"
#include "volumes/Paraboloid.h"
#include "benchmarking/Benchmarker.h"
#include "management/geo_manager.h"

using namespace vecgeom;

int main() {
    
    std::cout<<"Paraboloid Benchmark\n";
    UnplacedBox worldUnplaced = UnplacedBox(10., 10., 10.);
    UnplacedParaboloid paraboloidUnplaced = UnplacedParaboloid(3., 5., 7.); //rlo=3. - rhi=5. dz=7
    std::cout<<"Paraboloid created\n";
    LogicalVolume world = LogicalVolume("MBworld", &worldUnplaced);
    LogicalVolume paraboloid = LogicalVolume("paraboloid", &paraboloidUnplaced);
    world.PlaceDaughter(&paraboloid, &Transformation3D::kIdentity);
    VPlacedVolume *worldPlaced = world.Place();
    GeoManager::Instance().set_world(worldPlaced);
    std::cout<<"World set\n";
    
    
    Benchmarker tester(GeoManager::Instance().world());
    tester.SetVerbosity(3);
    tester.SetPointCount(1<<13);
    tester.SetRepetitions(1024);
    std::cout<<"Prepared to run benchmarker\n";
    tester.RunBenchmark();
    
    return 0;
}

