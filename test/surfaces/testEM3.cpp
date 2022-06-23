#include <iostream>
#include <string>

#include "test/benchmark/ArgParser.h"
#include <VecGeom/volumes/LogicalVolume.h>
#include <VecGeom/volumes/Box.h>
#include <VecGeom/management/GeoManager.h>
#include <VecGeom/base/RNG.h>

#include <VecGeom/surfaces/Model.h>
#include <VecGeom/surfaces/BrepHelper.h>
#include <VecGeom/surfaces/Navigator.h>
#include <VecGeom/navigation/GlobalLocator.h>
#include <VecGeom/navigation/NewSimpleNavigator.h>
#include <VecGeom/volumes/utilities/VolumeUtilities.h>
#include <VecGeom/navigation/NavStatePool.h>
#include "VecGeom/base/Stopwatch.h"

using namespace vecgeom;
// Forwards
void CreateVecGeomWorld(int, int);
bool ValidateNavigation(int, int, vgbrep::SurfData<vecgeom::Precision> const &);
void TestPerformance(int, int, vgbrep::SurfData<vecgeom::Precision> const &);
double PropagateRay(vecgeom::Vector3D<vecgeom::Precision> const &, vecgeom::Vector3D<vecgeom::Precision> const &,
                    vgbrep::SurfData<vecgeom::Precision> const &);

int main(int argc, char *argv[])
{
  using BrepHelper = vgbrep::BrepHelper<vecgeom::Precision>;

  OPTION_INT(layers, 50);
  OPTION_INT(absorbers, 2.);
  OPTION_INT(nvalidate, 10000);
  OPTION_INT(nbench, 1000000);
  OPTION_INT(verbose, 0);

  CreateVecGeomWorld(layers, absorbers);

  BrepHelper::Instance().SetVerbosity(verbose);

  if (!BrepHelper::Instance().CreateLocalSurfaces()) return 1;
  if (!BrepHelper::Instance().CreateCommonSurfacesFlatTop()) return 2;

  ValidateNavigation(nvalidate, layers, BrepHelper::Instance().GetSurfData());

  TestPerformance(nbench, layers, BrepHelper::Instance().GetSurfData());
  /*
  // Attempt some basic navigation
//  vecgeom::Vector3D<vecgeom::Precision> start(-8 * 1.1 * layers, 0, 0);
  vecgeom::Vector3D<vecgeom::Precision> start(-1.7, 0, -30);
   for (int i = 0; i < nrays; ++i) {
    vecgeom::Vector3D<vecgeom::Precision> direction;
    direction[2] = 1;
    direction[0] = 0.01 * (1. - 2 * vecgeom::RNG::Instance().uniform());
    direction[1] = 0.01 * (1. - 2 * vecgeom::RNG::Instance().uniform());
    direction.Normalize();
    PropagateRay(start, direction, BrepHelper::Instance().GetSurfData());
  }
  */
  // Test clearing surface data
  BrepHelper::Instance().ClearData();
}

void CreateVecGeomWorld(int NbOfLayers, int NbOfAbsorbers)
{
  const double CalorSizeYZ       = 40;
  const double GapThickness      = 2.3;
  const double AbsorberThickness = 5.7;

  const double LayerThickness = GapThickness + AbsorberThickness;
  const double CalorThickness = NbOfLayers * LayerThickness;

  const double WorldSizeX  = 1.2 * CalorThickness + 1000;
  const double WorldSizeYZ = 1.2 * CalorSizeYZ + 1000;

  auto worldSolid = new vecgeom::UnplacedBox(0.5 * WorldSizeX, 0.5 * WorldSizeYZ, 0.5 * WorldSizeYZ);
  auto worldLogic = new vecgeom::LogicalVolume("World", worldSolid);
  vecgeom::VPlacedVolume *worldPlaced = worldLogic->Place();

  //
  // Calorimeter
  //
  auto calorSolid = new vecgeom::UnplacedBox(0.5 * CalorThickness, 0.5 * CalorSizeYZ, 0.5 * CalorSizeYZ);
  auto calorLogic = new vecgeom::LogicalVolume("Calorimeter", calorSolid);
  vecgeom::Transformation3D origin;
  worldLogic->PlaceDaughter(calorLogic, &origin);

  //
  // Layers
  //
  auto layerSolid = new vecgeom::UnplacedBox(0.5 * LayerThickness, 0.5 * CalorSizeYZ, 0.5 * CalorSizeYZ);

  //
  // Absorbers
  //
  auto gapSolid = new vecgeom::UnplacedBox(0.5 * GapThickness, 0.5 * CalorSizeYZ, 0.5 * CalorSizeYZ);
  auto gapLogic = new vecgeom::LogicalVolume("Gap", gapSolid);
  vecgeom::Transformation3D gapPlacement(-0.5 * LayerThickness + 0.5 * GapThickness, 0, 0);

  auto absorberSolid = new vecgeom::UnplacedBox(0.5 * AbsorberThickness, 0.5 * CalorSizeYZ, 0.5 * CalorSizeYZ);
  auto absorberLogic = new vecgeom::LogicalVolume("Absorber", absorberSolid);
  vecgeom::Transformation3D absorberPlacement(0.5 * LayerThickness - 0.5 * AbsorberThickness, 0, 0);

  // Create a new LogicalVolume per layer, we need unique IDs for scoring.
  double xCenter = -0.5 * CalorThickness + 0.5 * LayerThickness;
  for (int i = 0; i < NbOfLayers; i++) {
    std::string name("Layer_");
    name += std::to_string(i);
    auto layerLogic = new vecgeom::LogicalVolume("Layer", layerSolid);
    vecgeom::Transformation3D placement(xCenter, 0, 0);
    calorLogic->PlaceDaughter(name.c_str(), layerLogic, &placement);

    layerLogic->PlaceDaughter(gapLogic, &gapPlacement);
    layerLogic->PlaceDaughter(absorberLogic, &absorberPlacement);

    xCenter += LayerThickness;
  }

  vecgeom::GeoManager::Instance().SetWorldAndClose(worldPlaced);
}

double PropagateRay(vecgeom::Vector3D<vecgeom::Precision> const &point,
                    vecgeom::Vector3D<vecgeom::Precision> const &direction,
                    vgbrep::SurfData<vecgeom::Precision> const &surfdata)
{
  // Locate the start point. This is not yet implemented in the surface model
  NavStateIndex in_state, out_state;
  int exit_surf   = 0;
  double dist_tot = 0;
  GlobalLocator::LocateGlobalPoint(GeoManager::Instance().GetWorld(), point, in_state, true);
  auto pt = point;
  printf("start: ");
  in_state.Print();
  do {
    auto distance = vgbrep::protonav::ComputeStepAndHit(pt, direction, in_state, out_state, surfdata, exit_surf);
    if (exit_surf != 0) {
      dist_tot += distance;
      pt += distance * direction;
    }
    printf("surface %d at dist = %g: ", exit_surf, distance);
    out_state.Print();
    in_state = out_state;
  } while (!out_state.IsOutside());

  return dist_tot;
}

bool ValidateNavigation(int npoints, int nbLayers, vgbrep::SurfData<vecgeom::Precision> const &surfdata)
{
  // prepare tracks to be used for benchmarking
  constexpr double tolerance = 10 * vecgeom::kTolerance;
  int num_errors             = 0;
  SOA3D<Precision> points(npoints);
  SOA3D<Precision> dirs(npoints);

  Vector3D<Precision> samplingVolume(nbLayers * 8, 20, 20);
  vecgeom::volumeUtilities::FillRandomPoints(samplingVolume, points);
  vecgeom::volumeUtilities::FillRandomDirections(dirs);

  // now setup all the navigation states
  int ndeep = GeoManager::Instance().getMaxDepth();
  NavStatePool origStates(npoints, ndeep);
  NavStatePool outputStates(npoints, ndeep);

  vecgeom::Precision *refSteps = new Precision[npoints];
  memset(refSteps, 0, sizeof(Precision) * npoints);

  auto *nav = vecgeom::NewSimpleNavigator<>::Instance();
  for (int i = 0; i < npoints; ++i) {
    Vector3D<Precision> const &pos = points[i];
    Vector3D<Precision> const &dir = dirs[i];
    GlobalLocator::LocateGlobalPoint(GeoManager::Instance().GetWorld(), pos, *origStates[i], true);
    nav->FindNextBoundaryAndStep(pos, dir, *origStates[i], *outputStates[i], vecgeom::kInfLength, refSteps[i]);
    // printf("dist: %g  ", refSteps[i]); origStates[i]->Print();
    // outputStates[i]->Print();

    // shoot the same ray in the surface model
    int exit_surf = 0;
    NavStateIndex out_state;
    auto distance = vgbrep::protonav::ComputeStepAndHit(pos, dir, *origStates[i], out_state, surfdata, exit_surf);
    if (out_state.GetNavIndex() != outputStates[i]->GetNavIndex() || std::abs(distance - refSteps[i]) > tolerance)
      num_errors++;
  }

  printf("=== Validation: num_erros = %d / %d\n", num_errors, npoints);
  delete[] refSteps;

  return num_errors == 0;
}

void TestPerformance(int npoints, int nbLayers, vgbrep::SurfData<vecgeom::Precision> const &surfdata)
{
  SOA3D<Precision> points(npoints);
  SOA3D<Precision> dirs(npoints);

  Vector3D<Precision> samplingVolume(nbLayers * 8, 20, 20);
  vecgeom::volumeUtilities::FillRandomPoints(samplingVolume, points);
  vecgeom::volumeUtilities::FillRandomDirections(dirs);

  // now setup all the navigation states
  int ndeep = GeoManager::Instance().getMaxDepth();
  NavStatePool origStates(npoints, ndeep);
  NavStateIndex out_state;
  vecgeom::Precision distance = 0;
  auto *nav                   = vecgeom::NewSimpleNavigator<>::Instance();

  // Locate all input points, without timing this operation
  for (int i = 0; i < npoints; ++i) {
    Vector3D<Precision> const &pos = points[i];
    GlobalLocator::LocateGlobalPoint(GeoManager::Instance().GetWorld(), pos, *origStates[i], true);
  }

  // Benchamrk primitive-based NewSimpleNavigator
  Stopwatch timer;
  timer.Start();
  for (int i = 0; i < npoints; ++i) {
    Vector3D<Precision> const &pos = points[i];
    Vector3D<Precision> const &dir = dirs[i];
    nav->FindNextBoundaryAndStep(pos, dir, *origStates[i], out_state, vecgeom::kInfLength, distance);
  }
  Precision time_prim = timer.Stop();

  int exit_surf = 0;
  Stopwatch timer1;
  timer1.Start();
  for (int i = 0; i < npoints; ++i) {
    Vector3D<Precision> const &pos = points[i];
    Vector3D<Precision> const &dir = dirs[i];
    distance = vgbrep::protonav::ComputeStepAndHit(pos, dir, *origStates[i], out_state, surfdata, exit_surf);
  }
  Precision time_surf = timer1.Stop();

  printf("Time for %d points: NewSimpleNavigator = %f [s]  vgbrep::protonav = %f\n", npoints, time_prim, time_surf);
}