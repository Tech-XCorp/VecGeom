/*
 * TestExportToROOT.cpp
 *
 *  Created on: 28.10.2014
 *      Author: swenzel
 */

#include "volumes/PlacedVolume.h"
#include "base/Global.h"
#include "base/Transformation3D.h"
#include "volumes/UnplacedBox.h"
#include "volumes/UnplacedTube.h"
#include "volumes/UnplacedTrd.h"
#include "volumes/UnplacedOrb.h"
#include "volumes/UnplacedParaboloid.h"
#include "volumes/UnplacedParallelepiped.h"
#include "management/RootGeoManager.h"
#include "management/GeoManager.h"
#include "TGeoManager.h"
#include <iostream>

using namespace vecgeom;


// create a VecGeom geometry

VPlacedVolume* SetupGeometry() {
  UnplacedBox *worldUnplaced = new UnplacedBox(10, 10, 10);
  UnplacedBox *boxUnplaced = new UnplacedBox(2.5, 2.5, 2.5);

  UnplacedTube *tube1Unplaced = new UnplacedTube( 0.5, 1., 0.5, 0., kTwoPi);
  UnplacedTube *tube2Unplaced = new UnplacedTube( 0.5, 1., 0.5, 0., kPi);

  UnplacedTrd *trdUnplaced = new UnplacedTrd( 0.1, 0.2, 0.15, 0.05 );

  UnplacedOrb *orbUnplaced = new UnplacedOrb( 0.1 );
  UnplacedParaboloid *paraUnplaced = new UnplacedParaboloid( 0.1, 0.2, 0.1 );
  UnplacedParallelepiped *epipedUnplaced =  new UnplacedParallelepiped( 0.1, 0.05, 0.05, 0.2, 0.4, 0.1 );

  Transformation3D *placement1 = new Transformation3D( 5,  5,  5,  0,  0,  0);
  Transformation3D *placement2 = new Transformation3D(-5,  5,  5, 45,  0,  0);
  Transformation3D *placement3 = new Transformation3D( 5, -5,  5,  0, 45,  0);
  Transformation3D *placement4 = new Transformation3D( 5,  5, -5,  0,  0, 45);
  Transformation3D *placement5 = new Transformation3D(-5, -5,  5, 45, 45,  0);
  Transformation3D *placement6 = new Transformation3D(-5,  5, -5, 45,  0, 45);
  Transformation3D *placement7 = new Transformation3D( 5, -5, -5,  0, 45, 45);
  Transformation3D *placement8 = new Transformation3D(-5, -5, -5, 45, 45, 45);

  Transformation3D *placement9 = new Transformation3D(-0.5,-0.5,-0.5,0,0,0);
  Transformation3D *placement10 = new Transformation3D(0.5,0.5,0.5,0,45,0);
  Transformation3D *idendity    = new Transformation3D();

  LogicalVolume *world = new LogicalVolume("world",worldUnplaced);
  LogicalVolume *box =   new LogicalVolume("lbox1",boxUnplaced);
  LogicalVolume *tube1 = new LogicalVolume("ltube1", tube1Unplaced);
  LogicalVolume *tube2 = new LogicalVolume("ltube2", tube2Unplaced);
  LogicalVolume *trd1  = new LogicalVolume("ltrd", trdUnplaced);

  LogicalVolume *orb1 = new LogicalVolume("lorb1", orbUnplaced);
  LogicalVolume *parab1 = new LogicalVolume("lparab1", paraUnplaced);
  LogicalVolume *epip1 = new LogicalVolume("lepip1", epipedUnplaced);

  world->PlaceDaughter(orb1, idendity);
  trd1->PlaceDaughter(parab1, idendity);
  world->PlaceDaughter(epip1, idendity);

  tube1->PlaceDaughter( trd1, idendity );
  box->PlaceDaughter( tube1, placement9 );
  box->PlaceDaughter( tube2, placement10 );

  world->PlaceDaughter(box, placement1);
  world->PlaceDaughter(box, placement2);
  world->PlaceDaughter(box, placement3);
  world->PlaceDaughter(box, placement4);
  world->PlaceDaughter(box, placement5);
  world->PlaceDaughter(box, placement6);
  world->PlaceDaughter(box, placement7);
  world->PlaceDaughter(box, placement8);
  return world->Place();
}


int main()
{
    VPlacedVolume const * world = SetupGeometry();
    GeoManager::Instance().SetWorld(world);
    GeoManager::Instance().CloseGeometry();
    int md1 = GeoManager::Instance().getMaxDepth();
    int mpv1 = GeoManager::Instance().GetPlacedVolumesCount();
    int mlv1 = GeoManager::Instance().GetLogicalVolumesCount();
    int ntotalnodes1 = GeoManager::Instance().GetTotalNodeCount();

    // exporting to ROOT file
    RootGeoManager::Instance().ExportToROOTGeometry( world, "geom1.root" );

    assert( ::gGeoManager->GetNNodes() == ntotalnodes1 );
    assert( ::gGeoManager->GetListOfVolumes()->GetEntries() == mlv1 );

    //
    RootGeoManager::Instance().Clear();

    // now try to read back in
    RootGeoManager::Instance().set_verbose(1);
    RootGeoManager::Instance().LoadRootGeometry("geom1.root");

    //// see if everything was restored
    // RootGeoManager::Instance().world()->logical_volume()->PrintContent(0);

    int md2 = GeoManager::Instance().getMaxDepth();
    int mpv2 = GeoManager::Instance().GetPlacedVolumesCount();
    int mlv2 = GeoManager::Instance().GetLogicalVolumesCount();
    int ntotalnodes2 = GeoManager::Instance().GetTotalNodeCount();

    assert( md2 == md1 );
    assert( mpv2 == mpv1 );
    assert( mlv2 == mlv1 );
    assert( mpv2 > 0);
    assert( mlv2 > 0);
    assert( ntotalnodes1 == ntotalnodes2 );

    return 0;
}
