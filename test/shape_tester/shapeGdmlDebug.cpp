//
// File:    shapeGdmlDebug.cpp
//
// Purpose: Loads a GDML geometry, and shapes are accessible and queried through their names.
//  A point (x,y,z) and a direction (vx,vy,vz), and the shape, point and
//  track is drawn using ROOT, from original point to next intersection with shape,
//  and all distances and safeties are compared with ROOT.
//
//  Note: ROOT is required for visualization.
//        Geant4 is also used when available, but it is not mandatory.

#include "VecGeom/management/RootGeoManager.h"
#include "VecGeom/volumes/LogicalVolume.h"
#include "VecGeom/volumes/PlacedVolume.h"
#include "VecGeom/volumes/UnplacedBox.h"
#include "VecGeom/benchmarking/Benchmarker.h"
#include "VecGeom/management/GeoManager.h"
#include "VecGeom/volumes/UnplacedBox.h"
#include "test/core/SetupGDMLGeometry.h"

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoBBox.h"
#include "TGeoTube.h"
#include "utilities/Visualizer.h"
#include <string>
#include <cmath>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#ifdef VECGEOM_GEANT4
#include "G4ThreeVector.hh"
#include "G4VSolid.hh"
#endif

using namespace vecgeom;

void showDaughterNames(LogicalVolume const* logmom) {
  const Vector<VPlacedVolume const*>& daughters = logmom->GetDaughters();
  cout << "Volume: <"<< logmom->GetName() <<">, # daughters="<< daughters.size() << endl;
  int i = 0;
  for (Vector<VPlacedVolume const*>::const_iterator idau = daughters.begin(); idau != daughters.end(); ++idau, ++i) {
    VPlacedVolume const* dau = *idau;
    auto *logdau = dau->GetLogicalVolume();
    cout << i <<' '<< dau->GetName() <<" of type "<< logdau->GetName() <<", #daus="<< logdau->GetDaughtersp()->size() << endl;
  }
}

void checkOverlaps(LogicalVolume const* logmom) {
  const Vector<VPlacedVolume const*>& daughters = logmom->GetDaughters();
  cout << "Volume: <"<< logmom->GetName() <<">, # daughters="<< daughters.size() << endl;
  int i = 0;
  for (Vector<VPlacedVolume const*>::const_iterator idau = daughters.begin(); idau != daughters.end(); ++idau, ++i) {
    VPlacedVolume const* dau = *idau;
  }
}

int main(int argc, char *argv[])
{
  if (argc < 9) {
    std::cerr << "Fixed shape in source code - user needs to give a local point + local dir\n";
    std::cerr << "example: " << argv[0] << " cms2021.gdml HBModule 10.0 0.8 -3.5 1 0 0\n";
    return 1;
  }

  // Load geometry from GDML file into gGeoManagerr
  //bool loaded = SetupGDMLGeometry(argv[1]);
  bool loaded = TGeoManager::Import(argv[1]);

  // from ROOT into VecGeom geometry
  RootGeoManager::Instance().LoadRootGeometry();
  double px   = atof(argv[3]);
  double py   = atof(argv[4]);
  double pz   = atof(argv[5]);
  double dirx = atof(argv[6]);
  double diry = atof(argv[7]);
  double dirz = atof(argv[8]);


  GeoManager& geoMgr = GeoManager::Instance();
  cout <<" # logical volumes: "<< geoMgr.GetRegisteredVolumesCount() << endl;
  cout <<" # placed volumes: "<< geoMgr.GetPlacedVolumesCount() << endl;
  cout <<" # nodes: "<< geoMgr.GetTotalNodeCount() << endl;
  cout <<" MaxDepth: "<< geoMgr.getMaxDepth() << endl;

  //geoMgr.GetWorld()->PrintContent();

  LogicalVolume const* vglogvol = geoMgr.FindLogicalVolume(argv[2]);
  cerr<<"spot 1, pvol="<< vglogvol <<"\n";
  showDaughterNames(vglogvol);

  checkOverlaps(vglogvol);

  VPlacedVolume const* vgpvol = geoMgr.FindPlacedVolume(argv[2]);
  cerr<<"spot 2, pvol="<< vgpvol <<"\n";

  VUnplacedVolume const* vgshape = vglogvol->GetUnplacedVolume();
  cerr<<"spot 3, pvol="<< vglogvol <<"\n";
  (void)vgshape;

  auto *rootShape = gGeoManager->FindVolumeFast(argv[2])->GetShape();
  cerr<<"spot 4, rootVol="<< rootShape <<"\n";

  // now get the VecGeom shape and benchmark it
  if (loaded && rootShape) {
    Vector3D<Precision> point(px, py, pz);
    Vector3D<Precision> dir(dirx, diry, dirz);
    // normalize direction vector
    dir = dir.Unit();

    SOA3D<Precision> pointcontainer(4);
    pointcontainer.resize(4);
    SOA3D<Precision> dircontainer(4);
    dircontainer.resize(4);
    Precision *output = new Precision[4];
    Precision *steps  = new Precision[4];

    if (argc > 9) {
      pointcontainer.set(0, point);
      dircontainer.set(0, dir.x(), dir.y(), dir.z());
      double px2   = atof(argv[9]);
      double py2   = atof(argv[10]);
      double pz2   = atof(argv[11]);
      double dirx2 = atof(argv[12]);
      double diry2 = atof(argv[13]);
      double dirz2 = atof(argv[14]);
      pointcontainer.set(1, px2, py2, pz2);
      dircontainer.set(1, dirx2, diry2, dirz2);
      pointcontainer.set(2, point);
      dircontainer.set(2, dir.x(), dir.y(), dir.z());
      pointcontainer.set(3, px2, py2, pz2);
      dircontainer.set(3, dirx2, diry2, dirz2);

      for (auto i = 0; i < 4; ++i) {
        steps[i] = vecgeom::kInfLength;
      }

    } else {
      for (auto i = 0; i < 4; ++i) {
        pointcontainer.set(i, point);
        dircontainer.set(i, dir.x(), dir.y(), dir.z());
        steps[i] = vecgeom::kInfLength;
      }
    }
    if (!dir.IsNormalized()) {
      std::cerr << "** Attention: Direction is not normalized **\n";
      std::cerr << "** Direction differs from 1 by "
                << std::sqrt(dir.x() * dir.x() + dir.y() * dir.y() + dir.z() * dir.z()) - 1. << "\n";
    }
    double dist;
    //std::cout << "VecGeom Capacity " << const_cast<VUnplacedVolume*>(vgpvol)->Capacity() << "\n";
    std::cout << "VecGeom CONTAINS " << vgpvol->Contains(point) << "\n";
    std::cout << "VecGeom INSIDE " << vgpvol->Inside(point) << "\n";
    dist = vgpvol->DistanceToIn(point, dir);
    std::cout << "VecGeom DI " << dist << "\n";
    if (dist >= vecgeom::kTolerance && dist < vecgeom::kInfLength) {
      std::cout << "VecGeom INSIDE(p=p+dist*dir) " << vgpvol->Inside(point + dir * dist) << "\n";
      if (vgpvol->Inside(point + dir * dist) == vecgeom::kOutside)
        std::cout << "VecGeom Distance seems to be too big  DI(p=p+dist*dir,-dir) "
                  << vgpvol->DistanceToIn(point + dir * dist, -dir) << "\n";
      if (vgpvol->Inside(point + dir * dist) == vecgeom::kInside)
        std::cout << "VecGeom Distance seems to be too small DO(p=p+dist*dir,dir) "
                  << vgpvol->DistanceToOut(point + dir * dist, dir) << "\n";
    }
    vgpvol->DistanceToIn(pointcontainer, dircontainer, steps, output);
    std::cout << "VecGeom DI-V " << output[0] << "\n";
    std::cout << "VecGeom DO " << vgpvol->DistanceToOut(point, dir) << "\n";
    vgpvol->DistanceToOut(pointcontainer, dircontainer, steps, output);
    std::cout << "VecGeom DO-V " << output[0] << "\n";

    std::cout << "VecGeom SI " << vgpvol->SafetyToIn(point) << "\n";
    std::cout << "VecGeom SO " << vgpvol->SafetyToOut(point) << "\n";

    std::cout << "ROOT Capacity " << rootShape->Capacity() << "\n";
    std::cout << "ROOT CONTAINS " << rootShape->Contains(&point[0]) << "\n";
    std::cout << "ROOT DI " << rootShape->DistFromOutside(&point[0], &dir[0]) << "\n";
    std::cout << "ROOT DO " << rootShape->DistFromInside(&point[0], &dir[0]) << "\n";
    std::cout << "ROOT SI " << rootShape->Safety(&point[0], false) << "\n";
    std::cout << "ROOT SO " << rootShape->Safety(&point[0], true) << "\n";

    TGeoShape const *rootback = vgpvol->ConvertToRoot();
    if (rootback) {
      std::cout << "ROOTBACKCONV CONTAINS " << rootback->Contains(&point[0]) << "\n";
      std::cout << "ROOTBACKCONV DI " << rootback->DistFromOutside(&point[0], &dir[0]) << "\n";
      std::cout << "ROOTBACKCONV DO " << rootback->DistFromInside(&point[0], &dir[0]) << "\n";
      std::cout << "ROOTBACKCONV SI " << rootback->Safety(&point[0], false) << "\n";
      std::cout << "ROOTBACKCONV SO " << rootback->Safety(&point[0], true) << "\n";
    } else {
      std::cerr << "ROOT backconversion failed\n";
    }

#ifdef VECGEOM_GEANT4
    G4ThreeVector g4p(point.x(), point.y(), point.z());
    G4ThreeVector g4d(dir.x(), dir.y(), dir.z());

    G4VSolid const *g4solid = vgpvol->ConvertToGeant4();
    if (g4solid != NULL) {
      std::cout << "G4 CONTAINS " << g4solid->Inside(g4p) << "\n";
      std::cout << "G4 DI " << g4solid->DistanceToIn(g4p, g4d) << "\n";
      std::cout << "G4 DO " << g4solid->DistanceToOut(g4p, g4d) << "\n";
      std::cout << "G4 SI " << g4solid->DistanceToIn(g4p) << "\n";
      std::cout << "G4 SO " << g4solid->DistanceToOut(g4p) << "\n";
    } else {
      std::cerr << "G4 conversion failed\n";
    }
#endif

    double step = 0;
    //    if( rootShape->Contains( &point[0] ) ){
    //      step = rootShape->DistFromInside( &point[0], &dir[0] );
    //    }
    //    else {
    //      step = rootShape->DistFromOutside( &point[0], &dir[0] );
    //    }

    // modified to show problem in DistanceToIn()
    step = rootShape->DistFromOutside(&point[0], &dir[0]);
    Visualizer visualizer;
    visualizer.AddVolume(*vgpvol);
    visualizer.AddPoint(point);
    visualizer.AddLine(point, point + step * dir);
    visualizer.Show();
  } else {
    std::cerr << " Error: problems constructing volume [" << rootShape << "] ... EXITING\n";
    return 1;
  }
}
