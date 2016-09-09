#include "volumes/LogicalVolume.h"
#include "volumes/Box.h"
#include "volumes/Trapezoid.h"
#include "benchmarking/Benchmarker.h"
#include "management/GeoManager.h"
#include "ArgParser.h"
#include <fstream>

using namespace vecgeom;

int main(int argc, char *argv[])
{
  OPTION_INT(npoints, 32768);
  OPTION_INT(nrep, 4);
  OPTION_INT(nstats, 1);

#if 0
  // temporary alert to deprecated use of npoints
  OPTION_INT(npoints, 0);
  if(npoints) {
    printf("\n***** ERROR: -npoints is now deprecated.  Please use -npoints instead.\n");
    std::exit(-1);
  }
#endif

  //===== Build a geometry with trapezoid(s)

  // world volume is a box, or a box-like trapezoid
  // UnplacedBox worldUnplaced = UnplacedBox(20., 20., 20.);
  UnplacedTrapezoid worldUnplaced = UnplacedTrapezoid(20., 0., 0., 20., 20., 20., 0., 20., 20., 20., 0.);

  //-- and here is for an internal trapezoid

  // validate construtor for input corner points -- add an xy-offset for non-zero theta,phi
  TrapCorners xyz;
  Precision xoffset = 9;
  Precision yoffset = -6;

  // define corner points
  // convention: p0(---); p1(+--); p2(-+-); p3(++-); p4(--+); p5(+-+); p6(-++); p7(+++)
  xyz[0] = Vector3D<Precision>(-2 + xoffset, -5 + yoffset, -15);
  xyz[1] = Vector3D<Precision>(2 + xoffset, -5 + yoffset, -15);
  xyz[2] = Vector3D<Precision>(-3 + xoffset, 5 + yoffset, -15);
  xyz[3] = Vector3D<Precision>(3 + xoffset, 5 + yoffset, -15);
  xyz[4] = Vector3D<Precision>(-4 - xoffset, -10 - yoffset, 15);
  xyz[5] = Vector3D<Precision>(4 - xoffset, -10 - yoffset, 15);
  xyz[6] = Vector3D<Precision>(-6 - xoffset, 10 - yoffset, 15);
  xyz[7] = Vector3D<Precision>(6 - xoffset, 10 - yoffset, 15);

  // create trapezoid
  UnplacedTrapezoid trapUnplaced(xyz);

  //.. and here is for a secon internal trapezoid
  // UnplacedTrapezoid trapUnplaced2(10, 0, 0, 10, 10, 10, 0, 10, 10, 10, 0);

  LogicalVolume world("world", &worldUnplaced);
  LogicalVolume trap("trap", &trapUnplaced);

  Transformation3D *transf = new Transformation3D(5, 2, 3, 15, 30, 45);
  world.PlaceDaughter(&trap, transf);
  // world.PlaceDaughter(&trap, &Transformation3D::kIdentity);

  VPlacedVolume *worldPlaced = world.Place();

  GeoManager::Instance().SetWorldAndClose(worldPlaced);

  Benchmarker tester(GeoManager::Instance().GetWorld());
  tester.SetPointCount(npoints);
  tester.SetRepetitions(nrep);
  tester.SetTolerance(1.e-9);
  tester.SetPoolMultiplier(1);

  //=== Here is for the validation + one perf data point displayed on screen
  tester.SetVerbosity(3);
  tester.SetMeasurementCount(1);
  auto errcode = tester.RunBenchmark();

  // clear benchmark results, so previous measurements won't be written out into the output .csv file
  tester.ClearResults();

  // Now run to collect statistics for performance plots - written to the .csv output file only
  if (nstats > 1) {
    // Idea is to start at npoints=2, and then increase it by x2 at a time until maxNpoints is reached
    tester.SetVerbosity(0);
    tester.SetMeasurementCount(nstats);
    npoints        = 2;
    int maxNpoints = 2048;
#ifdef VECGEOM_CUDA
    maxNpoints = 1048576;
#endif
    while (npoints <= maxNpoints) {
      tester.SetPointCount(npoints);
      tester.RunBenchmark();
      npoints *= 2;
    }

    // Save statistics data to a text file
    std::list<BenchmarkResult> results = tester.PopResults();
    std::ofstream outStream;
    outStream.open("trapBenchmarkData.csv", std::fstream::app);
    BenchmarkResult::WriteCsvHeader(outStream);
    for (auto i = results.begin(), iEnd = results.end(); i != iEnd; ++i) {
      i->WriteToCsv(outStream);
    }
    outStream.close();
  }

  // cleanup
  if (transf) delete transf;
  return errcode;
}
