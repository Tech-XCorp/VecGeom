/*
 * ReducedPolycone.cpp
 *
 *  Created on: Sep 20, 2018
 *      Author: ramansehgal
 */

#include "volumes/ReducedPolycone.h"
#include <iostream>
#include "base/Vector.h"

namespace vecgeom {

VECCORE_ATT_HOST_DEVICE
ReducedPolycone::ReducedPolycone(Vector<Vector2D<Precision>> rzVect)
{
  SetRZ(rzVect);
}

VECCORE_ATT_HOST_DEVICE
void ReducedPolycone::SetRMax()
{
  fRMax = fRZVect[0].x();
  for (unsigned int i = 1; i < fRZVect.size(); i++) {
    if (fRMax < fRZVect[i].x()) fRMax = fRZVect[i].x();
  }
  fRMax += 10.;
}

VECCORE_ATT_HOST_DEVICE
void ReducedPolycone::SetRZ(Vector<Vector2D<Precision>> rzVect)
{
  fRZVect.clear();
  for (unsigned int i = 0; i < rzVect.size(); i++) {
    fRZVect.push_back(rzVect[i]);
  }
  SetRMax();
}

VECCORE_ATT_HOST_DEVICE
ReducedPolycone::~ReducedPolycone() {}

VECCORE_ATT_HOST_DEVICE
void ReducedPolycone::ConvertToUniqueVector(Vector<Precision> &vect)
{
  /* Logic to sort : Currently using bubble sort, but if required
   * then some more efficient sorting algorithm can be used.
   *
   * Ideally this algo should be the part of Vector.h
   */
  for (unsigned int i = 0; i < vect.size(); i++) {
    for (unsigned int j = 0; j < vect.size() - 1; j++) {
      if (vect[j] >= vect[j + 1]) {
        Precision temp = vect[j];
        vect[j]        = vect[j + 1];
        vect[j + 1]    = temp;
      }
    }
  }

  // Copy the sorted vector to temporary vector to remove the duplicate
  Vector<Precision> tempVect;
  tempVect.push_back(vect[0]);
  for (unsigned int i = 1; i < vect.size(); i++) {
    if (tempVect[tempVect.size() - 1] != vect[i]) tempVect.push_back(vect[i]);
  }

  // Copying temporary vector back to original vector
  vect.clear();
  for (unsigned int i = 0; i < tempVect.size(); i++) {
    vect.push_back(tempVect[i]);
  }
}

VECCORE_ATT_HOST_DEVICE
Vector<Precision> ReducedPolycone::GetUniqueZVector()
{
  Vector<Precision> z;
  for (unsigned int i = 0; i < fRZVect.size(); i++) {
    z.push_back(fRZVect[i].y());
  }
  ConvertToUniqueVector(z);
  return z;
}

VECCORE_ATT_HOST_DEVICE
bool ReducedPolycone::GetLineIntersection(Precision p0_x, Precision p0_y, Precision p1_x, Precision p1_y,
                                          Precision p2_x, Precision p2_y, Precision p3_x, Precision p3_y,
                                          Precision *i_x, Precision *i_y)
{

  Precision s1_x, s1_y, s2_x, s2_y;
  s1_x = p1_x - p0_x;
  s1_y = p1_y - p0_y;
  s2_x = p3_x - p2_x;
  s2_y = p3_y - p2_y;

  if (s1_y == 0. && s2_y == 0.) return false;

  Precision s, t;
  Precision deno = (-s2_x * s1_y + s1_x * s2_y);
  if (deno == 0.) return false;
  s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / deno;
  t = (s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

  if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
    // Collision detected
    if (i_x != NULL) *i_x = p0_x + (t * s1_x);
    if (i_y != NULL) *i_y = p0_y + (t * s1_y);

    if (p0_x == p1_x && p2_y == p3_y) {
      *i_x = p0_x;
      *i_y = p2_y;
    }
    return true;
  }

  return false; // No collision
}
VECCORE_ATT_HOST_DEVICE
bool ReducedPolycone::GetLineIntersection(Line2D l1, Line2D l2)
{
  Vector2D<Precision> poi(0., 0.);
  return GetLineIntersection(l1.p1.x(), l1.p1.y(), l1.p2.x(), l1.p2.y(), l2.p1.x(), l2.p1.y(), l2.p2.x(), l2.p2.y(),
                             &poi.x(), &poi.y());
}
VECCORE_ATT_HOST_DEVICE
Vector<Line2D> ReducedPolycone::GetLineVector()
{
  Vector<Line2D> lineVect;
  for (unsigned int i = 0; i < fRZVect.size(); i++) {
    if (i == (fRZVect.size() - 1))
      lineVect.push_back(Line2D(fRZVect[i], fRZVect[0]));
    else
      lineVect.push_back(Line2D(fRZVect[i], fRZVect[i + 1]));
  }
  return lineVect;
}

VECCORE_ATT_HOST_DEVICE
void ReducedPolycone::CalcPoIVectorFor2DPolygon(Vector<Vector2D<Precision>> &poiVect, Vector<Precision> z)
{

  for (unsigned int i = 0; i < fRZVect.size(); i++) {
    for (unsigned int j = 0; j < z.size(); j++) {
      Vector2D<Precision> poi;
      bool valid = false;
      if (i == (fRZVect.size() - 1)) {
        valid = GetLineIntersection(Line2D(fRZVect[i], fRZVect[0]),
                                    Line2D(Vector2D<Precision>(0., z[j]), Vector2D<Precision>(fRMax, z[j])), poi);
      } else {
        valid = GetLineIntersection(Line2D(fRZVect[i], fRZVect[i + 1]),
                                    Line2D(Vector2D<Precision>(0., z[j]), Vector2D<Precision>(fRMax, z[j])), poi);
      }
      if (valid) {
        poiVect.push_back(poi);
      }
    }
  }
}
VECCORE_ATT_HOST_DEVICE
bool ReducedPolycone::Contour(Vector<Precision> z)
{
  bool contour = ContourCheck(z);
  // std::cout << "Ctonour : "  << contour << std::endl;
  if (!contour) {
#ifndef VECCORE_CUDA
    std::cerr << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
              << "@@@@ Polycone CAN'T handle contours of specified type @@@@ \n"
              << "@@@@ Kindly use GenericPolycone of Geant4             @@@@\n"
              << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
#endif
    return contour;
  }

  // Getting vector of all the line
  Vector<Line2D> lineVect = GetLineVector();

  for (int unsigned i = 2; i < lineVect.size(); i++) {
    for (unsigned int j = 0; j <= (i - 2); j++) {
      if (i == (lineVect.size() - 1)) {
        if (j == 0) {
          Vector2D<Precision> poi;
          contour &= GetLineIntersection(lineVect[i], lineVect[j], poi);
          contour &= (poi.x() == lineVect[0].p1.x()) && (poi.y() == lineVect[0].p1.y());
        } else {
          bool test = GetLineIntersection(lineVect[i], lineVect[j]);
          contour &= !test;
        }
      } else
        contour &= !GetLineIntersection(lineVect[i], lineVect[j]);
    }
  }
  return contour;
}
VECCORE_ATT_HOST_DEVICE
bool ReducedPolycone::ContourCheck(Vector<Precision> z)
{

  // Getting vector of all the line
  Vector<Line2D> lineVect = GetLineVector();

  // Creating vector of checkerLines
  Vector<Line2D> checkerLineVect;
  for (unsigned int i = 0; i < z.size() - 1; i++) {
    Precision zval = (z[i] + z[i + 1]) / 2.;
    checkerLineVect.push_back(Line2D(Vector2D<Precision>(0., zval), Vector2D<Precision>(fRMax, zval)));
  }

  /* Now iterating over above two line vectors to find whether the contour is
   * suitable to get converted to a polycone or not. */
  Vector<Vector2D<Precision>> poiVect;
  bool check = true;
  for (unsigned int i = 0; i < checkerLineVect.size(); i++) {
    poiVect.clear();
    for (unsigned int j = 0; j < lineVect.size(); j++) {
      Vector2D<Precision> poi;
      bool valid = false;
      valid      = GetLineIntersection(checkerLineVect[i], lineVect[j], poi);
      if (valid) {
        poiVect.push_back(poi);
      }
    }
    check &= (poiVect.size() == 2);
    if (!check) {
#ifndef VECCORE_CUDA
      //    std::cerr << "Not proper contour detected by checkerLine index : " << i << std::endl;
#endif
      break;
    }
  }
  return check;
}

VECCORE_ATT_HOST_DEVICE
bool ReducedPolycone::PointExist(Vector2D<Precision> pt)
{
  bool exist = false;
  for (unsigned int i = 0; i < fRZVect.size(); i++) {
    exist = ((fRZVect[i].x() == pt.x()) && (fRZVect[i].y() == pt.y()));
    if (exist) break;
  }
  return exist;
}

VECCORE_ATT_HOST_DEVICE
void ReducedPolycone::CreateNewContour()
{
  Vector<Precision> z;
  for (unsigned int i = 0; i < fRZVect.size(); i++) {
    z.push_back(fRZVect[i].y());
  }
  ConvertToUniqueVector(z);
  int numOfIterationsForContourModification = z.size();
  for (int i = 0; i < numOfIterationsForContourModification; i++) {
    Vector<Line2D> lineVect = GetLineVector();
    Vector<Vector2D<Precision>> modifiedRZ;
    Vector<Vector2D<Precision>> poiVect;
    for (unsigned int j = 0; j < lineVect.size(); j++) {
      Vector2D<Precision> poi;
      bool valid = GetLineIntersection(lineVect[j], Line2D(fRMax, z[i]), poi);
      if (valid) {
        if (!PointExist(poi)) {
          /* Modify Contour and add the PoI in the proper sequence.
           * Algo : Scan the lineVect and find the index of line that corresponds to
           *        line which gives new PoI, and then insert PoI and lineVect[j].p2
           */
          modifiedRZ.push_back(poi);
          modifiedRZ.push_back(lineVect[j].p2);
        } else {
          modifiedRZ.push_back(lineVect[j].p2);
        }
      } else {
        modifiedRZ.push_back(lineVect[j].p2);
      }
    }
    // Change the contour
    fRZVect.clear();
    fRZVect = modifiedRZ;
    modifiedRZ.clear();
  }
}

VECCORE_ATT_HOST_DEVICE
Vector<Vector<Precision>> ReducedPolycone::GetRandZVectorAtDiffZ(Vector<Vector2D<Precision>> poiVect,
                                                                 Vector<Precision> &dz)
{

  Vector<Precision> zVect;
  Vector<Precision> rVect;
  for (unsigned int i = 0; i < poiVect.size(); i++) {
    rVect.push_back(poiVect[i].x());
    zVect.push_back(poiVect[i].y());
  }
  ConvertToUniqueVector(zVect);

  Vector<Vector<Precision>> supVect;
  for (unsigned int i = 0; i < zVect.size(); i++) {
    Vector<Precision> rVect;
    for (unsigned int j = 0; j < poiVect.size(); j++) {
      if (poiVect[j].y() == zVect[i]) rVect.push_back(poiVect[j].x());
    }
    supVect.push_back(rVect);
    dz.push_back(zVect[i]);
  }
  return supVect;
}

VECCORE_ATT_HOST_DEVICE
void ReducedPolycone::ProcessContour(Vector<Precision> z)
{

  Vector<Precision> zV;
  Vector<Vector2D<Precision>> poiVect;
  CalcPoIVectorFor2DPolygon(poiVect, z);
  Vector<Vector<Precision>> supVect = GetRandZVectorAtDiffZ(poiVect, zV);
  // Make Unique R vector at different Z
  for (unsigned int i = 0; i < supVect.size(); i++)
    ConvertToUniqueVector(supVect[i]);

  assert(supVect.size() == z.size() && "Inconsistent vector sizes");
  Vector<Line2D> lineVect = GetLineVector();
  for (unsigned int i = 0; i < z.size() - 1; i++) {
    Vector<Line2D> sectionLine;
    int counter = 0;
    for (unsigned int j = 0; j < supVect[i].size(); j++) {

      Precision rval = supVect[i][j];
      Precision zval = z[i];
      for (unsigned int k = 0; k < lineVect.size(); k++) {
        Vector2D<Precision> poi;
        bool valid = false;
        Line2D line;
        bool cond = (lineVect[k].p1.x() == rval && lineVect[k].p1.y() == zval) ||
                    (lineVect[k].p2.x() == rval && lineVect[k].p2.y() == zval);
        if (cond) {
          line  = lineVect[k];
          valid = GetLineIntersection(
              line, Line2D(Vector2D<Precision>(0., z[i + 1]), Vector2D<Precision>(fRMax, z[i + 1])), poi);
          if (valid) {
            counter++;
            sectionLine.push_back(line);
          }
        }
      }
    }
    assert(sectionLine.size() == 2 && "Got more than two lines for a section");

    CreateSectionFromTwoLines(sectionLine[0], sectionLine[1]);
  }
}

VECCORE_ATT_HOST_DEVICE
void ReducedPolycone::CreateSectionFromTwoLines(Line2D l1, Line2D l2)
{
  Precision rmin1 = 0., rmin2 = 0., rmax1 = 0., rmax2 = 0., z1 = 0., z2 = 0.;

  if (l1.p1.y() < l1.p2.y()) {
    rmin1 = l1.p1.x();
    rmax1 = l2.p2.x();
    z1    = l1.p1.y();
    rmin2 = l1.p2.x();
    rmax2 = l2.p1.x();
    z2    = l1.p2.y();
  } else {
    rmin1 = l1.p2.x();
    rmax1 = l2.p1.x();
    z1    = l1.p2.y();
    rmin2 = l1.p1.x();
    rmax2 = l2.p2.x();
    z2    = l1.p1.y();
  }

  if (rmin2 > rmax2) Swap(rmin2, rmax2);
  if (rmin1 > rmax1) Swap(rmin1, rmax1);

  Section s(rmin1, rmax1, z1, rmin2, rmax2, z2);
  fSectionVect.push_back(s);
}

VECCORE_ATT_HOST_DEVICE
void ReducedPolycone::Swap(Precision &a, Precision &b)
{
  Precision temp;
  temp = a;
  a    = b;
  b    = temp;
}

VECCORE_ATT_HOST_DEVICE
bool ReducedPolycone::Check()
{
  CreateNewContour();
  Vector<Precision> zVect = GetUniqueZVector();
  bool contour            = Contour(zVect);
  if (contour) {
    ProcessContour(zVect);
  }
  return contour;
}
VECCORE_ATT_HOST_DEVICE
bool ReducedPolycone::GetPolyconeParameters(Vector<Precision> &rmin, Vector<Precision> &rmax, Vector<Precision> &z)
{
  bool contour = Check();
  if (contour) {
    // ProcessContour(rzVect,zVect);
    for (unsigned int i = 0; i < fSectionVect.size(); i++) {
      rmin.push_back(fSectionVect[i].rMin1);
      rmin.push_back(fSectionVect[i].rMin2);
      rmax.push_back(fSectionVect[i].rMax1);
      rmax.push_back(fSectionVect[i].rMax2);
      z.push_back(fSectionVect[i].z1);
      z.push_back(fSectionVect[i].z2);
    }
  }
  return contour;
}

} /* namespace vecgeom */