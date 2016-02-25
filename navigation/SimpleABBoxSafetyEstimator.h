/*
 * SimpleABBoxSafetyEstimator.h
 *
 *  Created on: Aug 28, 2015
 *      Author: swenzel
 */

#ifndef NAVIGATION_SIMPLEABBOXSAFETYESTIMATOR_H_
#define NAVIGATION_SIMPLEABBOXSAFETYESTIMATOR_H_

#include "navigation/VSafetyEstimator.h"
#include "management/ABBoxManager.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

//! a safety estimator using a (vectorized) search through bounding boxes to exclude certain daughter volumes
//! to talk to
class SimpleABBoxSafetyEstimator : public VSafetyEstimatorHelper<SimpleABBoxSafetyEstimator> {

private:
  // we keep a reference to the ABBoxManager ( avoids calling Instance() on this guy all the time )
  ABBoxManager &fABBoxManager;

  SimpleABBoxSafetyEstimator() : VSafetyEstimatorHelper<SimpleABBoxSafetyEstimator>(), fABBoxManager(ABBoxManager::Instance()) {}

  // convert index to physical daugher
  VPlacedVolume const *LookupDaughter(LogicalVolume const *lvol, int id) const {
    assert(id >= 0 && "access with negative index");
    assert(size_t(id) < lvol->GetDaughtersp()->size() && "access beyond size of daughterlist ");
    return lvol->GetDaughtersp()->operator[](id);
  }

  // helper function calculating some candidate volumes
  size_t GetSafetyCandidates_v(Vector3D<Precision> const &point, ABBoxManager::ABBoxContainer_v const &corners,
                               size_t size, ABBoxManager::BoxIdDistancePair_t *boxsafetypairs,
                               Precision upper_squared_limit) const {
    size_t count = 0;
#ifdef VECGEOM_VC
    Vector3D<float> pointfloat((float)point.x(), (float)point.y(), (float)point.z());
    size_t vecsize = size;
    for (size_t box = 0; box < vecsize; ++box) {
      ABBoxManager::Real_v safetytoboxsqr = ABBoxImplementation::ABBoxSafetySqr<kVcFloat, ABBoxManager::Real_t>(
          corners[2 * box], corners[2 * box + 1], pointfloat);

      ABBoxManager::Bool_v hit = safetytoboxsqr < ABBoxManager::Real_t(upper_squared_limit);
      if (Any(hit)) {
        for (size_t i = 0; i < kVcFloat::precision_v::Size; ++i) {
          if (hit[i]) {
            boxsafetypairs[count] =
                ABBoxManager::BoxIdDistancePair_t(box * kVcFloat::precision_v::Size + i, safetytoboxsqr[i]);
            count++;
          }
        }
      }
    }
#else
#pragma message("implementation for GetSafetyCandidates for scalar backend is missing")
#endif
    return count;
  }

private:
  VECGEOM_INLINE
  Precision TreatSafetyToIn(Vector3D<Precision> const &localpoint, VPlacedVolume const *pvol, Precision outsafety ) const {
    // a stack based workspace array
    static __thread ABBoxManager::BoxIdDistancePair_t boxsafetylist[VECGEOM_MAXDAUGHTERS] = {};

    double safety = outsafety; // we use the outsafety estimate as starting point
    double safetysqr = safety * safety;

    // safety to bounding boxes
    LogicalVolume const *lvol = pvol->GetLogicalVolume();
    if (safety > 0. && lvol->GetDaughtersp()->size() > 0) {
      int size;

      ABBoxManager::ABBoxContainer_v bboxes = fABBoxManager.GetABBoxes_v(lvol, size);
      // calculate squared bounding box safeties in vectorized way
      auto ncandidates = GetSafetyCandidates_v(localpoint, bboxes, size, boxsafetylist, safetysqr);
      // not sorting the candidate list ( which one could do )
      for (unsigned int candidate = 0; candidate < ncandidates; ++candidate) {
        auto boxsafetypair = boxsafetylist[candidate];
        if (boxsafetypair.second < safetysqr) {
          VPlacedVolume const *candidate = LookupDaughter(lvol, boxsafetypair.first);
          if (boxsafetypair.first > lvol->GetDaughtersp()->size())
            break;
          auto candidatesafety = candidate->SafetyToIn(localpoint);
#ifdef VERBOSE
          if (candidatesafety * candidatesafety > boxsafetypair.second && boxsafetypair.second > 0)
            std::cerr << "real safety smaller than boxsafety \n";
#endif
          if (candidatesafety < safety) {
            safety = candidatesafety;
            safetysqr = safety * safety;
          }
        }
      }
    }
    return safety;
  }

public:
  static constexpr const char *gClassNameString = "SimpleABBoxSafetyEstimator";

  VECGEOM_INLINE
  virtual Precision ComputeSafetyForLocalPoint(Vector3D<Precision> const &localpoint,
                                               VPlacedVolume const *pvol) const override {


    // safety to mother
    double safety = pvol->SafetyToOut(localpoint);
    return TreatSafetyToIn(localpoint,pvol,safety);
  }

#ifdef VECGEOM_BACKEND_PRECISION_NOT_SCALAR
  VECGEOM_INLINE
   virtual VECGEOM_BACKEND_PRECISION_TYPE ComputeSafetyForLocalPoint(Vector3D<VECGEOM_BACKEND_PRECISION_TYPE> const &localpoint,
                                                VPlacedVolume const *pvol, VECGEOM_BACKEND_PRECISION_TYPE::Mask m) const override {
     VECGEOM_BACKEND_PRECISION_TYPE safety(0.);
     if (Any(m)) {
       // SIMD safety to mother
       auto safety = pvol->SafetyToOut(localpoint);

       // now loop over the voxelized treatment of safety to in
       for (unsigned int i = 0; i < VECGEOM_BACKEND_PRECISION_TYPE::Size; ++i) {
         if (m[i]) {
           safety[i] = TreatSafetyToIn(Vector3D<Precision>(localpoint.x()[i], localpoint.y()[i], localpoint.z()[i]),
                                       pvol, safety[i]);
         }
         else{
           safety[i] = 0;
         }
       }
     }
     return safety;
  }
#endif


  // vector interface
  VECGEOM_INLINE
  virtual void ComputeSafetyForLocalPoints(SOA3D<Precision> const &localpoints, VPlacedVolume const *pvol,
                                           Precision *safeties) const override {
    // a stack based workspace array
    static __thread ABBoxManager::BoxIdDistancePair_t boxsafetylist[VECGEOM_MAXDAUGHTERS] = {};

    // safety to mother -- using vector interface
    pvol->SafetyToOut(localpoints, safeties);

    // safety to bounding boxes
    LogicalVolume const *lvol = pvol->GetLogicalVolume();
    if (!(lvol->GetDaughtersp()->size() > 0))
      return;

    // get bounding boxes (they are the same for all tracks)
    int numberofboxes;
    auto bboxes = fABBoxManager.GetABBoxes_v(lvol, numberofboxes);

    // now loop over particles
    for (int i = 0, ntracks = localpoints.size(); i < ntracks; ++i) {
      double safety = safeties[i];
      if (safeties[i] > 0.) {
        double safetysqr = safeties[i] * safeties[i];
        auto lpoint = localpoints[i];
        // vectorized search through bounding boxes -- quickly excluding many candidates
        auto ncandidates = GetSafetyCandidates_v(lpoint, bboxes, numberofboxes, boxsafetylist, safetysqr);
        // loop over remaining candidates
        for (unsigned int candidate = 0; candidate < ncandidates; ++candidate) {
          auto boxsafetypair = boxsafetylist[candidate];
          if (boxsafetypair.second < safetysqr) {
            VPlacedVolume const *candidate = LookupDaughter(lvol, boxsafetypair.first);
            if (boxsafetypair.first > lvol->GetDaughtersp()->size())
              break;
            auto candidatesafety = candidate->SafetyToIn(lpoint);
#ifdef VERBOSE
            if (candidatesafety * candidatesafety > boxsafetypair.second && boxsafetypair.second > 0)
              std::cerr << "real safety smaller than boxsafety \n";
#endif
            if (candidatesafety < safety) {
              safety = candidatesafety;
              safetysqr = safety * safety;
            }
          }
        }
      }
      // write back result
      safeties[i] = safety;
    }
  }

  static VSafetyEstimator *Instance() {
    static SimpleABBoxSafetyEstimator instance;
    return &instance;
  }

}; // end class

}} // end namespace


#endif /* NAVIGATION_SIMPLEABBOXSAFETYESTIMATOR_H_ */
