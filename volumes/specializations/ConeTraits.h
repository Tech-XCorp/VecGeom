/*
 * ConeTraits.h
 *
 *  Created on: May 14, 2014
 *      Author: swenzel
 */

#ifndef VECGEOM_VOLUMES_SPECIALIZATIONS_CONETRAITS_H_
#define VECGEOM_VOLUMES_SPECIALIZATIONS_CONETRAITS_H_

#include <string>
#include <csignal>
#include "volumes/UnplacedCone.h"

namespace VECGEOM_NAMESPACE {

namespace ConeTraits {

#define DEFINE_TRAIT_TYPE(name) \
    struct name { \
      static std::string toString() { \
        return #name; \
      } \
    } \

// A cone that encompasses all cases - not specialized and
// will do extra checks at runtime
DEFINE_TRAIT_TYPE(UniversalCone);

// A cone not having rmin or phi sector
DEFINE_TRAIT_TYPE(NonHollowCone);
// A cone without rmin but with a phi sector smaller than pi
DEFINE_TRAIT_TYPE(NonHollowConeWithSmallerThanPiSector);
// A cone without rmin but with a phi sector greater than pi
DEFINE_TRAIT_TYPE(NonHollowConeWithBiggerThanPiSector);
// A cone without rmin but with a phi sector equal to pi
DEFINE_TRAIT_TYPE(NonHollowConeWithPiSector);

// A cone with rmin and no phi sector
DEFINE_TRAIT_TYPE(HollowCone);
// A cone with rmin and a phi sector smaller than pi
DEFINE_TRAIT_TYPE(HollowConeWithSmallerThanPiSector);
// A cone with rmin and a phi sector greater than pi
DEFINE_TRAIT_TYPE(HollowConeWithBiggerThanPiSector);
// A cone with rmin and a phi sector equal to pi
DEFINE_TRAIT_TYPE(HollowConeWithPiSector);

#undef DEFINE_TRAIT_TYPE

// Mapping of cone types to certain characteristics
enum TreatmentType {
  YES = 0,
  NO,
  UNKNOWN
};


// asking for phi treatment
template <typename T>
struct NeedsPhiTreatment {
  static const TreatmentType value=YES;
};
template <>
struct NeedsPhiTreatment<NonHollowCone> {
  static const TreatmentType value=NO;
};
template <>
struct NeedsPhiTreatment<HollowCone> {
  static const TreatmentType value=NO;
};
template <>
struct NeedsPhiTreatment<UniversalCone> {
  static const TreatmentType value=UNKNOWN;
};


template<typename T>
VECGEOM_INLINE
bool checkPhiTreatment(const UnplacedCone& cone) {
  if(NeedsPhiTreatment<T>::value != UNKNOWN)
    return NeedsPhiTreatment<T>::value == YES;
  else
    // could use a direct constant for 2*M_PI here
    return cone.GetDPhi() < 2.*M_PI;
}

// asking for rmin treatment
template <typename T>
struct NeedsRminTreatment
{
  static const TreatmentType value=YES;
};
template <>
struct NeedsRminTreatment<NonHollowCone>
{
  static const TreatmentType value=NO;
};
template <>
struct NeedsRminTreatment<NonHollowConeWithSmallerThanPiSector>
{
  static const TreatmentType value=NO;
};
template <>
struct NeedsRminTreatment<NonHollowConeWithBiggerThanPiSector>
{
  static const TreatmentType value=NO;
};
template <>
struct NeedsRminTreatment<NonHollowConeWithPiSector>
{
  static const TreatmentType value=NO;
};
template <>
struct NeedsRminTreatment<UniversalCone>
{
  static const TreatmentType value=UNKNOWN;
};


template<typename T>
VECGEOM_INLINE
bool checkRminTreatment(const UnplacedCone& cone) {
  if(NeedsRminTreatment<T>::value != UNKNOWN)
    return NeedsRminTreatment<T>::value == YES;
  else
    return cone.GetRmin1() > 0 || cone.GetRmin2() >0;
}


// sector size
enum AngleType
{
  NOANGLE = 0,
  SMALLER_THAN_PI,
  ONE_PI,
  BIGGER_THAN_PI,
  UNKNOWN_AT_COMPILE_TIME
};

template<typename T>
struct SectorType {
  static const AngleType value=NOANGLE;
};

template<>
struct SectorType<UniversalCone> {
  static const AngleType value=UNKNOWN_AT_COMPILE_TIME;
};

template<>
struct SectorType<NonHollowConeWithSmallerThanPiSector> {
  static const AngleType value=SMALLER_THAN_PI;
};

template<>
struct SectorType<NonHollowConeWithPiSector> {
  static const AngleType value=ONE_PI;
};

template<>
struct SectorType<NonHollowConeWithBiggerThanPiSector> {
  static const AngleType value=BIGGER_THAN_PI;
};
template<>
struct SectorType<HollowConeWithSmallerThanPiSector> {
  static const AngleType value=SMALLER_THAN_PI;
};

template<>
struct SectorType<HollowConeWithPiSector> {
  static const AngleType value=ONE_PI;
};

template<>
struct SectorType<HollowConeWithBiggerThanPiSector> {
  static const AngleType value=BIGGER_THAN_PI;
};

} // end conetrait namespace

} // End global namespace

#endif /* VECGEOM_VOLUMES_SPECIALIZATIONS_CONETRAITS_H_ */
