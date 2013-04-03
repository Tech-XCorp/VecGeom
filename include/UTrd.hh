/**
* @file
* @author  John Doe <jdoe@example.com>
* @version 0.1
*
* @section LICENSE
*
* @section DESCRIPTION
*
* A simple Orb defined by half-lengths on the three axis. The center of the Orb matches the origin of the local reference frame.

See: http://geant4.web.cern.ch/geant4/G4UsersDocuments/UsersGuides/ForApplicationDeveloper/html/Detector/geomSolids.html
*/

#ifndef USOLIDS_UTrd
#define USOLIDS_UTrd

#include "VUSolid.hh"
#include "UUtils.hh"

class UTrd : public VUSolid
{
	enum ESide {kUndefined,kPX,kMX,kPY,kMY,kPZ,kMZ};
public:
	UTrd() : VUSolid(), fDx1(0), fDx2(0), fDy1(0), fDy2(0), fDz(0) {}
	UTrd(const std::string &pName, double pdx1, double pdx2, double pdy1, double pdy2, double pdz);
	virtual ~UTrd() {}

	// Navigation methods
	EnumInside     Inside (const UVector3 &aPoint) const;   

	virtual double SafetyFromInside ( const UVector3 &aPoint, bool aAccurate=false) const;

	inline double SafetyFromInsideAccurate ( const UVector3 &aPoint) const;

	virtual double SafetyFromOutside( const UVector3 &aPoint, bool aAccurate=false) const;

	inline double SafetyFromOutsideAccurate( const UVector3 &aPoint) const;

	virtual double DistanceToIn     ( const UVector3 &aPoint, 
		const UVector3 &aDirection,
		// UVector3       &aNormalVector,
		double aPstep = UUtils::kInfinity) const;          

	inline double  DistanceToInRoot     ( const UVector3 &aPoint, 
		const UVector3 &aDirection,
		// UVector3       &aNormalVector,
		double aPstep = UUtils::kInfinity) const;          


	virtual double DistanceToOut     ( const UVector3 &aPoint,
		const UVector3 &aDirection,
		UVector3       &aNormalVector, 
		bool           &aConvex,
		double aPstep = UUtils::kInfinity) const;

	inline double DistanceToOutRoot     ( const UVector3 &aPoint,
		const UVector3 &aDirection,
		UVector3       &aNormalVector, 
		bool           &aConvex,
		double aPstep = UUtils::kInfinity) const;

	virtual bool Normal ( const UVector3& aPoint, UVector3 &aNormal ) const; 

	inline bool NormalRoot ( const UVector3& aPoint, UVector3 &aNormal ) const; 
	inline bool NormalGeant4 ( const UVector3& aPoint, UVector3 &aNormal ) const; 

//	virtual void Extent ( EAxisType aAxis, double &aMin, double &aMax ) const;
	void Extent (UVector3 &aMin, UVector3 &aMax) const; 
	virtual double Capacity();
	virtual double SurfaceArea();
	inline VUSolid* Clone() const
	{
		return new UTrd(GetName(), fDx1, fDx2, fDy1, fDy2, fDz);
	}
	virtual UGeometryType GetEntityType() const { return "Trd";}
	virtual void ComputeBBox(UBBox * /*aBox*/, bool /*aStore = false*/) {}

	//G4Visualisation
	virtual void GetParametersList(int /*aNumber*/,double * /*aArray*/) const{} 
	virtual UPolyhedron* GetPolyhedron() const{return CreatePolyhedron();}

	std::ostream& StreamInfo( std::ostream& os ) const;

	UVector3 GetPointOnSurface() const;

	UPolyhedron* CreatePolyhedron () const;

private:  
	inline UVector3 ApproxSurfaceNormal( const UVector3& p ) const;
	double fDx1,fDx2,fDy1,fDy2,fDz;
};
#endif
