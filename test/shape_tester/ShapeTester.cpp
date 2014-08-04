//
// Implementation of the batch solid  test
//


#include <iomanip>
#include <sstream>
#include <ctime>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "ShapeTester.h"
#include "VUSolid.hh"

#include "base/Vector3D.h"
#include "volumes/Box.h"

#ifdef VECGEOM_ROOT
#include "TGeoShape.h"
#include "TGeoParaboloid.h"
#include "TGeoBBox.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoParaboloid.h"
#include "TGeoVolume.h"
#include "TPolyMarker3D.h"
#include "TRandom3.h"
#include "TColor.h"
#include "TROOT.h"
#include "TAttMarker.h"
#include "TH1D.h"
#endif
    
using namespace std;

ShapeTester::ShapeTester()
{
	SetDefaults();
}

ShapeTester::~ShapeTester()
{

}

void ShapeTester::SetDefaults()
{
	numCheckPoints = 10;
	maxPoints = 10000;
        fVerbose = 1 ;
	repeat = 1000;
	insidePercent = 100.0/3;
	outsidePercent = 100.0/3;
        edgePercent = 0;
        
	outsideMaxRadiusMultiple = 10;
	outsideRandomDirectionPercent = 50;
	differenceTolerance = 0.01;
        ifSaveAllData = true;//false;//true;
        ifMoreTests = true;//false;//true;
        ifDifUSolids = true;
        minDifference = VUSolid::Tolerance();
        difPointsInside = 0;
        difPointsSurface = 0;
        difPointsOutside = 0;

        definedNormal = false;
        ifException = false;
        maxErrorBreak =1000;

	method = "all";
	perftab = perflabels = NULL;
        volumeUSolids = NULL;

  //
  // Zero error list
  //
  errorList = 0;
      
}

UVector3 ShapeTester::GetRandomDirection() 
{
	double phi = 2.*UUtils::kPi*UUtils::Random();
	double theta = UUtils::ACos(1.-2.*UUtils::Random());
	double vx = std::sin(theta)*std::cos(phi);
	double vy = std::sin(theta)*std::sin(phi);
	double vz = std::cos(theta);
	UVector3 vec(vx,vy,vz);
	vec.Normalize();

	return vec;
} 

UVector3 ShapeTester::GetPointOnOrb(double r) 
{
	double phi = 2.*UUtils::kPi*UUtils::Random();
	double theta = UUtils::ACos(1.-2.*UUtils::Random());
	double vx = std::sin(theta)*std::cos(phi);
	double vy = std::sin(theta)*std::sin(phi);
	double vz = std::cos(theta);
	UVector3 vec(vx,vy,vz);
	vec.Normalize();
        vec=vec*r;
	return vec;
} 

// DONE: all set point methods are performance equivalent


void ShapeTester::TestConsistencySolids()
{
    	std::cout<<"% Performing CONSISTENCY TESTS: ConsistencyTests for Inside, Outside and Surface points " <<std::endl;

        TestInsidePoint();
        TestOutsidePoint();
        TestSurfacePoint();

        if (ifSaveAllData){
	  UVector3 point;
	  for (int i = 0; i < maxPoints; i++)
	  {
	   GetVectorUSolids(point, points, i);
           VUSolid::EnumInside inside = volumeUSolids->Inside(point);
           resultDoubleUSolids[i] = (double) inside;
	   }
           SaveResultsToFile("Inside");
        }
 
}
void ShapeTester::ShapeNormal()
{
  int nError =0;
  ClearErrors();
  int i;
  int numTrials =1000;
#ifdef VECGEOM_ROOT
  //Visualisation
   TPolyMarker3D *pm2 = 0;
    pm2 = new TPolyMarker3D();
    pm2->SetMarkerSize(0.2);
    pm2->SetMarkerColor(kBlue);
#endif

  for ( i = 0; i < maxPointsInside; i++)
  {
   UVector3 point = points[i+offsetInside];
   UVector3 dir = directions[i+offsetInside];    
   UVector3 norm;
   bool convex;
      
       VUSolid::EnumInside inside;
       int count = 0;
       double dist=volumeUSolids->DistanceToOut(point, dir ,norm,convex);
       point = point + dist*dir;
       for (int j = 0; j < numTrials; j++)
       { 
         UVector3 dir_new;
         do{
           dir_new=GetRandomDirection();
	   inside = volumeUSolids->Inside(point+dir_new*0.0000001);
	   count++;
         }while((inside!=vecgeom::EInside::kInside)&&(count < 1000));
           
         if(count>=1000)
         {ReportError( &nError,point, dir_new, 0, "SN: Can not find direction pointing Inside after 1000 trials");
            break;
         }
	 count = 0;
         dist=volumeUSolids->DistanceToOut(point, dir_new ,norm,convex);
         if ( dist <  VUSolid::Tolerance() ) {
	   if(inside == vecgeom::EInside::kInside)
           ReportError( &nError,point, dir_new, dist, "SN: DistanceToOut has to be  bigger than tolerance for point Inside");
         }
         double dot=norm.Dot(dir_new);
         if ( dot < 0. )
         {
           ReportError( &nError,point, dir_new, dot, "SN: Wrong direction of Normal calculated by DistanceToOut");
         }
         if(definedNormal)
         {
           UVector3 normal;
           bool valid =volumeUSolids->Normal(point,normal);
           if ( ! valid ) ReportError( &nError,point, dir_new, 0, "SN: Normal has to be valid for point on the Surface");
           dot = normal.Dot(dir_new);
            if ( dot < 0. )
            {
             ReportError( &nError,point, dir_new, dot, "SN: Wrong direction of Normal calculated by Normal");
             }
     
         }
         point = point + dist*dir_new;
#ifdef VECGEOM_ROOT
         //visualisation
         pm2->SetNextPoint(point.x(),point.y(),point.z());
#endif
         if(volumeUSolids->Inside(point)==vecgeom::EInside::kOutside)
         {
           ReportError( &nError,point, dir_new, 0, "SN: DistanceToOut is overshooting,  new point must be on the Surface"); break;
         }
         double safFromIn = volumeUSolids->SafetyFromInside (point);
         double safFromOut = volumeUSolids->SafetyFromOutside (point);
         if ( safFromIn > VUSolid::Tolerance()) ReportError( &nError,point, dir_new, safFromIn, "SN: SafetyFromInside must be zero on Surface ");
         if ( safFromOut > VUSolid::Tolerance()) ReportError( &nError,point, dir_new, safFromOut, "SN: SafetyFromOutside must be zero on Surface");
        
      }
      
  }   
 
#ifdef VECGEOM_ROOT
    //visualisation
    new TCanvas("shape03", "ShapeNormals", 1000, 800);
    pm2->Draw();
#endif
   std::cout<<"% "<<std::endl;    
   std::cout << "% TestShapeNormal reported = " << CountErrors() << " errors"<<std::endl; 
   std::cout<<"% "<<std::endl; 
}

void ShapeTester::ShapeDistances()
{
  int i;
  int nError = 0 ;
  ClearErrors();
  double maxDifOut=0, maxDifIn=0., delta =0.,tolerance=VUSolid::Tolerance();
  bool convex,convex2;
  UVector3 norm;
  UVector3 minExtent,maxExtent;

  volumeUSolids->Extent(minExtent,maxExtent);
  double maxX=std::max(std::fabs(maxExtent.x()),std::fabs(minExtent.x()));
  double maxY=std::max(std::fabs(maxExtent.y()),std::fabs(minExtent.y()));
  double maxZ=std::max(std::fabs(maxExtent.z()),std::fabs(minExtent.z()));
  double maxXYZ=2*std::sqrt(maxX*maxX+maxY*maxY+maxZ*maxZ);
  double dmove = maxXYZ;

 #ifdef VECGEOM_ROOT 
  //Histograms
   TH1D *hist1 = new TH1D("hTest1", "Residual DistancetoIn/Out",200,-20, 0);
     hist1->GetXaxis()->SetTitle("delta[mm] - first bin=overflow");
     hist1->GetYaxis()->SetTitle("count");
     hist1->SetMarkerStyle(kFullCircle);
   TH1D *hist2 = new TH1D("hTest2", "Accuracy distanceToIn for points near Surface",200,-20, 0);
     hist2->GetXaxis()->SetTitle("delta[mm] - first bin=overflow");
     hist2->GetYaxis()->SetTitle("count");
     hist2->SetMarkerStyle(kFullCircle);
  TH1D *hist3 = new TH1D("hTest3", "Accuracy distanceToOut for points near Surface",200,-20, 0);
     hist3->GetXaxis()->SetTitle("delta[mm] - first bin=overflow");
     hist3->GetYaxis()->SetTitle("count");
     hist3->SetMarkerStyle(kFullCircle);
#endif
      
  for ( i = 0; i < maxPointsInside; i++)
  {
   UVector3 point = points[i+offsetInside];
   UVector3 dir = directions[i+offsetInside];    
   double DistanceOut2 = volumeUSolids->DistanceToOut(point, dir ,norm,convex2);
   
   UVector3 pointIn = point+dir*DistanceOut2*(1.-10*tolerance);
   double DistanceOut =  volumeUSolids->DistanceToOut(pointIn, dir ,norm,convex);
   UVector3 pointOut = point+dir*DistanceOut2*(1.+10*tolerance);
   double DistanceIn =  volumeUSolids-> DistanceToIn(pointOut, -dir );
   //Calculate distances for convex or notconvex case 
   double DistanceToInSurf = volumeUSolids->DistanceToIn(point+dir*DistanceOut2,dir);
   if(DistanceToInSurf >= UUtils::kInfinity )
   {
     dmove = maxXYZ;
     if( !convex2 )ReportError( &nError,point, dir, DistanceToInSurf, "SD: Error in convexity, must be convex"); 
   }
   else
   {//reentering solid, it is not convex
     dmove = DistanceToInSurf*0.5;
     if( convex2 )ReportError( &nError,point, dir, DistanceToInSurf, "SD: Error in convexity, must be NOT convex"); 
   }
    double DistanceToIn2 =  volumeUSolids-> DistanceToIn(point+dir*dmove, -dir );
   //double Dif=maxXYZ-DistanceIn-DistanceOut2;
   //std::cout<<"Diff="<<Dif<<std::endl;
   
   if(DistanceOut > 1000.*tolerance)
   ReportError( &nError,pointIn, dir, DistanceOut, "SD: DistanceToOut is not precise");
   if(DistanceIn > 1000.*tolerance)
   ReportError( &nError,pointOut, dir, DistanceIn, "SD: DistanceToIn is not precise "); 

   if( maxDifOut < DistanceOut ) { maxDifOut = DistanceOut;}
   if( ( volumeUSolids->Inside(pointOut-dir*DistanceIn)!=vecgeom::EInside::kOutside )&&(maxDifIn < DistanceIn ))   
   { maxDifIn = DistanceIn;}
   
    double difDelta = dmove-DistanceOut2 - DistanceToIn2 ;
    if(difDelta > 1000.*tolerance)
    ReportError( &nError,point, dir, difDelta, "SD: Distances calculation is not precise");  
    if ( difDelta > delta) delta=std::fabs (difDelta) ; 
    
#ifdef VECGEOM_ROOT
    //Hstograms
    if(std::fabs(difDelta) < 1E-20) difDelta = 1E-30;
    if(std::fabs(DistanceIn) < 1E-20) difDelta = 1E-30;
    if(std::fabs(DistanceOut) < 1E-20) difDelta = 1E-30;
    hist1->Fill(std::max(0.5*std::log(std::fabs(difDelta)),-20.)); 
    hist2->Fill(std::max(0.5*std::log(std::fabs(DistanceIn)),-20.)); 
    hist3->Fill(std::max(0.5*std::log(std::fabs(DistanceOut)),-20.)); 
#endif
   

  }
  if(fVerbose){
   std::cout<<"% TestShapeDistances:: Accuracy max for DistanceToOut="<<maxDifOut<<" from asked accuracy eps="<<10*tolerance<<std::endl;
   std::cout<<"% TestShapeDistances:: Accuracy max for DistanceToIn="<<maxDifIn<<" from asked accuracy eps="<<10*tolerance<<std::endl;
   std::cout<<"% TestShapeDistances:: Accuracy max for Delta="<<delta<<std::endl;
  }
   std::cout<<"% "<<std::endl; 
   std::cout << "% TestShapeDistances reported = " << CountErrors() << " errors"<<std::endl; 
   std::cout<<"% "<<std::endl; 
#ifdef VECGEOM_ROOT
   //Histograms
   TCanvas *c4=new TCanvas("c4", "Residuals DistancsToIn/Out", 800, 600);
   c4->Update();
    hist1->Draw();
    TCanvas *c5=new TCanvas("c5", "Residuals DistancsToIn", 800, 600);
   c5->Update();
    hist2->Draw();
 TCanvas *c6=new TCanvas("c6", "Residuals DistancsToOut", 800, 600);
   c6->Update();
    hist3->Draw();
#endif
  
}

void ShapeTester::TestNormalSolids()
{
   UVector3 point, normal;

   for (int i = 0; i < maxPoints; i++)
   {
     GetVectorUSolids(point, points, i);
     bool valid = volumeUSolids->Normal(point, normal);
     if (ifSaveAllData) 
     {
      resultBoolUSolids[i] = valid;
      SetVectorUSolids(normal, resultVectorUSolids, i);
     }
   }

   SaveResultsToFile("Normal");
}

void ShapeTester::TestSafetyFromOutsideSolids()
{
    std::cout<<"% Performing SAFETYFromOUTSIDE TESTS: ShapeSafetyFromOutside " <<std::endl;
    ShapeSafetyFromOutside(1000);

    if (ifSaveAllData){
	UVector3 point;
	for (int i = 0; i < maxPoints; i++)
	{
           GetVectorUSolids(point, points, i);
	   double res = volumeUSolids->SafetyFromOutside(point,true);
           resultDoubleUSolids[i] = res;
        }
      SaveResultsToFile("SafetyFromOutside");
    }

}

void ShapeTester::TestSafetyFromInsideSolids()
{
    std::cout<<"% Performing SAFETYFromINSIDE TESTS: ShapeSafetyFromInside " <<std::endl;
    ShapeSafetyFromInside(1000);

    if (ifSaveAllData){
	UVector3 point;

	for (int i = 0; i < maxPoints; i++)
	{
 	 GetVectorUSolids(point, points, i);
	 double res = volumeUSolids->SafetyFromInside(point);
         resultDoubleUSolids[i] = res;
	}
  
        SaveResultsToFile("SafetyFromInside");
      }
}

void ShapeTester::PropagatedNormalU(const UVector3 &point, const UVector3 &direction, double distance, UVector3 &normal)
{
  normal.Set(0);
  if (distance < UUtils::kInfinity)
  {
    UVector3 shift = distance * direction;
    UVector3 surfacePoint = point + shift;
    volumeUSolids->Normal(surfacePoint, normal);
    VUSolid::EnumInside e = volumeUSolids->Inside(surfacePoint);
    if (e != vecgeom::EInside::kSurface)
        e = e;
  }
}

void ShapeTester::TestDistanceToInSolids()
{
  std::cout<<"% Performing DISTANCEtoIn TESTS: ShapeDistances, TestsAccuracyDistanceToIn and TestFarAwayPoint " <<std::endl;
  ShapeDistances();           
  TestAccuracyDistanceToIn(1000.);
  TestFarAwayPoint();
 
  if (ifSaveAllData) {  
      UVector3 point, direction;
      for (int i = 0; i < maxPoints; i++)
      {
        GetVectorUSolids(point, points, i);
	GetVectorUSolids(direction, directions, i);
	double res = volumeUSolids->DistanceToIn(point, direction);
        resultDoubleUSolids[i] = res;

        UVector3 normal;
        PropagatedNormalU(point, direction, res, normal);
        SetVectorUSolids(normal, resultVectorUSolids, i);
		
      }
      SaveResultsToFile("DistanceToIn");
}

}

void ShapeTester::TestDistanceToOutSolids()
{
    std::cout<<"% Performing DISTANCEtoOUT TESTS: Shape Normals " <<std::endl;
    ShapeNormal();

    if (ifSaveAllData){

        UVector3 point,normal,direction;
	bool convex;

	for (int i = 0; i < maxPoints; i++)
	{
          GetVectorUSolids(point, points, i);
          GetVectorUSolids(direction, directions, i);
          normal.Set(0);
          double res = volumeUSolids->DistanceToOut(point, direction, normal, convex);

	  resultDoubleUSolids[i] = res;
          resultBoolUSolids[i] = convex;
          SetVectorUSolids(normal, resultVectorUSolids, i);
	}
   }
   SaveResultsToFile("DistanceToOut");
}

void ShapeTester::TestFarAwayPoint()
{
  UVector3 point,point1,vec, direction, normal;
  int icount=0, icount1=0, nError = 0;
  double distIn,diff, difMax=0., maxDistIn =0.;
  double tolerance = VUSolid::Tolerance();
   ClearErrors();
   for ( int j=0; j<maxPointsSurface+maxPointsEdge; j++)
   {
    point = points[j+offsetSurface];
    vec = GetRandomDirection();
    if(volumeUSolids->DistanceToIn(point,vec) < UUtils::kInfinity)continue;
    point1= point;
   
    for (int i=0; i<10000; i++)
    {
          point1 = point1+vec*10000;
    }
    distIn =  volumeUSolids-> DistanceToIn(point1,-vec);
    if( (distIn < UUtils::kInfinity) && (distIn > maxDistIn )) maxDistIn = distIn;
    diff = std::fabs ( (point1 - point). Mag() - distIn );
    if( diff > tolerance ) icount++;
    if( diff >=  UUtils::kInfinity)
    {icount1++;
       ReportError( &nError,point, point1, diff, "TFA:  Point missed Solid (DistanceToIn = Infinity)");
    }
    else{ if ( diff > difMax ) difMax = diff; }
  }
   if(fVerbose)
   {
   std::cout<<"% TestFarAwayPoints:: number of points with big difference (( DistanceToIn- Dist) ) >  tolerance ="<<icount<<std::endl;
   std::cout <<"%  Maxdif = "<<difMax<<" from MaxDist="<<maxDistIn<<" Number of points missing Solid (DistanceToIn = Infinity) = "<<icount1<<std::endl;
   }
   std::cout<<"% "<<std::endl; 
   std::cout << "% TestFarAwayPoints reported = " << CountErrors() << " errors"<<std::endl;
   std::cout<<"% "<<std::endl; 
}

void ShapeTester::TestSurfacePoint()
{
  UVector3 point, pointSurf,vec,direction, normal;
  bool convex;
  int icount=0, icount1=0;
  double distIn,distOut;
  int iIn=0,iInNoSurf=0,iOut=0,iOutNoSurf=0;
  double tolerance = VUSolid::Tolerance();
  int nError=0;
  ClearErrors();
 
   for (int i = 0; i < maxPointsSurface+maxPointsEdge; i++)
     { //test GetPointOnSurface()
       point = points[offsetSurface+i];
       if(volumeUSolids->Inside(point) !=  vecgeom::EInside::kSurface)
       {icount++;   
	 UVector3 v(0,0,0);
         ReportError( &nError,point, v, 0, "TS:  Point on not on the Surface");}
         //test if for point on Surface distIn=distOut=0  
         UVector3 v = GetRandomDirection();
          distIn  = volumeUSolids->DistanceToIn(point,v);
          distOut = volumeUSolids->DistanceToOut(point,v, normal, convex);

          if( distIn == 0. && distOut == 0. )
	    { icount1++;
	      ReportError( &nError,point, v, 0, "TS: DistanceToIn=DistanceToOut=0 for point on Surface");
	  }
        //test Accuracy distance for points near Surface
        pointSurf=point+v*10*tolerance;
        VUSolid::EnumInside inside=volumeUSolids->Inside(pointSurf);
        if(inside != vecgeom::EInside::kSurface)
        {
           if(inside == vecgeom::EInside::kOutside)
           {
            for(int j = 0; j < 1000; j++ )
            {
             vec = GetRandomDirection();
             distIn   = volumeUSolids->DistanceToIn(pointSurf,vec);
	     if(distIn < UUtils::kInfinity)
	    {
              iIn++;
        
              VUSolid::EnumInside surfaceP = volumeUSolids->Inside(pointSurf + distIn*vec);
              if(surfaceP != vecgeom::EInside::kSurface )
	      {
                iInNoSurf++;
	        ReportError( &nError,pointSurf, vec, distIn, "TS: Wrong DistanceToIn for point near Surface");
		        
	      }
            }
	  }
        }
        else
        {
          for(int j = 0; j < 1000; j++ )
          {
            iOut++;
            vec = GetRandomDirection();
            distOut  = volumeUSolids->DistanceToOut(pointSurf, vec,normal, convex); 
            VUSolid::EnumInside surfaceP = volumeUSolids->Inside(pointSurf + distOut*vec);

            if(surfaceP != vecgeom::EInside::kSurface )
	    {
              iOutNoSurf++;
	      ReportError( &nError, pointSurf, vec, distOut, "TS: Wrong DistanceToOut for point near Surface" );
	    }
          }
        }
      

    }
 
  }
   if(fVerbose){
     std::cout<<"% TestSurfacePoints GetPointOnSurface() for Solid  "<<volumeUSolids->GetName()<<" had "<<icount<<" errors"<<std::endl;
     std::cout<<"% TestSurfacePoints both  DistanceToIN and DistanceToOut ==0 for "<<volumeUSolids->GetName()<<" had "<<icount1<<" errors"<<std::endl;
     std::cout<<"% TestSurfacePoints new moved point is not on Surface::iInNoSurf = "<<iInNoSurf<<";    iOutNoSurf = "<<iOutNoSurf<<std::endl;
  
   }
   std::cout<<"% "<<std::endl; 
   std::cout << "% Test Surface Point reported = " << CountErrors() << " errors"<<std::endl;	
   std::cout<<"% "<<std::endl; 
}

void ShapeTester::TestInsidePoint()
{
  int i, n = maxPointsOutside;
  int nError = 0;
  ClearErrors();
 
  UVector3 minExtent,maxExtent;
  volumeUSolids->Extent(minExtent,maxExtent);
  double maxX=std::max(std::fabs(maxExtent.x()),std::fabs(minExtent.x()));
  double maxY=std::max(std::fabs(maxExtent.y()),std::fabs(minExtent.y()));
  double maxZ=std::max(std::fabs(maxExtent.z()),std::fabs(minExtent.z()));
  double maxXYZ=2*std::sqrt(maxX*maxX+maxY*maxY+maxZ*maxZ);

 for (int j = 0; j < maxPointsInside; j++)
 {
    //Check values of Safety
    UVector3 point=points[j+offsetInside];
    double safeDistance = volumeUSolids->SafetyFromInside(point );
    if (safeDistance <= 0.0) {
	UVector3 zero(0);
	ReportError( &nError, point, zero, safeDistance, "TI: SafetyFromInside(p) <= 0");
    return;
    }
    double safeDistanceFromOut = volumeUSolids->SafetyFromOutside(point );
    if (safeDistanceFromOut != 0.0) {
	UVector3 zero(0);
	ReportError(  &nError, point, zero, safeDistance, "TI: SafetyFromOutside(p) not 0 for Point Inside" );
	//continue;
    }
   
    //Check values of Extent
    
    if (point.x() < minExtent.x() ||
	point.x() > maxExtent.x() || 
        point.y() < minExtent.y() || 
        point.y() > maxExtent.y() || 
        point.z() < minExtent.z() || 
	point.z() > maxExtent.z() ) {
	 UVector3 zero(0);
         ReportError(  &nError, point, zero, safeDistance, "TI: Point is outside Extent");
          } 

     //Check values with points and directions to outside points
     for( i=0; i < n; i++ ) 
     {
       UVector3 vr = points[i+offsetOutside] - point;
       UVector3 v = vr.Unit();
       bool valid,convex;
       valid=false;
       UVector3 norm;

       double dist = volumeUSolids->DistanceToOut( point, v, norm,convex);
       double NormalDist ;

       NormalDist = volumeUSolids->SafetyFromInside( point );
     
      if (dist > maxXYZ) {
        ReportError( &nError, point, v, dist, "TI: DistanceToOut(p,v) > Solid's Extent  dist = ");
      continue;
      }
      if (dist <= 0) {
        ReportError( &nError, point, v, NormalDist, "TI: DistanceToOut(p,v) <= 0  Normal Dist = ");
      continue;
      }
      if (dist >= UUtils::kInfinity) {
        ReportError( &nError, point, v, safeDistance, "TI: DistanceToOut(p,v) == kInfinity" );
      continue;
      }
      if (dist < safeDistance-1E-10) {
       ReportError( &nError, point, v, safeDistance, "TI: DistanceToOut(p,v) < DistanceToIn(p)");
      continue;
      }

      if (valid) {
      if (norm.Dot(v) < 0) {
         ReportError( &nError, point, v, safeDistance, "TI: Outgoing normal incorrect" );
	continue;
       }
      }
      //Check DistanceToIn, 0 for now, has to be -1 in future
      double distIn = volumeUSolids->DistanceToIn( point, v);
      if (distIn > 0.) {
	//ReportError( nError, point, v, distIn, "TI: DistanceToIn(p,v) has to be 0 or negative");
	//std::cout<<"distIn="<<distIn<<std::endl;
      continue;
      }
      //Move to the boundary and check
       UVector3 p = point + v*dist;
    
       VUSolid::EnumInside insideOrNot = volumeUSolids->Inside(p);
       if (insideOrNot == vecgeom::EInside::kInside) {
        ReportError( &nError, point, v, safeDistance, "TI: DistanceToOut(p,v) undershoots" );
      continue;
       }
       if (insideOrNot == vecgeom::EInside::kOutside) {
        ReportError( &nError, point, v, safeDistance, "TI: DistanceToOut(p,v) overshoots" );
       continue;
       }
       UVector3 norm1;
        valid = volumeUSolids->Normal( p ,norm1);
        if (norm1.Dot(v) < 0) {
	  if (volumeUSolids->DistanceToIn(p,v) != 0){
           ReportError( &nError, p, v, safeDistance, "TI: SurfaceNormal is incorrect" );
	  }
    }//End Check points and directions
  }
 }
    std::cout<<"% "<<std::endl; 
    std::cout << "% TestInsidePoint reported = " << CountErrors() << " errors"<<std::endl;
    std::cout<<"% "<<std::endl; 
}

void ShapeTester::TestOutsidePoint( )
{
  int i, n = maxPointsInside;
  int nError=0;
  ClearErrors();

  for( int j=0; j < maxPointsOutside; j++ ) {
    //std::cout<<"ConsistencyOutside check"<<j<<std::endl;
  UVector3 point = points[j+offsetOutside];
  double safeDistance = volumeUSolids->SafetyFromOutside( point );
  
  if (safeDistance <= 0.0) {
    UVector3 zero(0);
    ReportError(  &nError, point, zero, safeDistance,"T0: SafetyFromOutside(p) <= 0");
    return;
  }

   double safeDistanceFromInside = volumeUSolids->SafetyFromInside( point );
  
  if (safeDistanceFromInside != 0.0) {
	UVector3 zero(0);
    ReportError(  &nError, point, zero, safeDistance,"T0: SafetyFromInside(p) not 0 for point Outside");
    //continue;
  }
   
   for( i=0; i < n; i++ ) {
    UVector3 vr = points[i+offsetInside] - point;
    UVector3 v = vr.Unit();

    double dist = volumeUSolids->DistanceToIn( point, v );
    if (dist <= 0) {
      ReportError(  &nError, point, v, safeDistance, "T0: DistanceToIn(p,v) <= 0" );
      continue;
    }
    if (dist >= UUtils::kInfinity) {
      ReportError(0, point, v, safeDistance, "T0: DistanceToIn(p,v) == kInfinity" );
      continue;
    }
    if (dist < safeDistance-1E-10) {
      ReportError(  &nError, point, v, safeDistance, "T0: DistanceToIn(p,v) < DistanceToIn(p)" );
      continue;
    }

    UVector3 p = point + dist*v;
    VUSolid::EnumInside insideOrNot = volumeUSolids->Inside( p );
    if (insideOrNot == vecgeom::EInside::kOutside) {
      ReportError(  &nError, point, v, safeDistance, "T0: DistanceToIn(p,v) undershoots");
      continue;
    }
    if (insideOrNot == vecgeom::EInside::kInside) {
      ReportError(  &nError, point, v, safeDistance, "TO: DistanceToIn(p,v) overshoots" );
      continue;
    }

    dist = volumeUSolids->SafetyFromOutside( p );

    //if (dist != 0) {
    if (dist > VUSolid::Tolerance()) {
      ReportError(  &nError, p, v, safeDistance, "T02: DistanceToIn(p) should be zero" );
      // logger << "Dist != 0 : " << dist << endl;
      continue;
    }

    dist = volumeUSolids->SafetyFromInside( p );
    //if (dist != 0) {
    if (dist > VUSolid::Tolerance()) {
      ReportError( 0, p, v, safeDistance, "T02: DistanceToOut(p) should be zero" );
      continue;
    }

    dist = volumeUSolids->DistanceToIn( p, v );
    safeDistance = volumeUSolids->SafetyFromOutside( p );
    //
    // Beware! We might expect dist to be precisely zero, but this may not
    // be true at corners due to roundoff of the calculation of p = point + dist*v.
    // It should, however, *not* be infinity.
    //
    //if (dist != UUtils::kInfinity) {
     if (dist >= UUtils::kInfinity) {
      ReportError(  &nError, p, v, safeDistance, "T02: DistanceToIn(p,v) == kInfinity" );
      continue;
    }	

    bool valid,convex,convex1;
    valid=false;
    UVector3 norm;

    dist = volumeUSolids->DistanceToOut( p, v, norm, convex );
    if (dist == 0) continue;

    if (dist >= UUtils::kInfinity) {
      ReportError(  &nError, p, v, safeDistance, "T02: DistanceToOut(p,v) == kInfinity" );
      continue;
    }
    else if (dist < 0) {
      ReportError(  &nError, p, v, safeDistance, "T02: DistanceToOut(p,v) < 0");
      continue;
    }

    if (valid) {
      if (norm.Dot(v) < 0) {
	ReportError(  &nError, p, v, safeDistance, "T02: Outgoing normal incorrect" );
	continue;
      }
    }

    UVector3 norm1;
    valid = volumeUSolids->Normal( p,norm1 );
    if (norm1.Dot(v) > 0) {
      ReportError(  &nError, p, v, safeDistance, "T02: Ingoing surfaceNormal is incorrect" );
    }


    UVector3 p2 = p + v*dist;

    insideOrNot = volumeUSolids->Inside(p2);
    if (insideOrNot == vecgeom::EInside::kInside) {
      ReportError(  &nError, p, v, safeDistance, "T02: DistanceToOut(p,v) undershoots" );
      continue;
    }
    if (insideOrNot == vecgeom::EInside::kOutside) {
      ReportError(  &nError, p, v, safeDistance, "TO2: DistanceToOut(p,v) overshoots" );
      continue;
    }

    UVector3 norm2, norm3 ;
      valid = volumeUSolids->Normal( p2 , norm2);
    if (norm2.Dot(v) < 0) {
      if (volumeUSolids->DistanceToIn(p2,v) != 0)
	ReportError(  &nError, p2, v, safeDistance, "T02: Outgoing surfaceNormal is incorrect" );
    }
    if (convex) {
      if (norm.Dot(norm2) < 0.0) {
	ReportError(  &nError, p2, v, safeDistance, "T02: SurfaceNormal and DistanceToOut disagree on normal" );
      }
    }

    if (convex) {
      dist = volumeUSolids->DistanceToIn(p2,v);
      if (dist == 0) {
	//
	// We may have grazed a corner, which is a problem of design.
	// Check distance out
	//
	if (volumeUSolids->DistanceToOut(p2,v,norm3,convex1) != 0) {
	  ReportError(  &nError, p, v, safeDistance, "TO2: DistanceToOut incorrectly returns validNorm==true (line of sight)(c)");
	  continue;
	}
      }
      else if (dist != UUtils::kInfinity) {
	//ReportError(  &nError, p, v, safeDistance, "TO2: DistanceToOut incorrectly returns validNorm==true (line of sight)" );
	continue;
      }

      int k;
      //for( k=0; k < n; k++ ) {
        for( k=0; k < 10; k++ ) {
	UVector3 p2top = points[k+offsetInside] - p2;

	if (p2top.Dot(norm) > 0) {
	  ReportError(  &nError, p, v,safeDistance, 
		       "T02: DistanceToOut incorrectly returns validNorm==true (horizon)" );
	  continue;
	}
      }
    } // if valid normal
  } // Loop over inside points
 
   n = maxPointsOutside;

  for(int l=0; l < n; l++ ) {
    UVector3 vr =  points[l+offsetOutside] - point;
    if (vr.Mag2() < DBL_MIN) continue;

    UVector3 v = vr.Unit();

    double dist = volumeUSolids->DistanceToIn( point, v );

    if (dist <= 0) {
      ReportError(  &nError, point, v, safeDistance, "T03: DistanceToIn(p,v) <= 0" );
      continue;
    }
    if (dist >= UUtils::kInfinity) {
      //G4cout << "dist == kInfinity" << G4endl ;
      continue;
    }
    if (dist < safeDistance-1E-10) {
      ReportError(  &nError, point, v, safeDistance, "T03: DistanceToIn(p,v) < DistanceToIn(p)" );
      continue;
    }
    UVector3 p = point + dist*v;

     VUSolid::EnumInside insideOrNot = volumeUSolids->Inside( p );
     if (insideOrNot == vecgeom::EInside::kOutside) {
      ReportError(  &nError, point, v, safeDistance, "T03: DistanceToIn(p,v) undershoots" );
      continue;
    }
     if (insideOrNot == vecgeom::EInside::kInside) {
      ReportError(  &nError, point, v, safeDistance, "TO3: DistanceToIn(p,v) overshoots");
      continue;
    }
  } // Loop over outside points

  }
   std::cout<<"% "<<std::endl; 
   std::cout<< "% TestOutsidePoint reported = " << CountErrors() << " errors"<<std::endl;
   std::cout<<"% "<<std::endl; 
}
//
//Surface Checker 
//
void ShapeTester::TestAccuracyDistanceToIn(double dist)
{
  UVector3 point,pointSurf,v, direction, normal;
  bool convex;
  double distIn,distOut;
  double maxDistIn=0,diff=0,difMax=0;
  int nError=0;
  ClearErrors();
  int  iIn=0,iInNoSurf=0,iOut=0,iOutNoSurf=0,iWrongSideIn=0,iWrongSideOut=0,
        iOutInf=0,iOutZero=0,iInInf=0,iInZero=0,iSafIn=0,iSafOut=0;
  double tolerance = VUSolid::Tolerance();

#ifdef VECGEOM_ROOT
  //Histograms
   TH1D *hist1 = new TH1D("hTest1", "Accuracy DistancetoIn",200,-20, 0);
     hist1->GetXaxis()->SetTitle("delta[mm] - first bin=overflow");
     hist1->GetYaxis()->SetTitle("count");
     hist1->SetMarkerStyle(kFullCircle);
#endif

 //test Accuracy distance 
     for (int i = 0; i < maxPointsSurface+maxPointsEdge; i++)
     { 
      
      //test GetPointOnSurface
      pointSurf = points[i+offsetSurface];
      UVector3 vec = GetRandomDirection();

      point=pointSurf+vec*dist; 

      VUSolid::EnumInside inside=volumeUSolids->Inside(point);
      if(inside !=  vecgeom::EInside::kSurface)
      {
           if(inside ==  vecgeom::EInside::kOutside)
           {
            distIn   = volumeUSolids->DistanceToIn(pointSurf,vec);
            if(distIn >= UUtils::kInfinity){ 
             // Accuracy Test for convex part 
              distIn   = volumeUSolids->DistanceToIn(point,-vec);
              if(maxDistIn < distIn)maxDistIn = distIn;
              diff = ( (pointSurf-point).Mag() - distIn);
              if(diff > difMax) difMax = diff;
               if(std::fabs(diff) < 1E-20) diff = 1E-30;
#ifdef VECGEOM_ROOT
               hist1->Fill(std::max(0.5*std::log(std::fabs(diff)),-20.)); 
#endif
	    }

            for(int j = 0; j < 1000; j++ )
            {
             vec = GetRandomDirection();
             distIn   = volumeUSolids->DistanceToIn(point,vec);
             distOut = volumeUSolids->DistanceToOut(point,vec, normal, convex);
             iWrongSideOut++;
             if(distOut>=UUtils::kInfinity){
              iOutInf++;
              ReportError(  &nError, point, vec, distOut, "TAD: DistanceToOut is Infinity for point Outside");
             }
             if(std::fabs(distOut) < -tolerance){
              iOutZero++;//std::cout<<"WrongSide reply dOut="<<distOut<<std::endl;
              ReportError(  &nError, point, vec, distOut, "TAD: DistanceToOut is < tolerance for point Outside");
             }
             double safFromIn = volumeUSolids->SafetyFromInside(point);
             if(safFromIn > tolerance){
	      iSafOut++;
              ReportError(  &nError, point, vec,safFromIn , "TAD: SafetyFromInside is > tolerance for point Outside");
	     }

	     if(distIn <  UUtils::kInfinity)
	    {
              iIn++;
              
	      //Surface Test 
	      VUSolid::EnumInside surfaceP = volumeUSolids->Inside(point + distIn*vec);
              if(surfaceP !=  vecgeom::EInside::kSurface )
	      {
                iInNoSurf++;
	        ReportError(  &nError, point, vec, 0. , "TAD: Moved to Solid point is not on Surface");
	      }
            }
	  }
        }
        else
        {
          for(int j = 0; j < 1000; j++ )
          {
            iOut++;
            vec = GetRandomDirection();
                   
            distOut  = volumeUSolids->DistanceToOut(point,vec, normal, convex); 
	    VUSolid::EnumInside surfaceP = volumeUSolids->Inside(point + distOut*vec);
             distIn   = volumeUSolids->DistanceToIn(point,vec);
             iWrongSideIn++;
             if(distOut>=UUtils::kInfinity){
              iInInf++;
              ReportError(  &nError, point, vec,distOut , "TAD: Distance ToOut is Infinity  for point Inside");
             }
             if(std::fabs(distOut)<tolerance){
              iInZero++;
              ReportError(  &nError, point, vec,distOut , "TAD: Distance ToOut < tolerance  for point Inside");
             }
             double safFromOut = volumeUSolids->SafetyFromOutside(point);
             if(safFromOut > tolerance){
              iSafIn++;
              ReportError(  &nError, point, vec,safFromOut , "TAD: SafetyFromOutside > tolerance  for point Inside");
	     }
             if(surfaceP != vecgeom::EInside::kSurface )
	     {
              iOutNoSurf++;
	      ReportError(  &nError, point, vec, 0. , "TAD: Moved to Surface point is not on Surface");
	    }
          }
        }
      
     }
    }
   if(fVerbose){
     std::cout<<"TestAccuracyDistanceToIn::Errors  ::iInNoSurf = "<<iInNoSurf<<";    iOutNoSurf = "<<iOutNoSurf<<std::endl;
     std::cout<<"TestAccuracyDistanceToIn::Errors SolidUSolid ::From total number of points iIn = "<<iIn<<";     iOut = "<<iOut<<std::endl;
     std::cout<<"TestAccuracyDistanceToIn::WrongSide Infinity n="<<iInInf<<"+"<<iInZero<<std::endl;
     std::cout<<"TestAccuracyDistanceToIn::WrongSide Zero n="<<iOutInf<<"+"<<iOutZero<<std::endl;
     std::cout<<"TestAccuracyDistanceToIn::WrongSide Saf In and Out  n="<<iSafIn<<"+"<<iSafOut<<std::endl;
     std::cout << "TestForWrongSide:: Point is Outside DistanceToOut not zero = " << iWrongSideOut-iOutZero << "\n";
     std::cout<< "TestForWrongSide:: Point is Outside SafetyFromInside not zero = " << iSafOut << "\n";
     std::cout<< "TestForWrongSide:: Point is Inside DistanceToIn not zero = " << iWrongSideIn-iInZero << "\n";
     std::cout<< "TestForWrongSide:: Point is Inside SafetyFromOutside not zero = " << iSafIn << "\n";
   
   }
#ifdef VECGEOM_ROOT
    TCanvas *c7=new TCanvas("c7", "Accuracy DistancsToIn", 800, 600);
    c7->Update();
    hist1->Draw();
#endif
    std::cout<<"% "<<std::endl; 
    std::cout << "% TestAccuracyDistanceToIn reported = " << CountErrors() << " errors"<<std::endl;
    std::cout<<"% "<<std::endl; 
}

void  ShapeTester::ShapeSafetyFromInside(int max)
{
  UVector3 point,dir,pointSphere,norm;
  bool convex;
  int count=0, count1=0;
  int nError =0;
  ClearErrors();
#ifdef VECGEOM_ROOT
  //visualisation
   TPolyMarker3D *pm3 = 0;
    pm3 = new TPolyMarker3D();
    pm3->SetMarkerSize(0.2);
    pm3->SetMarkerColor(kBlue);
#endif

  if( max > maxPoints )max=maxPoints;
  for (int i = 0; i < max; i++)
  {
   GetVectorUSolids(point, points, i);
   double res = volumeUSolids->SafetyFromInside(point);
   for (int j=0;j<1000;j++)
   { 
     dir=GetRandomDirection();
     pointSphere=point+res*dir;
#ifdef VECGEOM_ROOT
     //visualisation
     pm3->SetNextPoint(pointSphere.x(),pointSphere.y(),pointSphere.z());
#endif
     double distOut=volumeUSolids->DistanceToOut(point,dir,norm,convex);
     if(distOut < res) {count1++;
        ReportError(  &nError, pointSphere, dir, distOut, "SSFI: DistanceToOut is underestimated,  less that Safety" );
     }
     if( volumeUSolids->Inside(pointSphere) == vecgeom::EInside::kOutside)
     { 
       ReportError(  &nError, pointSphere, dir, res, "SSFI: Safety is not safe, point on the SafetySphere is Outside" );
       double error=volumeUSolids->DistanceToIn(pointSphere,-dir);
       if(error>100*VUSolid::Tolerance())
       {   
	count++;
        
       }   
      }
     }
   }
   if(fVerbose){
     std::cout<<"% "<<std::endl; 
     std::cout<<"% ShapeSafetyFromInside ::  number of Points Outside Safety="<<count<<" number of points with  distance smaller that safety="<<count1<<std::endl;
     std::cout<<"% "<<std::endl; 
   }
#ifdef VECGEOM_ROOT
     //visualisation
    new TCanvas("shape", "ShapeSafetyFromInside", 1000, 800);
    pm3->Draw();
#endif
    std::cout<<"% "<<std::endl; 
    std::cout<< "% TestShapeSafetyFromInside reported = " << CountErrors() << " errors"<<std::endl;
    std::cout<<"% "<<std::endl; 
}

void  ShapeTester::ShapeSafetyFromOutside(int max)
{
  UVector3 point,temp,dir,pointSphere,normal;
  double res,error;
  int count=0, count1=0;
  int nError;
  ClearErrors();
#ifdef VECGEOM_ROOT
  //visualisation
   TPolyMarker3D *pm4 = 0;
    pm4 = new TPolyMarker3D();
    pm4->SetMarkerSize(0.2);
    pm4->SetMarkerColor(kBlue);
#endif

  UVector3 minExtent,maxExtent;
  volumeUSolids->Extent(minExtent,maxExtent);
  //double maxX=std::max(std::fabs(maxExtent.x()),std::fabs(minExtent.x()));
  //double maxY=std::max(std::fabs(maxExtent.y()),std::fabs(minExtent.y()));
  //double maxZ=std::max(std::fabs(maxExtent.z()),std::fabs(minExtent.z()));
  //double maxXYZ= std::sqrt(maxX*maxX+maxY*maxY+maxZ*maxZ);
  if( max > maxPointsOutside )max=maxPointsOutside;
  for (int i = 0; i < max; i++)
  {
    //GetVectorUSolids(point, points, i);
    point=points[i+offsetOutside];
    res = volumeUSolids->SafetyFromOutside(point);
    if(res>0)
    {     //Safety Sphere test
     bool convex;
     int numTrials = 1000;
     //if(res > maxXYZ) 
     //{
     // int dummy = (int)(std::pow((maxXYZ/res),2));
     // numTrials = numTrials*dummy;
     //}
     for (int j=0;j<numTrials;j++)
     { dir=GetRandomDirection();
       double distIn=volumeUSolids->DistanceToIn(point,dir);
       if(distIn < res){count1++;
        ReportError(  &nError, point, dir, distIn, "SSFO: DistanceToIn is underestimated,  less that Safety" );
       }
       pointSphere=point+res*dir;
       //std::cout<<"SFO "<<pointSphere<<std::endl;
#ifdef VECGEOM_ROOT
       //visualisation
       pm4->SetNextPoint(pointSphere.x(),pointSphere.y(),pointSphere.z());
#endif
       if( volumeUSolids->Inside(pointSphere) == vecgeom::EInside::kInside)
       { 
             ReportError(  &nError, pointSphere, dir, res, "SSFO: Safety is not safe, point on the SafetySphere is Inside" );
	error=volumeUSolids->DistanceToOut(pointSphere,-dir,normal,convex);
        if(error>100*VUSolid::Tolerance())
        {   
         count++;
	 
        }   
       }
      }
	  
    }                
 }
 if(fVerbose){
  std::cout<<"% "<<std::endl; 
  std::cout<<"% TestShapeSafetyFromOutside::  number of points Inside Safety Sphere ="<<count<<" number of points with Distance smaller that Safety="<<count1<<std::endl;
  std::cout<<"% "<<std::endl; 
 }
#ifdef VECGEOM_ROOT
    //visualisation
    new TCanvas("shapeTest", "ShapeSafetyFromOutside", 1000, 800);
    pm4->Draw();
#endif
   std::cout<<"% "<<std::endl; 
   std::cout<< "% TestShapeSafetyFromOutside reported = " << CountErrors() << " errors"<<std::endl;
   std::cout<<"% "<<std::endl; 
}

void ShapeTester::FlushSS(stringstream &ss)
{
	string s = ss.str();
	cout << s;
	*log << s;
	ss.str("");
}

void ShapeTester::Flush(const string &s)
{
	cout << s;
	*log << s;
}


// NEW: results written normalized to nano seconds per operation
double ShapeTester::NormalizeToNanoseconds(double time)
{
	double res = ((time * (double) 1e+9) / ((double)repeat * (double)maxPoints));
	return res;
}

double ShapeTester::MeasureTest (void (ShapeTester::*funcPtr)(int), const string &amethod)
{
	Flush("Measuring performance of method "+amethod+"\n");

	// NEW: storing data phase, timer is off
	// See: http://www.newty.de/fpt/fpt.html , The Function Pointer Tutorials
	(*this.*funcPtr)(0);

         clock_t start_t, end_t;

        start_t = clock();
	// performance phase, timer is on
	for (int i = 1; i <= repeat; i++)
	{
		(*this.*funcPtr)(i); 
	}

	end_t=clock();
        double realTime = end_t-start_t;

	stringstream ss;
	// NEW: write time per operation / bunch of operations
	ss << "Time elapsed: " << realTime << "s\n";
	ss << "Time per one repeat: " << realTime / repeat << "s\n";
	ss << "Nanoseconds per one method call: " << NormalizeToNanoseconds(realTime) << "\n";
	FlushSS(ss);

	return realTime;
}

void ShapeTester::CreatePointsAndDirectionsSurface()
{
  	UVector3 norm, point;   
	for (int i = 0; i < maxPointsSurface; i++)
	{
	
         UVector3 pointU;
         int retry = 100;
         do 
	   { bool surfaceExist=true;
	   if(surfaceExist) {pointU = volumeUSolids->GetPointOnSurface();}
           else {
                UVector3 dir = GetRandomDirection(), norm;
                bool convex;
                double random=UUtils::Random();
                int index = (int)maxPointsInside*random;
                double dist = volumeUSolids->DistanceToOut(points[index],dir,norm,convex);
                pointU = points[index]+dir*dist ;

           }
           if (retry-- == 0) break;
         }
         while (volumeUSolids->Inside(pointU) != vecgeom::EInside::kSurface);

	UVector3 vec = GetRandomDirection();
	directions[i] = vec;
  	point.Set(pointU.x(), pointU.y(), pointU.z());
        points[i+offsetSurface] = point;
        
	}

}
void ShapeTester::CreatePointsAndDirectionsEdge()
{
	UVector3 norm, point; 
       
	for (int i = 0; i < maxPointsEdge; i++)
	{
	 UVector3 pointU;
         int retry = 100;
         do 
         {
	  volumeUSolids->SamplePointsOnEdge(1,&pointU);
          if (retry-- == 0) break;
         }
         while (volumeUSolids->Inside(pointU) != vecgeom::EInside::kSurface);
 	 UVector3 vec = GetRandomDirection();
	 directions[i] = vec;
   
	point.Set(pointU.x(), pointU.y(), pointU.z());
        points[i+offsetEdge] = point;
        
	}
     
}

void ShapeTester::CreatePointsAndDirectionsOutside()
{

	UVector3 minExtent,maxExtent;
        volumeUSolids->Extent(minExtent,maxExtent);
	double maxX=std::max(std::fabs(maxExtent.x()),std::fabs(minExtent.x()));
	double maxY=std::max(std::fabs(maxExtent.y()),std::fabs(minExtent.y()));
	double maxZ=std::max(std::fabs(maxExtent.z()),std::fabs(minExtent.z()));
        double rOut=std::sqrt(maxX*maxX+maxY*maxY+maxZ*maxZ);
        
        for (int i = 0; i < maxPointsOutside; i++)
	{
	          
	   UVector3 vec, point;
           do
	   {
	    point.x() =  -1 + 2 * UUtils::Random();
	    point.y() = -1 + 2 * UUtils::Random(); 
	    point.z() = -1 + 2 * UUtils::Random();
            point *= rOut * outsideMaxRadiusMultiple;
	   }
	   while (volumeUSolids->Inside(point)!=vecgeom::EInside::kOutside);

  	   double random = UUtils::Random();
	   if (random <= outsideRandomDirectionPercent/100.) 
	   {
		vec = GetRandomDirection();
	   }
	   else
	   {
		UVector3 pointSurface= volumeUSolids->GetPointOnSurface();
		vec = pointSurface - point;
		vec.Normalize();
	   }
		
	   points[i+offsetOutside] = point;
	   directions[i+offsetOutside] = vec;
	}

}

// DONE: inside points generation uses random points inside bounding box
void ShapeTester::CreatePointsAndDirectionsInside()
{       
  UVector3 minExtent,maxExtent;
  volumeUSolids->Extent(minExtent,maxExtent);
  int i = 0; 
  while (i < maxPointsInside)
  {
   double x = RandomRange(minExtent.x(), maxExtent.x());
   double y = RandomRange(minExtent.y(), maxExtent.y());
   if (minExtent.y() == maxExtent.y())
   y = RandomRange(-1000, +1000);
   double z = RandomRange(minExtent.z(), maxExtent.z());
   UVector3 point0(x, y, z);
   if (volumeUSolids->Inside(point0)==vecgeom::EInside::kInside)
   {               
    UVector3 point(x, y, z);
    UVector3 vec = GetRandomDirection();
    points[i+offsetInside] = point;
    directions[i+offsetInside] = vec;
    i++;
   }
  }
}

void ShapeTester::CreatePointsAndDirections()
{ 
  maxPointsInside = (int) (maxPoints * (insidePercent/100));
  maxPointsOutside = (int) (maxPoints * (outsidePercent/100));
  maxPointsEdge = (int) (maxPoints * (edgePercent/100));
  maxPointsSurface = maxPoints - maxPointsInside - maxPointsOutside-maxPointsEdge;
      
  offsetInside = 0;
  offsetSurface = maxPointsInside;
  offsetEdge = offsetSurface + maxPointsSurface;
  offsetOutside = offsetEdge+maxPointsEdge;

  points.resize(maxPoints);
  directions.resize(maxPoints);
  resultDoubleDifference.resize(maxPoints);
  resultBoolUSolids.resize(maxPoints);
  resultDoubleUSolids.resize(maxPoints);

  resultVectorDifference.resize(maxPoints);
  resultVectorUSolids.resize(maxPoints);

  CreatePointsAndDirectionsOutside();
  CreatePointsAndDirectionsInside();
  CreatePointsAndDirectionsSurface();
 
}


#include <sys/types.h>  // For stat().
#include <sys/stat.h>   // For stat().


int directoryExists (string s)
{
  {
   struct stat status;
   stat(s.c_str(), &status);
   return (status.st_mode & S_IFDIR);
  }
  return false;
}


void ShapeTester::PrintCoordinates (stringstream &ss, const UVector3 &vec, const string &delimiter, int precision)
{ 
	ss.precision(precision);
	ss << vec.x() << delimiter << vec.y() << delimiter << vec.z();
}

string ShapeTester::PrintCoordinates (const UVector3 &vec, const string &delimiter, int precision)
{
	static stringstream ss;
	PrintCoordinates(ss, vec, delimiter, precision);
	string res(ss.str());
	ss.str("");
	return res;
}

string ShapeTester::PrintCoordinates (const UVector3 &vec, const char *delimiter, int precision)
{
	string d(delimiter);
	return PrintCoordinates(vec, d, precision);
}

void ShapeTester::PrintCoordinates (stringstream &ss, const UVector3 &vec, const char *delimiter, int precision)
{
	string d(delimiter);
	return PrintCoordinates(ss, vec, d, precision);
}


int ShapeTester::CountDoubleDifferences(const vector<double> &differences, const vector<double> &values1, const vector<double> &values2)	 
{
	int countOfDifferences = 0;
	stringstream ss;

	for (int i = 0; i < maxPoints; i++) 
	{
		double value1 = values1[i];
		double value2 = values2[i];
		double dif = differences[i];
		double difference = std::abs (dif);
		if (difference > std::abs (differenceTolerance*value1))
		{
			if (++countOfDifferences <= 10) ss << "Different point found: index " << i << 
				"; point coordinates:" << PrintCoordinates(points[i], ",") << 
				"; direction coordinates:" << PrintCoordinates(directions[i], ",") <<
				"; difference=" << difference << ")" << 
				"; value2 =" << value2 <<
				"; value1 = " << value1 << "\n";
		}
	}
	ss << "Number of differences is " << countOfDifferences << "\n";
	FlushSS(ss);
	return countOfDifferences;
}

int ShapeTester::CountDoubleDifferences(const vector<double> &differences)
{
	int countOfDifferences = 0;

	for (int i = 0; i < maxPoints; i++) 
	{
		double difference = std::abs (differences[i]);
		if (difference > differenceTolerance) countOfDifferences++;
	}
	stringstream ss;
	ss << "Number of differences is " << countOfDifferences << "\n";
	FlushSS(ss); 
	return countOfDifferences;
}

// NEW: output values precision setprecision (16)
// NEW: for each method, one file

// NEW: print also different point coordinates

void ShapeTester::VectorToDouble(const vector<UVector3> &vectorUVector, vector<double> &vectorDouble)
{
	UVector3 vec;

	int size = vectorUVector.size();
	for (int i = 0; i < size; i++)
	{
		vec = vectorUVector[i];
		double mag = vec.Mag();
		if (mag > 1.1) 
			mag = 1;
		vectorDouble[i] = mag;
	}
}

void ShapeTester::BoolToDouble(const std::vector<bool> &vectorBool, std::vector<double> &vectorDouble)
{
  int size = vectorBool.size();
  for (int i = 0; i < size; i++)
    vectorDouble[i] = (double) vectorBool[i];
}

int ShapeTester::SaveResultsToFile(const string &method1)
{
        string name=volumeUSolids->GetName();
	string filename1(folder+name+"_"+method1+".dat");
	std::cout<<"Saving all results to " <<filename1 <<std::endl;
	ofstream file(filename1.c_str());
	bool saveVectors = (method1 == "Normal");
	int prec = 16;
	if (file.is_open())
	{
		file.precision(prec);
		file << volumeString << "\n";
		string spacer("\t");
		for (int i = 0; i < maxPoints; i++)
		{
			
			file << PrintCoordinates(points[i], spacer, prec) << spacer << PrintCoordinates(directions[i], spacer, prec) << spacer; 
			if (saveVectors) file << PrintCoordinates(resultVectorUSolids[i], spacer, prec) << "\n";
			else file <<  resultDoubleUSolids[i] << "\n";
		}
		return 0;
	}
	std::cout<<"Unable to create file "<<filename1<<std::endl;
	return 1;
}

void ShapeTester::TestMethod(void (ShapeTester::*funcPtr)())
{

        std::cout<< "========================================================= "<<std::endl;
	std::cout<< "% Creating " <<  maxPoints << " points and directions for method =" <<method<<std::endl;

	CreatePointsAndDirections();
           cout.precision(20);

	std::cout<< "% Statistics: points=" << maxPoints << ",\n";

	std::cout << "%             ";
	std::cout << "surface=" << maxPointsSurface << ", inside=" << maxPointsInside << ", outside=" <<   
              maxPointsOutside << "\n";
	std::cout << "%     "<<std::endl;


	(*this.*funcPtr)();
     std::cout<< "========================================================= "<<std::endl;

}

//will run all tests. in this case, one file stream will be used
void ShapeTester::TestMethodAll()
{       
    method = "Consistency";
    TestMethod(&ShapeTester::TestConsistencySolids);
    if(definedNormal)TestMethod(&ShapeTester::TestNormalSolids);
    method = "SafetyFromInside";
    TestMethod(&ShapeTester::TestSafetyFromInsideSolids);
    method = "SafetyFromOutside";
    TestMethod(&ShapeTester::TestSafetyFromOutsideSolids);
    method = "DistanceToIn";
    TestMethod(&ShapeTester::TestDistanceToInSolids);
    method = "DistanceToOut";
    TestMethod(&ShapeTester::TestDistanceToOutSolids);
    method = "all";
}

void ShapeTester::SetFolder(const string &newFolder)
{
   cout << "Checking for existance of " << newFolder << endl;
      
   if (!directoryExists(newFolder))
   {
	string command;
	#ifdef WIN32
		_mkdir(newFolder.c_str());
	#else
		std::cout<<"try to create dir for "<<std::endl;
		mkdir(newFolder.c_str(), 0777);
	#endif
	if (!directoryExists(newFolder))
	{
		cout << "Directory "+newFolder+" does not exist, it must be created first\n";
		exit(1);
	}

        }
	folder = newFolder+"/";

}

void ShapeTester::Run(VUSolid *testVolume)
{
	stringstream ss;

	void (ShapeTester::*funcPtr)()=NULL;

	volumeUSolids= testVolume;
        std::ofstream logger("/log/box");
	log = &logger;
       
        SetFolder("log");
       
	if (method == "") method = "all";
	string name = testVolume->GetName();
	std::cout<< "\n\n";
	std::cout << "===============================================================================\n";
	std::cout << "Invoking test for method " << method << " on " << name << " ..." << "\nFolder is " << folder << std::endl;
	std::cout<< "===============================================================================\n";
	std::cout<< "\n";
	
	if (method == "Consistency") funcPtr = &ShapeTester::TestConsistencySolids;
	if (method == "Normal") funcPtr = &ShapeTester::TestNormalSolids;
	if (method == "SafetyFromInside") funcPtr = &ShapeTester::TestSafetyFromInsideSolids;
	if (method == "SafetyFromOutside") funcPtr = &ShapeTester::TestSafetyFromOutsideSolids;
	if (method == "DistanceToIn") funcPtr = &ShapeTester::TestDistanceToInSolids;
	if (method == "DistanceToOut") funcPtr = &ShapeTester::TestDistanceToOutSolids;

	if (method == "all") TestMethodAll();
	else if (funcPtr) TestMethod(funcPtr);
	else std::cout<< "Method " << method << " is not supported" << std::endl;

        ClearErrors();
        method = "all";
}
void ShapeTester::RunMethod(VUSolid *testVolume, std::string method1)
{
	stringstream ss;

	void (ShapeTester::*funcPtr)()=NULL;

	volumeUSolids= testVolume;
        std::ofstream logger("/log/box");
	log = &logger;
       
        SetFolder("log");
 
        method = method1;
 
	if (method == "") method = "all";
	string name = testVolume->GetName();
        
	std::cout<< "\n\n";
	std::cout << "===============================================================================\n";
	std::cout << "Invoking test for method " << method << " on " << name << " ..." << "\nFolder is " << folder << std::endl;
	std::cout<< "===============================================================================\n";
	std::cout<< "\n";
	
	if (method == "Consistency") funcPtr = &ShapeTester::TestConsistencySolids;
	if (method == "Normal") funcPtr = &ShapeTester::TestNormalSolids;
	if (method == "SafetyFromInside") funcPtr = &ShapeTester::TestSafetyFromInsideSolids;
	if (method == "SafetyFromOutside") funcPtr = &ShapeTester::TestSafetyFromOutsideSolids;
	if (method == "DistanceToIn") funcPtr = &ShapeTester::TestDistanceToInSolids;
	if (method == "DistanceToOut") funcPtr = &ShapeTester::TestDistanceToOutSolids;

	if (method == "all") TestMethodAll();
	else if (funcPtr) TestMethod(funcPtr);
	else std::cout<< "Method " << method << " is not supported" << std::endl;

	ClearErrors();
        method = "all";
}
//
// ReportError
//
// Report the specified error message, but only if it has not been reported a zillion
// times already.
//
void ShapeTester::ReportError( int *nError,  UVector3 &p, 
			   UVector3 &v, double distance,
			       std::string comment)//, std::ostream &logger )
{
  
  ShapeTesterErrorList *last=0, *errors = errorList;
  while( errors ) {
    
    if (errors->message == comment) {
      if ( ++errors->nUsed > 5 ) return;
      break;
    }
    last = errors;
    errors = errors->next;
  }

  if (errors == 0) {
    //
    // New error: add it the end of our list
    //
    errors = new ShapeTesterErrorList;
    errors->message = comment;
    errors->nUsed = 1;
    errors->next = 0;
    if (errorList) 
      last->next = errors;
    else
      errorList = errors;
  }

  //
  // Output the message
  //	
 
  std::cout << "% " << comment;
  if (errors->nUsed == 5) std::cout << " (any further such errors suppressed)";
  std::cout << " Distance = " << distance ; 
  std::cout << std::endl;

  std::cout << ++(*nError) << " " << p.x() << " " << p.y() << " " << p.z() 
	    << " " << v.x() << " " << v.y() << " " << v.z() << std::endl;
  
  //
  // if debugging mode we have to exit now
  //
  if(ifException){
     std::ostringstream text;
     text << "Abborting due to Debugging mode in solid: " << volumeUSolids->GetName();
     UUtils::Exception("ShapeTester", "Debugging mode", FatalErrorInArguments, 1, text.str().c_str());
  }
}
//
// ClearErrors
// Reset list of errors (and clear memory)
//
void ShapeTester::ClearErrors()
{
  ShapeTesterErrorList *here, *next;

  here = errorList;
  while( here ) {
    next = here->next;
    delete here;
    here = next;
  }
  errorList = 0;
}
//
// CountErrors
//
int ShapeTester::CountErrors() const
{
  ShapeTesterErrorList *here;
  int answer = 0;

  here = errorList;
  while( here ) {
    answer += here->nUsed;
    here = here->next;
  }

  return answer;
}
