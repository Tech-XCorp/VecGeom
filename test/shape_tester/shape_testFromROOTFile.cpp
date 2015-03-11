#include "ShapeTester.h"
#include "VUSolid.hh"
#include "management/RootGeoManager.h"
#include "volumes/LogicalVolume.h"
#include "management/GeoManager.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "base/Vector3D.h"

#ifdef VECGEOM_ROOT
#include "TApplication.h"
#endif
#include "stdlib.h"

using namespace vecgeom;


// benchmarking any available shape (logical volume) found in a ROOT file
// usage: BenchmarkShapeFromROOTFile detector.root logicalvolumename
// logicalvolumename should not contain trailing pointer information

int main(  int argc,char *argv[]) {

  if( argc < 3 )
  {
    std::cerr << "need to give root geometry file and logical volume name";
  }

  TGeoManager::Import( argv[1] );
  std::string testvolume( argv[2] );

  int found = 0;
  TGeoVolume * foundvolume = NULL;
  // now try to find shape with logical volume name given on the command line
  TObjArray *vlist = gGeoManager->GetListOfVolumes( );
  for( auto i = 0; i < vlist->GetEntries(); ++i )
  {
    TGeoVolume * vol = reinterpret_cast<TGeoVolume*>(vlist->At( i ));
    std::string fullname(vol->GetName());
    
    // strip off pointer information
    std::string strippedname(fullname, 0, fullname.length()-4);
    
    std::size_t founds = strippedname.compare(testvolume);
    if (founds == 0){
      found++;
      foundvolume = vol;
      
      std::cerr << "found matching volume " << fullname 
		<< " of type " << vol->GetShape()->ClassName() << "\n";
    }
  }

  std::cerr << "volume found " << found << " times \n";
  foundvolume->GetShape()->InspectShape();
  std::cerr << "volume capacity " 
	    << foundvolume->GetShape()->Capacity() << "\n";


  // now get the VecGeom shape and benchmark it
  if( foundvolume )
  {
    LogicalVolume * converted = RootGeoManager::Instance().Convert( foundvolume );
    VUSolid* shape= converted->Place()->ConvertToUSolids()->Clone();
    ShapeTester tester;

    std::cout << "\n==============Shape StreamInfo ========= \n";
    shape->StreamInfo(std::cout);

    if(argc>3)
    {
      if(strcmp(argv[3],"vis")==0)
      {
      #ifdef VECGEOM_ROOT
	TApplication theApp("App",0,0);
	tester.Run(shape);
	theApp.Run();
      #endif
      }
    }
    else
    {
      tester.Run(shape);
    }
    return 1;
  }
  else
  {
    std::cerr << " NO SUCH VOLUME [" 
	      <<  testvolume << "] FOUND ... EXITING \n";
    return 1;
  }
}
