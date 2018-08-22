// classes for dictionary generation

#include "base/Transformation3D.h"
#include "base/Vector.h"
#include "management/GeoManager.h"

#include "navigation/VNavigator.h"
#include "navigation/VLevelLocator.h"
#include "navigation/VSafetyEstimator.h"
#include "navigation/SimpleLevelLocator.h"

#include "volumes/Box.h"
#include "volumes/Paraboloid.h"
#include "volumes/Parallelepiped.h"
#include "volumes/Sphere.h"
#include "volumes/Trd.h"
#include "volumes/GenTrap.h"
#include "volumes/Trapezoid.h"
#include "volumes/Hype.h"
#include "volumes/Orb.h"
#include "volumes/CutTube.h"
#include "volumes/PlacedAssembly.h"
#include "volumes/Tube.h"
#include "volumes/BooleanVolume.h"
#include "volumes/Cone.h"
#include "volumes/Extruded.h"
#include "volumes/MultiUnion.h"
#include "volumes/Polycone.h"
#include "volumes/Polyhedron.h"
#include "volumes/ScaledShape.h"
#include "volumes/SExtru.h"
#include "volumes/Tessellated.h"
#include "volumes/Tet.h"
#include "volumes/Torus2.h"