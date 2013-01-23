#
# GEANT4 SBT Script to test G4ExtrudedSolid
# I.Hrivnacova, IPN Orsay 28/01/2008 

#
/test/maxPoints 1000
#
# --- extrudedSolid.a1.log
# Extruded solid with triangular polygon
#
/solid/G4ExtrudedSolid 3 (-0.3,0.0,0.3) (-0.3,0.3,-0.3) 2 (-0.3,0.3) (0.0,0.0) (0.0,0.0) (1.0,1.0)
/test/gridSizes 0.1 0.1 0.1 m
/test/errorFileName  log/extrudedSolid.a1.log
/test/run
/voxel/errorFileName log/extrudedSolidv.a1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/extrudedSolid.a2.log
/test/run
#
# --- extrudedSolid.b1.log
# Box defined as Extruded solid
#
/solid/G4ExtrudedSolid 4 (-0.3,-0.3,0.3,0.3) (-0.3,0.3,0.3,-0.3) 2 (-0.3,0.3) (0.0,0.0) (0.0,0.0) (1.0,1.0)
/test/gridSizes 0.1 0.1 0.1 m
/test/errorFileName  log/extrudedSolid.b1.log
/test/run
/voxel/errorFileName log/extrudedSolidv.b1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/extrudedSolid.b2.log
/test/run
#
# --- extrudedSolid.c1.log
# Extruded solid with 4 z-sections
#
/solid/G4ExtrudedSolid 8 (-0.3,-0.3,0.3,0.3,0.15,0.15,-0.15,-0.15) (-0.3,0.3,0.3,-0.3,-0.3,0.15,0.15,-0.3) 4 (-0.4,0.1,0.15,0.4) (-0.2,0.0,0.0,0.2) (0.1,0.0,0.0,0.2) (1.5,0.5,0.7,0.9)
/test/gridSizes 0.1 0.1 0.1 m
/test/errorFileName  log/extrudedSolid.c1.log
/test/run
/voxel/errorFileName log/extrudedSolidv.c1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/extrudedSolid.c2.log
/test/run
#
# --- extrudedSolid.d1.log
# Another extruded solid, where polygon decomposition was failing
# in Geant4 9.1
#
/solid/G4ExtrudedSolid 8 (-0.2,-0.2,0.1,0.1,0.2,0.2,-0.1,-0.1) (0.1,0.25,0.25,-0.1,-0.1,-0.25,-0.25,0.1) 2 (-0.2,0.2) (0.0,0.0) (0.0,0.0) (1.0,1.0)
/test/gridSizes 0.1 0.1 0.1 m
/test/errorFileName  log/extrudedSolid.d1.log
/test/run
/voxel/errorFileName log/extrudedSolidv.d1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/extrudedSolid.d2.log
/test/run
#
exit
