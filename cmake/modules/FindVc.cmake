# - Finds Vc installation ( the wrapper library to SIMD intrinsics )
# This module sets up Vc information 
# It defines:
# VC_FOUND          If the ROOT is found
# VC_INCLUDE_DIR    PATH to the include directory
# VC_LIBRARIES      Most common libraries
# VC_LIBRARY_DIR    PATH to the library directory 

# look if an environment variable VCROOT exists
set(VCROOT $ENV{VCROOT})

if( NOT VCROOT )
       	set(VC_FOUND FALSE)	
	message(STATUS "Vc library not found; try to set a VCROOT environment variable to the base installation path or add -DVCROOT= to the cmake command")	
else()	
	find_library(VC_LIBRARIES libVc.a PATHS ${VCROOT}/lib )
	if( VC_LIBRARIES  )	
		set(VC_FOUND TRUE)	
		set(VC_INCLUDE_DIR ${VCROOT}/include)
		message(STATUS "Found Vc library in ${VC_LIBRARIES}")		
	else()
		message(STATUS "No Vc library found in ${VCROOT}/lib")		
	endif()
endif()


