cmake_minimum_required( VERSION 2.6 FATAL_ERROR )
SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/../CMakeModules )
SET(Boost_ADDITIONAL_VERSIONS "1.34" "1.34.1" "1.35" "1.35.1" "1.36" "1.36.1" )

# General stuff
PROJECT(HancockReal CXX Fortran)
SET(CMAKE_PROJECT_NAME "Hancock REAL" )

# find errors/eclipse auto-incl-dir discovery
SET(CMAKE_VERBOSE_MAKEFILE ON)

# Packaging
SET(CPACK_PACKAGE_NAME                hancock_real )
SET(CPACK_PACKAGE_CONTACT             "Hannes Schulz <mail@hannes-schulz.de>" )
SET(CPACK_PACKAGE_VERSION             0.1 )
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Serializes Graphs for Use with REAL" )
SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE i386 )
SET(CPACK_DEBIAN_PACKAGE_PRIORITY     optional )
INCLUDE(CPack)

# No linebreaks in errors (gcc) (easier parsing)
IF(CMAKE_COMPILER_IS_GNUCC)
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0 ")
ENDIF(CMAKE_COMPILER_IS_GNUCC)
IF(CMAKE_COMPILER_IS_GNUCXX)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0 ")
	MESSAGE(STATUS ${CMAKE_CXX_FLAGS})
ENDIF(CMAKE_COMPILER_IS_GNUCXX)


# Precompiled Headers
find_package(PCHSupport QUIET)

# Boost
find_package( Boost 1.35 COMPONENTS program_options )
link_directories ( ${Boost_LIBRARY_DIRS} ${LAPACK_LIBRARIES} )
include_directories ( ${Boost_INCLUDE_DIRS} configuration )

# uBlas Backend
FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)

# Nana
find_package( Nana )
include_directories ( ${NANA_INCLUDE_DIRS} )

# general third party
include_directories ( ${CMAKE_SOURCE_DIR}/third_party )

# SDP
SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${CMAKE_SOURCE_DIR}/third_party/Csdp-6.0.1/include )
SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${CMAKE_SOURCE_DIR}/third_party/Csdp-6.0.1/lib )
find_package( CSDP )
IF(CSDP_FOUND)
	include_directories ( ${CSDP_INCLUDE_DIRS} )
	link_directories    ( ${CSDP_LIBRARY_DIRS} )
	ADD_DEFINITIONS(-DHAVE_CSDP)
ENDIF(CSDP_FOUND)
find_package( DSDP )
IF(DSDP_FOUND)
	include_directories ( ${DSDP_INCLUDE_DIRS} )
	link_directories ( ${DSDP_LIBRARY_DIRS} )
	ADD_DEFINITIONS(-DHAVE_DSDP)
ENDIF(DSDP_FOUND)

# subdirs: sub-libs
ADD_SUBDIRECTORY( seriation_gen )
ADD_SUBDIRECTORY( configuration )
ADD_SUBDIRECTORY( sdp_wrappers )
include_directories( 
  ${CMAKE_SOURCE_DIR}/seriation_gen 
  ${CMAKE_SOURCE_DIR}/configuration 
  ${CMAKE_SOURCE_DIR}/sdp_wrappers )

# subdirs: extensions
ADD_SUBDIRECTORY( python_bindings EXCLUDE_FROM_ALL ) 
ADD_SUBDIRECTORY( tests EXCLUDE_FROM_ALL ) 

ADD_CUSTOM_TARGET( test DEPENDS tests)
ADD_DEPENDENCIES( test tests )

MESSAGE(STATUS "Boost Libs: ${Boost_LIBRARIES} ${Boost_LIBRARY_DIRS}" )
add_executable( hrmain main.cpp )
target_link_libraries(hrmain 
	${Boost_LIBRARIES}  
	${NANA_LIBRARIES} 
	hrconfig 
	seriation_gen 
	sdp_wrappers )

#add_library( hello SHARED hello.cpp )

# generate precompiled header for the library above
#include_directories( ${CMAKE_CURRENT_BINARY_DIR})
#add_precompiled_header( hello ${PROJECT_SOURCE_DIR}/hello.hpp )
