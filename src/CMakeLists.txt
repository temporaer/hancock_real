cmake_minimum_required( VERSION 2.6 FATAL_ERROR )
SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/../CMakeModules )
SET(Boost_ADDITIONAL_VERSIONS "1.34" "1.34.1" "1.35" "1.35.1" "1.36" "1.36.1" )

# General stuff
PROJECT(HancockReal CXX Fortran)
SET(CMAKE_PROJECT_NAME "Hancock REAL" )

# find errors/eclipse auto-incl-dir discovery
#SET(CMAKE_VERBOSE_MAKEFILE ON)

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

IF(CMAKE_COMPILER_IS_GNUCC)
	IF(CMAKE_BUILD_TYPE EQUALS "Debug")
		# Debug mode: Use Safe Containers
		ADD_DEFINITIONS(-D_GLIBCXX_DEBUG)
	ELSE()
		# Other modes: Use Parallel Algorithms
		# turn off at single uses like so:
		# std::sort(v.begin(), v.end(), __gnu_parallel::sequential_tag());
		ADD_DEFINITIONS(-D_GLIBCXX_PARALLEL)
	ENDIF()
ENDIF(CMAKE_COMPILER_IS_GNUCC)


# Precompiled Headers
find_package(PCHSupport QUIET)

# Boost
find_package( Boost 1.35 COMPONENTS program_options )
link_directories ( ${Boost_LIBRARY_DIRS} )
include_directories ( ${Boost_INCLUDE_DIRS} configuration )

# uBlas Backend
FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)
include_directories ( ${LAPACK_INCLUDE_DIRS} )
link_directories    ( ${LAPACK_LIBRARY_DIRS} )
FIND_PATH(BOOST_BINDINGS_PATH boost/numeric/bindings/lapack/lapack.hpp PATHS /usr/local/include/boost-numeric-bindings )
include_directories( ${BOOST_BINDINGS_PATH} )

# Nana
find_package( Nana )
include_directories ( ${NANA_INCLUDE_DIRS} )

# general third party
include_directories ( ${CMAKE_SOURCE_DIR}/third_party )

# SDP-Solvers
SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${CMAKE_SOURCE_DIR}/third_party/Csdp-6.0.1/include )
SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${CMAKE_SOURCE_DIR}/third_party/Csdp-6.0.1/lib )
find_package( CSDP )
IF(CSDP_FOUND)
	include_directories ( ${CSDP_INCLUDE_DIRS} )
	link_directories    ( ${CSDP_LIBRARY_DIRS} )
	FIND_PROGRAM(CSDP_BINARY csdp PATHS ${CMAKE_SOURCE_DIR}/third_party/Csdp-6.0.1/solver/csdp )
ENDIF(CSDP_FOUND)
find_package( DSDP )
IF(DSDP_FOUND)
	include_directories ( ${DSDP_INCLUDE_DIRS} )
	link_directories ( ${DSDP_LIBRARY_DIRS} )
	FIND_PROGRAM(DSDP_BINARY dsdp5)
ENDIF(DSDP_FOUND)
FIND_PROGRAM(SDPA_BINARY sdpa PATHS ${CMAKE_SOURCE_DIR}/third_party/sdpa-7.1.1)
IF(SDPA_BINARY)
	GET_FILENAME_COMPONENT(SDPA_PATH ${SDPA_BINARY} PATH)
	SET(SDPA_FOUND true)
ENDIF(SDPA_BINARY)
FIND_PROGRAM(SDPLR_BINARY sdplr PATHS ${CMAKE_SOURCE_DIR}/third_party/SDPLR-1.02)
IF(SDPLR_BINARY)
	GET_FILENAME_COMPONENT(SDPLR_PATH ${SDPLR_BINARY} PATH)
	SET(SDPLR_FOUND true)
ENDIF(SDPLR_BINARY)

# subdirs: sub-libs
ADD_SUBDIRECTORY( tools )
ADD_SUBDIRECTORY( actions )
ADD_SUBDIRECTORY( adjmat_gen )
ADD_SUBDIRECTORY( seriation_gen )
ADD_SUBDIRECTORY( configuration )
ADD_SUBDIRECTORY( sdp_wrappers )
include_directories( 
  ${CMAKE_SOURCE_DIR}/tools
  ${CMAKE_SOURCE_DIR}/adjmat_gen 
  ${CMAKE_SOURCE_DIR}/actions
  ${CMAKE_SOURCE_DIR}/seriation_gen 
  ${CMAKE_SOURCE_DIR}/configuration 
  ${CMAKE_SOURCE_DIR}/sdp_wrappers )

# subdirs: optional extensions
ADD_SUBDIRECTORY( python_bindings EXCLUDE_FROM_ALL ) 
ADD_SUBDIRECTORY( tests EXCLUDE_FROM_ALL ) 

# docs
find_package(Doxygen QUIET)
if(DOXYGEN_EXECUTABLE)
	CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/../doc/Doxyfile.in
			${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
			IMMEDIATE)
	ADD_CUSTOM_TARGET(api-docs ${DOXYGEN_EXECUTABLE} doc/Doxyfile)
endif(DOXYGEN_EXECUTABLE)



MESSAGE(STATUS "Boost Libs: ${Boost_LIBRARIES} ${Boost_LIBRARY_DIRS}" )
add_executable( hrmain main.cpp )
target_link_libraries(hrmain 
	${Boost_LIBRARIES}  
	${NANA_LIBRARIES} 
	${LAPACK_LIBRARIES} 
	dl
	hrconfig 
	hrtools
	adjmat_gen
	seriation_gen
	)
