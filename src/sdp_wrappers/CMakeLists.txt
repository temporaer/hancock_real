SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/../CMakeModules )
SET(Boost_ADDITIONAL_VERSIONS "1.34" "1.34.1" "1.35" "1.35.1" "1.36" "1.36.1" )

# const-char* -> char* conversion, needed for Nana
IF( CMAKE_COMPILER_IS_GNUCXX )
	SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -Wno-write-strings\ -Wall )
ENDIF()

# Boost
find_package( Boost 1.35 COMPONENTS ublas )
link_directories ( ${Boost_LIBRARY_DIRS} )
include_directories ( ${Boost_INCLUDE_DIRS} . )

# Nana
find_package( Nana )
include_directories ( ${NANA_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/rand_walk_gen )
link_directories ( ${NANA_LIBRARY_DIRS} )

# DSDP
find_package( DSDP )
include_directories ( ${DSDP_INCLUDE_DIRS} )
link_directories ( ${DSDP_LIBRARY_DIRS} )

# CSDP
find_package( CSDP )
include_directories ( ${CSDP_INCLUDE_DIRS} )
link_directories    ( ${CSDP_LIBRARY_DIRS} )

SET(SRC  )
IF(DSDP_FOUND)
  MESSAGE(STATUS "Found DSDP, wrapping it..." )
  SET(SRC ${SRC} dsdp_wrapper.cpp)
ENDIF()
IF(CSDP_FOUND)
  MESSAGE(STATUS "Found CSDP, wrapping it..." )
  SET(SRC ${SRC} csdp_wrapper.cpp)
ENDIF()

add_library( sdp_wrappers SHARED ${SRC} )
target_link_libraries(sdp_wrappers  ${Boost_LIBRARIES} ${NANA_LIBRARIES} ${DSDP_LIBRARIES} ${CSDP_LIBRARIES} )
