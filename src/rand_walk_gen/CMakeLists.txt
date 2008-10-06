SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/../CMakeModules )
#SET(Boost_ADDITIONAL_VERSIONS "1.34" "1.34.1" "1.35" "1.35.1" "1.36" "1.36.1" )

# const-char* -> char* conversion, needed for Nana
IF( CMAKE_COMPILER_IS_GNUCXX )
	SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -Wno-write-strings\ -Wall )
ENDIF()

# Boost
#find_package( Boost 1.35 COMPONENTS ublas )
#link_directories ( ${Boost_LIBRARY_DIRS} )
#include_directories ( ${Boost_INCLUDE_DIRS} . )
include_directories ( . )

# Nana
#find_package( Nana )
#include_directories ( ${NANA_INCLUDE_DIRS} )


add_library( rand_walk_gen SHARED sdp_rand_walk_gen.cpp sdp_rand_walk_prob_gen.cpp )
target_link_libraries(rand_walk_gen  ${Boost_LIBRARIES} ${NANA_LIBRARIES} )
