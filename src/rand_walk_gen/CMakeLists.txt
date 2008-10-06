SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/../CMakeModules )

# const-char* -> char* conversion, needed for Nana
IF( CMAKE_COMPILER_IS_GNUCXX )
	SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -Wno-write-strings\ -Wall )
ENDIF()

include_directories ( . )

add_library( rand_walk_gen SHARED 
	sdp_rand_walk_gen.cpp 
	sdp_rand_walk_prob_gen.cpp )
target_link_libraries(rand_walk_gen  ${Boost_LIBRARIES} ${NANA_LIBRARIES} )
