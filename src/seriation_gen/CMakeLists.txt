SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/../CMakeModules )

# const-char* -> char* conversion, needed for Nana
IF( CMAKE_COMPILER_IS_GNUCXX )
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings -Wall" )
ENDIF()

include_directories ( 
	. 
	${CMAKE_SOURCE_DIR}/tools
	${CMAKE_SOURCE_DIR}/sdp_wrappers 
	)

add_library( seriation_gen SHARED 
	seriation_gen_cfg.cpp      seriation_gen.hpp
	sdp_seriation_gen.cpp      sdp_seriation_gen.hpp
	sdp_seriation_prob_gen.cpp sdp_seriation_prob_gen.hpp 
	sdp_prob.hpp
	)
target_link_libraries(seriation_gen  ${Boost_LIBRARIES} ${NANA_LIBRARIES} )
