SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/../CMakeModules )

# const-char* -> char* conversion, needed for Nana
IF( CMAKE_COMPILER_IS_GNUCXX )
	SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -Wno-write-strings\ -Wall )
ENDIF()

include_directories( 
	. 
	${CMAKE_SOURCE_DIR}/configuration
	${CMAKE_SOURCE_DIR}/seriation_gen
	)

SET(SRC ${SRC}  )
IF(DSDP_FOUND)
  MESSAGE(STATUS "Found DSDP, wrapping it..." )
  SET(SRC ${SRC} dsdp_wrapper.cpp)
ENDIF(DSDP_FOUND)
IF(CSDP_FOUND)
  MESSAGE(STATUS "Found CSDP, wrapping it..." )
  SET(SRC ${SRC} csdp_wrapper.cpp)
ENDIF(CSDP_FOUND)
SET(SRC ${SRC} sdp_wrapper_cfg.cpp sdp_wrapper.cpp )

add_library( sdp_wrappers SHARED ${SRC} )
