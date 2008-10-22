ENABLE_TESTING()

SET(Boost_USE_STATIC_LIBS ON)
find_package( Boost 1.35 COMPONENTS unit_test_framework program_options )

LINK_DIRECTORIES ( ${Boost_LIBRARY_DIRS}   )
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/configuration ${Boost_INCLUDE_DIRS} )

# Test Configuration.hpp
ADD_EXECUTABLE(test_cfg config_test.cpp )
TARGET_LINK_LIBRARIES ( test_cfg ${Boost_LIBRARIES} hrconfig )
ADD_TEST(test_cfg ./test_cfg configuration)

# Test SDP-Prob-Generation
ADD_EXECUTABLE(test_sdp_ser_prob_gen sdp_ser_prob_gen.cpp)
TARGET_LINK_LIBRARIES ( test_sdp_ser_prob_gen ${Boost_LIBRARIES} ${NANA_LIBRARIES} hrconfig seriation_gen )
ADD_TEST(test_sdp_ser_prob_gen ./test_sdp_ser_prob_gen)

# Test SDPA-Wrapper
# Note: Boost has a bug which results in error when child returns with 0. 
# Therefore we _need_ the cmdline param below.
IF(SDPA_FOUND)
	ADD_EXECUTABLE(test_sdpa_wrapper sdpa_wrapper_test.cpp)
	TARGET_LINK_LIBRARIES ( test_sdpa_wrapper ${Boost_LIBRARIES} ${NANA_LIBRARIES} hrconfig sdp_wrappers )
	ADD_TEST(test_sdpa_wrapper ./test_sdpa_wrapper --catch_system_errors=no )
ENDIF(SDPA_FOUND)

# Test SDPLR-Wrapper
# Note: Boost has a bug which results in error when child returns with 0. 
# Therefore we _need_ the cmdline param below.
IF(SDPLR_FOUND)
	ADD_EXECUTABLE(test_sdplr_wrapper sdplr_wrapper_test.cpp)
	TARGET_LINK_LIBRARIES ( test_sdplr_wrapper ${Boost_LIBRARIES} ${NANA_LIBRARIES} hrconfig sdp_wrappers )
	ADD_TEST(test_sdplr_wrapper ./test_sdplr_wrapper --catch_system_errors=no )
ENDIF(SDPLR_FOUND)
# Test DSDP-Wrapper
IF(DSDP_FOUND)
	ADD_EXECUTABLE(test_dsdp_wrapper dsdp_wrapper_test.cpp)
	TARGET_LINK_LIBRARIES ( test_dsdp_wrapper ${Boost_LIBRARIES} ${NANA_LIBRARIES} hrconfig sdp_wrappers )
	ADD_TEST(test_dsdp_wrapper ./test_dsdp_wrapper --catch_system_errors=no )
ENDIF(DSDP_FOUND)
