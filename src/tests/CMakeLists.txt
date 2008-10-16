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
