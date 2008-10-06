# Precompiled Headers
find_package(PCHSupport REQUIRED)

# Boost
find_package( Boost 1.34 COMPONENTS python REQUIRED )
link_directories ( ${Boost_LIBRARY_DIRS} )
include_directories ( ${Boost_INCLUDE_DIRS})

# Python-Dev (For Boost.Python)
find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_PATH})

# the library with python exports 
ADD_LIBRARY(python_bindings SHARED bindings.cpp)

# generate precompiled header for the library above
#include_directories( ${CMAKE_CURRENT_BINARY_DIR} )
#add_precompiled_header( python_bindings ${PROJECT_SOURCE_DIR}/hello.hpp )

# link the library above with boost.python
target_link_libraries ( python_bindings ${Boost_LIBRARIES} )

install(TARGETS python_bindings
        RUNTIME DESTINATION bin
		LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib/static )
