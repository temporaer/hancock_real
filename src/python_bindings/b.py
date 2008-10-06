import os
from pyplusplus import module_builder

#Creating an instance of class that will help you to expose your declarations
mb = module_builder.module_builder_t( [r"/home/sarx/prog/boost.python/pypp/src/hello.hpp"]
                                      , gccxml_path=r"" 
                                      , working_directory=r"/home/sarx/prog/boost.python/pypp/build/python"
                                      , include_paths=['/home/sarx/prog/boost.python/pypp/src']
                                      , define_symbols=[] )

#Creating code creator. After this step you should not modify/customize declarations.
mb.build_code_creator( module_name='libhello' )

#Writing code to file.
mb.write_module( './bindings.cpp' )
