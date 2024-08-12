cmake_minimum_required(VERSION 3.12)
project(fjpeg 
LANGUAGES C CXX
VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 11)

# Add the include directory to the project
include_directories(src)

# Add the source file(s) to the project
list(APPEND SOURCE_FILES src/fjpeg.cpp src/fjpeg_transquant.cpp src/fjpeg_huffman.cpp )
list(APPEND SOURCE_FILES_CLI src/fjpeg_cli.cpp)

# Create a static library
if(BUILD_SHARED_LIBS)
  add_library(fjpeg SHARED ${SOURCE_FILES})
else()
  add_library(fjpeg STATIC ${SOURCE_FILES})
  if(MSVC)
    set_target_properties(fjpeg PROPERTIES OUTPUT_NAME libfjpeg)
  endif()  
endif()

if(MSVC)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS) # Disable MSVC warnings
endif()

# Add the source file(s) to the project and link the library
add_executable(fjpeg-cli ${SOURCE_FILES_CLI})

target_link_libraries(fjpeg-cli PUBLIC fjpeg)

# Make the cli binary output name fjpeg
set_target_properties(fjpeg-cli PROPERTIES OUTPUT_NAME fjpeg)
set_target_properties(fjpeg-cli PROPERTIES RUNTIME_OUTPUT_NAME fjpeg)