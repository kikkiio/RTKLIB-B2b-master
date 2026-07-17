cmake_minimum_required(VERSION 3.10)

project(RTKLIB-B2b VERSION 1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Keep generated binaries and libraries under the repository root.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")

if(WIN32)
    add_compile_definitions(WIN32 _WIN32 WIN32_LEAN_AND_MEAN)
endif()

# Define the common library before applications that link to it.
add_subdirectory(src)
add_subdirectory(app/postppp)
add_subdirectory(app/rtppp)
add_subdirectory(app/postdecoder)
add_subdirectory(app/rtdecoder)
