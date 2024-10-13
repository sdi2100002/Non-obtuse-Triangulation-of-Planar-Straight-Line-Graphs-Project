# CMake version requirement
cmake_minimum_required(VERSION 3.1...3.23)

# Project name and language
project(test LANGUAGES CXX)

# Find and include CGAL
set(CGAL_DIR ../lib/CGAL-5.6.1/build)

# Find and include CGAL with Qt5 support
find_package(CGAL QUIET COMPONENTS Qt5)
if (NOT CGAL_FOUND)
    message(STATUS "CGAL or CGAL Qt5 not found. Compilation will not proceed.")
    return()
endif()

# Find Boost
find_package(Boost REQUIRED COMPONENTS system filesystem json)
if (NOT Boost_FOUND)
    message(STATUS "Boost not found. Compilation will not proceed.")
    return()
endif()

# Find Qt5 components (including Charts)
find_package(Qt5 COMPONENTS Widgets OpenGL Charts REQUIRED)
if (NOT Qt5_FOUND)
    message(STATUS "Qt5 not found. Compilation will not proceed.")
    return()
endif()

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# Source files
set(SOURCES
    firstp.cpp
    ./triangulation/triangulation.cpp
    # ./utils/utils.cpp
    ./graphics/graphics.cpp
    # ./utilsAlgo/utilsAlgo.cpp
)

# Create the executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Link the executable to CGAL and third-party libraries
target_link_libraries(${PROJECT_NAME} PRIVATE CGAL::CGAL CGAL::CGAL_Qt5)

# Link Qt5 libraries
target_link_libraries(${PROJECT_NAME} PRIVATE Qt5::Widgets Qt5::OpenGL Qt5::Charts)

# Enable CGAL viewer functionality
add_definitions(-DCGAL_USE_BASIC_VIEWER)
message(STATUS "Qt5 found, enabling CGAL Qt5 viewer.")

