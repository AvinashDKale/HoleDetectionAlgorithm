# CMakeList.txt : CMake project for HoleDetectionAlgorithmProject, include source and define

cmake_minimum_required (VERSION 3.8)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_BUILD_TYPE Debug)

#set(CMAKE_DEBUG_POSTFIX ON)

set(CMAKE_TOOLCHAIN_FILE $ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake)

project ("HoleDetectionAlgorithmProject" C CXX)



find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)
include_directories(${CGAL_INCLUDE_DIRS})

link_directories(${CGAL_LIBRARY_DIRS})

add_definitions(${CGAL_DEFINITIONS})


find_package(PCL  REQUIRED)
find_package(OpenCV CONFIG REQUIRED )


include_directories(${PCL_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})

link_directories(${PCL_LIBRARY_DIRS})

add_definitions(${PCL_DEFINITIONS})

include_directories(source include)
#ShapeReconstructionPointCloud.cpp
set(SRC main.cpp src/PointCloudPreprocessing.cpp src/FindHole.cpp)
#ShapeReconstructionPointCloud.h 
set(HDR  include/PointCloudPreprocessing.h include/FindHole.h)


# Add source to this project's executable.
add_executable(${PROJECT_NAME} ${SRC} ${HDR})
#target_link_libraries("HoleDetectionAlgorithmProject" PUBLIC CGAL::CGAL_Basic_viewer )
#target_link_libraries("HoleDetectionAlgorithmProject" PRIVATE ${PCL_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${PCL_LIBRARIES} ${CGAL_LIBRARIES} ${OpenCV_LIBS} opencv_ml opencv_dnn opencv_core opencv_flann )
# TODO: Add tests and install targets if needed.
