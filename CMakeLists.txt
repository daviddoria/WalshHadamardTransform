cmake_minimum_required(VERSION 2.6)

Project(WalshHadamardTransform)

include_directories(/media/portable/src/Eigen)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

ADD_EXECUTABLE(Test Test.cpp)
target_link_libraries(Test ${ITK_LIBRARIES})
