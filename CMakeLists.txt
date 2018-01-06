cmake_minimum_required(VERSION 3.5)
project("calabash")
set(CMAKE_CXX_STANDARD 14)

find_package(Eigen3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIRS}
    ${CMAKE_CURRENT_LIST_DIR}/calabash)

if(DEBUG_DIAMOND)
    add_definitions(-D_DEBUG_DIAMOND_)
endif()

file(GLOB src ${CMAKE_CURRENT_LIST_DIR}/calabash/*.cpp)

add_executable(${PROJECT_NAME} ${src})

