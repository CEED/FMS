cmake_minimum_required(VERSION 3.1)

project(include-fms)

find_package(fms REQUIRED
  HINTS ${FMS_DIR})

add_executable(main main.cpp)
target_link_libraries(main PRIVATE fms::fms)
