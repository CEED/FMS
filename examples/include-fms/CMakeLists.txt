cmake_minimum_required(VERSION 3.1...3.18)

project(include-fms)

find_package(fms REQUIRED)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE fms::fms)
