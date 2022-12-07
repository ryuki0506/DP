cmake_minimum_required(VERSION 3.2)
project(DP CXX)

add_subdirectory(src)

add_executable(main main.cpp)
target_link_libraries(main greetings)
