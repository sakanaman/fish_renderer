cmake_minimum_required(VERSION 3.15.2)
project(fish_renderer)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wextra -ferror-limit=100")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3")


add_subdirectory(ext)

add_subdirectory(lib)

add_subdirectory(src)

add_subdirectory(projects)

add_subdirectory(test)