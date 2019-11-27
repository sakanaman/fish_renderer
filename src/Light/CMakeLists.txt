add_library(IBL IBL.cpp)
target_include_directories(IBL
    PRIVATE ${CMAKE_SOURCE_DIR}/ext/stb)