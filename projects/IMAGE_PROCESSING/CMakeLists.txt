add_executable(environmentmap2cdf environmentmap2cdf.cpp)
target_include_directories(environmentmap2cdf
    PRIVATE ${CMAKE_SOURCE_DIR}/ext/stb
    PRIVATE ${CMAKE_SOURCE_DIR}/lib/RNG)
target_link_libraries(environmentmap2cdf stb_image stb_image_write random)