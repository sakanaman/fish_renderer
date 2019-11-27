add_executable(path_tracing_test main.cpp)
target_include_directories(path_tracing_test
    PRIVATE ${CMAKE_SOURCE_DIR}/lib/Camera
    PRIVATE ${CMAKE_SOURCE_DIR}/lib/DataStracture
    PRIVATE ${CMAKE_SOURCE_DIR}/lib/RNG
    PRIVATE ${CMAKE_SOURCE_DIR}/src/Parallel
    PRIVATE ${CMAKE_SOURCE_DIR}/src/Light
    PRIVATE ${CMAKE_SOURCE_DIR}/ext/stb
    PRIVATE ${CMAKE_SOURCE_DIR}/ext/tinyobjloader
    PRIVATE ${CMAKE_SOURCE_DIR}/lib/Math)
    target_link_libraries(path_tracing_test Camera 
                                            random
                                            IBL 
                                            parallel 
                                            Accelerator 
                                            Vec 
                                            stb_image_write 
                                            tinyobjloader 
                                            stb_image)