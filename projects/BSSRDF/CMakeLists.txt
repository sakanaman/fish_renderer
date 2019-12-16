add_executable(local_onb test_local_onb.cpp)
target_include_directories(local_onb
    PRIVATE ${CMAKE_SOURCE_DIR}/lib/Camera
    PRIVATE ${CMAKE_SOURCE_DIR}/lib/DataStracture
    PRIVATE ${CMAKE_SOURCE_DIR}/lib/RNG
    PRIVATE ${CMAKE_SOURCE_DIR}/src/Parallel
    PRIVATE ${CMAKE_SOURCE_DIR}/src/Light
    PRIVATE ${CMAKE_SOURCE_DIR}/ext/stb
    PRIVATE ${CMAKE_SOURCE_DIR}/ext/tinyobjloader
    PRIVATE ${CMAKE_SOURCE_DIR}/lib/Math)
    target_link_libraries(local_onb         Camera 
                                            random
                                            IBL 
                                            parallel 
                                            Accelerator 
                                            Vec 
                                            stb_image_write 
                                            tinyobjloader 
                                            stb_image)

add_library(bssrdf bssrdf.cpp)
target_include_directories(bssrdf 
PRIVATE ${CMAKE_SOURCE_DIR}/lib/DataStracture
PRIVATE ${CMAKE_SOURCE_DIR}/lib/Math
)                                            

add_executable(test_bssrdf test_bssrdf.cpp)
target_include_directories(test_bssrdf
    PRIVATE ${CMAKE_SOURCE_DIR}/lib/Camera
    PRIVATE ${CMAKE_SOURCE_DIR}/lib/DataStracture
    PRIVATE ${CMAKE_SOURCE_DIR}/lib/RNG
    PRIVATE ${CMAKE_SOURCE_DIR}/src/Parallel
    PRIVATE ${CMAKE_SOURCE_DIR}/src/Light
    PRIVATE ${CMAKE_SOURCE_DIR}/ext/stb
    PRIVATE ${CMAKE_SOURCE_DIR}/ext/tinyobjloader
    PRIVATE ${CMAKE_SOURCE_DIR}/lib/Math)
    target_link_libraries(test_bssrdf       bssrdf
                                            Camera 
                                            random
                                            IBL 
                                            parallel 
                                            Accelerator 
                                            Vec 
                                            stb_image_write 
                                            tinyobjloader 
                                            stb_image)