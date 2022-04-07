cmake_minimum_required(VERSION 3.21)
project(ConeCaster CXX  CUDA)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CUDA_STANDARD 20)

add_compile_options("/Oy-")

add_subdirectory("./source")
add_subdirectory("./include")
target_link_libraries(App PUBLIC Lib)
