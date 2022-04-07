cmake_minimum_required(VERSION 3.21)

add_library(Lib INTERFACE)
target_include_directories(Lib INTERFACE .)

target_precompile_headers(Lib INTERFACE
    Lib/core/.pch.hpp
    # Lib/file_management/.pch.hpp
    Lib/opengl/.pch.hpp
    Lib/display/.pch.hpp
)

find_package(glbinding CONFIG REQUIRED)
find_package(glfw3 CONFIG REQUIRED)
find_package(glm CONFIG REQUIRED)
find_package(imgui CONFIG REQUIRED)
find_package(RapidJSON CONFIG REQUIRED)
find_path(STB_INCLUDE_DIRS "stb.h")

target_link_libraries(Lib INTERFACE
    glbinding::glbinding glbinding::glbinding-aux
    glfw
    glm::glm
    imgui::imgui
    rapidjson)

target_include_directories(Lib INTERFACE
    ${STB_INCLUDE_DIRS})
