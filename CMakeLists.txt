cmake_minimum_required(VERSION 3.15)
project(hello_roadrunner)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Path to RoadRunner base
set(RR_ROOT "${CMAKE_SOURCE_DIR}/externals/libroadrunner-release")

# Executable
add_executable(hello_roadrunner agnostic_rr.cpp)

# Suppress warnings for MSVC
if (MSVC)
    target_link_options(hello_roadrunner PRIVATE /ignore:4217)
endif()

# Include headers
target_include_directories(hello_roadrunner PRIVATE
    "${RR_ROOT}/include"
)

# Link RoadRunner and dependencies using build configuration (Debug vs Release)
target_link_directories(hello_roadrunner PRIVATE
    "${RR_ROOT}/lib"
)

target_link_libraries(hello_roadrunner PRIVATE
    roadrunner
    libsbml-staticmd
    zlibstatic
    libexpatMD
)

# Auto-copy all .dlls from correct bin folder based on config
add_custom_command(TARGET hello_roadrunner POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${RR_ROOT}/bin"
        "$<TARGET_FILE_DIR:hello_roadrunner>"
)

# Auto-copy model files (you can add more later)
foreach(XML glycolysis.xml tca_cycle.xml oxphos.xml ppp.xml)
    add_custom_command(TARGET hello_roadrunner POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_SOURCE_DIR}/${XML}"
            "$<TARGET_FILE_DIR:hello_roadrunner>"
    )
endforeach()
