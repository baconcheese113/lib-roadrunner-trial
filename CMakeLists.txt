cmake_minimum_required(VERSION 3.15)
project(HelloRoadRunner)

set(CMAKE_CXX_STANDARD 17)

# Correct include paths
include_directories(${CMAKE_SOURCE_DIR}/externals/libroadrunner/include)
link_directories(${CMAKE_SOURCE_DIR}/externals/libroadrunner/lib)

add_definitions(-DLIBSBML_STATIC)

# Create the executable
add_executable(hello_roadrunner hello_roadrunner.cpp)

message(STATUS "Linking libraries from: ${CMAKE_SOURCE_DIR}/externals/libroadrunner/lib")

# Link RoadRunner + libSBML depending on build type
target_link_libraries(hello_roadrunner
    roadrunnerd
    libsbml-staticmdd
    zlibstaticd
    libexpatdMD
)


# Copy DLL after build
add_custom_command(TARGET hello_roadrunner POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_SOURCE_DIR}/externals/libroadrunner/bin/roadrunnerd.dll
        $<TARGET_FILE_DIR:hello_roadrunner>
)

add_custom_target(run ALL
    # Copy all DLLs in bin/
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/externals/libroadrunner/bin
        $<TARGET_FILE_DIR:hello_roadrunner>

    # Copy SBML model
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_SOURCE_DIR}/glycolysis.xml
        $<TARGET_FILE_DIR:hello_roadrunner>

    # Run program
    # COMMAND $<TARGET_FILE:hello_roadrunner>
    # DEPENDS hello_roadrunner
    # WORKING_DIRECTORY $<TARGET_FILE_DIR:hello_roadrunner>
)


