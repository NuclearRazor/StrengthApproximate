# Author: Ivan Blagopoluchnyy 2018 (c)

cmake_minimum_required(VERSION 3.2.1)

#Use C++11 standart methods in code
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# set config for project
set (PROJECT_NAME "StrengthApproximate" CXX)
project (${PROJECT_NAME})
set (PROJECT_SOURCE_DECLARATION_DIRECTORY ${PROJECT_SOURCE_DIR}/include)
set (PROJECT_SOURCE_DEFINITION_DIRECTORY ${PROJECT_SOURCE_DIR}/src)
set (MAIN_FILE strength_approximate.cpp)

# set options used by "CMake"
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

add_executable (${PROJECT_NAME} ${MAIN_FILE})

set_target_properties(${TARGET} PROPERTIES LINKER_LANGUAGE Cxx)


