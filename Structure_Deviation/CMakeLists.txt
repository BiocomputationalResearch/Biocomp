cmake_minimum_required(VERSION 2.6)
set(CMAKE_CXX_STANDARD 11)
cmake_policy(SET CMP0037 NEW)
project("Structure_Deviation")
include_directories("add")
file(GLOB SOURCES "src/*.cpp")
add_executable("Structure_Deviation" ${SOURCES})
target_link_libraries("Structure_Deviation" "openbabel")
