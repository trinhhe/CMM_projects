cmake_minimum_required(VERSION 3.5)

project(optLib CXX)

set (CMAKE_CXX_STANDARD 17)
FILE(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_library(${PROJECT_NAME} 
    ${SOURCES}
)
target_link_libraries(${PROJECT_NAME} 
    helpersLib
)

set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_LIST_DIR})