cmake_minimum_required(VERSION 3.14)

project(testSNARK LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_SOURCES
    snark.h
    circut.h circut.cpp
    setup.h setup.cpp
    funciton.h funciton.cpp
    main.cpp
)

add_executable(testSNARK
    ${PROJECT_SOURCES})
