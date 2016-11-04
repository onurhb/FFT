cmake_minimum_required(VERSION 3.6)
project(FFT)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")

set(SOURCE_FILES main.cpp FFTRecursive.cpp FFTRecursive.h FFT.cpp FFT.h)
add_executable(FFT ${SOURCE_FILES})