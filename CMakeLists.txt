cmake_minimum_required(VERSION 3.14)
project(fft1d C)

set(CMAKE_C_STANDARD 99)
FILE(GLOB MyCSources "${CMAKE_CURRENT_SOURCE_DIR}/libfixmath/*.c")

add_executable(fft1d ${MyCSources} main.c data_in.h FFTLib.c FFTLib.h)
