cmake_minimum_required(VERSION 3.14)
project(fft1d C)

set(CMAKE_C_STANDARD 99)

add_executable(fft1d main.c data_in.h FFTLib.c FFTLib.h)