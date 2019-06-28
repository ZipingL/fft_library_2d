cmake_minimum_required(VERSION 3.14)
project(fft2d_natural_c C)

set(CMAKE_C_STANDARD 99)

FILE(GLOB MyCSources "${CMAKE_CURRENT_SOURCE_DIR}/libfixmath/*.c")

#add paths for compiler to look for header files in
include_directories(./)

#create execu main.c with required .c and .h files
add_executable(fft2d_natural_c ${MyCSources}
        main.c

        data_in.h
        FFT2d_cn.c
        FFT2d_cn.h)
set_target_properties(fft2d_natural_c PROPERTIES LINK_FLAGS -Wl,--stack,100000000)
