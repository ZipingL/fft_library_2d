//
// Created by Ziping Liu on 6/20/2019.
//

#ifndef FFT1D_FFTLIB_H
#define FFT1D_FFTLIB_H

#include "complex.h"
#include "stdio.h"
#include "math.h"
#include <stdlib.h>
#include "stdbool.h"


typedef double complex cmplx;


cmplx* FFT(cmplx* samples, size_t size);

cmplx** FFT2(cmplx** samples, size_t num_rows, size_t num_cols);


// Requires preallocation of buffers size of the signal _size
void FFT_preallocation_expected(cmplx* samples, size_t size, cmplx* buff_0, cmplx* buff_1,
                                cmplx* output_buffer_even, cmplx* output_buffer_odd, bool output_type, bool flip);

void FFT2_preallocation_expected(cmplx** samples, size_t num_rows, size_t num_cols, cmplx** output_pre_transpose,
                                 cmplx* buf_0, cmplx* buf_1, cmplx* output_buff_even, cmplx* output_buff_odd);

#endif //FFT1D_FFTLIB_H
