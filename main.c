#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdint.h>
#include "data_in.h"
#include <x86intrin.h>
#include "FFTLib.h"
#define __x86_64__ 1

#ifdef __i386__
#  define RDTSC_DIRTY "%eax", "%ebx", "%ecx", "%edx"
#elif __x86_64__
#  define RDTSC_DIRTY "%rax", "%rbx", "%rcx", "%rdx"
#else
# error unknown platform
#endif

#define RDTSC_START(cycles)                                \
    do {                                                   \
        register unsigned cyc_high, cyc_low;               \
        asm volatile("CPUID\n\t"                           \
                     "RDTSC\n\t"                           \
                     "mov %%edx, %0\n\t"                   \
                     "mov %%eax, %1\n\t"                   \
                     : "=r" (cyc_high), "=r" (cyc_low)     \
                     :: RDTSC_DIRTY);                      \
        (cycles) = ((uint64_t)cyc_high << 32) | cyc_low;   \
    } while (0)

#define RDTSC_STOP(cycles)                                 \
    do {                                                   \
        register unsigned cyc_high, cyc_low;               \
        asm volatile("RDTSCP\n\t"                          \
                     "mov %%edx, %0\n\t"                   \
                     "mov %%eax, %1\n\t"                   \
                     "CPUID\n\t"                           \
                     : "=r" (cyc_high), "=r" (cyc_low)     \
                     :: RDTSC_DIRTY);                      \
        (cycles) = ((uint64_t)cyc_high << 32) | cyc_low;   \
    } while(0)




/*Author: Ziping Liu*/

typedef double complex cmplx;





//#define TEST1DFFT 1
//#define TEST1DFFT_PREALLOCATE 1
#define TEST2DFFT
#define TEST2DFFT_preallocate
int main()
{

#ifdef TEST1DFFT
    int N = 4;
    cmplx* buf_0 = malloc(sizeof(cmplx)*N);
    cmplx* buf_1 = malloc(sizeof(cmplx)*N);
    cmplx* buf_samples = malloc(sizeof(cmplx)*N);
    cmplx* buf_out_even = malloc(sizeof(cmplx)*N);
    cmplx* buf_out_odd = malloc(sizeof(cmplx)*N);

    buf_samples[0] = 0;
    buf_samples[1] = 1;
    buf_samples[2] = 0;
    buf_samples[3] = -1;

    cmplx* output = NULL;
#ifdef TEST1DFFT_PREALLOCATE
    FFT_preallocation_expected(buf_samples, N, buf_0, buf_1, buf_out_even, buf_out_odd, false, false);
    output = buf_out_even;
#else
    output = FFT(buf_samples, N);
#endif

    show("", output, N);


    free(output);

#endif


#ifdef TEST2DFFT
    // Test data using 32x32 input from data_in.h
    size_t N = 32;
    cmplx** x_n_2d_malloc = NULL;
    x_n_2d_malloc = malloc(sizeof(cmplx*)*N);
    for(int i = 0; i < N; i++)
    {
        x_n_2d_malloc[i] = malloc(sizeof(cmplx)*N);
        for(int j = 0; j< N; j++)
        {
            x_n_2d_malloc[i][j] = data_in[i][j];
        }
    }
    uint64_t cycles1 = 0, cycles2 = 0;
    RDTSC_START(cycles1);


#ifdef TEST2DFFT_preallocate
    cmplx** output_2d_prealloc = malloc(sizeof(cmplx*)*N);
    for(int i = 0; i < N; i++)
    {
        x_n_2d_malloc[i] = malloc(sizeof(cmplx)*N);
        for(int j = 0; j< N; j++)
        {
            x_n_2d_malloc[i][j] = data_in[i][j];
        }
    }

    cmplx** buf_0 = malloc(sizeof(cmplx*)*N);
    cmplx** buf_1 = malloc(sizeof(cmplx*)*N);
    cmplx** buf_samples = malloc(sizeof(cmplx*)*N);
    cmplx** buf_out_even = malloc(sizeof(cmplx*)*N);
    cmplx** buf_out_odd = malloc(sizeof(cmplx*)*N);


    for(int i = 0; i < N; i++) {

        buf_0[i] = malloc(sizeof(cmplx) * N);
        buf_1[i] = malloc(sizeof(cmplx) * N);
        buf_samples[i] = malloc(sizeof(cmplx) * N);
        buf_out_even[i] = malloc(sizeof(cmplx) * N);
        buf_out_odd[i] = malloc(sizeof(cmplx) * N);

    }

    /*void FFT2_preallocation_expected(cmplx** samples, size_t num_rows, size_t num_cols,
     * cmplx** output_pre_transpose,
       cmplx* buf_0, cmplx* buf_1, cmplx* output_buff_even, cmplx* output_buff_odd); */
    FFT2_preallocation_expected(x_n_2d_malloc, N, N, output_2d_prealloc, buf_0, buf_1, buf_out_even, buf_out_odd);
    cmplx** output_2d = output_2d_prealloc;
#else
    cmplx** output_2d = FFT2(x_n_2d_malloc,N, N);
#endif
    RDTSC_STOP(cycles2);

    printf("fft2 output for %dx%d input\n", N,N);

    for(int i = 0; i < N; i++)
    {
        printf("Row: %d: ",i);
        show("", output_2d[i], N);
    }

    //Compare C code results with Octave Results
    cmplx** E = malloc(sizeof(cmplx*)*N);
    double max = 0;
    for(int i = 0; i < N; i++)
    {
        E[i] = malloc(sizeof(cmplx)*N);
        for(int j = 0; j < N; j++)
        {
            E[i][j] = ( output_2d[i][j]) - (octave_output[i][j]);

            max = max < cabs(E[i][j]) ? cabs(E[i][j]) : max;
        }
    }
    printf("Max Err %.3f\n", max);
    for(int i = 0; i < N; i++)
    {
        free(x_n_2d_malloc[i]);
        //free(output_2d[i]);
    }
    //free(output_2d);
    free(x_n_2d_malloc);
    printf("fft2 cycle start: %10u \nfft2 cycle  stop: %10u \n", cycles1, cycles2);
#endif


    return 0;
}