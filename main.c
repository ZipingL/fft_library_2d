/*Author: Ziping Liu*/

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdint.h>
#include "data_in.h"
#include <x86intrin.h>
#include "FFT2d_cn.h"
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
    } while(0)s

typedef double complex cmplx;

//#define TEST1DFFT_FLOAT
//#define TEST1DFFT_FLOAT_PREALLOCATE
//#define TEST2DFFT_FLOAT
//#define TEST2DFFT_FLOAT_preallocate
#define TESTFFT_FIXED
//#define TEST2DFFT_FIXED
int main()
{

#ifdef TESTFFT_FIXED
    printf("Starting 2dFFt cn test\n");
    int N = 32;
    complex_number_t* buf_samples = malloc(sizeof(complex_number_t)*N);
    complex_number_t* twiddles = NULL;

    for(int i = 0; i < N; i ++)
    {
        buf_samples[i] .real = creal(first_col[i]) * pow(2, 15);
        buf_samples[i] .imag = cimag(first_col[i]) * pow(2, 15);

    }

    generateTwiddles(N, &twiddles);
    complex_number_t* output = NULL;
    output = FFT_fixed(buf_samples, N, twiddles );
    show_fixed("", output, N);

    show_fixed("Twiddles\n", twiddles, N);
    free(twiddles);
#endif

    //FFT2d_with_1dBatched_cn();
/* Test the float version of FFT */
#ifdef TEST1DFFT_FLOAT
    int N = 32;
    cmplx* buf_samples = malloc(sizeof(cmplx)*N);
    for(int i = 0; i < N; i ++)
    {
       buf_samples[i] = first_col[i];

    }
    cmplx* output = NULL;
        #ifdef TEST1DFFT_FLOAT_PREALLOCATE
            BinaryTree * output_tree = FFT_preallocate_memory(N);
            FFT_preallocation_expected(buf_samples, N, output_tree);
            output = output_tree->output_buf;
            show("", output, N);
            FFT_free_tree(output_tree);

        #else
        output = FFT(buf_samples, N);
        show("", output, N);
        free(output);

        #endif

    free(buf_samples);
#endif

/* Test the float version of FFT2D */
#ifdef TEST2DFFT_FLOAT
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


#ifdef TEST2DFFT_FLOAT_preallocate
    cmplx** output_2d_prealloc = malloc(sizeof(cmplx*)*N);
    for(int i = 0; i < N; i++)
    {
        x_n_2d_malloc[i] = malloc(sizeof(cmplx)*N);
        for(int j = 0; j< N; j++)
        {
            x_n_2d_malloc[i][j] = data_in[i][j];
        }
    }

    BinaryTree** tree_vector = malloc(sizeof(BinaryTree*)*N);

    for(int i = 0; i < N; i++) {

        tree_vector[i] = FFT_preallocate_memory(N);
    }


    FFT2_preallocation_expected(x_n_2d_malloc, N, N, output_2d_prealloc, tree_vector);
    cmplx** output_2d = output_2d_prealloc;
#else
    cmplx** output_2d = FFT2(x_n_2d_malloc,N, N);
#endif

    printf("\n\nfft2 output for %dx%d input\n", N,N);

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

#ifdef TEST2DFFT_FLOAT_preallocate
        free(output_2d[i]);
        FFT_free_tree(tree_vector[i]);
#endif
    }

#ifdef TEST2DFFT_FLOAT_preallocate
    free(output_2d);
    free(tree_vector);
#endif
    free(x_n_2d_malloc);
#endif


#ifdef TEST2DFFT_FIXED
    // Test data using 32x32 input from data_in.h
    size_t N = 32;
    complex_number_t** x_n_2d_malloc = NULL;
    x_n_2d_malloc = malloc(sizeof(complex_number_t*)*N);
    for(int i = 0; i < N; i++)
    {
        x_n_2d_malloc[i] = malloc(sizeof(complex_number_t)*N);
        for(int j = 0; j< N; j++)
        {
            x_n_2d_malloc[i][j].real = creal(data_in[i][j]) * pow(2,15);
            x_n_2d_malloc[i][j].imag = cimag(data_in[i][j]) * pow(2,15);

        }
    }

    complex_number_t ** output_2d = FFT2_fixed(x_n_2d_malloc, N, N);
    for(int i = 0; i < N; i++)
    {
        printf("Row: %d: ",i);
        show_fixed("", output_2d[i], N);
    }

#endif

    return 0;
}
