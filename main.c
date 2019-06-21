/*Author: Ziping Liu*/

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





typedef double complex cmplx;




#define TEST1DFFT
//#define TEST1DFFT_PREALLOCATE
 //#define TEST2DFFT
//#define TEST2DFFT_preallocate
int main()
{

#ifdef TEST1DFFT
    int N = 4;
    cmplx* buf_samples = malloc(sizeof(cmplx)*N);
    buf_samples[0] = 0;
    buf_samples[1] = 1;
    buf_samples[2] = 0;
    buf_samples[3] = -1;

    fixed_complex* buf_samples_fixed = malloc(sizeof(fixed_complex)*N);

    for(int i = 0; i < N; i ++)
    {
        buf_samples_fixed[i].real = fix16_from_dbl(buf_samples[i]);
        buf_samples_fixed[i].imaginary = fix16_to_dbl(0.0);
    }

    cmplx* output = NULL;
    #ifdef TEST1DFFT_PREALLOCATE
        BinaryTree * output_tree = FFT_preallocate_memory(N);
        FFT_preallocation_expected(buf_samples, N, output_tree);
        output = output_tree->output_buf;
        show("", output, N);
        FFT_free_tree(output_tree);

    #else
    fixed_complex* fixed_output = FFT_fixed(buf_samples_fixed, N);
    show_fixed("", fixed_output, N);
    free(output);

    #endif


    free(buf_samples);

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

    BinaryTree** tree_vector = malloc(sizeof(BinaryTree*)*N);

    for(int i = 0; i < N; i++) {

        tree_vector[i] = FFT_preallocate_memory(N);
    }


    FFT2_preallocation_expected(x_n_2d_malloc, N, N, output_2d_prealloc, tree_vector);
    cmplx** output_2d = output_2d_prealloc;
#else
    cmplx** output_2d = FFT2(x_n_2d_malloc,N, N);
#endif
    RDTSC_STOP(cycles2);

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

#ifdef TEST2DFFT_preallocate
        free(output_2d[i]);
        FFT_free_tree(tree_vector[i]);
#endif
    }

#ifdef TEST2DFFT_preallocate
    free(output_2d);
    free(tree_vector);
#endif
    free(x_n_2d_malloc);
    printf("fft2 cycle start: %10u \nfft2 cycle  stop: %10u \n", cycles1, cycles2);
#endif


    return 0;
}