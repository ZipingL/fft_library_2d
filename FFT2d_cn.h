//
// Created by Ziping Liu on 6/20/2019.
//

#ifndef FFT1D_FFT2D_CN_H
#define FFT1D_FFT2D_CN_H

#include "complex.h"
#include "stdio.h"
#include "math.h"
#include <stdlib.h>
#include "stdbool.h"
#include "array.h"

struct complex_number
{
    int32_t real;
    int32_t imag;
}; typedef struct complex_number complex_number_t;

void generateTwiddles(size_t size, complex_number_t ** fixed_out);
complex_number_t * FFT_fixed(complex_number_t* samples, size_t size, complex_number_t* twiddles);
complex_number_t** FFT2_fixed(complex_number_t** samples, size_t num_rows, size_t num_cols);



typedef double complex cmplx;

struct BinaryTree
{
    cmplx* output_buf;
    struct BinaryTree *parent;
    struct BinaryTree *left;
    struct BinaryTree *right;
    cmplx* samples_even;
    cmplx* samples_odd;
};

//TODO: free tree
typedef struct BinaryTree BinaryTree;

void show(const char * s, cmplx buf[], size_t N);
void show_fixed(const char * s, complex_number_t *buf, size_t N);

/*
 * Returns an array of complex numbers that is the FFT of the given array of samples with specified size.
 * Creates an output buffer using malloc that user needs to free
 * Params: samples, size
 * Returns: cmplx array
 */
cmplx* FFT(cmplx* samples, size_t size);

/*
 * Takes in a 2D array of complex signals and returns the 2-D FFT
 * Creates a matrix output buffer using malloc that user needs to free
 */
cmplx** FFT2(cmplx** samples, size_t num_rows, size_t num_cols);

/*
 * Works exactly the same as FFT function but requires all buffers to be allocated by user before calling the function
 * Please use FFT_preallocate_memory passing into output
 */
BinaryTree* FFT_preallocation_expected(cmplx* samples, size_t size, BinaryTree *output);


/*
 * Works exactly the same as FFT function but requires all buffers to be allocated by user before calling the function
 * Please create a vector of trees generated by FFT_preallocate memory, passed in through output[]
 */
void FFT2_preallocation_expected(cmplx** samples, size_t num_rows, size_t num_cols, cmplx** output_pre_transpose,
                               BinaryTree* output[]);

/*
 * Creates a binary tree. Each nodes stores the need buffers required by each recursive FFT call.
 * For example, a tree of 7 nodes are created for a signal size of 4
 * Since for a FFT on a signal of length 4, the FFT is recursively called 6 times (not including itself)
 *            0
 *          /   \
 *         0     0
 *        / \   / \
 *       0   0 0   0
 * Each node containing all the memory need by each recursive function call
 */
BinaryTree* FFT_preallocate_memory(size_t size);
void FFT_free_tree(BinaryTree* tree);


#endif //FFT1D_FFT2D_CN_H