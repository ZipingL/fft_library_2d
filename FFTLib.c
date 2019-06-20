//
// Created by a0232832 on 6/20/2019.
//

#include "FFTLib.h"
#define swap(a, b, temp)        \
        temp = a;               \
        a = b;                  \
        b = temp;





cmplx* FFT(cmplx* samples, size_t size)
{

    // base case, when the signal is split and transformed on to a point that it is length of 1
    // the confusing this to remember is when the signal is split to length of one,
    // the FFT returns the samples, not the transformed samples.
    // so in the case that size is say 2, this if statement is not hit.
    // that is when the samples[] actually get transformed to frequency space
    if(size == 1)
    {
        return samples;
    }

    // Split signal into even and odd indices
    int size_halved = size / 2;

    cmplx x_even[size_halved];
    cmplx x_odd[size_halved];

    for(int i = 0; i < size_halved; i++)
    {
        x_even[i] = samples[2*i];
        x_odd[i] = samples[2*i+1];
    }

    cmplx* f_even = NULL;
    f_even = FFT(x_even, size_halved);
    cmplx* f_odd = NULL;
    f_odd =  FFT(x_odd, size_halved) ;

    cmplx* output = (cmplx*) malloc(sizeof(cmplx)*size);
    for(int i = 0; i < size_halved; i++)
    {
        // for k values 0 to k/2
        // 0 to pi
        cmplx exponential = cexp(-2*I*M_PI*i/size) * f_odd[i];
        output[i] = f_even[i] + exponential;

        // for k values k/2 to size (N)
        // exploiting symmetry of sinusoidal
        // exponential is negative due to the pi to 2pi,
        // thus the addition of a constant pi to the exponential,
        // but instead of recalculating the exponential,
        // adding a negative in front does the same job
        output[i + size_halved] = f_even[i] - exponential;
    }
    if(f_even != NULL)
        free(f_even);
    if(f_even != NULL)
        free(f_odd);
    return output;
}




cmplx** FFT2(cmplx** samples, size_t num_rows, size_t num_cols)
{


    /* Do a set of FFT's on the columns of the input image "samples" */
    // Time complexity = num_col*num_rows*log(num_rows)
    cmplx* column = malloc(sizeof(cmplx)*num_rows);

    cmplx ** output_pre_transpose = malloc(sizeof(cmplx*)*num_cols);
    for(int i = 0; i < num_cols; i ++)
    {
        output_pre_transpose[i] = malloc(sizeof(cmplx)*num_rows);

        for(int j = 0; j < num_rows; j++)
            column[j] = samples[j][i];

        output_pre_transpose[i] = FFT(column, num_rows);
    }

    // Transpose the output //
    // Use the samples parameter to store the transposed matrix
    for(int i = 0; i < num_cols; i++)
        for(int j = 0; j < num_rows; j++)
            samples[j][i] = output_pre_transpose[i][j];

    /* Do a set of FFT's on the rows of the input image "samples" */
    // Time complexity = row_length*col_length*log(col_length)

    for(int i = 0; i < num_rows; i++)
    {
        output_pre_transpose[i] = FFT(samples[i], num_rows);
    }

    /* Total Time complexity == row_length*col_length*log(row_len + col_len) */





    return output_pre_transpose;
}


void FFT_preallocation_expected(cmplx* samples, size_t size, cmplx* buff_0, cmplx* buff_1,
        cmplx* output_buffer_even, cmplx* output_buffer_odd, bool output_type, bool flip)
{

    // base case, when the signal is split and transformed on to a point that it is length of 1
    // the confusing this to remember is when the signal is split to length of one,
    // the FFT returns the samples, not the transformed samples.
    // so in the case that size is say 2, this if statement is not hit.
    // that is when the samples[] actually get transformed to frequency space
    if(size == 1)
    {
        if(output_type == false)
        output_buffer_even[0] = samples[0];

        else
            output_buffer_odd[0] = samples[0];
        return;
    }

    // Split signal into even and odd indices
    int size_halved = size / 2;

    cmplx x_even[size_halved];
    cmplx x_odd[size_halved];

    for(int i = 0; i < size_halved; i++)
    {
        x_even[i] = samples[2*i];
        x_odd[i] = samples[2*i+1];
    }
    complex* temp = NULL;



    cmplx* f_even = NULL;
    FFT_preallocation_expected(x_even, size_halved, output_buffer_even, output_buffer_odd, buff_0, buff_1,false ,!flip);
    cmplx* f_odd = NULL;
    FFT_preallocation_expected(x_odd, size_halved, output_buffer_even, output_buffer_odd, buff_0, buff_1,true, !flip);


    f_even = buff_0;
    f_odd = buff_1;


    cmplx* output = output_type == false ? output_buffer_even : output_buffer_odd;

    for(int i = 0; i < size_halved; i++)
    {
        // for k values 0 to k/2
        // 0 to pi
        cmplx exponential = cexp(-2*I*M_PI*i/size) * f_odd[i];
        output[i] = f_even[i] + exponential;

        // for k values k/2 to size (N)
        // exploiting symmetry of sinusoidal
        // exponential is negative due to the pi to 2pi,
        // thus the addition of a constant pi to the exponential,
        // but instead of recalculating the exponential,
        // adding a negative in front does the same job
        output[i + size_halved] = f_even[i] - exponential;
    }

}


void FFT2_preallocation_expected(cmplx** samples, size_t num_rows, size_t num_cols, cmplx** output_pre_transpose,
                                    cmplx* buf_0, cmplx* buf_1, cmplx* output_buff_even, cmplx* output_buff_odd)
{


    /* Do a set of FFT's on the columns of the input image "samples" */
    // Time complexity = num_col*num_rows*log(num_rows)
    cmplx* column = malloc(sizeof(cmplx)*num_rows);

    //cmplx ** output_pre_transpose = malloc(sizeof(cmplx*)*num_rows);
    for(int i = 0; i < num_cols; i ++)
    {
        output_pre_transpose[i] = malloc(sizeof(cmplx)*num_cols);

        for(int j = 0; j < num_rows; j++)
            column[j] = samples[j][i];

        FFT_preallocation_expected(column, num_cols, buf_0, buf_1, output_buff_even, output_buff_odd, false, false);
        output_pre_transpose[i] = output_buff_even;
    }

    // Transpose the output //
    // Use the samples parameter to store the transposed matrix
    for(int i = 0; i < num_cols; i++)
        for(int j = 0; j < num_rows; j++)
            samples[j][i] = output_pre_transpose[i][j];

    /* Do a set of FFT's on the rows of the input image "samples" */
    // Time complexity = row_length*col_length*log(col_length)

    for(int i = 0; i < num_rows; i++)
    {
        FFT_preallocation_expected(samples[i], num_rows, buf_0, buf_1, output_buff_even, output_buff_odd, false, false);
        output_pre_transpose[i] = output_buff_even;
    }

    /* Total Time complexity == row_length*col_length*log(row_len + col_len) */


}