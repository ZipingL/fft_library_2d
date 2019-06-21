//
// Created by Ziping Liu on 6/20/2019.
//

#include "FFTLib.h"
/* declared in header file
struct BinaryTree
{
    cmplx* output_buf;
    struct BinaryTree *parent;
    struct BinaryTree *left;
    struct BinaryTree *right;
};

//TODO: free tree
typedef struct BinaryTree BinaryTree; */

void free_node(BinaryTree* tree)
{
    if(tree->output_buf)
        free(tree->output_buf);
    if(tree->samples_even)
        free(tree->samples_even);
    if(tree->samples_odd)
        free(tree->samples_odd);
}


BinaryTree* createNode(size_t size)
{
    BinaryTree* tree = malloc(sizeof(BinaryTree));

    tree->output_buf = malloc(sizeof(cmplx)*size);

    if(size >= 2)
    tree->samples_even = malloc(sizeof(cmplx)*(size/2));
    else
        tree->samples_even = NULL;

    if(size >= 2)
    tree->samples_odd = malloc(sizeof(cmplx)*(size/2));
        tree->samples_odd = NULL;


    tree->left = NULL;
    tree->right = NULL;
}

void show(const char * s, cmplx buf[], size_t N) {
    printf("%s", s);
    for (int i = 0; i < N; i++)
        printf("(%g, %gj) %s ", creal(buf[i]), cimag(buf[i]), i == (N-1) ? "":"|");

    printf("\n");
}


void show_fixed(const char * s, fixed_complex buf[], size_t N)
{
    printf("%s", s);
    for (int i = 0; i < N; i++)
        printf("(%f, %fj) %s ", fix16_to_dbl(buf[i].real),fix16_to_dbl(buf[i].imaginary), i == (N-1) ? "":"|");

    printf("\n");
}




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
    if(f_even != NULL && size_halved != 1)
        free(f_even);
    if(f_even != NULL && size_halved != 1)
        free(f_odd);
    return output;
}




cmplx** FFT2(cmplx** samples, size_t num_rows, size_t num_cols)
{


    /* Do a set of FFT's on the columns of the input image "samples" */
    // Time complexity = num_col*num_rows*log(num_rows)
    cmplx* column = malloc(sizeof(cmplx)*num_rows);

    cmplx ** output_pre_transpose = malloc(sizeof(cmplx*)*num_cols);
    for(int i = 0; i < num_rows; i ++)
    {
        output_pre_transpose[i] = malloc(sizeof(cmplx)*num_rows);

        for(int j = 0; j < num_rows; j++)
            column[j] = samples[j][i];

        output_pre_transpose[i] = FFT(column, num_rows);

        if(i == 0)
            show("first FFT on first col: ", output_pre_transpose[i], num_rows);
    }

    // Transpose the output //
    // Use the samples parameter to store the transposed matrix
    for(int i = 0; i < num_rows; i++)
        for(int j = 0; j < num_cols; j++)
            samples[j][i] = output_pre_transpose[i][j];

    /* Do a set of FFT's on the rows of the input image "samples" */
    // Time complexity = row_length*col_length*log(col_length)

    for(int i = 0; i < num_cols; i++)
    {
        output_pre_transpose[i] = FFT(samples[i], num_rows);

        if(i == 0)
            show("first FFT on first row: ", output_pre_transpose[i], num_rows);
    }

    /* Total Time complexity == row_length*col_length*log(row_len + col_len) */





    return output_pre_transpose;
}

BinaryTree* FFT_preallocate_memory(size_t size)
{
    if(size == 1)
    {
        return createNode(size);
    }


    size_t size_halved = size / 2;

    BinaryTree* right = FFT_preallocate_memory(size_halved);
    BinaryTree* left = FFT_preallocate_memory(size_halved);

    BinaryTree* root = createNode(size);
    root->right = right;
    root->left = left;

    return root;
}


void FFT_free_tree(BinaryTree* tree)
{
    if(tree == NULL)
        return;

    FFT_free_tree(tree->left);
    FFT_free_tree(tree->right);

    free_node(tree);
}

BinaryTree* FFT_preallocation_expected(cmplx* samples, size_t size, BinaryTree *output)
{

    // base case, when the signal is split and transformed on to a point that it is length of 1
    // the confusing this to remember is when the signal is split to length of one,
    // the FFT returns the samples, not the transformed samples.
    // so in the case that size is say 2, this if statement is not hit.
    // that is when the samples[] actually get transformed to frequency space
    if(size == 1)
    {
        output->output_buf[0] = samples[0];
        return output;
    }

    // Split signal into even and odd indices
    size_t size_halved = size / 2;

    cmplx x_even[size_halved];
    cmplx x_odd[size_halved];

    for(int i = 0; i < size_halved; i++)
    {
        x_even[i] = samples[2*i];
        x_odd[i] = samples[2*i+1];
    }

    BinaryTree* f_even = FFT_preallocation_expected(x_even, size_halved, output->right);
    BinaryTree* f_odd = FFT_preallocation_expected(x_odd, size_halved, output->left);

    for(int i = 0; i < size_halved; i++)
    {
        // for k values 0 to k/2
        // 0 to pi
        cmplx exponential = cexp(-2*I*M_PI*i/size) * f_odd->output_buf[i];
        output->output_buf[i] = f_even->output_buf[i] + exponential;

        // for k values k/2 to size (N)
        // exploiting symmetry of sinusoidal
        // exponential is negative due to the pi to 2pi,
        // thus the addition of a constant pi to the exponential,
        // but instead of recalculating the exponential,
        // adding a negative in front does the same job
        output->output_buf[i + size_halved] = f_even->output_buf[i] - exponential;
    }

    return output;

}


void FFT2_preallocation_expected(cmplx** samples, size_t num_rows, size_t num_cols, cmplx** output_pre_transpose,
                                   BinaryTree* output[])
{


    /* Do a set of FFT's on the columns of the input image "samples" */
    // Time complexity = num_col*num_rows*log(num_rows)
    cmplx* column = malloc(sizeof(cmplx)*num_rows);

    //cmplx ** output_pre_transpose = malloc(sizeof(cmplx*)*num_rows);
    for(int i = 0; i < num_cols; i ++)
    {
        output_pre_transpose[i] = malloc(sizeof(cmplx)*num_rows);

        for(int j = 0; j < num_rows; j++)
            column[j] = samples[j][i];

        FFT_preallocation_expected(column, num_rows, output[i]);
        output_pre_transpose[i] = output[i]->output_buf;

        if(i == 0)
            show("first FFT on first col: ", output_pre_transpose[i], num_rows);
    }

    // Transpose the output //
    // Use the samples parameter to store the transposed matrix
    for(int i = 0; i < num_rows; i++)
        for(int j = 0; j < num_cols; j++)
            samples[j][i] = output_pre_transpose[i][j];

    /* Do a set of FFT's on the rows of the input image "samples" */
    // Time complexity = row_length*col_length*log(col_length)

    for(int i = 0; i < num_cols; i++)
    {
        FFT_preallocation_expected(samples[i], num_cols,output[i]);
        output_pre_transpose[i] = output[i]->output_buf;

        if(i == 0)
            show("first FFT on first row: ", output_pre_transpose[i], num_rows);
    }

    /* Total Time complexity == row_length*col_length*log(row_len + col_len) */


}

fixed_complex* FFT_fixed(fixed_complex* samples, size_t size)
{
    // base case, when the signal is split and transformed on to a point that it is length of 1
    // the confusing this to remember is when the signal is split to length of one,
    // the FFT returns the samples, not the transformed samples.
    // so in the case that size is say 2, this if statement is not hit.
    // that is when the samples[] actually get transformed to frequency space
    if(size == 1)
    {
        fixed_complex* output = (fixed_complex*) malloc(sizeof(fixed_complex)*size);
        output[0] = samples[0];
        return output;
    }

    // Split signal into even and odd indices
    int size_halved = size / 2;

    fixed_complex x_even[size_halved];
    fixed_complex x_odd[size_halved];

    for(int i = 0; i < size_halved; i++)
    {
        x_even[i] = samples[2*i];
        x_odd[i] = samples[2*i+1];
    }

    fixed_complex* f_even = NULL;
    f_even = FFT_fixed(x_even, size_halved);
    fixed_complex* f_odd = NULL;
    f_odd =  FFT_fixed(x_odd, size_halved) ;

    fixed_complex* output = (fixed_complex*) malloc(sizeof(fixed_complex)*size);
    for(int i = 0; i < size_halved; i++)
    {
        // for k values 0 to k/2
        // 0 to pi
        /*         cmplx exponential = cexp(-2*I*M_PI*i/size) * f_odd[i]; */

        fix16_t a = fix16_smul(fix16_from_int(-2), fix16_pi);
        fix16_t b = fix16_sdiv(fix16_from_dbl(i),fix16_from_int((int)size));
        fix16_t cosine_part = fix16_cos(fix16_smul(a,b));

        a = fix16_smul(fix16_from_int(-2), fix16_pi);
        b = fix16_sdiv(fix16_from_dbl(i),fix16_from_int((int)size));
        fix16_t sine_part = fix16_sin(fix16_smul(a,b));


        fix16_t exp_real  = fix16_smul(  cosine_part, f_odd[i].real);
        fix16_t exp_imag  = fix16_smul( sine_part , f_odd[i].imaginary);


        /*         output[i] = f_even[i] + exponential; */
        output[i].real = fix16_sadd(f_even[i].real , exp_real);
        output[i].imaginary = fix16_sadd(f_even[i].imaginary , exp_imag);

        // for k values k/2 to size (N)
        // exploiting symmetry of sinusoidal
        // exponential is negative due to the pi to 2pi,
        // thus the addition of a constant pi to the exponential,
        // but instead of recalculating the exponential,
        // adding a negative in front does the same job
        // output[i + size_halved] = f_even[i] - exponential;
        fix16_t negate = fix16_from_int(-1);
        output[i + size_halved].real = fix16_sadd(f_even[i].real , fix16_smul(negate, exp_real));
        output[i + size_halved].imaginary = fix16_sadd(f_even[i].imaginary , fix16_smul(negate, exp_imag));
    }

    if(f_even != NULL  )
        free(f_even);
    if(f_even != NULL  )
        free(f_odd);
    return output;
}

cmplx** FFT2_fixed(cmplx** samples, size_t num_rows, size_t num_cols)
{


    /* Do a set of FFT's on the columns of the input image "samples" */
    // Time complexity = num_col*num_rows*log(num_rows)
    cmplx* column = malloc(sizeof(cmplx)*num_rows);

    cmplx ** output_pre_transpose = malloc(sizeof(cmplx*)*num_cols);
    for(int i = 0; i < num_rows; i ++)
    {
        output_pre_transpose[i] = malloc(sizeof(cmplx)*num_rows);

        for(int j = 0; j < num_rows; j++)
            column[j] = samples[j][i];

        output_pre_transpose[i] = FFT(column, num_rows);

        if(i == 0)
            show("first FFT on first col: ", output_pre_transpose[i], num_rows);
    }

    // Transpose the output //
    // Use the samples parameter to store the transposed matrix
    for(int i = 0; i < num_rows; i++)
        for(int j = 0; j < num_cols; j++)
            samples[j][i] = output_pre_transpose[i][j];

    /* Do a set of FFT's on the rows of the input image "samples" */
    // Time complexity = row_length*col_length*log(col_length)

    for(int i = 0; i < num_cols; i++)
    {
        output_pre_transpose[i] = FFT(samples[i], num_rows);

        if(i == 0)
            show("first FFT on first row: ", output_pre_transpose[i], num_rows);
    }

    /* Total Time complexity == row_length*col_length*log(row_len + col_len) */





    return output_pre_transpose;
}