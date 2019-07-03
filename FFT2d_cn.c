//
// Created by Ziping Liu on 6/20/2019.
//

#include "FFT2d_cn.h"

void generateTwiddles(size_t size, complex_number_t ** fixed_out)
{
    int32_t size_temp = size;
    int32_t twiddle_size = 0;

    while(size_temp > 1)
    {
        size_temp = floor(size_temp  /2);
        twiddle_size += size_temp;
    }


    *fixed_out = (complex_number_t*)malloc(sizeof(complex_number_t)*twiddle_size);

    size_temp = size;
    int j = 0;

    int stage = 1;
    while(size_temp > 1)
    {


        for(int i = 0; i < floor(size_temp / 2 ); i ++,j++)
        {

            (*fixed_out)[j].real =cos(-2*M_PI*i / floor(size_temp)) * pow(2, 15);
            (*fixed_out)[j].imag =sin(-2*M_PI*i / floor(size_temp)) * pow(2,15) ;
        }

        stage++;
        size_temp = floor(size_temp / 2);
    }
}

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

void show(const char * s, cmplx *buf, size_t N) {
    printf("%s", s);
    for (int i = 0; i < N; i++)
        printf("(%g, %gj) %s ", creal(buf[i]), cimag(buf[i]), i == (N-1) ? "":"|");

    printf("\n");
}
void show_fixed(const char * s, complex_number_t buf[], size_t N) {
    printf("%s", s);
    for (int i = 0; i < N; i++) {
        printf("(%.3f, %.3fj) %s ", ((float) buf[i].real) / pow(2, 15), ((float) buf[i].imag) / pow(2, 15),
             "" );

        //printf("(%.3f, %.3fj) %s ", ((float) buf[i].real), ((float) buf[i].imag) ,
           //    i == (N - 1) ? "" : "|");
    }
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

        printf("DEBUG FFT: (%g, %gj) \n"
               "EXP (%g, %gj) \n"
               "A: (%g, %gj) \n",
               creal(output[i]) , cimag(output[i]),
               creal(exponential), cimag(exponential),
               creal(cexp(-2*I*M_PI*i/size)), cimag( cexp(-2*I*M_PI*i/size)) );

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

complex_number_t * FFT_fixed(complex_number_t* samples, size_t size, complex_number_t* twiddles)
{

    // base case, when the signal is split and transformed on to a point that it is length of 1
    // the confusing this to remember is when the signal is split to length of one,
    // the FFT returns the samples, not the transformed samples.
    // so in the case that size is say 2, this if statement is not hit.
    // that is when the samples[] actually get transformed to frequency space
    if(size == 1)
    {
        complex_number_t* out = malloc(sizeof(complex_number_t));
        out[0].real = samples[0].real;
        out[0].imag = samples[0].imag;
        return out;
    }

    // Split signal into even and odd indices
    int size_halved = size / 2;

   // complex_number_t x_even[size_halved];
   // complex_number_t x_odd[size_halved];

    complex_number_t* x_even = malloc(sizeof(complex_number_t)*size_halved);
    complex_number_t* x_odd = malloc(sizeof(complex_number_t)*size_halved);

    for(int i = 0; i < size_halved; i++)
    {
        x_even[i] .real = samples[2*i].real;
        x_even[i] .imag = samples[2*i].imag;
        x_odd[i].real = samples[2*i+1].real;
        x_odd[i].imag = samples[2*i+1].imag;

    }

    complex_number_t* f_even = NULL;
    f_even = FFT_fixed(x_even, size_halved, size_halved == 1 ? NULL : twiddles + size_halved);
    complex_number_t* f_odd = NULL;
    f_odd =  FFT_fixed(x_odd, size_halved, size_halved == 1 ? NULL : twiddles + size_halved) ;

    complex_number_t* output =  malloc(sizeof(complex_number_t)*size);
    for(int i = 0; i < size_halved; i++)
    {
        // for k values 0 to k/2
        // 0 to pi
       // cmplx exponential = cexp(-2*I*M_PI*i/size) * f_odd[i];
       // output[i] = f_even[i] + exponential;

       // Multiply complex exp, e.g.
       // (3 + 2i)(1 + 4i) = 3 + 12i + 2i + 8i^2.
        complex_number_t exponential;
       int64_t exp_temp_r = (f_odd[i].real )*(twiddles[i].real ) ;
       int64_t exp_temp_i = (f_odd[i].imag*twiddles[i].imag);

       int64_t complex_product_real = exp_temp_r - exp_temp_i;



       exp_temp_r = f_odd[i].real*twiddles[i].imag;
       exp_temp_i = f_odd[i].imag*twiddles[i].real;

       int64_t complex_product_imag = exp_temp_r + exp_temp_i;

        exponential.real  =(complex_product_real >> 15);
        exponential.imag = (complex_product_imag >> 15);

        output[i].real = f_even[i].real +  exponential.real;
        output[i].imag = f_even[i].imag  + exponential.imag;


        // for k values k/2 to size (N)
        // exploiting symmetry of sinusoidal
        // exponential is negative due to the pi to 2pi,
        // thus the addition of a constant pi to the exponential,
        // but instead of recalculating the exponential,
        // adding a negative in front does the same job


        //output[i + size_halved] = f_even[i] - exponential;

        output[i+ size_halved].real = f_even[i].real -  exponential.real;
        output[i+ size_halved].imag = f_even[i].imag  -exponential.imag;


       /* printf("DEBUG FFT %d %d: (%g, %gj) \n"
               "EXP (%g, %gj) \n"
               "Twiddle: (%g, %gj) \n"
              // "Even FFT: (%g, %gj)\n"
               //"Odd FFT: (%g, %gj)\n"
               ,

            i, size,
            output[i].real / pow(2, 15)*1000, output[i].imag / pow(2,15)*1000,
               exponential.real / pow(2, 15)*1000, exponential .imag/ pow(2,15)*1000,
               twiddles[i].real / pow(2,15), twiddles[i].imag / pow(2,15)
               //,
               //f_even[i].real / pow(2, 15)*1000, f_even[i].imag/ pow(2,15)*1000,
             // f_odd[i].real / pow(2, 15)*1000, f_odd[i].imag/ pow(2,15)*1000
             );*/

     //   printf( "(%g, %gj)\n",
             //   twiddles[i].real / pow(2,15), twiddles[i].imag / pow(2,15));

    }

    free(x_even);
    free(x_odd);
    if(f_even != NULL && size_halved != 1)
        free(f_even);
    if(f_even != NULL && size_halved != 1)
        free(f_odd);
    return output;
}
void scale_buffer(complex_number_t * buf, size_t size, int scale_by)
{
    for (int i = 0; i < size; i ++)
    {
        buf[i].real = buf[i].real / (int) (pow(2, 15)*scale_by);
        buf[i].imag = buf[i].imag / (int) (pow(2, 15)*scale_by);

    }
}


complex_number_t** FFT2_fixed(complex_number_t** samples, size_t num_rows, size_t num_cols)
{
    complex_number_t * twiddles = malloc(sizeof(complex_number_t)*num_rows);

    generateTwiddles(num_cols, &twiddles);

    /* Do a set of FFT's on the columns of the input image "samples" */
    // Time complexity = num_col*num_rows*log(num_rows)
    complex_number_t* column = malloc(sizeof(complex_number_t)*num_rows);

    complex_number_t ** output_pre_transpose = malloc(sizeof(complex_number_t*)*num_cols);
    for(int i = 0; i < num_rows; i ++)
    {
       // output_pre_transpose[i] = malloc(sizeof(complex_number_t)*num_rows);

        for(int j = 0; j < num_rows; j++) {
            column[j].real = samples[j][i].real;
            column[j].imag = samples[j][i].imag;

        }
        output_pre_transpose[i] = FFT_fixed(column, num_rows, twiddles);

        if(i == 31)
            show_fixed("first FFT on first col: ", output_pre_transpose[i], num_rows);
    }


    // Scale and Transpose the output //
    // Use the samples parameter to store the transposed matrix
    for(int i = 0; i < num_rows; i++)
        for(int j = 0; j < num_cols; j++) {

            int64_t real = output_pre_transpose[i][j].real * pow(2,15);
            real = real / (10*pow(2,15));
            int64_t imag = output_pre_transpose[i][j].imag * pow(2,15);
            imag = imag / (10*pow(2,15));

            samples[j][i].real = real;
            samples[j][i].imag =imag;

        }
    /* Do a set of FFT's on the rows of the input image "samples" */
    // Time complexity = row_length*col_length*log(col_length)

    for(int i = 0; i < num_cols; i++)
    {
        output_pre_transpose[i] = FFT_fixed(samples[i], num_rows, twiddles);

        if(i == 0)
            show_fixed("first FFT on first row: ", output_pre_transpose[i], num_rows);
    }

    /* Total Time complexity == row_length*col_length*log(row_len + col_len) */

    return output_pre_transpose;
}
