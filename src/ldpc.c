#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "ldpc.h"
#include "list.h"
#include "random.h"
#include "basic.h"

// Create a ldpc matrix as described in Gallager's paper from 1963
h_matrix * create_base(size_t n, size_t j, size_t k)
{
    if (n % k)
    {
        return NULL;
    }

    size_t m = (n * j) / k;
    h_matrix *res = chm(m, n);
    i_list *perm = cil(n, n);
    // Fill the first horizontal part of the matrix
    for (size_t i = 0; i < n; i++)
    {
        shm(res, 1, (i / k), i);
    }

    // Fill the j-1 other bands
    for (size_t i = 1; i < j; i++)
    {
        permutation(perm);
        for (size_t x = 0; x < n; x++)
        {
            shm(res, 1, ((perm->list[x] + i * n) / k), x);
        }
    }
    return res;
}


// Create the encoder matrix with the base matrix
h_matrix * create_generator_matrix_h(h_matrix *mat)
{
    return juxtapose_h(mat, 1);
}


// Create the encoder matrix with the base matrix
a_matrix * create_generator_matrix_a(h_matrix *mat)
{
    return juxtapose_a(mat, 1);
}


// Create the decoder matrix with the base matrix
h_matrix * create_decoder_matrix_h(h_matrix *mat)
{
    return juxtapose_h(mat, 0);
}


// Create the decoder matrix with the base matrix
a_matrix * create_decoder_matrix_a(h_matrix *mat)
{
    return juxtapose_a(mat, 0);
}


// Encode a message mes with the generator matrix gen
h_list * encode_ldpc_h(h_matrix *gen, h_list *mes)
{
    return product_h(gen, mes);
}


// Encode a message mes with the generator matrix gen
h_list * encode_ldpc_a(a_matrix *gen, h_list *mes)
{
    return product_a(gen, mes);
}


// Decode the received message res with the decoding matrix mat
// Returns 0 if the decoding is succesful or the number of iterations if not
int decode_ldpc_a_basic(a_matrix *mat, h_list *mes, size_t nb_max)
{
    h_list *verif = product_a(mat, mes);
    i_list *count = cil(mes->n, mes->n);
    i_list *max_errors = cil(mes->n, mes->n);
    size_t iter = 0;

    char correct = is_all_nil(verif);


    while (!correct && (iter < nb_max))
    {
        set_all_i_list(count, 0);

        // Count the number of errors for each bit
        for (size_t i = 0; i < verif->n; i++)
        {
            if (verif->list[i])
            {
                for (size_t j = 0; j < mat->list_n[i]->n; j++)
                {
                    count->list[mat->list_n[i]->list[j]] ++;
                }
            }
        }

        // Flip the bits with the most errors
        max_i_list(count, max_errors);
        for (size_t i = 0; i < max_errors->n; i++)
        {
            mes->list[max_errors->list[i]] ^= 1;
        }

        iter ++;

        product_a_in_place(mat, mes, verif);
        correct = is_all_nil(verif);
    }

    // Free used lists
    free_h_list(verif);
    free_i_list(count);
    free_i_list(max_errors);

    return iter;
}


int demo_ldpc_basic(char *file_name, double s, size_t i_max)
{
    char debug = 1;
    int err;


    char path[1024];        // The source file
    char path_n[1024];      // The noisy file
    char path_d[1024];      // The decoded file

    getcwd(path, 1024);
    strcat(path, "/../files/bytes/");
    strcat(path, file_name);

    strcpy(path_n, path);
    strcpy(path_d, path);

    strcat(path, ".bt");
    strcat(path_n, "_ldpc_n.bt");
    strcat(path_d, "_ldpc_d.bt");

    FILE *fp = fopen(path, "r");
    FILE *fn = fopen(path_n, "w");
    FILE *fd = fopen(path_d, "w");

    if (debug)
    {
        printf("Opened files\n");
    }

    hh_list *mes_d = read_bit_file_h(fp, 8920);
    hh_list *mes_n = deep_copy(mes_d);

    if (debug)
    {
        printf("Created hh_lists\n");
    }

    h_matrix *base = create_base(8920, 40, 20);
    a_matrix *gen = cgm_a(base);
    a_matrix *dec = cdm_a(base);

    if (debug)
    {
        printf("Created LDPC matrices\n");
    }

    printf("Number of slices : %d\n", mes_d->n);

    for (size_t i = 0; i < mes_d->n; i++)
    {
        printf("\nSlice number %d :\n", i);

        h_list *res = encode_ldpc_a(gen, mes_d->list[i]);
        printf("\t- Encoded\n");

        s_list *noisy = add_noise(res, s);
        printf("\t- Added gaussian noise\n");

        h_list *noisy_h = decode_h_basic(noisy);
        if (debug)
        {
            printf("\t- Performed basic decoding\n");
        }

        copy_h(noisy_h, mes_n->list[i], 8920, 1);
        if (debug)
        {
            printf("\t- Copied basic decoding\n");
        }

        err = decode_ldpc_a_basic(dec, noisy_h, i_max);
        printf("\t- Decoded\n");
        if (debug)
        {
            printf("\t\t- Iterations : %d\n", err);
        }


        copy_h(noisy_h, mes_d->list[i], 8920, 1);
        if (debug)
        {
            printf("\t- Copied decoding\n");
        }

        free_h_list(noisy_h);
        free_s_list(noisy);
        free_h_list(res);
    }

    write_bit_file_h(fn, mes_n);
    write_bit_file_h(fd, mes_d);

    fclose(fp);
    fclose(fn);
    fclose(fd);

    free_hh_list(mes_d);
    free_hh_list(mes_n);

    free_h_matrix(base);
    free_a_matrix(gen);
    free_a_matrix(dec);

    return 0;
}
