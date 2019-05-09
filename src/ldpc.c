#include <stdio.h>
#include <stdlib.h>

#include "list.h"
#include "random.h"


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

    print_h_list(verif);

    char correct = is_all_nil(verif);


    while (!correct && (iter < nb_max))
    {
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
        print_i_list(count);
        for (size_t i = 0; i < max_errors->n; i++)
        {
            mes->list[max_errors->list[i]] ^= 1;
        }

        print_h_list(mes);

        iter ++;

        product_a_in_place(mat, mes, verif);
        print_h_list(verif);
        printf("\n");
        set_all_i_list(count, 0);
        correct = is_all_nil(verif);
    }

    if (correct)
    {
        return 0;
    }

    return iter;
}
