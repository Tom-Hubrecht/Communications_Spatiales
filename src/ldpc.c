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
// Returns the number of iterations
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


int decode_ldpc_proba(a_matrix *mat, s_list *mes, double s, size_t nb_max)
{
    // Used to store the sign of the llr
    h_list *alpha = chl(mes->n, mes->n);

    // Used to store the absolute value of the llr
    s_list *beta = csl(mes->n, mes->n);

    // To store f(beta)
    s_list *f_beta = csl(mes->n, mes->n);

    // To store the temporary results
    s_list *tmp = csl(mes->n, mes->n);

    // To store the sum of f(beta_{i,l})
    s_list *f_sum = csl(mat->n, mat->n);

    // Initialize alpha and beta
    double l;

    for (size_t d = 0; d < mes->n; d++)
    {
        l = log(p_trans(mes->list[d], 1, s) / p_trans(mes->list[d], 0, s));
        alpha->list[d] = (l > 0 ? 1 : -1);
        beta->list[d] = fabs(l);
    }

    size_t iter = 0;
    double m = min_s(beta);

    char a;
    double b;

    while (iter < nb_max && m < 5000000)
    {
        // Fill f_beta
        for (size_t i = 0; i < beta->n; i++)
        {
            f_beta->list[i] = f(beta->list[i]);
        }

        // Compute the right sums
        product_s_ip(mat, f_beta, f_sum);

        for (size_t d = 0; d < mes->n; d++)
        {
            tmp->list[d] = alpha->list[d] * beta->list[d];

            // Compute f of the sum and the product of alpha_{i,l}
            for (size_t i = 0; i < mat->list_m[d]->n; i++)
            {
                a = alpha->list[d];
                for (size_t k = 0; k < mat->list_n[i]->n; k++)
                {
                    a *= alpha->list[mat->list_n[i]->list[k]];
                }

                tmp->list[d] += a * f(f_sum->list[0]);
            }


        }



    }



    return 0;
}
