#include <stdio.h>
#include <stdlib.h>

#include "list.h"
#include "random.h"


// Create a ldpc matrix as described in Gallager's paper from 1963
h_matrix * create_ldpc(size_t n, size_t j, size_t k)
{
    if (n % k)
    {
        return NULL;
    }

    size_t m = (n * j) / k;
    h_matrix *res = chm(m, n);
    i_list *perm = chi(n);
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


// Create the encoder matrix with the ldpc matrix
h_matrix * create_generator_matrix(h_matrix *mat)
{
    return juxtapose(mat, 1);
}


// Encode a message mes with the generator matrix gen
h_list * encode_ldpc(h_matrix *gen, h_list *mes)
{
    return product_h(gen, mes);
}


// Decode the received message res with the decoding matrix mat