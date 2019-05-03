#include <stdio.h>

#ifndef LDPC_H
#define LDPC_H

h_matrix * create_ldpc(size_t n, size_t j, size_t k);
h_matrix * create_generator_matrix(h_matrix *mat);
h_list * encode_ldpc(h_matrix *gen, h_list *mes);

#define cldpc create_ldpc
#define cgm create_generator_matrix

#endif
