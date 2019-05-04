#include <stdio.h>

#ifndef LDPC_H
#define LDPC_H

h_matrix * create_base(size_t n, size_t j, size_t k);
h_matrix * create_generator_matrix_h(h_matrix *mat);
a_matrix * create_generator_matrix_a(h_matrix *mat);
h_list * encode_ldpc_h(h_matrix *gen, h_list *mes);
h_list * encode_ldpc_a(a_matrix *gen, h_list *mes);

#define cldpc create_ldpc
#define cgm_h create_generator_matrix_h

#endif
