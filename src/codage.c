#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "../lib/BitArray/bit_array.h"
#include "../lib/BitArray/bar.h"

#include "turbocode.h"
#include "ldpc.h"
#include "list.h"
#include "random.h"
#include "basic.h"
#include "demo.h"


int main(void)
{
    // Initialise the random generator
    srand((unsigned int)time(NULL));
    size_t x;
    double y;
    double z;

    //demo_ldpc_basic("junob", 0.3, 200);
    demo_turbo_iter("junob", 0.5, 10);
//    printf("Noise sigma : %f\n", s);
/*
    h_matrix *test = create_base(2000, 35, 20);
    a_matrix *gen_a = cgm_a(test);
    a_matrix *dec_a = cdm_a(test);
    h_list *mes = chl(2000, 2000);

    set_all_h_list(mes, 1);
    //print_h_list(mes);
    //printf("\n");


    h_list *res_a = encode_ldpc_a(gen_a, mes);
    s_list *noisy = add_noise(res_a, 0.6);
    h_list *noisy_h = decode_h_basic(noisy);

    int r = nb_errors(res_a, noisy_h);

    printf("%d -- %f \n", r, r / (float) 2000);

    printf("%d\n", decode_ldpc_a_basic(dec_a, noisy_h, 2000));

    //print_h_list(noisy_h);

    r = nb_errors(res_a, noisy_h);
    printf("%d -- %f \n", r, r / (float) 2000);

    free_s_list(noisy);
    free_h_list(noisy_h);
    free_h_matrix(test);
    free_a_matrix(gen_a);
    free_a_matrix(dec_a);
    free_h_list(mes);
    free_h_list(res_a);
*/

    return 0;
}
