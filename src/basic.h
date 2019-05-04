#include <stdio.h>

#include "list.h"

#ifndef BASIC_H
#define BASIC_H

h_list * decode_h_basic(s_list *mes);
s_list * add_noise(h_list *mes, double s);
int nb_errors(h_list *og_mes, h_list *mes);

#endif
