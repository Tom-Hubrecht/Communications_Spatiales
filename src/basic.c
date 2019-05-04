#include <stdio.h>

#include "basic.h"
#include "random.h"
//#include "list.h"


// Decode the s_list by comparing its values to 0
h_list * decode_h_basic(s_list *mes)
{
    h_list *res = chl(mes->n);
    for (size_t i = 0; i < mes->n; i++)
    {
        res->list[i] = (mes->list[i] > 0);
    }
    return res;
}


// Add a white noise with variance s^2 to the map of mes where
// 0 -> -1
// 1 -> 1
s_list * add_noise(h_list *mes, double s)
{
    s_list *res = csl(mes->n);
    for (size_t i = 0; i < res->n; i++)
    {
        res->list[i] = (2.0 * mes->list[i]) - 1.0 + box_muller(0.0, s);
    }
    return res;
}


// Computes the number of differences between og_mes and mes
int nb_errors(h_list *og_mes, h_list *mes)
{
    if (og_mes->n != mes->n)
    {
        return -1;
    }
    int d;
    for (size_t i = 0; i < mes->n; i++)
    {
        d += (og_mes->list[i] != mes->list[i]);
    }
    return d;
}

