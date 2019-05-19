#include <stdio.h>
#include <string.h>

#include "../lib/BitArray/bit_array.h"
#include "../lib/BitArray/bar.h"
#include "basic.h"
#include "list.h"

#ifndef TURBOCODE_H
#define TURBOCODE_H

#define k1 8
#define k2 (223 * 5)

#define pTrans pTransition

extern const size_t p[];

typedef bar buffer;

size_t pi(size_t s);
char yieldEncode(char d, bar *memState);
bar * initMemState(size_t n);
h_list * encode_turbo(h_list *buf);
int decode_part(s_list *X, s_list *Y, s_list *llr, double s);

h_list * recreate(s_list *mes, char f);
int split_s(s_list *buf, s_list *X, s_list *Y1, s_list *Y2);
int interleave(s_list *X, s_list *Y);
int deinterleave(s_list *X, s_list *Y);
void min_max(double mM[2], s_list *X);
void max_min(double mM[2], s_list *X);
h_list * decode_turbo_basic(s_list *buf, double s);
h_list * decode_turbo_iter(s_list *buf, double s, size_t i_max);

#endif
