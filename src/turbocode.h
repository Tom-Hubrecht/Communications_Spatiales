#include <stdio.h>
#include <string.h>

#include "bit_array.h"
#include "bar.h"

typedef unsigned int uint;

#ifndef TURBOCODE_H
#define TURBOCODE_H

#define k1 8
#define k2 (223 * 5)

#define pTrans pTransition

extern const uint p[];

typedef bar buffer;

uint pi(uint s);
char yieldEncode(char d, bar *memState);
bar * initMemState(uint n);
bar * encode_turbo(buffer *buf);
double * decode(double *X, double *Y, double s);
int decodePart(double *X, double *Y, double *llr, double s);
bar * recreate(double *mes, char f);
bar * decodeStreamOnce(double *buf, double s);
bar * decodeStreamParallel(double *buf, double s, size_t q);
void split(double *buf, double *X, double *Y1, double *Y2, size_t n);
int difference(bar *m1, bar *m2);
void interleave(double *X, double *Y);
void deinterleave(double *X, double *Y);

#endif
