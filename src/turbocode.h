#include <stdio.h>
#include <string.h>

#include "bit_array.h"
#include "bar.h"

typedef unsigned int uint;

#ifndef TURBOCODE_H
#define TURBOCODE_H

#define k1 8
#define k2 (223 * 5)

extern const uint p[];

typedef bar buffer;

uint pi(uint s);
char yieldEncode(char d, bar *memState);
bar * initMemState(uint n);
bar * encode(buffer *buf);
#endif

