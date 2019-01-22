#include <stdio.h>
#include <string.h>

#include "bit_array.h"
#include "bar.h"

#ifndef RANDOM_H
#define RANDOM_H

typedef unsigned int uint;
typedef bar generator;


bar * initGenerator(void);
char yield(generator *gen);
bar * sequence(generator *gen, uint n);
bar * combine(generator *gen, bar *message);
void resetGenerator(generator *gen);
void freeGenerator(generator *gen);

#endif