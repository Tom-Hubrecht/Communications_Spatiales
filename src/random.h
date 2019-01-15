#include <stdio.h>
#include <string.h>

#include "bit_array.h"
#include "bar.h"

typedef unsigned int uint;

bar * initGenerator();
char yield(bar *generator);
bar * sequence(bar *generator, uint n);
bar * combine(bar *generator, bar *message);
void resetGenerator(bar *generator);
void freeGenerator(bar *generator);
