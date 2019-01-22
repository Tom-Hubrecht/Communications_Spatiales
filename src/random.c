#include <stdio.h>
#include <string.h>

#include "bit_array.h"
#include "bar.h"

#include "random.h"

generator * initGenerator(void)
{
    bar *generator = barcreate(8);
    barfill(generator);

    return generator;
}


char yield(generator *gen)
{
    char res = barget(gen, 0);
    char a = res;
    a = a ^ barget(gen, 3);
    a = a ^ barget(gen, 5);
    a = a ^ barget(gen, 7);
    barshr(gen, 1, a);

    return res;
}


bar * sequence(generator *gen, uint n)
{
    bar *res = barcreate(n);
    
    for(uint i = 0; i < n; i++)
    {
        barmake(res, i, yield(gen));
    }
    
    return res;
}


bar * combine(generator *gen, bar *message)
{
    resetGenerator(gen);
    int n = barlen(message);
    bar *seq = sequence(gen, n);
    bar *res = barcreate(n);
    barxor(res, seq, message);

    return res;
}


void resetGenerator(generator *gen)
{
    barfill(gen);
}


void freeGenerator(generator *gen)
{
    bardestroy(gen);
}