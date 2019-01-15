#include <stdio.h>
#include <string.h>

#include "bit_array.h"
#include "bar.h"

#include "random.h"

typedef unsigned int uint;

bar * initGenerator()
{
    bar *generator = barcreate(8);
    barfill(generator);

    return generator;
}


char yield(bar *generator)
{
    char res = barget(generator, 0);
    char a = res;
    a = a ^ barget(generator, 3);
    a = a ^ barget(generator, 5);
    a = a ^ barget(generator, 7);
    barshr(generator, 1, a);

    return res;
}


bar * sequence(bar *generator, uint n)
{
    bar *res = barcreate(n);
    
    for(uint i = 0; i < n; i++)
    {
        barmake(res, i, yield(generator));
    }
    
    return res;
}


bar * combine(bar *generator, bar *message)
{
    resetGenerator(generator);
    int n = barlen(message);
    bar *seq = sequence(generator, n);
    bar *res = barcreate(n);
    barxor(res, seq, message);

    return res;
}


void resetGenerator(bar *generator)
{
    barfill(generator);
}


void freeGenerator(bar *generator)
{
    bardestroy(generator);
}