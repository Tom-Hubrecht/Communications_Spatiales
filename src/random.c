#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "bit_array.h"
#include "bar.h"

#include "random.h"


// Add a White Gaussian Noise with mean 0 and standard deviation s to encoded data
void addNoise(bar *message, double s, double noisy[])
{
    uint n = barlen(message);

    for(uint i = 0; i < n; i++)
    {
        noisy[i] = (2 * barget(message, i) - 1.0) + box_muller(0.0, s);
    }
}


// Compute the transition probability of the channel
double pTransition(double x, uint d, double s)
{
    char mu = 2 * d - 1; // -1 if d = 0, 1 if d = 1

    return exp(- pow((x - mu) / s, 2) / 2.0) / (s * sqrt(2.0 * M_PI));
}


// Create the pseudorandom generator
generator * initGenerator(void)
{
    bar *generator = barcreate(8);
    barfill(generator);

    return generator;
}


// Return a random value and update the memory state
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


// Generate a sequence of n random values
bar * sequence(generator *gen, uint n)
{
    bar *res = barcreate(n);

    for(uint i = 0; i < n; i++)
    {
        barmake(res, i, yield(gen));
    }

    return res;
}


// XOR a bit array with a random sequence
bar * combine(generator *gen, bar *message)
{
    resetGenerator(gen);
    int n = barlen(message);
    bar *seq = sequence(gen, n);
    bar *res = barcreate(n);
    barxor(res, seq, message);

    return res;
}


// Reset the memory state of the generator
void resetGenerator(generator *gen)
{
    barfill(gen);
}


void freeGenerator(generator *gen)
{
    bardestroy(gen);
}