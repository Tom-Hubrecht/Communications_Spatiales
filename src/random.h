#include <stdio.h>
#include <string.h>

#include "bit_array.h"
#include "bar.h"

#ifndef RANDOM_H
#define RANDOM_H

#define M_PI 3.14159265358979323846	/* pi */
#define inf INFINITY

typedef unsigned int uint;
typedef bar generator;


double box_muller(double m, double s);
void addNoise(bar *message, double s, double noisy[]);
double pTransition(double x, uint d, double s);
double normal(double x, double m, double s);
double mean(double *Z, size_t n);
double variance(double *Z, size_t n);


bar * initGenerator(void);
char yield(generator *gen);
bar * sequence(generator *gen, uint n);
bar * combine(generator *gen, bar *message);
void resetGenerator(generator *gen);
void freeGenerator(generator *gen);

#endif