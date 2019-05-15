#include <stdio.h>
#include <string.h>
#include <math.h>

#include "../lib/BitArray/bit_array.h"
#include "../lib/BitArray/bar.h"

#include "list.h"

#ifndef RANDOM_H
#define RANDOM_H

#define M_PI 3.14159265358979323846	/* pi */
#define inf INFINITY

typedef bar generator;


void permutation(i_list *list_i);
double box_muller(double m, double s);
void addNoise(bar *message, double s, double noisy[]);
double pTransition(double x, size_t d, double s);
double normal(double x, double m, double s);
double mean(double *Z, size_t n);
double variance(double *Z, size_t n);


bar * initGenerator(void);
char yield(generator *gen);
bar * sequence(generator *gen, size_t n);
bar * combine(generator *gen, bar *message);
void resetGenerator(generator *gen);
void freeGenerator(generator *gen);

#endif