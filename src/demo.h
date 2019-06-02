#include <stdio.h>

#include "list.h"
#include "ldpc.h"
#include "turbocode.h"

#ifndef DEMO_H
#define DEMO_H

int demo_turbo_basic(char *file_name, double s);
int demo_turbo_iter(char *file_name, double s, size_t i_max);

int demo_ldpc_basic(char *file_name, double s, size_t i_max);
int demo_ldpc_proba(char *file_name, double s, size_t i_max);

int demo_turbo_graph(size_t nb_iter, double p);

#endif
