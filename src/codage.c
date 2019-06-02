#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "../lib/BitArray/bit_array.h"
#include "../lib/BitArray/bar.h"

#include "turbocode.h"
#include "ldpc.h"
#include "list.h"
#include "random.h"
#include "basic.h"
#include "demo.h"


int main(int argc, char *argv[])
{
    // We only accept at most three arguments
    if (argc == 1 || argc > 4)
    {
        return 1;
    }

    int err = 0;
    int demos = 0;
    double s = 0.0;

    // Incorrect type for an argument
    sscanf(argv[1], "%d", &demos);
    sscanf(argv[3], "%lf", &s);

    char *file = argv[2];

    // Initialise the random generator
    srand((unsigned int)time(NULL));

    if (demos & 1)
    {
        demo_turbo_iter(file, s, 10);
    }
    if (demos & 2)
    {
        demo_ldpc_basic(file, s, 500);
    }
    if (demos & 4)
    {
        demo_ldpc_proba(file, s, 10);
    }
    if (demos & 8)
    {
        demo_turbo_graph(20, s);
    }


    return 0;
}
