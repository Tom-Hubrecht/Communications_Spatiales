#include "tests.h"
#include "turbocode.h"
#include "random.h"


int interleaver(void)
{
    int count[k1 * k2] = {0};
    int d;

    for(size_t i = 0; i < k1*k2; i++)
    {
        count[pi(i)] += 1;
    }

    for(size_t i = 0; i < k1*k2; i++)
    {
        if (count[i] != 1)
        {
            d += 1;
            printf("%d\n", count[i]);
        }
    }

    return d;
}

int main(void)
{
    int i = interleaver();
    if (i == 0)
    {
        printf("Interleaver function :\tPass\n");
    }
    else
    {
        printf("Interleaver function :\tFail\n");
    }

}