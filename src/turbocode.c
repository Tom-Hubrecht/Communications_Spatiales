/*
    Implementation of turbocode encoding as defined in CCSDS 131.0-B-3
*/

#include <stdio.h>
#include <string.h>

#include "bit_array.h"
#include "bar.h"

#include "turbocode.h"

const uint p[] = {0, 31, 37, 43, 47, 53, 59, 61, 67};


/* Tests functions */

void print2( bar *mes)
{
    int n = barlen(mes);
    for (int i = 0; i < n; ++i)
    {
        printf("%d", barget(mes, i));
    }

    printf("\n");

}

/*-----------------*/

uint pi(uint s)
{
    uint m = (s - 1) % 2;
    uint i = (s - 1) / (2 * k2);
    uint j = ((s - 1) / 2) - i * k2;
    uint t = (19 * i + 1) % (k1 / 2);
    uint q = (t % 8) + 1;
    uint c = (p[q] * j + 21 * m) % k2;

    return 2 * (t + c * (k1 / 2) + 1) - m;
}

// We use the connection vector G1 = 11011 for the message
// and G0 = 10011 for the component code

char yieldEncode(char d, bar *memState)
{
    char a = d;
    char g = d;

    // Encode the bit d for output
    a = a ^ barget(memState, 0);
    a = a ^ barget(memState, 2);
    a = a ^ barget(memState, 3);

    // Encode the bit d to update the memState
    g = a ^ barget(memState, 0);

    barshr(memState, 1, g);

    return a;
}


bar * initMemState(uint n)
{
    return barcreate(n);
}


bar * encode(buffer *buf)
{
    uint k = barlen(buf);
    // We use the rate 1/3 for convenience
    bar *encodedBuffer = barcreate(3 * (k + 4));

    // Create the shift-registers
    bar *registerA = initMemState(4);
    bar *registerB = initMemState(4);

    char a;
    char b;
    char d;

    // Encode the buffer
    for(uint i = 0; i < k; i++)
    {
        d = barget(buf, i);
        barmake(encodedBuffer, 3 * i, d);

        a = yieldEncode(d, registerA);
        barmake(encodedBuffer, 3 * i + 1, a);

        d = barget(buf, pi(i));
        b = yieldEncode(d, registerB);
        barmake(encodedBuffer, 3 * i + 2, b);
    }

    // Clean the registers
    for(uint i = k; i < (k + 4); i++)
    {
        d = barget(registerA, 3) ^ barget(registerA, 2);
        barmake(encodedBuffer, 3 * i, d);

        a = yieldEncode(d, registerA);
        barmake(encodedBuffer, 3 * i + 1, a);

        d = barget(registerB, 3) ^ barget(registerB, 2);
        b = yieldEncode(d, registerB);
        barmake(encodedBuffer, 3 * i + 2, b);
    }
    
    bardestroy(registerA);
    bardestroy(registerB);

    return encodedBuffer;
}