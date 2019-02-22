/*
    Implementation of turbocode encoding as defined in CCSDS 131.0-B-3
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "bit_array.h"
#include "bar.h"

#include "random.h"
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


// Decode the rate 1/3 code with constraint length 4
double * decode(double *buf, double s)
{

    uint n = k1 * k2 + 4; // The original message length plus the padding bits
                          // at the end
    double alpha[2 * (n + 1) * 16];         // alpha(i, k, m)
                                                // = alpha[i + 2*m + 32*k]
    double beta[(n + 1) * 16];              // beta(k, m)
                                                // = beta[m + 16*k]
    double gamma[2 * (n + 1) * 16];         // gamma(i, R_k, m', m)
                                                // = gamma[i + 2*m + 32*k]
                                                // We only have one choice for
                                                // m' knowing i and m
    double lambda[2 * (n + 1) * 16];        // lambda(i, k, m)
                                                // = lambda[i + 2*m + 32*k]
    double a[n + 1];                        // Used for normalization

    double *logLikelihood = malloc(sizeof(double) * (n + 1));
    if (logLikelihood == NULL)
    {
        return NULL;
    }

    double tmp[2];

    uint d; // The value of the k-th bit
    uint b; // The value of the k-th encoded bit
    uint m; // The previous state of the register
    uint i;
    double x;
    double y;

    // Compute recursively alpha, gamma, beta and lambda for y1 and y2
    // For rate 1/3, the received bits are x_0, y1_0, y2_0, x_1, ...

    // Initialize alpha
    alpha[0] = 1.0;
    alpha[1] = 1.0;

    // Initialize beta
    beta[16*n] = 1.0;

    for(uint k = 1; k <= n; k++) // k-th bit of the message
    {
        x = buf[3 * (k - 1)];
        y = buf[3 * (k - 1) + 1];

        for(uint S = 0; S < 16; S++) // Register state of the encoder
        {
            // Knowing d_k = i and S_k = S,
                // m = S_{k-1} = S/2 + (S & 8) ^ 8*(i ^ (S & 1))
            // b = i ^ (m & 1) ^ (m & 4)/4 ^ (m & 8)/8

            // d_k = 0
            m = S/2 + (S & 8) ^ 8*(S & 1);
            d = 0;
            b = (m & 1) ^ (m & 4)/4 ^ (m & 8)/8;
            gamma[2*S + 32*k] = pTransition(x, d, s) * pTransition(y, b, s) * 0.5;
            i = 2*m + 32*(k-1);
            alpha[2*S + 32*k] = gamma[2*S + 32*k] * (alpha[i] + alpha[1 + i]);

            // d_k = 1
            m = S/2 + (S & 8) ^ 8*((S & 1) ^ 1);
            d = 1;
            b = 1 - b;
            gamma[1 + 2*S + 32*k] = pTransition(x, d, s) * pTransition(y, b, s) * 0.5;
            i = 2*m + 32*(k-1);
            alpha[1 + 2*S + 32*k] = gamma[1 + 2*S + 32*k] * (alpha[i] + alpha[1 + i]);
            a[k] = a[k] + alpha[2*S + 32*k] + alpha[1 + 2*S + 32*k];
        }

        for(uint S = 0; S < 16; S++)
        {
            // Normalize alpha
            alpha[2*S + 32*k] = alpha[2*S + 32*k] / a[k];
            alpha[1 + 2*S + 32*k] = alpha[1 + 2*S + 32*k] / a[k];
        }
    }

    for(int k = (n - 1); k >= 0; k--)
    {
        for(uint S = 0; S < 16; S++)
        {
            // Knowing d_k = i and S_{k-1} = S,
                // m = S_k = (2*S & 15) + (S & 8)/8 ^ (i ^ (S & 4)/4)

            // d_k = 0
            m = (2*S & 15) + ((S & 8)/8) ^ ((S & 4)/4);
            beta[S + 16*k] = beta[m + 16*(k + 1)] * gamma[2*m + 32*(k + 1)];

            // d_k = 1
            m = (2*S & 15) + ((S & 8)/8) ^ (1 - (S & 4)/4);
            i = 2*m + 32*(k + 1);
            beta[S + 16*k] = beta[S + 16*k] + beta[m + 16*(k + 1)] * gamma[1 + i];
        }

        for(uint S = 0; S < 16; S++)
        {
            // Normalize beta
            beta[S + 16*k] = beta[S + 16*k] / a[k + 1];

            // Compute lambda
            lambda[2*S + 32*k] = alpha[2*S + 32*(k + 1)] * beta[S + 16*k];
            lambda[1 + 2*S + 32*k] = alpha[1 + 2*S + 32*(k + 1)] * beta[S + 16*k];
        }

        tmp[0] = 0.0;
        tmp[1] = 0.0;

        for(uint S = 0; S < 16; S++)
        {
            tmp[0] = tmp[0] + lambda[2*S + 32*k];
            tmp[1] = tmp[1] + lambda[1 + 2*S + 32*k];
        }

        logLikelihood[k] = log(tmp[1] / tmp[0]);

    }

    return logLikelihood;
}
