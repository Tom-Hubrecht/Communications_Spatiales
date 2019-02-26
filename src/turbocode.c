#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "random.h"
#include "turbocode.h"
#include "tests.h"


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

// In the CCSDS Standard, bits are numbered from 1 but we really start at 0
// so we need to replace pi(s) by pi(s + 1) - 1
// Instead, we return pi(s) - 1 and call pi(k + 1)
uint pi(uint s)
{
    uint m = s % 2;
    uint i = s / (2 * k2);
    uint j = (s / 2) - i * k2;
    uint t = (19 * i + 1) % (k1 / 2);
    uint q = (t % 8) + 1;
    uint c = (p[q] * j + 21 * m) % k2;

    return 2 * (t + c * (k1 / 2) + 1) - m - 1;
}


// We use the connection vector G1 = 11011 for the message
// and G0 = 10011 for the component code
char yieldEncode(char d, bar *m)
{
    char a = d ^ barget(m, 0) ^ barget(m, 2) ^ barget(m, 3);
    char g = d ^ barget(m, 2) ^ barget(m, 3);

    barshl(m, 1, g);

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


// Computes the llr given the arrays X and Y
void decodePart(double *X, double *Y, double *llr, double s)
{
    // Same as in decodeStream
    size_t t = sizeof(double);
    size_t n = k1 * k2 + 4; // The original message length plus the
                            // padding bits at the end

    double *alpha = malloc(2 * (n + 1) * 16 * t);   // alpha(i, k, m)
                                                // = alpha[i + 2*m + 32*k]

    double *beta = malloc((n + 1) * 16 * t);        // beta(k, m)
                                                // = beta[m + 16*k]

    double *gamma = malloc(2 * (n + 1) * 16 * t);   // gamma(i, R_k, m', m)
                                                // = gamma[i + 2*m + 32*k]
                                                // We only have one choice for
                                                // m' knowing i and m

    double *lambda = malloc(2 * (n + 1) * 16 * t);  // lambda(i, k, m)
                                                // = lambda[i + 2*m + 32*k]

    double *a = malloc((n + 1) * t);                // Used for normalization

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

    // Initialize the norm of alpha
    a[0] = 1.0;

    // Initialize beta
    beta[16*n] = 1.0;

    for(uint k = 1; k <= n; k++) // k-th bit of the message
    {
        x = X[k - 1];
        y = Y[k - 1];

        for(uint S = 0; S < 16; S++) // Register state of the encoder
        {
            // Knowing d_k = i and S_k = S,
                // m = S_{k-1} = S/2 + (S & 8) ^ 8*(i ^ (S & 1))
                // b = i ^ (m & 1) ^ (m & 4)/4 ^ (m & 8)/8

            // d_k = 0
            d = 0;
            m = S/2 + (S & 8) ^ 8*((S & 1) ^ d);
            b = d ^ (m & 1) ^ (m & 4)/4 ^ (m & 8)/8;
            gamma[2*S + 32*k] = pTrans(x, d, s) * pTrans(y, b, s) * 0.5;
            i = 2*m + 32*(k-1);
            alpha[2*S + 32*k] = gamma[2*S + 32*k] * (alpha[i] + alpha[1 + i]);

            // d_k = 1
            d = 1;
            m = S/2 + (S & 8) ^ 8*((S & 1) ^ d);
            b = d ^ (m & 1) ^ (m & 4)/4 ^ (m & 8)/8;
            gamma[1 + 2*S + 32*k] = pTrans(x, d, s) * pTrans(y, b, s) * 0.5;
            i = 2*m + 32*(k-1);
            alpha[1 + 2*S + 32*k] = gamma[1 + 2*S + 32*k] * (alpha[i] + alpha[1 + i]);
            a[k] += alpha[2*S + 32*k] + alpha[1 + 2*S + 32*k];
        }

        for(uint S = 0; S < 16; S++)
        {
            // Normalize alpha
            alpha[2*S + 32*k] /=  a[k];
            alpha[1 + 2*S + 32*k] /= a[k];
        }
    }

    // lambda(i, n, 0) = alpha(i, n, 0)
    // lambda(i, n, m) = 0 if m != 0
    lambda[32 * n] = alpha[32 * n];
    lambda[1 + 32 * n] = alpha [1 + 32 * n];

    for(int k = (n - 1); k > 0; k--)    // Compute the probabilities beta
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
            beta[S + 16*k] += beta[i / 2] * gamma[1 + i];
        }

        for(uint S = 0; S < 16; S++)
        {
            // Normalize beta
            beta[S + 16*k] /= a[k + 1];

            // Compute lambda
            lambda[2*S + 32*k] = alpha[2*S + 32*k] * beta[S + 16*k];
            lambda[1 + 2*S + 32*k] = alpha[1 + 2*S + 32*k] * beta[S + 16*k];
        }

        tmp[0] = 0.0;
        tmp[1] = 0.0;

        for(uint S = 0; S < 16; S++)
        {
            tmp[0] += lambda[2*S + 32*k];
            tmp[1] += lambda[1 + 2*S + 32*k];
        }

        llr[k - 1] = log(tmp[1] / tmp[0]);
    }

    // Free the arrays
    free(alpha);
    free(beta);
    free(gamma);
    free(lambda);
    free(a);
}


// Computes the number of bits that differ from both arrays
int difference(bar *m1, bar *m2)
{
    int n = barlen(m1);
    if (n != barlen(m2))
    {
        return -1;
    }

    int d = 0;
    for(int k = 0; k < n; k++)
    {
        d = d + (barget(m1, k) != barget(m2, k));
    }

    return d;
}


// If  f = 1 then we must deinterleave the array
bar * recreate(double *mes, char f)
{
    if (f != 0 && f != 1)
    {
        return NULL;
    }
    size_t i;

    bar *res = barcreate(k1 * k2);
    for(size_t k = 0; k < k1 * k2; k++)
    {
        i = f ? pi(k) : k;
        barmake(res, i, (mes[k] > 0));
    }

    return res;
}


// Divides the rate 1/3 stream into three arrays
void split(double *buf, double *X, double *Y1, double *Y2, size_t n)
{
    for(size_t i = 0; i < n; i++)
    {
        X[i] = buf[3 * i];
        Y1[i] = buf[3*i + 1];
        Y2[i] = buf[3*i + 2];
    }
}


// Executes a single pass through both decoders on the stream
bar * decodeStreamOnce(double *buf, double s)
{
    size_t d = sizeof(double);
    size_t n = k1 * k2 + 4; // The original message length plus the
                            // padding bits at the end

    double *llr = malloc(n * d);    // The 0-th bit is not considered

    double *X1 = malloc(n * d);
    double *X2 = malloc(n * d);
    double *Y1 = malloc(n * d);
    double *Y2 = malloc(n * d);

    split(buf, X1, Y1, Y2, n);
    decodePart(X1, Y1, llr, s);

    // We need to interleave the llr to match the pattern of Y2
    for(size_t i = 0; i < k1 * k2; i++)
    {
        X2[i] = llr[pi(i)];
    }

    decodePart(X2, Y2, llr, s);

    bar *res = recreate(llr, 1);

    // Free all the arrays
    free(llr);

    free(X1);
    free(X2);
    free(Y1);
    free(Y2);

    return res;
}


bar * decodeStreamParallel(double *buf, double s, size_t q)
{
    size_t d = sizeof(double);
    size_t n = k1 * k2 + 4; // The original message length plus the
                            // padding bits at the end

    double *llr = malloc(n * d);    // The 0-th bit is not considered

    double *X1 = malloc(n * d);
    double *X2 = malloc(n * d);
    double *Y1 = malloc(n * d);
    double *Y2 = malloc(n * d);
    double *Z  = malloc(n * d);

    split(buf, X1, Y1, Y2, n);

    for(size_t i = 0; i < q; i++)
    {
        decodePart(X1, Y1, llr, s);

        // We need to interleave the llr to match the pattern of Y2
        for(size_t i = 0; i < k1 * k2; i++)
        {
            X2[i] = llr[pi(i)];
        }

        decodePart(X2, Y2, llr, s);

    }


    bar *res = recreate(llr, 1);

    // Free all the arrays
    free(llr);

    free(X1);
    free(X2);
    free(Y1);
    free(Y2);

    return res;
}