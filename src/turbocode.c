#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include "../lib/BitArray/bit_array.h"
#include "../lib/BitArray/bar.h"
#include "random.h"
#include "turbocode.h"
#include "list.h"
#include "tests.h"


const size_t p[] = {0, 31, 37, 43, 47, 53, 59, 61, 67};


// In the CCSDS Standard, bits are numbered from 1 but we really start at 0
// so we need to replace pi(s) by pi(s + 1) - 1
// Instead, we return pi(s) - 1 and call pi(k + 1)
size_t pi(size_t s)
{
    size_t m = s % 2;
    size_t i = s / (2 * k2);
    size_t j = (s / 2) - i * k2;
    size_t t = (19 * i + 1) % (k1 / 2);
    size_t q = (t % 8) + 1;
    size_t c = (p[q] * j + 21 * m) % k2;

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


bar * initMemState(size_t n)
{
    return barcreate(n);
}


h_list * encode_turbo(h_list *buf)
{
    // We use the rate 1/3 for convenience
    h_list *buf_e = chl(3 * (buf->n + 4), 3 * (buf->n + 4));

    // Create the shift-registers
    bar *registerA = initMemState(4);
    bar *registerB = initMemState(4);

    char a;
    char b;
    char d;

    // Encode the buffer
    for(size_t i = 0; i < buf->n; i++)
    {
        d = buf->list[i];
        buf_e->list[3 * i] = d;

        a = yieldEncode(d, registerA);
        buf_e->list[3 * i + 1] = a;

        d = buf->list[pi(i)];
        b = yieldEncode(d, registerB);
        buf_e->list[3 * i + 2] = b;
    }

    // Clean the registers
    for(size_t i = buf->n; i < (buf->n + 4); i++)
    {
        d = barget(registerA, 3) ^ barget(registerA, 2);
        buf_e->list[3 * i] = d;

        a = yieldEncode(d, registerA);
        buf_e->list[3 * i + 1] = a;

        d = barget(registerB, 3) ^ barget(registerB, 2);
        b = yieldEncode(d, registerB);
        buf_e->list[3 * i + 2] = b;
    }

    bardestroy(registerA);
    bardestroy(registerB);

    return buf_e;
}


// Computes the llr given the arrays X and Y
int decode_part(s_list *X, s_list *Y, s_list *llr, double s)
{
    // Same as in decodeStream
    if (X->n != Y->n || X->n != llr->n)
    {
        return 1;
    }

    size_t n = X->n; // The original message length plus the
                            // padding bits at the end

    s_list *alpha = csl(32*(n + 1), 32*(n + 1));    // alpha(i, k, m)
                                                    // = alpha[i + 2*m + 32*k]

    s_list *beta = csl(16*(n + 1), 16*(n + 1));     // beta(k, m)
                                                    // = beta[m + 16*k]

    s_list *gamma = csl(32*(n + 1), 32*(n + 1));    // gamma(i, R_k, m', m)
                                                    // = gamma[i + 2*m + 32*k]
                                                // We only have one choice for
                                                // m' knowing i and m

    s_list *lambda = csl(32*(n + 1), 32*(n + 1));   // lambda(i, k, m)
                                                    // = lambda[i + 2*m + 32*k]

    s_list *a = csl(n + 1, n + 1);                // Used for normalization

    if (alpha == NULL || beta == NULL || gamma == NULL || lambda == NULL ||
        a == NULL)
    {
        return 1;
    }

    double tmp[2];

    size_t d; // The value of the k-th bit
    size_t b; // The value of the k-th encoded bit
    size_t m; // The previous state of the register
    size_t i;
    double x;
    double y;

    // Compute recursively alpha, gamma, beta and lambda for y1 and y2
    // For rate 1/3, the received bits are x_0, y1_0, y2_0, x_1, ...

    // Initialize alpha
    alpha->list[0] = 1.0;
    alpha->list[1] = 1.0;

    // Initialize the norm of alpha
    a->list[0] = 1.0;
    a->list[n] = 1.0;

    // Initialize beta
    beta->list[16*n] = 1.0;

    for(size_t k = 1; k <= n; k++) // k-th bit of the message
    {
        x = X->list[k - 1];
        y = Y->list[k - 1];

        for(size_t S = 0; S < 16; S++) // Register state of the encoder
        {
            // Knowing d_k = i and S_k = S,
                // m = S_{k-1} = S/2 + (S & 8) ^ 8*(i ^ (S & 1))
                // b = i ^ (m & 1) ^ (m & 4)/4 ^ (m & 8)/8

            // d_k = 0
            d = 0;
            m = S/2 + (S & 8) ^ 8*((S & 1) ^ d);
            b = d ^ (m & 1) ^ (m & 4)/4 ^ (m & 8)/8;
            gamma->list[2*S + 32*k] = pTrans(x, d, s) * pTrans(y, b, s);
            i = 2*m + 32*(k-1);
            alpha->list[2*S + 32*k] = gamma->list[2*S + 32*k] *
                                        (alpha->list[i] + alpha->list[1 + i]);

            // d_k = 1
            d = 1;
            m = S/2 + (S & 8) ^ 8*((S & 1) ^ d);
            b = d ^ (m & 1) ^ (m & 4)/4 ^ (m & 8)/8;
            gamma->list[1 + 2*S + 32*k] = pTrans(x, d, s) * pTrans(y, b, s);
            i = 2*m + 32*(k-1);
            alpha->list[1 + 2*S + 32*k] = gamma->list[1 + 2*S + 32*k] *
                                        (alpha->list[i] + alpha->list[1 + i]);
            a->list[k] += alpha->list[2*S + 32*k] + alpha->list[1 + 2*S + 32*k];
        }

        if (isnan(a->list[k]))
        {
            return 2;
        }

        for(size_t S = 0; S < 16; S++)
        {
            // Normalize alpha
            alpha->list[2*S + 32*k] /=  a->list[k];
            alpha->list[1 + 2*S + 32*k] /= a->list[k];
        }
    }

    // lambda(i, n, 0) = alpha(i, n, 0)
    // lambda(i, n, m) = 0 if m != 0
    lambda->list[32 * n] = alpha->list[32 * n];
    lambda->list[1 + 32 * n] = alpha->list[1 + 32 * n];

    for(int k = (n - 1); k > 0; k--)    // Compute the probabilities beta
    {
        for(size_t S = 0; S < 16; S++)
        {
            // Knowing d_k = i and S_{k-1} = S,
                // m = S_k = (2*S & 15) + (S & 8)/8 ^ (i ^ (S & 4)/4)

            // d_k = 0
            m = (2*S & 15) + ((S & 8)/8) ^ ((S & 4)/4);
            beta->list[S + 16*k] = beta->list[m + 16*(k + 1)] *
                                                gamma->list[2*m + 32*(k + 1)];

            // d_k = 1
            m = (2*S & 15) + ((S & 8)/8) ^ (1 - (S & 4)/4);
            i = 2*m + 32*(k + 1);
            beta->list[S + 16*k] += beta->list[i / 2] * gamma->list[1 + i];
        }

        for(size_t S = 0; S < 16; S++)
        {
            // Normalize beta
            beta->list[S + 16*k] /= a->list[k + 1];

            // Compute lambda
            lambda->list[2*S + 32*k] = alpha->list[2*S + 32*k] *
                                                        beta->list[S + 16*k];
            lambda->list[1 + 2*S + 32*k] = alpha->list[1 + 2*S + 32*k] *
                                                        beta->list[S + 16*k];
        }

        tmp[0] = 0.0;
        tmp[1] = 0.0;

        for(size_t S = 0; S < 16; S++)
        {
            tmp[0] += lambda->list[2*S + 32*k];
            tmp[1] += lambda->list[1 + 2*S + 32*k];
        }

        llr->list[k - 1] = log(tmp[1] / tmp[0]);
    }

    // Free the arrays
    free_s_list(alpha);
    free_s_list(beta);
    free_s_list(gamma);
    free_s_list(lambda);
    free_s_list(a);

    return 0;
}


// If  f = 1 then we must deinterleave the array
h_list * recreate(s_list *mes, char f)
{
    char debug = 0;

    if (f != 0 && f != 1)
    {
        if (debug)
        {
            printf("Incorrect value for f\n");
        }
        return NULL;
    }
    size_t i;

    h_list *res = chl(mes->n, mes->n);
    if (debug)
    {
        printf("Created recipient\n");
        printf("Sizes : %d, %d\n", mes->n, res->n);
    }

    for(size_t k = 0; k < mes->n; k++)
    {
        i = f ? pi(k) : k;
        if (debug)
        {
            //printf("%d, %d\n", i, k);
        }
        res->list[i] = (mes->list[k] > 0);
    }
    if (debug)
    {
        printf("Copied data\n");
    }
    return res;
}


// Divides the rate 1/3 stream into three arrays
int split_s(s_list *buf, s_list *X, s_list *Y1, s_list *Y2)
{
    if (buf->n % 3 || X->n != Y1->n || X->n != Y2->n || X->n != (buf->n / 3))
    {
        return 1;
    }

    for(size_t i = 0; i < X->n; i++)
    {
        X->list[i] = buf->list[3 * i];
        Y1->list[i] = buf->list[3*i + 1];
        Y2->list[i] = buf->list[3*i + 2];
    }
    return 0;
}


// Interleave the array Y and puts the result in X
int interleave(s_list *X, s_list *Y)
{
    if (Y->n > X->m_s)
    {
        return 1;
    }

    X->n = Y->n;
    double t;

    for(size_t i = 0; i < k1 * k2; i++)
    {
        // If the result is inf or -inf we replace it to prevent errors
        t = Y->list[pi(i)];
        if (isfinite(t))
        {
            X->list[i] = t;
        }
        if (isinf(t))
        {
            X->list[i] = 100 * (t == inf ? 1.0 : -1.0);
        }
    }
    return 0;
}


// Deinterleave the array Y and puts the result in X
int deinterleave(s_list *X, s_list *Y)
{
    if (Y->n > X->m_s)
    {
        return 1;
    }

    double t;

    for(size_t i = 0; i < k1 * k2; i++)
    {
        // If the result is inf or -inf we replace it to prevent errors
        t = Y->list[i];
        if (isfinite(t))
        {
            X->list[pi(i)] = t;
        }
        if (isinf(t))
        {
            X->list[pi(i)] = 100 * (t == inf ? 1.0 : -1.0);
        }
    }
    return 0;
}


// Compute the max negative and min positive number of a s_list
void min_max(double mM[2], s_list *X)
{
    double x;
    double m = -inf;
    double M = inf;
    for(size_t i = 0; i < X->n; i++)
    {
        x = X->list[i];
        if (x <= 0 && x > m)
        {
            m = x;
        }

        if (x >= 0 && x < M)
        {
            M = x;
        }
    }
    mM[0] = m;
    mM[1] = M;
}


// Compute the min negative and max positive number of an array
void max_min(double mM[2], s_list *X)
{
    double x;
    double m = 0;
    double M = 0;
    for(size_t i = 0; i < X->n; i++)
    {
        x = X->list[i];
        if (x < m)
        {
            m = x;
        }

        if (x > M)
        {
            M = x;
        }
    }
    mM[0] = m;
    mM[1] = M;
}


// Executes a single pass through both decoders on the stream
h_list * decode_turbo_basic(s_list *buf, double s)
{
    char debug = 0;

    char interleaved = 0;
    size_t n = buf->n / 3;

    s_list *llr = csl(n, n);    // The 0-th bit is not considered
    s_list *X1 = csl(n, n);
    s_list *X2 = csl(n, buf->n);
    s_list *Y1 = csl(n, n);
    s_list *Y2 = csl(n, n);

    if (debug)
    {
        printf("Created the s_lists\n");
        printf("Sizes :\n");
        printf("\t- buf : %d\n", buf->n);
        printf("\t- llr : %d\n", llr->n);
        printf("\t- X1 : %d\n", X1->n);
        printf("\t- X2 : %d\n", X2->n);
        printf("\t- Y1 : %d\n", Y1->n);
        printf("\t- Y2 : %d\n", Y2->n);
    }

    double t;
    int k;

    split_s(buf, X1, Y1, Y2);
    if (debug)
    {
        printf("Split incomming message\n");
    }

    k = decode_part(X1, Y1, llr, s);
    if (debug)
    {
        printf("First pass on the decoder\n");
        printf("\t- Error code : %d\n", k);
    }

    // If the values in llr are sufficiently big, there is no need to do a
    // second pass, thus we compute the min of positives and max of negatives
    // to check this situation
    double mM[2];
    double Mm[2];
    min_max(mM, llr);
    max_min(Mm, llr);

    if (Mm[0] > -26 || Mm[1] < 26)
    {
        interleaved = 1;

        // We need to interleave the llr to match the pattern of Y2
        interleave(X2, llr);
        if (debug)
        {
            printf("Interleaved the data\n");
        }

        k = decode_part(X2, Y2, llr, s);
        if (debug)
        {
            printf("Second pass on the decoder\n");
            printf("\t- Error code : %d\n", k);
        }
    }



    h_list *res = recreate(llr, interleaved);

    // Free all the arrays
    free_s_list(llr);
    free_s_list(X1);
    free_s_list(X2);
    free_s_list(Y1);
    free_s_list(Y2);
    if (debug)
    {
        printf("Freed structures\n");
    }

    return res;
}


// Executes at most i_max passes through both decoders on the stream
h_list * decode_turbo_iter(s_list *buf, double s, size_t i_max)
{
    char debug = 1;

    char done = 0;
    char interleaved = 0;
    size_t iter = 0;

    double mM[2];
    double Mm[2];

    size_t n = buf->n / 3;

    s_list *llr = csl(n, n);    // The 0-th bit is not considered
    s_list *X1 = csl(n, n);
    s_list *X2 = csl(n, buf->n);
    s_list *Y1 = csl(n, n);
    s_list *Y2 = csl(n, n);

    if (debug)
    {
        printf("Created the s_lists\n");
        printf("Sizes :\n");
        printf("\t- buf : %d\n", buf->n);
        printf("\t- llr : %d\n", llr->n);
        printf("\t- X1 : %d\n", X1->n);
        printf("\t- X2 : %d\n", X2->n);
        printf("\t- Y1 : %d\n", Y1->n);
        printf("\t- Y2 : %d\n", Y2->n);
    }

    double t;
    int k;

    split_s(buf, X1, Y1, Y2);
    if (debug)
    {
        printf("Split incomming message\n");
    }

    while (!done && iter < i_max)
    {
        iter ++;

        if (debug)
        {
            printf("Iteration %d / %d\n", iter, i_max);
        }

        interleaved = 0;

        k = decode_part(X1, Y1, llr, s);
        substract_s(llr, X1);
        if (debug)
        {
            printf("\tFirst pass on the decoder\n");
            printf("\t\t- Error code : %d\n", k);
        }

        // If the values in llr are sufficiently big, there is no need to do a
        // second pass, thus we compute the min of positives and max of
        // negatives to check for this situation
        max_min(Mm, llr);

        if (Mm[0] > -5000000 || Mm[1] < 5000000)
        {
            interleaved = 1;

            // We need to interleave the llr to match the pattern of Y2
            interleave(X2, llr);
            if (debug)
            {
                printf("\tInterleaved the data\n");
            }

            k = decode_part(X2, Y2, llr, s);
            substract_s(llr, X2);
            if (debug)
            {
                printf("\tSecond pass on the decoder\n");
                printf("\t\t- Error code : %d\n", k);
            }

            // Update X1 with the new values
            deinterleave(X1, llr);
        }

        if (Mm[0] < -5000000 && Mm[1] > 5000000)
        {
            done = 1;
            if (debug)
            {
                printf("Decoding complete\n");
            }
        }
    }

    h_list *res = recreate(llr, interleaved);

    // Free all the arrays
    free_s_list(llr);
    free_s_list(X1);
    free_s_list(X2);
    free_s_list(Y1);
    free_s_list(Y2);
    if (debug)
    {
        printf("Freed structures\n");
    }

    return res;
}
