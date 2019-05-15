#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "random.h"
#include "turbocode.h"
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


bar * encode_turbo(buffer *buf)
{
    size_t k = barlen(buf);
    // We use the rate 1/3 for convenience
    bar *encodedBuffer = barcreate(3 * (k + 4));

    // Create the shift-registers
    bar *registerA = initMemState(4);
    bar *registerB = initMemState(4);

    char a;
    char b;
    char d;

    // Encode the buffer
    for(size_t i = 0; i < k; i++)
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
    for(size_t i = k; i < (k + 4); i++)
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
int decodePart(double *X, double *Y, double *llr, double s)
{
    // Same as in decodeStream
    size_t t = sizeof(double);
    size_t n = k1 * k2 + 4; // The original message length plus the
                            // padding bits at the end

    double *alpha = calloc(2 * (n + 1) * 16, t);   // alpha(i, k, m)
                                                // = alpha[i + 2*m + 32*k]

    double *beta = calloc((n + 1) * 16, t);        // beta(k, m)
                                                // = beta[m + 16*k]

    double *gamma = calloc(2 * (n + 1) * 16, t);   // gamma(i, R_k, m', m)
                                                // = gamma[i + 2*m + 32*k]
                                                // We only have one choice for
                                                // m' knowing i and m

    double *lambda = calloc(2 * (n + 1) * 16, t);  // lambda(i, k, m)
                                                // = lambda[i + 2*m + 32*k]

    double *a = calloc((n + 1), t);                // Used for normalization

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
    alpha[0] = 1.0;
    alpha[1] = 1.0;

    // Initialize the norm of alpha
    a[0] = 1.0;
    a[n] = 1.0;

    // Initialize beta
    beta[16*n] = 1.0;

    for(size_t k = 1; k <= n; k++) // k-th bit of the message
    {
        x = X[k - 1];
        y = Y[k - 1];

        for(size_t S = 0; S < 16; S++) // Register state of the encoder
        {
            // Knowing d_k = i and S_k = S,
                // m = S_{k-1} = S/2 + (S & 8) ^ 8*(i ^ (S & 1))
                // b = i ^ (m & 1) ^ (m & 4)/4 ^ (m & 8)/8

            // d_k = 0
            d = 0;
            m = S/2 + (S & 8) ^ 8*((S & 1) ^ d);
            b = d ^ (m & 1) ^ (m & 4)/4 ^ (m & 8)/8;
            gamma[2*S + 32*k] = pTrans(x, d, s) * pTrans(y, b, s);
            i = 2*m + 32*(k-1);
            alpha[2*S + 32*k] = gamma[2*S + 32*k] * (alpha[i] + alpha[1 + i]);

            // d_k = 1
            d = 1;
            m = S/2 + (S & 8) ^ 8*((S & 1) ^ d);
            b = d ^ (m & 1) ^ (m & 4)/4 ^ (m & 8)/8;
            gamma[1 + 2*S + 32*k] = pTrans(x, d, s) * pTrans(y, b, s);
            i = 2*m + 32*(k-1);
            alpha[1 + 2*S + 32*k] = gamma[1 + 2*S + 32*k] *
                                        (alpha[i] + alpha[1 + i]);
            a[k] += alpha[2*S + 32*k] + alpha[1 + 2*S + 32*k];
        }

        if (isnan(a[k]))
        {
            return 2;
        }

        for(size_t S = 0; S < 16; S++)
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
        for(size_t S = 0; S < 16; S++)
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

        for(size_t S = 0; S < 16; S++)
        {
            // Normalize beta
            beta[S + 16*k] /= a[k + 1];

            // Compute lambda
            lambda[2*S + 32*k] = alpha[2*S + 32*k] * beta[S + 16*k];
            lambda[1 + 2*S + 32*k] = alpha[1 + 2*S + 32*k] * beta[S + 16*k];
        }

        tmp[0] = 0.0;
        tmp[1] = 0.0;

        for(size_t S = 0; S < 16; S++)
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

    return 0;
}


// Computes the llr given the arrays X and Y
int decodePartPar(double *X, double *Y, double *Z, double *llr, double s)
{
    // Same as in decodeStream
    size_t t = sizeof(double);
    size_t n = k1 * k2 + 4; // The original message length plus the
                            // padding bits at the end

    double *alpha = calloc(2 * (n + 1) * 16, t);   // alpha(i, k, m)
                                                // = alpha[i + 2*m + 32*k]

    double *beta = calloc((n + 1) * 16, t);        // beta(k, m)
                                                // = beta[m + 16*k]

    double *gamma = calloc(2 * (n + 1) * 16, t);   // gamma(i, R_k, m', m)
                                                // = gamma[i + 2*m + 32*k]
                                                // We only have one choice for
                                                // m' knowing i and m

    double *lambda = calloc(2 * (n + 1) * 16, t);  // lambda(i, k, m)
                                                // = lambda[i + 2*m + 32*k]

    double *a = calloc((n + 1), t);                // Used for normalization

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
    alpha[0] = 1.0;
    alpha[1] = 1.0;

    // Initialize the norm of alpha
    a[0] = 1.0;
    a[n] = 1.0;

    // Initialize beta
    beta[16*n] = 1.0;

    for(size_t k = 1; k <= n; k++) // k-th bit of the message
    {
        x = X[k - 1];
        y = Y[k - 1];

        for(size_t S = 0; S < 16; S++) // Register state of the encoder
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
            alpha[1 + 2*S + 32*k] = gamma[1 + 2*S + 32*k] *
                                        (alpha[i] + alpha[1 + i]);
            a[k] += alpha[2*S + 32*k] + alpha[1 + 2*S + 32*k];
        }

        if (isnan(a[k]))
        {
            return 2;
        }

        for(size_t S = 0; S < 16; S++)
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
        for(size_t S = 0; S < 16; S++)
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

        for(size_t S = 0; S < 16; S++)
        {
            // Normalize beta
            beta[S + 16*k] /= a[k + 1];

            // Compute lambda
            lambda[2*S + 32*k] = alpha[2*S + 32*k] * beta[S + 16*k];
            lambda[1 + 2*S + 32*k] = alpha[1 + 2*S + 32*k] * beta[S + 16*k];
        }

        tmp[0] = 0.0;
        tmp[1] = 0.0;

        for(size_t S = 0; S < 16; S++)
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

    return 0;
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


// Interleave the array Y and puts the result in X
void interleave(double *X, double *Y)
{
    double t;

    for(size_t i = 0; i < k1 * k2; i++)
    {
        // If the result is inf or -inf we replace it to prevent errors
        t = Y[pi(i)];
        if (isfinite(t))
        {
            X[i] = t;
        }
        if (isinf(t))
        {
            X[i] = t == inf ? 1.0 : -1.0;
        }
    }
}


// Deinterleave the array Y and puts the result in X
void deinterleave(double *X, double *Y)
{
    double t;

    for(size_t i = 0; i < k1 * k2; i++)
    {
        // If the result is inf or -inf we replace it to prevent errors
        t = Y[i];
        if (isfinite(t))
        {
            X[pi(i)] = t;
        }
        if (isinf(t))
        {
            X[pi(i)] = t == inf ? 1.0 : -1.0;
        }
    }
}


// Compute the max negative and min positive number of an array
void minMax(double *mM, double *X, size_t n)
{
    double x;
    double m = -inf;
    double M = inf;
    for(size_t i = 0; i < n; i++)
    {
        x = X[i];
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
void maxMin(double *mM, double *X, size_t n)
{
    double x;
    double m = 0;
    double M = 0;
    for(size_t i = 0; i < n; i++)
    {
        x = X[i];
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
bar * decodeStreamOnce(double *buf, double s)
{
    size_t d = sizeof(double);
    size_t n = k1 * k2 + 4; // The original message length plus the
                            // padding bits at the end

    double *llr = calloc(n, d);    // The 0-th bit is not considered

    double *X1 = calloc(n, d);
    double *X2 = calloc(n, d);
    double *Y1 = calloc(n, d);
    double *Y2 = calloc(n, d);

    double t;
    int k;

    split(buf, X1, Y1, Y2, n);

    decodePart(X1, Y1, llr, s);

    // If the values in llr are sufficiently big, there is no need to do a
    // second pass, thus we compute the min of positives and max of negatives
    // to check this situation
    double mM[2];
    double Mm[2];
    minMax(mM, llr, k1 * k2);
    maxMin(Mm, llr, k1 * k2);

    if (Mm[0] < -26 && Mm[1] > 26)
    {
        return recreate(llr, 0);
    }

    // We need to interleave the llr to match the pattern of Y2
    interleave(X2, llr);

    for(size_t i = 0; i < 4; i++)
    {
        printf("%f\t%f\n", llr[n - 5 + i], Y1[n - 4 + i]);
    }

    k = decodePart(X2, Y2, llr, s);

    printf("%d\n", k);


    // Free all the arrays
    free(llr);

    free(X1);
    free(X2);
    free(Y1);
    free(Y2);

    return recreate(llr, 1);
}


// Decode iteratively the stream
bar * decodeStreamParallel(double *buf, double s, size_t q)
{
    size_t d = sizeof(double);
    size_t n = k1 * k2 + 4; // The original message length plus the
                            // padding bits at the end

    double *llr = calloc(n, d);    // The 0-th bit is not considered

    double *X1 = calloc(n, d);
    double *X2 = calloc(n, d);
    double *Y1 = calloc(n, d);
    double *Y2 = calloc(n, d);
    double *Z  = calloc(n, d);

    split(buf, X1, Y1, Y2, n);

    for(size_t i = 0; i < q; i++)
    {
        decodePart(X1, Y1, llr, s);

        // We need to interleave the llr to match the pattern of Y2
        interleave(X2, llr);

        decodePart(X2, Y2, llr, s);

        deinterleave(X1, llr);

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