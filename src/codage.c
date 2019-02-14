#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "bit_array.h"
#include "bar.h"

#include "random.h"
#include "turbocode.h"

int getFileSize(FILE *f)
{
    fseek(f, 0, SEEK_END);
    int size = ftell(f);
    fseek(f, 0, SEEK_SET);

    return size;
}


void padding( bar *message, int k, int n)
{
    int m = barlen(message);

    if (m <= (n - 2 * k))
    {
        barsize(message, n);
        barshl(message, k, 0);
     }
}


bar * bitList(int n, int k)
{
    bar *res = barcreate(k);

    for (int i = 0; i < k; ++i)
    {
        barmake(res, i, (n >> i) & 1);
    }

    return res;
}


bar * fileToBitArray(FILE *fp)
{
    int fileSize = getFileSize(fp);
    unsigned char buffer[fileSize];
    bar *bitFile = barcreate(8 * fileSize); // getFileSize returns the length
                                            // in bytes

    int n;
    n = fread(buffer, fileSize, 1, fp);

    for (int i = 0; i < fileSize; ++i)
    {
        barsw8(bitFile, 8 * i, buffer[i]);
    }

    return bitFile;

}


/**void chunkBitFile(bar *bitFile, int k, int n)
{
    int m = barlen(bitFile);
    int b = n - 2 * k; // Block length
    int
}

**/

bar * loadFile(char *fileName, int k, int n)
{
    char cwd[1024];
    getcwd(cwd, 1025);
    strcat(cwd, "../files/");
    strcat(cwd, fileName);

    FILE *fp;
    fp = fopen(cwd, "rb");

    bar *bitFile = fileToBitArray(fp);

    fclose(fp);

    return bitFile;

}


bar * rsc(bar *mes, int k, int gen)
{
    int n = barlen(mes);
    bar *g = bitList(gen, k - 1);
    bar *A = barcreate(k - 1);
    bar *res = barcreate(n);
    char d;
    char a;
    for (int i = 0; i < n; ++i)
    {
        //d = barget(mes, i);
        a = barget(mes, i);
        for (int j = 0; j < k - 1; ++j)
        {
            a = a ^ (barget(g, j) & barget(A, j));
        }
        barmake(res, i, a);
        barshl(A, 1, a);
    }

    return res;
}


void print( bar *mes)
{
    int n = barlen(mes);
    for (int i = 0; i < n; ++i)
    {
        printf("%d", barget(mes, i));
    }

    printf("\n");

}

int main(void) {

    srand((unsigned int)time(NULL));
    uint x;
    double y;
    double z;

    bar *mes = barcreate(8920);

    barset(mes, 3);

    bar *gen = initGenerator();
    bar *res = combine(gen, encode(mes));

    uint n = barlen(res) - 1;
    double noisy[n + 1];

    addNoise(res, 0.5, noisy);

    double *logLikelihood = decode(noisy, 0.5);

    for(uint i = 0; i < 100; i++)
    {
        x = barget(res, 3 * i);
        y = noisy[3 * i];
        z = logLikelihood[3*i];
        printf("%d\t%f\t%f\n", x, y, z);
    }

/*    loadFile("test.txt", 2, 3);

    padding(mes, 5, 65536);
    bar *res = rsc(mes, 5, 31);

    bardestroy(mes);
    bardestroy(res);


    //print(generator);
    //barshl(mes, 3, 1);
    //print(mes);

    bardestroy(mes);
    bardestroy(res);
*/

    free(logLikelihood);

    return 0;
}