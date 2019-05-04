#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "bit_array.h"
#include "bar.h"

#include "random.h"
#include "turbocode.h"
#include "list.h"

#include "ldpc.h"

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


int main(void)
{
    // Initialise the random generator
    srand((unsigned int)time(NULL));
    uint x;
    double y;
    double z;

//    double s = 0.95;
//    printf("Noise sigma : %f\n", s);

    h_matrix *test = cldpc(20, 3, 4);
    h_matrix *gen_h = cgm_h(test);
    a_matrix *gen_a = convert_h(gen_h);
    h_list *mes = chl(20, 20);


    set_all_h_list(mes, 1);
    print_h_list(mes);
    printf("\n");
    print_h_matrix(gen_h);
    printf("\n");
    print_a_matrix(gen_a);
    printf("\n");


    h_list *res_h = encode_ldpc(gen_h, mes);
    h_list *res_a = product_a(gen_a, mes);
    print_h_list(res_h);
    printf("\n");
    print_h_list(res_a);


    free_h_matrix(test);
    free_h_matrix(gen_h);
    free_a_matrix(gen_a);
    free_h_list(mes);
    free_h_list(res_h);
    free_h_list(res_a);

/*
    bar *mes = barcreate(8920);

    bit_array_random(mes, 1.0);

    bar *gen = initGenerator();
    //bar *res = combine(gen, encode(mes));
    bar *res = encode(mes);

    uint n = barlen(res);
    double noisy[n];

    addNoise(res, s, noisy);

    bar *toto = decodeStreamParallel(noisy, s, 100);

    print(toto);

    printf("Number of error while decoding : %d\n", difference(mes, toto));


    bardestroy(toto);
    bardestroy(mes);
    bardestroy(res);
*/
    return 0;
}
