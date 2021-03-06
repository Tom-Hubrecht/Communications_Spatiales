#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "demo.h"
#include "list.h"
#include "ldpc.h"
#include "turbocode.h"
#include "random.h"
#include "basic.h"


int demo_turbo_basic(char *file_name, double s)
{
    char debug = 0;

    char path[1024];        // The source file
    char path_n[1024];      // The noisy file
    char path_d[1024];      // The decoded file

    getcwd(path, 1024);
    strcat(path, "/../files/bytes/");
    strcat(path, file_name);

    strcpy(path_n, path);
    strcpy(path_d, path);

    strcat(path, ".bt");
    strcat(path_n, "_turbo_n.bt");
    strcat(path_d, "_turbo_d.bt");

    FILE *fp = fopen(path, "r");
    FILE *fn = fopen(path_n, "w");
    FILE *fd = fopen(path_d, "w");

    if (debug)
    {
        printf("Opened files\n");
    }

    hh_list *mes_d = read_bit_file_h(fp, 8920);
    hh_list *mes_n = deep_copy(mes_d);

    if (debug)
    {
        printf("Created hh_lists\n");
    }

    printf("Number of slices : %d\n", mes_d->n);

    for (size_t i = 0; i < mes_d->n; i++)
    {
        printf("\nSlice number %d :\n", i);

        h_list *res = encode_turbo(mes_d->list[i]);
        printf("\t- Encoded\n");

        s_list *noisy = add_noise(res, s);
        printf("\t- Added gaussian noise\n");

        h_list *noisy_h = decode_h_basic(noisy);
        if (debug)
        {
            printf("\t- Performed basic decoding\n");
        }

        copy_h(noisy_h, mes_n->list[i], 8920, 3);
        if (debug)
        {
            printf("\t- Copied basic decoding\n");
        }

        h_list *dec = decode_turbo_basic(noisy, s);
        printf("\t- Decoded\n");

        copy_h(dec, mes_d->list[i], 8920, 1);
        if (debug)
        {
            printf("\t- Copied decoding\n");
        }

        free_h_list(dec);
        free_h_list(noisy_h);
        free_s_list(noisy);
        free_h_list(res);
    }

    write_bit_file_h(fn, mes_n);
    write_bit_file_h(fd, mes_d);

    fclose(fp);
    fclose(fn);
    fclose(fd);

    free_hh_list(mes_d);
    free_hh_list(mes_n);

    return 0;
}


int demo_turbo_iter(char *file_name, double s, size_t i_max)
{
    char debug = 1;

    char path[1024];        // The source file
    char path_n[1024];      // The noisy file
    char path_d[1024];      // The decoded file

    getcwd(path, 1024);
    strcat(path, "/../files/bytes/");
    strcat(path, file_name);

    strcpy(path_n, path);
    strcpy(path_d, path);

    strcat(path, ".bt");
    strcat(path_n, "_turbo_n.bt");
    strcat(path_d, "_turbo_d.bt");

    FILE *fp = fopen(path, "r");
    FILE *fn = fopen(path_n, "w");
    FILE *fd = fopen(path_d, "w");

    if (debug)
    {
        printf("Opened files\n");
    }

    hh_list *mes_d = read_bit_file_h(fp, 8920);
    hh_list *mes_n = deep_copy(mes_d);

    if (debug)
    {
        printf("Created hh_lists\n");
    }

    printf("Number of slices : %d\n", mes_d->n);

    for (size_t i = 0; i < mes_d->n; i++)
    {
        printf("\nSlice number %d :\n", i);

        h_list *res = encode_turbo(mes_d->list[i]);
        printf("\t- Encoded\n");

        s_list *noisy = add_noise(res, s);
        printf("\t- Added gaussian noise with s = %f\n", s);

        h_list *noisy_h = decode_h_basic(noisy);
        if (debug)
        {
            printf("\t- Performed basic decoding\n");
        }

        copy_h(noisy_h, mes_n->list[i], 8920, 3);
        if (debug)
        {
            printf("\t- Copied basic decoding\n");
        }

        h_list *dec = decode_turbo_iter(noisy, s, i_max);
        printf("\t- Decoded\n");

        dec->n = 8920;
        printf("\t\tNumber of errors : %d\n", nb_errors(mes_d->list[i], dec));

        copy_h(dec, mes_d->list[i], 8920, 1);
        if (debug)
        {
            printf("\t- Copied decoding\n");
        }

        free_h_list(dec);
        free_h_list(noisy_h);
        free_s_list(noisy);
        free_h_list(res);
    }

    write_bit_file_h(fn, mes_n);
    write_bit_file_h(fd, mes_d);

    fclose(fp);
    fclose(fn);
    fclose(fd);

    free_hh_list(mes_d);
    free_hh_list(mes_n);

    return 0;
}


int demo_ldpc_basic(char *file_name, double s, size_t i_max)
{
    char debug = 1;
    int err;

    char path[1024];        // The source file
    char path_n[1024];      // The noisy file
    char path_d[1024];      // The decoded file

    getcwd(path, 1024);
    strcat(path, "/../files/bytes/");
    strcat(path, file_name);

    strcpy(path_n, path);
    strcpy(path_d, path);

    strcat(path, ".bt");
    strcat(path_n, "_ldpc_basic_n.bt");
    strcat(path_d, "_ldpc_basic_d.bt");

    FILE *fp = fopen(path, "r");
    FILE *fn = fopen(path_n, "w");
    FILE *fd = fopen(path_d, "w");
    if (debug)
    {
        printf("Opened files\n");
    }

    hh_list *mes_d = read_bit_file_h(fp, 8920);
    hh_list *mes_n = deep_copy(mes_d);
    if (debug)
    {
        printf("Created hh_lists\n");
    }

    h_matrix *base = create_base(8920, 40, 20);
    a_matrix *gen = cgm_a(base);
    a_matrix *dec = cdm_a(base);
    if (debug)
    {
        printf("Created LDPC matrices\n");
    }

    printf("Number of slices : %d\n", mes_d->n);

    for (size_t i = 0; i < mes_d->n; i++)
    {
        printf("\nSlice number %d :\n", i);

        h_list *res = encode_ldpc_a(gen, mes_d->list[i]);
        printf("\t- Encoded\n");

        s_list *noisy = add_noise(res, s);
        printf("\t- Added gaussian noise\n");

        h_list *noisy_h = decode_h_basic(noisy);
        if (debug)
        {
            printf("\t- Performed basic decoding\n");
        }

        copy_h(noisy_h, mes_n->list[i], 8920, 1);
        if (debug)
        {
            printf("\t- Copied basic decoding\n");
        }

        err = decode_ldpc_a_basic(dec, noisy_h, i_max);
        printf("\t- Decoded\n");
        if (debug)
        {
            printf("\t\t- Iterations : %d\n", err);
        }

        copy_h(noisy_h, mes_d->list[i], 8920, 1);
        if (debug)
        {
            printf("\t- Copied decoding\n");
        }

        free_h_list(noisy_h);
        free_s_list(noisy);
        free_h_list(res);
    }

    write_bit_file_h(fn, mes_n);
    write_bit_file_h(fd, mes_d);

    fclose(fp);
    fclose(fn);
    fclose(fd);

    free_hh_list(mes_d);
    free_hh_list(mes_n);

    free_h_matrix(base);
    free_a_matrix(gen);
    free_a_matrix(dec);

    return 0;
}


int demo_ldpc_proba(char *file_name, double s, size_t i_max)
{
    char debug = 1;
    int err;

    char path[1024];        // The source file
    char path_n[1024];      // The noisy file
    char path_d[1024];      // The decoded file

    getcwd(path, 1024);
    strcat(path, "/../files/bytes/");
    strcat(path, file_name);

    strcpy(path_n, path);
    strcpy(path_d, path);

    strcat(path, ".bt");
    strcat(path_n, "_ldpc_proba_n.bt");
    strcat(path_d, "_ldpc_proba_d.bt");

    FILE *fp = fopen(path, "r");
    FILE *fn = fopen(path_n, "w");
    FILE *fd = fopen(path_d, "w");
    if (debug)
    {
        printf("Opened files\n");
    }

    hh_list *mes_d = read_bit_file_h(fp, 8920);
    hh_list *mes_n = deep_copy(mes_d);
    if (debug)
    {
        printf("Created hh_lists\n");
    }

    h_matrix *base = create_base(8920, 40, 20);
    a_matrix *gen = cgm_a(base);
    a_matrix *dec = cdm_a(base);
    if (debug)
    {
        printf("Created LDPC matrices\n");
    }

    printf("Number of slices : %d\n", mes_d->n);

    for (size_t i = 0; i < mes_d->n; i++)
    {
        printf("\nSlice number %d :\n", i);

        h_list *res = encode_ldpc_a(gen, mes_d->list[i]);
        printf("\t- Encoded\n");

        s_list *noisy = add_noise(res, s);
        printf("\t- Added gaussian noise\n");

        h_list *noisy_h = decode_h_basic(noisy);
        if (debug)
        {
            printf("\t- Performed basic decoding\n");
        }

        copy_h(noisy_h, mes_n->list[i], 8920, 1);
        if (debug)
        {
            printf("\t- Copied basic decoding\n");
        }

        h_list *dem = decode_ldpc_proba(dec, noisy, s, i_max);
        printf("\t- Decoded\n");

        copy_h(dem, mes_d->list[i], 8920, 1);
        if (debug)
        {
            printf("\t- Copied decoding\n");
        }

        free_h_list(noisy_h);
        free_s_list(noisy);
        free_h_list(res);
        free_h_list(dem);
    }

    write_bit_file_h(fn, mes_n);
    write_bit_file_h(fd, mes_d);

    fclose(fp);
    fclose(fn);
    fclose(fd);

    free_hh_list(mes_d);
    free_hh_list(mes_n);

    free_h_matrix(base);
    free_a_matrix(gen);
    free_a_matrix(dec);

    return 0;
}


int demo_turbo_graph(size_t nb_iter, double p)
{
    char path[1024];

    getcwd(path, 1024);
    strcat(path, "/../files/graph/turbo.grp");

    FILE *fp = fopen(path, "w");

    double noise[] = {2, 1, 0.9, 0.8, 0.7, 0.6, 0.5,
                                0.4, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01};
    double s;
    size_t nb_err;

    double ber;

    h_list *mes = chl(8920, 8920);
    random_h(mes, p);

    h_list *enc = encode_turbo(mes);

    for (size_t i = 0; i < 14; i++)
    {
        nb_err = 0;
        s = pow(10, -0.1*noise[i]);

        for (size_t k = 0; k < nb_iter; k++)
        {
            s_list *sig = add_noise(enc, s);
            h_list *res = decode_turbo_iter(sig, s, 20);
            res->n = 8920;

            nb_err += nb_errors(res, mes);

            free_s_list(sig);
            free_h_list(res);
        }
        printf("%f -- %d / %d\n", s, nb_err, 8920*nb_iter);

        ber = ((double) nb_err) / (8920.0 * nb_iter);
        fprintf(fp, "%f %f\n", noise[i], ber);
    }


    fclose(fp);

    free_h_list(mes);
}


int demo_ldpc_graph(size_t nb_iter, double p)
{
    char path[1024];

    getcwd(path, 1024);
    strcat(path, "/../files/graph/ldpc.grp");

    FILE *fp = fopen(path, "w");

    double noise[] = {5, 4, 3, 2, 1, 0.9, 0.8, 0.7, 0.6, 0.5,
                                0.4, 0.3, 0.2, 0.1};
    double s;
    size_t nb_err;

    double ber;

    h_list *mes = chl(8920, 8920);
    random_h(mes, p);

    h_matrix *base = create_base(8920, 40, 20);
    a_matrix *gen = cgm_a(base);
    a_matrix *dec = cdm_a(base);

    h_list *enc = encode_ldpc_a(gen, mes);

    for (size_t i = 0; i < 14; i++)
    {
        nb_err = 0;
        s = pow(10, -0.1*noise[i]);

        for (size_t k = 0; k < nb_iter; k++)
        {
            s_list *sig = add_noise(enc, s);
            h_list *sig_h = decode_h_basic(sig);
            decode_ldpc_a_basic(dec, sig_h, 100);
            sig_h->n = 8920;

            nb_err += nb_errors(sig_h, mes);

            free_s_list(sig);
            free_h_list(sig_h);
        }
        printf("%f -- %d / %d\n", s, nb_err, 8920*nb_iter);

        ber = ((double) nb_err) / (8920.0 * nb_iter);
        fprintf(fp, "%f %f\n", noise[i], ber);
    }


    fclose(fp);

    free_h_list(mes);
    free_h_matrix(base);
    free_a_matrix(gen);
    free_a_matrix(dec);

}


int demo_base_graph(size_t nb_iter, double p)
{
    char path[1024];

    getcwd(path, 1024);
    strcat(path, "/../files/graph/base.grp");

    FILE *fp = fopen(path, "w");

    double noise[] = {5, 4, 3, 2, 1, 0.9, 0.8, 0.7, 0.6, 0.5,
                                0.4, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01};
    double s;
    size_t nb_err;

    double ber;

    h_list *mes = chl(8920, 8920);
    random_h(mes, p);

    for (size_t i = 0; i < 17; i++)
    {
        nb_err = 0;
        s = pow(10, -0.1*noise[i]);

        for (size_t k = 0; k < nb_iter; k++)
        {
            s_list *sig = add_noise(mes, s);
            h_list *res = decode_h_basic(sig);

            nb_err += nb_errors(res, mes);

            free_s_list(sig);
            free_h_list(res);
        }
        printf("%f -- %d / %d\n", s, nb_err, 8920*nb_iter);

        ber = ((double) nb_err) / (8920.0 * nb_iter);
        fprintf(fp, "%f %f\n", noise[i], ber);
    }


    fclose(fp);

    free_h_list(mes);

}
