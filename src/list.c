#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "list.h"
#include "basic.h"


h_list * create_h_list(size_t n, size_t m_s)
{
    h_list *res = malloc(sizeof(h_list));
    res->n = n;
    res->m_s = m_s;
    res->list = calloc(m_s, sizeof(char));

    return res;
}


i_list * create_i_list(size_t n, size_t m_s)
{
    i_list *res = malloc(sizeof(i_list));
    res->n = n;
    res->m_s = m_s;
    res->list = calloc(m_s, sizeof(int));

    return res;
}


s_list * create_s_list(size_t n, size_t m_s)
{
    s_list *res = malloc(sizeof(s_list));
    res->n = n;
    res->m_s = m_s;
    res->list = calloc(m_s, sizeof(double));

    return res;
}


h_matrix * create_h_matrix(size_t n, size_t m)
{
    h_matrix *res = malloc(sizeof(h_matrix));
    res->n = n;
    res->m = m;
    res->mat = calloc(n * m, sizeof(char));

    return res;
}


a_matrix * create_a_matrix(size_t n, size_t m)
{
    a_matrix *res = malloc(sizeof(a_matrix));
    res-> n = n;
    res->m = m;
    res->list_m = malloc(m * sizeof(i_list *));
    res->list_n = malloc(n * sizeof(i_list *));

    return res;
}


hh_list * create_hh_list(size_t n, size_t m_e)
{
    hh_list *res = malloc(sizeof(hh_list));
    res->n = n;
    res->m_e = m_e;
    res->list = malloc(n * sizeof(h_list *));

    return res;
}


char get_h_list(h_list *list_h, size_t i)
{
    return list_h->list[i];
}


double get_s_list(s_list *list_s, size_t i)
{
    return list_s->list[i];
}


char get_h_matrix(h_matrix *mat_h, size_t i, size_t j)
{
    return mat_h->mat[(i * (mat_h->m)) + j];
}


void set_h_list(h_list *list_h, char x, size_t i)
{
    list_h->list[i] = x;
}


void set_s_list(s_list *list_s, double x, size_t i)
{
    list_s->list[i] = x;
}


void set_h_matrix(h_matrix *mat_h, char x, size_t i, size_t j)
{
    mat_h->mat[i * (mat_h->m) + j] = x;
}


void write_char_h(h_list *list_h, char x, size_t p)
{
    unsigned char t;
    if (p < (list_h->n - 7))
    {
        t = x;
        for (size_t i = 0; i < 8; i++)
        {
            list_h->list[p + i] = t % 2;
            t = t / 2;
        }
    }
}


char read_char_h(h_list *list_h, size_t p)
{
    if (p < (list_h->n - 7))
    {
        unsigned char x = 0;
        for (size_t i = 8; i > 0; i--)
        {
            x *= 2;
            x += list_h->list[p + i - 1];
        }
        return x;
    }
}


void write_bit_h(h_list *list_h, unsigned char x, size_t p)
{
    if (p < list_h->n)
    {
        list_h->list[p] = x % 2;
    }
}


char read_bit_h(h_list *list_h, size_t p)
{
    if (p < list_h->n)
    {
        return 255 * list_h->list[p];
    }
}


void set_all_h_list(h_list *list_h, char x)
{
    for (size_t i = 0; i < list_h->n; i++)
    {
        list_h->list[i] = x;
    }
}


void set_all_i_list(i_list *list_i, int x)
{
    for (size_t i = 0; i < list_i->n; i++)
    {
        list_i->list[i] = x;
    }
}


void set_all_s_list(s_list *list_s, double x)
{
    for (size_t i = 0; i < list_s->n; i++)
    {
        list_s->list[i] = x;
    }
}


// Determines if the h_list is only 0
char is_all_nil(h_list *list_h)
{
    char ok = 1;
    for (size_t i = 0; i < list_h->n; i++)
    {
        if (list_h->list[i])
        {
            ok = 0;
        }
    }
    return ok;
}


double min_s(s_list *list_s)
{
    double s = INFINITY;

    for (size_t i = 0; i < list_s->n; i++)
    {
        if (fabs(list_s->list[i]) < s)
        {
            s = fabs(list_s->list[i]);
        }
    }

    return s;
}


int append_i(i_list *list_i, int x)
{
    if (list_i->n == list_i->m_s)
    {
        // The list is maxed out
        return 1;
    }

    list_i->list[list_i->n] = x;
    list_i->n ++;

    return 0;
}


// Shift the list l elements to the right
int shift_i(i_list *list_i, int l)
{
    if (0 >= list_i->n + l >= list_i->m_s)
    {
        return 1;
    }

    list_i->n += l;
    for (size_t i = (list_i->n - 1); i >=0; i--)
    {
        list_i->list[i + list_i->n] = list_i->list[i];
    }
    return 0;
}


int substract_s(s_list *list_a, s_list *list_b)
{
    if (list_a->n != list_b->n)
    {
        return 1;
    }

    for (size_t i = 0; i < list_a->n; i++)
    {
        list_a->list[i] -= list_b->list[i];
    }
    return 0;
}


// Copy n elements of list_a into list_b if possible
int copy_h(h_list *list_a, h_list *list_b, size_t n, size_t s)
{
    if (n > list_b->m_s || n * s >list_a->n)
    {
        return 1;
    }

    list_b->n = n;
    for (size_t i = 0; i < n; i++)
    {
        list_b->list[i] = list_a->list[i * s];
    }
    return 0;
}


hh_list *deep_copy(hh_list *list_hh)
{
    hh_list *res = chhl(list_hh->n, list_hh->m_e);
    for (size_t i = 0; i < res->n; i++)
    {
        res->list[i] = chl(list_hh->list[i]->n, list_hh->list[i]->m_s);
        copy_h(list_hh->list[i], res->list[i], res->list[i]->n, 1);
    }
    return res;
}


void max_i_list(i_list *list_i, i_list *res)
{
    int m = list_i->list[0];
    res->n = 0;

    for (size_t i = 0; i < list_i->n; i++)
    {
        if (list_i->list[i] == m)
        {
            append_i(res, i);
        }

        if (list_i->list[i] > m)
        {
            m = list_i->list[i];
            res->n = 0;
            append_i(res, i);
        }
    }
}


h_list * product_h(h_matrix *mat, h_list *vect)
{
    if (mat->m != vect->n)
    {
        return NULL;
    }

    h_list *res = chl(mat->n, mat->n);
    for (size_t k = 0; k < mat->m; k++)
    {
        for (size_t i = 0; i < mat->n; i++)
        {
            if (ghm(mat, i, k) && ghl(vect, k))
            {
                res->list[i] = 1 - res->list[k];
            }
        }
    }
    return res;
}


h_list * product_a(a_matrix *mat, h_list *vect)
{
    if (mat->m != vect->n)
    {
        return NULL;
    }

    h_list *res = chl(mat->n, mat->n);
    for (size_t i = 0; i < mat->n; i++)
    {
        for (size_t k = 0; k < mat->list_n[i]->n; k++)
        {
            if (vect->list[mat->list_n[i]->list[k]])
            {
                res->list[i] ^= 1;
            }
        }
    }
    return res;
}


int product_a_in_place(a_matrix *mat, h_list *vect, h_list *res)
{
    if (mat->m != vect->n || mat->n > res->m_s)
    {
        return 1;
    }

    res->n = mat->n;
    set_all_h_list(res, 0);

    for (size_t i = 0; i < mat->n; i++)
    {
        for (size_t k = 0; k < mat->list_n[i]->n; k++)
        {
            if (vect->list[mat->list_n[i]->list[k]])
            {
                res->list[i] ^= 1;
            }
        }
    }
    return 0;
}


int product_s_ip(a_matrix *mat, s_list *vect, s_list *res)
{
    if (mat->m != vect->n || mat->n > res->m_s)
    {
        return 1;
    }

    res->n = mat->n;
    set_all_s_list(res, 0.0);

    for (size_t i = 0; i < mat->n; i++)
    {
        for (size_t k = 0; k < mat->list_n[i]->n; k++)
        {
            res->list[i] += vect->list[mat->list_n[i]->list[k]];
        }
    }
    return 0;
}


// If dir = 0 -> [A | I]
// If dir = 1 -> [I
//                A]
h_matrix * juxtapose_h(h_matrix *mat, char dir)
{
    int n = mat->n;
    int m = mat->m;
    h_matrix *res;

    if (dir)
    {
        res = chm(n + m, m);
        // Copy mat into res and add identity
        for (size_t j = 0; j < mat->m; j++)
        {
            shm(res, 1, j, j);
            for (size_t i = 0; i < mat->n; i++)
            {
                shm(res, ghm(mat, i, j), i + m, j);
            }
        }
    }
    else
    {
        res = chm(n, m + n);
        // Copy mat into res and add identity
        for (size_t i = 0; i < mat->n; i++)
        {
            shm(res, 1, i, i + m);
            for (size_t j = 0; j < mat->m; j++)
            {
                shm(res, ghm(mat, i, j), i, j);
            }
        }
    }
    return res;
}


// If dir = 0 -> [A | I]
// If dir = 1 -> [I
//                A]
a_matrix * juxtapose_a(h_matrix *mat, char dir)
{
    int n = mat->n;
    int m = mat->m;
    a_matrix *tmp = convert_h(mat);
    a_matrix *res;

    if (dir)
    {
        res = cam(n + m, m);

        // Copy tmp into res and add identity
        for (size_t i = 0; i < mat->m; i++)
        {
            res->list_n[i] = cil(1, 1);
            res->list_n[i]->list[0] = i;

            res->list_m[i] = cil(tmp->list_m[i]->n + 1, tmp->list_m[i]->n + 1);
            res->list_m[i]->list[0] = i;
            for (size_t k = 1; k < res->list_m[i]->n; k++)
            {
                res->list_m[i]->list[k] = tmp->list_m[i]->list[k -1] + m;
            }
        }

        for (size_t i = 0; i < mat->n; i++)
        {
            res->list_n[i + m] = cil(tmp->list_n[i]->n, tmp->list_n[i]->n);
            for (size_t k = 0; k < res->list_n[i + m]->n; k++)
            {
                res->list_n[i + m]->list[k] = tmp->list_n[i]->list[k];
            }
        }
    }
    else
    {
        res = cam(n, m + n);

        // Copy tmp into res and add identity
        for (size_t i = 0; i < mat->n; i++)
        {
            res->list_m[i + m] = cil(1, 1);
            res->list_m[i + m]->list[0] = i;

            res->list_n[i] = cil(tmp->list_n[i]->n, tmp->list_n[i]->n + 1);
            for (size_t k = 0; k < res->list_n[i]->n; k++)
            {
                res->list_n[i]->list[k] = tmp->list_n[i]->list[k];
            }
            append_i(res->list_n[i], m + i);
        }

        for (size_t i = 0; i < mat->m; i++)
        {
            res->list_m[i] = cil(tmp->list_m[i]->n, tmp->list_m[i]->n);
            for (size_t k = 0; k < res->list_m[i]->n; k++)
            {
                res->list_m[i]->list[k] = tmp->list_m[i]->list[k];
            }
        }
    }
    free_a_matrix(tmp);
    return res;
}


a_matrix * convert_h(h_matrix *mat)
{
    a_matrix *res = cam(mat->n, mat->m);

    // Initialise list_m and list_n
    for (size_t i = 0; i < mat->m; i++)
    {
        res->list_m[i] = cil(0, mat->n);
    }

    for (size_t i = 0; i < mat->n; i++)
    {
        res->list_n[i] = cil(0, mat->m);
    }

    // Fill the matrix
    for (size_t i = 0; i < mat->n; i++)
    {
        for (size_t j = 0; j < mat->m; j++)
        {
            if (ghm(mat, i, j))
            {
                append_i(res->list_m[j], i);
                append_i(res->list_n[i], j);
            }
        }
    }
    return res;
}


hh_list * read_file_h(FILE *fp, size_t n)
{
    if (n % 8)
    {
        return NULL;
    }

    size_t s = file_size(fp);
    size_t m = (s % n) ? (s / n) + 1 : (s / n); // Size in bytes
    m *= 8;                                     // Size in bits

    hh_list *res = chhl(m, s);
    for (size_t i = 0; i < m; i++)
    {
        res->list[i] = chl(n, n);
    }

    unsigned char *buffer = malloc(s * sizeof(char));
    fread(buffer, sizeof(char), s, fp);

    for (size_t i = 0; i < s; i++)
    {
        write_char_h(res->list[8 * i / n], buffer[i], (8 * i) % n);
    }

    return res;
}


int write_file_h(FILE *fp, hh_list *list_hh)
{
    size_t c = 0;
    unsigned char x;

    for (size_t i = 0; i < list_hh->n; i++)
    {
        if (list_hh->list[i]->n % 8)
        {
            return 1;
        }

        for (size_t j = 0; j < (list_hh->list[i]->n / 8); j++)
        {
            if (c < list_hh->m_e)
            {
                x = read_char_h(list_hh->list[i], 8 * j);
                fwrite(&x, sizeof(char), 1, fp);
            }
            c ++;
        }
    }
    return 0;
}


hh_list * read_bit_file_h(FILE *fp, size_t n)
{

    size_t s = file_size(fp);
    size_t m = (s % n) ? (s / n) + 1 : (s / n); // Number of lists to create

    hh_list *res = chhl(m, s);
    for (size_t i = 0; i < m; i++)
    {
        res->list[i] = chl(n, n);
    }

    unsigned char *buffer = malloc(s * sizeof(char));
    fread(buffer, sizeof(char), s, fp);

    for (size_t i = 0; i < s; i++)
    {
        write_bit_h(res->list[i / n], buffer[i], i % n);
    }

    return res;
}


int write_bit_file_h(FILE *fp, hh_list *list_hh)
{
    size_t c = 0;
    unsigned char x;

    for (size_t i = 0; i < list_hh->n; i++)
    {
        for (size_t j = 0; j < (list_hh->list[i]->n); j++)
        {
            if (c < list_hh->m_e)
            {
                x = read_bit_h(list_hh->list[i], j);
                fwrite(&x, sizeof(char), 1, fp);
            }
            c ++;
        }
    }
    return 0;
}


void print_h_list(h_list *list_h)
{
    printf("[");
    for (size_t i = 0; (i + 1) < list_h->n; i++)
    {
        printf("%d, ", list_h->list[i]);
    }
    printf("%d]\n", list_h->list[list_h->n - 1]);
}


void print_i_list(i_list *list_i)
{
    int n = list_i->n;
    printf("[");
    for (size_t i = 0; i < (n - 1); i++)
    {
        printf("%d, ", (list_i->list[i]));
    }
    printf("%d]\n", list_i->list[list_i->n - 1]);
}


void print_s_list(s_list *list_s)
{
    printf("[");
    for (size_t i = 0; (i + 1) < list_s->n; i++)
    {
        printf("%f, ", list_s->list[i]);
    }
    printf("%f]\n", list_s->list[list_s->n - 1]);
}


void print_h_matrix(h_matrix *mat_h)
{
    int n = mat_h->n;
    int m = mat_h->m;

    printf("[");
    for (size_t i = 0; (i + 1) < n; i++)
    {
        printf("[");
        for (size_t j = 0; (j + 1) < m; j++)
        {
            printf("%d, ", mat_h->mat[(i * m) + j]);
        }
        printf("%d],\n", mat_h->mat[(i * m) + m - 1]);
    }

    printf("[");
    for (size_t j = 0; (j + 1) < m; j++)
    {
        printf("%d, ", mat_h->mat[((n - 1) * m) + j]);
    }
    printf("%d]]\n", mat_h->mat[(n * m) - 1]);
}


void print_a_matrix(a_matrix *mat_a)
{
    // Print the m vertical lists
    for (size_t i = 0; i < mat_a->m; i++)
    {
        print_i_list(mat_a->list_m[i]);
    }
    printf("\n");

    // Print the n horizontal lists
    for (size_t i = 0; i < mat_a->n; i++)
    {
        print_i_list(mat_a->list_n[i]);
    }
}


void print_hh_list(hh_list *list_hh)
{
    for (size_t i = 0; i < list_hh->n; i++)
    {
        print_h_list(list_hh->list[i]);
    }
}


void free_h_list(h_list *list_h)
{
    free(list_h->list);
    free(list_h);
}


void free_i_list(i_list *list_i)
{
    free(list_i->list);
    free(list_i);
}


void free_s_list(s_list *list_s)
{
    free(list_s->list);
    free(list_s);
}


void free_h_matrix(h_matrix *mat_h)
{
    free(mat_h->mat);
    free(mat_h);
}


void free_a_matrix(a_matrix *mat_a)
{
    for (size_t i = 0; i < mat_a->m; i++)
    {
        free_i_list(mat_a->list_m[i]);
    }
    free(mat_a->list_m);

    for (size_t i = 0; i < mat_a->n; i++)
    {
        free_i_list(mat_a->list_n[i]);
    }
    free(mat_a->list_n);

    free(mat_a);
}


void free_hh_list(hh_list *list_hh)
{
    for (size_t i = 0; i < list_hh->n; i++)
    {
        free_h_list(list_hh->list[i]);
    }
    free(list_hh->list);

    free(list_hh);
}
