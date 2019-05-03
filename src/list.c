#include <stdio.h>
#include <stdlib.h>

#include "list.h"


h_list * create_h_list(size_t n)
{
    h_list *res = malloc(sizeof(h_list));
    res->n = n;
    res->list = calloc(n, sizeof(char));

    return res;
}


i_list * create_i_list(size_t n)
{
    i_list *res = malloc(sizeof(i_list));
    res->n = n;
    res->list = calloc(n, sizeof(int));

    return res;
}


s_list * create_s_list(size_t n)
{
    s_list *res = malloc(sizeof(s_list));
    res->n = n;
    res->list = calloc(n, sizeof(double));

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


void set_all_h_list(h_list *list_h, char x)
{
    for (size_t i = 0; i < list_h->n; i++)
    {
        list_h->list[i] = x;
    }
}


h_list * product_h(h_matrix *mat, h_list *vect)
{
    if (mat->m != vect->n)
    {
        return NULL;
    }

    h_list *res = chl(mat->n);
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


// If dir = 0 -> [A | I]
// If dir = 1 -> [I
//                A]
h_matrix * juxtapose(h_matrix *mat, char dir)
{
    int n = mat->n;
    int m = mat-> m;
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