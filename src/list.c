#include <stdio.h>
#include <stdlib.h>

#include "list.h"


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
    if (list_i->n + l >= list_i->m_s)
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
            res->list_m[i + n] = cil(1, 1);
            res->list_m[i + n]->list[0] = i;

            res->list_n[i] = cil(tmp->list_n[i]->n, tmp->list_n[i]->n + 1);
            for (size_t k = 0; k < res->list_m[i]->n; k++)
            {
                res->list_n[i]->list[k] = tmp->list_n[i]->list[k];
            }
            append_i(res->list_n[i], n + i);
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
