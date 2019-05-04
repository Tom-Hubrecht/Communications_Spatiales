#include <stdio.h>

#ifndef LIST_H
#define LIST_H

typedef struct h_list {
  int n;      // Number of elements
  int m_s;    // Maximum number of elements
  char *list;
} h_list;

typedef struct s_list {
  int n;      // Number of elements
  int m_s;    // Maximum number of elements
  double *list;
} s_list;

typedef struct i_list {
  int n;      // Number of elements
  int m_s;    // Maximum number of elements
  int *list;
} i_list;

typedef struct h_matrix {
  int n;  // n lignes
  int m;  // m colonnes
  char *mat;
} h_matrix;

// Use a sparse matrix structure similar to that of David MacKay's
typedef struct a_matrix {
  int n;            // n lignes
  int m;            // m colonnes
  i_list **list_m;    // Liste des coordonées verticales non nulles
  i_list **list_n;    // Liste des coordonées horizontales non nulles

} a_matrix ;

h_list * create_h_list(size_t n, size_t m_s);
i_list * create_i_list(size_t n, size_t m_s);
s_list * create_s_list(size_t n, size_t m_s);
h_matrix * create_h_matrix(size_t n, size_t m);

char get_h_list(h_list *list_h, size_t i);
double get_s_list(s_list *list_s, size_t i);
char get_h_matrix(h_matrix *mat_h, size_t i, size_t j);

void set_h_list(h_list *list_h, char x, size_t i);
void set_s_list(s_list *list_s, double x, size_t i);
void set_h_matrix(h_matrix *mat_h, char x, size_t i, size_t j);

void set_all_h_list(h_list *list_h, char x);

int append_i(i_list *list_i, int x);

h_list * product_h(h_matrix *mat, h_list *vect);
h_list * product_a(a_matrix *mat, h_list *vect);
h_matrix * juxtapose_h(h_matrix *mat, char dir);
a_matrix * juxtapose_a(h_matrix *mat, char dir);

a_matrix * convert_h(h_matrix *mat);

void print_h_list(h_list *list_h);
void print_i_list(i_list *list_i);
void print_s_list(s_list *list_s);
void print_h_matrix(h_matrix *mat_h);
void print_a_matrix(a_matrix *mat_a);

void free_h_list(h_list *list_h);
void free_i_list(i_list *list_i);
void free_s_list(s_list *list_s);
void free_h_matrix(h_matrix *mat_h);
void free_a_matrix(a_matrix *mat_a);

// Abréviations
#define chl create_h_list
#define cil create_i_list
#define csl create_s_list
#define chm create_h_matrix
#define cam create_a_matrix

#define ghl get_h_list
#define gsl get_s_list
#define ghm get_h_matrix

#define shl set_h_list
#define ssl set_s_list
#define shm set_h_matrix

#define ph product_h

#endif
