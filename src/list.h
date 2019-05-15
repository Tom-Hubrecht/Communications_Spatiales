#include <stdio.h>

#ifndef LIST_H
#define LIST_H

typedef struct h_list {
  size_t n;      // Number of elements
  size_t m_s;    // Maximum number of elements
  char *list;
} h_list;

typedef struct s_list {
  size_t n;      // Number of elements
  size_t m_s;    // Maximum number of elements
  double *list;
} s_list;

typedef struct i_list {
  size_t n;      // Number of elements
  size_t m_s;    // Maximum number of elements
  int *list;
} i_list;

typedef struct h_matrix {
  size_t n;  // n lignes
  size_t m;  // m colonnes
  char *mat;
} h_matrix;

// Use a sparse matrix structure similar to that of David MacKay's
typedef struct a_matrix {
  size_t n;            // n lignes
  size_t m;            // m colonnes
  i_list **list_m;    // Liste des coordonées verticales non nulles
  i_list **list_n;    // Liste des coordonées horizontales non nulles

} a_matrix ;

typedef struct hh_list {
  size_t n;
  size_t m_e;         // Number of elements to write into file
  h_list **list;
} hh_list;

h_list * create_h_list(size_t n, size_t m_s);
i_list * create_i_list(size_t n, size_t m_s);
s_list * create_s_list(size_t n, size_t m_s);
h_matrix * create_h_matrix(size_t n, size_t m);
a_matrix * create_a_matrix(size_t n, size_t m);
hh_list * create_hh_list(size_t n, size_t m_e);

char get_h_list(h_list *list_h, size_t i);
double get_s_list(s_list *list_s, size_t i);
char get_h_matrix(h_matrix *mat_h, size_t i, size_t j);

void set_h_list(h_list *list_h, char x, size_t i);
void set_s_list(s_list *list_s, double x, size_t i);
void set_h_matrix(h_matrix *mat_h, char x, size_t i, size_t j);

void write_char_h(h_list *list_h, char x, size_t p);
char read_char_h(h_list *list_h, size_t p);

void write_bit_h(h_list *list_h, unsigned char x, size_t p);
char read_bit_h(h_list *list_h, size_t p);

void set_all_h_list(h_list *list_h, char x);
void set_all_i_list(i_list *list_i, int x);

char is_all_nil(h_list *list_h);

int append_i(i_list *list_i, int x);
int shift_i(i_list *list_i, int l);

int copy_h(h_list *list_a, h_list *list_b, size_t n, size_t s);
hh_list *deep_copy(hh_list *list_hh);

void max_i_list(i_list *list_i, i_list *res);

h_list * product_h(h_matrix *mat, h_list *vect);
h_list * product_a(a_matrix *mat, h_list *vect);
int product_a_in_place(a_matrix *mat, h_list *vect, h_list *res);
h_matrix * juxtapose_h(h_matrix *mat, char dir);
a_matrix * juxtapose_a(h_matrix *mat, char dir);

a_matrix * convert_h(h_matrix *mat);

hh_list * read_file_h(FILE *fp, size_t n);
int write_file_h(FILE *fp, hh_list *list_hh);

hh_list * read_bit_file_h(FILE *fp, size_t n);
int write_bit_file_h(FILE *fp, hh_list *list_hh);

void print_h_list(h_list *list_h);
void print_i_list(i_list *list_i);
void print_s_list(s_list *list_s);
void print_h_matrix(h_matrix *mat_h);
void print_a_matrix(a_matrix *mat_a);
void print_hh_list(hh_list *list_hh);

void free_h_list(h_list *list_h);
void free_i_list(i_list *list_i);
void free_s_list(s_list *list_s);
void free_h_matrix(h_matrix *mat_h);
void free_a_matrix(a_matrix *mat_a);
void free_hh_list(hh_list *list_hh);

// Abréviations
#define chl create_h_list
#define cil create_i_list
#define csl create_s_list
#define chm create_h_matrix
#define cam create_a_matrix
#define chhl create_hh_list

#define ghl get_h_list
#define gsl get_s_list
#define ghm get_h_matrix

#define shl set_h_list
#define ssl set_s_list
#define shm set_h_matrix

#define ph product_h

#endif
