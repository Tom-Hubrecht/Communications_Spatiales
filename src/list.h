#include <stdio.h>

#ifndef LIST_H
#define LIST_H

typedef struct {
  int n;
  char *list;
} h_list;

typedef struct {
  int n;
  double *list;
} s_list;

typedef struct {
  int n;
  int *list;
} i_list;

typedef struct {
  int n;  // n lignes
  int m;  // m colonnes
  char *mat;
} h_matrix;

h_list * create_h_list(size_t n);
i_list * create_i_list(size_t n);
s_list * create_s_list(size_t n);
h_matrix * create_h_matrix(size_t n, size_t m);

char get_h_list(h_list *list_h, size_t i);
double get_s_list(s_list *list_s, size_t i);
char get_h_matrix(h_matrix *mat_h, size_t i, size_t j);

void set_h_list(h_list *list_h, char x, size_t i);
void set_s_list(s_list *list_s, double x, size_t i);
void set_h_matrix(h_matrix *mat_h, char x, size_t i, size_t j);

void set_all_h_list(h_list *list_h, char x);

h_list * product_h(h_matrix *mat, h_list *vect);
h_matrix * juxtapose(h_matrix *mat, char dir);


void print_h_list(h_list *list_h);
void print_i_list(i_list *list_i);
void print_s_list(s_list *list_s);
void print_h_matrix(h_matrix *mat_h);

void free_h_list(h_list *list_h);
void free_i_list(i_list *list_i);
void free_s_list(s_list *list_s);
void free_h_matrix(h_matrix *mat_h);

// Abréviations
#define chl create_h_list
#define chi create_i_list
#define csl create_s_list
#define chm create_h_matrix

#define ghl get_h_list
#define gsl get_s_list
#define ghm get_h_matrix

#define shl set_h_list
#define ssl set_s_list
#define shm set_h_matrix

#define ph product_h

#endif
