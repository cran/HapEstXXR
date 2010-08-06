  /* matalloc.c */
  /**************/

#include <stdlib.h>
#include "longhap.h"
#include "message.h"
#include "matalloc.h"
#include "string.h"


/*****************************************************************************/
/*                                                                           */
/* Multidimensional arrays for type char                                     */
/*                                                                           */
/*****************************************************************************/


/*** array of char ***********************************************************/


char *c_array(unsigned n)
{
  char *a = (char*) malloc(n * sizeof(char));

  if (a == NULL)
    fatal("Memory block of %d bytes unavailable", n * sizeof(char));
  return a;
}


char *init_c_array(char *a, unsigned n, char val)
{
  unsigned i;

  if (a == NULL)
    a = c_array(n);
  for (i = 0; i < n; i++)
    a[i] = val;
  return a;
}


char *assign_c_array(char *a, unsigned n, ui2d_mapfun f)
{
  unsigned i;

  if (a && f)
    for (i = 0; i < n; i++)
      a[i] = f(i);
  return a;
}


void destroy_c_array(char *a)
{
  if (a)
    free(a);
}


/*** 2-dimensional array of char *********************************************/


char **c_array2(unsigned x, unsigned y)
{
  char   **a;
  unsigned   i;

  if ((a = (char**) malloc((x + 1) * sizeof(char*))) == NULL)
    fatal("Memory block of %d bytes unavailable", ((x + 1) * sizeof(char*)));
  if ((a[0] = (char*) malloc(x * y *  sizeof(char))) == NULL)
    fatal("Memory block of %d bytes unavailable", x * y *  sizeof(char));
  a[x] = NULL;
  for (i = 1; i < x; i++)
    a[i] = &a[0][i * y];
  return a;
}


char **init_c_array2(char **a, unsigned x, unsigned y, char val)
{
  unsigned   i, j;

  if (a == NULL)
    a = c_array2(x, y);
  for (i = 0; i < x; i++)
    for (j = 0; j < y; j++)
      a[i][j] = val;
  return a;
}


void destroy_c_array2(char **a)
{
  if (a) {
    free(a[0]);
    free(a);
  }
}


/*** 3-dimensional array of char *********************************************/


char ***c_array3(unsigned x, unsigned y, unsigned z)
{
  char   ***a;
  unsigned    i;

  if ((a = (char***) malloc((x + 1) * sizeof(char**))) == NULL)
    fatal("Memory block of %d bytes unavailable", (x + 1) * sizeof(char**));
  if ((a[0] = (char**) malloc((x + 1) * y * sizeof(char*))) == NULL)
    fatal("Memory block of %d bytes unavailable", (x + 1) * y * sizeof(char*));
  if ((a[0][0] = (char*) malloc((x * y * z) * sizeof(char))) == NULL)
    fatal("Memory block of %d bytes unavailable", (x * y * z) * sizeof(char));
  a[x] = NULL;
  for (i = 1; i < x; i++)
    a[i] = &a[0][y * i];
  for (i = 1; i < (x * y); i++)
    a[0][i] = &a[0][0][z * i];
#ifdef DebugInfo
  printf("Consuming %d bytes\n",   ((x + 1)       * sizeof(char**))
                                 + (((x + 1) * y) * sizeof(char*))
                                 + ((x * y * z)   * sizeof(char)));
#endif
  return a;
}


char ***init_c_array3(char ***a, unsigned x, unsigned y, unsigned z, char val)
{
  unsigned i, j, k;

  if (a == NULL)
    a = c_array3(x, y, z);
  for (i = 0; i < x; i++)
    for (j = 0; j < y; j++)
      for (k = 0; k < z; k++)
        a[i][j][k] = val;
  return a;
}


void destroy_c_array3(char ***a)
{
  if (a) {
    free(a[0][0]);
    free(a[0]);
    free(a);
  }
}


/*****************************************************************************/
/*                                                                           */
/* Multidimensional arrays for type double                                   */
/*                                                                           */
/*****************************************************************************/


/*** array of double *********************************************************/

#define N          1000
#define Efficiency

double *d_array(unsigned n)
{
  double *a = (double*) malloc(n * sizeof(double));

  if (a == NULL)
    fatal("Memory block of %d bytes unavailable", n * sizeof(double));
  return a;
}


double *init_d_array(double *a, unsigned n, double val)
{
#ifdef Efficiency
  static double *zero = NULL;
#endif

  int i;

#ifdef Efficiency
  if (zero == NULL) {
    zero = (double*) malloc(N * sizeof(double));
    for (i = 0; i < N; i++)
      zero[i] = 0.0;
  }
#endif
  if (a == NULL)
    a = d_array(n);
#ifdef Efficiency
  if (val == 0.0) {
    for (i = n; i > 0; i -= N)
      memcpy(&a[n - i], zero, sizeof(double) * ((i > N) ? N : i));
  } else {
#endif
    for (i = 0; i < n; i++)
      a[i] = val;
#ifdef Efficiency
  }
#endif
  return a;
}


double *dup_d_array(double *a, double *b, unsigned n)
{
  if (b) {
    if (a == NULL)
      a = d_array(n);
    memcpy(a, b, n * sizeof(double));
  }
  return a;
}


double *assign_d_array(double *a, unsigned n, ui2d_mapfun f)
{
  unsigned i;

  if (a && f)
    for (i = 0; i < n; i++)
      a[i] = f(i);
  return a;
}


void destroy_d_array(double *a)
{
  if (a)
    free(a);
}


/*** array of double *********************************************************/


double **d_array2(unsigned x, unsigned y)
{
  double   **a;
  unsigned   i;

  if ((a = (double**) malloc((x + 1) * sizeof(double*))) == NULL)
    fatal("Memory block of %d bytes unavailable", (x + 1) * sizeof(double*));
  if ((a[0] = (double*) malloc(x * y *  sizeof(double))) == NULL)
    fatal("Memory block of %d bytes unavailable", x * y *  sizeof(double));
  a[x] = NULL;
  for (i = 1; i < x; i++)
    a[i] = &a[0][i * y];
  return a;
}


double **init_d_array2(double **a, unsigned x, unsigned y, double val)
{
  unsigned   i, j;

  if (a == NULL)
    a = d_array2(x, y);
  for (i = 0; i < x; i++)
    for (j = 0; j < y; j++)
      a[i][j] = val;
  return a;
}


void destroy_d_array2(double **a)
{
  if (a) {
    free(a[0]);
    free(a);
  }
}


/*** 3-dimensional array of double *******************************************/


double ***d_array3(unsigned x, unsigned y, unsigned z)
{
  double   ***a;
  unsigned    i;

  if ((a = (double***) malloc((x + 1) * sizeof(double**))) == NULL)
    fatal("Memory block of %d bytes unavailable", (x + 1) * sizeof(double**));
  if ((a[0] = (double**) malloc(((x + 1) * y) * sizeof(double*))) == NULL)
    fatal("Memory block of %d bytes unavailable", (x + 1) * y * sizeof(double*));
  if ((a[0][0] = (double*) malloc(x * y * z * sizeof(double))) == NULL)
    fatal("Memory block of %d bytes unavailable", x * y * z * sizeof(double));
  a[x] = NULL;
  for (i = 1; i < x; i++)
    a[i] = &a[0][i * y];
  for (i = 1; i < (x * y); i++)
    a[0][i] = &a[0][0][z * i];
#ifdef DebugInfo
  printf("Consuming %d bytes\n",   ((x + 1)       * sizeof(double**))
                                 + (((x + 1) * y) * sizeof(double*))
                                 + ((x * y * z)   * sizeof(double)));
#endif
  return a;
}


double ***init_d_array3(double ***a, unsigned x, unsigned y, unsigned z, double val)
{
  unsigned i, j, k;

  if (a == NULL)
    a = d_array3(x, y, z);
  for (i = 0; i < x; i++)
    for (j = 0; j < y; j++)
      for (k = 0; k < z; k++)
        a[i][j][k] = val;
  return a;
}


double ***assign_d_array3(double ***a, unsigned x, unsigned y, unsigned z)
{
  unsigned i, j, k;

  if (a == NULL)
    a = d_array3(x, y, z);
  for (i = 0; i < x; i++)
    for (j = 0; j < y; j++)
      for (k = 0; k < z; k++)
        a[i][j][k] = (double) (i + j + k);
  return a;
}


void destroy_d_array3(double ***a)
{
  if (a) {
    free(a[0][0]);
    free(a[0]);
    free(a);
  }
}


/*****************************************************************************/
/*                                                                           */
/* Multidimensional arrays for type int                                      */
/*                                                                           */
/*****************************************************************************/


/*** array of int ************************************************************/


int *i_array(unsigned n)
{
  int *a = (int*) malloc(n * sizeof(int));

  if (a == NULL)
    fatal("Memory block of %d bytes unavailable", n * sizeof(int));
  return a;
}


int *init_i_array(int *a, unsigned n, int val)
{
  unsigned i;

  if (a == NULL)
    a = i_array(n);
  for (i = 0; i < n; i++)
    a[i] = val;
  return a;
}


int *assign_i_array(int *a, unsigned n, ui2d_mapfun f)
{
  unsigned i;

  if (a && f)
    for (i = 0; i < n; i++)
      a[i] = f(i);
  return a;
}


void destroy_i_array(int *a)
{
  if (a)
    free(a);
}


/*** array of int ************************************************************/


int **i_array2(unsigned x, unsigned y)
{
  int   **a;
  unsigned   i;

  if ((a = (int**) malloc((x + 1) * sizeof(int*))) == NULL)
    fatal("Memory block of %d bytes unavailable", (x + 1) * sizeof(int*));
  if ((a[0] = (int*) malloc(x * y *  sizeof(int))) == NULL)
    fatal("Memory block of %d bytes unavailable", x * y *  sizeof(int));
  a[x] = NULL;
  for (i = 1; i < x; i++)
    a[i] = &a[0][i * y];
  return a;
}


int **init_i_array2(int **a, unsigned x, unsigned y, int val)
{
  unsigned   i, j;

  if (a == NULL)
    a = i_array2(x, y);
  for (i = 0; i < x; i++)
    for (j = 0; j < y; j++)
      a[i][j] = val;
  return a;
}


void destroy_i_array2(int **a)
{
  if (a) {
    free(a[0]);
    free(a);
  }
}


/*** 3-dimensional array of int **********************************************/


int ***i_array3(unsigned x, unsigned y, unsigned z)
{
  int   ***a;
  unsigned    i;

  if ((a = (int***) malloc((x + 1) * sizeof(int**))) == NULL)
    fatal("Memory block of %d bytes unavailable", (x + 1) * sizeof(int**));
  if ((a[0] = (int**) malloc(((x + 1) * y) * sizeof(int*))) == NULL)
    fatal("Memory block of %d bytes unavailable", (x + 1) * y * sizeof(int*));
  if ((a[0][0] = (int*) malloc((x * y * z) * sizeof(int))) == NULL)
    fatal("Memory block of %d bytes unavailable", x * y * z * sizeof(int));
  a[x] = NULL;
  for (i = 1; i < x; i++)
    a[i] = &a[0][y * i];
  for (i = 1; i < (x * y); i++)
    a[0][i] = &a[0][0][z * i];
#ifdef DebugInfo
  printf("Consuming %d bytes\n",   ((x + 1)       * sizeof(int**))
                                 + (((x + 1) * y) * sizeof(int*))
                                 + ((x * y * z)   * sizeof(int)));
#endif
  return a;
}


int ***init_i_array3(int ***a, unsigned x, unsigned y, unsigned z, int val)
{
  unsigned i, j, k;

  if (a == NULL)
    a = i_array3(x, y, z);
  for (i = 0; i < x; i++)
    for (j = 0; j < y; j++)
      for (k = 0; k < z; k++)
        a[i][j][k] = val;
  return a;
}


void destroy_i_array3(int ***a)
{
  if (a) {
    free(a[0][0]);
    free(a[0]);
    free(a);
  }
}


/*****************************************************************************/
/*                                                                           */
/* Multidimensional arrays for type unsigned                                 */
/*                                                                           */
/*****************************************************************************/


/*** array of unsigned *******************************************************/


unsigned *u_array(unsigned n)
{
  unsigned *a = (unsigned*) malloc(n * sizeof(unsigned));

  if (a == NULL)
    fatal("Memory block of %d bytes unavailable", n * sizeof(unsigned));
  return a;
}


unsigned *init_u_array(unsigned *a, unsigned n, unsigned val)
{
  unsigned i;

  if (a == NULL)
    a = u_array(n);
  for (i = 0; i < n; i++)
    a[i] = val;
  return a;
}


void destroy_u_array(unsigned *a)
{
  if (a)
    free(a);
}


/*** array of unsigned *******************************************************/


unsigned **u_array2(unsigned x, unsigned y)
{
  unsigned   **a;
  unsigned   i;

  if ((a = (unsigned**) malloc((x + 1) * sizeof(unsigned*))) == NULL)
    fatal("Memory block of %d bytes unavailable", (x + 1) * sizeof(unsigned*));
  if ((a[0] = (unsigned*) malloc(x * y *  sizeof(unsigned))) == NULL)
    fatal("Memory block of %d bytes unavailable", x * y *  sizeof(unsigned));
  a[x] = NULL;
  for (i = 1; i < x; i++)
    a[i] = &a[0][i * y];
  return a;
}


unsigned **init_u_array2(unsigned **a, unsigned x, unsigned y, unsigned val)
{
  unsigned   i, j;

  if (a == NULL)
    a = u_array2(x, y);
  for (i = 0; i < x; i++)
    for (j = 0; j < y; j++)
      a[i][j] = val;
  return a;
}


void destroy_u_array2(unsigned **a)
{
  if (a) {
    free(a[0]);
    free(a);
  }
}


/*** 3-dimensional array of unsigned *****************************************/


unsigned ***u_array3(unsigned x, unsigned y, unsigned z)
{
  unsigned   ***a;
  unsigned    i;

  if ((a = (unsigned***) malloc((x + 1) * sizeof(unsigned**))) == NULL)
    fatal("Memory block of %d bytes unavailable", (x + 1) * sizeof(unsigned**));
  if ((a[0] = (unsigned**) malloc(((x + 1) * y) * sizeof(unsigned*))) == NULL)
    fatal("Memory block of %d bytes unavailable", (x + 1) * y * sizeof(unsigned*));
  if ((a[0][0] = (unsigned*) malloc((x * y * z) * sizeof(unsigned))) == NULL)
    fatal("Memory block of %d bytes unavailable", x * y * z * sizeof(unsigned));
  a[x] = NULL;
  for (i = 1; i < x; i++)
    a[i] = &a[0][y * i];
  for (i = 1; i < (x * y); i++)
    a[0][i] = &a[0][0][z * i];
#ifdef DebugInfo
  printf("Consuming %d bytes\n",   ((x + 1)       * sizeof(unsigned**))
                                 + (((x + 1) * y) * sizeof(unsigned*))
                                 + ((x * y * z)   * sizeof(unsigned)));
#endif
  return a;
}


unsigned ***init_u_array3(unsigned ***a, unsigned x, unsigned y, unsigned z, unsigned val)
{
  unsigned i, j, k;

  if (a == NULL)
    a = u_array3(x, y, z);
  for (i = 0; i < x; i++)
    for (j = 0; j < y; j++)
      for (k = 0; k < z; k++)
        a[i][j][k] = val;
  return a;
}


void destroy_u_array3(unsigned ***a)
{
  if (a) {
    free(a[0][0]);
    free(a[0]);
    free(a);
  }
}
