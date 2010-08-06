  /* matalloc.h */
  /**************/


#define cmatrix(r, c)                       c_array2(r, c)
#define init_cmatrix(array, val, r, c)      init_c_array2(array, r, c, val)
#define destroy_cmatrix(mat)                destroy_c_array2(array)

#define dmatrix(r, c)                       d_array2(r, c)
#define init_dmatrix(array, val, r, c)      init_d_array2(array, r, c, val)
#define destroy_dmatrix(array)              destroy_d_array2(array)

#define imatrix(r, c)                       i_array2(r, c)
#define init_imatrix(array, val, c, r)      init_i_array2(array, r, c, val)
#define destroy_imatrix(r, c)               destroy_i_array2(array)

#define uimatrix(r, c)                      u_array2(r, c)
#define init_uimatrix(array, val, c, r)     init_u_array2(array, r, c, val)
#define destroy_uimatrix(array)             destroy_u_array2(array)


#define dvector(n)                          d_array(n)
#define init_dvector(array, val, n)         init_d_array(array, n, val)
#define is_in_dvector(array, val, n)        destroy_d_array(array)

#define ivector(n)                          i_array(n)
#define init_ivector(array, val, n)         init_i_array(array, n, val)
#define destroy_ivector(array)              destroy_i_array(array)

#define uivector(n)                         u_array(n)
#define init_uivector(array, val, n)        init_u_array(array, n, val)
#define destroy_uivector(array, val, n)     destroy_u_array(array)

/*****************************************************************************/
/*                                                                           */
/* Multidimensional arrays for type char                                     */
/*                                                                           */
/*****************************************************************************/

/*** array of char ***********************************************************/

char *c_array(unsigned n);

char *init_c_array(char *a, unsigned n, char val);

void destroy_c_array(char *a);

/*** 2-dimensional array of char *********************************************/

char **c_array2(unsigned x, unsigned y);

char **init_c_array2(char **a, unsigned x, unsigned y, char val);

void destroy_c_array2(char **a);

/*** 3-dimensional array of char *********************************************/

char ***c_array3(unsigned x, unsigned y, unsigned z);

char ***init_c_array3(char ***a, unsigned x, unsigned y, unsigned z, char val);

void destroy_c_array3(char ***a);

/*****************************************************************************/
/*                                                                           */
/* Multidimensional arrays for type double                                   */
/*                                                                           */
/*****************************************************************************/

/*** array of double *********************************************************/

typedef double (*ui2d_mapfun)(unsigned);

double *d_array(unsigned n);

double *init_d_array(double *a, unsigned n, double val);

double *dup_d_array(double *a, double *b, unsigned n);

double *assign_d_array(double *a, unsigned n, ui2d_mapfun f);

void destroy_d_array(double *a);

/*** 2-dimensional array of double *******************************************/

double **d_array2(unsigned x, unsigned y);

double **init_d_array2(double **a, unsigned x, unsigned y, double val);

void destroy_d_array2(double **a);

/*** 3-dimensional array of double *******************************************/

double ***d_array3(unsigned x, unsigned y, unsigned z);

double ***init_d_array3(double ***a, unsigned x, unsigned y, unsigned z, double val);

double ***assign_d_array3(double ***a, unsigned x, unsigned y, unsigned z);

void destroy_d_array3(double ***a);

/*****************************************************************************/
/*                                                                           */
/* Multidimensional arrays for type int                                      */
/*                                                                           */
/*****************************************************************************/

/*** array of int ************************************************************/

int *i_array(unsigned n);

int *init_i_array(int *a, unsigned n, int val);

int *assign_i_array(int *a, unsigned n, ui2d_mapfun f);

void destroy_i_array(int *a);

/*** 2-dimensional array of int **********************************************/

int **i_array2(unsigned x, unsigned y);

int **init_i_array2(int **a, unsigned x, unsigned y, int val);

void destroy_i_array2(int **a);

/*** 3-dimensional array of int **********************************************/

int ***i_array3(unsigned x, unsigned y, unsigned z);

int ***init_i_array3(int ***a, unsigned x, unsigned y, unsigned z, int val);

void destroy_i_array3(int ***a);

/*****************************************************************************/
/*                                                                           */
/* Multidimensional arrays for type unsigned                                 */
/*                                                                           */
/*****************************************************************************/

/*** array of unsigned *******************************************************/

unsigned *u_array(unsigned n);

unsigned *init_u_array(unsigned *a, unsigned n, unsigned val);

void destroy_u_array(unsigned *a);

/*** 2-dimensional array of unsigned *****************************************/

unsigned **u_array2(unsigned x, unsigned y);

unsigned **init_u_array2(unsigned **a, unsigned x, unsigned y, unsigned val);

void destroy_u_array2(unsigned **a);

/*** 3-dimensional array of unsigned *****************************************/

unsigned ***u_array3(unsigned x, unsigned y, unsigned z);

unsigned ***init_u_array3(unsigned ***a, unsigned x, unsigned y, unsigned z, unsigned val);

void destroy_u_array3(unsigned ***a);

