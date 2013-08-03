/*******************************************

  matrix.c

********************************************/


#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include "matrix.h"

/*  'i'   integer
    length
    ' '   format 'c' character , 'f' double , 'i' integer
    array
    'e'   end of input
    example:
    print_vec  ( 'i' , nb , 'i' , b , 'e' ) ;
*/
void print_vec ( char format , ... ) {      
  va_list p ;                                            
  int i,n ;
  va_start ( p , format ) ;
  format = n = va_arg ( p , int ) ;
  format = va_arg ( p , int ) ;
  Rprintf ( "\n" ) ;
  while ( format != 'e' ) {
    switch ( tolower(format) ) {
      case 'i': {
         int *a = va_arg ( p , int* ) ;
         for ( i=0 ; i<n ; i++ ) {
           Rprintf ( "%i  " , a[i] ) ;
           if ( ((i+1)%10) == 0 ) Rprintf("\n");
         }
         break;
      }
      case 'f': {
        double *a = va_arg ( p , double* ) ;
         for ( i=0 ; i<n ; i++ ) {
           Rprintf ( "%f " , a[i] ) ;
           if ( ((i+1)%10) == 0 ) Rprintf("\n");
         }
         break;
      }
      case 'c': {
        char *a = va_arg ( p , char* ) ;
         for ( i=0 ; i<n ; i++ ) {
           Rprintf ( "%c" , a[i] ) ;
           if ( ((i+1)%10) == 0 ) Rprintf("\n");
         }
         break;
      }
    }
    format = va_arg(p,int) ;

  }
  Rprintf ( "\n" ) ;
  va_end ( p ) ;

} /** End of print_vec **/

/*******************************************************************
 *  Functions                                                      *
 *                                                                 *
 *******************************************************************/
char** c_matrix2 ( int x , int y ) {
  char   **a;
  int   i;

  if ( (a = (char**) malloc((x) * sizeof(char*))) == NULL )
    error("Error in c_matrix2: Memory block of %d bytes unavailable.", ((x)*sizeof(char*)));

   for (i = 0; i < x; i++) {
       if ((a[i] = (char*) malloc( y *  sizeof(char))) == NULL)
         error("Error in c_matrix2: Memory block of %d bytes unavailable.", y*sizeof(char));
   }

  return a;

} /** End of c_matrix2 **/

void init_c_matrix2 ( char **a , int x , int y , char value ) {

  int   i,j ;
  for ( i=0;i<x;i++ ) { for ( j=0;j<y;j++ ) a[i][j] = value ;}

} /** End of init_c_matrix2 **/

void free_c_matrix2 ( char **a , int x  ) {
    int   i ;
    for (i = 0; i < x; i++) {
      free(a[i]) ;
   }
   free ((char**)a) ;


} /** End of free_c_matrix2  **/

void print_c_matrix2 ( char **a , int x , int y ) {
  int   i,j ;
  Rprintf("\n");
  for ( i=0;i<x;i++ ) {
    for ( j=0;j<y;j++ )  {
      Rprintf ( "%c" , a[i][j] ) ;
    }
    Rprintf("\n");
  }
  Rprintf("\n");
} /** End of print_c_matrix2 **/


// Convert binary-to-decimal
// e.g. '222'=0    '211' = 3     '111' = 7
// a = '222'
// res = 0
int Bin2Dec ( char *a , int len ) {
  int res=0 , i ;
  for ( i=0 ; i<len ; i++) { if ( a[len-i-1]=='1' ) {res +=  (int) pow (2,i ) ; } }
  return ( res ) ;
} /** End of Bin2Dec **/

