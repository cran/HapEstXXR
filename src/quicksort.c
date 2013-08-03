/******************************************************
 quicksort.c
 2007-08-22
******************************************************/

#include <string.h>
#include <stdio.h>
#include "R.h"

void d_swap ( double *a , int i, int j )
{
  double tmp ;
  tmp = a[i] ;
  a[i] = a[j] ;
  a[j] = tmp ;
}  /* End d_swap */

void c_swap ( char **a , int i, int j )
{
  char *tmp ;
  tmp = a[i] ;
  a[i] = a[j] ;
  a[j] = tmp ;
}  /* End c_swap */


void quicksort ( double *a , char **b  , char **c ,  int links , int rechts )
{
   double pivot ;
   int l, r;
   if (links < rechts )
   {
      l = links; r = rechts;
      pivot = a[ ( (int) ( (links+rechts)/2 ) ) ];
      do
      {
         while  ( a[l] < pivot ) l++;
         while  ( a[r] > pivot ) r--;
         if (l <= r )
         {
            d_swap ( a , l , r ) ;
	       		c_swap ( b , l , r ) ;
            c_swap ( c , l , r ) ;
            l++; r--;
         }
      }
      while ( l <= r );
      quicksort(a,b,c, links, r);
      quicksort(a,b,c, l, rechts);
   }
} /* End quicksort */

