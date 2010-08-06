
#include <stdio.h>
#include <stdlib.h>
#include "R.h"


#define  source(row,col,nrow)   source[ (row) + ((col)*(nrow)) ]

void exist_pattern ( int * source  , int *ndim ,
                     int * pattern , int *nlen  ,
                     double * exist , int *search_to ) {

  int i , j , k , count , find , in ,
      nr = ndim[0] ,
      nc = ndim[1] ;

  if ( nc != (*nlen) ) {
     error("\n Error in exist_pattern: dimension of source (nc=%i) clashed with length of pattern (nlen=%i) \n" , nc , *nlen );
  }
  *exist = 0 ;
  i=0 ;
  while ( (search_to[0]>=0) & (i<search_to[0]) & ((*exist)!=1) ) {
     count = 0 ;
     
     j=0 ;
     in = 1 ;
     while ( j<(*nlen) && (in==1) ) {
       // find pattern in source matrix
       k = 0 ;
       find = 0 ;
       while ( (k<(*nlen)) & (find==0) ) {
          if ( source(i,j,nr) == pattern[k] ) {
            find=1 ;
            count++;
          }
          k++ ;
       } // while k
       if (find==0) {
         in = 0;
       }
       j++;
     } // while j
     if ( count == (*nlen) ) {
       (*exist) = 1 ;
     }
    i++ ;
  } // while i          
} // end of exist_pattern

