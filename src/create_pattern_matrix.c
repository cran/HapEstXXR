
#include <stdio.h>
#include <stdlib.h>
#include "R.h"

#define  pattern(row,col,nrow)   pattern[ (row) + ((col)*(nrow)) ]
#define  newpat(row,col,nrow)    newpat[ (row) + ((col)*(nrow)) ]
#define a(i,j,nrow)  a[(i)+((j)*(nrow))]

void  print_array (int *a , int *nrow , int *ncol ) {

  int i,j;
  Rprintf ("\n");
  for ( i=0 ; i<*nrow ; i++) {
    for ( j=0 ; j<*ncol ; j++ ) {
      Rprintf("%i ",a(i,j,*nrow));
    }
    Rprintf ("\n") ;
  }
} // end of print_array


void i_quicksort(int l[],int n) ;
void exist_pattern ( int * source  , int *ndim ,
                     int * pattern , int *nlen  ,
                     double * exist , int *search_to ) ;

void create_pattern_matrix ( int *pattern  , int *patterndim ,
                             int *snps     , int *snplen  ,
                             int * newpat  , int *newpatdim,
                             int * len ) {
  int i , ii, j , k ,  find , tmp[newpatdim[1]] , it ,
    search_to[1];
  double exist[1] ;

  // copy from pattern newpat ;
  ii = 0 ;
  for ( i=0 ; i<patterndim[0] ; i++ ) {
    for ( j=0 ; j<snplen[0] ; j++ ) {
      k = 0 ;
      find = 0 ;
      while ( (k<patterndim[1]) && (find==0) ) {
        // check if value exist in pattern ;
        if (pattern(i,k,patterndim[0])==snps[j]) {
          find=1 ;
        }
        k++;
      } // while k      
      if ( find==0 ) { 
        k=0;
        while ( k<patterndim[1] ) {
          newpat(ii,k,newpatdim[0]) = pattern(i,k,patterndim[0]) ;    
          k++;
        } // while k      
      
        newpat(ii,newpatdim[1]-1,newpatdim[0]) = snps[j] ;
        
        // exist ?? ---------------------------------------
        for ( it=0 ; it<newpatdim[1] ; it++ ) {
           tmp[it] = newpat(ii,it,newpatdim[0]) ;        

        }
        exist[0]=0 ;
        it=0 ;
         search_to[0] =  ii ;
          exist_pattern ( newpat  , newpatdim ,
                     tmp ,  &newpatdim[1] ,
                     exist , search_to ) ;
    
        if ( *exist==0 ) {
          // sort newpattern ;
          i_quicksort(tmp,newpatdim[1]-1);
          // save sorted pattern
          for ( it=0 ; it<newpatdim[1] ; it++ ) {
             newpat(ii,it,newpatdim[0])= tmp[it] ;        
          }   
          ii++ ;
        }
      } // if find ;      
    } // for j ;
  } // for i ;
  *len = ii ;
} // create_pattern_matrix



