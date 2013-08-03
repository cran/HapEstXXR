/*******************************************************************************

  Create powersets of a given vector
  April, 1, 2010

*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "R.h"
#include "Rmath.h"
#include "Rinternals.h"
#include <R_ext/Parse.h>
#include <R_ext/PrtUtil.h>

SEXP powerset( SEXP a )
{
  int i,ii=0,j,jj,k,nprotect = 0 , nset=0 ;
  int temp[8000],tmp=0,dec ;

  SEXP ans,value2 ;
  nset = pow(2,length(a)) ;
  PROTECT(ans = allocVector(VECSXP, nset-1));
  nprotect++;
  for(i = 0; i < nset; i++) {
    /* search for */
    k=0 ;
    j=0;
    dec = i ;
    while ( dec>0 ) {
      tmp = dec % 2 ;
      dec=dec/2 ;
      if ( tmp==1 ) {
        temp[k++]=j;
      }
      j++;
    }
    if ( i!=0 ) {
      switch(TYPEOF(a)) {
        case INTSXP: PROTECT (value2 = allocVector(INTSXP,k)) ;   break;
        case REALSXP: PROTECT (value2 = allocVector(REALSXP,k)) ; break;
        case STRSXP: PROTECT (value2 = allocVector(STRSXP,k));  break;
        default: error("type of parameter x are not supported.");
      }
      nprotect++ ;
      for ( jj=0 ; jj<k ; jj++ ) {
        switch(TYPEOF(a)) {
          case INTSXP:  INTEGER(value2)[jj] = INTEGER(a)[temp[jj]] ; break;
          case REALSXP: REAL(value2)[jj] = REAL(a)[temp[jj]] ; break;
          case STRSXP:  SET_STRING_ELT(value2, jj,  STRING_ELT(a, temp[jj] )  ); break;
          default: error("type of parameter are not supported");
        }
      }
      SET_VECTOR_ELT( ans, ii++ , value2 );
    }
  }
  UNPROTECT(nprotect);
  return(ans) ;
}


SEXP powerset_file( SEXP a , SEXP filename )
{
  int i,ii=0,j,jj,k,nprotect = 0 , nset=0 ;
  int temp[8000],tmp=0,dec ;
  FILE *fp;
  fp=fopen(CHAR(STRING_ELT(filename,0)), "w");

  SEXP ans,value2 ;
  nset = pow(2,length(a)) ;
  PROTECT(ans = allocVector(VECSXP, nset-1));
  nprotect++;
  for(i = 0; i < nset; i++) {
    /* search for */
    k=0 ;
    j=0;
    dec = i ;
    while ( dec>0 ) {
      tmp = dec % 2 ;
      dec=dec/2 ;
      if ( tmp==1 ) {
        temp[k++]=j;
      }
      j++;
    }

    if ( i!=0 ) {
      switch(TYPEOF(a)) {
        case INTSXP: PROTECT (value2 = allocVector(INTSXP,k)) ;   break;
        case REALSXP: PROTECT (value2 = allocVector(REALSXP,k)) ; break;
        case STRSXP: PROTECT (value2 = allocVector(STRSXP,k));  break;
        default: error("type of parameter x are not supported.");
      }
      nprotect++ ;
      for ( jj=0 ; jj<k ; jj++ ) {
        switch(TYPEOF(a)) {
          case INTSXP:  INTEGER(value2)[jj] = INTEGER(a)[temp[jj]]             ;
                             fprintf (fp, "%i " , INTEGER(value2)[jj] ) ; break;
          case REALSXP: REAL(value2)[jj] = REAL(a)[temp[jj]]                   ;
                                fprintf (fp, "%f " , REAL(value2)[jj] ) ; break;
          case STRSXP:  SET_STRING_ELT(value2, jj,  STRING_ELT(a, temp[jj] )  );
                  fprintf (fp, "%s " , CHAR(STRING_ELT(a, temp[jj] )) ) ; break;
          default: error("type of parameter are not supported");
        }
      }
      SET_VECTOR_ELT( ans, ii++ , value2 );
    }
    fprintf(fp,"\n")  ;
  }

  fclose(fp) ;
  UNPROTECT(nprotect);

  return(ans) ;
}


void powerset_only_file( SEXP a , SEXP filename )
{
  int i,j,jj,k , nset=0 ;
  int temp[8000],tmp=0,dec ;
  FILE *fp;
  fp=fopen(CHAR(STRING_ELT(filename,0)), "w");
  nset = pow(2,length(a)) ;
  for(i = 0; i < nset; i++) {
    /* search for */
    k=0 ;
    j=0;
    dec = i ;
    while ( dec>0 ) {
      tmp = dec % 2 ;
      dec=dec/2 ;
      if ( tmp==1 ) {
        temp[k++]=j;
      }
      j++;
    }

    if ( i!=0 ) {
      for ( jj=0 ; jj<k ; jj++ ) {
        switch(TYPEOF(a)) {
          case INTSXP:  fprintf (fp, "%i " , INTEGER(a)[temp[jj]] ) ; break;
          case REALSXP: fprintf (fp, "%f " , REAL(a)[temp[jj]] ) ; break;
          case STRSXP:  fprintf (fp, "%s " , CHAR(STRING_ELT(a, temp[jj] )) ) ; break;
          default: error("type of parameter are not supported");
        }
      }
    }
    fprintf(fp,"\n")  ;
  }

  fclose(fp) ;
}


