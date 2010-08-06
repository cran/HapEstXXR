/*********************************************************************/
/*                                                                   */
/*    itegeppXXR.c (author Klaus Rohde) is haplotype estimation       */
/*    routine for samples of individual genotypes (EM-algorithm)     */
/*    The program calculates the weight for all haplotype pair       */
/*    configurations of the individuals and outputs a design matrix  */
/*    for the substantial (not lower weighted) haplotype pairs.      */
/*    If *tog is set to 1                                            */
/*    the program prints out a design matrix for each individual     */
/*    by the weights of the possible haplotype pairs                 */
/*    If *tog is set to 0                                            */
/*    the program prints out a design matrix over haplotypes         */
/*    Parameter lim limits the essential haplotypes to freq>=lim     */
/*    Input may be simulated via prodtdt (prodtdt.c)                 */
/*                Version 01.04.2006                                 */
/*                                                                   */
/*********************************************************************/
/* Version history:

0.1-1: June 07, 2007
	- Maximum naumber of SNPs is limited to 16.
	- If des=0 (haplotypes) then will be maximum number of haplotypes limited to 100.
	- If des=1 (haplotype pairs) then will be maximum number of haplotypes limited to 60.

0.1-2: 2008
	- Change in likea():
	  Computation of the likelihood is now: f += ((double) mg[i]) * log(prob[i]) ;
  - Change in call of itegeppXR  !!!
    alternate char **likeres to double *likeres 
  - Änderung Ausgabeformat: desres auf 8.6f  und  hapres auf 9.6f
  - Änderung Designmatrix für tog=0
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fundamental.h"
#include "bool.h"
#include "matalloc.h"
#include "message.h"
#include "readin.h"
#include "R.h"

/* 2^16=65536  --> max. 16 SNPs */
#define Hapco                 65537
#define SignificanceLevel      0.001
#define LoopPrecision          0.00001
#define LineLength            80

#define Epsilon                0.000001


enum EstimationMethod { BestState = 0 , Iteration = 1 };


char    **geno;           /* Auflistung tatsaechlch vorkommender Genotypen  */
double   *hap,            /* Prob. tats. vorkomm. Haplotypen                */
         *haptmp,         /* temp. Prob. tats. vorkomm. Haplotypen          */
 	 *hapnew,         /* temp. Prob. tats. vorkomm. Haplotypen          */
 	 *prob,           /* Prob. tats. vorkomm. Genotypen                 */
         *genoProb,       /* Prob. tats. vorkomm. Genotypen equal prior pr. */
// 	 *p,              /* probabilities used for simulation              */
         *max_prob,       /* maximal probability of each genotype           */
//         *qtrait,         /* quantitative traits                            */
	  p_value,        /* p-value of anova                               */
        **pgen;
int       first,          /* flaf for simulations                           */
         *genoId,         /* list of sorted genotypes                       */
          nall,           /* number of alleles in sample                    */
          np,             /* number of probands in the input file           */
          nh,             /* number of tats. vorkomm. Haplotypen            */
         *nstate,         /* number of possible states per genotype         */
//         *instate,        /* number of possible states per individual       */
 	 *mstate,         /* nu of possible states per genotype best state  */
          ng,             /* number of unique haplotypes which we extracted */
         *merke,          /* number of heterozygous loci per genotype       */
         *nulmer,         /* number of unknown loci per genotype            */
          len,            /* string length of the phenotype description     */
         *mg,             /* number of tats. vorkomm. Genotypes             */
         *hc,             /* number of tats. vorkomm. Haplotypes            */
 	 *ge,             /* pointer from person to genotype                */
  	  Loci,
          People,
          jjx;            /* design output count                            */
uint   ***state,          /* list haplotypes pairs per genotype             */
        **htpp,
	  *po;            /* list of 2 to power                             */


void pspace(int n, char *hpres)
{
	int i;
	char hp[40];
	sprintf(hp,"   ***");                /*** ROHDE 11.09.2001 ***/
	for(i=0;i<n;i++)sprintf(hp+6+i," ");
	hp[n+6] = '\0';
	strcat(hpres,hp);
}

void swap_i(int *u, int *v)
{
  int w = *u;   *u = *v;   *v = w;
}


void swap_ui(uint *u, uint *v)
{
  uint w = *u;   *u = *v;   *v = w;
}


void swap_d(double *u, double *v)
{
  double w = *u;   *u = *v;   *v = w;
}


int fill(int n)
{
  int i;

  for (i=0; i<nh; i++)
    if ( hc[i] == n )
      return i;
  hc[nh++] = n;
  return nh-1;
}


bool compatible(char *g, int *ac)
{
  int i;

  for (i=0; i<len; i++)
    switch ( g[i] ) {
      case '1':
        if ((ac[0] & po[i]) || (ac[1] & po[i]))
          return false;
        break;
      case '2':
        if ( ! (ac[0] & po[i]) || ! (ac[1] & po[i]))
          return false;
        break;
      case '3':
        if ((ac[0] & po[i]) == (ac[1] & po[i]))
          return false;  /* '0' is always compatible */
    }
  return true;
}


void addon(int i)
{
  int    j, k;
  double weight = (double) mg[i] / (double) nstate[i] / (double) nall;

  for (j=0; j<nstate[i]; j++)
    for (k=0; k<2; k++)
      hap[state[i][j][k]] += weight;
}


void recprob(int i)
{
  int         j, k, zz;
  double      prod;
  uint        ac[2];

  prob[i] = 0.0;
  zz      = nstate[i];
  for (j=0; j<zz; j++) {
    ac[0] = state[i][j][0];
    ac[1] = state[i][j][1];

    prod     = hap[ac[0]] * hap[ac[1]];
    prob[i] += prod;
    /*** begin computation of mean probabilities *****************************/
    if ( ac[0] == ac[1] )
      haptmp[ac[0]] += 2.0 * prod;
    else
      for (k=0; k<2; k++)
        haptmp[ac[k]] += prod;
  }
}


void selprob(int n)
{
  int    i, j,
         maxi,
         k         = 0,
         zz;
  double prodmax   = 0,
         prodfirst = 0,
         prod;

  zz        = nstate[n];
  mstate[n] = 0;
  for (i=0; i<zz; i++) {
    maxi    = i;
    prodmax = hap[state[n][i][0]] * hap[state[n][i][1]];

    for (j=i+1; j<zz; j++) {
      prod = hap[state[n][j][0]]*hap[state[n][j][1]];
      if (prod > prodmax) {
        maxi    = j;
        prodmax = prod;
      }
    }     /* end j-loop */

    if ( k == 0 ) {
      prodfirst = prodmax;
      k         = 1;
    }
    else
      if ( fabs(prodfirst - prodmax)/prodfirst > Epsilon )  break;

    swap_ui(&state[n][maxi][0], &state[n][i][0]);
    swap_ui(&state[n][maxi][1], &state[n][i][1]);
    mstate[n] ++;
  }     /* end i-loop */
  max_prob[n] = prodfirst;
}


double likea(void)
{
  int    i, j;
  double f = 0.0;

  for (i=0; i<ng; i++) {
    init_dvector(haptmp, 0.0 , nh);
    recprob(i);
    if (first == 0)
    
    // change computation of Likelihood   SKn
    
      for (j=0; j<nh; j++)
        if(prob[i] > 0.0) hapnew[j] += haptmp[j] * (double)mg[i] / prob[i];
    if ( prob[i] <= 0.0 )
      f += 10000.0;
    else
      f += ((double) mg[i]) * log(prob[i]) ;
  }
  return f;
}


void rechap(int i, int k1, int k2)
{
  static unsigned int m, ac[2];
  int                 j, k, hapi;
  double              vv;

  k = k1;
  if (k == 0) {
    ac[0] = ac[1] = 0;
    m     = 0;
  }
  if ( k == k2 ) {
    vv                = nulmer[i] ? (1.0/pow(4.0,nulmer[i])) : 1.0;
    vv               *= (double)mg[i]/(double)po[merke[i]]/(double)nall;
    hapi              = fill(ac[0]);
    state[i][m][0]    = hapi;
    hap[hapi]        += vv;
    hapi              = fill(ac[1]);
    state[i][m++][1]  = hapi;
    hap[hapi]        += vv;
    nstate[i]         = m;
  }
  else {
    if ((geno[i][k] == '1') || (geno[i][k] == '0')) {
      if ( ac[0] & po[k] )         ac[0] -= po[k];
      if ( ac[1] & po[k] )         ac[1] -= po[k];
      rechap(i, k+1, k2);
    }
    if ((geno[i][k] == '2') || (geno[i][k] == '0')) {
      if ( !(ac[0] & po[k]) )      ac[0] += po[k];
      if ( !(ac[1] & po[k]) )      ac[1] += po[k];
      rechap(i, k+1, k2);
    }
    if ((geno[i][k] == '3') || (geno[i][k] == '0')) {
      for (j=0; j<2; j++) {
        if (j == 0) {
          if ( ac[0] & po[k] )     ac[0] -= po[k];
          if ( !(ac[1] & po[k]) )  ac[1] += po[k];
        }
        else {
          if ( !(ac[0] & po[k]) )  ac[0] += po[k];
          if ( ac[1] & po[k] )     ac[1] -= po[k];
        }
        rechap(i, k+1, k2);
      }
    }
  }
}


void printHaplotype(int haplo, int len, char *hpres)
{
  int i;
  char hp[len];

  for (i=0; i<len-1; i++)
    sprintf(hp+i,"%c", ((haplo & po[i]) ? '2' : '1'));
    *(hp+len-1)='\0';
    strcat(hpres,hp);
}


void sortByProb(double *prob, int *geno, int n)
{
  int    i, j;

  for (i=1; i<n; i++)
    for (j=i; (j > 0) && (prob[j] > prob[j-1]); j--) {
      swap_d(&prob[j], &prob[j-1]);
      swap_i(&geno[j], &geno[j-1]);
    }
}



void itegeppXXR(int *tog, double *lim, char **gent, double *qtrait,
int *xnp, double *likeres, char **freqres, char **hapres, char **desres)
{
  char      lino[10000], lin[10000];
  double    likold,
            pe, pex,                 /* 10.3. 2000 ROHDE */
	    *p2max,
            gsum;                    /* 10.3. 2000 ROHDE */
  int       i, inp, it, j, k, ki, kj, h, s, glev, non,
            ac[2],
	    drei,
	    null,
	    df,
	    combinations,
	    nz,
	    iqual,
            nhap,
           *hlist,
	    **pimax,
            h1x, h2x;
  uint      iterations, h1, h2;
  bool      loop;

  // new for create design matrix (tog=0)
  double  pehh,
          *peh;
  

/* Max. 16 SNPs */
 if ( strlen(gent[0]) > 16 ) error ("Number of SNPs should smaller than 17.") ;

  np       = *xnp;
  len      =  strlen(gent[0]) + 1;

  mg       = ivector(np);
  merke    = ivector(np);
  nulmer   = ivector(np);
  ge       = ivector(np);
  hlist    = ivector(np);
  po       = uivector(len);

  geno     = cmatrix(np, len);
  max_prob = init_dvector(NULL, 0.0, np);
  prob	   = init_dvector(NULL, 0.0, np);


  hap	   = init_dvector (NULL,  0.0, Hapco);
  hc	   = init_ivector (NULL, -1,Hapco);

  po[0]=1;
  for(i=1;i<len;i++)po[i] = 2*po[i-1];
  combinations = po[len-1];


    init_dvector(hap,  0.0, Hapco);
    init_ivector (hc, -1,Hapco);

    ng = 0;

    /* read input data */
    for(inp=0;inp<np;inp++){
    drei = 0;
    null = 0;
    for (i=0; i<len-1; i++) {
    if(i < len-1 && (gent[inp][i] < 48 || gent[inp][i] > 51) ){
    Rprintf("%d %d %d\n",inp, i, gent[inp][i]);
    //Rprintf("\n Error in data person %d\n",inp+1); /* ROHDE 15.03.2000 */
    error("\n Error in data person %d\n",inp+1);;
    }
    if ( gent[inp][i] == '3' )  drei  ++;
    if ( gent[inp][i] == '0' )  null  ++;
    }
    gent[inp][len-1] = '\0';

    it = 1;
    for (i=0; i<ng; i++) {
    if ( strncmp (geno[i], gent[inp], len) == 0 ) {

    /*** a certain genotype was found more than just once ***/
    ge[inp] = i;
    mg[i] ++;
    it       = 0;
    merke[i] = drei;
    break;
    }
    }
    if (it) {
    /*** a certain genotype was encountered the first time ***/
    strcpy (geno[ng], gent[inp]);
    ge[inp] = ng;
    mg[ng]    = 1;
    merke[ng] = drei;
    nulmer[ng] = null;
    ng ++;
    }
    }   /* end while */
    People 	= np;
    Loci 	= len-1;
    /* end of reading sample data */

    nall = 2 * np;
    nstate  = init_ivector (NULL, 0, ng);
    mstate  = init_ivector (NULL, 0, ng);

    state = (uint***) calloc(ng , sizeof(uint**));
    for (i=0; i<ng; i++) {
      nz       = po[merke[i]] * po[nulmer[i]] * po[nulmer[i]];
      state[i] = uimatrix(nz, 2);
    }
    /*** sort genotypes by weights *******************************************/
    genoProb = dvector(ng);
    genoId   = ivector(ng);
    for (i=0; i<ng; i++) {
      genoId[i]   = i;
      genoProb[i] = ((double)mg[i])/((double) po[merke[i]])/pow(4.0,nulmer[i]);
    }
    sortByProb(genoProb, genoId, ng);

    glev=0;
    for(i=0;i<ng;i++)if(genoProb[i] >= SignificanceLevel)glev++;

    /*** process sorted genotypes ********************************************/

    nh = 0;

    for (i=0; i<glev; i++) {
/*    printf("\n ng: %d glev: %d  i: %d",ng,glev+1,i+1);    */
      rechap(genoId[i], 0, len-1);
/*    printf("\n %s >> %d\n",geno[genoId[i]],mg[genoId[i]]);
      for(k=0;k<16;k++)printf("%2d:%g ",hc[k],hap[k]);
      printf("\n");                                         */
    }

    for (i=glev; i<ng; i++) {
      s = 0;
/*    printf("\n ng: %d glev: %d  i: %d",ng,glev+1,i+1);    */
      for (j=0; j<nh; j++) {
        ac[0] = hc[j];
        for (k=j; k<nh; k++) {
          ac[1] = hc[k];
	  if ( compatible(geno[genoId[i]], ac) ) {
            state[genoId[i]][s][0] = j;
            state[genoId[i]][s][1] = k;
            s ++;
            if ( j != k ) {
              state[genoId[i]][s][0] = k;
              state[genoId[i]][s][1] = j;
              s ++;
	    }
	  }
        }
      }
      nstate[genoId[i]] = s;
    }
    for (i=glev; i<ng; i++) {
      addon(genoId[i]);
    }

/*  printf("\n");
    printf("\ngloop: %d ng: %d glev: %d nh: %d\n",gloop,ng,glev,nh);  */

    /*** now comes the output that does not need simulated annealing *********/

    first  = 1;           /*** start likelihood outside of annealing loops ***/

    df = nh;

    hapnew = init_dvector(NULL, 0.0, nh );
    haptmp = init_dvector(NULL, 0.0, nh );

    for (i=0; i<ng; i++)selprob(i);

    likold = likea();

    /* Continue computation of mean probabilities */

    for(i=0;i<ng;i++)
      for(j=0;j<mstate[i];j++) {
        double pp = 0.0;

        if ( nstate[i] > 1 ) {
          h          = state[i][j][0];
          hapnew[h] += (double)mg[i] / (double)mstate[i];
          h          = state[i][j][1];
          pp        += hapnew[h];
          hapnew[h] += ((double) mg[i]) / ((double) mstate[i]);
          pp        += hapnew[h];
        }
        else {
          h          = state[i][j][0];
          hapnew[h] += 2.0 * ((double) mg[i]) / ((double) mstate[i]);
          pp        += hapnew[h];
        }
      }


    non = 0;
    for (i=0; i<nh; i++) {
    if(hapnew[i]==0.0)non++;
    else hapnew[i] /= (double) nall;
    }

    for (i=0; i<nh; i++) {
    if(hapnew[i]==0.0)hapnew[i] = 0.0001/(double)non;
    else hapnew[i] *= 0.9999;
    }


	iterations = 0;
    first = 0;

    do {
      loop = 0;
      iterations ++;
/*    printf("gloop:%3d  count: %d\n",gloop,iterations);  */
      /* Recompute mean probabilities */
      for (i=0; i<nh; i++) {
        if ( fabs(hap[i] - hapnew[i]) > LoopPrecision )
          loop = 1;
        hap[i] = hapnew[i];
      }

      init_dvector(prob,   0.0, np);
      init_dvector(haptmp, 0.0, nh);
      init_dvector(hapnew, 0.0, nh);
      likold = likea();
      for (i=0; i<nh; i++)
        hapnew[i] /= (double) nall;

    } while (loop);

/*  Rprintf("\n");
    Rprintf("  Results Ensemble means: \n\n"); */
    nhap = 0;
    j    = 0;
    for (i=0; i<nh; i++)
      {
      if ( hapnew[i] >= *lim ) {
    /* 07.06.2007  S.Knüppel > Beschränken der geschätzten Haplotypen. */
    if ( (*tog==0) &&  ((nhap+1) > 100) ) {
     error ("Error in itegeppXXR: Too much estimated haplotypes. Increase option lim.") ;
    }
    if ( (*tog==1) &&  ((nhap+1) > 60) ) {
     error ("Error in itegeppXXR: Too much estimated haplotypes. Increase option lim.") ;
    }



  	sprintf(lino,"\0");
        printHaplotype(hc[i], len, lino);
	/*
        printf("    hapnew[%8d] = %7.4f  (%7.4f)\n",
	    hc[i], hapnew[i], hap[i]);
        */
	sprintf(lin,"%9.6f\0", hapnew[i]);
	strcat(lino,lin);
	strcpy(freqres[j],lino);
	j++;
        hlist[nhap++] = i;



      }
      }
      k = 0;
      htpp = init_uimatrix(NULL,0,nhap+1,nhap+1);
      for(i=0;i<nhap+1;i++)
	      for(j=i;j<nhap+1;j++)htpp[i][j]=k++;

      pgen = init_dmatrix(NULL,0.0,ng,(nhap+1)*(nhap+2)/2);

      /* start find best states after MLE   10.3.2000 ROHDE  */

   	  pimax  = imatrix(ng,10);               /* ROHDE  10.3.2000 */
          p2max = init_dvector(NULL,0.0,ng);     /* ROHDE  10.3.2000 */
	  for(i=0;i<ng;i++)max_prob[i] = 0.0;    /* ROHDE  10.3.2000 */
      for (i=0;i<ng;i++){
      for(j=0;j<10;j++)pimax[genoId[i]][j] = -1;
          iqual=1;
      for (j=0;j<nstate[genoId[i]];j++){
      pe = hapnew[state[genoId[i]][j][0]] * hapnew[state[genoId[i]][j][1]];
	  if( state[genoId[i]][j][0] != state[genoId[i]][j][1] ) pe += pe;

	  if(pe > p2max[genoId[i]]){

      if (pe > max_prob[genoId[i]]){
		p2max[genoId[i]] = max_prob[genoId[i]];
      	max_prob[genoId[i]] = pe;
      	pimax[genoId[i]][0]=j;
		for(k=1;k<10;k++)pimax[genoId[i]][k]=-1;
		iqual = 1;                 /***   ROHDE  04.09.2001 ***/
      	}

	  else{
	  if (pe == max_prob[genoId[i]] && iqual < 9){
	  	for(k=0;k<iqual;k++) if(state[genoId[i]][j][0] ==
						state[genoId[i]][pimax[genoId[i]][k]][1]) pe=0.0;
	  	if(pe > 0.0)pimax[genoId[i]][iqual++]=j;
	  	}
	    else p2max[genoId[i]] = pe;
		}
	  }

      }
      }    /* end of maximum state search */

/*    Rprintf("\n Haplotypes after MLE\n");           */
      jjx = 0;
      for(i=0;i<np;i++){
        sprintf(lino,"%i %s >> \0",i, geno[ge[i]]);
        for(k=0;k<10;k++){
          j = pimax[ge[i]][k];
          if(j > -1){
            if(k>0)pspace(len+3,lino);   /*** ROHDE  11.09.2001 ***/
            printHaplotype(hc[state[ge[i]][j][0]],len,lino);
            strcat(lino," <> \0");
            printHaplotype(hc[state[ge[i]][j][1]],len,lino);
            sprintf(lin,"  P>> %9.7f D>> %9.7f",
      	                 max_prob[ge[i]],max_prob[ge[i]]-p2max[ge[i]]);
            strcat(lino,lin);
          } else break;
        }
        strcpy(hapres[jjx++],lino);
      }

      /* endfind best states after MLE   10.3.2000 ROHDE  */


/*    Rprintf("\n\n       Likelihood = %f\n", likold);
      Rprintf("\n");                                      */
      sprintf(lino,"Likelihood = %f\0", likold);
     // strcpy(likeres[0],lino);
      (*likeres) = likold ;
      
/*  Sample over states for each genotype ***********************************/

      for(i=0;i<ng;i++){
      gsum = 0.0;
      for(j=0;j<nstate[genoId[i]];j++){
      h1 = state[genoId[i]][j][0];
      h2 = state[genoId[i]][j][1];
      h1x = h2x = 0;
      for(ki=1;ki<=nhap;ki++) if( h1 == hlist[ki-1] ) h1x=ki;
      for(kj=1;kj<=nhap;kj++) if( h2 == hlist[kj-1] ) h2x=kj;
      if(h1x>0 && h2x>0){
	      if(h2x < h1x){ k=h1x; h1x=h2x; h2x=k;}
	      pgen[genoId[i]][htpp[h1x-1][h2x-1]] += hapnew[h1]*hapnew[h2];
	                                     gsum += hapnew[h1]*hapnew[h2];
	}
	else{ pgen[genoId[i]][htpp[nhap][nhap]] += hapnew[h1]*hapnew[h2];
		                           gsum += hapnew[h1]*hapnew[h2];
	}
      }
      for(k=0;k<(nhap+1)*(nhap+2)/2;k++)pgen[genoId[i]][k] /= gsum;
      }

/*    for(i=0;i<ng;i++){
      Rprintf("i:%2d  %s\t",i,geno[genoId[i]]);
      for(ki=0;ki<nhap+1;ki++){
	for(kj=ki;kj<nhap+1;kj++)
               Rprintf("%4.2f ",pgen[genoId[i]][htpp[ki][kj]]);
	if(kj<ki)printf("0.000\t ");
	else	printf("%4.2f\t",pgen[genoId[i]][htpp[ki][kj]]);
      Rprintf("\t");
      }
      Rprintf("\n");
      }
*/

      jjx = 0;

      if (*tog == 1){
      for(i=0;i<np;i++){
      /*
      printf("\n%4s %s %4.2f >> ",pid[i],geno[ge[i]],qtrait[i]);
      */
      strcpy(lino,"\0");
      for(ki=0;ki<nhap;ki++){  /*  each haplotype alone   */
        for(kj=ki;kj<nhap;kj++){
          sprintf(lin,"%8.6f \0",pgen[ge[i]][htpp[ki][kj]]);
	        strcat(lino,lin);
	      }

      }
      sprintf(lin,"%8.6f\0",pgen[ge[i]][htpp[nhap][nhap]]);
      strcat(lino,lin);
      strcpy(desres[jjx],lino);
      jjx++;
      }
      }
      
 /* geändert nach Klaus; Bildung Designmatrix
    16.09.2008 */
 
 /*     if(*tog == 0){
      for(i=0;i<np;i++){
      //
      //printf("\n%4s %s %4.2f >> ",id[i],geno[ge[i]],qtrait[i]);
      //
      strcpy(lino,"\0");
          pex = 0.0;
      for(j=0;j<nhap;j++){
          pe = 0.0;
      for(ki=0;ki<nhap;ki++){    // over all haplotype pairs  
      for(kj=ki;kj<nhap;kj++){
      if(ki==j && kj==j && pgen[ge[i]][htpp[ki][kj]] > 0.0)
	      pe +=2.0*pgen[ge[i]][htpp[ki][kj]];
      else if ((ki==j || kj==j) && pgen[ge[i]][htpp[ki][kj]] > 0.0)
                                    pe += pgen[ge[i]][htpp[ki][kj]];
      }
      }
      pex += pe;
      sprintf(lin,"%8.6f \0",pe);
      strcat(lino,lin);
      }
      sprintf(lin,"%8.6f\0",2.0-pex);
      strcat(lino,lin);
      strcpy(desres[jjx],lino);
      jjx++;
      }
      }
*/
/* new: nach Klaus; 17.09.2008 */
        if(*tog == 0){
        
      peh = init_dvector(NULL, 0.0, nhap+1);
      for(i=0;i<np;i++){
      /*
      printf("\n%4s %s %4.2f >> ",id[i],geno[ge[i]],qtrait[i]);
      */
      strcpy(lino,"\0");
      for(j=0;j<nhap;j++){
      gsum  = 0.0;
      /*
      for(ki=0;ki<nhap;ki++){    * over all haplotype pairs  *
      for(kj=ki;kj<nhap;kj++){
      if(ki==j && kj==j && pgen[ge[i]][htpp[ki][kj]] > 0.0)
	      pe +=2.0*pgen[ge[i]][htpp[ki][kj]];
      else if ((ki==j || kj==j) && pgen[ge[i]][htpp[ki][kj]] > 0.0)
                                    pe += pgen[ge[i]][htpp[ki][kj]];
      }
      }
      */
      for(ki=0;ki<nstate[ge[i]];ki++){
      h1 = state[ge[i]][ki][0];
      h2 = state[ge[i]][ki][1];
      h  = hlist[j];
      pex = hapnew[h1]*hapnew[h2];
      gsum += 2*pex;
      if((h == h1) && (h == h2))peh[j] += 2*pex;
      else if((h == h1) || (h == h2))peh[j] += pex;
      }  /* end nstate */
      }  /* end nhap   */
      pehh = 0.0;
      for(j=0;j<nhap;j++){
      pehh += 2*peh[j];
      sprintf(lin,"%8.6f \0",2*peh[j]/gsum);
      strcat(lino,lin);
      }  /* end print */
      sprintf(lin,"%8.6f\0",2.0-pehh/gsum);
      strcat(lino,lin);
      strcpy(desres[jjx],lino);
      jjx++;
      init_dvector(peh, 0.0, nhap+1);
      }  /*  end np  */
      destroy_d_array(peh);
      }

    
    destroy_c_array2(geno);
    destroy_u_array(po);
    for ( i=0;i<ng;i++) { destroy_u_array2(state[i]) ; }
    free((uint***)state);
    destroy_u_array2(htpp);
    destroy_i_array(nstate);
    destroy_i_array(mstate);
    destroy_i_array(genoId);
    destroy_i_array(mg);
    destroy_i_array(merke);
    destroy_i_array(nulmer);
    destroy_i_array(ge);
    destroy_i_array(hlist);
    destroy_d_array(prob);
    destroy_d_array(max_prob);
    destroy_d_array(hapnew);
    destroy_d_array(haptmp);
    destroy_d_array(hap);
    destroy_d_array(genoProb);
    destroy_d_array2(pgen);
    destroy_i_array(hc);
    destroy_d_array(p2max);
    destroy_i_array2(pimax);
}
