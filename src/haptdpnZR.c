/*********************************************************************/
/*                                                                   */
/*    haptdpnR.c (author Klaus Rohde) is haplotype estimation        */
/*    routine with inclusion of nuclear family information,          */
/*    and subsequent TDT analysis according Lazzaroni/Lange          */
/*    The program is written as R-function for inclusion             */
/*    The program calculates the weight for all haplotype pair       */
/*    configurations of a family and carries out a subsequent        */
/*    TDT statistics taking the weights into account                 */
/*    A parent pair may have no or up to ten children. Every         */
/*    pair of unrelated people may be taken as as parent pair        */
/*    without children. Unknowns (0) in the genotypes are not        */
/*    allowed. Best states are found according to MLE.               */
/*    Parameter lim limits the essential haplotypes to freq>=lim     */
/*    Additionally the program prints out a design matrix for        */
/*    each individual by the weights of the possible haplotype       */
/*    pairs. Input may be simulated via prodtdt (prodtdt.c)          */
/*                Version 01.04.2006                                 */
/*                                                                   */
/*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "readin.h"
#include "fundamental.h"
#include "bool.h"
#include "matalloc.h"
#include "message.h"
#include "rng.h"
#include "R.h"

#define NoSeparateRNGs
#define Loops                100
#define RandomSeed         12345
#define Hapco              40000
#define SignificanceLevel      0.001
#define LoopPrecision          0.001
#define LineLength            80

#define Epsilon                0.001
#define perLoop(x)             ((x) / ((double) Loops))
#define InputFile               "platz.dat"


enum EstimationMethod { BestState = 0 , Iteration = 1 };


char    **geno;           /* Auflistung tatsaechlch vorkommender Genotypen  */
double   *hap,            /* Prob. tats. vorkomm. Haplotypen                */
         *haptmp,         /* temp. Prob. tats. vorkomm. Haplotypen          */
	 *hapnew,         /* temp. Prob. tats. vorkomm. Haplotypen          */
	 *prob,           /* Prob. tats. vorkomm. Genotypen                 */
         *genoProb,       /* Prob. tats. vorkomm. Genotypen equal prior pr. */
	**pep,            /* probabilities of states for each family        */
         *max_prob;       /* maximal probability of each genotype           */
int       first,          /* flaf for simulations                           */
	  nall,           /* number of alleles in sample                    */
          np,             /* number of probands in the input file           */
	  nf,             /* number of families in input file               */
          nh,             /* number of tats. vorkomm. Haplotypen            */
         *nstate,         /* number of possible states per genotype         */
	 *nfstate,        /* number of states per family                    */
	 *mstate,         /* nu of possible states per genotype best state  */
          ng,             /* number of unique haplotypes which we extracted */
         *merke,          /* number of heterozygous loci per genotype       */
         *nulmer,         /* number of unknown loci per genotype            */
          len,            /* string length of the phenotype description     */
         *mg,             /* number of tats. vorkomm. Genotypes             */
         *hc,             /* number of tats. vorkomm. Haplotypes            */
	 *ge,             /* pointer from person to genotype                */
	  itp, itn,
	  InitialHaplos,
	  Loci,
	  People,
	  jjx,
	 *xsk;
uint   ***state,          /* list haplotypes pairs per genotype             */
       ***fstate,         /* list states of families                        */
	 *po;             /* list of 2 to power                             */

typedef struct {
	   	int fa;
	   	int mo;
		int nchi;
	   	int c[10];
		char name[10];
		char id[12][10];
	  	} ped;

typedef struct {
	        int fhet;
		int mhet;
		int mult[10];
		int chap[10][4][2];
		int fqc;
		int mqc;
		int qc[10];
                } tdt;

ped      *fam;
tdt      *fdat;

typedef struct {
  unsigned          hap[4];
  struct state_rec *next;
} state_rec;


state_rec *new_state_rec(unsigned a, unsigned b,
                         unsigned c, unsigned d, state_rec *n)
{
  state_rec *sr = calloc(1, sizeof(state_rec));

  if (sr) {
    sr->hap[0] = a;
    sr->hap[1] = b;
    sr->hap[2] = c;
    sr->hap[3] = d;
    sr->next = (struct state_rec * ) n;
  }
  return sr;
}


void delete_state_rec(state_rec *sr)
{
  state_rec *tmp = sr;

  if (sr)
    while (sr) {
      sr = (state_rec * ) sr->next;
      free(tmp);
      tmp = sr;
    }
}

void pspace_haptdpnZR(char *geline, int n)
{
        int i;                      /* 01.11.2006 ROHDE  */

        for(i=0;i<n;i++)geline[i] = 32;
        geline[n] = '\0';
}


void swap_i_haptdpnZR(int *u, int *v)
{
  int w = *u;   *u = *v;   *v = w;
}


void swap_ui_haptdpnZR(uint *u, uint *v)
{
  uint w = *u;   *u = *v;   *v = w;
}


void swap_d_haptdpnZR(double *u, double *v)
{
  double w = *u;   *u = *v;   *v = w;
}


bool repeating(int *d, int n)
{
  int i;

  for (i=0; i<n; i++)
    if ( d[i] == d[n] )
      return true;
  return false;
}



int comp(char *a, char *b)
{               /* in order to allow single unknowns(0) in children */
   int i=0;
   while(a[i]){
   if(b[i] != 48 && a[i] != b[i])return 1;
   i++;
   }
   return 0;
}



uint *initPosition(uint *pos, uint n)
{
  uint i;

  for (i=0; i<=n; i++)
    pos[i] = (uint) pow(2.0, ((double) i));
  return pos;
}


int drawHaplo(double *prob, double chance)
{
  int    i;

  for (i=0; chance > prob[i]; i++)
    chance -= prob[i];
  return i;
}


int fill_haptdpnZR(int n)
{
  int i;

  for (i=0; i<nh; i++)
    if ( hc[i] == n )
      return i;
  hc[nh++] = n;
  return nh-1;
}



void recprob_haptdpnZR(int i)
{
  int         j, k, zz, mult;
  double      prod;
  uint        ac[4];

  prob[i] = 0.0;
  zz      = nfstate[i];
  for (j=0; j<zz; j++) {
    mult = 1;
    ac[0] = fstate[i][j][0];
    ac[1] = fstate[i][j][1];
    ac[2] = fstate[i][j][2];
    ac[3] = fstate[i][j][3];

/*  if(ac[0] != ac[1])mult *= 2;
    if(ac[2] != ac[3])mult *= 2;                                      */

    prod     = mult * hap[ac[0]] * hap[ac[1]] * hap[ac[2]] * hap[ac[3]];
    prob[i] += prod;

    if ( ac[0] == ac[1] && ac[2] == ac[3] && ac[0] == ac[2] )
      haptmp[ac[0]] += 4.0 * prod;         
    else if( ac[0] == ac[1] ){
	    haptmp[ac[0]] += 2 * prod;
	    for(k=2;k<4;k++)haptmp[ac[k]] +=prod;
    }
    else if( ac[2] == ac[3] ){
	    haptmp[ac[2]] += 2 * prod;
	    for(k=0;k<2;k++)haptmp[ac[k]] +=prod;
    }
    else for (k=0; k<4; k++) haptmp[ac[k]] += prod; 
  }
}


void selprob_haptdpnZR(int n)
{
  int    i, j,
         maxi,
	 mult,
         k         = 0,
         zz;
  double prodmax   = 0,
         prodfirst = 0,
         prod;

  zz        = nfstate[n];
  mstate[n] = 0;
  for (i=0; i<zz; i++) {  
    mult = 1;
    maxi    = i;
/*  if(fstate[n][i][0] != fstate[n][i][1] ) mult *= 2;
    if(fstate[n][i][2] != fstate[n][i][3] ) mult *= 2;           */
    prodmax = hap[fstate[n][i][0]] * hap[fstate[n][i][1]];
    prodmax *= hap[fstate[n][i][2]] * hap[fstate[n][i][3]];
    prodmax *= mult;

    for (j=i+1; j<zz; j++) {  
      mult = 1;
/*    if(fstate[n][j][0] != fstate[n][j][1] ) mult *= 2;
      if(fstate[n][j][2] != fstate[n][j][3] ) mult *= 2;         */
      prod = hap[fstate[n][j][0]]*hap[fstate[n][j][1]];
      prod *= hap[fstate[n][j][2]]*hap[fstate[n][j][3]];
/*    prod *= mult;                                              */
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
	/*
      if ( fabs(prodfirst - prodmax)/prodfirst > Epsilon )  break;
        */
      if ( fabs(prodfirst - prodmax)/prodfirst > Epsilon )  break;

    swap_ui_haptdpnZR(&fstate[n][maxi][0], &fstate[n][i][0]);
    swap_ui_haptdpnZR(&fstate[n][maxi][1], &fstate[n][i][1]);
    swap_ui_haptdpnZR(&fstate[n][maxi][2], &fstate[n][i][2]);
    swap_ui_haptdpnZR(&fstate[n][maxi][3], &fstate[n][i][3]);
    mstate[n] ++;
  }     /* end i-loop */
  max_prob[n] = prodfirst;
}


double likea_haptdpnZR(void)
{
  int    i, j;
  double f = 0.0;

  for (i=0; i<nf; i++) if(xsk[i]){
    init_dvector(haptmp, 0.0 , nh);
    recprob_haptdpnZR(i);
    if (first == 0)
      for (j=0; j<nh; j++)
        if(prob[i] > 0.0)hapnew[j] += haptmp[j] / prob[i];
    if ( prob[i] <= 0.0 )
      f += 10000.0;
    else
      f -= log(prob[i]);
  }
  return f;
}


void rechap_haptdpnZR(int i, int k1, int k2)
{
  static unsigned int m, ac[2];
  int                 j, k;

  k = k1;
  if (k == 0) {
    ac[0] = ac[1] = 0;
    m     = 0;
  }
  if ( k == k2 ) {
/*  if( ac[0] > ac[1] )return;   */
    state[i][m][0]    = ac[0];
    state[i][m++][1]  = ac[1];

    nstate[i]         = m;
  }
  else {
    if ((geno[i][k] == '1') || (geno[i][k] == '0')) {
      if ( ac[0] & po[k] )         ac[0] -= po[k];
      if ( ac[1] & po[k] )         ac[1] -= po[k];
      rechap_haptdpnZR(i, k+1, k2);
    }
    if ((geno[i][k] == '2') || (geno[i][k] == '0')) {
      if ( !(ac[0] & po[k]) )      ac[0] += po[k];
      if ( !(ac[1] & po[k]) )      ac[1] += po[k];
      rechap_haptdpnZR(i, k+1, k2);
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
        rechap_haptdpnZR(i, k+1, k2);
      }
    }
  }
}


double log_L(int *mg, double *max_prob)
{
  double x = 0.0;
  int    i;

  for (i=0; i<ng; i++)
    x -= ((double) mg[i]) * log(max_prob[i]);
  return x;
}

void printHaplotype_haptdpnZR(int haplo, int len, char *hpres)
{
  int i;
  char hp[len];
	      
  for (i=0; i<len-1; i++)
    sprintf(hp+i,"%c", ((haplo & po[i]) ? '2' : '1'));
    *(hp+len-1)='\0';
    strcat(hpres,hp);
}


void sortByProb_haptdpnZR(double *prob, int *geno, int n)
{
  int    i, j;

  for (i=1; i<n; i++)
    for (j=i; (j > 0) && (prob[j] > prob[j-1]); j--) {
      swap_d_haptdpnZR(&prob[j], &prob[j-1]);
      swap_i_haptdpnZR(&geno[j], &geno[j-1]);
    }
}



void states(int ii, ped* fam, char **gent)
{
  int i, i1, i2, j, k, k1, k2, k3, nfa, nfm, hf[2], hm[2];
  int it1[10], itest, itt, k1e, k2e;
  char gev[len];

  nfa = nstate[fam[ii].fa];
  nfm = nstate[fam[ii].mo];

  for(i=0;i<nfa;i++){
     hf[0] = state[fam[ii].fa][i][0];
     hf[1] = state[fam[ii].fa][i][1];
     if(hf[0] == hf[1])k1e=1;else k1e=2;
     for(j=0;j<nfm;j++){
        hm[0] = state[fam[ii].mo][j][0];
        hm[1] = state[fam[ii].mo][j][1];
        if(hm[0] == hm[1])k2e=1;else k2e=2;
	for(k=0;k<10;k++)it1[k] = 0;
	for(k1=0;k1<k1e;k1++)
	for(k2=0;k2<k2e;k2++){
	for(k3=0;k3<len-1;k3++){
	   i1  =  hf[k1] & po[k3]; 
	   i2  =  hm[k2] & po[k3]; 
	   if( i1 == i2 ) gev[k3] = i1 ? '2' : '1';
	   else gev[k3] = '3';
	   }
	   gev[len-1]='\0';
	   for(itt=0;itt<fam[ii].nchi-1;itt++)
     if( comp(gev,gent[fam[ii].c[itt]]) == 0 ) it1[itt] = 1;
	   for(k=0,itest=0;k<fam[ii].nchi-1;k++)itest += it1[k];

	   if( itest == fam[ii].nchi - 1 )
	      {
	      fstate[ii][nfstate[ii]][0]=hf[0];
	      fstate[ii][nfstate[ii]][1]=hf[1];
	      fstate[ii][nfstate[ii]][2]=hm[0];
	      fstate[ii][nfstate[ii]++][3]=hm[1];
	      for(k=0;k<10;k++)it1[k] = 0;
	      }
	   }
        }
     }
}

void hapcount(int ii)
{
  int i, j, hapi, mult;

  for(i=0;i<nfstate[ii];i++){
          mult = 1;
/*        if(fstate[ii][i][0] != fstate[ii][i][1]) mult *= 2; 
          if(fstate[ii][i][2] != fstate[ii][i][3]) mult *= 2;   */
	  for(j=0;j<4;j++){
		  hapi = fill_haptdpnZR(fstate[ii][i][j]);
		  hap[hapi] += (double)mult/(double)nfstate[ii];
		  fstate[ii][i][j] = hapi;
	  }
  }
}

void printChild(int ii, int jj, int kkk, char **gent, int prin, char *hres)
{
  int i1, i2, k1, k1e, k2, k2e, k3, hf[2], hm[2];
  int kk[4][2], mult;
  char gev[len];

  hf[0] = hc[fstate[ii][kkk][0]];
  hf[1] = hc[fstate[ii][kkk][1]];
  hm[0] = hc[fstate[ii][kkk][2]];
  hm[1] = hc[fstate[ii][kkk][3]];

  if(hf[0] == hf[1])k1e=1;else k1e=2;
  if(hm[0] == hm[1])k2e=1;else k2e=2;

  fdat[ii].fhet = (k1e == 2) ? 1 : 0;
  fdat[ii].mhet = (k2e == 2) ? 1 : 0;

  mult = 0;
  for(k1=0;k1<k1e;k1++)
    for(k2=0;k2<k2e;k2++){
       for(k3=0;k3<len-1;k3++){
	   i1  =  hf[k1] & po[k3]; 
	   i2  =  hm[k2] & po[k3]; 
	   if( i1 == i2 ) gev[k3] = i1 ? '2' : '1';
	   else gev[k3] = '3';
	   }
	   gev[len-1]='\0';
	   if( comp(gev,gent[fam[ii].c[jj]]) == 0 ){
		   kk[mult][0] = k1;
		   kk[mult][1] = k2;
		   mult++;
	   }
    }

           /* ROHDE 14.1.2003 identical heterozygous parents */
           fdat[ii].mult[jj] = mult;
	   /* ROHDE 14.1.2003 identical heterozygous parents */
           if( mult == 2 && (
	   (hf[0] == hm[0] && hf[1] == hm[1]) ||
	   (hf[0] == hm[1] && hf[1] == hm[0])))mult=1;

	   if(mult == 1){
		   k1 = kk[0][0];
		   k2 = kk[0][1];
		   if(prin)printHaplotype_haptdpnZR(hf[k1],len,hres);
		   fdat[ii].chap[jj][0][0] = hf[k1];
		   fdat[ii].chap[jj][2][0] = k1 ? hf[0] : hf[1];
		   if(prin)strcat(hres," <> \0");
		   if(prin)printHaplotype_haptdpnZR(hm[k2],len,hres);
		   fdat[ii].chap[jj][1][0] = hm[k2];
		   fdat[ii].chap[jj][3][0] = k2 ? hm[0] : hm[1];
                   fdat[ii].mult[jj] = mult;
           }
	   else{
		   k1 = kk[0][0];
		   k2 = kk[0][1];
		   if(prin)printHaplotype_haptdpnZR(hf[k1],len,hres);
		   fdat[ii].chap[jj][0][0] = hf[k1];
		   fdat[ii].chap[jj][2][0] = k1 ? hf[0] : hf[1];
		   fdat[ii].chap[jj][0][1] = fdat[ii].chap[jj][2][0];
		   fdat[ii].chap[jj][2][1] = fdat[ii].chap[jj][0][0];
		   if(prin)strcat(hres," <> \0");
		   if(prin)printHaplotype_haptdpnZR(hm[k2],len,hres);
		   fdat[ii].chap[jj][1][0] = hm[k2];
		   fdat[ii].chap[jj][3][0] = k2 ? hm[0] : hm[1];
		   fdat[ii].chap[jj][1][1] = fdat[ii].chap[jj][3][0];
		   fdat[ii].chap[jj][3][1] = fdat[ii].chap[jj][1][0];
		   if(prin){
		   strcat(hres," >> \0");
		   printHaplotype_haptdpnZR(fdat[ii].chap[jj][0][1],len,hres);
		   strcat(hres," <> \0");
		   printHaplotype_haptdpnZR(fdat[ii].chap[jj][1][1],len,hres);
		   }
                   }
}

void tdtpchi(int nfam, int*hcc, int *hl, int nhap, double **om, char **gent, char **dres)
{
  int i, j, ks, ii, jj, h1, h2, mult; 
  char lino[400], lin[400];
  char* CharNull = "\0";
  void prom();

  for(i=0;i<=nhap;i++)
  for(j=i;j<=nhap;j++) om[i][j]=0.0;
  om[nhap][nhap] = 0.0;

  /*
  printf("\n%4s %4s %s >> ",fam[nfam].name,fam[nfam].id[0],geno[fam[nfam].fa]);
  sprintf(lino,"%4s%4s  %4.3f\0",fam[nfam].name,fam[nfam].id[0],fdat[nfam].fqc);
  */
  for(i=0;i<nfstate[nfam];i++){
  h1 = hc[fstate[nfam][i][0]];     /* father */
  h2 = hc[fstate[nfam][i][1]];
  prom(h1,h2,nfam,i,hcc,hl,nhap,1,om);
  }
  
  strcpy(lino,"\0");
  for(i=0;i<nhap;i++)
  for(j=i;j<nhap;j++){
  /* sprintf(lin," %4.3f \0",om[i][j]); 06.11.2014/SKn */
     sprintf(lin," %4.3f %s",om[i][j], CharNull);
  strcat(lino,lin);
  om[i][j] = 0.0;
  }
  /* sprintf(lin," %4.3f \0",om[nhap][nhap]); 06.11.2014/SKn */
     sprintf(lin," %4.3f %s",om[nhap][nhap], CharNull);
  om[nhap][nhap] = 0.0;
  strcat(lino,lin);
  strcpy(dres[jjx],lino);
  jjx++;

  /*
  printf("\n%4s %4s %s >> ",fam[nfam].name,fam[nfam].id[1],geno[fam[nfam].mo]);
  sprintf(lino,"%4s%4s  %4.3f\0",fam[nfam].name,fam[nfam].id[1],fdat[nfam].mqc);
  */
  for(i=0;i<nfstate[nfam];i++){
  h1 = hc[fstate[nfam][i][2]];     /* mother */
  h2 = hc[fstate[nfam][i][3]];
  prom(h1,h2,nfam,i,hcc,hl,nhap,1,om);
  }
  
  strcpy(lino,"\0");
  for(i=0;i<nhap;i++)
  for(j=i;j<nhap;j++){
  /* sprintf(lin," %4.3f \0",om[i][j]); 06.11.2014/SKn */
     sprintf(lin," %4.3f %s",om[i][j], CharNull);
  strcat(lino,lin);
  }
  /* sprintf(lin," %4.3f \0",om[nhap][nhap]); 06.11.2014/SKn */
  sprintf(lin," %4.3f %s",om[nhap][nhap], CharNull);
  strcat(lino,lin);
  strcpy(dres[jjx],lino);
  jjx++;

  for(j=0;j<fam[nfam].nchi-1;j++){ /* children */

  for(ii=0;ii<=nhap;ii++)
  for(jj=ii;jj<=nhap;jj++) om[ii][jj]=0.0;

  /*
  printf("\n%4s %4s %s >> ",
		  fam[nfam].name,fam[nfam].id[j+2],gent[fam[nfam].c[j]]);
  printf("   %4.3f", fdat[nfam].qc[j]);
  sprintf(lino,"%4s%4s  %4.3f\0",fam[nfam].name,fam[nfam].id[j+2],fdat[nfam].qc[j]);
  */
  for(i=0;i<nfstate[nfam];i++){
  printChild(nfam,j,i,gent,0,lino);
  mult = fdat[nfam].mult[j];
  for(ks=0;ks<mult;ks++){
  h1 = fdat[nfam].chap[j][0][ks];
  h2 = fdat[nfam].chap[j][1][ks];
  prom(h1,h2,nfam,i,hcc,hl,nhap,mult,om);
  }
  }
  strcpy(lino,"\0");
  for(ii=0;ii<nhap;ii++)
  for(jj=ii;jj<nhap;jj++){
  /* sprintf(lin," %4.3f \0",om[ii][jj]); 06.11.2014/SKn */
     sprintf(lin," %4.3f %s",om[ii][jj], CharNull);
  strcat(lino,lin);
  }
  /* sprintf(lin," %4.3f \0",om[nhap][nhap]); 06.11.2014/SKn */
  sprintf(lin," %4.3f %s",om[nhap][nhap], CharNull);
  strcat(lino,lin);
  strcpy(dres[jjx],lino);
  jjx++;
  }
}
  

void prom(int h1,int h2,int nfam,int ist,int *hcc,int *hl,int nhap,int mult, double **om)
{
  int i, j, found;
	  
  found = 0;
  for(i=0;i<nhap;i++)
  for(j=i;j<nhap;j++){
    if((hcc[hl[i]] == h1 && hcc[hl[j]] == h2) ||
       (hcc[hl[i]] == h2 && hcc[hl[j]] == h1))
    {
    om[i][j] += pep[nfam][ist]/mult;
    found = 1;
    }
  }
  if(found==0)om[nhap][nhap] += pep[nfam][ist]/mult;

}


double tdtmax_haptdpnZR(RNG rng,int nfam, int nhap, int *hcc, int *hl, int *imax, int length, char **gent, int sim, char **tres)
{
  char      lino[400], lin[400], spa[20];
  char*     CharNull = "\0";
  int 	    i, j, js, k, ks, t = 0 /*SKn*/, n = 0 /*SKn*/, y, ifa, ifb, mult;
  double    s, smax, mt, mut, *sis, **tmat;
  void      printHaplotype_haptdpnZR(), pspace_haptdpnZR();
  void      printChild();

  smax = 0.0;
  jjx  = 0;

  tmat = (double**) calloc(2,sizeof(double*));
  for(i=0;i<2;i++) tmat[i] = (double*)calloc(nhap+1,sizeof(double));
  sis = (double *)calloc(nhap+1,sizeof(double));


  for(i=0;i<=nhap;i++){  /* all essental haplotypes */
	  tmat[0][i]=tmat[1][i]=0.0;
	  sis[i] = 0.0;
  }

  s = 0.0; t = 0; n = 0; mt= 0.0; mut = 0.0;
  for(j=0;j<nfam;j++){             /* all families */
  for(k=0;k<fam[j].nchi - 1; k++){ /*all children */
  if(sim){ y = (Chance(rng) < 0.5)? 1 : 0;
  if(y==0)itp++;  else itn++;}
  for(js=0;js<nfstate[j];js++){    /* all states per family */
  printChild(j,k,js,gent,0,lino);
  mult = fdat[j].mult[k];
  for(ks=0;ks<mult;ks++){
	if( fdat[j].fhet ){
	ifa=0;
        for(i=0;i<nhap;i++){
		if( hcc[hl[i]] == fdat[j].chap[k][0][ks]){
			if(sim && y==1){
			       if(fdat[j].qc[k]== 1)tmat[1][i] +=pep[j][js]/mult;
			       if(fdat[j].qc[k]== 2)tmat[0][i] +=pep[j][js]/mult;
			}
			else{
			       if(fdat[j].qc[k]== 1)tmat[0][i] +=pep[j][js]/mult;
			       if(fdat[j].qc[k]== 2)tmat[1][i] +=pep[j][js]/mult;
			}
			ifa=1;
		}
	}
	ifb=0;
        for(i=0;i<nhap;i++){
		if( hcc[hl[i]] == fdat[j].chap[k][2][ks]){
			if(sim && y==1){
			       if(fdat[j].qc[k]== 1)tmat[0][i] +=pep[j][js]/mult;
			       if(fdat[j].qc[k]== 2)tmat[1][i] +=pep[j][js]/mult;
			}
			else{
			       if(fdat[j].qc[k]== 1)tmat[1][i] +=pep[j][js]/mult;
			       if(fdat[j].qc[k]== 2)tmat[0][i] +=pep[j][js]/mult;
			}
			ifb=1;
		}
	}
	if(ifa==0 && ifb==1){if(sim && y==1){
			       if(fdat[j].qc[k]== 1)tmat[1][nhap] +=pep[j][js]/mult;
			       if(fdat[j].qc[k]== 2)tmat[0][nhap] +=pep[j][js]/mult;
	                     }
			     else{
			       if(fdat[j].qc[k]== 1)tmat[0][nhap] +=pep[j][js]/mult;
			       if(fdat[j].qc[k]== 2)tmat[1][nhap] +=pep[j][js]/mult;
	                     }
	}
	if(ifa==1 && ifb==0){if(sim && y==1){
	                       if(fdat[j].qc[k]== 1)tmat[0][nhap] +=pep[j][js]/mult;
	                       if(fdat[j].qc[k]== 2)tmat[1][nhap] +=pep[j][js]/mult;
	                     }
			     else{
			       if(fdat[j].qc[k]== 1)tmat[1][nhap] +=pep[j][js]/mult;
			       if(fdat[j].qc[k]== 2)tmat[0][nhap] +=pep[j][js]/mult;
                             }

	}
	}

	if( fdat[j].mhet ){
	ifa=0;
        for(i=0;i<nhap;i++){
		if( hcc[hl[i]] == fdat[j].chap[k][1][ks]){
			if(sim && y==1){
			      if(fdat[j].qc[k]== 1)tmat[1][i] +=pep[j][js]/mult;
		              if(fdat[j].qc[k]== 2)tmat[0][i] +=pep[j][js]/mult;
			}
			else{
			      if(fdat[j].qc[k]== 1)tmat[0][i] +=pep[j][js]/mult;
			      if(fdat[j].qc[k]== 2)tmat[1][i] +=pep[j][js]/mult;
			}
			ifa=1;
		}
	}
	ifb=0;
        for(i=0;i<nhap;i++){
		if( hcc[hl[i]] == fdat[j].chap[k][3][ks]){
			if(sim && y==1){
			      if(fdat[j].qc[k]== 1)tmat[0][i] +=pep[j][js]/mult;
			      if(fdat[j].qc[k]== 2)tmat[1][i] +=pep[j][js]/mult;
			}
			else{
			      if(fdat[j].qc[k]== 1)tmat[1][i] +=pep[j][js]/mult;
			      if(fdat[j].qc[k]== 2)tmat[0][i] +=pep[j][js]/mult;
			}
			ifb=1;
		}
	}
	if(ifa==0 && ifb==1){if(sim && y==1){
		              if(fdat[j].qc[k]== 1)tmat[1][nhap] +=pep[j][js]/mult;
		              if(fdat[j].qc[k]== 2)tmat[0][nhap] +=pep[j][js]/mult;
	                     }
			     else{
			      if(fdat[j].qc[k]== 1)tmat[0][nhap] +=pep[j][js]/mult;
			      if(fdat[j].qc[k]== 2)tmat[1][nhap] +=pep[j][js]/mult;
	                     }
	}

	if(ifa==1 && ifb==0){if(sim && y==1){
		              if(fdat[j].qc[k]== 1)tmat[0][nhap] +=pep[j][js]/mult;
		              if(fdat[j].qc[k]== 2)tmat[1][nhap] +=pep[j][js]/mult;
	                     }
			     else{
			      if(fdat[j].qc[k]== 1)tmat[1][nhap] +=pep[j][js]/mult;
			      if(fdat[j].qc[k]== 2)tmat[0][nhap] +=pep[j][js]/mult;
	                     }
	}
	}
  }
  }
  }
  }


        smax = 0;
        for(i=0;i<nhap;i++){   /* not assigned haplotypes i=nhap skipped */
		mut = tmat[0][i] - tmat[1][i];
                mt  = tmat[0][i] + tmat[1][i];
		if(mt)sis[i]  = mut*mut/mt;
		if( sis[i] >= smax){ /* 2006-12-12 ROHDE */
			smax = sis[i];
		        if(!sim)*imax = hl[i];
		}
		s  += sis[i];

        }
        if(!sim){
	  /*
	  for(i=0;i<nhap;i++){
	  printf("\n%4d ",hl[i]);
	  printHaplotype_haptdpnZR(hcc[hl[i]],length);
	  printf(" >> tr:%8.4f  ntr:%8.4f  s1: %8.4f  s: %8.4f\n",
	  tmat[0][i],tmat[1][i],sis[i],s);
          }
	  printf("\n%4d ",hl[nhap]);
          pspace_haptdpnZR(length-1);
	  printf(" >> tr:%8.4f  ntr:%8.4f  s1: %8.4f  s: %8.4f\n",
	  tmat[0][nhap],tmat[1][nhap],sis[nhap],s);
	  */

	  for(i=0;i<nhap;i++){
/*	  Rprintf("\n%4d ",hcc[hl[i]]); *******************************************/
	  /* sprintf(lino,"%4d \0",hcc[hl[i]]); 06.11.2014/SKn */
	     sprintf(lino,"%4d %s",hcc[hl[i]], CharNull);
    printHaplotype_haptdpnZR(hcc[hl[i]],length,lino);
/*	  Rprintf(" >> tr:%8.4f  ntr:%8.4f  s1: %8.4f  s: %8.4f\n",
			            tmat[0][i],tmat[1][i],sis[i],s);   **************************/
	  /* sprintf(lin," >> tr:%8.4f  ntr:%8.4f  s1: %8.4f  s: %8.4f\0",
			            tmat[0][i],tmat[1][i],sis[i],s); 06.11.2014/SKn */
       sprintf(lin," >> tr:%8.4f  ntr:%8.4f  s1: %8.4f  s: %8.4f%s",
			            tmat[0][i],tmat[1][i],sis[i],s, CharNull);
                  
	  strcat(lino,lin);
	  strcpy(tres[jjx++],lino);
          }

	  pspace_haptdpnZR(spa,length);

/*	  Rprintf("\n%s ",spa);        rest of haplotypes collected */
	  strcpy(lino,spa);
	  strcat(lino,"    \0");
/*	  Rprintf(" >> tr:%8.4f  ntr:%8.4f  s1: %8.4f  s: %8.4f\n",
			            tmat[0][nhap],tmat[1][nhap],sis[nhap],s);  ******************/
	  /* sprintf(lin," >> tr:%8.4f  ntr:%8.4f  s1: %8.4f  s: %8.4f\0",
			            tmat[0][nhap],tmat[1][nhap],sis[nhap],s); 06.11.201/SKn */
       sprintf(lin," >> tr:%8.4f  ntr:%8.4f  s1: %8.4f  s: %8.4f%s",
			            tmat[0][nhap],tmat[1][nhap],sis[nhap],s, CharNull);
    strcat(lino,lin);
	  strcpy(tres[jjx++],lino);

	}

	for(i=0;i<2;i++)free( (double*)tmat[i] );
	free( (double**)tmat );
	free( (double*)sis );
  return s;
}

/*
int main(int argc, char *argv[])
{
  char     *filename = InputFile,  * name of the input file              *
            lino[10],
	  **mgent,
	  **mfamid,
	  **mpid,
	  **substr;
  int       inp,
	   *mtrait;
  uint      nsub;
  FILE      *fp;
  void      haptdpn();

  if ( argc > 1 )
  filename = argv[1];
  fileLength (filename, &np, &nf, &len, lino);

  if ( (fp = fopen(filename, "r")) == NULL )
  fatal ("Cannot reopen input file %s", filename);

  mgent     = cmatrix(np, len);
  mfamid    = cmatrix(np,10);
  mpid      = cmatrix(np,10);
  mtrait    = ivector(np);

    inp = 0;
    while((substr=readline(fp,&nsub,0)) != NULL) {
          strcpy(mfamid[inp],substr[0]);
          strcpy(mpid[inp],substr[1]);
          strcpy(mgent[inp],substr[2]);
          sscanf( substr[3], "%d",mtrait+inp);
	  inp++;
	  
    }
    if(inp != np){
	          printf("There is something wrong!\n");
	          exit(1);
    }
    fclose (fp);
    haptdpn(mfamid,mpid,mgent,mtrait,np);
    return 0;
}
*/

void haptdpnZR(char **famid, char **pid, char **gent, int *qtrait, int *xnp,
	     double *lim, char **likeres, char **freqres, char **hapres,
	     char **tdtres, char **pvres)
{
  bool      loop;
 /* char     *filename = InputFile;   name of the input file              */
  char*     CharNull = "\0";
  double    likold,
            pe,
	    smax,
	    /*ssim,*/
	   *p2max;
  int     **pimax,
         npimax=1000 ,
           *hlist,
            i, it, j, k, h, non,
	    drei,
      null,
	    df,
	    combinations = 0 /*SKn*/,
	    imaxhap, /*nmin,*/
	    nz, fz, ex, iqual, nhap, inp;
  uint      iterations;
  ulong     globalIterationCount = 0;
  RNG       rng;

  np = *xnp;

/*  (int)(pow(np,2)+1)
  (int)(pow(np,2)+1)  */
  char lino[ 40000  ], lin[ 40000 ];
  
  nf = 1;
  for(i=0;i<np-1;i++)if(strcmp(famid[i],famid[i+1]))nf++;
  len      = strlen(gent[0])+1;

  mg       = ivector(np);
  merke    = ivector(np);
  nulmer   = ivector(np);
  ge       = ivector(np);
  hlist    = ivector(np);

  po       = uivector(len);

  geno     = cmatrix(np, len);
  max_prob = init_dvector(NULL, 0.0, np);
  prob	   = init_dvector(NULL, 0.0, np);
       
  fam      = (ped *)calloc(nf+1,sizeof(ped));
  fdat     = (tdt *)calloc(nf+1,sizeof(tdt));
    
  hap	   = init_dvector (NULL,  0.0, Hapco);
  hc	   = init_ivector (NULL, -1,Hapco);
	       
  po[0]=1;
  for(i=1;i<len;i++)po[i] = 2*po[i-1];
  combinations = po[len-1];
		   
  /*
  if ( (fp = fopen(filename, "r")) == NULL )
  fatal ("Cannot reopen input file %s", filename);
  */
					     

  InitRNG(-10396, rng);

  iterations = 0;

    /*
    printf("\n Estimation of Multilocus Haplotypes \n");
    printf("*************************************\n");
    printf("\n\nInput file: %s\n\n",filename);
    */
	   
    ng = fz = nf = 0;
    strcpy(fam[0].name,famid[0]);
    
    /* read input data */
    for(inp=0;inp<np;inp++){

    if(inp==0)strcpy(lino,famid[0]);
    if( inp && strcmp(lino,famid[inp]) ){   /* new family */
	nf++; fz = 0;
	strcpy(lino,famid[inp]);
	strcpy(fam[nf].name,lino);
	}

    drei = 0;
    null = 0;
    for (i=0; i<len; i++) {
    if ( gent[inp][i] == '3' )  drei  ++;
    if ( gent[inp][i] == '0' )  null  ++;
    }
    gent[inp][len-1] = '\0';


    if(fz < 2)               /*** parents ***/
    {
    it = 1;
    for (i=0; i<ng; i++) {
    if ( strncmp (geno[i], gent[inp], len) == 0 ) {
									   
    /*
    Rprintf("%2d << %s :::: %3d << %s  %d\n",i,geno[i],inp,gent[inp],ng);
    */

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

    /*
    Rprintf("%2d >> %s :::: %3d >> %s\n",ng,geno[ng],inp,gent[inp]);
    */
    
    ge[inp]    = ng;
    mg[ng]    = 1;
    merke[ng] = drei;
    nulmer[ng] = null;
    ng ++;
    }
    }

    switch (fz){
    case 0   : {fam[nf].fa = ge[inp];
		strcpy(fam[nf].id[0],pid[inp]);
		fdat[nf].fqc = qtrait[inp];
		break;}
    case 1   : {fam[nf].mo = ge[inp];
		strcpy(fam[nf].id[1],pid[inp]);
		fdat[nf].mqc = qtrait[inp];
	        break;}
    default  : {fam[nf].c[fz-2] = inp;
		fam[nf].nchi = fz;
		strcpy(fam[nf].id[fz],pid[inp]);
		fdat[nf].qc[fz-2] = qtrait[inp];
		break;}
    }
    fz++;

    }   /* end while */

    People 	= 2*(nf+1);
    Loci 	= len-1;
    nf += 1;
    /* end of reading sample data */

    nall = 4*nf;
    nstate  = init_ivector (NULL, 0, ng);
    nfstate = init_ivector (NULL, 0, nf);
    mstate  = init_ivector (NULL, 0, nf);
    xsk     = init_ivector (NULL, 1, nf);

    state = (uint***) calloc(ng , sizeof(uint**));                
    fstate = (uint***) calloc(nf, sizeof(uint**));
    for (i=0; i<ng; i++) {
      nz       = po[merke[i]] * po[nulmer[i]] * po[nulmer[i]];
      state[i] = uimatrix(nz, 2);
    }

    
    nh = 0;

    for(i=0;i<ng;i++)rechap_haptdpnZR(i, 0, len-1);

    /* Hier kommt Analyse ueber Familien rein */

    for(i=0;i<nf;i++)
	 fstate[i] = uimatrix(4*nstate[fam[i].fa]*nstate[fam[i].mo],4);

    for(i=0;i<nf;i++)states(i ,fam,gent);

    ex = 0;
    for(i=0;i<nf;i++){
	    if( nfstate[i] == 0 ){
	    xsk[i] = 0;
	    nall  -= 4;  /* skip this family completely */
	    Rprintf("Contradiction in family %s\n",fam[i].name);
	    ex = 1;
	    }
    }
  /*  Rprintf("\n");*/

/*
     if(ex){
     Rprintf("\n");
     error("error in haplotypes");
    }

    destroy_u_array3(state);
*/
    
    for(i=0;i<ng;i++)destroy_u_array2(state[i]);
    free( (uint***)state );
    
    pep = (double**)malloc(nf*sizeof(double*));
    for(i=0;i<nf;i++) pep[i] = init_dvector(NULL,0.0,nfstate[i]);

    for(i=0;i<nf;i++)if(xsk[i]==1)hapcount(i);

    for(i=0;i<nh;i++)hap[i] /= nall;
    first  = 1;           /*** start likelihood outside of annealing loops ***/

    df = nh;

    hapnew = init_dvector(NULL, 0.0, nh );
    haptmp = init_dvector(NULL, 0.0, nh );

    for (i=0; i<nf; i++)if(xsk[i]==1)selprob_haptdpnZR(i);


    likold = likea_haptdpnZR();

    /* Continue computation of mean probabilities */
    

    for(i=0;i<nf;i++)if(xsk[i]){
      for(j=0;j<mstate[i];j++) {
	for(k=0;k<4;k++){
		h = fstate[i][j][k];
		hapnew[h] += 1.0 / mstate[i];
	}
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


    first = 0;

    do {
      loop = 0;
      iterations ++;
      globalIterationCount ++;
      /* Recompute mean probabilities */
      for (i=0; i<nh; i++) {
        if ( fabs(hap[i] - hapnew[i]) > LoopPrecision )
          loop = 1;
        hap[i] = hapnew[i];
      }
      /************ prints best state starting point
      if ( iterations == 1 ) {
        printf("  Iteration %d (Best State) \n", iterations);
        printf("\n");
        for (i=0; i<nh; i++)
          if ( hapnew[i] >= 0.01 ) {
            printf("       hapnew[%8d] = %9.7f    ", hc[i], hapnew[i]);
            printHaplotype_haptdpnZR(hc[i], len);
            printf("\n");
          }
        printf("\n");
	*/

        for (i=0; i<nf; i++)if(xsk[i]==1)selprob_haptdpnZR(i);

        /*	
        for(i=0;i<nf;i++){
        printf("\n%s %s %s >> ",fam[i].name,fam[i].id[0],geno[fam[i].fa]);
        printHaplotype_haptdpnZR(hc[fstate[i][0][0]],len);
        printf(" <> ");
        printHaplotype_haptdpnZR(hc[fstate[i][0][1]],len);
	
        printf("\n%s %s %s >> ",fam[i].name,fam[i].id[1],geno[fam[i].mo]);
        printHaplotype_haptdpnZR(hc[fstate[i][0][2]],len);
        printf(" <> ");
        printHaplotype_haptdpnZR(hc[fstate[i][0][3]],len);

	for(k=0;k<fam[i].nchi-1;k++){
        printf("\n%s %s %s >> ",fam[i].name,fam[i].id[k+2],gent[fam[i].c[k]]);
	printChild(i, k,0);
	}
        }
        printf("\n\n");
      }
      ********** end print best state starting point       *********/

      init_dvector(prob,   0.0, np);
      init_dvector(haptmp, 0.0, nh);
      init_dvector(hapnew, 0.0, nh);
      likold = likea_haptdpnZR();
      for (i=0; i<nh; i++)
        hapnew[i] /= (double) nall;
     
      /*
        printf("       Likelihood = %f\n", likold);
      */

    } while (loop);

    /*
    printf("\n");
    printf("  Results Ensemble means: \n\n");
    */
    nhap = 0;
    jjx = 0;
    for (i=0; i<nh; i++)
      {
      
      if ( hapnew[i] >= *lim ) {
	/* sprintf(lino,"%2d >> \0",i+1); 06.11.2014/SKn */
     sprintf(lino,"%2d >> %s",i+1, CharNull);
        printHaplotype_haptdpnZR(hc[i], len,lino);
        /* sprintf(lin," >> %8d >> %7.4f\0",
	hc[i], hapnew[i]  ); 06.11.2014/SKn */ /* , hap[i]); */
        sprintf(lin," >> %8d >> %7.4f%s",
  hc[i], hapnew[i], CharNull);/* , hap[i]); */
strcat(lino,lin);
	strcpy(freqres[jjx++],lino);
	hlist[nhap++] = i;
      }
      }

      /* sprintf(lino,"Likelihood = %f\0", likold); 06.11.2014/SKn */
      sprintf(lino,"Likelihood = %f%s", likold, CharNull);
      strcpy(likeres[0],lino);


      /* start find best states after MLE   10.8.2000 ROHDE  */
      pimax  = imatrix(nf,npimax);               /* ROHDE  10.8.2000 */
      p2max = init_dvector(NULL,0.0,nf);     /* ROHDE  10.8.2000 */
      for(i=0;i<nf;i++)if(xsk[i]==1)max_prob[i] = 0.0;    /* ROHDE  10.8.2000 */
      for (i=0;i<nf;i++)if(xsk[i]==1){
      for(j=0;j<npimax;j++)pimax[i][j] = -1;
      iqual=1;
      for (j=0;j<nfstate[i];j++){
      pe = hapnew[fstate[i][j][0]] * hapnew[fstate[i][j][1]];
      pe *= hapnew[fstate[i][j][2]] * hapnew[fstate[i][j][3]];
	  if( fstate[i][j][0] != fstate[i][j][1] ) pe += pe;
	  if( fstate[i][j][2] != fstate[i][j][3] ) pe += pe;

	  if(pe > p2max[i]){
	  
      if (pe > max_prob[i]){
	p2max[i] = max_prob[i];
      	max_prob[i] = pe;
      	pimax[i][0]=j;
		for(k=1;k<npimax;k++)pimax[i][k]=-1;
		iqual = 1;                      /* 10.08.2000 ROHDE  */
      	}
	  
	  else{
	  if (pe == max_prob[i] && iqual < 9) pimax[i][iqual++]=j;
	  else p2max[i] = pe;
	  }
      }
      }    /* end of maximum state search */
      }


      /*
      printf("\n Haplotypes after MLE\n");
      */
      jjx=0;
      for(i=0;i<nf;i++)if(xsk[i]==1){
/*    printf("\n%3d  %s >> ",i+1,geno[fam[i].fa]);      */
    /*  sprintf(lino,"%8s %8s %s >> \0",
		      fam[i].name,fam[i].id[0],geno[fam[i].fa]); 06.11.2014/SKn */
        sprintf(lino,"%8s %8s %s >> %s",
  	      fam[i].name,fam[i].id[0],geno[fam[i].fa], CharNull);
          
      printHaplotype_haptdpnZR(hc[fstate[i][pimax[i][0]][0]],len,lino);
      strcat(lino," <> \0");
      printHaplotype_haptdpnZR(hc[fstate[i][pimax[i][0]][1]],len,lino);

      /* sprintf(lin,"  P>> %9.7f D>> %9.7f\0",
	                 max_prob[i],max_prob[i]-p2max[i]); 06.11.2014/SKn */
         sprintf(lin,"  P>> %9.7f D>> %9.7f%s",
                   max_prob[i],max_prob[i]-p2max[i], CharNull);                   
      strcat(lino,lin);
      strcpy(hapres[jjx++],lino);
/*    printf("\n%3d  %s >> ",i+1,geno[fam[i].mo]);      */
      /* sprintf(lino,"%8s %8s %s >> \0",
		      fam[i].name,fam[i].id[1],geno[fam[i].mo]); 06.11.2014/SKn */
      sprintf(lino,"%8s %8s %s >> %s",
  	      fam[i].name,fam[i].id[1],geno[fam[i].mo], CharNull);        
      printHaplotype_haptdpnZR(hc[fstate[i][pimax[i][0]][2]],len,lino);
      strcat(lino," <> \0");
      printHaplotype_haptdpnZR(hc[fstate[i][pimax[i][0]][3]],len,lino);
      strcpy(hapres[jjx++],lino);

      for(k=0;k<fam[i].nchi-1;k++){
/*    printf("\n%3d  %s >> ",i+1,gent[fam[i].c[k]]);     */
      /* sprintf(lino,"%8s %8s %s >> \0",
		      fam[i].name,fam[i].id[k+2],gent[fam[i].c[k]]); */
      sprintf(lino,"%8s %8s %s >> %s",
  	      fam[i].name,fam[i].id[k+2],gent[fam[i].c[k]], CharNull);      
      printChild(i, k,pimax[i][0],gent,1,lino);
      /*
      Rprintf("%s\n",lino);
      */
      strcpy(hapres[jjx++],lino);
      }
		      
      }

      /* endfind best states after MLE   10.3.2000 ROHDE  */

      for(i=0;i<nf;i++)if(xsk[i]){
      for(j=0;j<nfstate[i];j++){
      pe  = hapnew[fstate[i][j][0]]*hapnew[fstate[i][j][1]];
      pe *= hapnew[fstate[i][j][2]]*hapnew[fstate[i][j][3]];
      pep[i][j] = pe;
      }
      }

      for(i=0;i<nf;i++)if(xsk[i]){
      pe = 0.0;
      for(j=0;j<nfstate[i];j++) pe += pep[i][j];
      for(j=0;j<nfstate[i];j++) pep[i][j] /= pe;
      }

      /*
      printf("\n Family State Probabilities:\n\n");
      for(i=0;i<nf;i++){
      printf("%3d ",i+1);
      for(j=0;j<nfstate[i];j++)printf("%5.3lf ",pep[i][j]);
      printf("\n");
      }
      * here start omega calculation */

      /*
      printf("\n Individuals by Haplotype Pairs: \n");
      

      sprintf(lino,"\n fid pid  trait \0");
      for(i=1;i<=nhap;i++)
      for(j=i;j<=nhap;j++){
      sprintf(lin,"%2d,%2d \0",i,j);
      strcat(lino,lin);
      }
      sprintf(lin," X, Y\0");
      strcat(lino,lin);
      strcpy(desres[0],lino);

      om = init_dmatrix(NULL,0.0,nhap+1,nhap+1);
      jjx = 0;
      for(i=0;i<nf;i++)if(xsk[i])tdtpchi(i, hc, hlist, nhap, om, gent, desres);
      */

      /* here starts tdt statistics                       */

/*      Rprintf("\n Weighted TDT-Statistics: \n"); ****************************/

      /*
      nc = 0; qmean = 0.0;
      for(i=0;i<nf;i++)
      for(j=0;j<fam[i].nchi-1;j++){
      nc++;
      qmean += fdat[i].qc[j];
      }

      qmean /= nc;

      for(i=0;i<nf;i++)
      for(j=0;j<fam[i].nchi-1;j++) fdat[i].qc[j] -= qmean;
      */
      
    // smax = tdtmax_haptdpnZR(rng,nf,nhap,hc,hlist,&imaxhap,len,gent,0,tdtres);

      /* start simulations

    nmin = 0; itn=0; itp=0;
    for(i=0;i<10000;i++){
    ssim = tdtmax_haptdpnZR(rng,nf,nhap,hc,hlist,&imaxhap,len,gent,1,tdtres);
    if( ssim < smax ) nmin++;
    }
    

    Rprintf("\n\n TDT - statistics itn: %d itp: %d\n",itn,itp);
    Rprintf("\n S=%10.4e  >> %d  SIG=%10.4e\n",
		    smax, hc[imaxhap], 1.0 - nmin/10000.0);

    sprintf(lino," S=%g  >> %d  SIG=%7.4f\0",
		    smax, hc[imaxhap], 1.0 - nmin/10000.0);
		    
		sprintf(lino,"%10.4e>>%10.4e\0",
		    smax, 1.0 - nmin/10000.0);
    strcpy(pvres[0],lino);  */
    
    // copy from HapRshort.c //
    smax = tdtmax_haptdpnZR(rng,nf,nhap,hc,hlist,&imaxhap,len,gent,0,tdtres);
    /* sprintf(lino,"%g \0",smax); 06.11.2014/SKn*/
    sprintf(lino,"%g %s",smax, CharNull);
    lin[0]='\0';
    printHaplotype_haptdpnZR(hc[imaxhap],len,lin);
    strcat(lino,lin);
    /* sprintf(lin," nhap: %3d\0",nhap); 06.11.2014/SKn */
    sprintf(lin," nhap: %3d%s",nhap, CharNull);
    strcat(lino,lin);
    strcpy(pvres[0],lino);
    // End of copy //

    destroy_i_array(mg);
    destroy_i_array(merke);
    destroy_i_array(nulmer);
    destroy_i_array(hlist);
    destroy_i_array(ge);
    destroy_i_array2(pimax);
    destroy_i_array(mstate);
    destroy_i_array(nstate);
    destroy_i_array(hc);
    destroy_i_array(nfstate);
    destroy_i_array(xsk);
/*  destroy_u_array3(fstate);  */
    for(i=0;i<nf;i++)destroy_u_array2(fstate[i]);
    free( (uint***)fstate);
    destroy_u_array(po);
    destroy_c_array2(geno);
    destroy_d_array(hap);
    destroy_d_array(hapnew);
    destroy_d_array(haptmp);
    destroy_d_array(prob);
    destroy_d_array(p2max);
    destroy_d_array(max_prob);
    for(i=0;i<nf;i++)destroy_d_array(pep[i]);
    free( (double**)pep);
    //destroy_d_array2(pep);
    //destroy_d_array2(om);
    free( (ped*)fam );
    free( (tdt*)fdat );

}
