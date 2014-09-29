   /* rng.c */
    /*********/



#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "bool.h"
#include "rng.h"
#include "message.h"
#include "R.h"

#define Sqr(x) ((x)*(x))


/* void   time(long int *); */




int SomeNumber(void)
{
  time_t now;

  time(&now);
  return (int) now;
}




void InitRNG(int seed, RNG rng)
{
  static bool used = false;

  /*06.11.2014/SKn
  if (!used) {
    used = true;
    if (seed < 0)
      srand(SomeNumber()); 
      
    else
      srand(seed);
  }*/
  
  /*rng[0] = rand() % 32767;
  rng[1] = rand() % 32767;
  rng[2] = rand() % 32767;*/
  /* Use random number generator of R; 06.11.2014/SKn */
  rng[0] = (int)(unif_rand() * 32767) % 32767;
  rng[1] = (int)(unif_rand() * 32767) % 32767;
  rng[2] = (int)(unif_rand() * 32767) % 32767;
}



double Chance(RNG rng)
{
  double r;

  rng[0] = rng[0] % 177 * 171 - rng[0] / 177 * 2;
  rng[1] = rng[1] % 176 * 172 - rng[1] / 176 * 35;
  rng[2] = rng[2] % 178 * 170 - rng[2] / 178 * 63;
  if (rng[0] < 0)
    rng[0] += 30269;
  if (rng[1] < 0)
    rng[1] += 30307;
  if (rng[2] < 0)
    rng[2] += 30323;
  r = ((double) rng[0]) / 30269.0 +
      ((double) rng[1]) / 30307.0 +
      ((double) rng[2]) / 30323.0;
  return (r - floor(r));
}




double SND(void)
  /* This is "gasdev()" from "Numerical Recipes in C", chapter 7, page 217 */
{
  static int    iset = 0;
  static double gset;
  static RNG    rng;
  static bool   used = false;


  double fac, r, v1, v2;


  if (!used) {
    used = true;
    InitRNG(0, rng);
  }
  if (iset == 0) {
    do {
      v1 = 2.0 * Chance(rng) - 1.0;
      v2 = 2.0 * Chance(rng) - 1.0;
      r  = Sqr(v1) + Sqr(v2);
    } while ((r >= 1.0) || (r == 0.0));
    fac = sqrt(-2.0 * log(r) / r);
    /* Now make the Box-Muller transformation to get two normal deviates.   *
     * Return one and save the other for next time.                         */
    gset = v1 * fac;
    iset = 1;
    return v2 * fac;
  } else {
    iset = 0;
    return gset;
  }
}
