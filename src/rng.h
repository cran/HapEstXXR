    /* rng.h */
    /*********/


#ifndef __RNG_H__
#define __RNG_H__


typedef long     RNG[3];


int SomeNumber();

void InitRNG(int seed, RNG rng);

double Chance(RNG rng);

double SND(void);
  /* This is "gasdev()" from "Numerical Recipes in C", chapter 7, page 217 */


#endif /* __RNG_H__ */
