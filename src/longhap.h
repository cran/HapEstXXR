  /* longhap.h */
  /*************/

typedef enum { BIT, BYTE, WORD } Unit;


void initLonghap(unsigned loci);


unsigned longhapLength(Unit u);


unsigned bit(unsigned n);


unsigned hapcodeLength(unsigned loci);


unsigned * newLonghap();

 
unsigned *clearLonghap(unsigned *h);


unsigned querySNP(unsigned pos, unsigned hap[]);


unsigned * setSNP(unsigned pos, unsigned hap[]);


unsigned * unsetSNP(unsigned pos, unsigned hap[]);


unsigned * copyLonghap(unsigned dest[], unsigned src[]);


unsigned * printLonghap(unsigned hap[], char appended);

