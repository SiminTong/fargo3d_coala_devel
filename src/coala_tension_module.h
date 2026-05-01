#include "coala_precision.h"

void TSPSI(unsigned int N, flt* COALA_X, flt* COALA_Y, flt* YP, flt* SIGMA, unsigned int IER, unsigned int SigErr);
void SIGS(unsigned int N,flt* COALA_X,flt* COALA_Y,flt* YP,flt* SIGMA, unsigned int IER,unsigned int SigErr);
void SNHCSH(flt COALA_X, flt* SINHM,flt* COSHM, flt* COSHMM);
void YPC1(unsigned int N, flt* COALA_X, flt* COALA_Y, flt* YP, unsigned int IER);
flt HVAL(flt T, unsigned int N, flt* COALA_X, flt* COALA_Y, flt* YP, flt* SIGMA, unsigned int IER);
flt HPVAL(flt T, unsigned int N, flt* COALA_X, flt* COALA_Y, flt* YP, flt* SIGMA, unsigned int IER);
int INTRVL (flt T,unsigned int N,flt* COALA_X);







