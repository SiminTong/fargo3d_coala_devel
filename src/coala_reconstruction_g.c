//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************


#include "coala_precision.h"
#include <math.h>
#include "coala_polynomials_legendre.h"
#include "coala_reconstruction_g.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


//-------------------------------------------
//+
//  Reconstruction of polynomial gh for given x in bin j, Eq.14
//+
//-------------------------------------------


/** @brief Reconstruct the approximation of g using the components gij and the Legendre polynomial basis
 *
 * @param    nbins             number of dust bins
 * @param    kpol              degree of polynomials for approximation
 * @param    massgrid          grid of masses, borders value of mass bins
 * @param    massbins          arithmetic mean value of massgrid for each mass bins
 * @param    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
 * @param    gij               components of g on the polynomial basis
 * @param    j                 indice of bin
 * @param    x                 value to evalute the polynomial
 * @return   res               reconstruction of g in bin j evaluated in x
 */
flt grecons(unsigned int nbins,unsigned int kpol,flt* massgrid,flt* massbins,flt* mat_coeffs_leg,flt* gij,unsigned int j,flt x){

   flt xjgridl = massgrid[j];
   flt xjgridr = massgrid[j+1];
   flt hj = xjgridr-xjgridl;
   flt xj = massbins[j];
   flt xij = ((flt)2)*(x-xj)/hj;
   flt LegP;

   flt ai[kpol+1];
   memset(ai, 0, (kpol + 1) * sizeof(flt));

   flt res = 0;
   for (unsigned int i=0;i<=kpol;i++){
   	// Safe bounds check
		if ( i > kpol) {
		   fprintf(stderr, "Polynomial order exceeds kpol=%u (i=%u)\n",kpol,i);
		}

	   for (unsigned int p=0; p<=i;p++ ){
	      ai[p] =  mat_coeffs_leg[p + i*(i+1)];
	   }    
      LegP = phi_pol(i,ai,xij);
      res += gij[i+j*(kpol+1)]*LegP;
	}
	return res;
}


// //-------------------------------------------
// //+
// // Coefficient of polynomials gh in bin j up to order 3
// //+
// //-------------------------------------------
// void compute_coeff_grecons(unsigned int kflux,unsigned int nbins,flt* massgrid,flt* massbins,flt* gij, flt* coeff_gh){
// 	flt coeff_Leg[4];
// 	flt a00; flt a01; flt a02; flt a03;
// 	flt a10; flt a11; flt a12; flt a13;
// 	flt a20; flt a21; flt a22; flt a23;
// 	flt a30; flt a31; flt a32; flt a33;
// 	flt g0j; flt g1j; flt g2j; flt g3j;

//    for (unsigned int j=0;j<=nbins-1;j++){
//       flt hj = massgrid[j+1]-massgrid[j];
//       flt xj = massbins[j];
  
//       if (kflux==0){
// 			g0j = gij[0+j*(kflux+1)];  g1j = ((flt)0);  g2j = ((flt)0);  g3j = ((flt)0);
	
// 			compute_coeff_Leg(0,coeff_Leg);
// 			a00 = coeff_Leg[0];  a01 = coeff_Leg[1];  a02 = coeff_Leg[2];  a03 = coeff_Leg[3];
// 			a10=((flt)0);  a11=((flt)0);  a12=((flt)0);  a13=((flt)0);
// 			a20=((flt)0);  a21=((flt)0);  a22=((flt)0);  a23=((flt)0);
// 			a30=((flt)0);  a31=((flt)0);  a32=((flt)0);  a33=((flt)0);
// 		}else if (kflux==1){	
// 			g0j = gij[0+j*(kflux+1)]; g1j = gij[1+j*(kflux+1)]; g2j = ((flt)0); g3j = ((flt)0);
			

// 			compute_coeff_Leg(0,coeff_Leg);
// 			a00 = coeff_Leg[0];  a01 = coeff_Leg[1];  a02 = coeff_Leg[2];  a03 = coeff_Leg[3];
// 			compute_coeff_Leg(1,coeff_Leg);
// 			a10 = coeff_Leg[0];  a11 = coeff_Leg[1];  a12 = coeff_Leg[2];  a13 = coeff_Leg[3];
// 			a20=((flt)0);  a21=((flt)0);  a22=((flt)0);  a23=((flt)0);
// 			a30=((flt)0);  a31=((flt)0);  a32=((flt)0);  a33=((flt)0);
// 		}else if (kflux==2){
// 			g0j = gij[0+j*(kflux+1)];  g1j = gij[1+j*(kflux+1)];  g2j = gij[2+j*(kflux+1)];  g3j = ((flt)0);
			
// 			compute_coeff_Leg(0,coeff_Leg);
// 			a00 = coeff_Leg[0];  a01 = coeff_Leg[1];  a02 = coeff_Leg[2];  a03 = coeff_Leg[3];
// 			compute_coeff_Leg(1,coeff_Leg);
// 			a10 = coeff_Leg[0];  a11 = coeff_Leg[1];  a12 = coeff_Leg[2];  a13 = coeff_Leg[3];
// 			compute_coeff_Leg(2,coeff_Leg);
// 			a20 = coeff_Leg[0];  a21 = coeff_Leg[1];  a22 = coeff_Leg[2];  a23 = coeff_Leg[3];
// 			a30=((flt)0);  a31=((flt)0);  a32=((flt)0);  a33=((flt)0);
// 		}else if (kflux==3){
// 			g0j = gij[0+j*(kflux+1)];  g1j = gij[1+j*(kflux+1)];  g2j = gij[2+j*(kflux+1)];  g3j = gij[3+j*(kflux+1)];
			
// 			compute_coeff_Leg(0,coeff_Leg);
// 			a00 = coeff_Leg[0];  a01 = coeff_Leg[1];  a02 = coeff_Leg[2];  a03 = coeff_Leg[3];
// 			compute_coeff_Leg(1,coeff_Leg);
// 			a10 = coeff_Leg[0];  a11 = coeff_Leg[1];  a12 = coeff_Leg[2];  a13 = coeff_Leg[3];
// 			compute_coeff_Leg(2,coeff_Leg);
// 			a20 = coeff_Leg[0];  a21 = coeff_Leg[1];  a22 = coeff_Leg[2];  a23 = coeff_Leg[3];
// 			compute_coeff_Leg(3,coeff_Leg);
// 			a30 = coeff_Leg[0];  a31 = coeff_Leg[1];  a32 = coeff_Leg[2];  a33 = coeff_Leg[3];
// 		}else{
// 			printf("limiter.c -> coeff_grecons -> Wrong order \n");
// 			exit(-1);
// 		}
      
      
    
// 		coeff_gh[0+j*4] = (a00*g0j*pow(hj,3) + a10*g1j*pow(hj,3) + a20*g2j*pow(hj,3) + a30*g3j*pow(hj,3) - ((flt)2)*a01*g0j*pow(hj,2)*xj - ((flt)2)*a11*g1j*pow(hj,2)*xj - ((flt)2)*a21*g2j*pow(hj,2)*xj - ((flt)2)*a31*g3j*pow(hj,2)*xj + ((flt)4)*a02*g0j*hj*pow(xj,2) + ((flt)4)*a12*g1j*hj*pow(xj,2) + ((flt)4)*a22*g2j*hj*pow(xj,2) + ((flt)4)*a32*g3j*hj*pow(xj,2) - ((flt)8)*a03*g0j*pow(xj,3) - ((flt)8)*a13*g1j*pow(xj,3) - ((flt)8)*a23*g2j*pow(xj,3) - ((flt)8)*a33*g3j*pow(xj,3))/pow(hj,3);
    
// 		coeff_gh[1+j*4] = (((flt)2)*a01*g0j*pow(hj,2) + ((flt)2)*a11*g1j*pow(hj,2) + ((flt)2)*a21*g2j*pow(hj,2) + ((flt)2)*a31*g3j*pow(hj,2) - ((flt)8)*a02*g0j*hj*xj - ((flt)8)*a12*g1j*hj*xj - ((flt)8)*a22*g2j*hj*xj - ((flt)8)*a32*g3j*hj*xj + ((flt)24)*a03*g0j*pow(xj,2) + ((flt)24)*a13*g1j*pow(xj,2) + ((flt)24)*a23*g2j*pow(xj,2) + ((flt)24)*a33*g3j*pow(xj,2))/pow(hj,3);
    
// 		coeff_gh[2+j*4] = (((flt)4)*a02*g0j*hj + ((flt)4)*a12*g1j*hj + ((flt)4)*a22*g2j*hj + ((flt)4)*a32*g3j*hj - ((flt)24)*a03*g0j*xj - ((flt)24)*a13*g1j*xj - ((flt)24)*a23*g2j*xj - ((flt)24)*a33*g3j*xj)/pow(hj,3);
      
//       coeff_gh[3+j*4] = (((flt)8)*a03*g0j + ((flt)8)*a13*g1j + ((flt)8)*a23*g2j + ((flt)8)*a33*g3j)/pow(hj,3);
// 	}

// }



