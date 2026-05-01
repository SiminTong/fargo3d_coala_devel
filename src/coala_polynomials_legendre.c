//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

//-------------------------------------------
// Functions for Legendre polynomials
// Equations refere to Lombart & Laibe (2020) (doi:10.1093/mnras/staa3682)
//-------------------------------------------
#include "coala_precision.h"
#include "coala_polynomials_legendre.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief Normalisation coefficient of Legendre polynomials
 *
 * @param   i   order of polynomials
 * @return      normalisation coefficient
 */
flt coeffnorm(unsigned int i){
	return ((flt)2)/(((flt)2)*i+((flt)1));
}

/**
 * @brief Coefficients of Legendre polynomial of order k 
 *
 * @param    k            order of polynomials
 * @return   coeffs_Leg   polynomials coefficients, sorted from low to high order 
 */
void compute_coeff_Leg(unsigned int k,flt* coeff_Leg, unsigned int size_coeff_leg){

	if (k >= size_coeff_leg) {
        printf("Error: k=%u exceeds array size coeff leg %zu\n", k, size_coeff_leg);
        exit(-1);
   }

	if (k==0){
		coeff_Leg[0] = ((flt)1);
	}else if (k==1){
		coeff_Leg[0] = ((flt)0); coeff_Leg[1] = ((flt)1);
	}else if (k==2){
		coeff_Leg[0] = -((flt)5e-1); coeff_Leg[1] = ((flt)0); coeff_Leg[2] = ((flt)15e-1);
	}else if (k==3){
		coeff_Leg[0] = ((flt)0); coeff_Leg[1] = -((flt)15e-1); coeff_Leg[2] = ((flt)0); coeff_Leg[3] = ((flt)25e-1);
	}else if (k==4){
		coeff_Leg[0] = ((flt)3)/((flt)8); coeff_Leg[1] = ((flt)0); coeff_Leg[2] = -((flt)30)/((flt)8); coeff_Leg[3] = ((flt)0); coeff_Leg[4] = ((flt)35)/((flt)8);
	}else if (k==5){
		coeff_Leg[0] = ((flt)0); coeff_Leg[1] = ((flt)15)/((flt)8); coeff_Leg[2] = ((flt)0); coeff_Leg[3] = -((flt)70)/((flt)8); coeff_Leg[4] = ((flt)0); coeff_Leg[5] = ((flt)63)/((flt)8);
	}else if (k==6){
		coeff_Leg[0] = -((flt)5)/((flt)16); coeff_Leg[1] = ((flt)0); coeff_Leg[2] = ((flt)105)/((flt)16); coeff_Leg[3] = ((flt)0); coeff_Leg[4] = -((flt)315)/((flt)16); coeff_Leg[5] = ((flt)0); coeff_Leg[6] = ((flt)231)/((flt)16);
	}else if (k==7){
		coeff_Leg[0] = ((flt)0); coeff_Leg[1] = -((flt)35)/((flt)16); coeff_Leg[2] = ((flt)0); coeff_Leg[3] = ((flt)315)/((flt)16); coeff_Leg[4] = ((flt)0); coeff_Leg[5] = -((flt)693)/((flt)16); coeff_Leg[6] = ((flt)0); coeff_Leg[7] = ((flt)429)/((flt)16);
	}else if (k==8){
		coeff_Leg[0] = ((flt)35)/((flt)128); coeff_Leg[1] = ((flt)0); coeff_Leg[2] = -((flt)1260)/((flt)128); coeff_Leg[3] = ((flt)0); coeff_Leg[4] = ((flt)6930)/((flt)128); coeff_Leg[5] = ((flt)0); coeff_Leg[6] = -((flt)12012)/((flt)128); coeff_Leg[7] = ((flt)0),((flt)6435)/((flt)128);
	}else if (k==9){
		coeff_Leg[0] = ((flt)0); coeff_Leg[1] = ((flt)315)/((flt)128); coeff_Leg[2] = ((flt)0); coeff_Leg[3] = -((flt)4620)/((flt)128); coeff_Leg[4] = ((flt)0); coeff_Leg[5] = ((flt)18018)/((flt)128); coeff_Leg[6] = ((flt)0); coeff_Leg[7] = -((flt)25740)/((flt)128); coeff_Leg[8] = ((flt)0); coeff_Leg[9] = ((flt)12155)/((flt)128);
	}else if (k==10){
		coeff_Leg[0] = -((flt)63)/((flt)256); coeff_Leg[1] = ((flt)0); coeff_Leg[2] = ((flt)3465)/((flt)256); coeff_Leg[3] = ((flt)0); coeff_Leg[4] = -((flt)30030)/((flt)256); coeff_Leg[5] = ((flt)0); coeff_Leg[6] = ((flt)90090)/((flt)256); coeff_Leg[7] = ((flt)0); coeff_Leg[8] = -((flt)109395)/((flt)256); coeff_Leg[9] = ((flt)0); coeff_Leg[10] = ((flt)46189)/((flt)256);
	}else{
		printf("coala_polynomials_legendre.c ->  coeff_Leg -> Wrong order, need kpol <= 10\n");
		exit(-1);
	}
		
}

/**
 * @brief Matrix of coefficients of Legendre polynomial of order up to k 
 *
 * @param   k            order of polynomials
 * @return  mat_coeffs   polynomials coefficients, sorted from low to high order 
 */
void compute_mat_coeffs(unsigned int k,flt* mat_coeffs){

	flt coeffs[k+1];
	memset(coeffs, 0, (k+1)*sizeof(flt));

   for (unsigned int i=0; i<=k; i++){
      unsigned int dim = i+1;
   	compute_coeff_Leg(i,coeffs,dim);
      
      for (unsigned int p=0; p<=i; p++){
         mat_coeffs[p + i*(i+1)] = coeffs[p];
      }

   }

}


/**
 * @brief Evaluate polynomial sum_{i=0}^{k} a_i x^i by Horner's method
 *
 * @param    i    order of polynomials
 * @param    ai   polynomials coefficients, sorted from low to high order
 * @return        evaluation of polynomial of order i at x
 */
flt phi_pol(unsigned int i,flt* ai,flt x){

	flt res = ai[i]; // Start with highest coefficient

   // Horner method need coefficients from high to low ordre => reverse do loop
   for (unsigned int j = i; j > 0; j--){
      res = res * x + ai[j-1];
	   
	}


   return res;
   
}

/**
 * @brief Coefficients of the derivative of Legendre polynomial of order k
 *
 * @param    k             order of polynomials
 * @param    pol_coeffs    polynomials coefficients, sorted from low to high order
 * @return   dpol_coeffs   coefficients of the derivative of polynomials, sorted from low to high order
 */
void polynomial_derivative_coeffs(unsigned int k, flt* pol_coeffs,flt* dpol_coeffs){



	if (k==0){
		dpol_coeffs[0] = ((flt)0);
		return;
	}
   
   for (unsigned int i=1; i<=k; i++){
      dpol_coeffs[i-1] = ((flt)i) * pol_coeffs[i];
   }
   
      
}

/**
 * @brief Derivative of P_k(xij) with respect to x, where xij = 2/hj*(x-xj)
 *
 * @param    k    order of polynomials
 * @param    ak   polynomials coefficients, sorted from low to high order
 * @param    hj   width of bin j
 * @param    xij  variable mapping the mass bin j in [-1,1], needed for Legendre polynomials
 * @return        evaluation at xij of the derivative of P_k(xij) with respect to x
 */
flt dphi_pol_k(unsigned int k,flt* ak,flt hj,flt xij){

	if (k == 0) return (flt)0;

   flt dpol_coeffs[k];  // degree k-1 derivative

   polynomial_derivative_coeffs(k, ak, dpol_coeffs);

   return phi_pol(k - 1, dpol_coeffs, xij) * ((flt)2) / hj;
   
}
	

