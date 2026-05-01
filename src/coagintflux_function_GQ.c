//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

//--------------------------
// Function to generate tabintflux with GQ method
// ballistic kernel (cross-section) + simple kernels
// Eq.(84) applied on Eq.(56) in Lombart et al., (2024) on arxiv
//--------------------------

#include "coala_polynomials_legendre.h"
#include "coala_functions_integrands_coag.h"
#include "coagintflux_function_GQ.h"
#include <stdlib.h>
#include <stdio.h>
#include "math.h"
#include <string.h>

/**
 * @brief Evaluate the triple integral for coagulation (term in DG scheme) depending only masses with Gauss-Legendre quadrature method.
 *
 * This function is used to calculate the array for the coagulation flux as precomputation.
 *
 * @param    kernel            select the collisional kernel function
 * @param    K0                constant value of the kernel function (used to adapt to code unit)
 * @param    Q                 number of points for Gauss-Legendre quadrature
 * @param    vecnodes          nodes of the Legendre polynomials
 * @param    vecweights        weights coefficients for the Gauss-Legendre polynomials
 * @param    nbins             number of dust bins
 * @param    kpol              degree of polynomials for approximation
 * @param    massgrid          grid of masses, borders value of mass bins
 * @param    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
 * @param    j                 index corresponding to the mass of the new formed grain
 * @param    k                 degree of polynomials in Legendre basis for approximation in bin j
 * @param    lp                index corresponding to the mass of one colliding grain
 * @param    l                 index corresponding to the mass of the second colliding grain
 * @param    ip                degree of polynomials in Legendre basis for approximation in bin lp
 * @param    i                 degree of polynomials in Legendre basis for approximation in bin l
 * @return                     triple integral for the coagulation flux evaluated at j,lp,l,ip,i
 */
flt coagintfluxfunction_GQ(unsigned int kernel,flt K0,unsigned int Q, flt* vecnodes, flt* vecweights,
                           unsigned int nbins,unsigned int kpol,flt* massgrid,flt* mat_coeffs_leg,
                           unsigned int j, unsigned int k,unsigned int lp,unsigned int l,
                           unsigned int ip, unsigned int i){

   flt xlgridl = massgrid[l];
   flt xlgridr = massgrid[l+1];
   flt hl = xlgridr - xlgridl;
   flt xl = ((flt)5e-1)*(massgrid[l+1]+massgrid[l]);

   flt xlpgridl = massgrid[lp];
   flt xlpgridr = massgrid[lp+1];
   flt hlp = xlpgridr - xlpgridl;
   flt xlp = ((flt)5e-1)*(massgrid[lp+1]+massgrid[lp]);

   flt xjgridr = massgrid[j+1];
   flt xjgridl = massgrid[j];
   flt hj = xjgridr - xjgridl;
   flt xj = ((flt)5e-1)*(massgrid[j+1]+massgrid[j]);
   
   flt xmin = massgrid[0];
   flt xmax = massgrid[nbins];


   flt ai[kpol+1];
   flt aip[kpol+1];
   flt ak[kpol+1];

   // Zero them immediately after
   memset(ai, 0, (kpol + 1) * sizeof(flt));
   memset(aip, 0, (kpol + 1) * sizeof(flt));
   memset(ak, 0, (kpol + 1) * sizeof(flt));


   // Safe bounds check
   if (i > kpol || ip > kpol || k > kpol) {
      fprintf(stderr, "Polynomial order exceeds kpol=%u (i=%u, ip=%u, k=%u)\n", 
             kpol, i, ip,k);
      return 0.0;
   }


   if (i > 0){
      for (unsigned int p=0; p<=i;p++ ){
         ai[p] =  mat_coeffs_leg[p + i*(i+1)];
      }
   }else{
      ai[i] = mat_coeffs_leg[0];
   }

   if (ip > 0){
      for (unsigned int p=0; p<=ip;p++ ){
         aip[p] =  mat_coeffs_leg[p + ip*(ip+1)];

      }
   } else{
      aip[ip] = mat_coeffs_leg[0];
   }


   if (k > 0){
      for (unsigned int p=0; p<=k;p++ ){
         ak[p] =  mat_coeffs_leg[p + k*(k+1)];

      }
   } else{
      ak[k] = mat_coeffs_leg[0];
   }


   flt xj_alpha; flt xij;
   flt a_ulp; flt b_ulp; flt ulp_alpha; flt xilp;
   flt a_vl; flt b_vl; flt vl_alpha; flt xil;

   flt res = ((flt) 0 );

   for (unsigned int alpha_x=0; alpha_x < Q; alpha_x++){
      xj_alpha = xj + ((flt)5e-1)*hj*vecnodes[alpha_x];
      xij = vecnodes[alpha_x];

      for (unsigned int alpha_u=0;alpha_u<Q;alpha_u++){
         a_ulp = fmax(xmin, xlpgridl);
         b_ulp = fmin(xj_alpha,xlpgridr);
         ulp_alpha = ((flt)5e-1)*(b_ulp+a_ulp) + ((flt)5e-1)*(b_ulp-a_ulp)*vecnodes[alpha_u];
         xilp = ((flt)2)*(ulp_alpha-xlp)/hlp;

         if (xj_alpha > xlpgridl){

            for (unsigned int alpha_v=0;alpha_v<Q;alpha_v++){
               a_vl = fmax(xj_alpha - ulp_alpha + xmin, xlgridl);
               b_vl = fmin(xmax - ulp_alpha + xmin, xlgridr);
               vl_alpha = ((flt)5e-1)*(b_vl + a_vl) + ((flt)5e-1)*(b_vl - a_vl)*vecnodes[alpha_v];
               xil = ((flt)2)*(vl_alpha-xl)/hl;

               if (xmax - ulp_alpha + xmin > xlgridl && xlgridr > xj_alpha - ulp_alpha + xmin){

                  res += ((flt)125e-3)*hj*(b_ulp-a_ulp)*(b_vl - a_vl)
                           *vecweights[alpha_x]*vecweights[alpha_u]*vecweights[alpha_v]
                           *func_coag_intflux(kernel,K0,k,i,ip,ak,ai,aip,ulp_alpha,vl_alpha,hj,xij,xilp,xil);


               }

            }
         }

      }


   }


   return res;
   

}