//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

//-------------------------------------------
// Generate tabflux and tabintflux 
// Equations refere to Lombart & Laibe (2020) (doi:10.1093/mnras/staa3682)
//-------------------------------------------
#include "coala_precision.h"
#include "coagflux_function_GQ.h"
#include "coagintflux_function_GQ.h"
#include <stdio.h>
#include <stdlib.h>
#include "coala_progress_bar.h"
#include "math.h"


/**
 * @brief Precompute array depending only on massgrid to evaluate the coagulation flux
 *        and array dimension to keep only non zero values
 *
 * DG scheme with piecewise constant approximation
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
 * @return   dim_tabflux       number of non zero values in tabflux_v0
 * @return   tabflux_v0        1D array partially filled with non zero values to evaluate coagulation flux 
 * @return   ind_tabflux_v0    1D array with indices corresponding to the non zero values
 */
void compute_dim_coagtabflux_GQ_k0(unsigned int kernel,flt K0,unsigned int Q, flt* vecnodes, flt* vecweights,
                                    unsigned int nbins,unsigned int kpol,flt* massgrid,flt* mat_coeffs_leg,
                                    int* dim_tabflux,flt* tabflux_v0, int* ind_tabflux_v0){
   flt res;
   int u = 0;
   unsigned int iprogress = 0;
   for (unsigned int j=0; j<nbins; j++){
      for (unsigned int lp=0; lp<=j; lp++){
         for (unsigned int l=0; l<nbins; l++){
            res = coagfluxfunction_GQ(kernel,K0,Q,vecnodes,vecweights,
                                       nbins,kpol,massgrid,mat_coeffs_leg,
                                       j,lp,l,0,0);

            if (isnan(res)){
               printf("for k=0, NAN in tabflux GQ\n");
               exit(-1);
            }

            if (res != ((flt)0)){
               tabflux_v0[u] = res;
               ind_tabflux_v0[u*3 + 0] = j;
               ind_tabflux_v0[u*3 + 1] = lp;
               ind_tabflux_v0[u*3 + 2] = l;
               u += 1;
            }

         }
      }

      progressbar(nbins,iprogress);
      iprogress += 1;

   } 
   printf("\n");

   //get dim for tabflux
   *dim_tabflux = u;


}


/**
 * @brief Precompute array with only non zero values to evaluate the coagulation flux
 *
 * DG scheme with piecewise constant approximation
 *
 * @param   dim_tabflux       dimension for tabflux 1D array
 * @param   tabflux_v0        1D array partially filled with non zero values to evaluate coagulation flux 
 * @param   ind_tabflux_v0    1D array with indices corresponding to the non zero values in tabflux_v0
 * @return  tabflux           1D array with only non zero values to evaluate coagulation flux 
 * @return  ind_tabflux       1D array with indices corresponding to values in tabflux
 */
void compute_coagtabflux_GQ_k0(unsigned int dim_tabflux,flt* tabflux_v0, int* ind_tabflux_v0, flt* tabflux, int* ind_tabflux){
   for (unsigned int j=0; j<dim_tabflux; j++){
      tabflux[j] = tabflux_v0[j];
      ind_tabflux[j*3 + 0] = ind_tabflux_v0[j*3 + 0];
      ind_tabflux[j*3 + 1] = ind_tabflux_v0[j*3 + 1];
      ind_tabflux[j*3 + 2] = ind_tabflux_v0[j*3 + 2];
   }
}


/**
 * @brief Precompute arrays depending only on massgrid to evaluate the coagulation flux and term including the integral of the flux
 *        and dimensions of arrays to keep only non zero values
 *
 * DG scheme with piecewise polynomial approximation
 *
 * @param    kernel               select the collisional kernel function
 * @param    K0                   constant value of the kernel function (used to adapt to code unit)
 * @param    Q                    number of points for Gauss-Legendre quadrature
 * @param    vecnodes             nodes of the Legendre polynomials
 * @param    vecweights           weights coefficients for the Gauss-Legendre polynomials
 * @param    nbins                number of dust bins
 * @param    kpol                 degree of polynomials for approximation
 * @param    massgrid             grid of masses, borders value of mass bins
 * @param    mat_coeffs_leg       array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
 * @return   dim_tabflux          number of non zero values in tabflux_v0
 * @return   tabflux_v0           1D array partially filled with non zero values to evaluate coagulation flux 
 * @return   ind_tabflux_v0       1D array with indices corresponding to the non zero values in tabflux_v0
 * @return   dim_tabintflux       number of non zero values in tabintflux_v0
 * @return   tabintflux_v0        1D array partially filled with non zero values to evaluate term including the integral of coagulation flux 
 * @return   ind_tabintflux_v0    1D array with indices corresponding to the non zero values in tabintflux_v0
 */
void compute_dim_coagtabflux_coagtabintflux_GQ(unsigned int kernel,flt K0,unsigned int Q, flt* vecnodes, flt* vecweights,
                                                unsigned int nbins,unsigned int kpol, flt* massgrid, flt* mat_coeffs_leg,
                                                int* dim_tabflux,flt* tabflux_v0,int* ind_tabflux_v0,
                                                int* dim_tabintflux,flt* tabintflux_v0,int* ind_tabintflux_v0){


   flt res;
   int u_flux = 0;
   int u_intflux = 0;

   unsigned int iprogress = 0;
   for (unsigned int j=0; j<nbins; j++){
      for (unsigned int lp=0; lp<=j; lp++){
         for (unsigned int l=0; l<nbins; l++){
            for (unsigned int ip=0; ip<=kpol; ip++){
               for (unsigned int i=0; i<=kpol; i++){
                  //flux
                  res = coagfluxfunction_GQ(kernel,K0,Q,vecnodes,vecweights,
                                             nbins,kpol,massgrid,mat_coeffs_leg,
                                             j,lp,l,ip,i);



                  if (isnan(res)){
                     printf(" k > 0, NAN in tabflux \n");
                     exit(-1);
                  }

                  if (res != ((flt)0)){
                     tabflux_v0[u_flux] = ((flt)res);
                     ind_tabflux_v0[u_flux*5 + 0] = j;
                     ind_tabflux_v0[u_flux*5 + 1] = lp;
                     ind_tabflux_v0[u_flux*5 + 2] = l;
                     ind_tabflux_v0[u_flux*5 + 3] = ip;
                     ind_tabflux_v0[u_flux*5 + 4] = i;
                     u_flux += 1;
                  }

                  //intflux
                  for (unsigned int k=1; k<=kpol; k++){


                     
                     //intflux
                     res = coagintfluxfunction_GQ(kernel,K0,Q,vecnodes,vecweights,
                                                   nbins,kpol,massgrid,mat_coeffs_leg,
                                                   j,k,lp,l,ip,i);

                     if (isnan(res)){
                        printf(" k>0, NAN in tabintflux \n");
                        exit(-1);
                     }

                     if (res != ((flt)0)){
                        tabintflux_v0[u_intflux] = ((flt)res);
                        ind_tabintflux_v0[u_intflux*6 + 0] = j;
                        ind_tabintflux_v0[u_intflux*6 + 1] = k;
                        ind_tabintflux_v0[u_intflux*6 + 2] = lp;
                        ind_tabintflux_v0[u_intflux*6 + 3] = l;
                        ind_tabintflux_v0[u_intflux*6 + 4] = ip;
                        ind_tabintflux_v0[u_intflux*6 + 5] = i;
                        u_intflux += 1;
                     }
                  }
               }
            }
         }
      }

      progressbar(nbins,iprogress);
      iprogress += 1;
   } 
   printf("\n");


   //get dim for tablux and tabintflux
   *dim_tabflux = u_flux;
   *dim_tabintflux = u_intflux;

}



/**
 * @brief Precompute arrays with only non zero values to evaluate the coagulation flux and term including the integral of the flux
 *
 * DG scheme with piecewise polynomial approximation
 *
 * @param   dim_tabflux          dimension for non tabflux 1D array
 * @param   tabflux_v0           1D array partially filled with non zero values to evaluate coagulation flux 
 * @param   ind_tabflux_v0       1D array with indices corresponding to the non zero values in tabflux_v0
 * @param   dim_tabintflux       dimension for tabintflux 1D array
 * @param   tabintflux_v0        1D array partially filled with non zero values to evaluate term including the integral of coagulation flux 
 * @param   ind_tabintflux_v0    1D array with indices corresponding to the non zero values in tabintflux_v0
 * @return  tabflux              1D array filled with non zero values to evaluate coagulation flux 
 * @return  ind_tabflux          1D array with indices corresponding to values in tabflux
 * @return  tabintflux           1D array partially filled with non zero values to evaluate term including the integral of coagulation flux 
 * @return  ind_tabintflux       1D array with indices corresponding to values in tabintflux
 */
void compute_coagtabflux_coagtabintflux_GQ(unsigned int dim_tabflux,flt* tabflux_v0,int* ind_tabflux_v0,
                                             unsigned int dim_tabintflux,flt* tabintflux_v0,int* ind_tabintflux_v0,
                                             flt* tabflux,int* ind_tabflux,flt* tabintflux,int* ind_tabintflux){
   for (unsigned int u=0; u<dim_tabintflux; u++){

      tabintflux[u] = tabintflux_v0[u];
      ind_tabintflux[u*6 + 0] = ind_tabintflux_v0[u*6 + 0];
      ind_tabintflux[u*6 + 1] = ind_tabintflux_v0[u*6 + 1];
      ind_tabintflux[u*6 + 2] = ind_tabintflux_v0[u*6 + 2];
      ind_tabintflux[u*6 + 3] = ind_tabintflux_v0[u*6 + 3];
      ind_tabintflux[u*6 + 4] = ind_tabintflux_v0[u*6 + 4];
      ind_tabintflux[u*6 + 5] = ind_tabintflux_v0[u*6 + 5];

      if (u<dim_tabflux){
         tabflux[u] = tabflux_v0[u];
         ind_tabflux[u*5 + 0] = ind_tabflux_v0[u*5 + 0];
         ind_tabflux[u*5 + 1] = ind_tabflux_v0[u*5 + 1];
         ind_tabflux[u*5 + 2] = ind_tabflux_v0[u*5 + 2];
         ind_tabflux[u*5 + 3] = ind_tabflux_v0[u*5 + 3];
         ind_tabflux[u*5 + 4] = ind_tabflux_v0[u*5 + 4];

      }
   }
}



