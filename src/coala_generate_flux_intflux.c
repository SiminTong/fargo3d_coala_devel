//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

//-------------------------------------------
// Generate flux and intflux 
// Equations refere to Lombart & Laibe (2020) (doi:10.1093/mnras/staa3682)
//-------------------------------------------
#include "coala_precision.h"
#include "coala_generate_flux_intflux.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief compute the approximation of the coagulation flux with DG scheme k=0, for simple kernels
 *
 * Coagulation flux is defined at the right boundary of mass bins, i.e. flux(j) ~ F(m_{j+1/2})
 *
 * @param    nbins          number of dust bins
 * @param    gij            components of g on the polynomial basis
 * @param    dim_tabflux    dimension of 1D array tabflux
 * @param    tabflux        1D array to evaluate coagulation flux
 * @param    ind_tabflux    1D array of indices corresponding to values in tabflux
 * @return   flux           1d array with approximation of the coagulation flux in each bin
 */
void compute_flux_k0(unsigned int nbins,flt* gij,unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,flt* flux){
   
   unsigned int j,lp,l;
   for (unsigned int j=0; j<nbins; j++){
      flux[j] = ((flt)0);
   }

   for (unsigned int u=0; u<dim_tabflux; u++){
      j  = ind_tabflux[u*3 + 0];
      lp = ind_tabflux[u*3 + 1];
      l  = ind_tabflux[u*3 + 2];

      flux[j] += gij[lp]*gij[l]*tabflux[u];

   }

   //mass conservation
   flux[nbins-1] = ((flt)0);

}




/**
 * @brief compute the approximation of the coagulation flux with DG scheme k=0, for simple kernels
 *
 * Coagulation flux is defined at the right boundary of mass bins, i.e. flux(j) ~ F(m_{j+1/2})
 *
 * @param    nbins          number of dust bins
 * @param    gij            components of g on the polynomial basis
 * @param    dim_tabflux    dimension of 1D array tabflux
 * @param    tabflux        1D array to evaluate coagulation flux
 * @param    ind_tabflux    1D array of indices corresponding to value in tabflux
 * @param    dv             1D array of the differential velocity between grains in bins lp and l
 * @return   flux           1D array with approximation of the coagulation flux in each bin
 */
void compute_flux_k0_kdv(unsigned int nbins,flt* gij,unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,flt* dv,flt* flux){
   
   unsigned int j,lp,l,p,pp;
   for (unsigned int j=0; j<nbins; j++){
      flux[j] = ((flt)0);
   }
   for (unsigned int u=0; u<dim_tabflux; u++){
      j  = ind_tabflux[u*3 + 0];
      lp = ind_tabflux[u*3 + 1];
      l  = ind_tabflux[u*3 + 2];

      flux[j] += gij[lp]*gij[l]* tabflux[u]* dv[lp + l*nbins] ;

   }

   //mass conservation
   flux[nbins-1] = ((flt)0);

}


/**
 * @brief compute the approximation of the coagulation flux with DG scheme k > 0, for simple kernels
 *
 * Coagulation flux is defined at the right boundary of mass bins, i.e. flux(j) ~ F(m_{j+1/2})
 *
 * @param    nbins          number of dust bins
 * @param    kpol           degree of polynomials for approximation
 * @param    gij            components of g on the polynomial basis
 * @param    dim_tabflux    dimension of 1D array tabflux
 * @param    tabflux        1D array to evaluate coagulation flux
 * @param    ind_tabflux    1D array of indices corresponding to value in tabflux
 * @return   flux           1D array with approximation of the coagulation flux in each bin
 */
void compute_flux(unsigned int nbins,unsigned int kpol,flt* gij,
                     unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
                     flt* flux){
   
   int j,lp,l,ip,i;

   for (unsigned int j=0; j<nbins; j++){
      flux[j] = ((flt)0);
   }

   for (unsigned int u=0; u<dim_tabflux; u++){
      j  = ind_tabflux[u*5 + 0];
      lp = ind_tabflux[u*5 + 1];
      l  = ind_tabflux[u*5 + 2];
      ip = ind_tabflux[u*5 + 3];
      i  = ind_tabflux[u*5 + 4];


      flux[j] += gij[ip + lp*(kpol+1)] * gij[i + l*(kpol+1)]* tabflux[u];

   }

   //mass conservation
   flux[nbins-1] = ((flt)0);

}



/**
 * @brief compute the approximation of the coagulation flux with DG scheme k > 0, for ballistic kernel
 *
 * Coagulation flux is defined at the right boundary of mass bins, i.e. flux(j) ~ F(m_{j+1/2})
 *
 * @param    nbins          number of dust bins
 * @param    kpol           degree of polynomials for approximation
 * @param    gij            components of g on the polynomial basis
 * @param    dim_tabflux    dimension of 1D array tabflux
 * @param    tabflux        1D array to evaluate coagulation flux
 * @param    ind_tabflux    1D array of indices corresponding to value in tabflux
 * @param    dv             1D array of the differential velocity between grains in bins lp and l
 * @return   flux           1D array with approximation of the coagulation flux in each bin
 */
void compute_flux_kdv(unsigned int nbins,unsigned int kpol,flt* gij, 
                        unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
                        flt* dv, flt* flux){
   
   unsigned int j,lp,l,ip,i;

   for (unsigned int u=0; u<dim_tabflux; u++){
      j  = ind_tabflux[u*5 + 0];
      lp = ind_tabflux[u*5 + 1];
      l  = ind_tabflux[u*5 + 2];
      ip = ind_tabflux[u*5 + 3];
      i  = ind_tabflux[u*5 + 4];

      flux[j] += gij[ip + lp*(kpol+1)] * gij[i + l*(kpol+1)]
                        * tabflux[u]
                        * dv[l + lp*nbins];

   }

   //mass conservation
   flux[nbins-1] = ((flt)0);
}


/**
 * @brief compute the approximation of the coagulation flux and the term including the integral of the flux with DG scheme k > 0, for simple kernels
 *
 * Coagulation flux is defined at the right boundary of mass bins, i.e. flux(j) ~ F(m_{j+1/2})
 *
 * @param    nbins          number of dust bins
 * @param    kpol           degree of polynomials for approximation
 * @param    gij            components of g on the polynomial basis
 * @param    dim_tabflux    dimension of 1D array tabflux
 * @param    tabflux        1D array to evaluate coagulation flux
 * @param    ind_tabflux    1D array of indices corresponding to value in tabflux
 * @param    dim_tabintflux dimension of 1D array tabintflux
 * @param    tabintflux     1D array partially filled with non zero values to evaluate term including the integral of coagulation flux 
 * @param    ind_tabintflux 1D array of indices corresponding to the non zero values
 * @return   flux           1D array with approximation of the coagulation flux in each bin
 * @return   intflux        1D array with approximation of the term including the integral of coagulation flux in each bin
 */
void compute_flux_intflux(unsigned int nbins,unsigned int kpol,flt* gij,
                           unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,
                           unsigned int dim_tabintflux,flt* tabintflux,int* ind_tabintflux,
                           flt* flux,flt* intflux){
   
   unsigned int j_f,lp_f,l_f,ip_f,i_f;

   unsigned int j_intf,k_intf,lp_intf,l_intf,ip_intf,i_intf;

   for (unsigned int j=0; j<nbins; j++){
      flux[j] = ((flt)0);
      for (unsigned int k=0; k<=kpol;k++){
         intflux[k+j*(kpol+1)] = ((flt)0);
      }
   }

   for (unsigned int u=0; u<dim_tabintflux; u++){
      //intflux
      j_intf  = ind_tabintflux[u*6 + 0];
      k_intf  = ind_tabintflux[u*6 + 1];
      lp_intf = ind_tabintflux[u*6 + 2];
      l_intf  = ind_tabintflux[u*6 + 3];
      ip_intf = ind_tabintflux[u*6 + 4];
      i_intf  = ind_tabintflux[u*6 + 5];

      // printf("j_intf=%d \n",j_intf);

      intflux[k_intf + j_intf*(kpol+1)] += gij[ip_intf + lp_intf*(kpol+1)] * gij[i_intf + l_intf*(kpol+1)] * tabintflux[u];
   
      //flux
      if (u<dim_tabflux){
         j_f  = ind_tabflux[u*5 + 0];
         lp_f = ind_tabflux[u*5 + 1];
         l_f  = ind_tabflux[u*5 + 2];
         ip_f = ind_tabflux[u*5 + 3];
         i_f  = ind_tabflux[u*5 + 4];

         flux[j_f] += gij[ip_f + lp_f*(kpol+1)] * gij[i_f + l_f*(kpol+1)] * tabflux[u];
      }
   }

   //mass conservation
   flux[nbins-1] = ((flt)0);
   

}


/**
 * @brief compute the approximation of the coagulation flux and the term including the integral of the flux with DG scheme k > 0, for ballistic kernel
 *
 * Coagulation flux is defined at the right boundary of mass bins, i.e. flux(j) ~ F(m_{j+1/2})
 *
 * @param    nbins          number of dust bins
 * @param    kpol           degree of polynomials for approximation
 * @param    gij            components of g on the polynomial basis
 * @param    dim_tabflux    dimension of 1D array tabflux
 * @param    tabflux        1D array to evaluate coagulation flux
 * @param    ind_tabflux    1D array of indices corresponding to value in tabflux
 * @param    dim_tabintflux dimension of 1D array tabintflux
 * @param    tabintflux     1D array partially filled with non zero values to evaluate term including the integral of coagulation flux 
 * @param    ind_tabintflux 1D array of indices corresponding to the non zero values
 * @param    dv             1D array of the differential velocity between grains in bins lp and l
 * @return   flux           1D array with approximation of the coagulation flux in each bin
 * @return   intflux        1D array with approximation of the term including the integral of coagulation flux in each bin
 */
void compute_flux_intflux_kdv(unsigned int nbins,unsigned int kpol,flt* gij,
                              unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,
                              unsigned int dim_tabintflux,flt* tabintflux,int* ind_tabintflux,
                              flt* dv, flt* flux,flt* intflux){
   
   unsigned int j_f,lp_f,l_f,ip_f,i_f;
   unsigned int j_intf,k_intf,lp_intf,l_intf,ip_intf,i_intf;

   for (unsigned int j=0; j<nbins; j++){
      flux[j] = ((flt)0);
      
      for (unsigned int k=0; k<=kpol;k++){
         intflux[k+j*(kpol+1)] = ((flt)0);
      }
   }

   for (unsigned int u=0; u<dim_tabintflux; u++){
      //intflux
      j_intf  = ind_tabintflux[u*6 + 0];
      k_intf  = ind_tabintflux[u*6 + 1];
      lp_intf = ind_tabintflux[u*6 + 2];
      l_intf  = ind_tabintflux[u*6 + 3];
      ip_intf = ind_tabintflux[u*6 + 4];
      i_intf  = ind_tabintflux[u*6 + 5];


      intflux[k_intf + j_intf*(kpol+1)] += gij[ip_intf + lp_intf*(kpol+1)] * gij[i_intf + l_intf*(kpol+1)]
                                                * tabintflux[u]
                                                * dv[l_intf + lp_intf*nbins];

      
      //flux
      if (u<dim_tabflux){
         j_f  = ind_tabflux[u*5 + 0];
         lp_f = ind_tabflux[u*5 + 1];
         l_f  = ind_tabflux[u*5 + 2];
         ip_f = ind_tabflux[u*5 + 3];
         i_f  = ind_tabflux[u*5 + 4];

         flux[j_f] += gij[ip_f + lp_f*(kpol+1)] * gij[i_f + l_f*(kpol+1)]
                           * tabflux[u]
                           * dv[l_f + lp_f*nbins];
      }
   }

   //mass conservation
   flux[nbins-1] = ((flt)0);

   
   
}





