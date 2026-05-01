//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

#include "coala_precision.h"

void compute_flux_k0(unsigned int nbins,flt* gij,unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,flt* flux);
void compute_flux_k0_kdv(unsigned int nbins,flt* gij,unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,flt* dv,flt* flux);

void compute_flux(unsigned int nbins,unsigned int kpol,flt* gij,
                           unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
                           flt* flux);
void compute_flux_kdv(unsigned int nbins,unsigned int kpol,flt* gij, 
                        unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
                        flt* dv, flt* flux);

void compute_flux_intflux(unsigned int nbins,unsigned int kpol,flt* gij,
                                 unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
                                 unsigned int dim_tabintflux,flt* tabintflux, int* ind_tabintflux,
                                 flt* flux,flt* intflux);
void compute_flux_intflux_kdv(unsigned int nbins,unsigned int kpol,flt* gij,
                                 unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
                                 unsigned int dim_tabintflux,flt* tabintflux, int* ind_tabintflux,
                                 flt* dv, flt* flux,flt* intflux);
                                                                                                   



                                                                                              