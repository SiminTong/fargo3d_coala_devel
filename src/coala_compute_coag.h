//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

#include "coala_precision.h"

void compute_coag_k0(flt eps, flt coeff_CFL,unsigned int nbins,flt* massgrid,flt* gij,
                            unsigned int dim_tabflux,flt* tabflux, int* ind_tabflux,
                            flt* flux, flt* tabdtCFL,flt dthydro,
                            flt* gij_new,int* nsub, int* ndt);

void compute_coag_k0_kdv(flt eps,flt coeff_CFL,unsigned int nbins,flt* massgrid,flt* gij,
                            unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,
                            flt* flux,flt* tabdtCFL,flt dthydro,flt* dv,
                            flt* gij_new,int* nsub, int* ndt);

void compute_coag(flt eps,flt coeff_CFL,unsigned int nbins,unsigned int kpol,
                    flt* massgrid, flt* massbins,flt* mat_coeffs_leg,flt* gij,
                    unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
                    unsigned int dim_tabintflux,flt* tabintflux, int* ind_tabintflux,
                    flt* flux,flt* intflux,flt* tabgamma,flt* tabdtCFL,flt dthydro,
                    flt* gij_new,int* nsub, int* ndt);


void compute_coag_kdv(flt eps, flt coeff_CFL,unsigned int nbins,unsigned int kpol,
                        flt* massgrid, flt* massbins,flt* mat_coeffs_leg,flt* gij,
                        unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
                        unsigned int dim_tabintflux,flt* tabintflux, int* ind_tabintflux,
                        flt* dv,flt* flux,flt* intflux,flt* tabgamma,flt* tabdtCFL,flt dthydro,
                        flt* gij_new,int* nsub, int* ndt);






