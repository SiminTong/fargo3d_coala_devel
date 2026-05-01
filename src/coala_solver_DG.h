//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

#include "coala_precision.h"


flt compute_CFL_k0(flt eps,unsigned int nbins,flt* massgrid,flt* gij,
                    unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,
                    flt* flux,flt* tabdtCFL);
void Lk0_func(unsigned int nbins,flt* massgrid,flt* gij,
                  unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,
                  flt* flux,flt* Lk0);
void solver_k0(flt eps,unsigned int nbins,flt* massgrid,flt* gij,
                  unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,
                  flt* flux, flt dt,flt* gij_new);


flt compute_CFL_k0_kdv(flt eps,unsigned int nbins,flt* massgrid,flt* gij,
                        unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,
                        flt* dv, flt* flux,flt* tabdtCFL);     
void Lk0_func_kdv(unsigned int nbins,flt* massgrid,flt* gij,
                    unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,
                    flt* dv, flt* flux,flt* Lk0);
void solver_k0_kdv(flt eps,unsigned int nbins,flt* massgrid,flt* gij,
                    unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,
                    flt* dv, flt* flux, flt dt,flt* gij_new);


flt compute_CFL(flt eps,unsigned int nbins,unsigned int kpol,flt* massgrid,flt* gij,
                    unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,  
                    flt* flux,flt* tabdtCFL);      
void Lk_func(unsigned int nbins,unsigned int kpol,flt* massgrid,flt* mat_coeffs_leg,flt* gij,
                    unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
                    unsigned int dim_tabintflux,flt* tabintflux,int* ind_tabintflux,
                    flt* flux,flt* intflux,flt* Lk);
void solver(flt eps,unsigned int nbins,unsigned int kpol,
            flt* massgrid, flt* massbins,flt* mat_coeffs_leg,flt* gij,
            unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
            unsigned int dim_tabintflux,flt* tabintflux, int* ind_tabintflux,
            flt* flux,flt* intflux,flt* tabgamma,flt dt,flt* gij_new);




flt compute_CFL_kdv(flt eps,unsigned int nbins,unsigned int kpol,flt* massgrid,flt* gij,
                    unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
                    flt* dv, flt* flux,flt* tabdtCFL);
void Lk_func_kdv(unsigned int nbins,unsigned int kpol,flt* massgrid,flt* mat_coeffs_leg,flt* gij,
                    unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
                    unsigned int dim_tabintflux,flt* tabintflux, int* ind_tabintflux,
                    flt* dv, flt* flux,flt* intflux,flt* Lk);
void solver_kdv(flt eps,unsigned int nbins,unsigned int kpol,flt* massgrid, flt* massbins,flt* mat_coeffs_leg,flt* gij,
                    unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
                    unsigned int dim_tabintflux,flt* tabintflux, int* ind_tabintflux,
                    flt* dv, flt* flux,flt* intflux,flt* tabgamma,flt dt,flt* gij_new);



                                                    

