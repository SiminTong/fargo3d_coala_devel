//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************
#include "coala_precision.h"

void compute_dim_coagtabflux_GQ_k0(unsigned int kernel,flt K0,unsigned int Q, flt* vecnodes, flt* vecweights,
                                    unsigned int nbins,unsigned int kpol,flt* massgrid,flt* mat_coeffs_leg,
                                    int* dim_tabflux,flt* tabflux_v0, int* ind_tabflux_v0);

void compute_coagtabflux_GQ_k0(unsigned int dim_tabflux,flt* tabflux_v0, int* ind_tabflux_v0, flt* tabflux, int* ind_tabflux);

void compute_dim_coagtabflux_coagtabintflux_GQ(unsigned int kernel,flt K0,unsigned int Q, flt* vecnodes, flt* vecweights,
                                                unsigned int nbins,unsigned int kpol, flt* massgrid, flt* mat_coeffs_leg,
                                                int* dim_tabflux,flt* tabflux_v0,int* ind_tabflux_v0,
                                                int* dim_tabintflux,flt* tabintflux_v0,int* ind_tabintflux_v0);

void compute_coagtabflux_coagtabintflux_GQ(unsigned int dim_tabflux,flt* tabflux_v0,int* ind_tabflux_v0,
                                             unsigned int dim_tabintflux,flt* tabintflux_v0,int* ind_tabintflux_v0,
                                             flt* tabflux,int* ind_tabflux,flt* tabintflux,int* ind_tabintflux);


                                                                 