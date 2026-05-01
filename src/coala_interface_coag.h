//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

#include "coala_precision.h"

void coala_coag_k0(unsigned int nbins,flt* massgrid,
                    unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,
                    flt* rhodust,flt eps_rhodust,flt* dv,flt dthydro,flt* new_rhodust);
void coala_coag(unsigned int nbins,unsigned int kpol,
                flt* massgrid,flt* massbins,flt* mat_coeffs_leg,
                unsigned int Q,flt* vecnodes,flt* vecweights,
                unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
                unsigned int dim_tabintflux,flt* tabintflux, int* ind_tabintflux,
                flt* rhodust,flt eps_rhodust,flt* dv,flt dthydro,flt* new_rhodust);