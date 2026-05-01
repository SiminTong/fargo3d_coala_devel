//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

#include "coala_precision.h"

void L2proj_GQ_k0(flt eps,unsigned int nbins,flt* massgrid, flt* massbins,unsigned int Q,flt* vecnodes,flt* vecweights, flt* gij);

void L2proj_GQ(flt eps,unsigned int nbins,unsigned int kpol,flt* massgrid, flt* massbins,flt* mat_coeffs_leg,unsigned int Q,flt* vecnodes,flt* vecweights, flt* gij);
void L2proj_gij_GQ(unsigned int nbins,unsigned int kpol,flt* massgrid, flt* massbins,flt* mat_coeffs_leg,unsigned int Q,flt* vecnodes,flt* vecweights,flt eps_rhodust, flt* rhodust, flt eps_gij, flt* gij);
