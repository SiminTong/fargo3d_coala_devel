//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

#include "coala_precision.h"

flt coagintfluxfunction_GQ(unsigned int kernel,flt K0,unsigned int Q, flt* vecnodes, flt* vecweights,
                           unsigned int nbins,unsigned int kpol,flt* massgrid,flt* mat_coeffs_leg,
                           unsigned int j, unsigned int k,unsigned int lp,unsigned int l,
                           unsigned int ip, unsigned int i);