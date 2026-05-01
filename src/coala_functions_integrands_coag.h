//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

#include "coala_precision.h"

flt func_coag_flux(unsigned int kernel, flt K0,unsigned int i, unsigned int ip, flt* ai, flt* aip, flt u, flt v, flt xilp, flt xil);
flt func_coag_intflux(unsigned int kernel, flt K0,unsigned int k, unsigned int i, unsigned int ip, flt* ak, flt* ai, flt* aip, flt u, flt v, flt hj, flt xij, flt xilp, flt xil);