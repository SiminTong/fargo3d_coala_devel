//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

#include "coala_precision.h"

flt kconst(flt K0,flt u,flt v);
flt kadd(flt K0,flt u,flt v);
flt kmul(flt K0,flt u,flt v);
flt cross_section(flt K0,flt u,flt v);
flt func_kernel(unsigned int kernel, flt K0,flt u, flt v);