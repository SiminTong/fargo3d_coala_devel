//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

#include "coala_precision.h"

void init_grid(unsigned int nbins,flt smax, flt smin, flt rhograin, 
                flt unit_m,flt* massgrid,flt* massbins);
void init_grid_test(unsigned int nbins, flt massmax, flt massmin, 
                    flt* massgrid,flt* massbins,flt* massmeanlog);