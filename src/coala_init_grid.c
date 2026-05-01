//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

#include "coala_precision.h"
#include <stdio.h>
#include "math.h"


/** @brief Generate massgrid and massbins in logscale from smin and smax
 *
 * @param    nbins        number of dust bins
 * @param    smax         maximum grain size
 * @param    smin         minimum grain size
 * @param    rhograin     intrinsic grain density
 * @param    unit_m       mass unit from hydro code
 * @return   massgrid     grid of masses, borders value of mass bins
 * @return   massbins     arithmetic mean value of massgrid for each mass bins
 */
void init_grid(unsigned int nbins,flt smax, flt smin, flt rhograin, 
                flt unit_m,flt* massgrid,flt* massbins){
    
    flt massmax,massmin,pi;

    pi = ((flt)4)*atan(((flt)1));

    massmax = ((flt)4)*pi*rhograin*pow(smax,3)/((flt)3)/unit_m;
    massmin = ((flt)4)*pi*rhograin*pow(smin,3)/((flt)3)/unit_m;

    flt r = pow(massmax/massmin,((flt)1)/((flt)nbins));
    massgrid[0] = massmin;
    for (unsigned int j=0;j<nbins;j++){
        massgrid[j+1] = r*massgrid[j];

        massbins[j] = (massgrid[j+1]+massgrid[j])/((flt)2);

    }
   
}

/** @brief Generate massgrid and massbins in logscale from massmin and massmax in dimensionless(for tests)
 *
 * @param    nbins        number of dust bins
 * @param    massmax      maximum grain mass
 * @param    masssmin     minimum grain mass
 * @return   massgrid     grid of masses, borders value of mass bins
 * @return   massbins     arithmetic mean value of massgrid for each mass bins
 * @return   massmeanlog  geometric mean value of massgrid for each mass bins
 */
void init_grid_test(unsigned int nbins, flt massmax, flt massmin, 
                    flt* massgrid,flt* massbins,flt* massmeanlog){
    
    flt r = pow(massmax/massmin,((flt)1)/((flt)nbins));
    massgrid[0] = massmin;
    for (unsigned int j=0;j<nbins;j++){
        massgrid[j+1] = r*massgrid[j];

        massbins[j] = (massgrid[j+1]+massgrid[j])/((flt)2);

        massmeanlog[j] =sqrt(massgrid[j+1]*massgrid[j]);

    }
   
}
