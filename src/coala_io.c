//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

#include "coala_precision.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h> 
#include "coala_io.h"



/** @brief Check negative values + limit values to eps for gij with k=0
 *
 *
 * @param    eps     minimum value for mass distribution approximation gij
 * @param    nbins   number of dust bins
 * @param    step    step in RK3 from time solver
 * @return   gij     checked and corrected components of g on the polynomial basis
 */
void check_gij_k0(flt eps,unsigned int nbins,unsigned int step,flt* gij){
    for (unsigned int j=0;j<=nbins-1;j++){
        if (gij[j] < ((flt)0)){
            for (unsigned int i=0;i<nbins;i++){
                printf("i=%d, gij[i]=%.15e\n",i,gij[i]);
            }

            printf("step %d, Negative value in j=%d, gij[j]=%.15e \n",step,j,gij[j]);
            exit(-1);

        }

        if (gij[j] < eps){
            gij[j] = eps;
        }

    }
}

/** @brief Check negative values + limit values to eps for gij with k>0
 *
 *
 * @param    eps     minimum value for mass distribution approximation gij
 * @param    nbins   number of dust bins
 * @param    kpol    degree of polynomials for approximation
 * @param    step    step in RK3 from time solver
 * @return   gij     checked and corrected components of g on the polynomial basis
 */

void check_gij(flt eps,unsigned int nbins,unsigned int kpol,unsigned int step,flt* gij){
    for (unsigned int j=0;j<=nbins-1;j++){
        if (gij[0+j*(kpol+1)] < ((flt)0)){
            for (unsigned int i=0;i<nbins;i++){
                printf("i=%d, gij=",i);
                for (unsigned int k=0;k<=kpol;k++){
                    printf("%.15e  ",gij[k+i*(kpol+1)]);

                }
                printf("\n");
            }

            printf("step %d, negative value in j=%d, gij1[0+j*(kpol+1)]=%.15e \n",step,j,gij[0+j*(kpol+1)]);
            exit(-1);

        }

        if (gij[0+j*(kpol+1)] < eps){
            gij[0+j*(kpol+1)] = eps;
            for (unsigned int k=1;k<=kpol;k++){
                gij[k+j*(kpol+1)] = ((flt)0);
            }


        }
        
    }

}
    