//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

//-------------------------------------------
// Function for Limiter for DG scheme (Zhang,Shu 2010)
// To compute gamma coefficient in Eq.28
// Equations refere to Lombart & Laibe (2020) (doi:10.1093/mnras/staa3682)
//-------------------------------------------
#include "coala_precision.h"
#include <math.h>
#include "coala_polynomials_legendre.h"
#include "coala_reconstruction_g.h"
#include "coala_limiter.h"
#include <stdlib.h>
#include <stdio.h>

// //-------------------------------------------
// //+
// //  Minimum value of a general polynomials in interval [xmin,xmax]
// //  which depends only on the coefficient of the polynomials
// //  ONLY for polynomials until order 1
// //  Generated from Mathematica
// //+
// //-------------------------------------------
// flt minvalpolk1(flt a,flt b,flt xmin,flt xmax){

//     if (xmin > ((flt)0) && xmax > xmin && b <= ((flt)0)) {
//         return a+b*xmax;

//     }else if (xmin > ((flt)0) && xmax > xmin && b > ((flt)0)) {
//         return a+b*xmin;

//     }else{
//         return 0;
//     }


// }



/** @brief Minimum value of the approximation of g in each bin
 *
 * @param    nbins                number of dust bins
 * @param    kpol                 degree of polynomials for approximation
 * @param    massgrid             grid of masses, borders value of mass bins
 * @param    massbins             arithmetic mean value of massgrid for each mass bins
 * @param    mat_coeffs_leg       array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
 * @param    gij                  components of g on the polynomial basis
 * @param    j                    indice of bin
 * @param    x                    value to evalute the polynomial
 * @return   tabminval_approx_g   reconstruction of g in bin j evaluated in x
 */
void minval_approx_g(unsigned int nbins,unsigned int kpol,flt* massgrid,flt* massbins,flt* mat_coeffs_leg,flt* gij,flt* tabminval_approx_g){

    if (kpol == 1){
        flt g_left; flt g_right;
        for (unsigned int j=0;j<nbins;j++){
            g_left  = grecons(nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,j,massgrid[j]);
            g_right = grecons(nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,j,massgrid[j+1]);
            tabminval_approx_g[j] = fmin(g_left,g_right);
        }

    }else if (kpol > 1){
        flt xjgridl; flt xjgridr; flt dx; flt xval; flt minval;
        int npoints = 400;
        flt* func_pol = malloc((npoints) * sizeof(flt)); 

        for (unsigned int j=0;j<nbins;j++){
            xjgridl = massgrid[j];
            xjgridr = massgrid[j+1];
            dx = (xjgridr-xjgridl)/(((flt)npoints)-((flt)1));

            func_pol[0] = grecons(nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,j,xjgridl);
            minval = func_pol[0];
            for (unsigned int i=1;i<npoints;i++){
                xval = xjgridl + i*dx;
                func_pol[i] = grecons(nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,j,xval);
                minval = fmin(minval,func_pol[i]);
            }

            tabminval_approx_g[j] = minval;

        }

        free(func_pol);

        // for (unsigned int j = 0; j < nbins; j++) {
        //     flt x_left = massgrid[j], x_right = massgrid[j+1];
        //     flt h = x_right - x_left;
        //     flt tol = 1e-6 * h;
            
        //     // Evaluate endpoints
        //     flt g_left = grecons(nbins, kpol, massgrid, massbins, mat_coeffs_leg, gij, j, x_left);
        //     flt g_right = grecons(nbins, kpol, massgrid, massbins, mat_coeffs_leg, gij, j, x_right);
        //     flt minval = fmin(g_left, g_right);
            
        //     // Adaptive sampling
        //     flt x_mid = 0.5 * (x_left + x_right);
        //     flt g_mid = grecons(nbins, kpol, massgrid, massbins, mat_coeffs_leg, gij, j, x_mid);
        //     minval = fmin(minval, g_mid);
            
        //     // Refine if needed (ternary search like)
        //     for (int iter = 0; iter < 10; iter++) {
        //         flt x1 = x_left + 0.25 * h, x2 = x_left + 0.75 * h;
        //         flt g1 = grecons(nbins, kpol, massgrid, massbins, mat_coeffs_leg, gij, j, x1);
        //         flt g2 = grecons(nbins, kpol, massgrid, massbins, mat_coeffs_leg, gij, j, x2);
        //         minval = fmin(minval, fmin(g1, g2));
                
        //         if (fabs(g1 - g2) < tol) break;  // Converged
        //         h *= 0.5;
        //     }
            
        //     tabminval_approx_g[j] = minval;
        // }


    }else{
        printf("Need correct kpol for minvalpol function in minval_approx_g subroutine \n");
        exit(-1);
    }  

}




/* @brief Limiter coefficient to ensure positivity of the numerical solution (Zhang and Shu 2010)
 *
 * @param    eps               minimum value for mass distribution approximation gij
 * @param    nbins             number of dust bins
 * @param    kpol              degree of polynomials for approximation
 * @param    massgrid          grid of masses, borders value of mass bins
 * @param    massbins          arithmetic mean value of massgrid for each mass bins
 * @param    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
 * @param    gij               components of g on the polynomial basis
 * @return   tabgamma         limiter coefficient in each bin
 */
void gammafunction(flt eps,unsigned int nbins,unsigned int kpol,flt* massgrid,flt* massbins,flt* mat_coeffs_leg,flt* gij,flt* tabgamma){
    flt* tabminval_approx_g = malloc(nbins*sizeof(flt));
    flt meangh;

    //Liu 2019
    if (kpol==0){
        for (unsigned int j=0;j<=nbins-1;j++){
            tabgamma[j] = ((flt)1);
        }
    }else{
        minval_approx_g(nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,tabminval_approx_g);

        
        for (unsigned int j=0;j<=nbins-1;j++){
            if (gij[0+j*(kpol+1)] == tabminval_approx_g[j]){

            }else{
                meangh = gij[0+j*(kpol+1)];

                tabgamma[j] = fmin(((flt)1),fabsl((meangh-eps)/(meangh-tabminval_approx_g[j])));
            }
        }
    }

}
