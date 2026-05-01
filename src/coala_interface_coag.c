//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

#include "coala_io.h"
#include "coala_L2proj_GQ.h"
#include "coala_limiter.h"
#include "coala_compute_coag.h"
#include "math.h"
#include <stdlib.h>
#include "coala_interface_coag.h"




/** @brief Compute coagulation solver in hydro code for 1 hydro time-step 
 *
 * DG scheme k=0, piecewise constant approximation
 * Function for interface with hydro code
 *
 * @param    nbins         number of dust bins
 * @param    massgrid      grid of masses, borders value of mass bins
 * @param    dim_tabflux   dimension of 1D array tabflux
 * @param    tabflux       1D array to evaluate coagulation flux
 * @param    ind_tabflux   1D array of indices corresponding to values in tabflux
 * @param    rhodust       1D array dust density for each grain size
 * @param    eps_rhodust   minimum value for dust density
 * @param    dv            1D array of the differential velocity between grains in bins lp and l
 * @param    dthydro       hydro timestep
 * @return   new_rhodust   evolved 1D array dust density for each grain size
 */
void coala_coag_k0(unsigned int nbins,flt* massgrid,
					unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,
					flt* rhodust,flt eps_rhodust,flt* dv,flt dthydro,flt* new_rhodust){

	//variable for DG schemels
	flt* flux     = NULL;
    flt* tabdtCFL = NULL;
	flux          = malloc((nbins) * sizeof(flt));
	tabdtCFL      = malloc((nbins) * sizeof(flt));


	//rhodust -> gij
	flt* gij     = NULL;
    flt* gij_new = NULL;
    gij          = malloc((nbins) * sizeof(flt));
    gij_new      = malloc((nbins) * sizeof(flt));

    flt eps_gij =  eps_rhodust/massgrid[nbins];

    for (unsigned int j=0; j<nbins; j++){
    	gij[j] = eps_gij;
    	if (rhodust[j] > eps_rhodust){
    		gij[j] = rhodust[j]/(massgrid[j+1]-massgrid[j]);
    	}
    }


    //CFL coeff for SSPRK-3
    flt coeff_CFL = ((flt)3e-1);

    //to get number of sub-cyling timestep
    int nsub, ndt; 

    //Solve coagulation for 1 dthydro
    compute_coag_k0_kdv(eps_gij,coeff_CFL,nbins,massgrid,gij,
                            dim_tabflux,tabflux,ind_tabflux,
                            flux,tabdtCFL,dthydro,dv,
                            gij_new,&nsub,&ndt);


    //gij -> rhodust
    for (unsigned int j=0; j<nbins; j++){
    	new_rhodust[j] = fmax(gij_new[j]*(massgrid[j+1]-massgrid[j]),eps_rhodust);
    }

    //free variable
    free(flux);
    free(tabdtCFL);
    free(gij);
    free(gij_new);

}


/** @brief Compute coagulation solver in hydro code for 1 hydro time-step 
 *
 * DG scheme k>0, piecewise polynomial approximation
 * Function for interface with hydro code
 *
 * @param    nbins             number of dust bins
 * @param    kpol              degree of polynomials for approximation
 * @param    massgrid          grid of masses, borders value of mass bins
 * @param    massbins          arithmetic mean value of massgrid for each mass bins
 * @param    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
 * @param    Q                 number of points for Gauss-Legendre quadrature
 * @param    vecnodes          nodes of the Legendre polynomials
 * @param    vecweights        weights coefficients for the Gauss-Legendre polynomials
 * @param    dim_tabflux       dimension of 1D array tabflux
 * @param    tabflux           1D array to evaluate coagulation flux
 * @param    ind_tabflux       1D array of indices corresponding to values in tabflux
 * @param    dim_tabintflux    dimension of 1D array tabintflux
 * @param    tabintflux        1D array to evaluate the term including the integral of coagulation flux
 * @param    ind_tabintflux    1D array of indices corresponding to values in tabintflux
 * @param    rhodust           1D array dust density for each grain size
 * @param    eps_rhodust       minimum value for dust density
 * @param    dv                2D array of the differential velocity between grains in bins lp and l
 * @param    dthydro           hydro timestep
 * @return   new_rhodust       evolved 1D array dust density for each grain size
 */
void coala_coag(unsigned int nbins,unsigned int kpol,
				flt* massgrid,flt* massbins,flt* mat_coeffs_leg,
				unsigned int Q,flt* vecnodes,flt* vecweights,
				unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
                unsigned int dim_tabintflux,flt* tabintflux, int* ind_tabintflux,
                flt* rhodust,flt eps_rhodust,flt* dv,flt dthydro,flt* new_rhodust){

	//variable for DG scheme
	flt* flux = NULL;
	flux = malloc((nbins) * sizeof(flt));
	flt* intflux = NULL;
	intflux = malloc((nbins*(kpol+1)) * sizeof(flt));
	flt* tabdtCFL = NULL;
	tabdtCFL = malloc((nbins) * sizeof(flt));
	flt* tabgamma = NULL;
	tabgamma = malloc((nbins) * sizeof(flt));

	//rhodust -> gij + enforcing mass conservation after interpolation
	flt* gij = NULL;
    gij = malloc((nbins) * sizeof(flt));

    flt eps_gij =  eps_rhodust/massgrid[nbins+1];

    L2proj_gij_GQ(nbins,kpol,massgrid,massbins,
    				mat_coeffs_leg,Q,vecnodes,vecweights,
    				eps_rhodust,rhodust,eps_gij,gij);

    //positivity of the polynomials approx
   	//apply scale limiter
   	gammafunction(eps_gij,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,tabgamma);
    for (unsigned int j=0;j<=nbins-1;j++){
        for (unsigned int k=1;k<=kpol;k++){
            gij[k+j*(kpol+1)] = tabgamma[j]*gij[k+j*(kpol+1)];
        }
    }

    //limit values to eps_gij
    check_gij(eps_gij,nbins,kpol,0,gij);



    //CFL coeff for SSPRK-3
    flt* gij_new = NULL;
    gij_new = malloc((nbins) * sizeof(flt));

    
    flt coeff_CFL = ((flt)3e-1);

    //to get number of sub-cyling timestep
    int nsub, ndt; 

    //Solve coagulation for 1 dthydro
    compute_coag_kdv(eps_gij,coeff_CFL,nbins,kpol,
    					massgrid,massbins,mat_coeffs_leg,gij,
						dim_tabflux,tabflux,ind_tabflux,
		                dim_tabintflux,tabintflux,ind_tabintflux,
						dv,flux,intflux,tabgamma,tabdtCFL,dthydro,
						gij_new,&nsub,&ndt);



    //gij -> rhodust
    for (unsigned int j=0; j<nbins; j++){
    	new_rhodust[j] = fmax(gij_new[0+j*(kpol+1)]*(massgrid[j+1]-massgrid[j]),eps_rhodust);
    }

    //free variable
    free(flux);
    free(intflux);
    free(tabdtCFL);
    free(tabgamma);
    free(gij);
    free(gij_new);

}


