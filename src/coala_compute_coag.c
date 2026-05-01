//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

#include "coala_precision.h"
#include "coala_L2proj_GQ.h"
#include "coala_solver_DG.h"
#include "coala_compute_coag.h"
#include "math.h"



/** @brief Compute coagulation solver (simple kernels) for 1 hydro time-step 
 *
 * DG scheme k=0, piecewise constant approximation
 *
 * @param    eps           minimum value for mass distribution approximation gij
 * @param    coeff_CFL     timestep coefficient for stability of the SSPRK order 3 scheme
 * @param    nbins         number of dust bins
 * @param    massgrid      grid of masses, borders value of mass bins
 * @param    gij           components of g on the polynomial basis
 * @param    dim_tabflux   dimension of 1D array tabflux
 * @param    tabflux       1D array to evaluate coagulation flux
 * @param    ind_tabflux   1D array of indices corresponding to values in tabflux
 * @param    flux          1D array with approximation of the coagulation flux in each bin
 * @param    tabdtCFL      1D array with CFL timestep values in each bin
 * @param    dthydro       hydro timestep
 * @return   gij_new       evolved components of g on the polynomial basis after dthydro
 * @return   nsub          number of subcycling coagulation timestep to reach dthydro
 * @return   ndt           number of hydro timestep, when coagulation CFL > dthydro
 */
void compute_coag_k0(flt eps, flt coeff_CFL,unsigned int nbins,flt* massgrid,flt* gij,
                            unsigned int dim_tabflux,flt* tabflux, int* ind_tabflux,
                            flt* flux, flt* tabdtCFL,flt dthydro,
                            flt* gij_new,int* nsub, int* ndt){

    *nsub = 0;
    *ndt  = 0;

    flt dtCFLsub,dt,dtsub,dtlast;
    //evaluate coagulation CFL
    dtCFLsub = compute_CFL_k0(eps,nbins,massgrid,gij,
                                dim_tabflux,tabflux,ind_tabflux,
                                flux,tabdtCFL);

    dtCFLsub = coeff_CFL*dtCFLsub;
    // printf("dtCFLsub=%.15e \n",dtCFLsub);

    //compare hydro timestep and coagulation CFL
    dt = fmin(dtCFLsub,dthydro);
    // dt = dthydro;

    
    //coagulation subcycling timesteps
    if (dt < dthydro){
        dtsub = ((flt)0);
        while (dtsub<dthydro && dthydro-dtsub>dtCFLsub){
            dtsub += dtCFLsub;

            *nsub += 1;
            
            solver_k0(eps,nbins,massgrid,gij,
                        dim_tabflux,tabflux,ind_tabflux,
                        flux,dtCFLsub,gij_new);
            
            dtCFLsub = compute_CFL_k0(eps,nbins,massgrid,gij_new,
                                        dim_tabflux,tabflux,ind_tabflux,
                                        flux,tabdtCFL);

            dtCFLsub = coeff_CFL*dtCFLsub;
            // printf("dtCFLsub=%.15e \n",dtCFLsub);

            gij = gij_new;

            
                
    
        }
    
        //last timestep to reach dthydro
        dtlast = dthydro-dtsub;
        *nsub += 1;
        solver_k0(eps,nbins,massgrid,gij,dim_tabflux,tabflux,ind_tabflux,flux,dtlast,gij_new);
    

    }else {
        //when coagulation CFL > hydro timstep
        *ndt += 1;

        solver_k0(eps,nbins,massgrid,gij,dim_tabflux,tabflux,ind_tabflux,flux,dt,gij_new);


        
    }


}

/** @brief Compute coagulation solver (physical kernel) for 1 hydro time-step 
 *
 * DG scheme k=0, piecewise constant approximation
 *
 * @param    eps           minimum value for mass distribution approximation gij
 * @param    coeff_CFL     timestep coefficient for stability of the SSPRK order 3 scheme
 * @param    nbins         number of dust bins
 * @param    massgrid      grid of masses, borders value of mass bins
 * @param    gij           components of g on the polynomial basis
 * @param    dim_tabflux   dimension of 1D array tabflux
 * @param    tabflux       1D array to evaluate coagulation flux
 * @param    ind_tabflux   1D array of indices corresponding to values in tabflux
 * @param    flux          1D array with approximation of the coagulation flux in each bin
 * @param    tabdtCFL      1D array with CFL timestep values in each bin
 * @param    dthydro       hydro timestep
 * @param    dv            1D array of the differential velocity between grains in bins lp and l
 * @return   gij_new       evolved components of g on the polynomial basis after dthydro
 * @return   nsub          number of subcycling coagulation timestep to reach dthydro
 * @return   ndt           number of hydro timestep, when coagulation CFL > dthydro
 */
void compute_coag_k0_kdv(flt eps,flt coeff_CFL,unsigned int nbins,flt* massgrid,flt* gij,
                            unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,
                            flt* flux,flt* tabdtCFL,flt dthydro,flt* dv,
                            flt* gij_new,int* nsub, int* ndt){

    *nsub = 0;
    *ndt  = 0;

    flt dtCFLsub,dt,dtsub,dtlast;
    //evaluate coagulation CFL
    dtCFLsub = compute_CFL_k0_kdv(eps,nbins,massgrid,gij,
                                    dim_tabflux,tabflux,ind_tabflux,
                                    dv,flux,tabdtCFL);

    dtCFLsub = coeff_CFL*dtCFLsub;

    //compare hydro timestep and coagulation CFL
    dt = fmin(dtCFLsub,dthydro);
    // printf("dt=%.15e \n",dt,eps);


    //coagulation subcycling timesteps
    if (dt < dthydro){
        dtsub = ((flt)0);
        while (dtsub<dthydro && dthydro-dtsub>dtCFLsub){
            dtsub += dtCFLsub;

            *nsub += 1;
            
            solver_k0_kdv(eps,nbins,massgrid,gij,
                            dim_tabflux,tabflux,ind_tabflux,
                            dv,flux,dtCFLsub,gij_new);
            
            dtCFLsub = compute_CFL_k0_kdv(eps,nbins,massgrid,gij_new,
                                            dim_tabflux,tabflux,ind_tabflux,
                                            dv,flux,tabdtCFL);

            dtCFLsub = coeff_CFL*dtCFLsub;
            // printf("dtCFLsub=%.15e \n",dtCFLsub);
            

            gij = gij_new;
            
            
        
        }
    
        //last timestep to reach dthydro
        dtlast = dthydro-dtsub;
        *nsub += 1;
        solver_k0_kdv(eps,nbins,massgrid,gij,
                        dim_tabflux,tabflux,ind_tabflux,
                        dv,flux,dtlast,gij_new);
    
    }else{
        //when coagulation CFL > hydro timstep
        *ndt += 1;
        solver_k0_kdv(eps,nbins,massgrid,gij,
                        dim_tabflux,tabflux,ind_tabflux,
                        dv,flux,dt,gij_new);

        
    }

}




/** @brief Compute coagulation solver (simple kernels) for 1 hydro time-step 
 *
 * DG scheme k>0, piecewise polynomial approximation
 *
 * @param    eps               minimum value for mass distribution approximation gij
 * @param    coeff_CFL         timestep coefficient for stability of the SSPRK order 3 scheme
 * @param    nbins             number of dust bins
 * @param    kpol              degree of polynomials for approximation
 * @param    massgrid          grid of masses, borders value of mass bins
 * @param    massbins          arithmetic mean value of massgrid for each mass bins
 * @param    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
 * @param    gij               components of g on the polynomial basis
 * @param    dim_tabflux       dimension of 1D array tabflux
 * @param    tabflux           1D array to evaluate coagulation flux
 * @param    ind_tabflux       1D array of indices corresponding to values in tabflux
 * @param    dim_tabintflux    dimension of 1D array tabintflux
 * @param    tabintflux        1D array to evaluate the term including the integral of coagulation flux
 * @param    ind_tabintflux    1D array of indices corresponding to values in tabintflux
 * @param    flux              1D array with approximation of the coagulation flux in each bin
 * @param    intflux           1D array with approximation of the term including intergral of the flux in each bin
 * @param    tabgamma          1D array with limiter coefficient for reconstruction in each bin to ensure positivity
 * @param    tabdtCFL          1D array with CFL timestep values in each bin
 * @param    dthydro           hydro timestep
 * @return   gij_new           evolved components of g on the polynomial basis after dthydro
 * @return   nsub              number of subcycling coagulation timestep to reach dthydro
 * @return   ndt               number of hydro timestep, when coagulation CFL > dthydro
 */
void compute_coag(flt eps,flt coeff_CFL,unsigned int nbins,unsigned int kpol,
                    flt* massgrid, flt* massbins,flt* mat_coeffs_leg,flt* gij,
                    unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
                    unsigned int dim_tabintflux,flt* tabintflux, int* ind_tabintflux,
                    flt* flux,flt* intflux,flt* tabgamma,flt* tabdtCFL,flt dthydro,
                    flt* gij_new,int* nsub, int* ndt){

    *nsub = 0;
    *ndt  = 0;

    flt dtCFLsub,dt,dtsub,dtlast;
    //evaluate coagulation CFL
    dtCFLsub = compute_CFL(eps,nbins,kpol,massgrid,gij,
                            dim_tabflux,tabflux,ind_tabflux,
                            flux,tabdtCFL);

    dtCFLsub = coeff_CFL*dtCFLsub;

    //compare hydro timestep and coagulation CFL
    dt = fmin(dtCFLsub,dthydro);
    // dt = dthydro;

    // printf("dtCFL=%.15e\n",dtCFL);


    //coagulation subcycling timesteps
    if (dt < dthydro){
        dtsub = ((flt)0);
        while (dtsub<dthydro && dthydro-dtsub>dtCFLsub){
            dtsub += dtCFLsub;

            *nsub += 1;

            solver(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,
                        dim_tabflux,tabflux,ind_tabflux,
                        dim_tabintflux,tabintflux,ind_tabintflux,
                        flux,intflux,tabgamma,dtCFLsub,gij_new);
            
            dtCFLsub = compute_CFL(eps,nbins,kpol,massgrid,gij_new,
                                dim_tabflux,tabflux,ind_tabflux,
                                flux,tabdtCFL);

            dtCFLsub = coeff_CFL*dtCFLsub;

            // printf("dtCFLsub=%.15e \n",dtCFLsub);
                
            for (unsigned int i = 0; i < nbins * (kpol + 1); i++) {
                gij[i] = gij_new[i];
            }
    
        }
    
        //last timestep to reach dthydro
        dtlast = dthydro-dtsub;
        *nsub += 1;

        solver(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,
                    dim_tabflux,tabflux,ind_tabflux,
                    dim_tabintflux,tabintflux,ind_tabintflux, 
                    flux,intflux,tabgamma,dtlast,gij_new);
    
    }else {
        //when coagulation CFL > hydro timstep
        *ndt += 1;
        solver(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,
                    dim_tabflux,tabflux,ind_tabflux,
                    dim_tabintflux,tabintflux,ind_tabintflux,
                    flux,intflux,tabgamma,dt,gij_new);
            
    }


}



/** @brief Compute coagulation solver (physical kernel) for 1 hydro time-step 
 *
 * DG scheme k>0, piecewise polynomial approximation
 *
 * @param    eps               minimum value for mass distribution approximation gij
 * @param    coeff_CFL         timestep coefficient for stability of the SSPRK order 3 scheme
 * @param    nbins             number of dust bins
 * @param    kpol              degree of polynomials for approximation
 * @param    massgrid          grid of masses, borders value of mass bins
 * @param    massbins          arithmetic mean value of massgrid for each mass bins
 * @param    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
 * @param    gij               components of g on the polynomial basis
 * @param    dim_tabflux       dimension of 1D array tabflux
 * @param    tabflux           1D array to evaluate coagulation flux
 * @param    ind_tabflux       1D array of indices corresponding to values in tabflux
 * @param    dim_tabintflux    dimension of 1D array tabintflux
 * @param    tabintflux        1D array to evaluate the term including the integral of coagulation flux
 * @param    ind_tabintflux    1D array of indices corresponding to values in tabintflux
 * @param    dv                1D array of the differential velocity between grains in bins lp and l
 * @param    flux              1D array with approximation of the coagulation flux in each bin
 * @param    intflux           1D array with approximation of the term including intergral of the flux in each bin
 * @param    tabgamma          1D array with limiter coefficient for reconstruction in each bin to ensure positivity
 * @param    tabdtCFL          1D array with CFL timestep values in each bin
 * @param    dthydro           hydro timestep
 * @return   gij_new           evolved components of g on the polynomial basis after dthydro
 * @return   nsub              number of subcycling coagulation timestep to reach dthydro
 * @return   ndt               number of hydro timestep, when coagulation CFL > dthydro
 */
void compute_coag_kdv(flt eps, flt coeff_CFL,unsigned int nbins,unsigned int kpol,
                        flt* massgrid, flt* massbins,flt* mat_coeffs_leg,flt* gij,
                        unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
                        unsigned int dim_tabintflux,flt* tabintflux, int* ind_tabintflux,
                        flt* dv,flt* flux,flt* intflux,flt* tabgamma,flt* tabdtCFL,flt dthydro,
                        flt* gij_new,int* nsub, int* ndt){

    *nsub = 0;
    *ndt  = 0;

    flt dtCFLsub,dt,dtsub,dtlast;

    //evaluate coagulation CFL
    dtCFLsub = compute_CFL_kdv(eps,nbins,kpol,massgrid,gij,
                                dim_tabflux,tabflux,ind_tabflux,
                                dv,flux,tabdtCFL);

    dtCFLsub = coeff_CFL*dtCFLsub;

    //compare hydro timestep and coagulation CFL
    dt = fmin(dtCFLsub,dthydro);


    //coagulation subcycling timesteps
    if (dt < dthydro){
        dtsub = ((flt)0);
        while (dtsub<dthydro && dthydro-dtsub>dtCFLsub){
            dtsub += dtCFLsub;

            *nsub += 1;

            solver_kdv(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,
                            dim_tabflux,tabflux,ind_tabflux,
                            dim_tabintflux,tabintflux,ind_tabintflux,
                            dv,flux,intflux,tabgamma,dtCFLsub,gij_new);
            
            dtCFLsub = compute_CFL_kdv(eps,nbins,kpol,massgrid,gij_new,
                                        dim_tabflux,tabflux,ind_tabflux,
                                        dv,flux,tabdtCFL);

            dtCFLsub = coeff_CFL*dtCFLsub;

            // printf("dtCFLsub 2 = %.15e\n",dtCFLsub);

            for (unsigned int i = 0; i < nbins * (kpol + 1); i++) {
                gij[i] = gij_new[i];
            }

        
                
    
        }
        //last timestep to reach dthydro
        dtlast = dthydro-dtsub;
        *nsub += 1;

        solver_kdv(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,
                    dim_tabflux,tabflux,ind_tabflux,
                    dim_tabintflux,tabintflux,ind_tabintflux,
                    dv,flux,intflux,tabgamma,dtlast,gij_new);
    
    }else {
        //when coagulation CFL > hydro timstep
        *ndt += 1;
        solver_kdv(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,
                        dim_tabflux,tabflux,ind_tabflux,
                        dim_tabintflux,tabintflux,ind_tabintflux,
                        dv,flux,intflux,tabgamma,dt,gij_new);
        
    }


}


