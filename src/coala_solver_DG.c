//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

//-------------------------------------------
// Functions for time solver DGGQ scheme 
// Equations refere to Lombart & Laibe (2020) (doi:10.1093/mnras/staa3682)
// GQ evaluation from Lombart et al., 2024
//-------------------------------------------
#include "coala_precision.h"
#include <math.h>
#include "coala_polynomials_legendre.h"
#include "coala_generate_flux_intflux.h"
#include "coala_limiter.h"
#include "coala_solver_DG.h"
#include "coala_io.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


/**
 * @brief Coagulation CFL condition for DG scheme k=0 piecewise constant approximation (simple kernels)
 *
 * CFL formulation from Filbet & Laurencot 2004, dt <= mean_g * dm/dF
 *
 * @param    eps            minimum value for mass distribution approximation gij
 * @param    nbins          number of dust bins
 * @param    massgrid       grid of masses, borders value of mass bins
 * @param    gij            components of g on the polynomial basis
 * @param    dim_tabflux    dimension of 1D array tabflux
 * @param    tabflux        1D array to evaluate coagulation flux
 * @param    ind_tabflux    1D array of indices corresponding to values in tabflux
 * @param    flux           1D array with approximation of the coagulation flux in each bin
 * @param    tabdtCFL       1D array with CFL timestep values in each bin
 * @return   dtCFL          CFL timestep restriction as min value of tabdtCFL, for coagulation, DG k=0
 */
flt compute_CFL_k0(flt eps,unsigned int nbins,flt* massgrid,flt* gij,
                    unsigned int dim_tabflux, flt* tabflux,int* ind_tabflux,
                    flt* flux,flt* tabdtCFL){ 

    compute_flux_k0(nbins,gij,dim_tabflux,tabflux,ind_tabflux,flux);
    // for (unsigned int j=0;j<nbins;j++){
    //     printf("j=%d, flux=%.15e\n",j,flux[j]);
    // }
    // exit(-1);

    //bin 0, where no mass entering the mass domain
    tabdtCFL[0] = fabs(gij[0]*(massgrid[1]-massgrid[0])/(flux[0]));

   
    for (unsigned int j=1;j<nbins;j++){
        if (gij[j] > eps){
            tabdtCFL[j] = fabs(gij[j]*(massgrid[j+1]-massgrid[j])/(flux[j]-flux[j-1]));
        }else{
            tabdtCFL[j] = ((flt)1e3);
        }
        
    }

    flt dtCFL = tabdtCFL[0];
    for(unsigned int i = 1; i < nbins; i++){
        dtCFL = fmin(dtCFL,tabdtCFL[i]);
    }

    if (dtCFL == ((flt)0)){
        printf("j,gij,hj,flux,tabdtCFL");
        for (unsigned int j=0; j<nbins; j++){
            printf("%d,%.15e,%.15e,%.15e,%.15e \n ",j,gij[j],massgrid[j+1]-massgrid[j],flux[j],tabdtCFL[j]);
        }
        printf("dtCFL = %.15e\n",dtCFL);
        printf("issue in CFL coagulation k=0 simple kernels \n");
        exit(-1);

    }
   
    return dtCFL;
   
}



/**
 * @brief compute the DG operator L for piecewise constant approximation and simple kernels (see Lombart et al., 2021) 
 *
 * It is used for the time solver
 *
 * @param    nbins          number of dust bins
 * @param    massgrid       grid of masses, borders value of mass bins
 * @param    gij            components of g on the polynomial basis
 * @param    dim_tabflux    dimension of 1D array tabflux
 * @param    tabflux        1D array to evaluate coagulation flux
 * @param    ind_tabflux    1D array of indices corresponding to values in tabflux
 * @param    flux           1D array with approximation of the coagulation flux in each bin
 * @return   Lk0            DG operator for piecewise constant approximation in each bin
 */
void Lk0_func(unsigned int nbins,flt* massgrid,flt* gij,
                unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,
                flt* flux,flt* Lk0){
    compute_flux_k0(nbins,gij,dim_tabflux,tabflux,ind_tabflux,flux);

    // for (unsigned int j=0;j<nbins;j++){
    //     printf("j=%d, flux=%.15e\n",j,flux[j]);
    // }
    // exit(-1);
        
    //bin 0 no incoming mass flux
    Lk0[0] = -flux[0]/(massgrid[1]-massgrid[0]);
    for (unsigned int j=1;j<=nbins-1;j++){
      Lk0[j] = (- (flux[j] - flux[j-1]))/(massgrid[j+1]-massgrid[j]);
    }
    
}


/**
 * @brief Function to compute SSPRK order 3 time solver with piecewise constant approximation for simple kernels
 *
 * See Zhang & Shu 2010 and Lombart et al., 2021
 *
 * @param    eps            minimum value for mass distribution approximation gij
 * @param    nbins          number of dust bins
 * @param    massgrid       grid of masses, borders value of mass bins
 * @param    gij            components of g on the polynomial basis
 * @param    dim_tabflux    dimension of 1D array tabflux
 * @param    tabflux        1D array to evaluate coagulation flux
 * @param    ind_tabflux    1D array of indices corresponding to values in tabflux
 * @param    flux           1D array with approximation of the coagulation flux in each bin
 * @param    dt             timestep
 * @return   gij_new        evolved components of g on the polynomial basis after 1 timestep
 */
void solver_k0(flt eps,unsigned int nbins,flt* massgrid,flt* gij,
                unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,
                flt* flux, flt dt, flt* gij_new){

    flt* gij1 = malloc(nbins * sizeof(flt));
    flt* gij2 = malloc(nbins * sizeof(flt));

    //step 1
    flt* Lk0 = malloc(nbins * sizeof(flt));
    Lk0_func(nbins,massgrid,gij,dim_tabflux,tabflux,ind_tabflux,flux,Lk0);
    for (unsigned int j=0;j<=nbins-1;j++){
        gij1[j] = gij[j] + dt*Lk0[j];
    }

    //check negative values + limit values to eps
    check_gij_k0(eps,nbins,1,gij1);

    
    //step 2
    flt* Lk0_1 = malloc(nbins * sizeof(flt));
    Lk0_func(nbins,massgrid,gij1,dim_tabflux,tabflux,ind_tabflux,flux,Lk0_1);
    for (unsigned int j=0;j<=nbins-1;j++){
        gij2[j] = ((flt)3)*gij[j]/((flt)4) + (gij1[j] + dt*Lk0_1[j])/((flt)4);
    }

    //check negative values + limit values to eps
    check_gij_k0(eps,nbins,2,gij2);
    
    //step 3
    flt* Lk0_2 = malloc(nbins * sizeof(flt));
    Lk0_func(nbins,massgrid,gij2,dim_tabflux,tabflux,ind_tabflux,flux,Lk0_2);
    for (unsigned int j=0;j<=nbins-1;j++){
        gij_new[j] = gij[j]/((flt)3) + ((flt)2)*(gij2[j] + dt*Lk0_2[j])/((flt)3);
    }

    //check negative values + limit values to eps
    check_gij_k0(eps,nbins,3,gij_new);

    //clear variables
    free(gij1);
    free(gij2);
    free(Lk0);
    free(Lk0_1);
    free(Lk0_2);
    gij1  = NULL;
    gij2  = NULL;
    Lk0   = NULL;
    Lk0_1 = NULL;
    Lk0_2 = NULL;
    
}




/**
 * @brief Coagulation CFL condition for DG scheme k=0 piecewise constant approximation (physical kernel)
 *
 * CFL formulation from Filbet & Laurencot 2004, dt <= mean_g * dm/dF
 *
 * @param    eps            minimum value for mass distribution approximation gij
 * @param    nbins          number of dust bins
 * @param    massgrid       grid of masses, borders value of mass bins
 * @param    gij            components of g on the polynomial basis
 * @param    dim_tabflux    dimension of 1D array tabflux
 * @param    tabflux        1D array to evaluate coagulation flux
 * @param    ind_tabflux    1D array of indices corresponding to values in tabflux
 * @param    dv             1D array of the differential velocity between grains in bins lp and l
 * @param    flux           1D array with approximation of the coagulation flux in each bin
 * @param    tabdtCFL       1D array with CFL timestep values in each bin
 * @return   dtCFL          CFL timestep restriction as min value of tabdtCFL, for coagulation, DG k=0
 */
flt compute_CFL_k0_kdv(flt eps,unsigned int nbins,flt* massgrid,flt* gij,
                        unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,
                        flt* dv, flt* flux,flt* tabdtCFL){  

    compute_flux_k0_kdv(nbins,gij,dim_tabflux,tabflux,ind_tabflux,dv,flux);
    // for (unsigned int j=0;j<nbins;j++){
    //     printf("j=%d, flux=%.15e\n",j,flux[j]);
    // }
    // exit(-1);

    //bin 0, where no mass entering the mass domain
    tabdtCFL[0] = fabs(gij[0]*(massgrid[1]-massgrid[0])/(flux[0]));
    


    for (unsigned int j=1;j<nbins;j++){
        if (gij[j] > eps){
            tabdtCFL[j] = fabs(gij[j]*(massgrid[j+1]-massgrid[j])/(flux[j]-flux[j-1]));
        }else{
            tabdtCFL[j] = ((flt)1e3);
        }
        
    }

    flt dtCFL = tabdtCFL[0];
    for(unsigned int i = 1; i < nbins; i++){
        dtCFL = fmin(dtCFL,tabdtCFL[i]);
    }

    if (dtCFL == ((flt)0)){
        printf("j,gij,hj,flux,tabdtCFL");
        for (unsigned int j=0; j<nbins; j++){
            printf("%d,%.15e,%.15e,%.15e,%.15e \n ",j,gij[j],massgrid[j+1]-massgrid[j],flux[j],tabdtCFL[j]);
        }
        printf("dtCFL = %.15e\n",dtCFL);
        printf("issue in CFL coagulation k=0 physical kernel \n");
        exit(-1);

    }
   
    return dtCFL;
   
}


/**
 * @brief compute the DG operator L for piecewise constant approximation and physical kernel (see Lombart et al., 2021) 
 *
 * It is used for the time solver
 *
 * @param    nbins          number of dust bins
 * @param    massgrid       grid of masses, borders value of mass bins
 * @param    gij            components of g on the polynomial basis
 * @param    dim_tabflux    dimension of 1D array tabflux
 * @param    tabflux        1D array to evaluate coagulation flux
 * @param    ind_tabflux    1D array of indices corresponding to values in tabflux
 * @param    dv             1D array of the differential velocity between grains in bins lp and l
 * @param    flux           1D array with approximation of the coagulation flux in each bin
 * @return   Lk0            DG operator for piecewise constant approximation in each bin
 */
void Lk0_func_kdv(unsigned int nbins,flt* massgrid,flt* gij,
                    unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,
                    flt* dv, flt* flux,flt* Lk0){

    compute_flux_k0_kdv(nbins,gij,dim_tabflux,tabflux,ind_tabflux,dv,flux);
        
    //bin 0 no incoming mass flux
    Lk0[0] = -flux[0]/(massgrid[1]-massgrid[0]);
    for (unsigned int j=1;j<=nbins-1;j++){
      Lk0[j] = (- (flux[j] - flux[j-1]))/(massgrid[j+1]-massgrid[j]);
    }
    
}


/**
 * @brief Function to compute SSPRK order 3 time solver with piecewise constant approximation for physical kernel
 *
 * See Zhang & Shu 2010 and Lombart et al., 2021
 *
 * @param    eps            minimum value for mass distribution approximation gij
 * @param    nbins          number of dust bins
 * @param    massgrid       grid of masses, borders value of mass bins
 * @param    gij            components of g on the polynomial basis
 * @param    dim_tabflux    dimension of 1D array tabflux
 * @param    tabflux        1D array to evaluate coagulation flux
 * @param    ind_tabflux    1D array of indices corresponding to values in tabflux
 * @param    dv             1D array of the differential velocity between grains in bins lp and l
 * @param    flux           1D array with approximation of the coagulation flux in each bin
 * @param    dt             timestep
 * @return   gij_new        evolved components of g on the polynomial basis after 1 timestep
 */
void solver_k0_kdv(flt eps,unsigned int nbins,flt* massgrid,flt* gij,
                    unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux,
                    flt* dv, flt* flux, flt dt, flt* gij_new){

    flt* gij1 = malloc(nbins * sizeof(flt));
    flt* gij2 = malloc(nbins * sizeof(flt));

    //step 1
    flt* Lk0 = malloc(nbins * sizeof(flt));
    Lk0_func_kdv(nbins,massgrid,gij,dim_tabflux,tabflux,ind_tabflux,dv,flux,Lk0);
    for (unsigned int j=0;j<=nbins-1;j++){
        gij1[j] = gij[j] + dt*Lk0[j];
    }

    //check negative values + limit values to eps
    check_gij_k0(eps,nbins,1,gij1);
    
    //step 2
    flt* Lk0_1 = malloc(nbins * sizeof(flt));
    Lk0_func_kdv(nbins,massgrid,gij1,dim_tabflux,tabflux,ind_tabflux,dv,flux,Lk0_1);
    for (unsigned int j=0;j<=nbins-1;j++){
        gij2[j] = ((flt)3)*gij[j]/((flt)4) + (gij1[j] + dt*Lk0_1[j])/((flt)4);
    }

    //check negative values + limit values to eps
    check_gij_k0(eps,nbins,2,gij2);
    
    //step 3
    flt* Lk0_2 = malloc(nbins * sizeof(flt));
    Lk0_func_kdv(nbins,massgrid,gij2,dim_tabflux,tabflux,ind_tabflux,dv,flux,Lk0_2);
    for (unsigned int j=0;j<=nbins-1;j++){
        gij_new[j] = gij[j]/((flt)3) + ((flt)2)*(gij2[j] + dt*Lk0_2[j])/((flt)3);
    }

    //check negative values + limit values to eps
    check_gij_k0(eps,nbins,3,gij_new);

    //clear variables
    free(gij1);
    free(gij2);
    free(Lk0);
    free(Lk0_1);
    free(Lk0_2);
    gij1  = NULL;
    gij2  = NULL;
    Lk0   = NULL;
    Lk0_1 = NULL;
    Lk0_2 = NULL;
    
}



/**
 * @brief Coagulation CFL condition for DG scheme k>0 piecewise polynomial approximation (simple kernels)
 *
 * CFL formulation from Filbet & Laurencot 2004, dt <= mean_g * dm/dF
 *
 * @param    eps            minimum value for mass distribution approximation gij
 * @param    nbins          number of dust bins
 * @param    kpol           degree of polynomials for approximation
 * @param    massgrid       grid of masses, borders value of mass bins
 * @param    gij            components of g on the polynomial basis
 * @param    dim_tabflux    dimension of 1D array tabflux
 * @param    tabflux        1D array to evaluate coagulation flux
 * @param    ind_tabflux    1D array of indices corresponding to values in tabflux
 * @param    flux           1D array with approximation of the coagulation flux in each bin
 * @param    tabdtCFL       1D array with CFL timestep values in each bin
 * @return   dtCFL          CFL timestep restriction as min value of tabdtCFL, for coagulation, DG k=0
 */
flt compute_CFL(flt eps,unsigned int nbins,unsigned int kpol,flt* massgrid,flt* gij,
                unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
                flt* flux,flt* tabdtCFL){

    
    compute_flux(nbins,kpol,gij,dim_tabflux,tabflux,ind_tabflux,flux);
    // for (unsigned int j=0;j<nbins;j++){
    //     printf("j=%d, flux=%.15e\n",j,flux[j]);
    // }
    // exit(-1);



    //bin 0, where no mass entering the mass domain
    tabdtCFL[0] = fabs(gij[0]*(massgrid[1]-massgrid[0])/(flux[0]));

   
    for (unsigned int j=1;j<nbins;j++){
        if (gij[0+j*(kpol+1)] > eps){
            tabdtCFL[j] = fabs(gij[0+j*(kpol+1)]*(massgrid[j+1]-massgrid[j])/(flux[j]-flux[j-1]));
        }else{
            tabdtCFL[j] = ((flt)1e3);
        }
        
    }

    flt dtCFL = tabdtCFL[0];
    for(unsigned int i = 1; i < nbins; i++){
        dtCFL = fmin(dtCFL,tabdtCFL[i]);
    }

    if (dtCFL == ((flt)0)){
        printf("j,gij,hj,flux,tabdtCFL");
        for (unsigned int j=0; j<nbins; j++){
            printf("%d,%.15e,%.15e,%.15e,%.15e \n ",j,gij[0+j*(kpol+1)],massgrid[j+1]-massgrid[j],flux[j],tabdtCFL[j]);
        }
        printf("dtCFL = %.15e\n",dtCFL);
        printf("issue in CFL coagulation k>0 simple kernels \n");
        exit(-1);

    }
   
    return dtCFL;
   
}



/**
 * @brief compute the DG operator L for piecewise polynomial approximation and simple kernels (see Lombart et al., 2021) 
 *
 * It is used for the time solver
 *
 * @param    nbins          number of dust bins
 * @param    kpol           degree of polynomials for approximation
 * @param    massgrid       grid of masses, borders value of mass bins
 * @param    mat_coeffs_leg array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
 * @param    gij            components of g on the polynomial basis
 * @param    dim_tabflux    dimension of 1D array tabflux
 * @param    tabflux        1D array to evaluate coagulation flux
 * @param    ind_tabflux    1D array of indices corresponding to values in tabflux
 * @param    dim_tabintflux dimension of 1D array tabintflux
 * @param    tabintflux     1D array to evaluate the term including the integral of coagulation flux
 * @param    ind_tabintflux 1D array of indices corresponding to values in tabintflux
 * @param    flux           1D array with approximation of the coagulation flux in each bin
 * @param    intflux        1D array with approximation of the term including intergral of the flux in each bin
 * @return   Lk             DG operator for piecewise polynomial approximation in each bin
 */
void Lk_func(unsigned int nbins,unsigned int kpol,flt* massgrid,flt* mat_coeffs_leg,flt* gij,
                unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
                unsigned int dim_tabintflux,flt* tabintflux, int* ind_tabintflux,
                flt* flux,flt* intflux,flt* Lk){


    compute_flux_intflux(nbins,kpol,gij,
                            dim_tabflux,tabflux,ind_tabflux,
                            dim_tabintflux,tabintflux,ind_tabintflux,
                            flux,intflux);

    // for (unsigned int j=0;j<nbins;j++){
    //     printf("j=%d, flux=%.15e\n",j,flux[j]);
    // }

    // for (unsigned int j=0;j<nbins;j++){
    //     printf("j=%d, intflux=",j);
    //     for (unsigned int k=0;k<=kpol;k++){
    //         printf("%.15e  ",intflux[k+j*(kpol+1)]);

    //     }
    //     printf("\n");
    // }
    // exit(-1);
        
    flt ak[kpol+1];
    memset(ak, 0, (kpol + 1) * sizeof(flt));
    flt hj; flt c; flt LegPleft; flt LegPright;

    for (unsigned int k=0;k<=kpol;k++){

        // Safe bounds check
        if ( k > kpol) {
            fprintf(stderr, "Polynomial order exceeds kpol=%u (k=%u)\n", 
                 kpol,k);
        }


        c = coeffnorm(k);
        
        
        for (unsigned int p=0; p<=k;p++ ){
            ak[p] =  mat_coeffs_leg[p + k*(k+1)];
        }

        LegPleft  = phi_pol(k,ak,-((flt)1));
        LegPright = phi_pol(k,ak,((flt)1));

        for (unsigned int j=0;j<=nbins-1;j++){
            hj = massgrid[j+1]-massgrid[j];
        
            if (j==0){
                Lk[k+j*(kpol+1)] = (((flt)2)*(intflux[k+j*(kpol+1)]-(flux[j]*LegPright)))/(c*hj);
            }else{
                Lk[k+j*(kpol+1)] = (((flt)2)*(intflux[k+j*(kpol+1)]-(flux[j]*LegPright - flux[j-1]*LegPleft)))/(c*hj);
            }

        }

    }
    
}


/**
 * @brief Function to compute SSPRK order 3 time solver with piecewise polynomial approximation for simple kernels
 *
 * See Zhang & Shu 2010 and Lombart et al., 2021
 *
 * @param    eps            minimum value for mass distribution approximation gij
 * @param    nbins          number of dust bins
 * @param    kpol           degree of polynomials for approximation
 * @param    massgrid       grid of masses, borders value of mass bins
 * @param    massbins       arithmetic mean value of massgrid for each mass bins
 * @param    mat_coeffs_leg array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
 * @param    gij            components of g on the polynomial basis
 * @param    dim_tabflux    dimension of 1D array tabflux
 * @param    tabflux        1D array to evaluate coagulation flux
 * @param    ind_tabflux    1D array of indices corresponding to values in tabflux
 * @param    dim_tabintflux dimension of 1D array tabintflux
 * @param    tabintflux     1D array to evaluate the term including the integral of coagulation flux
 * @param    ind_tabintflux 1D array of indices corresponding to values in tabintflux
 * @param    flux           1D array with approximation of the coagulation flux in each bin
 * @param    intflux        1D array with approximation of the term including intergral of the flux in each bin
 * @param    tabgamma       1D array with limiter coefficient for reconstruction in each bin to ensure positivity
 * @param    dt             timestep
 * @return   gij_new        evolved components of g on the polynomial basis after 1 timestep
 */
void solver(flt eps,unsigned int nbins,unsigned int kpol,
            flt* massgrid, flt* massbins,flt* mat_coeffs_leg,flt* gij,
            unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
            unsigned int dim_tabintflux,flt* tabintflux, int* ind_tabintflux,
            flt* flux,flt* intflux,flt* tabgamma,flt dt,flt* gij_new){

    flt* gij1 = malloc(nbins*(kpol+1) * sizeof(flt));
    flt* gij2 = malloc(nbins*(kpol+1) * sizeof(flt));

    //step 1
    flt* Lk = malloc(nbins*(kpol+1) * sizeof(flt));
    Lk_func(nbins,kpol,massgrid,mat_coeffs_leg,gij,
                dim_tabflux,tabflux,ind_tabflux,
                dim_tabintflux,tabintflux,ind_tabintflux,
                flux,intflux,Lk);

    for (unsigned int j=0;j<=nbins-1;j++){
        for (unsigned int k=0;k<=kpol;k++){
            gij1[k+j*(kpol+1)] = gij[k+j*(kpol+1)] + dt*Lk[k+j*(kpol+1)];
        }
    }

    //apply scaling limiter
    gammafunction(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij1,tabgamma);
    for (unsigned int j=0;j<=nbins-1;j++){
        for (unsigned int k=1;k<=kpol;k++){
            gij1[k+j*(kpol+1)] = tabgamma[j]*gij1[k+j*(kpol+1)]; 
        }
    }



    //check negative values + limit values to eps
    check_gij(eps,nbins,kpol,1,gij1);
    
    //step 2
    flt* Lk_1 = malloc(nbins*(kpol+1) * sizeof(flt));
    Lk_func(nbins,kpol,massgrid,mat_coeffs_leg,gij1,
                dim_tabflux,tabflux,ind_tabflux,
                dim_tabintflux,tabintflux,ind_tabintflux,
                flux,intflux,Lk_1);
    for (unsigned int j=0;j<=nbins-1;j++){
        for (unsigned int k=0;k<=kpol;k++){
            gij2[k+j*(kpol+1)] = ((flt)3)*gij[k+j*(kpol+1)]/((flt)4) + (gij1[k+j*(kpol+1)] + dt*Lk_1[k+j*(kpol+1)])/((flt)4);
        }
    }

    //apply scaling limiter
    gammafunction(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij2,tabgamma);
    for (unsigned int j=0;j<=nbins-1;j++){
        for (unsigned int k=1;k<=kpol;k++){
            gij2[k+j*(kpol+1)] = tabgamma[j]*gij2[k+j*(kpol+1)]; 
        }
    }

    //check negative values + limit values to eps
    check_gij(eps,nbins,kpol,2,gij2);
    
    //step 3
    flt* Lk_2 = malloc(nbins*(kpol+1) * sizeof(flt));
    Lk_func(nbins,kpol,massgrid,mat_coeffs_leg,gij2,
                dim_tabflux,tabflux,ind_tabflux,
                dim_tabintflux,tabintflux,ind_tabintflux,
                flux,intflux,Lk_2);
    for (unsigned int j=0;j<=nbins-1;j++){
        for (unsigned int k=0;k<=kpol;k++){
            gij_new[k+j*(kpol+1)] = gij[k+j*(kpol+1)]/((flt)3) + ((flt)2)*(gij2[k+j*(kpol+1)] + dt*Lk_2[k+j*(kpol+1)])/((flt)3);
        }
    }


    //apply scaling limiter
    gammafunction(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij_new,tabgamma);
    for (unsigned int j=0;j<=nbins-1;j++){
        for (unsigned int k=1;k<=kpol;k++){
            gij_new[k+j*(kpol+1)] = tabgamma[j]*gij_new[k+j*(kpol+1)]; 
        }
    }



    //check negative values + limit values to eps
    check_gij(eps,nbins,kpol,3,gij_new);


    //clear variables
    free(gij1);
    free(gij2);
    free(Lk);
    free(Lk_1);
    free(Lk_2);
    gij1 = NULL;
    gij2 = NULL;
    Lk   = NULL;
    Lk_1 = NULL;
    Lk_2 = NULL;
    
}


/**
 * @brief Coagulation CFL condition for DG scheme k>0 piecewise polynomial approximation (physical kernel)
 *
 * CFL formulation from Filbet & Laurencot 2004, dt <= mean_g * dm/dF
 *
 * @param    eps            minimum value for mass distribution approximation gij
 * @param    nbins          number of dust bins
 * @param    kpol           degree of polynomials for approximation
 * @param    massgrid       grid of masses, borders value of mass bins
 * @param    gij            components of g on the polynomial basis
 * @param    dim_tabflux    dimension of 1D array tabflux
 * @param    tabflux        1D array to evaluate coagulation flux
 * @param    ind_tabflux    1D array of indices corresponding to values in tabflux
 * @param    dv             1D array of the differential velocity between grains in bins lp and l
 * @param    flux           1D array with approximation of the coagulation flux in each bin
 * @param    tabdtCFL       1D array with CFL timestep values in each bin
 * @return   dtCFL          CFL timestep restriction as min value of tabdtCFL, for coagulation, DG k=0
 */
flt compute_CFL_kdv(flt eps,unsigned int nbins,unsigned int kpol,flt* massgrid,flt* gij,
                    unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
                    flt* dv, flt* flux,flt* tabdtCFL){


    // for(unsigned int j = 0 ; j < nbins;j++){
    //     for(unsigned int k=0; k<=kpol;k++){
    //         printf("%.15e   ", gij[k+j*(kpol+1)]);
            
    //     }
    //     printf("\n");
        
    // }
    
    compute_flux_kdv(nbins,kpol,gij,dim_tabflux,tabflux,ind_tabflux,dv,flux);
    // printf("CFL flux\n");
    // for (unsigned int j=0;j<nbins;j++){
    //     printf("j=%d, flux=%.15e\n",j,flux[j]);
    // }
    // exit(-1);

    //bin 0, where no mass entering the mass domain
    tabdtCFL[0] = fabs(gij[0]*(massgrid[1]-massgrid[0])/(flux[0]));

   
    for (unsigned int j=1;j<nbins;j++){
        if (gij[0+j*(kpol+1)] > eps){
            tabdtCFL[j] = fabs(gij[0+j*(kpol+1)]*(massgrid[j+1]-massgrid[j])/(flux[j]-flux[j-1]));
        }else{
            tabdtCFL[j] = ((flt)1e3);
        }
        
    }

    flt dtCFL = tabdtCFL[0];
    for(unsigned int i = 1; i < nbins; i++){
        dtCFL = fmin(dtCFL,tabdtCFL[i]);
    }

    if (dtCFL == ((flt)0)){
        printf("j,gij,hj,flux,tabdtCFL");
        for (unsigned int j=0; j<nbins; j++){
            printf("%d,%.15e,%.15e,%.15e,%.15e \n ",j,gij[0+j*(kpol+1)],massgrid[j+1]-massgrid[j],flux[j],tabdtCFL[j]);
        }
        printf("dtCFL = %.15e\n",dtCFL);
        printf("issue in CFL coagulation k>0 simple kernels \n");
        exit(-1);

    }
   
    return dtCFL;
   
}


/**
 * @brief compute the DG operator L for piecewise polynomial approximation and physical kernel (see Lombart et al., 2021) 
 *
 * It is used for the time solver
 *
 * @param    nbins          number of dust bins
 * @param    kpol           degree of polynomials for approximation
 * @param    massgrid       grid of masses, borders value of mass bins
 * @param    mat_coeffs_leg array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
 * @param    gij            components of g on the polynomial basis
 * @param    dim_tabflux    dimension of 1D array tabflux
 * @param    tabflux        1D array to evaluate coagulation flux
 * @param    ind_tabflux    1D array of indices corresponding to values in tabflux
 * @param    dim_tabintflux dimension of 1D array tabintflux
 * @param    tabintflux     1D array to evaluate the term including the integral of coagulation flux
 * @param    ind_tabintflux 1D array of indices corresponding to values in tabintflux
 * @param    dv             1D array of the differential velocity between grains in bins lp and l
 * @param    flux           1D array with approximation of the coagulation flux in each bin
 * @param    intflux        1D array with approximation of the term including intergral of the flux in each bin
 * @return   Lk             DG operator for piecewise polynomial approximation in each bin
 */
void Lk_func_kdv(unsigned int nbins,unsigned int kpol,flt* massgrid,flt* mat_coeffs_leg,flt* gij,
                    unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
                    unsigned int dim_tabintflux,flt* tabintflux, int* ind_tabintflux,
                    flt* dv, flt* flux,flt* intflux,flt* Lk){

    compute_flux_intflux_kdv(nbins,kpol,gij,
                                dim_tabflux,tabflux,ind_tabflux,
                                dim_tabintflux,tabintflux,ind_tabintflux,
                                dv,flux,intflux);

    // for (unsigned int j=0;j<nbins;j++){
    //     printf("j=%d, flux=%.15e\n",j,flux[j]);
    // }

    // for (unsigned int j=0;j<nbins;j++){
    //     printf("j=%d, intflux=",j);
    //     for (unsigned int k=0;k<=kpol;k++){
    //         printf("%.15e  ",intflux[k+j*(kpol+1)]);

    //     }
    //     printf("\n");
    // }
    // exit(-1);

        
    flt ak[kpol+1];
    memset(ak, 0, (kpol + 1) * sizeof(flt));

    flt hj; flt c; flt LegPleft; flt LegPright;

    for (unsigned int k=0;k<=kpol;k++){

        // Safe bounds check
        if ( k > kpol) {
            fprintf(stderr, "Polynomial order exceeds kpol=%u (k=%u)\n",kpol,k);
        }

        c = coeffnorm(k);
        
        
        for (unsigned int p=0; p<=k;p++ ){
            ak[p] =  mat_coeffs_leg[p + k*(k+1)];
        }

        LegPleft  = phi_pol(k,ak,-((flt)1));
        LegPright = phi_pol(k,ak,((flt)1));

        for (unsigned int j=0;j<=nbins-1;j++){
            hj = massgrid[j+1]-massgrid[j];
        
            if (j==0){
                Lk[k+j*(kpol+1)] = (((flt)2)*(intflux[k+j*(kpol+1)]-(flux[j]*LegPright)))/(c*hj);
            }else{
                Lk[k+j*(kpol+1)] = (((flt)2)*(intflux[k+j*(kpol+1)]-(flux[j]*LegPright - flux[j-1]*LegPleft)))/(c*hj);
            }

        }

    }
    
}

/**
 * @brief Function to compute SSPRK order 3 time solver with piecewise polynomial approximation for physical kernel
 *
 * See Zhang & Shu 2010 and Lombart et al., 2021
 *
 * @param    eps            minimum value for mass distribution approximation gij
 * @param    nbins          number of dust bins
 * @param    kpol           degree of polynomials for approximation
 * @param    massgrid       grid of masses, borders value of mass bins
 * @param    massbins       arithmetic mean value of massgrid for each mass bins
 * @param    mat_coeffs_leg array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
 * @param    gij            components of g on the polynomial basis
 * @param    dim_tabflux    dimension of 1D array tabflux
 * @param    tabflux        1D array to evaluate coagulation flux
 * @param    ind_tabflux    1D array of indices corresponding to values in tabflux
 * @param    dim_tabintflux dimension of 1D array tabintflux
 * @param    tabintflux     1D array to evaluate the term including the integral of coagulation flux
 * @param    ind_tabintflux 1D array of indices corresponding to values in tabintflux
 * @param    dv             1D array of the differential velocity between grains in bins lp and l
 * @param    flux           1D array with approximation of the coagulation flux in each bin
 * @param    intflux        1D array with approximation of the term including intergral of the flux in each bin
 * @param    tabgamma       1D array with limiter coefficient for reconstruction in each bin to ensure positivity
 * @param    dt             timestep
 * @return   gij_new        evolved components of g on the polynomial basis after 1 timestep
 */
void solver_kdv(flt eps,unsigned int nbins,unsigned int kpol,flt* massgrid, flt* massbins,flt* mat_coeffs_leg,flt* gij,
                    unsigned int dim_tabflux,flt* tabflux,int* ind_tabflux, 
                    unsigned int dim_tabintflux,flt* tabintflux, int* ind_tabintflux,
                    flt* dv, flt* flux,flt* intflux,flt* tabgamma,flt dt,flt* gij_new){

    

    flt* gij1 = malloc(nbins*(kpol+1) * sizeof(flt));
    flt* gij2 = malloc(nbins*(kpol+1) * sizeof(flt));

    // printf("gij = \n");
    // for(unsigned int j = 0 ; j < nbins;j++){
    //     for(unsigned int k=0; k<=kpol;k++){
    //         printf("%.15e   ", gij[k+j*(kpol+1)]);
    //     }
    //     printf("\n");

    // }

    //step 1
    flt* Lk = malloc(nbins*(kpol+1) * sizeof(flt));
    Lk_func_kdv(nbins,kpol,massgrid,mat_coeffs_leg,gij,
                    dim_tabflux,tabflux,ind_tabflux,
                    dim_tabintflux,tabintflux,ind_tabintflux,
                    dv, flux,intflux,Lk);


    // printf("Lk = \n");
    // for(unsigned int j = 0 ; j < nbins;j++){
    //     for(unsigned int k=0; k<=kpol;k++){
    //         printf("%.15e   ", Lk[k+j*(kpol+1)]);
    //     }
    //     printf("\n");

    // }



    for (unsigned int j=0;j<=nbins-1;j++){
        for (unsigned int k=0;k<=kpol;k++){
            gij1[k+j*(kpol+1)] = gij[k+j*(kpol+1)] + dt*Lk[k+j*(kpol+1)];
        }
    }

    // printf("gij1 before gamma = \n");
    // for(unsigned int j = 0 ; j < nbins;j++){
    //     for(unsigned int k=0; k<=kpol;k++){
    //         printf("%.15e   ", gij1[k+j*(kpol+1)]);
    //     }
    //     printf("\n");

    // }

    //apply scaling limiter
    gammafunction(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij1,tabgamma);
    for (unsigned int j=0;j<=nbins-1;j++){
        for (unsigned int k=1;k<=kpol;k++){
            gij1[k+j*(kpol+1)] = tabgamma[j]*gij1[k+j*(kpol+1)]; 
        }
    }

    // printf("gij1 = \n");
    // for(unsigned int j = 0 ; j < nbins;j++){
    //     for(unsigned int k=0; k<=kpol;k++){
    //         printf("%.15e   ", gij1[k+j*(kpol+1)]);
    //     }
    //     printf("\n");

    // }


    //check negative values + limit values to eps
    check_gij(eps,nbins,kpol,1,gij1);

    
    //step 2
    flt* Lk_1 = malloc(nbins*(kpol+1) * sizeof(flt));
    Lk_func_kdv(nbins,kpol,massgrid,mat_coeffs_leg,gij1,
                    dim_tabflux,tabflux,ind_tabflux,
                    dim_tabintflux,tabintflux,ind_tabintflux,
                    dv,flux,intflux,Lk_1);
    for (unsigned int j=0;j<=nbins-1;j++){
        for (unsigned int k=0;k<=kpol;k++){
            gij2[k+j*(kpol+1)] = ((flt)3)*gij[k+j*(kpol+1)]/((flt)4) + (gij1[k+j*(kpol+1)] + dt*Lk_1[k+j*(kpol+1)])/((flt)4);
        }
    }

    //apply scaling limiter
    gammafunction(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij2,tabgamma);
    for (unsigned int j=0;j<=nbins-1;j++){
        for (unsigned int k=1;k<=kpol;k++){
            gij2[k+j*(kpol+1)] = tabgamma[j]*gij2[k+j*(kpol+1)]; 
        }
    }


    //check negative values + limit values to eps
    check_gij(eps,nbins,kpol,2,gij2);

    
    //step 3
    flt* Lk_2 = malloc(nbins*(kpol+1) * sizeof(flt));
    Lk_func_kdv(nbins,kpol,massgrid,mat_coeffs_leg,gij2,
                    dim_tabflux,tabflux,ind_tabflux,
                    dim_tabintflux,tabintflux,ind_tabintflux,
                    dv,flux,intflux,Lk_2);


    for (unsigned int j=0;j<=nbins-1;j++){
        for (unsigned int k=0;k<=kpol;k++){
            gij_new[k+j*(kpol+1)] = gij[k+j*(kpol+1)]/((flt)3) + ((flt)2)*(gij2[k+j*(kpol+1)] + dt*Lk_2[k+j*(kpol+1)])/((flt)3);
        }
    }



    //apply scaling limiter
    gammafunction(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij_new,tabgamma);

    for (unsigned int j=0;j<=nbins-1;j++){
        for (unsigned int k=1;k<=kpol;k++){
            gij_new[k+j*(kpol+1)] = tabgamma[j]*gij_new[k+j*(kpol+1)]; 
        }
    }


    //check negative values + limit values to eps
    check_gij(eps,nbins,kpol,3,gij_new);


    //clear variables
    free(gij1);
    free(gij2);
    free(Lk);
    free(Lk_1);
    free(Lk_2);
    gij1 = NULL;
    gij2 = NULL;
    Lk = NULL;
    Lk_1 = NULL;
    Lk_2 = NULL;
    
}




