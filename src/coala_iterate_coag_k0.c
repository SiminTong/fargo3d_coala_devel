//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

#include "coala_precision.h"
#include "coala_io.h"
#include "coala_polynomials_legendre.h"
#include "coala_reconstruction_g.h"
#include "coala_generate_tabflux_tabintflux_GQ.h"
#include "coala_L2proj_GQ.h"
#include <stdio.h>
#include <time.h>
#include "coala_setup.h"
#include <stdlib.h>
#include "string.h"
#include "coala_limiter.h"
#include "coala_progress_bar.h"
#include "coala_compute_coag.h"



/** @brief Iterate coagulation solver to reach the time ndthydro x dthydro  
 *
 * DG scheme k=0, piecewise constant approximation
 *
 * @param    nbins         number of dust bins
 * @param    kpol          degree of polynomials for approximation
 * @param    Q             number of points for Gauss-Legendre quadrature
 * @param    vecnodes      nodes of the Legendre polynomials
 * @param    vecweights    weights coefficients for the Gauss-Legendre polynomials
 * @param    massgrid      grid of masses, borders value of mass bins
 * @param    massbins      arithmetic mean value of massgrid for each mass bins
 * @param    massmeanlog   geometric mean value of massgrid for each mass bins
 * @param    ndthydro      number of hydro timestep
 * @param    dthydro       hydro timestep
 * @return   gij_new       evolved gij at time=ndthydro x dthydro
 */
void iterate_coag_k0(unsigned int nbins,unsigned int kpol,
                        unsigned int Q,flt* vecnodes,flt* vecweights,
                        flt* massgrid,flt* massbins,flt* massmeanlog,
                        unsigned int ndthydro,flt dthydro){

    //variables for pre computing
    flt* tabflux_v0     = NULL;
    int* ind_tabflux_v0 = NULL;
    flt* tabflux        = NULL;
    int* ind_tabflux    = NULL;

    tabflux_v0  = malloc((nbins*nbins*nbins) * sizeof(flt));
    ind_tabflux_v0 = malloc((nbins*nbins*nbins)*3 * sizeof(int));

    
    // NULL check allocations for precomputing
    if (!tabflux_v0 || !ind_tabflux_v0) {
        fprintf(stderr, "Malloc failed in for precomputing v0 in iterate_coag_k0\n");
        exit(-1);
    }


    flt* mat_coeffs_leg = NULL;
    flt* flux           = NULL;
    flt* tabdtCFL       = NULL;
    flt* gij            = NULL;
    flt* gij_new        = NULL;

    mat_coeffs_leg  = malloc(((kpol+1)*(kpol+1)) * sizeof(flt));
    tabdtCFL        = malloc((nbins)*sizeof(flt));
    flux            = malloc((nbins)*sizeof(flt));
    gij             = malloc((nbins)*sizeof(flt));
    gij_new         = malloc((nbins)*sizeof(flt));

    // NULL check allocations for variables for coag solver
    if (!mat_coeffs_leg || !tabdtCFL || !flux || !gij || !gij_new) {
        fprintf(stderr, "Malloc failed for variables for coagulation solver in iterate_coag_k0\n");
        exit(-1);
    }

    
    
    // to write data
    FILE *fp;

    // Using time point and system_clock
    clock_t start_tot, end_tot;
    clock_t start, end;
    flt elapsed_seconds;

    
    //coefficients for Legendre polynomials
    compute_mat_coeffs(kpol,mat_coeffs_leg);

    //precomputation part: generate tabflux for DG scheme 
    int dim_tabflux;

    printf("Precomputing arrays ... \n");
    
    start = clock();

    
    
    compute_dim_coagtabflux_GQ_k0(kernel,K0,Q,vecnodes,vecweights,
                                    nbins,kpol,massgrid,mat_coeffs_leg,
                                    &dim_tabflux,tabflux_v0,ind_tabflux_v0);

    tabflux     = malloc(dim_tabflux * sizeof(flt));
    ind_tabflux = malloc(dim_tabflux*3 * sizeof(int));
    // NULL check allocations for precomputing
    if (!tabflux || !ind_tabflux) {
        fprintf(stderr, "Malloc failed in for precomputing in iterate_coag_k0\n");
        exit(-1);
    }


    compute_coagtabflux_GQ_k0(dim_tabflux,tabflux_v0,ind_tabflux_v0,tabflux,ind_tabflux);

    end = clock();
    elapsed_seconds = (flt)(end - start) / CLOCKS_PER_SEC;

    printf("Tabflux generated in %.3e s \n",elapsed_seconds);
    
    
    //write log simu
    if (save){
        fp = fopen(path_log,"a");
        fprintf(fp,"Tabflux generated in %.3e s \n",elapsed_seconds);
        fclose(fp);
    }


    start = clock();
    //generate gij (component on polynomial basis)
    L2proj_GQ_k0(eps,nbins,massgrid, massbins,Q,vecnodes,vecweights,gij);
    
    end = clock();
    elapsed_seconds = (flt)(end - start) / CLOCKS_PER_SEC;


    check_gij_k0(eps,nbins,0,gij);

    //total mass t0
    flt M1_t0 = (flt) 0.;
    for (unsigned int j=0;j<nbins;j++){
        M1_t0 += (massgrid[j+1]-massgrid[j])*gij[j];
    }
    
    
    printf("Gij start in %.3e s \n",elapsed_seconds);
    printf("gij start = \n");
    for(unsigned int j = 0 ; j < nbins;j++){
        printf("%.15e\n", gij[j]);
    }
    printf("\n");

    //write gij start
    if (save){
        fp = fopen(path_gij,"w");
        for(unsigned int j = 0 ; j < nbins;j++){
            fprintf(fp,"%.15e\n", gij[j]);
        }
        fclose(fp);
    }

    printf("\n");
    printf("\033[96m>>>Time solver<<<\033[0m \n");

    
    //iterate solver time
    //define variables for time solver
    flt time;
    unsigned int iprogress=0;  
    int nsub, ndt, tot_nsub, tot_ndt;  
    tot_nsub = 0;
    tot_ndt  = 0;
  
    
    //write time and gij
    time = ((flt)0);
    if (save){
        fp = fopen(path_time,"w");
        fprintf(fp,"%.15e\n",time);
        fclose(fp);

        fp = fopen(path_gij_t0,"w");
        for(unsigned int j = 0 ; j < nbins;j++){
            fprintf(fp,"%.15e\n", gij[j]);
        }
        fclose(fp);
    }


    start = clock();
    for (unsigned int it=0;it<ndthydro;it++){
        //compute coagulation for each hydro time step
        compute_coag_k0(eps,coeff_CFL,nbins,massgrid,gij,
                            dim_tabflux,tabflux,ind_tabflux,
                            flux,tabdtCFL,dthydro,
                            gij_new,&nsub,&ndt);

        for (unsigned int j=0; j<nbins; j++){
            gij[j] = gij_new[j];
        }

        tot_nsub += nsub;
        tot_ndt  += ndt;

        time += dthydro;
        if (save){
            fp = fopen(path_time,"a");
            fprintf(fp,"%.15e\n",time);
            fclose(fp);

            fp = fopen(path_gij,"a");
            for(unsigned int j = 0 ; j < nbins;j++){
                fprintf(fp,"%.15e\n", gij[j]);
            }
            fclose(fp);
        }
    
        progressbar(ndthydro,iprogress);
        iprogress += 1;

    }

    end = clock();
    elapsed_seconds = (flt)(end - start) / CLOCKS_PER_SEC;

    printf("\n");
    printf("gij end = \n");
    for(unsigned int j = 0 ; j < nbins;j++){
        printf("%.15e\n", gij[j]);
    }
    printf("\n");

    //write gij end
    if (save){
        fp = fopen(path_gij_tend,"w");
        for(unsigned int j = 0 ; j < nbins;j++){
            fprintf(fp,"%.15e\n", gij[j]);
        }
        fclose(fp);
        
    }


    //check mass conservation
    flt M1_tend = (flt) 0.;
    for (unsigned int j=0;j<nbins;j++){
        M1_tend += (massgrid[j+1]-massgrid[j])*gij[j];
    }
    printf("M1 t0 = %.3e \n",M1_t0);
    printf("M1 tend = %.3e \n",M1_tend);
    printf("diff M1 = %.3e \n",M1_tend-M1_t0);
    printf("\n");
    printf("Number of sub-timestep < dthydro for coagulation = %d \n",tot_nsub);
    printf("Number of timestep at dthydro = %d \n",tot_ndt);
    printf("Total number timesteps = %d \n",tot_ndt+tot_nsub);
    printf("Time solver in %.3e s \n",elapsed_seconds);
    printf("Time per timestep in %.3e s \n",elapsed_seconds/(tot_ndt+tot_nsub));


    //free pointers
    free(tabflux_v0);
    free(ind_tabflux_v0);
    free(tabflux);
    free(ind_tabflux);
    free(mat_coeffs_leg);
    free(tabdtCFL);
    free(gij);
    free(flux);

    tabflux_v0     = NULL;
    ind_tabflux_v0 = NULL;
    tabflux        = NULL;
    ind_tabflux    = NULL;
    mat_coeffs_leg = NULL;

    tabdtCFL = NULL;
    gij      = NULL;
    flux     = NULL;

    
}

