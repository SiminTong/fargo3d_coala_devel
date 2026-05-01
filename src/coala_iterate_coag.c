//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************


#include "coala_precision.h"
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
#include "coala_io.h"
#include "coala_compute_coag.h"




/** @brief Iterate coagulation solver to reach the time ndthydro x dthydro  
 *
 * DG scheme k>0, piecewise polynomial approximation
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
 * @return   gij_new	   evolved gij at time=ndthydro x dthydro
 */
void iterate_coag(unsigned int nbins,unsigned int kpol,
					unsigned int Q,flt* vecnodes,flt* vecweights,
					flt* massgrid,flt* massbins,flt* massmeanlog,
					unsigned int ndthydro,flt dthydro){

	//variables for pre computing
	flt* tabflux_v0        = NULL;
	flt* tabintflux_v0 	   = NULL;
	int* ind_tabflux_v0    = NULL;
	int* ind_tabintflux_v0 = NULL;
	flt* tabflux           = NULL;
	flt* tabintflux        = NULL;
	int* ind_tabflux       = NULL;
	int* ind_tabintflux    = NULL;

	tabflux_v0 = malloc((nbins*nbins*nbins*(kpol+1)*(kpol+1)) * sizeof(flt));
	tabintflux_v0 = malloc((nbins*nbins*nbins*(kpol+1)*(kpol+1)*(kpol+1)) * sizeof(flt));

	ind_tabflux_v0 = malloc((nbins*nbins*nbins*(kpol+1)*(kpol+1))*5 * sizeof(int));
	ind_tabintflux_v0 = malloc((nbins*nbins*nbins*(kpol+1)*(kpol+1)*(kpol+1))*6 * sizeof(int));

	
	// NULL check allocations for precomputing
	if (!tabflux_v0 || !ind_tabflux_v0 || !tabintflux_v0 || !ind_tabintflux_v0 ) {
	    fprintf(stderr, "Malloc failed in for precomputing v0 in iterate_coag\n");
	    exit(-1);
	}

	


	// variables for coagulation solver
	flt* mat_coeffs_leg = NULL;
    flt* flux           = NULL;
    flt* intflux        = NULL;
    flt* tabdtCFL       = NULL;
    flt* tabgamma       = NULL;
    flt* gij            = NULL;
    flt* gij_new        = NULL;

    mat_coeffs_leg  = malloc(((kpol+1)*(kpol+1)) * sizeof(flt));
    tabdtCFL 	    = malloc((nbins)*sizeof(flt));
    tabgamma 		= malloc((nbins) * sizeof(flt));
    flux     		= malloc((nbins)*sizeof(flt));
    intflux  		= malloc((nbins*(kpol+1))*sizeof(flt));
    gij      		= malloc((nbins*(kpol+1))*sizeof(flt));
    gij_new  		= malloc((nbins*(kpol+1))*sizeof(flt));

    // NULL check allocations for variables for coag solver
	if (!mat_coeffs_leg || !tabdtCFL || !tabgamma || !flux || !intflux || !gij || !gij_new) {
	    fprintf(stderr, "Malloc failed for variables for coagulation solver in iterate_coag\n");
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

    
    //precomputing part: generate tabflux tabintflux for DG scheme 
    int dim_tabflux;
	int dim_tabintflux; 
	

	printf("Precomputing arrays ... \n");

    start = clock(); 
    
	compute_dim_coagtabflux_coagtabintflux_GQ(kernel,K0,Q,vecnodes,vecweights,
                                                nbins,kpol,massgrid,mat_coeffs_leg,
                                                &dim_tabflux,tabflux_v0,ind_tabflux_v0,
                                                &dim_tabintflux,tabintflux_v0,ind_tabintflux_v0);

	
	tabflux = malloc(dim_tabflux * sizeof(flt));
	tabintflux = malloc(dim_tabintflux * sizeof(flt));

	ind_tabflux = malloc(dim_tabflux*5 * sizeof(int));
	ind_tabintflux = malloc(dim_tabintflux*6 * sizeof(int));	

	// NULL check allocations for intflux variables
	if (!tabflux || !ind_tabflux || !tabintflux || !ind_tabintflux) {
	    fprintf(stderr, "Malloc failed in for precomputing in iterate_coag\n");
	    exit(-1);
	}									
	

	compute_coagtabflux_coagtabintflux_GQ(dim_tabflux,tabflux_v0,ind_tabflux_v0,
											dim_tabintflux,tabintflux_v0,ind_tabintflux_v0,
											tabflux,ind_tabflux,tabintflux,ind_tabintflux);

	end = clock();
	elapsed_seconds = (flt)(end - start) / CLOCKS_PER_SEC;
	printf("Tabflux tabintflux generated in %.3e s \n",elapsed_seconds);


	//log simu
    if (save){
        fp = fopen(path_log,"a");
        fprintf(fp,"Tabflux tabintflux generated in %.3e s \n",elapsed_seconds);
        fclose(fp);
    }


    start = clock();
	//generate gij (component on polynomial basis)
	
	L2proj_GQ(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,Q,vecnodes,vecweights,gij);
	
	end = clock();
	elapsed_seconds = (flt)(end - start) / CLOCKS_PER_SEC;


	//apply scaling limiter
	
    gammafunction(eps,nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,tabgamma);
    for (unsigned int j=0;j<=nbins-1;j++){
        for (unsigned int k=1;k<=kpol;k++){
            gij[k+j*(kpol+1)] = tabgamma[j]*gij[k+j*(kpol+1)];
        }
    }

    // for(unsigned int j = 0 ; j < nbins;j++){
    //     printf("tabgamma = %.15e\n", tabgamma[j]);
    // }

    //limit values to eps
    check_gij(eps,nbins,kpol,1,gij);

    //total mass t0
    flt M1_t0 = (flt) 0.;
   	for (unsigned int j=0;j<nbins;j++){
     	M1_t0 += (massgrid[j+1]-massgrid[j])*gij[0+j*(kpol+1)];
   	}

	printf("Gij start in %.3e s \n",elapsed_seconds);
	printf("gij start = \n");
	for(unsigned int j = 0 ; j < nbins;j++){
		for (unsigned int k=0; k<=kpol; k++){
			printf("%.15e   ", gij[k+j*(kpol+1)]);
		}
		printf("\n");
		
	}

	printf("\n");


	//write gij start
	if (save){
		fp = fopen(path_gij_t0,"w");
		for (unsigned int j=0;j<nbins;j++){
			for (unsigned int k=0; k<=kpol; k++){
				fprintf(fp,"%.15e   ",gij[k+j*(kpol+1)]);
			}
			fprintf(fp,"\n");
			
		}
		fclose(fp);
	}

	
	//write gmassmeanlog start
	flt gmassmeanlog;
	if (save){

		fp = fopen(path_gt0_massmeanlog,"w");
		for (unsigned int j=0;j<nbins;j++){
			gmassmeanlog = grecons(nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,j,massmeanlog[j]);
			fprintf(fp,"%.15e\n",gmassmeanlog);
		}
		fclose(fp);


	}

	
	printf("\n");
	printf("\033[96m>>>Time solver<<<\033[0m \n");


	//iterate solver time    
	//define variables for time solver
    flt time;
    unsigned int iprogress = 0;
    int nsub, ndt, tot_nsub, tot_ndt;  
    tot_nsub = 0;
    tot_ndt  = 0;
    
    
	//write time and gij
	time = ((flt)0);
	if (save){
		fp = fopen(path_time,"w");
		fprintf(fp,"%.15e\n",time);
		fclose(fp);

		fp = fopen(path_gij,"w");
		for (unsigned int j=0;j<nbins;j++){
			for (unsigned int k=0; k<=kpol; k++){
				fprintf(fp,"%.15e   ",gij[k+j*(kpol+1)]);
			}
			fprintf(fp,"\n");
			
		}
		fclose(fp);
	}
	
	start = clock();
	for (unsigned int it=0;it<ndthydro;it++){

		compute_coag(eps,coeff_CFL,nbins,kpol,
                    massgrid,massbins,mat_coeffs_leg,gij,
                    dim_tabflux,tabflux,ind_tabflux, 
                    dim_tabintflux,tabintflux,ind_tabintflux,
                    flux,intflux,tabgamma,tabdtCFL,dthydro,
                    gij_new,&nsub,&ndt);

		for (unsigned int i = 0; i < nbins * (kpol + 1); i++) {
		    gij[i] = gij_new[i];
		}

		tot_nsub += nsub;
      	tot_ndt  += ndt;
		
		time += dthydro;
		if (save){
			fp = fopen(path_time,"a");
			fprintf(fp,"%.15e\n",time);
			fclose(fp);

			fp = fopen(path_gij,"a");
			for (unsigned int j=0;j<nbins;j++){
				for (unsigned int k=0; k<=kpol; k++){
					fprintf(fp,"%.15e   ",gij[k+j*(kpol+1)]);
				}
				fprintf(fp,"\n");
				
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
		for(unsigned int k=0; k<=kpol;k++){
			printf("%.15e   ", gij[k+j*(kpol+1)]);
		}
		printf("\n");
		
	}
	printf("\n");

	//write gij end
	if (save){
		fp = fopen(path_gij_tend,"w");
		for (unsigned int j=0;j<nbins;j++){
			for(unsigned int k=0; k<=kpol;k++){
				fprintf(fp,"%.15e   ", gij[k+j*(kpol+1)]);
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
		
	}
	
	//write gmassmeanlog end
	if (save){
		fp = fopen(path_gtend_massmeanlog,"w");
		for (unsigned int j=0;j<nbins;j++){
			gmassmeanlog = grecons(nbins,kpol,massgrid,massbins,mat_coeffs_leg,gij,j,massmeanlog[j]);
			fprintf(fp,"%.15e\n",gmassmeanlog);
		}
		fclose(fp);
		
	}


	//check mass conservation
	flt M1_tend = (flt) 0.;
   	for (unsigned int j=0;j<nbins;j++){
     	M1_tend += (massgrid[j+1]-massgrid[j])*gij[0+j*(kpol+1)];
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
	free(tabintflux_v0);
	free(ind_tabflux_v0);
	free(ind_tabintflux_v0);
	free(tabflux);
	free(tabintflux);
	free(ind_tabflux);
	free(ind_tabintflux);
	free(mat_coeffs_leg);
    free(tabdtCFL);
    free(tabgamma);
    free(flux);
    free(intflux);
    free(gij);
    free(gij_new);
	
	tabflux_v0        = NULL;
	tabintflux_v0     = NULL;
	ind_tabflux_v0    = NULL;
	ind_tabintflux_v0 = NULL;
	tabflux           = NULL;
	tabintflux        = NULL;
	ind_tabflux       = NULL;
	ind_tabintflux    = NULL;

	mat_coeffs_leg = NULL;
    tabdtCFL       = NULL;
    tabgamma       = NULL;
    flux           = NULL;
    intflux        = NULL;
    gij            = NULL;
    gij_new        = NULL;


	
}
