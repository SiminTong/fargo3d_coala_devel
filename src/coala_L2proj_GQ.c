//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

//-------------------------------------------
// L2 projection on Legendre polynomials basis 
// for initial condition in DG scheme Eq.19 with gauss legendre quadrature
// Equations refere to Lombart & Laibe (2020) (doi:10.1093/mnras/staa3682)
//-------------------------------------------

#include <math.h>
#include "coala_precision.h"
#include "coala_polynomials_legendre.h"
#include "coala_L2proj_GQ.h"
#include <stdlib.h>
#include <stdio.h>
#include "coala_tension_module.h"
#include <string.h>


/** @brief Compute the initial gij coefficient from the dimensionless initial function g(x)=x*exp(-x)
 *
 * DG scheme k=0, piecewise constant approximation
 *
 * @param    eps           minimum value for mass distribution approximation gij
 * @param    nbins         number of dust bins
 * @param    massgrid      grid of masses, borders value of mass bins
 * @param    massbins      arithmetic mean value of massgrid for each mass bins
 * @param    Q             number of points for Gauss-Legendre quadrature
 * @param    vecnodes      nodes of the Legendre polynomials
 * @param    vecweights    weights coefficients for the Gauss-Legendre polynomials
 * @return   gij           initial components of g on the polynomial basis
 */
void L2proj_GQ_k0(flt eps,unsigned int nbins,flt* massgrid, flt* massbins,unsigned int Q,flt* vecnodes,flt* vecweights, flt* gij){
	
	flt xj,hj,c,sum,LegP,xjalpha;
	for (unsigned int j=0;j<=nbins-1;j++){
		hj = massgrid[j+1]-massgrid[j];
		xj = massbins[j];

		
		c = coeffnorm(0);
		sum = ((flt)0);
			
		for (unsigned int alpha=0;alpha<Q;alpha++){
			xjalpha = xj + hj*vecnodes[alpha]/((flt)2);

			sum += vecweights[alpha]*xjalpha*exp(-xjalpha);
		
		}

		if (sum < ((flt)0)){
			printf("j=%d, gij[j]=%.15e\n",j,gij[j]);
			printf("Negative value in function L2proj_GQ_k0");
			exit(-1);
		}else{
			gij[j] = fmax(eps,sum/c);
		}
		
	}
}




/** @brief Compute the initial gij coefficient from the dimensionless initial function g(x)=x*exp(-x)
 *
 * DG scheme k>0, piecewise polynomial approximation
 *
 * @param    eps               minimum value for mass distribution approximation gij
 * @param    nbins             number of dust bins
 * @param    kpol              degree of polynomials for approximation
 * @param    massgrid      	   grid of masses, borders value of mass bins
 * @param    massbins      	   arithmetic mean value of massgrid for each mass bins
 * @param    mat_coeffs_leg    array containing on each line Legendre polynomial coefficients from degree 0 to kpol, on each line coefficients are ordered from low to high orders
 * @param    Q                 number of points for Gauss-Legendre quadrature
 * @param    vecnodes          nodes of the Legendre polynomials
 * @param    vecweights        weights coefficients for the Gauss-Legendre polynomials
 * @return   gij               initial components of g on the polynomial basis
 */
void L2proj_GQ(flt eps,unsigned int nbins,unsigned int kpol,flt* massgrid, flt* massbins,
				flt* mat_coeffs_leg,unsigned int Q,flt* vecnodes,flt* vecweights, flt* gij){
	
	flt ak[kpol+1];
	memset(ak, 0, (kpol + 1) * sizeof(flt));

	flt xj,hj,c,sum,LegP,xjalpha;
	for (unsigned int j=0;j<=nbins-1;j++){
		hj = massgrid[j+1]-massgrid[j];
		xj = massbins[j];

		for (unsigned int k=0;k<=kpol;k++){

		    // Safe bounds check
		    if ( k > kpol) {
		    	fprintf(stderr, "Polynomial order exceeds kpol=%u (k=%u)\n",kpol, k);
		   	} 


	        if (k > 0){
			    for (unsigned int p=0; p<=k;p++ ){
			        ak[p] =  mat_coeffs_leg[p + k*(k+1)];

			    }
			} else{
			      ak[k] = mat_coeffs_leg[0];
			}

			c = coeffnorm(k);
			sum = ((flt)0);


			for (unsigned int alpha=0;alpha<Q;alpha++){
				xjalpha = xj + hj*vecnodes[alpha]/((flt)2);
				LegP = phi_pol(k,ak,vecnodes[alpha]);
				sum += vecweights[alpha]*xjalpha*exp(-xjalpha)*LegP;
			
			}

			if (k==0){
				if (sum < ((flt)0)){
					printf("j=%d, gij[0+j*(kpol+1)]k=%.15e\n",j,gij[0+j*(kpol+1)]);
					printf("Negative value in function L2proj_GQ");
					exit(-1);
				}else if (sum/c < eps){
					gij[0+j*(kpol+1)] = eps;
					for (unsigned int ik=1;ik<=kpol;ik++){
						gij[ik+j*(kpol+1)] = ((flt)0);
					}
					
					break;
				}else{
					gij[k+j*(kpol+1)] = sum/c;
				}
			}else {
				gij[k+j*(kpol+1)] = sum/c;
			}



		}	
	}
}




/** @brief Compute the initial gij coefficient from interpolation of cumulative function of the dust density
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
 * @param    eps_rhodust       minimum value for dust density
 * @param    rhodust           1D array dust density for each grain size
 * @return   gij               initial components of g on the polynomial basis
 */
void L2proj_gij_GQ(unsigned int nbins,unsigned int kpol,flt* massgrid, flt* massbins,
					flt* mat_coeffs_leg,unsigned int Q,flt* vecnodes,flt* vecweights,
					flt eps_rhodust, flt* rhodust, flt eps_gij,flt* gij){
	
	//interpolation to get gij

	flt sum;
	flt Vfunc[nbins+1];
	Vfunc[0] = ((flt)0);
	for (unsigned int j=1;j<=nbins;j++){
		sum = ((flt)0);
		for (unsigned int i=0;i<j;i++){
			sum += rhodust[i];
		}
		Vfunc[j] = sum;
	}
	

	flt YP[nbins+1];
	flt SIGMA[nbins+1];
	flt log_massgrid[nbins+1];
	for (unsigned int j=0;j<nbins+1;j++){
		YP[j] = ((flt)0);
		SIGMA[j] = ((flt)0);
		log_massgrid[j] = log(massgrid[j]);
	}
	int IER;
	int SigErr;
	//interpolation Vfunc(log(x))
	TSPSI(nbins+1,log_massgrid,Vfunc,YP,SIGMA,IER,SigErr);

	// printf("TSPSI IER=%d, SigErr=%d \n ",IER,SigErr);
	for (unsigned int j=0;j<nbins;j++){
		for (unsigned int k=0;k<=kpol;k++){
			gij[k+j*(kpol+1)] = ((flt)0);
		}
	}


	flt ak[kpol+1];
	memset(ak, 0, (kpol + 1) * sizeof(flt));

	flt xj,hj,c,LegP,xjalpha,interp_func;

	
	for (unsigned int j=0;j<nbins;j++){

		if (rhodust[j] > eps_rhodust){

			hj = massgrid[j+1]-massgrid[j];
			xj = massbins[j];

			for (unsigned int k=0;k<=kpol;k++){

				// Safe bounds check
			    if ( k > kpol) {
			    	fprintf(stderr, "Polynomial order exceeds kpol=%u (k=%u)\n", 
			             	kpol,k);
			   	} 

		        for (unsigned int p=0; p<=k;p++ ){
		            ak[p] =  mat_coeffs_leg[p + k*(k+1)];
		        }

				c = coeffnorm(k);
				sum = ((flt)0);
				
				for (unsigned int alpha=0;alpha<Q;alpha++){
					xjalpha = xj + hj*vecnodes[alpha]/((flt)2);
					LegP = phi_pol(k,ak,vecnodes[alpha]);

					//interpolation V(log(x))
					interp_func = HPVAL(log(xjalpha),nbins+1,log_massgrid,Vfunc,YP,SIGMA,IER)/xjalpha;
					// exit(-1);
					// printf("j=%d, k=%d, alpha=%d, interp_func=%.15e \n",j,k,alpha,interp_func);
					// printf("IER=%d ",IER);

					sum += vecweights[alpha]*interp_func*LegP;
				
				}

				if (k==0){
					if (sum < ((flt)0)){
						printf("j=%d, gij[0+j*(kpol+1)]k=%.15e\n",j,gij[0+j*(kpol+1)]);
						printf("Negative value in function L2proj_gij_GQ");
						exit(-1);
					}else if (sum/c < eps_gij){
						gij[0+j*(kpol+1)] = eps_gij;
						for (unsigned int ik=1;ik<=kpol;ik++){
							gij[ik+j*(kpol+1)] = ((flt)0);
						}
						
						break;
					}else{
						gij[k+j*(kpol+1)] = sum/c;
					}
				}else {
					gij[k+j*(kpol+1)] = sum/c;
				}


			}	
		}

		//enforce mass conservation after interpolation + L2 projection
        gij[0+j*(kpol+1)] = rhodust[j]/(massgrid[j+1]-massgrid[j]);



	}

	
}


