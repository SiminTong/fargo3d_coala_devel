//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
#include "coala_precision.h"
#include "coala_generate_tabflux_tabintflux_GQ.h"
#include "coala_GQ_legendre_nodes_weights.h"
#include "coala_init_grid.h"
#include "coala_polynomials_legendre.h"
#include "coala_interface_coag.h"
//<\INCLUDES>


void Coala_init(){

    int nbins = NFLUIDS-1;

    massgrid    = malloc((nbins+1)*sizeof(flt));
    massbins    = malloc((nbins)*sizeof(flt));
    massmeanlog = malloc((nbins)*sizeof(flt));

    unsigned int kernel = 3;
    flt rhograin = RHOSOLID/(MSTAR_CGS/(R0_CGS*R0_CGS*R0_CGS))*(MSTAR/(R0*R0*R0));
    flt K0              = M_PI*pow((4./3.)*M_PI*rhograin,-2./3.);
    unsigned int Q      = 5;
    int kpol            = 0;
    // in K0, rhograin has to be given in hydro code unit or physical unit. It needs to be coherent with the mass unit

    // Need to define array for the Gauss-Quadrature method
    vecnodes   = malloc(Q*sizeof(flt));
    vecweights = malloc(Q*sizeof(flt));

    GQLeg_nodes(Q,vecnodes);
    GQLeg_weights(Q,vecweights);

    init_grid(nbins,SMAX*R0/R0_CGS,SMIN*R0/R0_CGS,rhograin,MSTAR,massgrid,massbins);

  
    tabflux_v0    = malloc((nbins*nbins*nbins) * sizeof(flt));
    ind_tabflux_v0 = malloc((nbins*nbins*nbins)*3 * sizeof(int));
    mat_coeffs_leg  = malloc(((kpol+1)*(kpol+1)) * sizeof(flt));

    compute_mat_coeffs(kpol,mat_coeffs_leg);
    compute_dim_coagtabflux_GQ_k0(kernel,K0,Q,vecnodes,vecweights,
                                nbins,kpol,massgrid,mat_coeffs_leg,
                                &dim_tabflux,tabflux_v0,ind_tabflux_v0);

    tabflux = malloc(dim_tabflux * sizeof(flt));
    ind_tabflux = malloc(dim_tabflux*3 * sizeof(int));
    compute_coagtabflux_GQ_k0(dim_tabflux,tabflux_v0,ind_tabflux_v0,tabflux,ind_tabflux);

}

void Coala_main(real dt){

//<USER_DEFINED>
   real *densgas;
   real *densities[NFLUIDS-1];
   real *velocitiesx[NFLUIDS-1];
   real *velocitiesy[NFLUIDS-1];
   real *velocitiesz[NFLUIDS-1];
   real *stokes[NFLUIDS-1];
   real eps_rhodust = 1.0e-10;

   int ii;

   densgas = Fluids[0]->Density->field_cpu;

   for (ii=1; ii<NFLUIDS; ii++) {

    INPUT(Fluids[ii]->Density);
    densities[ii-1]  = Fluids[ii]->Density->field_cpu;
    
    INPUT(Fluids[ii]->Energy);
    stokes[ii-1]  = Fluids[ii]->Energy->field_cpu;
    
#ifdef X
    INPUT(Fluids[ii]->Vx_temp);
    velocitiesx[ii-1] = Fluids[ii]->Vx_temp->field_cpu;
#endif

#ifdef Y
    INPUT(Fluids[ii]->Vy_temp);
    velocitiesy[ii-1] = Fluids[ii]->Vy_temp->field_cpu;
#endif

#ifdef Z
    INPUT(Fluids[ii]->Vz_temp);
    velocitiesz[ii-1] = Fluids[ii]->Vz_temp->field_cpu;
#endif
  }
//<\USER_DEFINED>

//<EXTERNAL>
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = XIP; 
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  int nbins = NFLUIDS-1;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  int o;
  int p;
  real hdust;
  real stokesp;
  real stokeso;
  real dvr;
  flt rhodust[NFLUIDS-1];
  flt invsize[NFLUIDS-1];
  flt rhos[NFLUIDS-1];
  flt new_rhodust[NFLUIDS-1];
  flt dv[(NFLUIDS-1)*(NFLUIDS-1)];  
//<\INTERNAL>

//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for(k=1; k<size_z; k++) {
#endif
#ifdef Y
    for(j=1; j<size_y; j++) {
#endif
#ifdef X
      for(i=XIM; i<size_x; i++) {
#endif
//<#>

	ll = l;

    //fill arrays with densities and relative velocity        
    for (o=0; o<NFLUIDS-1; o++) {
      
      //surface density to volumetric density
      hdust      = MIN(1.,sqrt(ALPHA/(ALPHA + stokes[o][ll])));
      rhodust[o] = densities[o][ll]/(sqrt(2*M_PI)*ASPECTRATIO*ymed(j)*hdust);    

        for (p=0; p<NFLUIDS-1; p++) {
            dvr = 0.;
    #ifdef X
            dvr += (velocitiesx[o][ll]-velocitiesx[p][ll])*(velocitiesx[o][ll]-velocitiesx[p][ll]);
    #endif
    #ifdef Y
            dvr += (velocitiesy[o][ll]-velocitiesy[p][ll])*(velocitiesy[o][ll]-velocitiesy[p][ll]);
    #endif
    #ifdef Z
            dvr += (velocitiesz[o][ll]-velocitiesz[p][ll])*(velocitiesz[o][ll]-velocitiesz[p][ll]);
    #endif
            dv[p+o*(NFLUIDS-1)] = pow(dvr, 0.5);
        }
    }

    // COALA KERNEL
    coala_coag_k0(nbins,massgrid,dim_tabflux,tabflux,ind_tabflux,rhodust,eps_rhodust,dv,dt,new_rhodust);

    // Update density with coala solution
    for (o=0; o<NFLUIDS-1; o++) {
        hdust      = MIN(1.,sqrt(ALPHA/(ALPHA + stokes[o][ll])));
       densities[o][ll] = new_rhodust[o]*(sqrt(2*M_PI)*ASPECTRATIO*ymed(j)*hdust);
	}
//<\#>
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
//<\MAIN_LOOP>
}

 