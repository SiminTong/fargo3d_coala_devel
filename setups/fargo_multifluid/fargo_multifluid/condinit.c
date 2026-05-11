#include "fargo3d.h"

void _CondInit(int id) {

  int i,j,k;
  real r, omega;
  
  real *rho  = Density->field_cpu;
  real *cs   = Energy->field_cpu;
  real *vphi = Vx->field_cpu;
  real *vr   = Vy->field_cpu;
  
  real rhog, rhod;
  real vk;
  real stokes_plus[NFLUIDS];
  real stokes[NFLUIDS-1];
  real epsilons[NFLUIDS-1];
  real sq = SQ;
  real slope;
  // Stokes numbers                                                                      
  real smax = TSMAX;
  real smin = TSMIN;
  real ds   = (log(smax)-log(smin))/(NFLUIDS-1);
  for(int n=0;n<NFLUIDS;n++){
    stokes_plus[n] = smin*exp(ds*n);
  }

  //Dust size-distribution                                                                 
  slope = 4-sq;
  for(int n=0; n<NFLUIDS-1; n++){
    if(slope != 0.) {
      epsilons[n]  = pow(stokes_plus[n+1],slope) - pow(stokes_plus[n],slope) ;
      epsilons[n] *= EPSILON/(pow(smax, slope) - pow(smin,slope));
    }
    else{
      epsilons[n]  = log(stokes_plus[n+1]/stokes_plus[n]);
      epsilons[n] *= EPSILON/log(smax/smin);
    }
    stokes[n] = stokes_plus[n+1];
    if( NFLUIDS == 2) stokes[n] = TSMAX;
  }

#ifdef DRAGFORCE
  if(id > 0) {
#ifdef STOKESNUMBER
   Coeffval[0]   = 1.0/stokes[id-1];
#endif
#ifdef DUSTSIZE
    Coeffval[1]   = 1.0/(stokes[id-1]*R0/R0_CGS);    
    Coeffval[2]   = RHOSOLID/(MSTAR_CGS/(R0_CGS*R0_CGS*R0_CGS))*(MSTAR/(R0*R0*R0));
#endif
    if(CPU_Master) printf("Ts %1.16f \t eps %1.16f \n", stokes[id-1], epsilons[id-1]);
  }
#endif
  i = j = k = 0;
  
  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      for (i=0; i<Nx+2*NGHX; i++) {
	
	r     = Ymed(j);
	omega = sqrt(G*MSTAR/r/r/r);                       //Keplerian frequency
	rhog  = SIGMA0*pow(r/R0,-SIGMASLOPE);              //Gas surface density
        rhod  = rhog*EPSILON;                              //Dust surface density

	if (Fluidtype == GAS) {
	  rho[l]   = rhog;
	  vphi[l]  = omega*r*sqrt(1.0 + pow(ASPECTRATIO,2)*pow(r/R0,2*FLARINGINDEX)*
				  (2.0*FLARINGINDEX - 1.0 - SIGMASLOPE));
	  vr[l]    = 0.0;
	  cs[l]    = ASPECTRATIO*pow(r/R0,FLARINGINDEX)*omega*r;
	}
	
	if (Fluidtype == DUST) {
	  rho[l]  = rhog*epsilons[id-1];
	  vphi[l] = omega*r;
	  vr[l]   = 0.0;
	  cs[l]   = 0.0;
	}
	
	vphi[l] -= OMEGAFRAME*r;
	
      }
    }
  }
}

void CondInit() {
  
  int id_gas = 0;
  int feedback = YES;
  //We first create the gaseous fluid and store it in the array Fluids[]
  Fluids[id_gas] = CreateFluid("gas",GAS);

  //We now select the fluid
  SelectFluid(id_gas);

  //and fill its fields
  _CondInit(0);

  //We repeat the process for the dust fluids
  char dust_name[MAXNAMELENGTH];
  int id_dust;

  for(id_dust = 1; id_dust<NFLUIDS; id_dust++) {
    sprintf(dust_name,"dust%d",id_dust); //We assign different names to the dust fluids

    Fluids[id_dust]  = CreateFluid(dust_name, DUST);
    SelectFluid(id_dust);
    _CondInit(id_dust);

  }

  #ifdef COLLISIONS
  /*We now fill the collision matrix (Feedback from dust included)
   Note: ColRate() moves the collision matrix to the device.
   If feedback=NO, gas does not feel the drag force.*/
  
  ColRate(INVSTOKES1, id_gas, 1, feedback);
  ColRate(INVSTOKES2, id_gas, 2, feedback);
  #endif
}
