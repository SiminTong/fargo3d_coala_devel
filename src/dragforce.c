//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>


void DragForce(real dt) {
  
  
  // Centered drag coefficient pre-factor
  FARGO_SAFE(DragForce_Coeff()); //store dragcoeff in Qs
  //-------------------------------------------------------------------------------------
  Reset_field(DensStar);  
  MULTIFLUID(DragForce_SumC(dt)); 

#ifdef X
  Reset_field(Mpx);
  MULTIFLUID(DragForce_SumCV(dt,0));
#endif

#ifdef Y  
  Reset_field(Mpy);
  MULTIFLUID(DragForce_SumCV(dt,1)); 
#endif

#ifdef Z  
  Reset_field(Mpz); 
  MULTIFLUID(DragForce_SumCV(dt,2)); 
#endif
  
#ifdef X
  MULTIFLUID(DragForce_UpdateVel(dt,0)); //update velocities
#endif

#ifdef Y
  MULTIFLUID(DragForce_UpdateVel(dt,1)); //update velocities
#endif
  
#ifdef Z  
  MULTIFLUID(DragForce_UpdateVel(dt,2)); //update velocities
#endif
}
