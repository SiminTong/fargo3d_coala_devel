//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

#include "coala_precision.h"

flt coeffnorm(unsigned int i);
void compute_coeff_Leg(unsigned int k,flt* coeff_Leg, unsigned int size_coeff_leg);
void compute_mat_coeffs(unsigned int k,flt* mat_coeffs);
flt phi_pol(unsigned int i,flt* ai,flt x);
void polynomial_derivative_coeffs(unsigned int k, flt* pol_coeffs,flt* dpol_coeffs);
flt dphi_pol_k(unsigned int k,flt* ak,flt hj,flt xij);
