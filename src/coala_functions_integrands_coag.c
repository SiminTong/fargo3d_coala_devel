//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

//-------------------------------------------
// integrand functions for coag flux
//-------------------------------------------

#include "coala_functions_integrands_coag.h"
#include "coala_collision_kernels.h"
#include "coala_polynomials_legendre.h"
#include "math.h"
#include <stdlib.h>
#include <stdio.h>

/**
 * @brief Evaluate the integrand of the coagulation flux, i.e. integrand of the double integral
 *
 * @param    kernel    select the collisional kernel function
 * @param    K0        constant value of the kernel function (used to adapt to code unit)
 * @param    i         order of polynomials
 * @param    ip        order of polynomials
 * @param    ai        coefficients of polynomial of degree i, sorted from low to high order
 * @param    aip       coefficients of polynomial of degree ip, sorted from low to high order
 * @param    u         mass variable (colliding grain of mass u)
 * @param    v         mass variable (colliding grain of mass v)
 * @param    xilp      variable mapping the mass bin lp in [-1,1], needed for Legendre polynomials
 * @param    xil       variable mapping the mass bin l in [-1,1], needed for Legendre polynomials
 * @return             integrand of cogulation flux evaluated at u,v,xilp,xil
 */
flt func_coag_flux(unsigned int kernel, flt K0,unsigned int i, unsigned int ip, flt* ai, flt* aip, flt u, flt v, flt xilp, flt xil){


   return func_kernel(kernel,K0,u,v)*phi_pol(ip,aip,xilp)*phi_pol(i,ai,xil)/v;
}

/**
 * @brief Evaluate the integrand of the term inclunding the integral of the coagulation flux, i.e. integrand of the triple integral
 *
 * @param    kernel    select the collisional kernel function
 * @param    K0        constant value of the kernel function (used to adapt to code unit)
 * @param    k         order of polynomials
 * @param    i         order of polynomials
 * @param    ip        order of polynomials
 * @param    ak        coefficients of polynomial of degree k, sorted from low to high order
 * @param    ai        coefficients of polynomial of degree i, sorted from low to high order
 * @param    aip       coefficients of polynomial of degree ip, sorted from low to high order
 * @param    hj        size of mass bin j
 * @param    u         mass variable (colliding grain of mass u)
 * @param    v         mass variable (colliding grain of mass v)
 * @param    xij       variable mapping the mass bin j in [-1,1], needed for Legendre polynomials
 * @param    xilp      variable mapping the mass bin lp in [-1,1], needed for Legendre polynomials
 * @param    xil       variable mapping the mass bin l in [-1,1], needed for Legendre polynomials
 * @return             integrand of the term including the integral of cogulation flux evaluated at u,v,xij,xilp,xil
 */
flt func_coag_intflux(unsigned int kernel, flt K0,unsigned int k, unsigned int i, unsigned int ip, flt* ak, flt* ai, flt* aip, flt u, flt v, flt hj, flt xij, flt xilp, flt xil){
   

   return dphi_pol_k(k,ak,hj,xij)*func_kernel(kernel,K0,u,v)*phi_pol(ip,aip,xilp)*phi_pol(i,ai,xil)/v;
}




