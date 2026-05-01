//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************


//-------------------------------------------
// Collision kernel functions
//-------------------------------------------

#include "coala_collision_kernels.h"
#include "math.h"
#include <stdlib.h>
#include <stdio.h>

/**
 * @brief Compute constant kernel 
 *
 * @param    K0    constant value of the kernel function (used to adapt to code unit)
 * @param    u     mass variable (colliding grain of mass u)
 * @param    v     mass variable (colliding grain of mass v)
 * @return         evaluate constant kernel
 */
flt kconst(flt K0,flt u,flt v){
   return K0;
}

/**
 * @brief Compute additive kernel 
 *
 * @param    K0    constant value of the kernel function (used to adapt to code unit)
 * @param    u     mass variable (colliding grain of mass u)
 * @param    v     mass variable (colliding grain of mass v)
 * @return         evaluate constant kernel
 */
flt kadd(flt K0,flt u,flt v){
   return K0*(u+v);
   
}

/**
 * @brief Compute multiplicative kernel 
 *
 * @param    K0    constant value of the kernel function (used to adapt to code unit)
 * @param    u     mass variable (colliding grain of mass u)
 * @param    v     mass variable (colliding grain of mass v)
 * @return         evaluate constant kernel
 */
flt kmul(flt K0, flt u,flt v){
   return K0*u*v;
   
}

/**
 * @brief Compute the cross-section term in the ballistic kernel K = sigma * dv
 *
 * @param    K0    constant value of the kernel function (used to adapt to code unit)
 * @param    u     mass variable (colliding grain of mass u)
 * @param    v     mass variable (colliding grain of mass v)
 * @return         evaluate constant kernel
 */
flt cross_section(flt K0,flt u,flt v){
   // return K0*(pow(u,((flt)2)/((flt)3)) + ((flt)2) * pow(u,((flt)1)/((flt)3))*pow(v,((flt)1)/((flt)3)) + pow(v,((flt)2)/((flt)3)));
   return K0*pow(pow(u,((flt)1)/((flt)3)) + pow(v,((flt)1)/((flt)3)),((flt)2)) ;
   
}

/**
 * @brief Compute the cross-section term in the ballistic kernel K = sigma * dv
 *
 * @param    kernel    select the collisional kernel function
 * @param    K0        constant value of the kernel function (used to adapt to code unit)
 * @param    u         mass variable (colliding grain of mass u)
 * @param    v         mass variable (colliding grain of mass v)
 * @return   res       evaluate constant kernel
 */
flt func_kernel(unsigned int kernel, flt K0,flt u, flt v){
   flt res;
   switch (kernel){
      case 0:
         res = kconst(K0,u,v);
         break;
      case 1:
         res = kadd(K0,u,v);
         break;
      case 2:
         res = kmul(K0,u,v);
         break;
      case 3:
         res = cross_section(K0,u,v);
         break;
      default :
         printf(" only simple kernels and cross section for physical kernel\n");
         exit(-1);
   }

   return res;
   
}