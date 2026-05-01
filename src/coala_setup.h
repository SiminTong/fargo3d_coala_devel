//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

#include "coala_precision.h"
#include <stdbool.h> 
#include <math.h>

extern unsigned int kernel;
extern flt K0;
extern unsigned int nbins;
extern unsigned int kpol;
extern unsigned int Q;
extern flt eps;
extern flt coeff_CFL;
extern flt dthydro;
extern unsigned int ndthydro;
extern flt massmin;
extern flt massmax;

extern bool save;

extern char path_data[100];
extern char path_massgrid[100];
extern char path_massbins[100];
extern char path_log[100];
extern char path_gtend_massmeanlog[100];
extern char path_massmeanlog[100];
extern char path_gt0_massmeanlog[100];
extern char path_gij[100];
extern char path_gij_t0[100];
extern char path_gij_tend[100];
extern char path_time[100];


void init_path_files(unsigned int nbins, unsigned int kpol, unsigned int kernel,char path_data[]);

