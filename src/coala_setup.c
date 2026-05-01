//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

//-------------------------------------------
//setup file to choose parameters
//-------------------------------------------
#include "coala_precision.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h> 

//-------------------------------------------
// kernel    -> choose among simple kernels
// K0        -> normalisation coefficient for simple kernels
// nbins     -> number of grain size bins
// kpol      -> order of polynomials
// Q         -> number of Gauss points for Gauss-quadrature method
// eps       -> minimum value for mass density distribution
// coeff_CFL -> coefficient for SSPRK 3 time solver
// dthydro   -> hydro timestep 
// ndthydro  -> number of dthydro
// massmax   -> maximal mass of grains to consider
// minmass   -> minimal mass of grains to consider
//-------------------------------------------

//-------------------------------------------
//choose kernel
//0->kconst
//1->kadd
//-------------------------------------------
unsigned int kernel=0;
flt K0 = ((flt)1);

unsigned int nbins  = 20;
unsigned int kpol   = 0;
unsigned int Q      = 10;
flt eps             = ((flt)1e-30);
flt coeff_CFL       = ((flt)3e-1);

//default setup simple kernels (value defined in main program coala.c)
flt dthydro;
unsigned int ndthydro;

flt massmin = ((flt)1e-3);
flt massmax = ((flt)1e6);



//-------------------------------------------
//Save gij, massgrid in data/ directory
//-------------------------------------------
bool save = true;


//-------------------------------------------
//path variables for data/
//-------------------------------------------
char path_data[100];
char path_massgrid[100];
char path_massbins[100];
char path_log[100];
char path_gtend_massmeanlog[100];
char path_massmeanlog[100];
char path_gt0_massmeanlog[100];
char path_gij[100];
char path_gij_t0[100];
char path_gij_tend[100];
char path_time[100];




//-------------------------------------------
//function to compute variable for path and create data. directory
//-------------------------------------------
void init_path_files(unsigned int nbins, unsigned int kpol, unsigned int kernel,char path_data[]){

    char strkernel[15];

    switch (kernel){
        case 0:
            strcpy(strkernel,"kconst");
            break;
        case 1:
            strcpy(strkernel,"kadd");
            break;
        default:
            printf("issue in choosing kernel, setup \n");
            exit(-1);
    }


    char strnbins[5];
    char strkpol[2]; 
    sprintf(strnbins, "%d", nbins);
    sprintf(strkpol, "%d", kpol);

    // printf("nbins =%s \n",strnbins);


    strcpy(path_data,"../data/");
    strcat(path_data,strkernel);
    strcat(path_data,"/nbins=");
    strcat(path_data,strnbins);
    strcat(path_data,"/kpol=");
    strcat(path_data,strkpol);
    strcat(path_data,"/");
    // printf("path_data=%s \n",path_data);
    char command[200];
    strcpy(command,"[ ! -f ");
    strcat(command,path_data);
    strcat(command," ] && mkdir -p ");
    strcat(command,path_data);

    // printf("command =%s \n",command);
    system(command);


}









