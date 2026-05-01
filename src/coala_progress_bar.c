//***********************************************************************************
// Coala code
// Copyright(C) Maxime Lombart <maxime.lombart@cea.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see LICENCE for more information
//***********************************************************************************

#include <stdio.h>

void progressbar(unsigned int nstep, unsigned int iprogress){
    //progress bar
    int dim_bar = 39;
	char bar[39] = {'[','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_','_',']'};

    // printf("len bar = %d\n ",strlen(bar));
    // printf("%c\n",bar[0]);
    // printf("%c\n",bar[38]);
    



	if (nstep==1){
		printf("progress -> 100%% \n");
	}else{

        // for (unsigned int j=1;j<iprogress*(strlen(bar)-2)/(nstep-1); j++){
		for (unsigned int j=1; j<dim_bar-1;j++){
			if (j<iprogress*(dim_bar-1)/(nstep-1)){
				bar[j] = '#';
	            // printf("j=%d\n",j);
	            // printf("%s \n",bar);
			}
            
            
        }

        // printf("%d\n",strlen(bar));
        // exit(-1);
		
		// printf("\r%s \n",bar);
		printf("\r%.39s  %3.f%%",bar,100.0*iprogress/(nstep-1));
		fflush(stdout);



	}

}