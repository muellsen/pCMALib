

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include "bbobStructures.h" /* Include all declarations for BBOB calls */


void init_bbob (int *dim,int *fun,char *output_path)
{
/* declare your own internal stuff */

  
    unsigned int instances[15] = {1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5};
    unsigned int idx_dim, ifun, idx_instances, seed;

    
    /**************************************************
     *          BBOB Mandatory initialization         *
     *************************************************/
    /* retrieve all default parameters of BBOB calls  */
    ParamStruct params = fgeneric_getDefaultPARAMS();

    /* modify the following parameters, choosing a different setting
     * for each new experiment */
    strcpy(params.dataPath, "PUT_MY_BBOB_DATA_PATH");  /* different folder for each experiment! */
    /* please beforehand run from the command-line 'python createfolders.py PUT_MY_BBOB_DATA_PATH'
     * to create the necessary folder structure to run an experiment. */
    strcpy(params.algName, "PUT ALGORITHM NAME");
    strcpy(params.comments, "PUT MORE DETAILED INFORMATION, PARAMETER SETTINGS ETC");

    seed = time(NULL);
    srand(seed); /* used by MY_OPTIMIZER */
 //   printf("MY_OPTIMIZER seed: %d\n", seed);

    /* To make the noise deterministic. */
    /* fgeneric_noiseseed(30); printf("seed for the noise set to: 30\n"); */


		/* set DIM, funcId, instanceId to initialize BBOB fgeneric */
		params.DIM = *dim;
		params.funcId = *fun;
		params.instanceId = 1;
                /* call the BBOB initialization */
        fgeneric_initialize_(params);
}
