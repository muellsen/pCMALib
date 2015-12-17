#include "stdafx.h"
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include "bbobStructures.h" /* Include all declarations for BBOB calls */


void BBOB_init()
{
	unsigned int dim[6] = {2, 3, 5, 10, 20, 40};
    unsigned int instances[15] = {1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5};
    unsigned int idx_dim, ifun, idx_instances, seed;    
    
    /**************************************************
     *          BBOB Mandatory initialization         *
     *************************************************/
    /* retrieve all default parameters of BBOB calls  */
    ParamStruct params = fgeneric_getDefaultPARAMS();
    /* modify the following parameters, choosing a different setting
     * for each new experiment */
    strcpy(params.dataPath, "test");  /* different folder for each experiment! */
    /* please beforehand run from the command-line 'python createfolders.py PUT_MY_BBOB_DATA_PATH'
     * to create the necessary folder structure to run an experiment. */
    strcpy(params.algName, " NAME");
    strcpy(params.comments, "PUT MORE DETAILED INFORMATION, PARAMETER SETTINGS ETC");

        /* To make the noise deterministic. */
    /* fgeneric_noiseseed(30); printf("seed for the noise set to: 30\n"); */

    
    params.DIM = dim[idx_dim];
    params.funcId = ifun;
    params.instanceId = instances[idx_instances];
    /* call the BBOB initialization */
    fgeneric_initialize(params);
}
