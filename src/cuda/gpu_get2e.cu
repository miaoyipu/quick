/*
 *  gpu_get2e.cpp
 *  new_quick
 *
 *  Created by Yipu Miao on 6/17/11.
 *  Copyright 2011 University of Florida.All rights reserved.
 *
 */

#include "gpu.h"


static 
#ifndef TEST
__constant__
#endif
 gpu_simulation_type devSim;

static
#ifndef TEST
__constant__
#endif
int devTrans[TRANSDIM*TRANSDIM*TRANSDIM];

static
#ifndef TEST
__constant__
#endif
int devMcal[MCALDIM*3];

static
#ifndef TEST
__constant__
#endif
int Sumindex[10]={0,0,1,4,10,20,35,56,84,120};


void upload_sim_to_constant(_gpu_type gpu){
    cudaError_t status;
#ifdef TEST
    memcpy(&devSim, &gpu->gpu_sim, sizeof(gpu_simulation_type));
#else    
    status = cudaMemcpyToSymbol(devSim, &gpu->gpu_sim, sizeof(gpu_simulation_type));
    PRINTERROR(status, " cudaMemcpyToSymbol, copy to constants failed")
#endif
}


void upload_para_to_const(){

    int trans[TRANSDIM*TRANSDIM*TRANSDIM] = {0};
    // Data to trans
    {
        LOC3(trans, 0, 0, 0, TRANSDIM, TRANSDIM, TRANSDIM) =   1;
        LOC3(trans, 0, 0, 1, TRANSDIM, TRANSDIM, TRANSDIM) =   4;
        LOC3(trans, 0, 0, 2, TRANSDIM, TRANSDIM, TRANSDIM) =  10;
        LOC3(trans, 0, 0, 3, TRANSDIM, TRANSDIM, TRANSDIM) =  20;
        LOC3(trans, 0, 0, 4, TRANSDIM, TRANSDIM, TRANSDIM) =  35;
        LOC3(trans, 0, 0, 5, TRANSDIM, TRANSDIM, TRANSDIM) =  56;
        LOC3(trans, 0, 0, 6, TRANSDIM, TRANSDIM, TRANSDIM) =  84;
        LOC3(trans, 0, 0, 7, TRANSDIM, TRANSDIM, TRANSDIM) = 120;
        LOC3(trans, 0, 1, 0, TRANSDIM, TRANSDIM, TRANSDIM) =   3;
        LOC3(trans, 0, 1, 1, TRANSDIM, TRANSDIM, TRANSDIM) =   6;
        LOC3(trans, 0, 1, 2, TRANSDIM, TRANSDIM, TRANSDIM) =  17;
        LOC3(trans, 0, 1, 3, TRANSDIM, TRANSDIM, TRANSDIM) =  32;
        LOC3(trans, 0, 1, 4, TRANSDIM, TRANSDIM, TRANSDIM) =  48;
        LOC3(trans, 0, 1, 5, TRANSDIM, TRANSDIM, TRANSDIM) =  67;
        LOC3(trans, 0, 1, 6, TRANSDIM, TRANSDIM, TRANSDIM) = 100;
        LOC3(trans, 0, 2, 0, TRANSDIM, TRANSDIM, TRANSDIM) =   9;
        LOC3(trans, 0, 2, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  16;
        LOC3(trans, 0, 2, 2, TRANSDIM, TRANSDIM, TRANSDIM) =  23;
        LOC3(trans, 0, 2, 3, TRANSDIM, TRANSDIM, TRANSDIM) =  42;
        LOC3(trans, 0, 2, 4, TRANSDIM, TRANSDIM, TRANSDIM) =  73;
        LOC3(trans, 0, 2, 5, TRANSDIM, TRANSDIM, TRANSDIM) = 106;
        LOC3(trans, 0, 3, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  19;
        LOC3(trans, 0, 3, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  31;
        LOC3(trans, 0, 3, 2, TRANSDIM, TRANSDIM, TRANSDIM) =  43;
        LOC3(trans, 0, 3, 3, TRANSDIM, TRANSDIM, TRANSDIM) =  79;
        LOC3(trans, 0, 3, 4, TRANSDIM, TRANSDIM, TRANSDIM) = 112;
        LOC3(trans, 0, 4, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  34;
        LOC3(trans, 0, 4, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  49;
        LOC3(trans, 0, 4, 2, TRANSDIM, TRANSDIM, TRANSDIM) =  74;
        LOC3(trans, 0, 4, 3, TRANSDIM, TRANSDIM, TRANSDIM) = 113;
        LOC3(trans, 0, 5, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  55;
        LOC3(trans, 0, 5, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  68;
        LOC3(trans, 0, 5, 2, TRANSDIM, TRANSDIM, TRANSDIM) = 107;
        LOC3(trans, 0, 6, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  83;
        LOC3(trans, 0, 6, 1, TRANSDIM, TRANSDIM, TRANSDIM) = 101;
        LOC3(trans, 0, 7, 0, TRANSDIM, TRANSDIM, TRANSDIM) = 119;
        LOC3(trans, 1, 0, 0, TRANSDIM, TRANSDIM, TRANSDIM) =   2;
        LOC3(trans, 1, 0, 1, TRANSDIM, TRANSDIM, TRANSDIM) =   7;
        LOC3(trans, 1, 0, 2, TRANSDIM, TRANSDIM, TRANSDIM) =  15;
        LOC3(trans, 1, 0, 3, TRANSDIM, TRANSDIM, TRANSDIM) =  28;
        LOC3(trans, 1, 0, 4, TRANSDIM, TRANSDIM, TRANSDIM) =  50;
        LOC3(trans, 1, 0, 5, TRANSDIM, TRANSDIM, TRANSDIM) =  69;
        LOC3(trans, 1, 0, 6, TRANSDIM, TRANSDIM, TRANSDIM) = 102;
        LOC3(trans, 1, 1, 0, TRANSDIM, TRANSDIM, TRANSDIM) =   5;
        LOC3(trans, 1, 1, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  11;
        LOC3(trans, 1, 1, 2, TRANSDIM, TRANSDIM, TRANSDIM) =  26;
        LOC3(trans, 1, 1, 3, TRANSDIM, TRANSDIM, TRANSDIM) =  41;
        LOC3(trans, 1, 1, 4, TRANSDIM, TRANSDIM, TRANSDIM) =  59;
        LOC3(trans, 1, 1, 5, TRANSDIM, TRANSDIM, TRANSDIM) =  87;
        LOC3(trans, 1, 2, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  13;
        LOC3(trans, 1, 2, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  25;
        LOC3(trans, 1, 2, 2, TRANSDIM, TRANSDIM, TRANSDIM) =  36;
        LOC3(trans, 1, 2, 3, TRANSDIM, TRANSDIM, TRANSDIM) =  60;
        LOC3(trans, 1, 2, 4, TRANSDIM, TRANSDIM, TRANSDIM) =  88;
        LOC3(trans, 1, 3, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  30;
        LOC3(trans, 1, 3, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  40;
        LOC3(trans, 1, 3, 2, TRANSDIM, TRANSDIM, TRANSDIM) =  61;
        LOC3(trans, 1, 3, 3, TRANSDIM, TRANSDIM, TRANSDIM) =  94;
        LOC3(trans, 1, 4, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  52;
        LOC3(trans, 1, 4, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  58;
        LOC3(trans, 1, 4, 2, TRANSDIM, TRANSDIM, TRANSDIM) =  89;
        LOC3(trans, 1, 5, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  71;
        LOC3(trans, 1, 5, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  86;
        LOC3(trans, 1, 6, 0, TRANSDIM, TRANSDIM, TRANSDIM) = 104;
        LOC3(trans, 2, 0, 0, TRANSDIM, TRANSDIM, TRANSDIM) =   8;
        LOC3(trans, 2, 0, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  14;
        LOC3(trans, 2, 0, 2, TRANSDIM, TRANSDIM, TRANSDIM) =  22;
        LOC3(trans, 2, 0, 3, TRANSDIM, TRANSDIM, TRANSDIM) =  44;
        LOC3(trans, 2, 0, 4, TRANSDIM, TRANSDIM, TRANSDIM) =  75;
        LOC3(trans, 2, 0, 5, TRANSDIM, TRANSDIM, TRANSDIM) = 108;
        LOC3(trans, 2, 1, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  12;
        LOC3(trans, 2, 1, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  24;
        LOC3(trans, 2, 1, 2, TRANSDIM, TRANSDIM, TRANSDIM) =  37;
        LOC3(trans, 2, 1, 3, TRANSDIM, TRANSDIM, TRANSDIM) =  62;
        LOC3(trans, 2, 1, 4, TRANSDIM, TRANSDIM, TRANSDIM) =  90;
        LOC3(trans, 2, 2, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  21;
        LOC3(trans, 2, 2, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  38;
        LOC3(trans, 2, 2, 2, TRANSDIM, TRANSDIM, TRANSDIM) =  66;
        LOC3(trans, 2, 2, 3, TRANSDIM, TRANSDIM, TRANSDIM) =  99;
        LOC3(trans, 2, 3, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  46;
        LOC3(trans, 2, 3, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  64;
        LOC3(trans, 2, 3, 2, TRANSDIM, TRANSDIM, TRANSDIM) =  98;
        LOC3(trans, 2, 4, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  77;
        LOC3(trans, 2, 4, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  92;
        LOC3(trans, 2, 5, 0, TRANSDIM, TRANSDIM, TRANSDIM) = 110;
        LOC3(trans, 3, 0, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  18;
        LOC3(trans, 3, 0, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  27;
        LOC3(trans, 3, 0, 2, TRANSDIM, TRANSDIM, TRANSDIM) =  45;
        LOC3(trans, 3, 0, 3, TRANSDIM, TRANSDIM, TRANSDIM) =  80;
        LOC3(trans, 3, 0, 4, TRANSDIM, TRANSDIM, TRANSDIM) = 114;
        LOC3(trans, 3, 1, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  29;
        LOC3(trans, 3, 1, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  39;
        LOC3(trans, 3, 1, 2, TRANSDIM, TRANSDIM, TRANSDIM) =  63;
        LOC3(trans, 3, 1, 3, TRANSDIM, TRANSDIM, TRANSDIM) =  95;
        LOC3(trans, 3, 2, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  47;
        LOC3(trans, 3, 2, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  65;
        LOC3(trans, 3, 2, 2, TRANSDIM, TRANSDIM, TRANSDIM) =  97;
        LOC3(trans, 3, 3, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  81;
        LOC3(trans, 3, 3, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  96;
        LOC3(trans, 3, 4, 0, TRANSDIM, TRANSDIM, TRANSDIM) = 116;
        LOC3(trans, 4, 0, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  33;
        LOC3(trans, 4, 0, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  51;
        LOC3(trans, 4, 0, 2, TRANSDIM, TRANSDIM, TRANSDIM) =  76;
        LOC3(trans, 4, 0, 3, TRANSDIM, TRANSDIM, TRANSDIM) = 115;
        LOC3(trans, 4, 1, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  53;
        LOC3(trans, 4, 1, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  57;
        LOC3(trans, 4, 1, 2, TRANSDIM, TRANSDIM, TRANSDIM) =  91;
        LOC3(trans, 4, 2, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  78;
        LOC3(trans, 4, 2, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  93;
        LOC3(trans, 4, 3, 0, TRANSDIM, TRANSDIM, TRANSDIM) = 117;
        LOC3(trans, 5, 0, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  54;
        LOC3(trans, 5, 0, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  70;
        LOC3(trans, 5, 0, 2, TRANSDIM, TRANSDIM, TRANSDIM) = 109;
        LOC3(trans, 5, 1, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  72;
        LOC3(trans, 5, 1, 1, TRANSDIM, TRANSDIM, TRANSDIM) =  85;
        LOC3(trans, 5, 2, 0, TRANSDIM, TRANSDIM, TRANSDIM) = 111;
        LOC3(trans, 6, 0, 0, TRANSDIM, TRANSDIM, TRANSDIM) =  82;
        LOC3(trans, 6, 0, 1, TRANSDIM, TRANSDIM, TRANSDIM) = 103;
        LOC3(trans, 6, 1, 0, TRANSDIM, TRANSDIM, TRANSDIM) = 105;
        LOC3(trans, 7, 0, 0, TRANSDIM, TRANSDIM, TRANSDIM) = 118;
    }
    
    int Mcal[3*MCALDIM] = {0};
    {
        LOC2(Mcal, 0,   0, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,   1, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,   2, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,   3, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,   4, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,   5, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,   6, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,   7, 3, MCALDIM) =   2;
        LOC2(Mcal, 0,   8, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,   9, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  10, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  11, 3, MCALDIM) =   2;
        LOC2(Mcal, 0,  12, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  13, 3, MCALDIM) =   2;
        LOC2(Mcal, 0,  14, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  15, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  16, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  17, 3, MCALDIM) =   3;
        LOC2(Mcal, 0,  18, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  19, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  20, 3, MCALDIM) =   2;
        LOC2(Mcal, 0,  21, 3, MCALDIM) =   2;
        LOC2(Mcal, 0,  22, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  23, 3, MCALDIM) =   2;
        LOC2(Mcal, 0,  24, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  25, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  26, 3, MCALDIM) =   3;
        LOC2(Mcal, 0,  27, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  28, 3, MCALDIM) =   3;
        LOC2(Mcal, 0,  29, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  30, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  31, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  32, 3, MCALDIM) =   4;
        LOC2(Mcal, 0,  33, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  34, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  35, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  36, 3, MCALDIM) =   2;
        LOC2(Mcal, 0,  37, 3, MCALDIM) =   2;
        LOC2(Mcal, 0,  38, 3, MCALDIM) =   3;
        LOC2(Mcal, 0,  39, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  40, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  41, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  42, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  43, 3, MCALDIM) =   2;
        LOC2(Mcal, 0,  44, 3, MCALDIM) =   3;
        LOC2(Mcal, 0,  45, 3, MCALDIM) =   2;
        LOC2(Mcal, 0,  46, 3, MCALDIM) =   3;
        LOC2(Mcal, 0,  47, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  48, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  49, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  50, 3, MCALDIM) =   4;
        LOC2(Mcal, 0,  51, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  52, 3, MCALDIM) =   4;
        LOC2(Mcal, 0,  53, 3, MCALDIM) =   5;
        LOC2(Mcal, 0,  54, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  55, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  56, 3, MCALDIM) =   4;
        LOC2(Mcal, 0,  57, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  58, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  59, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  60, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  61, 3, MCALDIM) =   2;
        LOC2(Mcal, 0,  62, 3, MCALDIM) =   3;
        LOC2(Mcal, 0,  63, 3, MCALDIM) =   2;
        LOC2(Mcal, 0,  64, 3, MCALDIM) =   3;
        LOC2(Mcal, 0,  65, 3, MCALDIM) =   2;
        LOC2(Mcal, 0,  66, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  67, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  68, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  69, 3, MCALDIM) =   5;
        LOC2(Mcal, 0,  70, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  71, 3, MCALDIM) =   5;
        LOC2(Mcal, 0,  72, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  73, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  74, 3, MCALDIM) =   2;
        LOC2(Mcal, 0,  75, 3, MCALDIM) =   4;
        LOC2(Mcal, 0,  76, 3, MCALDIM) =   2;
        LOC2(Mcal, 0,  77, 3, MCALDIM) =   4;
        LOC2(Mcal, 0,  78, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  79, 3, MCALDIM) =   3;
        LOC2(Mcal, 0,  80, 3, MCALDIM) =   3;
        LOC2(Mcal, 0,  81, 3, MCALDIM) =   6;
        LOC2(Mcal, 0,  82, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  83, 3, MCALDIM) =   0;
        LOC2(Mcal, 0,  84, 3, MCALDIM) =   5;
        LOC2(Mcal, 0,  85, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  86, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  87, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  88, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  89, 3, MCALDIM) =   2;
        LOC2(Mcal, 0,  90, 3, MCALDIM) =   4;
        LOC2(Mcal, 0,  91, 3, MCALDIM) =   2;
        LOC2(Mcal, 0,  92, 3, MCALDIM) =   4;
        LOC2(Mcal, 0,  93, 3, MCALDIM) =   1;
        LOC2(Mcal, 0,  94, 3, MCALDIM) =   3;
        LOC2(Mcal, 0,  95, 3, MCALDIM) =   3;
        LOC2(Mcal, 0,  96, 3, MCALDIM) =   3;
        LOC2(Mcal, 0,  97, 3, MCALDIM) =   2;
        LOC2(Mcal, 0,  98, 3, MCALDIM) =   2;
        LOC2(Mcal, 0,  99, 3, MCALDIM) =   0;
        LOC2(Mcal, 0, 100, 3, MCALDIM) =   0;
        LOC2(Mcal, 0, 101, 3, MCALDIM) =   1;
        LOC2(Mcal, 0, 102, 3, MCALDIM) =   6;
        LOC2(Mcal, 0, 103, 3, MCALDIM) =   1;
        LOC2(Mcal, 0, 104, 3, MCALDIM) =   6;
        LOC2(Mcal, 0, 105, 3, MCALDIM) =   0;
        LOC2(Mcal, 0, 106, 3, MCALDIM) =   0;
        LOC2(Mcal, 0, 107, 3, MCALDIM) =   2;
        LOC2(Mcal, 0, 108, 3, MCALDIM) =   5;
        LOC2(Mcal, 0, 109, 3, MCALDIM) =   2;
        LOC2(Mcal, 0, 110, 3, MCALDIM) =   5;
        LOC2(Mcal, 0, 111, 3, MCALDIM) =   0;
        LOC2(Mcal, 0, 112, 3, MCALDIM) =   0;
        LOC2(Mcal, 0, 113, 3, MCALDIM) =   3;
        LOC2(Mcal, 0, 114, 3, MCALDIM) =   4;
        LOC2(Mcal, 0, 115, 3, MCALDIM) =   3;
        LOC2(Mcal, 0, 116, 3, MCALDIM) =   4;
        LOC2(Mcal, 0, 117, 3, MCALDIM) =   7;
        LOC2(Mcal, 0, 118, 3, MCALDIM) =   0;
        LOC2(Mcal, 0, 119, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,   0, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,   1, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,   2, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,   3, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,   4, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,   5, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,   6, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,   7, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,   8, 3, MCALDIM) =   2;
        LOC2(Mcal, 1,   9, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  10, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  11, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  12, 3, MCALDIM) =   2;
        LOC2(Mcal, 1,  13, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  14, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  15, 3, MCALDIM) =   2;
        LOC2(Mcal, 1,  16, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  17, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  18, 3, MCALDIM) =   3;
        LOC2(Mcal, 1,  19, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  20, 3, MCALDIM) =   2;
        LOC2(Mcal, 1,  21, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  22, 3, MCALDIM) =   2;
        LOC2(Mcal, 1,  23, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  24, 3, MCALDIM) =   2;
        LOC2(Mcal, 1,  25, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  26, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  27, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  28, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  29, 3, MCALDIM) =   3;
        LOC2(Mcal, 1,  30, 3, MCALDIM) =   3;
        LOC2(Mcal, 1,  31, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  32, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  33, 3, MCALDIM) =   4;
        LOC2(Mcal, 1,  34, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  35, 3, MCALDIM) =   2;
        LOC2(Mcal, 1,  36, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  37, 3, MCALDIM) =   2;
        LOC2(Mcal, 1,  38, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  39, 3, MCALDIM) =   3;
        LOC2(Mcal, 1,  40, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  41, 3, MCALDIM) =   2;
        LOC2(Mcal, 1,  42, 3, MCALDIM) =   3;
        LOC2(Mcal, 1,  43, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  44, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  45, 3, MCALDIM) =   3;
        LOC2(Mcal, 1,  46, 3, MCALDIM) =   2;
        LOC2(Mcal, 1,  47, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  48, 3, MCALDIM) =   4;
        LOC2(Mcal, 1,  49, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  50, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  51, 3, MCALDIM) =   4;
        LOC2(Mcal, 1,  52, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  53, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  54, 3, MCALDIM) =   5;
        LOC2(Mcal, 1,  55, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  56, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  57, 3, MCALDIM) =   4;
        LOC2(Mcal, 1,  58, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  59, 3, MCALDIM) =   2;
        LOC2(Mcal, 1,  60, 3, MCALDIM) =   3;
        LOC2(Mcal, 1,  61, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  62, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  63, 3, MCALDIM) =   3;
        LOC2(Mcal, 1,  64, 3, MCALDIM) =   2;
        LOC2(Mcal, 1,  65, 3, MCALDIM) =   2;
        LOC2(Mcal, 1,  66, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  67, 3, MCALDIM) =   5;
        LOC2(Mcal, 1,  68, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  69, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  70, 3, MCALDIM) =   5;
        LOC2(Mcal, 1,  71, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  72, 3, MCALDIM) =   2;
        LOC2(Mcal, 1,  73, 3, MCALDIM) =   4;
        LOC2(Mcal, 1,  74, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  75, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  76, 3, MCALDIM) =   4;
        LOC2(Mcal, 1,  77, 3, MCALDIM) =   2;
        LOC2(Mcal, 1,  78, 3, MCALDIM) =   3;
        LOC2(Mcal, 1,  79, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  80, 3, MCALDIM) =   3;
        LOC2(Mcal, 1,  81, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  82, 3, MCALDIM) =   6;
        LOC2(Mcal, 1,  83, 3, MCALDIM) =   0;
        LOC2(Mcal, 1,  84, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  85, 3, MCALDIM) =   5;
        LOC2(Mcal, 1,  86, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  87, 3, MCALDIM) =   2;
        LOC2(Mcal, 1,  88, 3, MCALDIM) =   4;
        LOC2(Mcal, 1,  89, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  90, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  91, 3, MCALDIM) =   4;
        LOC2(Mcal, 1,  92, 3, MCALDIM) =   2;
        LOC2(Mcal, 1,  93, 3, MCALDIM) =   3;
        LOC2(Mcal, 1,  94, 3, MCALDIM) =   1;
        LOC2(Mcal, 1,  95, 3, MCALDIM) =   3;
        LOC2(Mcal, 1,  96, 3, MCALDIM) =   2;
        LOC2(Mcal, 1,  97, 3, MCALDIM) =   3;
        LOC2(Mcal, 1,  98, 3, MCALDIM) =   2;
        LOC2(Mcal, 1,  99, 3, MCALDIM) =   1;
        LOC2(Mcal, 1, 100, 3, MCALDIM) =   6;
        LOC2(Mcal, 1, 101, 3, MCALDIM) =   0;
        LOC2(Mcal, 1, 102, 3, MCALDIM) =   0;
        LOC2(Mcal, 1, 103, 3, MCALDIM) =   6;
        LOC2(Mcal, 1, 104, 3, MCALDIM) =   1;
        LOC2(Mcal, 1, 105, 3, MCALDIM) =   2;
        LOC2(Mcal, 1, 106, 3, MCALDIM) =   5;
        LOC2(Mcal, 1, 107, 3, MCALDIM) =   0;
        LOC2(Mcal, 1, 108, 3, MCALDIM) =   0;
        LOC2(Mcal, 1, 109, 3, MCALDIM) =   5;
        LOC2(Mcal, 1, 110, 3, MCALDIM) =   2;
        LOC2(Mcal, 1, 111, 3, MCALDIM) =   3;
        LOC2(Mcal, 1, 112, 3, MCALDIM) =   4;
        LOC2(Mcal, 1, 113, 3, MCALDIM) =   0;
        LOC2(Mcal, 1, 114, 3, MCALDIM) =   0;
        LOC2(Mcal, 1, 115, 3, MCALDIM) =   4;
        LOC2(Mcal, 1, 116, 3, MCALDIM) =   3;
        LOC2(Mcal, 1, 117, 3, MCALDIM) =   0;
        LOC2(Mcal, 1, 118, 3, MCALDIM) =   7;
        LOC2(Mcal, 1, 119, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,   0, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,   1, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,   2, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,   3, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,   4, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,   5, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,   6, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,   7, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,   8, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,   9, 3, MCALDIM) =   2;
        LOC2(Mcal, 2,  10, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  11, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  12, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  13, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  14, 3, MCALDIM) =   2;
        LOC2(Mcal, 2,  15, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  16, 3, MCALDIM) =   2;
        LOC2(Mcal, 2,  17, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  18, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  19, 3, MCALDIM) =   3;
        LOC2(Mcal, 2,  20, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  21, 3, MCALDIM) =   2;
        LOC2(Mcal, 2,  22, 3, MCALDIM) =   2;
        LOC2(Mcal, 2,  23, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  24, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  25, 3, MCALDIM) =   2;
        LOC2(Mcal, 2,  26, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  27, 3, MCALDIM) =   3;
        LOC2(Mcal, 2,  28, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  29, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  30, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  31, 3, MCALDIM) =   3;
        LOC2(Mcal, 2,  32, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  33, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  34, 3, MCALDIM) =   4;
        LOC2(Mcal, 2,  35, 3, MCALDIM) =   2;
        LOC2(Mcal, 2,  36, 3, MCALDIM) =   2;
        LOC2(Mcal, 2,  37, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  38, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  39, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  40, 3, MCALDIM) =   3;
        LOC2(Mcal, 2,  41, 3, MCALDIM) =   3;
        LOC2(Mcal, 2,  42, 3, MCALDIM) =   2;
        LOC2(Mcal, 2,  43, 3, MCALDIM) =   3;
        LOC2(Mcal, 2,  44, 3, MCALDIM) =   2;
        LOC2(Mcal, 2,  45, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  46, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  47, 3, MCALDIM) =   4;
        LOC2(Mcal, 2,  48, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  49, 3, MCALDIM) =   4;
        LOC2(Mcal, 2,  50, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  51, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  52, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  53, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  54, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  55, 3, MCALDIM) =   5;
        LOC2(Mcal, 2,  56, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  57, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  58, 3, MCALDIM) =   4;
        LOC2(Mcal, 2,  59, 3, MCALDIM) =   3;
        LOC2(Mcal, 2,  60, 3, MCALDIM) =   2;
        LOC2(Mcal, 2,  61, 3, MCALDIM) =   3;
        LOC2(Mcal, 2,  62, 3, MCALDIM) =   2;
        LOC2(Mcal, 2,  63, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  64, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  65, 3, MCALDIM) =   2;
        LOC2(Mcal, 2,  66, 3, MCALDIM) =   5;
        LOC2(Mcal, 2,  67, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  68, 3, MCALDIM) =   5;
        LOC2(Mcal, 2,  69, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  70, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  71, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  72, 3, MCALDIM) =   4;
        LOC2(Mcal, 2,  73, 3, MCALDIM) =   2;
        LOC2(Mcal, 2,  74, 3, MCALDIM) =   4;
        LOC2(Mcal, 2,  75, 3, MCALDIM) =   2;
        LOC2(Mcal, 2,  76, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  77, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  78, 3, MCALDIM) =   3;
        LOC2(Mcal, 2,  79, 3, MCALDIM) =   3;
        LOC2(Mcal, 2,  80, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  81, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  82, 3, MCALDIM) =   0;
        LOC2(Mcal, 2,  83, 3, MCALDIM) =   6;
        LOC2(Mcal, 2,  84, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  85, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  86, 3, MCALDIM) =   5;
        LOC2(Mcal, 2,  87, 3, MCALDIM) =   4;
        LOC2(Mcal, 2,  88, 3, MCALDIM) =   2;
        LOC2(Mcal, 2,  89, 3, MCALDIM) =   4;
        LOC2(Mcal, 2,  90, 3, MCALDIM) =   2;
        LOC2(Mcal, 2,  91, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  92, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  93, 3, MCALDIM) =   3;
        LOC2(Mcal, 2,  94, 3, MCALDIM) =   3;
        LOC2(Mcal, 2,  95, 3, MCALDIM) =   1;
        LOC2(Mcal, 2,  96, 3, MCALDIM) =   2;
        LOC2(Mcal, 2,  97, 3, MCALDIM) =   2;
        LOC2(Mcal, 2,  98, 3, MCALDIM) =   3;
        LOC2(Mcal, 2,  99, 3, MCALDIM) =   6;
        LOC2(Mcal, 2, 100, 3, MCALDIM) =   1;
        LOC2(Mcal, 2, 101, 3, MCALDIM) =   6;
        LOC2(Mcal, 2, 102, 3, MCALDIM) =   1;
        LOC2(Mcal, 2, 103, 3, MCALDIM) =   0;
        LOC2(Mcal, 2, 104, 3, MCALDIM) =   0;
        LOC2(Mcal, 2, 105, 3, MCALDIM) =   5;
        LOC2(Mcal, 2, 106, 3, MCALDIM) =   2;
        LOC2(Mcal, 2, 107, 3, MCALDIM) =   5;
        LOC2(Mcal, 2, 108, 3, MCALDIM) =   2;
        LOC2(Mcal, 2, 109, 3, MCALDIM) =   0;
        LOC2(Mcal, 2, 110, 3, MCALDIM) =   0;
        LOC2(Mcal, 2, 111, 3, MCALDIM) =   4;
        LOC2(Mcal, 2, 112, 3, MCALDIM) =   3;
        LOC2(Mcal, 2, 113, 3, MCALDIM) =   4;
        LOC2(Mcal, 2, 114, 3, MCALDIM) =   3;
        LOC2(Mcal, 2, 115, 3, MCALDIM) =   0;
        LOC2(Mcal, 2, 116, 3, MCALDIM) =   0;
        LOC2(Mcal, 2, 117, 3, MCALDIM) =   0;
        LOC2(Mcal, 2, 118, 3, MCALDIM) =   0;
        LOC2(Mcal, 2, 119, 3, MCALDIM) =   7;
    }
    
    // upload to trans device location
    cudaError_t status;
#ifdef TEST
    memcpy(devTrans, trans, sizeof(QUICKDouble)*TRANSDIM*TRANSDIM*TRANSDIM);
    memcpy(devMcal, Mcal, sizeof(QUICKDouble)*3*MCALDIM);
#else    
    status = cudaMemcpyToSymbol(devTrans, trans, sizeof(int)*TRANSDIM*TRANSDIM*TRANSDIM);
    PRINTERROR(status, " cudaMemcpyToSymbol, Trans copy to constants failed")
    
    status = cudaMemcpyToSymbol(devMcal, Mcal, sizeof(int)*3*MCALDIM);
    PRINTERROR(status, " cudaMemcpyToSymbol, Mcal copy to constants failed")
#endif
}

__global__ void get2e_kernel()
{
    const unsigned int tidx = threadIdx.x;
    const unsigned int tidy = threadIdx.y;
    const unsigned int bidx = blockIdx.x;
    const unsigned int bidy = blockIdx.y;
#ifndef TEST    
    gpu_shell(bidx,bidy,tidx,tidy);
#endif
}

#ifndef TEST
__device__
#endif
 void gpu_shell( unsigned int II, unsigned int JJ, unsigned int KK, unsigned int LL)
{

    QUICKDouble DNMax1 = MAX(4.0*LOC2(devSim.cutMatrix, II-1, JJ-1, devSim.nshell, devSim.nshell),
                             4.0*LOC2(devSim.cutMatrix, KK-1, LL-1, devSim.nshell, devSim.nshell));
    QUICKDouble DNMax2 = MAX(MAX(LOC2(devSim.cutMatrix, II-1, LL-1, devSim.nshell, devSim.nshell),
                                 LOC2(devSim.cutMatrix, II-1, KK-1, devSim.nshell, devSim.nshell)),
                             MAX(LOC2(devSim.cutMatrix, JJ-1, KK-1, devSim.nshell, devSim.nshell),
                                 LOC2(devSim.cutMatrix, JJ-1, LL-1, devSim.nshell, devSim.nshell)));
    QUICKDouble DNMax  = MAX(DNMax1, DNMax2);
    if ((LOC2(devSim.YCutoff, KK-1, LL-1, devSim.nshell, devSim.nshell)         > devSim.integralCutoff) &&
        (LOC2(devSim.YCutoff, KK-1, LL-1, devSim.nshell, devSim.nshell) * DNMax > devSim.integralCutoff)) {
        
            int NII1 = devSim.Qstart[II-1];
            int NII2 = devSim.Qfinal[II-1];
            int NJJ1 = devSim.Qstart[JJ-1];
            int NJJ2 = devSim.Qfinal[JJ-1];
            int NKK1 = devSim.Qstart[KK-1];
            int NKK2 = devSim.Qfinal[KK-1];
            int NLL1 = devSim.Qstart[LL-1];
            int NLL2 = devSim.Qfinal[LL-1];
            
            for (int i = NII1; i<= NII2; i++) {
                for (int j = NJJ1; j<=NJJ2; j++) {
                    for (int k = NKK1; k <= NKK2; k++) {
                        for (int l = NLL1; l<= NLL2; l++) {
                            iclass(i,j,k,l,II,JJ,KK,LL, DNMax);
                        }
                    }
                }
                
            }
        }
    
}


#ifndef TEST
__device__
#endif
void iclass(int I, int J, int K, int L, unsigned int II, unsigned int JJ, unsigned int KK, unsigned int LL, QUICKDouble DNMax)
{
    QUICKDouble RA[3],RB[3],RC[3],RD[3],AAtemp[3],P[3],Ptemp[3],Q[3],W[3];
    QUICKDouble Qtemp[3], WQtemp[3], WPtemp[3];
    QUICKDouble FM[14] = {0};
    QUICKDouble store[STOREDIM*STOREDIM]={0};
    int katomA,katomB,katomC,katomD;
    int NII2,NJJ2,NKK2,NLL2;
    int NNAB,NNCD,NNA,NNC,NABCD,NABCDTYPE;
    
    const QUICKDouble PI = 3.1415926535897932384626433832795;
    const QUICKDouble X0 = 2.0 * pow (PI, 2.5);
    
    NNA  = Sumindex[I+1]+1;
    NNAB = Sumindex[I+J+2];
    NNC  = Sumindex[K+1]+1;
    NNCD = Sumindex[K+L+2];
    
    
    /* 
     kAtom A, B, C ,D is the coresponding atom for shell ii, jj, kk, ll
     and be careful with the index difference between Fortran and C++, 
     Fortran starts array index with 1 and C++ starts 0.
     */
    katomA = devSim.katom[II-1];
    katomB = devSim.katom[JJ-1];
    katomC = devSim.katom[KK-1];
    katomD = devSim.katom[LL-1];
    
    /*
     NII1 is the starting angular momenta for shell i and NII2 is the ending
     angular momenta.So it is with other varibles
     */
    
    NII2   = devSim.Qfinal[II-1];
    NJJ2   = devSim.Qfinal[JJ-1];
    NKK2   = devSim.Qfinal[KK-1];
    NLL2   = devSim.Qfinal[LL-1];
    
    NABCDTYPE =(int) (NII2+NJJ2)*10u+NKK2+NLL2; // NABCDTYPE is used for hrr
    NABCD= (int) NII2+NJJ2+NKK2+NLL2;
    
    /*
     RA, RB, RC, and RD are the coordinates for atom katomA, katomB, katomC and katomD, 
     which means they are corrosponding coorinates for shell II, JJ, KK, and LL.
     */
    for (int i = 0; i<3; i++) {
        RA[i]=(QUICKDouble) LOC2(devSim.xyz, i , katomA-1, 3, devSim.natom);
        RB[i]=(QUICKDouble) LOC2(devSim.xyz, i , katomB-1, 3, devSim.natom);
        RC[i]=(QUICKDouble) LOC2(devSim.xyz, i , katomC-1, 3, devSim.natom); 
        RD[i]=(QUICKDouble) LOC2(devSim.xyz, i , katomD-1, 3, devSim.natom);               
    }
    
    
    for (int JJJ = 0; JJJ < devSim.kprim[JJ-1]; JJJ++) {
        for (int III = 0; III < devSim.kprim[II-1]; III++) {
            
            /* In the following comments, we have I, J, K, L denote the primitive gaussian function we use, and
             for example, expo(III, ksumtype(II)) stands for the expo for the IIIth primitive guassian function for II shell, 
             we use I to express the corresponding index.
             AA = expo(I)
             BB = expo(J)
             AB = expo(I)+expo(J)
             1
             ABtemp = -------------------
             2(expo(I) + expo(J))
             */
            QUICKDouble AA = (QUICKDouble) LOC2(devSim.gcexpo, III , devSim.Ksumtype[II-1]-1, 6, devSim.nbasis); 
            QUICKDouble BB = (QUICKDouble) LOC2(devSim.gcexpo, JJJ , devSim.Ksumtype[JJ-1]-1, 6, devSim.nbasis);
            QUICKDouble AB = AA + BB;
            QUICKDouble ABtemp = (QUICKDouble) 0.5 /AB;
            /*
             --->                --->
             ->     expo(I) * xyz (I) + expo(J) * xyz(J)
             P  = ---------------------------------------
             expo(I) + expo(J)
             
             -->            -->
             ----->         expo(I)*xyz(I)+expo(J)*xyz(J)                                 -->            -->
             AAtemp = ----------------------------------- * (expo(I) + expo(J)) = expo(I)*xyz(I)+expo(J)*xyz(J)
             expo(I) + expo(J)
             
             ----->   ->  ->
             Ptemp  = P - A
             */
            for (int i = 0; i<3; i++) {
                P[i] = (LOC2(devSim.xyz, i, devSim.katom[II-1]-1, 3, devSim.natom ) * AA +  
                        LOC2(devSim.xyz, i, devSim.katom[JJ-1]-1, 3, devSim.natom ) * BB) / (AA + BB) ;
                AAtemp[i] = AB * P[i];
                Ptemp[i] = P[i] - RA[i];
            }
            for (int LLL = 0 ; LLL < devSim.kprim[LL-1]; LLL++) {
                for (int KKK = 0; KKK < devSim.kprim[KK-1]; KKK++) {
                
                    QUICKDouble cutoffPrim = DNMax * LOC2(devSim.cutPrim, devSim.kstart[II-1]+III-1, devSim.kstart[JJ-1]+JJJ-1, devSim.jbasis, devSim.jbasis) \
                                                   * LOC2(devSim.cutPrim, devSim.kstart[KK-1]+KKK-1, devSim.kstart[LL-1]+LLL-1, devSim.jbasis, devSim.jbasis);
                    if (cutoffPrim > devSim.primLimit) {
                    /*
                     CC = expo(L)
                     DD = expo(K)
                     CD = expo(L)+expo(K)
                     1
                     CDtemp = ----------------------
                     2(expo(I) + expo(J))
                     ABCD = AB + CD = expo(I)+expo(J)+expo(K)+expo(L)
                     AB * CD      (expo(I)+expo(J))*(expo(K)+expo(L))
                     Rou(Greek Letter) =   ----------- = ------------------------------------
                     AB + CD         expo(I)+expo(J)+expo(K)+expo(L)
                     _______________________________
                     ABCDsqrt = \/expo(I)+expo(J)+expo(K)+expo(L)
                     expo(I)+expo(J)                        expo(K)+expo(L)
                     ABcom = --------------------------------  CDcom = --------------------------------
                     expo(I)+expo(J)+expo(K)+expo(L)           expo(I)+expo(J)+expo(K)+expo(L)
                     
                     ABCDtemp = 1/2(expo(I)+expo(J)+expo(K)+expo(L))                    
                     */
                    QUICKDouble CC = (QUICKDouble) LOC2(devSim.gcexpo, KKK, devSim.Ksumtype[KK-1]-1, 6, devSim.nbasis);
                    QUICKDouble DD = (QUICKDouble) LOC2(devSim.gcexpo, LLL, devSim.Ksumtype[LL-1]-1, 6, devSim.nbasis);  
                    QUICKDouble CD = (QUICKDouble) CC + DD;
                    QUICKDouble ABCD = (QUICKDouble) AB + CD;
                    QUICKDouble ROU = (QUICKDouble) AB*CD/ABCD;
                    QUICKDouble RPQ = (QUICKDouble) 0.0;
                    QUICKDouble ABCDsqrt = (QUICKDouble) sqrt(ABCD);
                    QUICKDouble CDtemp = (QUICKDouble) 0.5/ CD;
                    QUICKDouble ABcom = (QUICKDouble) AB/ABCD;
                    QUICKDouble CDcom = (QUICKDouble) CD/ABCD;
                    QUICKDouble ABCDtemp = (QUICKDouble) 0.5/ABCD;
                    /*
                     Q' is the weighting center of K and L
                     --->           --->
                     ->  ------>       expo(K)*xyz(K)+expo(L)*xyz(L)
                     Q = P'(K,L)  = ------------------------------
                     expo(K) + expo(L)
                     
                     W' is the weight center for I, J, K, L
                     
                     --->             --->             --->            --->
                     ->     expo(I)*xyz(I) + expo(J)*xyz(J) + expo(K)*xyz(K) +expo(L)*xyz(L)
                     W = -------------------------------------------------------------------
                     expo(I) + expo(J) + expo(K) + expo(L)                                            
                     ->  ->  2
                     RPQ =| P - Q | 
                     
                     */
                    for (int i = 0; i<3; i++) {
                        Q[i] = (LOC2(devSim.xyz, i, devSim.katom[KK-1]-1, 3, devSim.natom ) * CC +  \
                                LOC2(devSim.xyz, i, devSim.katom[LL-1]-1, 3, devSim.natom ) * DD) / (CC + DD) ;
                        W[i] = (AAtemp[i] + Q[i] * CD) / ABCD;
                        
                        RPQ  += (Q[i]-P[i])*(Q[i]-P[i]);
                        /*
                         ---->   ->  ->
                         Qtemp = Q - K
                         ----->   ->  ->
                         WQtemp = W - Q
                         ----->   ->  ->
                         WPtemp = W - P
                         */
                        Qtemp[i] = Q[i]-RC[i];
                        WQtemp[i] = W[i] - Q[i];
                        WPtemp[i] = W[i] - P[i];
                    }
                    /*
                     ->  -> 2
                     T = ROU * | P - Q|
                     */
                    QUICKDouble T = (QUICKDouble) RPQ * ROU;
                    FmT(NABCD, T, FM);
                    QUICKDouble YVerticalTemp[VDIM1*VDIM2*VDIM3] = {(QUICKDouble) 0.0};
                    /*
                     Starting point for Vertical Recursion.
                     */
                    for (int i = 0; i<= NABCD; i++) {
                        LOC3(YVerticalTemp, 0, 0, i,VDIM1, VDIM2, VDIM3) = FM[i]/ABCDsqrt;
                    }
                    vertical(NABCDTYPE, YVerticalTemp, Ptemp, WPtemp, Qtemp, WQtemp, ABCDtemp, ABtemp, CDtemp, ABcom, CDcom);
                    
                    /*
                     X2 is the multiplication of four indices normalized coeffecient
                     */
                    QUICKDouble DAB = (QUICKDouble) 0.0;
                    QUICKDouble DCD = (QUICKDouble) 0.0;
                    for (int i=0; i<3; i++) {
                        DAB += (QUICKDouble) pow((LOC2(devSim.xyz, i, devSim.katom[II-1]-1, 3, devSim.natom)-
                                    LOC2(devSim.xyz, i, devSim.katom[JJ-1]-1, 3, devSim.natom)),2);
                        
                        DCD += (QUICKDouble) pow((LOC2(devSim.xyz, i, devSim.katom[KK-1]-1, 3, devSim.natom)-
                                    LOC2(devSim.xyz, i, devSim.katom[LL-1]-1, 3, devSim.natom)),2);
                    }
                    QUICKDouble X2 = (QUICKDouble) exp(-AA*BB/(AA+BB)*DAB)/(AA+BB);
                    X2 = X2 * X0 * LOC2(devSim.gccoeff, III, devSim.Ksumtype[II-1] + I - 1, 6, devSim.nbasis) * \
                    LOC2(devSim.gccoeff, JJJ, devSim.Ksumtype[JJ-1] + J - 1, 6, devSim.nbasis);
                    
                    X2 = X2 * (exp(-CC*DD/(CC+DD)*DCD))/(CC+DD);
                    X2 = X2 * LOC2(devSim.gccoeff, LLL, devSim.Ksumtype[LL-1] + L - 1, 6, devSim.nbasis) *
                    LOC2(devSim.gccoeff, KKK, devSim.Ksumtype[KK-1] + K - 1, 6, devSim.nbasis);
                    for (int i = NNC; i<= NNCD; i++) {
                        for (int j = NNA; j<= NNAB; j++) {
                            LOC2(store, j-1, i-1, STOREDIM, STOREDIM) += X2 * LOC3(YVerticalTemp, j-1,i-1,0, VDIM1, VDIM2, VDIM3);
                        }
                    }
                    }
                }
            }
            
        }
    }
    
    int NBI1 = LOC2(devSim.Qsbasis, II-1, I, devSim.nshell, 4);
    int NBI2 = LOC2(devSim.Qfbasis, II-1, I, devSim.nshell, 4);
    int NBJ1 = LOC2(devSim.Qsbasis, JJ-1, J, devSim.nshell, 4);
    int NBJ2 = LOC2(devSim.Qfbasis, JJ-1, J, devSim.nshell, 4);
    int NBK1 = LOC2(devSim.Qsbasis, KK-1, K, devSim.nshell, 4);
    int NBK2 = LOC2(devSim.Qfbasis, KK-1, K, devSim.nshell, 4);
    int NBL1 = LOC2(devSim.Qsbasis, LL-1, L, devSim.nshell, 4);
    int NBL2 = LOC2(devSim.Qfbasis, LL-1, L, devSim.nshell, 4);
    
    
    // IJKLTYPE is the I, J, K,L type
    int IJKLTYPE = (int) (1000 * I + 100 *J + 10 * K + L);
    
    // maxIJKL is the max of I,J,K,L
    int maxIJKL = (int)MAX(MAX(I,J),MAX(K,L));
    
    if (((maxIJKL == 2)&&(J != 0 || L!=0)) || (maxIJKL >= 3)) {
        IJKLTYPE = 999;
    }
    
    int III1 = devSim.Ksumtype[II-1]+NBI1;
    int III2 = devSim.Ksumtype[II-1]+NBI2;
    int JJJ1 = devSim.Ksumtype[JJ-1]+NBJ1;
    int JJJ2 = devSim.Ksumtype[JJ-1]+NBJ2;
    int KKK1 = devSim.Ksumtype[KK-1]+NBK1;
    int KKK2 = devSim.Ksumtype[KK-1]+NBK2;
    int LLL1 = devSim.Ksumtype[LL-1]+NBL1;
    int LLL2 = devSim.Ksumtype[LL-1]+NBL2;
    
    if ((II < JJ) && (II < KK) && (KK < LL)) {
        for (int III = III1; III <= III2; III++) {
            for (int JJJ = JJJ1; JJJ <= JJJ2; JJJ++) {
                for (int KKK = KKK1; KKK <= KKK2 ; KKK++) {
                    for (int LLL = LLL1; LLL <= LLL2; LLL++) {
                        QUICKDouble Y = (QUICKDouble) hrrwhole(III, JJJ, KKK, LLL, IJKLTYPE, store, RA, RB, RC, RD);
//                                                printf("c%i  %i %e %i %i %i %i %i %i %i %i\n", IJKLTYPE, NABCDTYPE,Y, II,JJ,KK,LL, III,JJJ,KKK,LLL);
                        QUICKDouble DENSEKI = (QUICKDouble) LOC2(devSim.dense, KKK-1, III-1, devSim.nbasis, devSim.nbasis);
                        QUICKDouble DENSEKJ = (QUICKDouble) LOC2(devSim.dense, KKK-1, JJJ-1, devSim.nbasis, devSim.nbasis);
                        QUICKDouble DENSELJ = (QUICKDouble) LOC2(devSim.dense, LLL-1, JJJ-1, devSim.nbasis, devSim.nbasis);
                        QUICKDouble DENSELI = (QUICKDouble) LOC2(devSim.dense, LLL-1, III-1, devSim.nbasis, devSim.nbasis);
                        QUICKDouble DENSELK = (QUICKDouble) LOC2(devSim.dense, LLL-1, KKK-1, devSim.nbasis, devSim.nbasis);
                        QUICKDouble DENSEJI = (QUICKDouble) LOC2(devSim.dense, JJJ-1, III-1, devSim.nbasis, devSim.nbasis);
/*
printf("1, %f %f %f %f %f %f\n",DENSEKI, DENSEKJ, DENSELJ, DENSELI, DENSELK,DENSEJI);                                                        

printf("BEFORE\n");
QUICKDouble tmp1;
TESTCODE(JJJ,III)
TESTCODE(LLL,KKK)
TESTCODE(KKK,III)
TESTCODE(LLL,III)                        
TESTCODE(JJJ,KKK)
TESTCODE(JJJ,LLL)
TESTCODE(KKK,JJJ)
TESTCODE(LLL,JJJ)
*/                        // Find the (ij|kl) integrals where j>i, k>i, l>k, and k and j are equal.
                        QUICKULL val1 = (QUICKULL) (fabs(2.0*DENSELK*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                        if ( DENSELK*Y < (QUICKDouble)0.0)
                        val1 = 0ull - val1;
                        
                        QUICKULL val2 = (QUICKULL) (fabs(2.0*DENSEJI*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                        if ( DENSEJI*Y < (QUICKDouble)0.0)
                        val2 = 0ull - val2;
                        
                        QUICKULL val3 = (QUICKULL) (fabs(0.5*DENSELJ*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                        if ( DENSELJ*Y < (QUICKDouble)0.0)
                        val3 = 0ull - val3;
                        
                        QUICKULL val4 = (QUICKULL) (fabs(0.5*DENSEKJ*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                        if ( DENSEKJ*Y < (QUICKDouble)0.0)
                        val4 = 0ull - val4;
                        
                        QUICKULL val5 = (QUICKULL) (fabs(0.5*DENSELI*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                        if ( DENSELI*Y < (QUICKDouble)0.0)
                        val5 = 0ull - val5;
                        
                        QUICKULL val6 = (QUICKULL) (fabs(0.5*DENSEKI*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                        if ( DENSEKI*Y < (QUICKDouble)0.0)
                        val6 = 0ull - val6;
                        
                        QUICKADD(LOC2(devSim.oULL, JJJ-1, III-1, devSim.nbasis, devSim.nbasis), val1);
                        QUICKADD(LOC2(devSim.oULL, LLL-1, KKK-1, devSim.nbasis, devSim.nbasis), val2);
                        QUICKADD(LOC2(devSim.oULL, KKK-1, III-1, devSim.nbasis, devSim.nbasis), 0ull-val3);
                        QUICKADD(LOC2(devSim.oULL, LLL-1, III-1, devSim.nbasis, devSim.nbasis), 0ull-val4);
                        QUICKADD(LOC2(devSim.oULL, JJJ-1, KKK-1, devSim.nbasis, devSim.nbasis), 0ull-val5);
                        QUICKADD(LOC2(devSim.oULL, JJJ-1, LLL-1, devSim.nbasis, devSim.nbasis), 0ull-val6);
                        QUICKADD(LOC2(devSim.oULL, KKK-1, JJJ-1, devSim.nbasis, devSim.nbasis), 0ull-val5);
                        QUICKADD(LOC2(devSim.oULL, LLL-1, JJJ-1, devSim.nbasis, devSim.nbasis), 0ull-val6);
/*printf("AFTER\n");

TESTCODE(JJJ,III)
TESTCODE(LLL,KKK)
TESTCODE(KKK,III)
TESTCODE(LLL,III)                        
TESTCODE(JJJ,KKK)
TESTCODE(JJJ,LLL)
TESTCODE(KKK,JJJ)
TESTCODE(LLL,JJJ)
  */                      
                    }
                }
            }
        }
    }else {
        for (int III = III1; III <= III2; III++) {
            for (int JJJ = MAX(III,JJJ1); JJJ <= JJJ2; JJJ++) {
                for (int KKK = MAX(III,KKK1); KKK <= KKK2; KKK++) {
                    for (int LLL = MAX(KKK,LLL1); LLL <= LLL2; LLL++) {
                        if (III < KKK) {
                            QUICKDouble Y = (QUICKDouble) hrrwhole(III, JJJ, KKK, LLL, IJKLTYPE, store, RA, RB, RC, RD);

                            //                            printf("a%i  %i %e %i %i %i %i %i %i %i %i\n", IJKLTYPE, NABCDTYPE,Y, II,JJ,KK,LL, III,JJJ,KKK,LLL);
                            if ((III < JJJ)&&(KKK < LLL)) {
                                QUICKDouble DENSEKI = (QUICKDouble) LOC2(devSim.dense, KKK-1, III-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEKJ = (QUICKDouble) LOC2(devSim.dense, KKK-1, JJJ-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSELJ = (QUICKDouble) LOC2(devSim.dense, LLL-1, JJJ-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSELI = (QUICKDouble) LOC2(devSim.dense, LLL-1, III-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSELK = (QUICKDouble) LOC2(devSim.dense, LLL-1, KKK-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEJI = (QUICKDouble) LOC2(devSim.dense, JJJ-1, III-1, devSim.nbasis, devSim.nbasis);
    /*                            
printf("2, %f %f %f %f %f %f\n",DENSEKI, DENSEKJ, DENSELJ, DENSELI, DENSELK,DENSEJI);                                
printf("BEFORE\n");
QUICKDouble tmp1;
TESTCODE(JJJ,III)
TESTCODE(LLL,KKK)
TESTCODE(KKK,III)
TESTCODE(LLL,III)                        
TESTCODE(JJJ,KKK)
TESTCODE(JJJ,LLL)
TESTCODE(KKK,JJJ)
TESTCODE(LLL,JJJ)
*/
                                // Find the (ij|kl) integrals where j>i, k>i, l>k, and k and j are equal.
                                QUICKULL val1 = (QUICKULL) (fabs(2.0*DENSELK*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                if ( DENSELK*Y < (QUICKDouble)0.0)
                                val1 = 0ull - val1;
                                
                                QUICKULL val2 = (QUICKULL) (fabs(2.0*DENSEJI*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                if ( DENSEJI*Y < (QUICKDouble)0.0)
                                val2 = 0ull - val2;
//printf("DENSELJ*Y=%f\n",Y);                                
                                QUICKULL val3 = (QUICKULL) (fabs(0.5*DENSELJ*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                if ( DENSELJ*Y < (QUICKDouble)0.0)
                                val3 = 0ull - val3;
                                
                                QUICKULL val4 = (QUICKULL) (fabs(0.5*DENSEKJ*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                if ( DENSEKJ*Y < (QUICKDouble)0.0)
                                val4 = 0ull - val4;
                                
                                QUICKULL val5 = (QUICKULL) (fabs(0.5*DENSELI*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                if ( DENSELI*Y < (QUICKDouble)0.0)
                                val5 = 0ull - val5;
                                
                                QUICKULL val6 = (QUICKULL) (fabs(0.5*DENSEKI*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                if ( DENSEKI*Y < (QUICKDouble)0.0)
                                val6 = 0ull - val6;
                                
                                QUICKADD(LOC2(devSim.oULL, JJJ-1, III-1, devSim.nbasis, devSim.nbasis), val1);
                                QUICKADD(LOC2(devSim.oULL, LLL-1, KKK-1, devSim.nbasis, devSim.nbasis), val2);
                                QUICKADD(LOC2(devSim.oULL, KKK-1, III-1, devSim.nbasis, devSim.nbasis), 0ull-val3);
                                QUICKADD(LOC2(devSim.oULL, LLL-1, III-1, devSim.nbasis, devSim.nbasis), 0ull-val4);
                                QUICKADD(LOC2(devSim.oULL, JJJ-1, KKK-1, devSim.nbasis, devSim.nbasis), 0ull-val5);
                                QUICKADD(LOC2(devSim.oULL, JJJ-1, LLL-1, devSim.nbasis, devSim.nbasis), 0ull-val6);
                                QUICKADD(LOC2(devSim.oULL, KKK-1, JJJ-1, devSim.nbasis, devSim.nbasis), 0ull-val5);
                                QUICKADD(LOC2(devSim.oULL, LLL-1, JJJ-1, devSim.nbasis, devSim.nbasis), 0ull-val6);
/*printf("AFTER\n");
TESTCODE(JJJ,III)
TESTCODE(LLL,KKK)
TESTCODE(KKK,III)
TESTCODE(LLL,III)                        
TESTCODE(JJJ,KKK)
TESTCODE(JJJ,LLL)
TESTCODE(KKK,JJJ)
TESTCODE(LLL,JJJ)
  */                              
                                
                            }else if ((III == JJJ)&&(KKK == LLL)) {
                                
                                // Find  all the (ii|jj) integrals.
                                QUICKDouble DENSEJI = LOC2(devSim.dense, KKK-1, III-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEJJ = LOC2(devSim.dense, KKK-1, KKK-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEII = LOC2(devSim.dense, III-1, III-1, devSim.nbasis, devSim.nbasis);
/*
printf("3, %f %f %f\n",DENSEJI, DENSEJJ, DENSEII);                                                                
printf("BEFORE\n");
QUICKDouble tmp1;
TESTCODE(III,III)
TESTCODE(KKK,KKK)
TESTCODE(KKK,III)
*/
                                QUICKULL val1 = (QUICKULL) (fabs(DENSEJJ*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                if ( DENSEJJ*Y < (QUICKDouble)0.0)
                                val1 = 0ull - val1;                               
                                
                                QUICKULL val2 = (QUICKULL) (fabs(DENSEII*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                if ( DENSEII*Y < (QUICKDouble)0.0)
                                val2 = 0ull - val2;
                                
                                QUICKULL val3 = (QUICKULL) (fabs(0.5*DENSEJI*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                if ( DENSEJI*Y < (QUICKDouble)0.0)
                                val3 = 0ull - val3;
                                
                                QUICKADD(LOC2(devSim.oULL, III-1, III-1, devSim.nbasis, devSim.nbasis), val1);
                                QUICKADD(LOC2(devSim.oULL, KKK-1, KKK-1, devSim.nbasis, devSim.nbasis), val2);
                                QUICKADD(LOC2(devSim.oULL, KKK-1, III-1, devSim.nbasis, devSim.nbasis), 0ull-val3);
/*printf("after\n");
TESTCODE(III,III)
TESTCODE(KKK,KKK)
TESTCODE(KKK,III)                                
  */                              
                            }else if ((JJJ == KKK)&&(JJJ==LLL)) {

                                // Find all the (ij|jj) integrals.
                                QUICKDouble DENSEJI = (QUICKDouble) LOC2(devSim.dense, JJJ-1, III-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEJJ = (QUICKDouble) LOC2(devSim.dense, JJJ-1, JJJ-1, devSim.nbasis, devSim.nbasis);
/*printf("4, %f %f\n",DENSEJI, DENSEJJ);                                
printf("BEFORE\n");
QUICKDouble tmp1;
TESTCODE(JJJ,III)
TESTCODE(JJJ,JJJ)
  */                            QUICKULL val1 = (QUICKULL) (fabs(0.5*DENSEJJ*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                if ( DENSEJJ*Y < (QUICKDouble)0.0)
                                val1 = 0ull - val1;                               
                                
                                QUICKULL val2 = (QUICKULL) (fabs(DENSEJI*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                if ( DENSEJI*Y < (QUICKDouble)0.0)
                                val2 = 0ull - val2;
                                
                                QUICKADD(LOC2(devSim.oULL, JJJ-1, III-1, devSim.nbasis, devSim.nbasis), val1);
                                QUICKADD(LOC2(devSim.oULL, JJJ-1, JJJ-1, devSim.nbasis, devSim.nbasis), val2);
/*printf("after\n");
TESTCODE(JJJ,III)
TESTCODE(JJJ,JJJ)
  */                          }else if ((KKK == LLL)&&(III<JJJ)&&(JJJ!=KKK)) {
                                
                                //Find all the (ij|kk) integrals where j>i, k>j.
                                QUICKDouble DENSEKI = (QUICKDouble) LOC2(devSim.dense, KKK-1, III-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEKJ = (QUICKDouble) LOC2(devSim.dense, KKK-1, JJJ-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEKK = (QUICKDouble) LOC2(devSim.dense, KKK-1, KKK-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEJI = (QUICKDouble) LOC2(devSim.dense, JJJ-1, III-1, devSim.nbasis, devSim.nbasis);
//printf("5, %f %f %f %f\n",DENSEKI, DENSEKJ, DENSEKK,DENSEJI);                                                                
                                QUICKULL val1 = (QUICKULL) (fabs(DENSEKK*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                if ( DENSEKK*Y < (QUICKDouble)0.0)
                                val1 = 0ull - val1;                               
                                
                                QUICKULL val2 = (QUICKULL) (fabs(2.0*DENSEJI*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                if ( DENSEJI*Y < (QUICKDouble)0.0)
                                val2 = 0ull - val2;
                                
                                QUICKULL val3 = (QUICKULL) (fabs(0.5*DENSEKJ*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                if ( DENSEKJ*Y < (QUICKDouble)0.0)
                                val3 = 0ull - val3;
                                
                                QUICKULL val4 = (QUICKULL) (fabs(0.5*DENSEKI*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                if ( DENSEKI*Y < (QUICKDouble)0.0)
                                val4 = 0ull - val4;
/*printf("BEFORE\n");
QUICKDouble tmp1;
TESTCODE(JJJ,III)
TESTCODE(KKK,KKK)                                
TESTCODE(KKK,III)                                
TESTCODE(KKK,JJJ)
TESTCODE(JJJ,KKK)                                                                
  */                              
                                QUICKADD(LOC2(devSim.oULL, JJJ-1, III-1, devSim.nbasis, devSim.nbasis), val1);
                                QUICKADD(LOC2(devSim.oULL, KKK-1, KKK-1, devSim.nbasis, devSim.nbasis), val2);
                                QUICKADD(LOC2(devSim.oULL, KKK-1, III-1, devSim.nbasis, devSim.nbasis), 0ull-val3);
                                QUICKADD(LOC2(devSim.oULL, KKK-1, JJJ-1, devSim.nbasis, devSim.nbasis), 0ull-val4);
                                QUICKADD(LOC2(devSim.oULL, JJJ-1, KKK-1, devSim.nbasis, devSim.nbasis), 0ull-val4);
/*printf("after\n");
TESTCODE(JJJ,III)
TESTCODE(KKK,KKK)                                
TESTCODE(KKK,III)                                
TESTCODE(KKK,JJJ)
TESTCODE(JJJ,KKK)                                                                                                
  */                              
                            }else if ((III==JJJ)&&(KKK<LLL)) {
                                
                                //Find all the (ii|jk) integrals where j>i, k>j.
                                QUICKDouble DENSEII = (QUICKDouble) LOC2(devSim.dense, III-1, III-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEJI = (QUICKDouble) LOC2(devSim.dense, KKK-1, III-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEKI = (QUICKDouble) LOC2(devSim.dense, LLL-1, III-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEKJ = (QUICKDouble) LOC2(devSim.dense, LLL-1, KKK-1, devSim.nbasis, devSim.nbasis);
//printf("6, %f %f %f %f\n",DENSEII, DENSEJI, DENSEKI,DENSEKJ);                                                                
                                
                                QUICKULL val1 = (QUICKULL) (fabs(DENSEII*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                if ( DENSEII*Y < (QUICKDouble)0.0)
                                val1 = 0ull - val1;                               
                                
                                QUICKULL val2 = (QUICKULL) (fabs(2.0*DENSEKJ*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                if ( DENSEKJ*Y < (QUICKDouble)0.0)
                                val2 = 0ull - val2;
                                
                                QUICKULL val3 = (QUICKULL) (fabs(0.5*DENSEKI*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                if ( DENSEKI*Y < (QUICKDouble)0.0)
                                val3 = 0ull - val3;
                                
                                QUICKULL val4 = (QUICKULL) (fabs(0.5*DENSEJI*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                if ( DENSEJI*Y < (QUICKDouble)0.0)
                                val4 = 0ull - val4;
/*printf("BEFORE\n");
QUICKDouble tmp1;
TESTCODE(LLL,KKK)
TESTCODE(III,III)                                
TESTCODE(KKK,III)                                
TESTCODE(LLL,III)
  */                              
                                QUICKADD(LOC2(devSim.oULL, LLL-1, KKK-1, devSim.nbasis, devSim.nbasis), val1);
                                QUICKADD(LOC2(devSim.oULL, III-1, III-1, devSim.nbasis, devSim.nbasis), val2);
                                QUICKADD(LOC2(devSim.oULL, KKK-1, III-1, devSim.nbasis, devSim.nbasis), 0ull-val3);
                                QUICKADD(LOC2(devSim.oULL, LLL-1, III-1, devSim.nbasis, devSim.nbasis), 0ull-val4);
/*printf("after\n");
TESTCODE(LLL,KKK)
TESTCODE(III,III)                                
TESTCODE(KKK,III)                                
TESTCODE(LLL,III)                                
  */                          } 
                            
                        }else {
                        
                            if (JJJ <= LLL) {
                                QUICKDouble Y = (QUICKDouble) hrrwhole(III, JJJ, KKK, LLL, IJKLTYPE, store, RA, RB, RC, RD);
                            //    printf("b %i  %i %e %i %i %i %i %i %i %i %i\n", IJKLTYPE, NABCDTYPE, Y, II,JJ,KK,LL, III,JJJ,KKK,LLL);
                                if((III==JJJ)&&(III==KKK)&&(III==LLL)){
                                    // do all the (ii|ii) integrals
                                    QUICKDouble DENSEII = (QUICKDouble) LOC2(devSim.dense, III-1, III-1, devSim.nbasis, devSim.nbasis);
//printf("7, %f\n",DENSEII);                              
                                    QUICKULL val1 = (QUICKULL) (fabs(0.5*DENSEII*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                    if ( DENSEII*Y < (QUICKDouble)0.0)
                                    val1 = 0ull - val1;                               
/*printf("BEFORE\n");
QUICKDouble tmp1;
TESTCODE(III,III)                                    
  */                                  QUICKADD(LOC2(devSim.oULL, III-1, III-1, devSim.nbasis, devSim.nbasis), val1);
/*printf("after\n");
TESTCODE(III,III)                                    
  */                              }else if ((III == JJJ) && (III == KKK) && (III < LLL)){
                                    
                                    QUICKDouble DENSEJI = (QUICKDouble) LOC2(devSim.dense, LLL-1, III-1, devSim.nbasis, devSim.nbasis);
                                    QUICKDouble DENSEII = (QUICKDouble) LOC2(devSim.dense, III-1, III-1, devSim.nbasis, devSim.nbasis);
/*
printf("8, %f %f\n",DENSEJI, DENSEII);
QUICKDouble tmp1,tmp2;
printf("BEFORE\n");
TESTCODE(LLL,III)
TESTCODE(III,III)
*/
/*
            if (LOC2(devSim.oULL, LLL-1, III-1, devSim.nbasis, devSim.nbasis) >= 0x8000000000000000ull) {
                tmp1  = -(QUICKDouble)(LOC2(devSim.oULL, LLL-1, III-1, devSim.nbasis, devSim.nbasis) ^ 0xffffffffffffffffull);
            }
            else
            {
                tmp1  = (QUICKDouble) LOC2(devSim.oULL, LLL-1, III-1, devSim.nbasis, devSim.nbasis);
            }       
            
            if (LOC2(devSim.oULL, III-1, III-1, devSim.nbasis, devSim.nbasis) >= 0x8000000000000000ull) {
                tmp2  = -(QUICKDouble)(LOC2(devSim.oULL, III-1, III-1, devSim.nbasis, devSim.nbasis) ^ 0xffffffffffffffffull);
            }
            else
            {
                tmp2  = (QUICKDouble) LOC2(devSim.oULL, III-1, III-1, devSim.nbasis, devSim.nbasis);
            }       
            
            printf("before = %f %f\n",tmp1*ONEOVERENERGYSCALE, \
                             tmp2*ONEOVERENERGYSCALE);
                                                                      
*/                  

                                    QUICKULL val1 = (QUICKULL) (fabs(0.5*DENSEII*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                    if ( DENSEII*Y < (QUICKDouble)0.0)
                                    val1 = 0ull - val1;                               
                                    
                                    QUICKULL val2 = (QUICKULL) (fabs(DENSEJI*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                    if ( DENSEJI*Y < (QUICKDouble)0.0)
                                    val2 = 0ull - val2;                               
                                    
                                    QUICKADD(LOC2(devSim.oULL, LLL-1, III-1, devSim.nbasis, devSim.nbasis), val1);
                                    QUICKADD(LOC2(devSim.oULL, III-1, III-1, devSim.nbasis, devSim.nbasis), val2);

/*
printf("after\n");

TESTCODE(LLL,III)
TESTCODE(III,III)
*/
/*
            if (LOC2(devSim.oULL, LLL-1, III-1, devSim.nbasis, devSim.nbasis) >= 0x8000000000000000ull) {
                tmp1  = -(QUICKDouble)(LOC2(devSim.oULL, LLL-1, III-1, devSim.nbasis, devSim.nbasis) ^ 0xffffffffffffffffull);
            }
            else
            {
                tmp1  = (QUICKDouble) LOC2(devSim.oULL, LLL-1, III-1, devSim.nbasis, devSim.nbasis);
            }       
            
            if (LOC2(devSim.oULL, III-1, III-1, devSim.nbasis, devSim.nbasis) >= 0x8000000000000000ull) {
                tmp2  = -(QUICKDouble)(LOC2(devSim.oULL, III-1, III-1, devSim.nbasis, devSim.nbasis) ^ 0xffffffffffffffffull);
            }
            else
            {
                tmp2  = (QUICKDouble) LOC2(devSim.oULL, III-1, III-1, devSim.nbasis, devSim.nbasis);
            }       
            
            printf("after = %f %f\n",tmp1*ONEOVERENERGYSCALE,\
                             tmp2*ONEOVERENERGYSCALE);
  */                                                                    
                                }else if ((III == KKK) && (JJJ == LLL) && (III < JJJ)){
                                    
                                    QUICKDouble DENSEJI = (QUICKDouble) LOC2(devSim.dense, JJJ-1, III-1, devSim.nbasis, devSim.nbasis);
                                    QUICKDouble DENSEJJ = (QUICKDouble) LOC2(devSim.dense, JJJ-1, JJJ-1, devSim.nbasis, devSim.nbasis);
                                    QUICKDouble DENSEII = (QUICKDouble) LOC2(devSim.dense, III-1, III-1, devSim.nbasis, devSim.nbasis);
/*
printf("9, %f %f %f\n",DENSEJI, DENSEJJ,DENSEII);                                    
printf("BEFORE\n");
QUICKDouble tmp1;

TESTCODE(JJJ,III)                                    
TESTCODE(JJJ,JJJ)                                    
TESTCODE(III,III)                                    
*/
                                    QUICKULL val1 = (QUICKULL) (fabs(1.5*DENSEJI*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                    if ( DENSEJI*Y < (QUICKDouble)0.0)
                                    val1 = 0ull - val1;                               
                                    
                                    QUICKULL val2 = (QUICKULL) (fabs(0.5*DENSEII*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                    if ( DENSEII*Y < (QUICKDouble)0.0)
                                    val2 = 0ull - val2;                               
                                    
                                    QUICKULL val3 = (QUICKULL) (fabs(0.5*DENSEJJ*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                    if ( DENSEJJ*Y < (QUICKDouble)0.0)
                                    val3 = 0ull - val3;                               
                                    
                                    
                                    QUICKADD(LOC2(devSim.oULL, JJJ-1, III-1, devSim.nbasis, devSim.nbasis), val1);
                                    QUICKADD(LOC2(devSim.oULL, JJJ-1, JJJ-1, devSim.nbasis, devSim.nbasis), 0ull-val2);
                                    QUICKADD(LOC2(devSim.oULL, III-1, III-1, devSim.nbasis, devSim.nbasis), 0ull-val3);
/*
printf("after\n");
TESTCODE(JJJ,III)                                    
TESTCODE(JJJ,JJJ)                                    
TESTCODE(III,III)                                    
  */                                  
                                }else if ((III == KKK) && (III <  JJJ) && (JJJ < LLL)){
                                    
                                    QUICKDouble DENSEKI = (QUICKDouble) LOC2(devSim.dense, LLL-1, III-1, devSim.nbasis, devSim.nbasis);
                                    QUICKDouble DENSEKJ = (QUICKDouble) LOC2(devSim.dense, LLL-1, JJJ-1, devSim.nbasis, devSim.nbasis);
                                    QUICKDouble DENSEII = (QUICKDouble) LOC2(devSim.dense, III-1, III-1, devSim.nbasis, devSim.nbasis);
                                    QUICKDouble DENSEJI = (QUICKDouble) LOC2(devSim.dense, JJJ-1, III-1, devSim.nbasis, devSim.nbasis);
//printf("10, %f %f %f %f\n",DENSEKI,DENSEKJ, DENSEII, DENSEJI);                                    
                                    QUICKULL val1 = (QUICKULL) (fabs(1.5*DENSEKI*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                    if ( DENSEKI*Y < (QUICKDouble)0.0)
                                    val1 = 0ull - val1;                               
                                    
                                    QUICKULL val2 = (QUICKULL) (fabs(1.5*DENSEJI*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                    if ( DENSEJI*Y < (QUICKDouble)0.0)
                                    val2 = 0ull - val2;                               
                                    
                                    QUICKULL val3 = (QUICKULL) (fabs(1.0*DENSEKJ*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                    if ( DENSEKJ*Y < (QUICKDouble)0.0)
                                    val3 = 0ull - val3;                               
                                    
                                    QUICKULL val4 = (QUICKULL) (fabs(0.5*DENSEII*Y*ENERGYSCALE) + (QUICKDouble)0.5);
                                    if ( DENSEII*Y < (QUICKDouble)0.0)
                                    val4 = 0ull - val4;
/*printf("BEFORE\n");
QUICKDouble tmp1;
TESTCODE(JJJ,III)                                    
TESTCODE(LLL,III)                                    
TESTCODE(III,III)                                    
TESTCODE(LLL,JJJ)                                    
*/
                                    QUICKADD(LOC2(devSim.oULL, JJJ-1, III-1, devSim.nbasis, devSim.nbasis), val1);
                                    QUICKADD(LOC2(devSim.oULL, LLL-1, III-1, devSim.nbasis, devSim.nbasis), val2);
                                    QUICKADD(LOC2(devSim.oULL, III-1, III-1, devSim.nbasis, devSim.nbasis), 0ull-val3);
                                    QUICKADD(LOC2(devSim.oULL, LLL-1, JJJ-1, devSim.nbasis, devSim.nbasis), 0ull-val4);
/*printf("after\n");
TESTCODE(JJJ,III)                                    
TESTCODE(LLL,III)                                    
TESTCODE(III,III)                                    
TESTCODE(LLL,JJJ)
  */                              }
                            }
                        }
                    }
                }
            }
        }
    }
}


#ifndef TEST
__device__
#endif
  void FmT(int MaxM, QUICKDouble X, QUICKDouble* vals)
{
    const QUICKDouble PI = (QUICKDouble) 3.14159265358E0 ;
    const QUICKDouble PIE4 = (QUICKDouble) PI/4.0 ;
    const QUICKDouble XINV = (QUICKDouble) 1.0 /X;
    const QUICKDouble E = (QUICKDouble) exp(-X);
    QUICKDouble WW1;
    
    if (X > 5.0 ) {
        if (X>15.0 ) {
            if (X>33.0 ) {
                WW1 = sqrt(PIE4 * XINV);
            }else {
                WW1 = (( 1.9623264149430E-01 *XINV-4.9695241464490E-01 )*XINV - \
                       6.0156581186481E-05 )*E + sqrt(PIE4*XINV);
            }
        }else if (X>10.0 ) {
            WW1 = (((-1.8784686463512E-01 *XINV+2.2991849164985E-01 )*XINV - \
                    4.9893752514047E-01 )*XINV-2.1916512131607E-05 )*E + sqrt(PIE4*XINV);
        }else {
            WW1 = (((((( 4.6897511375022E-01  *XINV-6.9955602298985E-01 )*XINV + \
                       5.3689283271887E-01 )*XINV-3.2883030418398E-01 )*XINV + \
                     2.4645596956002E-01 )*XINV-4.9984072848436E-01 )*XINV - \
                   3.1501078774085E-06 )*E + sqrt(PIE4*XINV);
        }
    }else if (X >1.0 ) {
        if (X>3.0 ) {
            QUICKDouble Y = (QUICKDouble) X - 4.0 ;
            QUICKDouble F1 = ((((((((((-2.62453564772299E-11 *Y+3.24031041623823E-10  )*Y- \
                                      3.614965656163E-09 )*Y+3.760256799971E-08 )*Y- \
                                    3.553558319675E-07 )*Y+3.022556449731E-06 )*Y- \
                                  2.290098979647E-05 )*Y+1.526537461148E-04 )*Y- \
                                8.81947375894379E-04 )*Y+4.33207949514611E-03 )*Y- \
                              1.75257821619926E-02 )*Y+5.28406320615584E-02 ;
            WW1 = (X+X)*F1+E;
        }else {
            QUICKDouble Y = (QUICKDouble) X - 2.0 ;
            QUICKDouble F1 = ((((((((((-1.61702782425558E-10 *Y+1.96215250865776E-09  )*Y- \
                                      2.14234468198419E-08  )*Y+2.17216556336318E-07  )*Y- \
                                    1.98850171329371E-06  )*Y+1.62429321438911E-05  )*Y- \
                                  1.16740298039895E-04  )*Y+7.24888732052332E-04  )*Y- \
                                3.79490003707156E-03  )*Y+1.61723488664661E-02  )*Y- \
                              5.29428148329736E-02  )*Y+1.15702180856167E-01 ;
            WW1 = (X+X)*F1+E;
        }
        
    }else if (X > 3.0E-7 ) {
        QUICKDouble F1 =(((((((( -8.36313918003957E-08 *X+1.21222603512827E-06  )*X- \
                               1.15662609053481E-05  )*X+9.25197374512647E-05  )*X- \
                             6.40994113129432E-04  )*X+3.78787044215009E-03  )*X- \
                           1.85185172458485E-02  )*X+7.14285713298222E-02  )*X- \
                         1.99999999997023E-01  )*X+3.33333333333318E-01 ;
        WW1 = (X+X)*F1+E;
    }else {
        WW1 = (1.0 -X)/(QUICKDouble)(2.0 * MaxM+1);
    }
    
    if (X > 3.0E-7 ) {
        vals[0] = WW1;
        for (int m = 1; m<= MaxM; m++) {
            vals[m] = (((2*m-1)*vals[m-1])- E)*0.5*XINV;
        }
    }else {
        vals[MaxM] = WW1;
        for (int m = MaxM-1; m >=0; m--) {
            vals[m] = (2.0 * X * vals[m+1] + E) / (QUICKDouble)(m*2+1);
        }
    }
}

#ifndef TEST
__device__
#endif
 void vertical(int NABCDTYPE, QUICKDouble* YVerticalTemp, QUICKDouble* Ptemp, QUICKDouble* WPtemp, \
                QUICKDouble* Qtemp, QUICKDouble* WQtemp,QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom)
{
    switch (NABCDTYPE) {
        // SSSS oribital
        case 0:
        break;
        // PSSS orbital
        case 10:
        PSSS(0, YVerticalTemp, Ptemp, WPtemp);
        break;
        // SSPS orbital
        case 1:
        SSPS(0, YVerticalTemp, Qtemp, WQtemp);
        break;
        // PSPS orbital
        case 11:
        PSSS(0, YVerticalTemp, Ptemp, WPtemp);
        SSPS(0, YVerticalTemp, Qtemp, WQtemp);
        SSPS(1, YVerticalTemp, Qtemp, WQtemp);
        PSPS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        break;
        case 20:
        PSSS(0, YVerticalTemp, Ptemp, WPtemp);
        PSSS(1, YVerticalTemp, Ptemp, WPtemp);
        DSSS(0, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        break;
        case 2:
        SSPS(0, YVerticalTemp, Qtemp, WQtemp);
        SSPS(1, YVerticalTemp, Qtemp, WQtemp);
        SSDS(0, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        break;
        case 21:
        SSPS(0, YVerticalTemp, Qtemp, WQtemp);
        
        PSSS(0, YVerticalTemp, Ptemp, WPtemp);
        SSPS(1, YVerticalTemp, Qtemp, WQtemp);
        PSPS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        
        PSSS(1, YVerticalTemp, Ptemp, WPtemp);
        DSSS(0, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        
        PSSS(2, YVerticalTemp, Ptemp, WPtemp);
        DSSS(1, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSPS(0, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        break;
        case 12:
        SSPS(0, YVerticalTemp, Qtemp, WQtemp);
        
        PSSS(0, YVerticalTemp, Ptemp, WPtemp);
        SSPS(1, YVerticalTemp, Qtemp, WQtemp);
        PSPS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        
        SSDS(0, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        
        SSPS(2, YVerticalTemp, Qtemp, WQtemp);
        SSDS(1, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        PSDS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        break;
        case 22:
        SSPS(0, YVerticalTemp, Qtemp, WQtemp);
        
        PSSS(0, YVerticalTemp, Ptemp, WPtemp);
        SSPS(1, YVerticalTemp, Qtemp, WQtemp);
        PSPS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        
        SSDS(0, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSPS(2, YVerticalTemp, Qtemp, WQtemp);
        SSDS(1, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        PSDS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        
        PSSS(1, YVerticalTemp, Ptemp, WPtemp);
        DSSS(0, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        PSSS(2, YVerticalTemp, Ptemp, WPtemp);
        DSSS(1, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSPS(0, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        
        SSPS(3, YVerticalTemp, Qtemp, WQtemp);
        SSDS(2, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        PSDS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        
        PSPS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        DSDS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp, ABtemp, CDcom);
        break;
        
        case 30:
        PSSS(0, YVerticalTemp, Ptemp, WPtemp);
        PSSS(1, YVerticalTemp, Ptemp, WPtemp);
        DSSS(0, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        PSSS(2, YVerticalTemp, Ptemp, WPtemp);
        DSSS(1, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        //FSSS(0)
        break;
        
        
        case 3:
        SSPS(0, YVerticalTemp, Qtemp, WQtemp);
        SSPS(1, YVerticalTemp, Qtemp, WQtemp);
        SSDS(0, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSPS(2, YVerticalTemp, Qtemp, WQtemp);
        SSDS(1, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        //SSFS(0)
        break;
        
        case 40:
        PSSS(0, YVerticalTemp, Ptemp, WPtemp);
        PSSS(1, YVerticalTemp, Ptemp, WPtemp);
        DSSS(0, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        PSSS(2, YVerticalTemp, Ptemp, WPtemp);
        DSSS(1, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        //FSSS(0)
        PSSS(3, YVerticalTemp, Ptemp, WPtemp);
        DSSS(2, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        //FSSS(1)
        //GSSS(0)
        break;
        
        case 4:
        SSPS(0, YVerticalTemp, Qtemp, WQtemp);
        SSPS(1, YVerticalTemp, Qtemp, WQtemp);
        SSDS(0, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSPS(2, YVerticalTemp, Qtemp, WQtemp);
        SSDS(1, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        //SSFS(0)
        SSPS(3, YVerticalTemp, Qtemp, WQtemp);
        SSDS(2, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        //SSFS(1)
        //SSGS(0)
        break;
        
        
        case 31:
        PSSS(0, YVerticalTemp, Ptemp, WPtemp);
        PSSS(1, YVerticalTemp, Ptemp, WPtemp);
        DSSS(0, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        PSSS(2, YVerticalTemp, Ptemp, WPtemp);
        DSSS(1, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        //FSSS(0)
        PSSS(3, YVerticalTemp, Ptemp, WPtemp);
        DSSS(2, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        //FSSS(1)
        SSPS(0, YVerticalTemp, Qtemp, WQtemp);
        SSPS(1, YVerticalTemp, Qtemp, WQtemp);
        PSPS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        //DSPS(0, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        //FSPS(0)
        break;
        
        case 13:
        SSPS(0, YVerticalTemp, Qtemp, WQtemp);
        SSPS(1, YVerticalTemp, Qtemp, WQtemp);
        SSDS(0, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSPS(2, YVerticalTemp, Qtemp, WQtemp);
        SSDS(1, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        //SSFS(0)
        SSPS(3, YVerticalTemp, Qtemp, WQtemp);
        SSDS(2, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        //SSFS(1)
        PSSS(0, YVerticalTemp, Ptemp, WPtemp);
        PSPS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSDS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        //PSFS(0)
        break;
        
        case 41:
        PSSS(0, YVerticalTemp, Ptemp, WPtemp);
        PSSS(1, YVerticalTemp, Ptemp, WPtemp);
        DSSS(0, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        PSSS(2, YVerticalTemp, Ptemp, WPtemp);
        DSSS(1, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        //FSSS(0)
        PSSS(3, YVerticalTemp, Ptemp, WPtemp);
        DSSS(2, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        //FSSS(1)
        //GSSS(0)
        PSSS(4, YVerticalTemp, Ptemp, WPtemp);
        DSSS(3, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        //FSSS(2)
        //GSSS(1)
        SSPS(0, YVerticalTemp, Qtemp, WQtemp);
        SSPS(1, YVerticalTemp, Qtemp, WQtemp);
        PSPS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        //DSPS(0, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        //FSPS(0)
        //GSPS(0)
        break;
        
        case 14:
        SSPS(0, YVerticalTemp, Qtemp, WQtemp);
        SSPS(1, YVerticalTemp, Qtemp, WQtemp);
        SSDS(0, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSPS(2, YVerticalTemp, Qtemp, WQtemp);
        SSDS(1, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        //SSFS(0)
        SSPS(3, YVerticalTemp, Qtemp, WQtemp);
        SSDS(2, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        //SSFS(1)
        PSSS(0, YVerticalTemp, Ptemp, WPtemp);
        PSSS(1, YVerticalTemp, Ptemp, WPtemp);
        PSPS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        //SSGS(0)
        SSPS(4, YVerticalTemp, Qtemp, WQtemp);
        SSDS(3, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        //SSFS(2)
        //SSGS(1)
        PSDS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        //PSFS(0)
        //PSGS(0)
        break;
    
        case 32:
        PSSS(0, YVerticalTemp, Ptemp, WPtemp);
        PSSS(1, YVerticalTemp, Ptemp, WPtemp);
        DSSS(0, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        PSSS(2, YVerticalTemp, Ptemp, WPtemp);
        DSSS(1, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        // FSSS(0)
        PSSS(3, YVerticalTemp, Ptemp, WPtemp);
        DSSS(2, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        //FSSS(1)
        SSPS(0, YVerticalTemp, Qtemp, WQtemp);
        SSPS(1, YVerticalTemp, Qtemp, WQtemp);
        PSPS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        DSPS(0, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        //FSPS(0)
        PSSS(4, YVerticalTemp, Ptemp, WPtemp);
        DSSS(3, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        //FSSS(2)
        //FSPS(1)
        DSPS(1, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        SSDS(0, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSPS(2, YVerticalTemp, Qtemp, WQtemp);
        SSDS(1, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        PSDS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        SSPS(3, YVerticalTemp, Qtemp, WQtemp);
        SSDS(2, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        PSDS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSPS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        DSDS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp, ABtemp, CDcom);
        //FSDS(0)
        break;
        
        case 23:
        SSPS(0, YVerticalTemp, Qtemp, WQtemp);
        SSPS(1, YVerticalTemp, Qtemp, WQtemp);
        SSDS(0, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSPS(2, YVerticalTemp, Qtemp, WQtemp);
        SSDS(1, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        //SSFS(0)
        SSPS(3, YVerticalTemp, Qtemp, WQtemp);
        SSDS(2, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        //SSFS(1)
        PSSS(0, YVerticalTemp, Ptemp, WPtemp);
        PSSS(1, YVerticalTemp, Ptemp, WPtemp);
        PSPS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSDS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        //PSFS(0)
        SSPS(4, YVerticalTemp, Qtemp, WQtemp);
        SSDS(3, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        //SSFS(2)
        //PSFS(1)
        PSDS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        DSSS(0, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        PSSS(2, YVerticalTemp, Ptemp, WPtemp);
        DSSS(1, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSPS(0, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        PSPS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        DSDS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp, ABtemp, CDcom);
        //DSFS(0)
        break;
        
        case 42:
        PSSS(0, YVerticalTemp, Ptemp, WPtemp);
        PSSS(1, YVerticalTemp, Ptemp, WPtemp);
        DSSS(0, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        PSSS(2, YVerticalTemp, Ptemp, WPtemp);
        DSSS(1, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        //FSSS(0)
        PSSS(3, YVerticalTemp, Ptemp, WPtemp);
        DSSS(2, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        //FSSS(1)
        //GSSS(0)
        PSSS(4, YVerticalTemp, Ptemp, WPtemp);
        DSSS(3, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        //FSSS(2)
        //GSSS(1)
        DSPS(0, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        //FSPS(0)
        //GSPS(0)
        PSSS(5, YVerticalTemp, Ptemp, WPtemp);
        DSSS(4, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        //FSSS(3)
        //GSSS(2)
        DSPS(1, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        //FSPS(1)
        //GSPS(1)
        SSPS(0, YVerticalTemp, Qtemp, WQtemp);
        SSPS(1, YVerticalTemp, Qtemp, WQtemp);
        PSPS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        SSDS(0, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSPS(2, YVerticalTemp, Qtemp, WQtemp);
        SSDS(1, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        PSDS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        SSPS(3, YVerticalTemp, Qtemp, WQtemp);
        SSDS(2, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        PSDS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSPS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        DSDS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp, ABtemp, CDcom);
        //FSDS(0)
        //GSDS(0)
        break;
        
        case 24:
        SSPS(0, YVerticalTemp, Qtemp, WQtemp);
        SSPS(1, YVerticalTemp, Qtemp, WQtemp);
        SSDS(0, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSPS(2, YVerticalTemp, Qtemp, WQtemp);
        SSDS(1, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        //SSFS(0)
        SSPS(3, YVerticalTemp, Qtemp, WQtemp);
        SSDS(2, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        //SSFS(1)
        //SSGS(0)
        SSPS(4, YVerticalTemp, Qtemp, WQtemp);
        SSDS(3, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        //SSFS(2)
        //SSGS(1)
        PSDS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        //PSFS(0)
        //PSGS(0)
        SSPS(5, YVerticalTemp, Qtemp, WQtemp);
        SSDS(4, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        
        //SSFS(3)
        
        //SSGS(2)
        PSDS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        //PSFS(1)
        
        //PSGS(1)
        
        PSSS(0, YVerticalTemp, Ptemp, WPtemp);
        
        PSSS(1, YVerticalTemp, Ptemp, WPtemp);
        PSPS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        
        DSSS(0, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        PSSS(2, YVerticalTemp, Ptemp, WPtemp);
        DSSS(1, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSPS(0, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        
        PSSS(3, YVerticalTemp, Ptemp, WPtemp);
        DSSS(2, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSPS(1, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        
        PSPS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        DSDS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp, ABtemp, CDcom);
        
        //DSFS(0)
        //DSGS(0)
        break;
        
        
        case 33:
        PSSS(0, YVerticalTemp, Ptemp, WPtemp);
        PSSS(1, YVerticalTemp, Ptemp, WPtemp);
        PSSS(2, YVerticalTemp, Ptemp, WPtemp);
        PSSS(3, YVerticalTemp, Ptemp, WPtemp);
        PSSS(4, YVerticalTemp, Ptemp, WPtemp);
        PSSS(5, YVerticalTemp, Ptemp, WPtemp);
        
        SSPS(0, YVerticalTemp, Qtemp, WQtemp);
        SSPS(1, YVerticalTemp, Qtemp, WQtemp);
        SSPS(2, YVerticalTemp, Qtemp, WQtemp);
        SSPS(3, YVerticalTemp, Qtemp, WQtemp);
        SSPS(4, YVerticalTemp, Qtemp, WQtemp);
        
        PSPS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSPS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSPS(2, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSPS(3, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        
        DSSS(0, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSSS(1, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSSS(2, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSSS(3, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSSS(4, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        
        SSDS(0, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSDS(1, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSDS(2, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSDS(3, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        
        DSPS(0, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        DSPS(1, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        DSPS(2, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        DSPS(3, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        
        PSDS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSDS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSDS(2, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        
        DSDS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp, ABtemp, CDcom);
        DSDS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp, ABtemp, CDcom);
        
        //FSSS(0)
        //FSSS(1)
        //FSSS(2)
        //FSSS(3)
        
        //SSFS(0)
        //SSFS(1)
        //SSFS(2)


        //FSPS(0)
        //FSPS(1)
        //FSPS(2)
        
        //PSFS(0)
        //PSFS(1)
        
        //FSDS(0)
        //DSFS(0)
        //FSDS(1)
        //FSFS(0)
        break;
        
        case 43:
        PSSS(0, YVerticalTemp, Ptemp, WPtemp);
        PSSS(1, YVerticalTemp, Ptemp, WPtemp);
        PSSS(2, YVerticalTemp, Ptemp, WPtemp);
        PSSS(3, YVerticalTemp, Ptemp, WPtemp);
        PSSS(4, YVerticalTemp, Ptemp, WPtemp);
        PSSS(5, YVerticalTemp, Ptemp, WPtemp);
        PSSS(6, YVerticalTemp, Ptemp, WPtemp);
        
        SSPS(0, YVerticalTemp, Qtemp, WQtemp);
        SSPS(1, YVerticalTemp, Qtemp, WQtemp);
        SSPS(2, YVerticalTemp, Qtemp, WQtemp);
        SSPS(3, YVerticalTemp, Qtemp, WQtemp);
        SSPS(4, YVerticalTemp, Qtemp, WQtemp);
        SSPS(5, YVerticalTemp, Qtemp, WQtemp);
        
        PSPS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSPS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSPS(2, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSPS(3, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSPS(4, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        
        DSSS(0, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSSS(1, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSSS(2, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSSS(3, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSSS(4, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSSS(5, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        
        SSDS(0, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSDS(1, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSDS(2, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSDS(3, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSDS(4, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        
        DSPS(0, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        DSPS(1, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        DSPS(2, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        DSPS(3, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        DSPS(4, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        
        PSDS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSDS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSDS(2, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSDS(3, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        
        DSDS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp, ABtemp, CDcom);
        DSDS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp, ABtemp, CDcom);
        DSDS(2, YVerticalTemp, Ptemp, WPtemp, ABCDtemp, ABtemp, CDcom);
        
        //FSSS(0)
        //FSSS(1)
        //FSSS(2)
        //FSSS(3)
        //FSSS(4)
        
        //SSFS(0)
        //SSFS(1)
        //SSFS(2)
        //SSFS(3)
        
        //FSPS(0)
        //FSPS(1)
        //FSPS(2)
        //FSPS(3)
        
        //PSFS(0)
        //PSFS(1)
        //PSFS(2)
    
        //FSDS(0)
        //DSFS(0)
        //FSDS(1)
        //DSFS(1)
        //FSDS(2)
        //FSFS(0)
        
        //GSSS(0)
        //GSSS(1)
        //GSSS(2)
        //GSSS(3)
        
        //GSPS(0)
        //GSPS(1)
        //GSPS(2)
        
        //GSDS(0)
        //GSDS(1)
        //GSFS(0)
        break;
        
        case 34:
        SSPS(0, YVerticalTemp, Qtemp, WQtemp);
        SSPS(1, YVerticalTemp, Qtemp, WQtemp);
        SSPS(2, YVerticalTemp, Qtemp, WQtemp);
        SSPS(3, YVerticalTemp, Qtemp, WQtemp);
        SSPS(4, YVerticalTemp, Qtemp, WQtemp);
        SSPS(5, YVerticalTemp, Qtemp, WQtemp);
        SSPS(6, YVerticalTemp, Qtemp, WQtemp);
        
        PSSS(0, YVerticalTemp, Ptemp, WPtemp);
        PSSS(1, YVerticalTemp, Ptemp, WPtemp);
        PSSS(2, YVerticalTemp, Ptemp, WPtemp);
        PSSS(3, YVerticalTemp, Ptemp, WPtemp);
        PSSS(4, YVerticalTemp, Ptemp, WPtemp);
        PSSS(5, YVerticalTemp, Ptemp, WPtemp);
        
        PSPS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSPS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSPS(2, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSPS(3, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSPS(4, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        
        SSDS(0, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSDS(1, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSDS(2, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSDS(3, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSDS(4, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSDS(5, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        
        DSSS(0, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSSS(1, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSSS(2, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSSS(3, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSSS(4, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        
        PSDS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSDS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSDS(2, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSDS(3, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSDS(4, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        
        DSPS(0, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        DSPS(1, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        DSPS(2, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        DSPS(3, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        
        
        DSDS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp, ABtemp, CDcom);
        DSDS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp, ABtemp, CDcom);
        DSDS(2, YVerticalTemp, Ptemp, WPtemp, ABCDtemp, ABtemp, CDcom);
        
        
        //SSFS(0)
        //SSFS(1)
        //SSFS(2)
        //SSFS(3)
        //SSFS(4)
        
        //FSSS(0)
        //FSSS(1)
        //FSSS(2)
        //FSSS(3)
        
        //PSFS(0)
        //PSFS(1)
        //PSFS(2)
        //PSFS(3)
    
        //FSPS(0)
        //FSPS(1)
        //FSPS(2)
        
        //FSDS(0)
        //DSFS(0)
        //FSDS(1)
        //DSFS(1)
        //DSFS(2)
        //FSFS(0)
        
        //SSGS(0)
        //SSGS(1)
        //SSGS(2)
        //SSGS(3)
        
        //PSGS(0)
        //PSGS(1)
        //PSGS(2)
        
        //DSGS(0)
        //DSGS(1)
        
        //FSGS(0)
        
        break;
        case 44:
        
        PSSS(0, YVerticalTemp, Ptemp, WPtemp);
        PSSS(1, YVerticalTemp, Ptemp, WPtemp);
        PSSS(2, YVerticalTemp, Ptemp, WPtemp);
        PSSS(3, YVerticalTemp, Ptemp, WPtemp);
        PSSS(4, YVerticalTemp, Ptemp, WPtemp);
        PSSS(5, YVerticalTemp, Ptemp, WPtemp);
        PSSS(6, YVerticalTemp, Ptemp, WPtemp);
        PSSS(7, YVerticalTemp, Ptemp, WPtemp);
        
        SSPS(0, YVerticalTemp, Qtemp, WQtemp);
        SSPS(1, YVerticalTemp, Qtemp, WQtemp);
        SSPS(2, YVerticalTemp, Qtemp, WQtemp);
        SSPS(3, YVerticalTemp, Qtemp, WQtemp);
        SSPS(4, YVerticalTemp, Qtemp, WQtemp);
        SSPS(5, YVerticalTemp, Qtemp, WQtemp);
        SSPS(6, YVerticalTemp, Qtemp, WQtemp);
        
        PSPS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSPS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSPS(2, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSPS(3, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSPS(4, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSPS(5, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        
        DSSS(0, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSSS(1, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSSS(2, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSSS(3, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSSS(4, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSSS(5, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        DSSS(6, YVerticalTemp, Ptemp, WPtemp, ABtemp, CDcom);
        
        SSDS(0, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSDS(1, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSDS(2, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSDS(3, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSDS(4, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        SSDS(5, YVerticalTemp, Qtemp, WQtemp, CDtemp, ABcom);
        
        DSPS(0, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        DSPS(1, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        DSPS(2, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        DSPS(3, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        DSPS(4, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        DSPS(5, YVerticalTemp, Qtemp, WQtemp, ABCDtemp);
        
        PSDS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSDS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSDS(2, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSDS(3, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        PSDS(4, YVerticalTemp, Ptemp, WPtemp, ABCDtemp);
        
        DSDS(0, YVerticalTemp, Ptemp, WPtemp, ABCDtemp, ABtemp, CDcom);
        DSDS(1, YVerticalTemp, Ptemp, WPtemp, ABCDtemp, ABtemp, CDcom);
        DSDS(2, YVerticalTemp, Ptemp, WPtemp, ABCDtemp, ABtemp, CDcom);
        DSDS(3, YVerticalTemp, Ptemp, WPtemp, ABCDtemp, ABtemp, CDcom);
        
        //FSSS(0)
        //FSSS(1)
        //FSSS(2)
        //FSSS(3)
        //FSSS(4)
        //FSSS(5)
        
        //SSFS(0)
        //SSFS(1)
        //SSFS(2)
        //SSFS(3)
        //SSFS(4)
        
        //FSPS(0)
        //FSPS(1)
        //FSPS(2)
        //FSPS(3)
        //FSPS(4)
        
        //PSFS(0)
        //PSFS(1)
        //PSFS(2)
        //PSFS(3)
    
        //FSDS(0)
        //FSDS(1)
        //FSDS(2)
        //FSDS(3)
        
        //DSFS(0)
        //DSFS(1)
        //DSFS(2)
        
        //FSFS(0)
        //FSFS(1)
        
        //GSSS(0)
        //GSSS(1)
        //GSSS(2)
        //GSSS(3)
        //GSSS(4)
        
        //SSGS(0)
        //SSGS(1)
        //SSGS(2)
        //SSGS(3)
        
        //GSPS(0)
        //GSPS(1)
        //GSPS(2)
        //GSPS(3)
        
        //PSGS(0)
        //PSGS(1)
        //PSGS(2)
        
        //GSDS(0)
        //GSDS(1)
        //GSDS(2)
        
        //DSGS(0)
        //DSGS(1)
        
        //GSFS(0)
        //GSFS(1)
        
        //FSGS(0)
        //GSGS(0)
        break;
        
        default:
        break;
    }
    
}

#ifndef TEST
__device__
#endif
 void PSSS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble* Ptemp, QUICKDouble* WPtemp)
{
    for (int i = 0; i<3; i++) {
        LOC3(YVerticalTemp, i+1, 0,  mtemp, VDIM1, VDIM2, VDIM3) = Ptemp[i]  * LOC3( YVerticalTemp, 0, 0, mtemp, VDIM1, VDIM2, VDIM3) \
                                                                 + WPtemp[i] * LOC3( YVerticalTemp, 0, 0, mtemp+1, VDIM1, VDIM2, VDIM3);
    }
}

#ifndef TEST
__device__
#endif
  void SSPS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble* Qtemp, QUICKDouble* WQtemp)
{
    for (int i = 0; i<3; i++) {
        LOC3(YVerticalTemp, 0, i+1,  mtemp, VDIM1, VDIM2, VDIM3) = Qtemp[i]  * LOC3( YVerticalTemp, 0, 0, mtemp, VDIM1, VDIM2, VDIM3) \
                                                                 + WQtemp[i] * LOC3( YVerticalTemp, 0, 0, mtemp+1, VDIM1, VDIM2, VDIM3);
    }
}

#ifndef TEST
__device__
#endif
  void PSPS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble* Ptemp, QUICKDouble* WPtemp, QUICKDouble ABCDtemp)
{
    for (int i = 0; i<3; i++) {
        for (int j = 0; j<3; j++) {
            LOC3(YVerticalTemp, i+1 ,j+1, mtemp, VDIM1,VDIM2,VDIM3) = Ptemp[i] * LOC3( YVerticalTemp, 0, j+1, mtemp, VDIM1, VDIM2, VDIM3) \
                                                                   +  WPtemp[i]* LOC3( YVerticalTemp, 0, j+1, mtemp+1, VDIM1, VDIM2, VDIM3);
            if (i == j) {
                LOC3(YVerticalTemp, i+1, j+1, mtemp, VDIM1, VDIM2, VDIM3) += ABCDtemp * LOC3( YVerticalTemp, 0, 0, mtemp+1, VDIM1, VDIM2, VDIM3);
            }
        }
    }
}

#ifndef TEST
__device__
#endif
 void DSSS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble* Ptemp, QUICKDouble* WPtemp, QUICKDouble ABtemp, QUICKDouble CDcom)
{
    int B[3];
    for (int i = 4; i<10; i++) {
        for (int j = 0; j<3; j++) {
            B[j] = LOC2(devMcal, j, i, 3, MCALDIM);
        }
        for (int j = 0; j<3; j++) {
            if (B[j] != 0) {
                B[j] = LOC2(devMcal, j, i, 3, MCALDIM) - 1;
                int ii = (int) LOC3(devTrans, B[0], B[1], B[2], TRANSDIM, TRANSDIM, TRANSDIM);
                LOC3(YVerticalTemp, i, 0, mtemp, VDIM1, VDIM2, VDIM3) =  Ptemp[j] * LOC3( YVerticalTemp, ii-1, 0, mtemp, VDIM1, VDIM2, VDIM3) \
                                                                      + WPtemp[j] * LOC3( YVerticalTemp, ii-1, 0, mtemp+1,VDIM1, VDIM2, VDIM3);
                if (LOC2(devMcal, j, i, 3, MCALDIM) > 1) {
                    LOC3(YVerticalTemp, i, 0, mtemp, VDIM1, VDIM2, VDIM3) += ABtemp * (LOC3(YVerticalTemp, 0, 0, mtemp, VDIM1, VDIM2, VDIM3) \
                                                                             - CDcom * LOC3(YVerticalTemp, 0, 0, mtemp+1, VDIM1, VDIM2, VDIM3));
                }
                break;
            }
        }
    }
}

#ifndef TEST
__device__
#endif
 void SSDS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble* Qtemp, QUICKDouble* WQtemp, QUICKDouble CDtemp, QUICKDouble ABcom)
{
    int B[3];
    for (int i = 4; i<10; i++) {
        for (int j = 0; j<3; j++) {
            B[j] = LOC2(devMcal, j, i, 3, MCALDIM);
        }
        for (int j = 0; j<3; j++) {
            if (B[j] != 0) {
                B[j] = LOC2(devMcal, j, i, 3, MCALDIM) - 1;
                int ii = (int) LOC3(devTrans, B[0], B[1], B[2], TRANSDIM, TRANSDIM, TRANSDIM);
                LOC3(YVerticalTemp, 0, i, mtemp, VDIM1, VDIM2, VDIM3) =  Qtemp[j] * LOC3( YVerticalTemp, 0, ii-1, mtemp, VDIM1, VDIM2, VDIM3) \
                                                                      + WQtemp[j] * LOC3( YVerticalTemp, 0, ii-1, mtemp+1,VDIM1, VDIM2, VDIM3);
                if (LOC2(devMcal, j, i, 3, MCALDIM) > 1) {
                    LOC3(YVerticalTemp, 0, i, mtemp, VDIM1, VDIM2, VDIM3) +=  CDtemp * (LOC3(YVerticalTemp, 0, 0, mtemp, VDIM1, VDIM2, VDIM3) \
                                                                              - ABcom * LOC3(YVerticalTemp, 0, 0, mtemp+1, VDIM1, VDIM2, VDIM3));
                }
                break;
            }
        }
    }
}

#ifndef TEST
__device__
#endif
 void DSPS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble* Qtemp, QUICKDouble* WQtemp, QUICKDouble ABCDtemp)
{
    int B[3];
    for (int i = 4; i<10; i++) {
        for (int jtemp = 2; jtemp <= 4; jtemp++) {
            for (int j = 0; j<3; j++) {
                B[j] = LOC2(devMcal, j, i, 3, MCALDIM);
            }
            for (int j = 0; j<3; j++) {
                if (LOC2(devMcal, j, jtemp-1, 3, MCALDIM) != 0 ) {
                    LOC3(YVerticalTemp, i, jtemp-1, mtemp, VDIM1, VDIM2, VDIM3) =  Qtemp[j] * LOC3( YVerticalTemp, i, 0, mtemp, VDIM1, VDIM2, VDIM3) \
                                                                                + WQtemp[j] * LOC3( YVerticalTemp, i, 0, mtemp+1,VDIM1, VDIM2, VDIM3);
                    if (B[j] != 0) {
                        B[j] = LOC2(devMcal, j, i, 3, MCALDIM) -1;
                        int ii = (int) LOC3(devTrans, B[0], B[1], B[2], TRANSDIM, TRANSDIM, TRANSDIM);
                        LOC3(YVerticalTemp, i, jtemp-1, mtemp, VDIM1, VDIM2, VDIM3) += \
                        ABCDtemp * LOC3(YVerticalTemp, ii-1, 0, mtemp+1, VDIM1, VDIM2, VDIM3) * LOC2(devMcal, j, i, 3, MCALDIM);
                    }
                    break;
                }
            }
        }
    }
}
 
#ifndef TEST
__device__
#endif
 void PSDS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble* Ptemp, QUICKDouble* WPtemp, QUICKDouble ABCDtemp)
 {
    int B[3];
    for (int i = 4; i<10; i++) {
        for (int jtemp = 2; jtemp <= 4; jtemp++) {
            for (int j = 0; j<3; j++) {
                B[j] = LOC2(devMcal, j, i, 3, MCALDIM);
            }
            for (int j = 0; j<3; j++) {
                if (LOC2(devMcal, j, jtemp-1, 3, MCALDIM) != 0 ) {
                    LOC3(YVerticalTemp, jtemp-1, i, mtemp, VDIM1, VDIM2, VDIM3) =  Ptemp[j] * LOC3( YVerticalTemp, 0, i, mtemp, VDIM1, VDIM2, VDIM3) \
                                                                                + WPtemp[j] * LOC3( YVerticalTemp, 0, i, mtemp+1,VDIM1, VDIM2, VDIM3);
                    if (B[j] != 0) {
                        B[j] = LOC2(devMcal, j, i, 3, MCALDIM) -1;
                        int ii = (int) LOC3(devTrans, B[0], B[1], B[2], TRANSDIM, TRANSDIM, TRANSDIM);
                        LOC3(YVerticalTemp, jtemp-1, i, mtemp, VDIM1, VDIM2, VDIM3) += \
                        ABCDtemp * LOC3(YVerticalTemp, 0, ii-1, mtemp+1, VDIM1, VDIM2, VDIM3) * LOC2(devMcal, j, i, 3, MCALDIM);
                    }
                    break;
                }
            }
        }
    }
 }
 
#ifndef TEST
__device__
#endif
void DSDS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble* Ptemp, QUICKDouble* WPtemp, QUICKDouble ABCDtemp, QUICKDouble ABtemp, QUICKDouble CDcom)
{
    int A[3],B[3];
    for (int i = 4; i<10; i++) {
        for (int jtemp = 4; jtemp <10; jtemp++) {
            for (int j = 0; j<3; j++) {
                B[j] = LOC2(devMcal, j, i, 3, MCALDIM);
                A[j] = LOC2(devMcal, j, jtemp, 3, MCALDIM);
            } 
            for (int j = 0; j<3; j++) {
                if (LOC2(devMcal, j, i, 3, MCALDIM) != 0) {
                    B[j] = LOC2(devMcal, j, i, 3, MCALDIM) -1;
                    int ii = (int) LOC3(devTrans, B[0], B[1], B[2], TRANSDIM, TRANSDIM, TRANSDIM);
                    LOC3(YVerticalTemp, i, jtemp, mtemp, VDIM1, VDIM2, VDIM3) = Ptemp[j] * LOC3( YVerticalTemp, ii-1, jtemp, mtemp, VDIM1, VDIM2, VDIM3) + \
                    +WPtemp[j] * LOC3( YVerticalTemp, ii-1, jtemp, mtemp+1, VDIM1, VDIM2, VDIM3);
                    if (LOC2(devMcal, j, i, 3, MCALDIM) >= 2) {
                        LOC3(YVerticalTemp, i, jtemp, mtemp, VDIM1, VDIM2, VDIM3) += \
                          ABtemp * (LOC3(YVerticalTemp, 0, jtemp, mtemp, VDIM1, VDIM2, VDIM3)
                         - CDcom *  LOC3(YVerticalTemp, 0, jtemp, mtemp+1, VDIM1, VDIM2, VDIM3));
                    }
                    if (A[j] != 0) {
                        A[j] = A[j]-1;
                        int iii = (int) LOC3(devTrans,  A[0], A[1], A[2], TRANSDIM, TRANSDIM, TRANSDIM);
                        LOC3(YVerticalTemp, i, jtemp, mtemp, VDIM1, VDIM2, VDIM3) += \
                        ABCDtemp * LOC3(YVerticalTemp, ii-1, iii-1, mtemp+1, VDIM1, VDIM2, VDIM3) * LOC2(devMcal, j, jtemp, 3, MCALDIM);
                    }
                    break;
                }
            }
        }
    }
}
 
 


#ifndef TEST
__device__
#endif
QUICKDouble hrrwhole(int III, int JJJ, int KKK, int LLL, int IJKLTYPE, QUICKDouble* store, \
                     QUICKDouble* RA, QUICKDouble* RB, QUICKDouble* RC, QUICKDouble* RD)
{
    QUICKDouble Y;
    int NA[3],NB[3],NC[3],ND[3];
    int MA,MB,MC,MD;
    for (int i = 0; i<3; i++) {
        NA[i] = (int) LOC2(devSim.KLMN,i,III-1,3,devSim.nbasis);
        NB[i] = (int) LOC2(devSim.KLMN,i,JJJ-1,3,devSim.nbasis);
        NC[i] = (int) LOC2(devSim.KLMN,i,KKK-1,3,devSim.nbasis);
        ND[i] = (int) LOC2(devSim.KLMN,i,LLL-1,3,devSim.nbasis);

    }
    
    MA =  (int) LOC3(devTrans, NA[0], NA[1], NA[2], TRANSDIM, TRANSDIM, TRANSDIM);
    MB =  (int) LOC3(devTrans, NB[0], NB[1], NB[2], TRANSDIM, TRANSDIM, TRANSDIM);
    MC =  (int) LOC3(devTrans, NC[0], NC[1], NC[2], TRANSDIM, TRANSDIM, TRANSDIM);         
    MD =  (int) LOC3(devTrans, ND[0], ND[1], ND[2], TRANSDIM, TRANSDIM, TRANSDIM);
    

    switch (IJKLTYPE) {
        case 0:
        case 10:
        case 1000:
        case 1010:
        {
            Y = (QUICKDouble) LOC2(store, MA-1, MC-1, STOREDIM, STOREDIM);
            break;
        }
        case 2000:
        case 20:
        case 2010:
        case 1020:
        case 2020:
        {
            Y = (QUICKDouble) LOC2(store, MA-1, MC-1, STOREDIM, STOREDIM) * devSim.cons[III-1] * devSim.cons[JJJ-1] * devSim.cons[KKK-1] * devSim.cons[LLL-1];
            break;
        }
        case 100:
        {
            for (int i = 0; i< 3; i++) {
                if (NB[i] != 0) {
                    Y = (QUICKDouble) LOC2(store, MB-1, 0, STOREDIM, STOREDIM) + (RA[i]-RB[i])*LOC2(store, 0, 0, STOREDIM, STOREDIM);
                    break;
                }
            }
            break;
        }
        case 110:
        {
            for (int i = 0; i< 3; i++) {
                if ( NB[i] != 0) {
                    Y = (QUICKDouble) LOC2(store, MB-1, MC-1 , STOREDIM, STOREDIM) + (RA[i]-RB[i])*LOC2(store, 0, MC-1, STOREDIM, STOREDIM);
                    break;
                }
            }
            break;
        }
        case 101:
        {
            QUICKDouble Y1,Y2;
            for (int i = 0; i<3; i++) {
                if (ND[i] != 0) {
                    QUICKDouble c = (QUICKDouble) (RC[i] - RD[i]);
                    Y1 = (QUICKDouble) LOC2(store, MB-1, MD-1 , STOREDIM, STOREDIM) + c * LOC2(store, MB-1,  0, STOREDIM, STOREDIM);
                    Y2 = (QUICKDouble) LOC2(store,    0, MD-1 , STOREDIM, STOREDIM) + c * LOC2(store,    0,  0, STOREDIM, STOREDIM);
                    break;
                }
            }
            for (int i = 0; i<3; i++) {
                if (NB[i] != 0) {
                    Y = Y1 + (RA[i]-RB[i])*Y2;
                    break;
                }
            }
            break;
        }
        case 111:
        {
            QUICKDouble Y1,Y2;
            int MCD = (int) LOC3(devTrans, NC[0]+ND[0], NC[1]+ND[1], NC[2]+ND[2], TRANSDIM, TRANSDIM, TRANSDIM);
            for (int i = 0; i<3; i++) {
                if (ND[i] != 0) {
                    QUICKDouble c = (RC[i] - RD[i]);
                    Y1 = (QUICKDouble) LOC2(store, MB-1, MCD-1 , STOREDIM, STOREDIM) + c * LOC2(store, MB-1, MC-1 , STOREDIM, STOREDIM);
                    Y2 = (QUICKDouble) LOC2(store,    0, MCD-1 , STOREDIM, STOREDIM) + c * LOC2(store,    0, MC-1 , STOREDIM, STOREDIM);
                    break;
                }
            }
            for (int i = 0; i<3; i++) {
                if (NB[i] != 0) {
                    Y = Y1 + (RA[i]-RB[i])*Y2;
                    break;
                }
            }
            break;
        }
        case 1100:
        {
            int MAB = (int) LOC3(devTrans, NA[0]+NB[0], NA[1]+NB[1], NA[2]+NB[2], TRANSDIM, TRANSDIM, TRANSDIM);
            for (int i = 0; i<3; i++) {
                if (NB[i] != 0) {
                    Y = (QUICKDouble) LOC2(store, MAB-1, 0 , STOREDIM, STOREDIM) + (RA[i]-RB[i])*LOC2(store, MA-1, 0, STOREDIM, STOREDIM);
                    break;
                }
            }
            break;
        }
        case 1110:
        {
            int MAB = (int) LOC3(devTrans, NA[0]+NB[0], NA[1]+NB[1], NA[2]+NB[2], TRANSDIM, TRANSDIM, TRANSDIM);
            for (int i = 0; i<3; i++) {
                if (NB[i] != 0) {
                    Y = (QUICKDouble) LOC2(store, MAB-1, MC-1, STOREDIM, STOREDIM) + (RA[i]-RB[i])*LOC2(store, MA-1, MC-1, STOREDIM, STOREDIM);
                    break;
                }
            }
            break;
        }
        case 1101:
        {
            QUICKDouble Y1,Y2;
            int MAB = (int) LOC3(devTrans, NA[0]+NB[0], NA[1]+NB[1], NA[2]+NB[2], TRANSDIM, TRANSDIM, TRANSDIM);
            for (int i = 0; i<3; i++) {
                if (ND[i] != 0) {
                    QUICKDouble c = (QUICKDouble) (RC[i] - RD[i]);
                    Y1 = (QUICKDouble) LOC2(store, MAB-1, MD-1 , STOREDIM, STOREDIM) + c * LOC2(store, MAB-1, 0 , STOREDIM, STOREDIM);
                    Y2 = (QUICKDouble) LOC2(store,  MA-1, MD-1 , STOREDIM, STOREDIM) + c * LOC2(store,  MA-1, 0 , STOREDIM, STOREDIM);
                    break;
                }
            }
            for (int i = 0; i<3; i++) {
                if (NB[i] != 0) {
                    Y = Y1 + (RA[i] - RB[i])*Y2;
                    break;
                }
            }
            break;
        }
        case 1111:
        {
            QUICKDouble Y1,Y2;
            int MAB = (int) LOC3(devTrans, NA[0]+NB[0], NA[1]+NB[1], NA[2]+NB[2], TRANSDIM, TRANSDIM, TRANSDIM);
            int MCD = (int) LOC3(devTrans, NC[0]+ND[0], NC[1]+ND[1], NC[2]+ND[2], TRANSDIM, TRANSDIM, TRANSDIM);
            for (int i = 0; i<3; i++) {
                if (ND[i] != 0) {
                    QUICKDouble c = (QUICKDouble) (RC[i] - RD[i]);
                    Y1 = (QUICKDouble) LOC2(store, MAB-1, MCD-1 , STOREDIM, STOREDIM) + c * LOC2(store, MAB-1, MC-1 , STOREDIM, STOREDIM);
                    Y2 = (QUICKDouble) LOC2(store,  MA-1, MCD-1 , STOREDIM, STOREDIM) + c * LOC2(store,  MA-1, MC-1 , STOREDIM, STOREDIM);
                    break;
                }
            }
            
            for (int i = 0; i<3; i++) {
                if (NB[i] != 0) {
                    Y = Y1 + (RA[i] - RB[i])*Y2;
                    break;
                }
            }
            break;
        }
        case 1:
        {
            for (int i = 0; i< 3; i++) {
                if (ND[i] != 0) {
                    Y = (QUICKDouble) LOC2(store, 0, MD-1 , STOREDIM, STOREDIM) + (RC[i] - RD[i]) * LOC2(store, 0, 0, STOREDIM, STOREDIM);
                    break;
                }
            }
            break;
        }
        case 11:
        {
            int MCD = (int) LOC3(devTrans, NC[0]+ND[0], NC[1]+ND[1], NC[2]+ND[2], TRANSDIM, TRANSDIM, TRANSDIM);
            for (int i = 0; i<3; i++) {
                if (ND[i] != 0) {
                    Y = (QUICKDouble) LOC2(store, 0, MCD-1 , STOREDIM, STOREDIM) + (RC[i] - RD[i]) * LOC2(store, 0, MC-1, STOREDIM, STOREDIM);
                    break;
                }
            }
            break;
        }
        case 1001:
        {   
            for (int i = 0; i<3; i++) {
                if (ND[i] != 0) {
                    Y = (QUICKDouble) LOC2(store, MA-1, MD-1 , STOREDIM, STOREDIM) + (RC[i] - RD[i]) * LOC2(store, MA-1, 0, STOREDIM, STOREDIM);
                    break;
                }
            }

            break;
        }
        case 1011:
        {
            int MCD = (int) LOC3(devTrans, NC[0]+ND[0], NC[1]+ND[1], NC[2]+ND[2], TRANSDIM, TRANSDIM, TRANSDIM);
            for (int i = 0; i<3; i++) {
                if (ND[i] != 0) {
                    Y = (QUICKDouble) LOC2(store, MA-1, MCD-1 , STOREDIM, STOREDIM) + (RC[i] - RD[i]) * LOC2(store, MA-1, MC-1, STOREDIM, STOREDIM);
                    break;
                }
            }
            break;
        }
        default:
        {   
            int numAngularL, numAngularR;
            int angularL[20], angularR[20];
            QUICKDouble coefAngularL[20], coefAngularR[20];
            Y = (QUICKDouble) 0;
//            lefthrr
//            lefthrr
            numAngularL = 1;
            // delete the above line.
            for (int i = 0; i<numAngularL; i++) {
                for (int j = 0; j<numAngularR; j++) {
                    Y += coefAngularL[i] * coefAngularR[i] * LOC2(store, angularL[i]-1, angularR[i]-1 , STOREDIM, STOREDIM);
                }
            }
            
            Y = Y * devSim.cons[III-1] * devSim.cons[JJJ-1] * devSim.cons[KKK-1] * devSim.cons[LLL-1];
            break;
        }
    }
    return Y;
}









