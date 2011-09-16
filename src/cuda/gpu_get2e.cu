/*
 *  gpu_get2e.cpp
 *  new_quick
 *
 *  Created by Yipu Miao on 6/17/11.
 *  Copyright 2011 University of Florida.All rights reserved.
 *
 */

#include "gpu.h"
#include <cuda.h>

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
    status = cudaMemcpyToSymbol("devSim", &gpu->gpu_sim, sizeof(gpu_simulation_type), 0, cudaMemcpyHostToDevice);
    PRINTERROR(status, " cudaMemcpyToSymbol, sim copy to constants failed")
#endif
}

void upload_para_to_const(){

    int trans[TRANSDIM*TRANSDIM*TRANSDIM];
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
    
    int Mcal[3*MCALDIM];
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

#ifdef DEBUG
static float totTime;
#endif

void get2e(_gpu_type gpu)
{
//    dim3 blocks(64,64);
//    gpu->threadsPerBlock = 1;
#ifdef DEBUG
    cudaEvent_t start,end;
    cudaEventCreate(&start);
    cudaEventCreate(&end);
    cudaEventRecord(start, 0);
#endif
	get2e_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();

#ifdef DEBUG
    cudaEventRecord(end, 0);
    cudaEventSynchronize(end);
    float time;
    cudaEventElapsedTime(&time, start, end);
    totTime+=time;
    printf("this cycle:%f ms total time:%f ms\n", time, totTime);
    cudaEventDestroy(start);
    cudaEventDestroy(end);
#endif
    
}

__global__ void 
get2e_kernel()
{
    unsigned int offside = blockIdx.x*blockDim.x+threadIdx.x;
    int totalThreads = blockDim.x*gridDim.x;

    QUICKULL jshell   = (QUICKULL) devSim.sqrQshell;
    QUICKULL myInt    = (QUICKULL) jshell*jshell / totalThreads;

    if ((jshell*jshell - myInt*totalThreads)> offside) myInt++;

    for (QUICKULL i = 1; i<=myInt; i++) {

        QUICKULL currentInt = totalThreads * (i-1)+offside;        
        QUICKULL a = (QUICKULL) currentInt/jshell;
        QUICKULL b = (QUICKULL) (currentInt - a*jshell);
        
        /*
        QUICKULL a, b;
        double aa = (double)((currentInt+1)*1E-4);
        QUICKULL t = (QUICKULL)(sqrt(aa)*1E2);
        if ((currentInt+1)==t*t) {
            t--;
        }
        
        QUICKULL k = currentInt-t*t;
        if (k<=t) {
            a = k;
            b = t;
        }else {
            a = t;
            b = 2*t-k;
        }*/

        
        
        int II = devSim.sorted_YCutoffIJ[a].x;
        int JJ = devSim.sorted_YCutoffIJ[a].y;
        int KK = devSim.sorted_YCutoffIJ[b].x;
        int LL = devSim.sorted_YCutoffIJ[b].y;        
        
        int ii = devSim.sorted_Q[II];
        int jj = devSim.sorted_Q[JJ];
        int kk = devSim.sorted_Q[KK];
        int ll = devSim.sorted_Q[LL];
        
        if (ii<=kk){
            int nshell = devSim.nshell;
            QUICKDouble DNMax = MAX(MAX(4.0*LOC2(devSim.cutMatrix, ii, jj, nshell, nshell), 4.0*LOC2(devSim.cutMatrix, kk, ll, nshell, nshell)),
                                    MAX(MAX(LOC2(devSim.cutMatrix, ii, ll, nshell, nshell),     LOC2(devSim.cutMatrix, ii, kk, nshell, nshell)),
                                        MAX(LOC2(devSim.cutMatrix, jj, kk, nshell, nshell),     LOC2(devSim.cutMatrix, jj, ll, nshell, nshell))));
            
            if ((LOC2(devSim.YCutoff, kk, ll, nshell, nshell) * LOC2(devSim.YCutoff, ii, jj, nshell, nshell))> devSim.integralCutoff && \
                (LOC2(devSim.YCutoff, kk, ll, nshell, nshell) * LOC2(devSim.YCutoff, ii, jj, nshell, nshell) * DNMax) > devSim.integralCutoff) {
                
                int iii = devSim.sorted_Qnumber[II];
                int jjj = devSim.sorted_Qnumber[JJ];
                int kkk = devSim.sorted_Qnumber[KK];
                int lll = devSim.sorted_Qnumber[LL];
                
                iclass(iii, jjj, kkk, lll, ii, jj, kk, ll, DNMax);
                
            }
        }
    }
}

__device__ QUICKDouble quick_dsqr(QUICKDouble a)
{
    return a*a;
}



#ifndef TEST
__device__
#endif
void iclass(int I, int J, int K, int L, unsigned int II, unsigned int JJ, unsigned int KK, unsigned int LL, QUICKDouble DNMax)
{

    QUICKDouble RAx, RAy, RAz;
    QUICKDouble RBx, RBy, RBz;
    QUICKDouble RCx, RCy, RCz;
    QUICKDouble RDx, RDy, RDz;

    int kPrimI, kPrimJ, kPrimL, kPrimK;
    int kStartI, kStartJ, kStartK, kStartL;

    /* 
     kAtom A, B, C ,D is the coresponding atom for shell ii, jj, kk, ll
     and be careful with the index difference between Fortran and C++, 
     Fortran starts array index with 1 and C++ starts 0.
     */
    /*
     NII1 is the starting angular momenta for shell i and NII2 is the ending
     angular momenta.So it is with other varibles
     */
     
    /*RA, RB, RC, and RD are the coordinates for atom katomA, katomB, katomC and katomD, 
     which means they are corrosponding coorinates for shell II, JJ, KK, and LL.
     */
    RAx = LOC2(devSim.xyz, 0 , devSim.katom[II]-1, 3, devSim.natom);
    RAy = LOC2(devSim.xyz, 1 , devSim.katom[II]-1, 3, devSim.natom);
    RAz = LOC2(devSim.xyz, 2 , devSim.katom[II]-1, 3, devSim.natom);
    
    RCx = LOC2(devSim.xyz, 0 , devSim.katom[KK]-1, 3, devSim.natom);
    RCy = LOC2(devSim.xyz, 1 , devSim.katom[KK]-1, 3, devSim.natom);
    RCz = LOC2(devSim.xyz, 2 , devSim.katom[KK]-1, 3, devSim.natom);
    
    
    kPrimI = devSim.kprim[II];
    kPrimJ = devSim.kprim[JJ];
    kPrimK = devSim.kprim[KK];
    kPrimL = devSim.kprim[LL];
    
    kStartI = devSim.kstart[II];
    kStartJ = devSim.kstart[JJ];
    kStartK = devSim.kstart[KK];
    kStartL = devSim.kstart[LL];
    
    
    
    QUICKDouble store[STOREDIM*STOREDIM];
	memset(store, 0, STOREDIM*STOREDIM*sizeof(QUICKDouble));
    
    for (int JJJ = 0; JJJ < kPrimJ; JJJ++) {
        for (int III = 0; III < kPrimI; III++) {
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
            QUICKDouble AB = LOC4(devSim.expoSum, III, JJJ, II, JJ, 6, 6, devSim.jshell, devSim.jshell);
            
			/*
                              --->                --->
             ->     expo(I) * xyz (I) + expo(J) * xyz(J)
             P  = ---------------------------------------
                              expo(I) + expo(J)
             
                                    -->             -->
             ----->        expo(I)*xyz(I) + expo(J)*xyz(J)                                 -->            -->
             AAtemp = ----------------------------------- * (expo(I) + expo(J)) = expo(I)*xyz(I)+expo(J)*xyz(J)
                                  expo(I) + expo(J)
             
             ----->   ->  ->
             Ptemp  = P - A
             */            
             
             QUICKDouble cutoffPrim = DNMax * LOC2(devSim.cutPrim, kStartI+III-1, kStartJ+JJJ-1, devSim.jbasis, devSim.jbasis);
            
            QUICKDouble X1 = LOC4(devSim.Xcoeff, kStartI+III-1, kStartJ+JJJ-1, I, J, devSim.jbasis, devSim.jbasis, 4, 4);
            
            QUICKDouble Px = LOC4(devSim.weightedCenterX, III, JJJ, II, JJ, 6, 6, devSim.jshell, devSim.jshell);
            QUICKDouble Py = LOC4(devSim.weightedCenterY, III, JJJ, II, JJ, 6, 6, devSim.jshell, devSim.jshell);
            QUICKDouble Pz = LOC4(devSim.weightedCenterZ, III, JJJ, II, JJ, 6, 6, devSim.jshell, devSim.jshell);
            
                
            for (int LLL = 0 ; LLL < kPrimL; LLL++) {
                for (int KKK = 0; KKK < kPrimK; KKK++) {
                    if (cutoffPrim * LOC2(devSim.cutPrim, kStartK+KKK-1, kStartL+LLL-1, devSim.jbasis, devSim.jbasis) > devSim.primLimit) {
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

                                      expo(I)+expo(J)                        expo(K)+expo(L)
                         ABcom = --------------------------------  CDcom = --------------------------------
                                  expo(I)+expo(J)+expo(K)+expo(L)           expo(I)+expo(J)+expo(K)+expo(L)
                         
                         ABCDtemp = 1/2(expo(I)+expo(J)+expo(K)+expo(L))                    
                         */                        
                        QUICKDouble CD = LOC4(devSim.expoSum, KKK, LLL, KK, LL, 6, 6, devSim.jshell, devSim.jshell);
                        QUICKDouble ABCD = 1/(AB+CD);

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
                         
                        ---->   ->  ->
                         Qtemp = Q - K
                         ----->   ->  ->
                         WQtemp = W - Q
                         ----->   ->  ->
                         WPtemp = W - P

                         ->  -> 2
                         T = ROU * | P - Q|
                         */
                         
                        QUICKDouble Qx = LOC4(devSim.weightedCenterX, KKK, LLL, KK, LL, 6, 6, devSim.jshell, devSim.jshell);
                        QUICKDouble Qy = LOC4(devSim.weightedCenterY, KKK, LLL, KK, LL, 6, 6, devSim.jshell, devSim.jshell);
                        QUICKDouble Qz = LOC4(devSim.weightedCenterZ, KKK, LLL, KK, LL, 6, 6, devSim.jshell, devSim.jshell);
                                                                        
                        QUICKDouble T = AB * CD * ABCD * ( quick_dsqr(Px-Qx) + quick_dsqr(Py-Qy) + quick_dsqr(Pz-Qz));

                        QUICKDouble YVerticalTemp[VDIM1*VDIM2*VDIM3];
                        FmT(I+J+K+L, T, YVerticalTemp);
                        QUICKDouble X2 = sqrt(ABCD) * X1 * LOC4(devSim.Xcoeff, kStartK+KKK-1, kStartL+LLL-1, K, L, devSim.jbasis, devSim.jbasis, 4, 4);
                        
                        for (int i = 0; i<=I+J+K+L; i++) {
                            VY(0, 0, i) = VY(0, 0, i) * X2;
                        }
                        
                        /*
                         X2 is the multiplication of four indices normalized coeffecient
                         */
                        
                       /* 
                        if (I+J+K+L != 0) {
                       
                            QUICKDouble tempx = (Px*AB+Qx*CD)*ABCD;
                            QUICKDouble tempy = (Py*AB+Qy*CD)*ABCD;
                            QUICKDouble tempz = (Pz*AB+Qz*CD)*ABCD;
                       
                            ABCD = ABCD * 0.5;
                            QUICKDouble ABtemp, CDtemp;
                            
                            if (I+J>0) {
                                //PSSS(0)
                                VY( 1, 0, 0) = (Px-RAx) * VY( 0, 0, 0) + (tempx - Px) * VY( 0, 0, 1);
                                VY( 2, 0, 0) = (Py-RAy) * VY( 0, 0, 0) + (tempy - Py) * VY( 0, 0, 1);
                                VY( 3, 0, 0) = (Pz-RAz) * VY( 0, 0, 0) + (tempz - Pz) * VY( 0, 0, 1);
                            }
                            
                            if (K+L>0) {
                                //SSPS(0)
                                VY( 0, 1, 0) = (Qx-RCx) * VY( 0, 0, 0) + (tempx - Qx) * VY( 0, 0, 1);
                                VY( 0, 2, 0) = (Qy-RCy) * VY( 0, 0, 0) + (tempy - Qy) * VY( 0, 0, 1);
                                VY( 0, 3, 0) = (Qz-RCz) * VY( 0, 0, 0) + (tempz - Qz) * VY( 0, 0, 1);
                            }
                            
                            if ((I+J>0 && K+L>0) || K+L>1) {
                                //SSPS(1, YVerticalTemp, Qx-RCx, Qy-RCy, Qz-RCz, tempx - Qx, tempy - Qy, tempz - Qz);
                                VY( 0, 1, 1) = (Qx-RCx) * VY( 0, 0, 1) + (tempx - Qx) * VY( 0, 0, 2);
                                VY( 0, 2, 1) = (Qy-RCy) * VY( 0, 0, 1) + (tempy - Qy) * VY( 0, 0, 2);
                                VY( 0, 3, 1) = (Qz-RCz) * VY( 0, 0, 1) + (tempz - Qz) * VY( 0, 0, 2);
                            }
                            
                            if (I+J>0 && K+L>0) {
                                
                                //PSPS(0, YVerticalTemp, Px-RAx, Py-RAy, Pz-RAz, tempx - Px, tempy - Py, tempz - Pz, 0.5*ABCD);
                                VY( 1, 1, 0) = (Px-RAx) * VY( 0, 1, 0) + (tempx - Px) * VY( 0, 1, 1) + ABCD * VY( 0, 0, 1);
                                VY( 2, 1, 0) = (Py-RAy) * VY( 0, 1, 0) + (tempy - Py) * VY( 0, 1, 1);
                                VY( 3, 1, 0) = (Pz-RAz) * VY( 0, 1, 0) + (tempz - Pz) * VY( 0, 1, 1);

                                VY( 1, 2, 0) = (Px-RAx) * VY( 0, 2, 0) + (tempx - Px) * VY( 0, 2, 1);
                                VY( 2, 2, 0) = (Py-RAy) * VY( 0, 2, 0) + (tempy - Py) * VY( 0, 2, 1) + ABCD * VY( 0, 0, 1);
                                VY( 3, 2, 0) = (Pz-RAz) * VY( 0, 2, 0) + (tempz - Pz) * VY( 0, 2, 1);

                                VY( 1, 3, 0) = (Px-RAx) * VY( 0, 3, 0) + (tempx - Px) * VY( 0, 3, 1);
                                VY( 2, 3, 0) = (Py-RAy) * VY( 0, 3, 0) + (tempy - Py) * VY( 0, 3, 1);
                                VY( 3, 3, 0) = (Pz-RAz) * VY( 0, 3, 0) + (tempz - Pz) * VY( 0, 3, 1) + ABCD * VY( 0, 0, 1);
                            }
                            
                            if (I+J>1) {
                                ABtemp = 0.5/AB;
                                //PSSS(1, YVerticalTemp, Px-RAx, Py-RAy, Pz-RAz, tempx - Px, tempy - Py, tempz - Pz);
                                VY( 1, 0, 1) = (Px-RAx) * VY( 0, 0, 1) + (tempx - Px) * VY( 0, 0, 2);
                                VY( 2, 0, 1) = (Py-RAy) * VY( 0, 0, 1) + (tempy - Py) * VY( 0, 0, 2);
                                VY( 3, 0, 1) = (Pz-RAz) * VY( 0, 0, 1) + (tempz - Pz) * VY( 0, 0, 2);
                                
                                //DSSS(0, YVerticalTemp, Px-RAx, Py-RAy, Pz-RAz, tempx - Px, tempy - Py, tempz - Pz, 0.5/AB, CD*ABCD);
                                VY( 4, 0, 0) = (Px-RAx) * VY( 2, 0, 0) + (tempx - Px) * VY( 2, 0, 1);
                                VY( 5, 0, 0) = (Py-RAy) * VY( 3, 0, 0) + (tempy - Py) * VY( 3, 0, 1);
                                VY( 6, 0, 0) = (Px-RAx) * VY( 3, 0, 0) + (tempx - Px) * VY( 3, 0, 1);
                                
                                VY( 7, 0, 0) = (Px-RAx) * VY( 1, 0, 0) + (tempx - Px) * VY( 1, 0, 1) + ABtemp *(VY( 0, 0, 0) - CD * 2 * ABCD * VY( 0, 0, 1));
                                VY( 8, 0, 0) = (Py-RAy) * VY( 2, 0, 0) + (tempy - Py) * VY( 2, 0, 1) + ABtemp *(VY( 0, 0, 0) - CD * 2 * ABCD * VY( 0, 0, 1));
                                VY( 9, 0, 0) = (Pz-RAz) * VY( 3, 0, 0) + (tempz - Pz) * VY( 3, 0, 1) + ABtemp *(VY( 0, 0, 0) - CD * 2 * ABCD * VY( 0, 0, 1));
                            }
                            
                            if (K+L>1) {
                                CDtemp = 1/CD;
                                //SSPS(1, YVerticalTemp, Qx-RCx, Qy-RCy, Qz-RCz, tempx - Qx, tempy - Qy, tempz - Qz);
                                VY( 0, 1, 1) = (Qx-RCx) * VY( 0, 0, 1) + (tempx - Qx) * VY( 0, 0, 2);
                                VY( 0, 2, 1) = (Qy-RCy) * VY( 0, 0, 1) + (tempy - Qy) * VY( 0, 0, 2);
                                VY( 0, 3, 1) = (Qz-RCz) * VY( 0, 0, 1) + (tempz - Qz) * VY( 0, 0, 2);
                                
                                //SSDS(0, YVerticalTemp, Qx-RCx, Qy-RCy, Qz-RCz, tempx - Qx, tempy - Qy, tempz - Qz, 0.5/CD, AB*ABCD);
                                VY( 0, 4, 0) = (Qx-RCx) * VY( 0, 2, 0) + (tempx - Qx) * VY( 0, 2, 1);
                                VY( 0, 5, 0) = (Qy-RCy) * VY( 0, 3, 0) + (tempy - Qy) * VY( 0, 3, 1);
                                VY( 0, 6, 0) = (Qx-RCx) * VY( 0, 3, 0) + (tempx - Qx) * VY( 0, 3, 1);
                                
                                VY( 0, 7, 0) = (Qx-RCx) * VY( 0, 1, 0) + (tempx - Qx) * VY( 0, 1, 1) + 0.5 * CDtemp *(VY( 0, 0, 0) - AB*2 * ABCD * VY( 0, 0, 1)) ;
                                VY( 0, 8, 0) = (Qy-RCy) * VY( 0, 2, 0) + (tempy - Qy) * VY( 0, 2, 1) + 0.5 * CDtemp *(VY( 0, 0, 0) - AB*2 * ABCD * VY( 0, 0, 1)) ;
                                VY( 0, 9, 0) = (Qz-RCz) * VY( 0, 3, 0) + (tempz - Qz) * VY( 0, 3, 1) + 0.5 * CDtemp *(VY( 0, 0, 0) - AB*2 * ABCD * VY( 0, 0, 1)) ;
                            }
                            
                            if (I+J>1 && K+L>0) {
                                //PSSS(2, YVerticalTemp, Px-RAx, Py-RAy, Pz-RAz, tempx - Px, tempy - Py, tempz - Pz);
                                VY( 1, 0, 2) = (Px-RAx) * VY( 0, 0, 2) + (tempx - Px) * VY( 0, 0, 3);
                                VY( 2, 0, 2) = (Py-RAy) * VY( 0, 0, 2) + (tempy - Py) * VY( 0, 0, 3);
                                VY( 3, 0, 2) = (Pz-RAz) * VY( 0, 0, 2) + (tempz - Pz) * VY( 0, 0, 3);
                                
                                //DSSS(1, YVerticalTemp, Px-RAx, Py-RAy, Pz-RAz, tempx - Px, tempy - Py, tempz - Pz, 0.5/AB, CD*ABCD);
                                VY( 4, 0, 1) = (Px-RAx) * VY( 2, 0, 1) + (tempx - Px) * VY( 2, 0, 2);
                                VY( 5, 0, 1) = (Py-RAy) * VY( 3, 0, 1) + (tempy - Py) * VY( 3, 0, 2);
                                VY( 6, 0, 1) = (Px-RAx) * VY( 3, 0, 1) + (tempx - Px) * VY( 3, 0, 2);
                                
                                VY( 7, 0, 1) = (Px-RAx) * VY( 1, 0, 1) + (tempx - Px) * VY( 1, 0, 2) + ABtemp *( VY( 0, 0, 1) - 2 * CD * ABCD * VY( 0, 0, 2));
                                VY( 8, 0, 1) = (Py-RAy) * VY( 2, 0, 1) + (tempy - Py) * VY( 2, 0, 2) + ABtemp *( VY( 0, 0, 1) - 2 * CD * ABCD * VY( 0, 0, 2));
                                VY( 9, 0, 1) = (Pz-RAz) * VY( 3, 0, 1) + (tempz - Pz) * VY( 3, 0, 2) + ABtemp *( VY( 0, 0, 1) - 2 * CD * ABCD * VY( 0, 0, 2));
                                
                                //DSPS(0, YVerticalTemp, Qx-RCx, Qy-RCy, Qz-RCz, tempx - Qx, tempy - Qy, tempz - Qz, 0.5*ABCD);
                                VY( 4, 1, 0) = (Qx-RCx) * VY( 4, 0, 0) + (tempx - Qx) * VY( 4, 0, 1) + ABCD * VY( 2, 0, 1);
                                VY( 4, 2, 0) = (Qy-RCy) * VY( 4, 0, 0) + (tempy - Qy) * VY( 4, 0, 1) + ABCD * VY( 1, 0, 1);
                                VY( 4, 3, 0) = (Qz-RCz) * VY( 4, 0, 0) + (tempz - Qz) * VY( 4, 0, 1);
                                
                                VY( 5, 1, 0) = (Qx-RCx) * VY( 5, 0, 0) + (tempx - Qx) * VY( 5, 0, 1);
                                VY( 5, 2, 0) = (Qy-RCy) * VY( 5, 0, 0) + (tempy - Qy) * VY( 5, 0, 1) + ABCD * VY( 3, 0, 1);
                                VY( 5, 3, 0) = (Qz-RCz) * VY( 5, 0, 0) + (tempz - Qz) * VY( 5, 0, 1) + ABCD * VY( 2, 0, 1);
                                
                                VY( 6, 1, 0) = (Qx-RCx) * VY( 6, 0, 0) + (tempx - Qx) * VY( 6, 0, 1) + ABCD * VY( 3, 0, 1);
                                VY( 6, 2, 0) = (Qy-RCy) * VY( 6, 0, 0) + (tempy - Qy) * VY( 6, 0, 1);
                                VY( 6, 3, 0) = (Qz-RCz) * VY( 6, 0, 0) + (tempz - Qz) * VY( 6, 0, 1) + ABCD * VY( 1, 0, 1);
                                
                                VY( 7, 1, 0) = (Qx-RCx) * VY( 7, 0, 0) + (tempx - Qx) * VY( 7, 0, 1) + 2 * ABCD * VY( 1, 0, 1);
                                VY( 7, 2, 0) = (Qy-RCy) * VY( 7, 0, 0) + (tempy - Qy) * VY( 7, 0, 1);
                                VY( 7, 3, 0) = (Qz-RCz) * VY( 7, 0, 0) + (tempz - Qz) * VY( 7, 0, 1);
                                
                                VY( 8, 1, 0) = (Qx-RCx) * VY( 8, 0, 0) + (tempx - Qx) * VY( 8, 0, 1);
                                VY( 8, 2, 0) = (Qy-RCy) * VY( 8, 0, 0) + (tempy - Qy) * VY( 8, 0, 1) + 2 * ABCD * VY( 2, 0, 1);
                                VY( 8, 3, 0) = (Qz-RCz) * VY( 8, 0, 0) + (tempz - Qz) * VY( 8, 0, 1);
                                
                                VY( 9, 1, 0) = (Qx-RCx) * VY( 9, 0, 0) + (tempx - Qx) * VY( 9, 0, 1);
                                VY( 9, 2, 0) = (Qy-RCy) * VY( 9, 0, 0) + (tempy - Qy) * VY( 9, 0, 1);
                                VY( 9, 3, 0) = (Qz-RCz) * VY( 9, 0, 0) + (tempz - Qz) * VY( 9, 0, 1) + 2 * ABCD * VY( 3, 0, 1);
                            }
                            
                            if (I+J>0 && K+L>1) {
                            
                                //SSPS(2, YVerticalTemp, Qx-RCx, Qy-RCy, Qz-RCz, tempx - Qx, tempy - Qy, tempz - Qz);
                                VY( 0, 1, 2) = (Qx-RCx) * VY( 0, 0, 2) + (tempx - Qx) * VY( 0, 0, 3);
                                VY( 0, 2, 2) = (Qy-RCy) * VY( 0, 0, 2) + (tempy - Qy) * VY( 0, 0, 3);
                                VY( 0, 3, 2) = (Qz-RCz) * VY( 0, 0, 2) + (tempz - Qz) * VY( 0, 0, 3);
                                
                                //SSDS(1, YVerticalTemp, Qx-RCx, Qy-RCy, Qz-RCz, tempx - Qx, tempy - Qy, tempz - Qz, 0.5/CD, AB*ABCD);
                                VY( 0, 4, 1) = (Qx-RCx) * VY( 0, 2, 1) + (tempx - Qx) * VY( 0, 2, 2);
                                VY( 0, 5, 1) = (Qy-RCy) * VY( 0, 3, 1) + (tempy - Qy) * VY( 0, 3, 2);
                                VY( 0, 6, 1) = (Qx-RCx) * VY( 0, 3, 1) + (tempx - Qx) * VY( 0, 3, 2);
                                
                                VY( 0, 7, 1) = (Qx-RCx) * VY( 0, 1, 1) + (tempx - Qx) * VY( 0, 1, 2) + 1 * CDtemp *( 0.5 * VY( 0, 0, 1) - AB*ABCD * VY( 0, 0, 2)) ;
                                VY( 0, 8, 1) = (Qy-RCy) * VY( 0, 2, 1) + (tempy - Qy) * VY( 0, 2, 2) + 1 * CDtemp *( 0.5 * VY( 0, 0, 1) - AB*ABCD * VY( 0, 0, 2)) ;
                                VY( 0, 9, 1) = (Qz-RCz) * VY( 0, 3, 1) + (tempz - Qz) * VY( 0, 3, 2) + 1 * CDtemp *( 0.5 * VY( 0, 0, 1) - AB*ABCD * VY( 0, 0, 2)) ;
                                
                                //PSDS(0, YVerticalTemp, Px-RAx, Py-RAy, Pz-RAz, tempx - Px, tempy - Py, tempz - Pz, 0.5*ABCD);
                                VY( 1, 4, 0) = (Px-RAx) * VY( 0, 4, 0) + (tempx - Px) * VY( 0, 4, 1) + ABCD * VY( 0, 2, 1);
                                VY( 2, 4, 0) = (Py-RAy) * VY( 0, 4, 0) + (tempy - Py) * VY( 0, 4, 1) + ABCD * VY( 0, 1, 1);
                                VY( 3, 4, 0) = (Pz-RAz) * VY( 0, 4, 0) + (tempz - Pz) * VY( 0, 4, 1) ;
                                
                                VY( 1, 5, 0) = (Px-RAx) * VY( 0, 5, 0) + (tempx - Px) * VY( 0, 5, 1);
                                VY( 2, 5, 0) = (Py-RAy) * VY( 0, 5, 0) + (tempy - Py) * VY( 0, 5, 1) + ABCD * VY( 0, 3, 1);
                                VY( 3, 5, 0) = (Pz-RAz) * VY( 0, 5, 0) + (tempz - Pz) * VY( 0, 5, 1) + ABCD * VY( 0, 2, 1);
                                
                                VY( 1, 6, 0) = (Px-RAx) * VY( 0, 6, 0) + (tempx - Px) * VY( 0, 6, 1) + ABCD * VY( 0, 3, 1);
                                VY( 2, 6, 0) = (Py-RAy) * VY( 0, 6, 0) + (tempy - Py) * VY( 0, 6, 1);
                                VY( 3, 6, 0) = (Pz-RAz) * VY( 0, 6, 0) + (tempz - Pz) * VY( 0, 6, 1) + ABCD * VY( 0, 1, 1);
                                
                                VY( 1, 7, 0) = (Px-RAx) * VY( 0, 7, 0) + (tempx - Px) * VY( 0, 7, 1) + 2 * ABCD * VY( 0, 1, 1);
                                VY( 2, 7, 0) = (Py-RAy) * VY( 0, 7, 0) + (tempy - Py) * VY( 0, 7, 1);
                                VY( 3, 7, 0) = (Pz-RAz) * VY( 0, 7, 0) + (tempz - Pz) * VY( 0, 7, 1);

                                VY( 1, 8, 0) = (Px-RAx) * VY( 0, 8, 0) + (tempx - Px) * VY( 0, 8, 1);
                                VY( 2, 8, 0) = (Py-RAy) * VY( 0, 8, 0) + (tempy - Py) * VY( 0, 8, 1) + 2 * ABCD * VY( 0, 2, 1);
                                VY( 3, 8, 0) = (Pz-RAz) * VY( 0, 8, 0) + (tempz - Pz) * VY( 0, 8, 1);
                                
                                VY( 1, 9, 0) = (Px-RAx) * VY( 0, 9, 0) + (tempx - Px) * VY( 0, 9, 1);
                                VY( 2, 9, 0) = (Py-RAy) * VY( 0, 9, 0) + (tempy - Py) * VY( 0, 9, 1);
                                VY( 3, 9, 0) = (Pz-RAz) * VY( 0, 9, 0) + (tempz - Pz) * VY( 0, 9, 1) + 2 * ABCD * VY( 0, 3, 1);                                

                            }

                            if (I+J>1 && K+L>1) {
                                //SSPS(3, YVerticalTemp, Qx-RCx, Qy-RCy, Qz-RCz, tempx - Qx, tempy - Qy, tempz - Qz);
                                VY( 0, 1, 3) = (Qx-RCx) * VY( 0, 0, 3) + (tempx - Qx) * VY( 0, 0, 4);
                                VY( 0, 2, 3) = (Qy-RCy) * VY( 0, 0, 3) + (tempy - Qy) * VY( 0, 0, 4);
                                VY( 0, 3, 3) = (Qz-RCz) * VY( 0, 0, 3) + (tempz - Qz) * VY( 0, 0, 4);
                                
                                //SSDS(2, YVerticalTemp, Qx-RCx, Qy-RCy, Qz-RCz, tempx - Qx, tempy - Qy, tempz - Qz, 0.5/CD, AB*ABCD);
                                VY( 0, 4, 2) = (Qx-RCx) * VY( 0, 2, 2) + (tempx - Qx) * VY( 0, 2, 3);
                                VY( 0, 5, 2) = (Qy-RCy) * VY( 0, 3, 2) + (tempy - Qy) * VY( 0, 3, 3);
                                VY( 0, 6, 2) = (Qx-RCx) * VY( 0, 3, 2) + (tempx - Qx) * VY( 0, 3, 3);
                                
                                VY( 0, 7, 2) = (Qx-RCx) * VY( 0, 1, 2) + (tempx - Qx) * VY( 0, 1, 3) + CDtemp *( 0.5 * VY( 0, 0, 2) - AB*ABCD * VY( 0, 0, 3)) ;
                                VY( 0, 8, 2) = (Qy-RCy) * VY( 0, 2, 2) + (tempy - Qy) * VY( 0, 2, 3) + CDtemp *( 0.5 * VY( 0, 0, 2) - AB*ABCD * VY( 0, 0, 3)) ;
                                VY( 0, 9, 2) = (Qz-RCz) * VY( 0, 3, 2) + (tempz - Qz) * VY( 0, 3, 3) + CDtemp *( 0.5 * VY( 0, 0, 2) - AB*ABCD * VY( 0, 0, 3)) ;
                                
                                //PSDS(1, YVerticalTemp, Px-RAx, Py-RAy, Pz-RAz, tempx - Px, tempy - Py, tempz - Pz, 0.5*ABCD);
                                VY( 1, 4, 1) = (Px-RAx) * VY( 0, 4, 1) + (tempx - Px) * VY( 0, 4, 2) + ABCD * VY( 0, 2, 2);
                                VY( 2, 4, 1) = (Py-RAy) * VY( 0, 4, 1) + (tempy - Py) * VY( 0, 4, 2) + ABCD * VY( 0, 1, 2);
                                VY( 3, 4, 1) = (Pz-RAz) * VY( 0, 4, 1) + (tempz - Pz) * VY( 0, 4, 2) ;
                                
                                VY( 1, 5, 1) = (Px-RAx) * VY( 0, 5, 1) + (tempx - Px) * VY( 0, 5, 2);
                                VY( 2, 5, 1) = (Py-RAy) * VY( 0, 5, 1) + (tempy - Py) * VY( 0, 5, 2) + ABCD * VY( 0, 3, 2);
                                VY( 3, 5, 1) = (Pz-RAz) * VY( 0, 5, 1) + (tempz - Pz) * VY( 0, 5, 2) + ABCD * VY( 0, 2, 2);
                                
                                VY( 1, 6, 1) = (Px-RAx) * VY( 0, 6, 1) + (tempx - Px) * VY( 0, 6, 2) + ABCD * VY( 0, 3, 2);
                                VY( 2, 6, 1) = (Py-RAy) * VY( 0, 6, 1) + (tempy - Py) * VY( 0, 6, 2);
                                VY( 3, 6, 1) = (Pz-RAz) * VY( 0, 6, 1) + (tempz - Pz) * VY( 0, 6, 2) + ABCD * VY( 0, 1, 2);
                                
                                VY( 1, 7, 1) = (Px-RAx) * VY( 0, 7, 1) + (tempx - Px) * VY( 0, 7, 2) + 2 * ABCD * VY( 0, 1, 2);
                                VY( 2, 7, 1) = (Py-RAy) * VY( 0, 7, 1) + (tempy - Py) * VY( 0, 7, 2);
                                VY( 3, 7, 1) = (Pz-RAz) * VY( 0, 7, 1) + (tempz - Pz) * VY( 0, 7, 2);

                                VY( 1, 8, 1) = (Px-RAx) * VY( 0, 8, 1) + (tempx - Px) * VY( 0, 8, 2);
                                VY( 2, 8, 1) = (Py-RAy) * VY( 0, 8, 1) + (tempy - Py) * VY( 0, 8, 2) + 2 * ABCD * VY( 0, 2, 2);
                                VY( 3, 8, 1) = (Pz-RAz) * VY( 0, 8, 1) + (tempz - Pz) * VY( 0, 8, 2);
                                
                                VY( 1, 9, 1) = (Px-RAx) * VY( 0, 9, 1) + (tempx - Px) * VY( 0, 9, 2);
                                VY( 2, 9, 1) = (Py-RAy) * VY( 0, 9, 1) + (tempy - Py) * VY( 0, 9, 2);
                                VY( 3, 9, 1) = (Pz-RAz) * VY( 0, 9, 1) + (tempz - Pz) * VY( 0, 9, 2) + 2 * ABCD * VY( 0, 3, 2);
                                
                                //PSPS(1, YVerticalTemp, Px-RAx, Py-RAy, Pz-RAz, tempx - Px, tempy - Py, tempz - Pz, 0.5*ABCD);
                                VY( 1, 1, 1) = (Px-RAx) * VY( 0, 1, 1) + (tempx - Px) * VY( 0, 1, 2) + ABCD * VY( 0, 0, 2);
                                VY( 2, 1, 1) = (Py-RAy) * VY( 0, 1, 1) + (tempy - Py) * VY( 0, 1, 2);
                                VY( 3, 1, 1) = (Pz-RAz) * VY( 0, 1, 1) + (tempz - Pz) * VY( 0, 1, 2);

                                VY( 1, 2, 1) = (Px-RAx) * VY( 0, 2, 1) + (tempx - Px) * VY( 0, 2, 2);
                                VY( 2, 2, 1) = (Py-RAy) * VY( 0, 2, 1) + (tempy - Py) * VY( 0, 2, 2) + ABCD * VY( 0, 0, 2);
                                VY( 3, 2, 1) = (Pz-RAz) * VY( 0, 2, 1) + (tempz - Pz) * VY( 0, 2, 2);

                                VY( 1, 3, 1) = (Px-RAx) * VY( 0, 3, 1) + (tempx - Px) * VY( 0, 3, 2);
                                VY( 2, 3, 1) = (Py-RAy) * VY( 0, 3, 1) + (tempy - Py) * VY( 0, 3, 2);
                                VY( 3, 3, 1) = (Pz-RAz) * VY( 0, 3, 1) + (tempz - Pz) * VY( 0, 3, 2) + ABCD * VY( 0, 0, 2);
                                
                                //DSDS(0, YVerticalTemp, Px-RAx, Py-RAy, Pz-RAz, tempx - Px, tempy - Py, tempz - Pz, 0.5*ABCD, 0.5/AB, CD*ABCD); 
                                VY( 4, 4, 0) = (Px-RAx) * VY( 2, 4, 0) + (tempx - Px) * VY( 2, 4, 1) + ABCD * VY( 2, 2, 1);
                                VY( 4, 5, 0) = (Px-RAx) * VY( 2, 5, 0) + (tempx - Px) * VY( 2, 5, 1);
                                VY( 4, 6, 0) = (Px-RAx) * VY( 2, 6, 0) + (tempx - Px) * VY( 2, 6, 1) + ABCD * VY( 2, 3, 1);
                                VY( 4, 7, 0) = (Px-RAx) * VY( 2, 7, 0) + (tempx - Px) * VY( 2, 7, 1) + 2 * ABCD * VY( 2, 1, 1);
                                VY( 4, 8, 0) = (Px-RAx) * VY( 2, 8, 0) + (tempx - Px) * VY( 2, 8, 1);
                                VY( 4, 9, 0) = (Px-RAx) * VY( 2, 9, 0) + (tempx - Px) * VY( 2, 9, 1);
                                
                                VY( 5, 4, 0) = (Py-RAy) * VY( 3, 4, 0) + (tempy - Py) * VY( 3, 4, 1) + ABCD * VY( 3, 1, 1);
                                VY( 5, 5, 0) = (Py-RAy) * VY( 3, 5, 0) + (tempy - Py) * VY( 3, 5, 1) + ABCD * VY( 3, 3, 1);
                                VY( 5, 6, 0) = (Py-RAy) * VY( 3, 6, 0) + (tempy - Py) * VY( 3, 6, 1);
                                VY( 5, 7, 0) = (Py-RAy) * VY( 3, 7, 0) + (tempy - Py) * VY( 3, 7, 1);
                                VY( 5, 8, 0) = (Py-RAy) * VY( 3, 8, 0) + (tempy - Py) * VY( 3, 8, 1) + 2 * ABCD * VY( 3, 2, 1);
                                VY( 5, 9, 0) = (Py-RAy) * VY( 3, 9, 0) + (tempy - Py) * VY( 3, 9, 1);
                                
                                VY( 6, 4, 0) = (Px-RAx) * VY( 3, 4, 0) + (tempx - Px) * VY( 3, 4, 1) + ABCD * VY( 3, 2, 1);
                                VY( 6, 5, 0) = (Px-RAx) * VY( 3, 5, 0) + (tempx - Px) * VY( 3, 5, 1);
                                VY( 6, 6, 0) = (Px-RAx) * VY( 3, 6, 0) + (tempx - Px) * VY( 3, 6, 1) + ABCD * VY( 3, 3, 1);
                                VY( 6, 7, 0) = (Px-RAx) * VY( 3, 7, 0) + (tempx - Px) * VY( 3, 7, 1) + 2 * ABCD * VY( 3, 1, 1);
                                VY( 6, 8, 0) = (Px-RAx) * VY( 3, 8, 0) + (tempx - Px) * VY( 3, 8, 1);
                                VY( 6, 9, 0) = (Px-RAx) * VY( 3, 9, 0) + (tempx - Px) * VY( 3, 9, 1);
                                
                                VY( 7, 4, 0) = (Px-RAx) * VY( 1, 4, 0) + (tempx - Px) * VY( 1, 4, 1) +  ABtemp * (VY( 0, 4,0)-2 * CD*ABCD*VY( 0, 4,1)) + ABCD * VY( 1, 2, 1);
                                VY( 7, 5, 0) = (Px-RAx) * VY( 1, 5, 0) + (tempx - Px) * VY( 1, 5, 1) +  ABtemp * (VY( 0, 5,0)-2 * CD*ABCD*VY( 0, 5,1));
                                VY( 7, 6, 0) = (Px-RAx) * VY( 1, 6, 0) + (tempx - Px) * VY( 1, 6, 1) +  ABtemp * (VY( 0, 6,0)-2 * CD*ABCD*VY( 0, 6,1)) + ABCD * VY( 1, 3, 1);
                                VY( 7, 7, 0) = (Px-RAx) * VY( 1, 7, 0) + (tempx - Px) * VY( 1, 7, 1) +  ABtemp * (VY( 0, 7,0)-2 * CD*ABCD*VY( 0, 7,1)) + 2*ABCD * VY( 1, 1, 1);
                                VY( 7, 8, 0) = (Px-RAx) * VY( 1, 8, 0) + (tempx - Px) * VY( 1, 8, 1) +  ABtemp * (VY( 0, 8,0)-2 * CD*ABCD*VY( 0, 8,1));
                                VY( 7, 9, 0) = (Px-RAx) * VY( 1, 9, 0) + (tempx - Px) * VY( 1, 9, 1) +  ABtemp * (VY( 0, 9,0)-2 * CD*ABCD*VY( 0, 9,1));
                                
                                
                                VY( 8, 4, 0) = (Py-RAy) * VY( 2, 4, 0) + (tempy - Py) * VY( 2, 4, 1) +  ABtemp * (VY( 0, 4,0)-2 * CD*ABCD*VY( 0, 4,1)) + ABCD * VY( 2, 1, 1);
                                VY( 8, 5, 0) = (Py-RAy) * VY( 2, 5, 0) + (tempy - Py) * VY( 2, 5, 1) +  ABtemp * (VY( 0, 5,0)-2 * CD*ABCD*VY( 0, 5,1)) + ABCD * VY( 2, 3, 1);
                                VY( 8, 6, 0) = (Py-RAy) * VY( 2, 6, 0) + (tempy - Py) * VY( 2, 6, 1) +  ABtemp * (VY( 0, 6,0)-2 * CD*ABCD*VY( 0, 6,1));
                                VY( 8, 7, 0) = (Py-RAy) * VY( 2, 7, 0) + (tempy - Py) * VY( 2, 7, 1) +  ABtemp * (VY( 0, 7,0)-2 * CD*ABCD*VY( 0, 7,1));
                                VY( 8, 8, 0) = (Py-RAy) * VY( 2, 8, 0) + (tempy - Py) * VY( 2, 8, 1) +  ABtemp * (VY( 0, 8,0)-2 * CD*ABCD*VY( 0, 8,1)) + 2*ABCD * VY( 2, 2, 1);
                                VY( 8, 9, 0) = (Py-RAy) * VY( 2, 9, 0) + (tempy - Py) * VY( 2, 9, 1) +  ABtemp * (VY( 0, 9,0)-2 * CD*ABCD*VY( 0, 9,1));
                                
                                VY( 9, 4, 0) = (Pz-RAz) * VY( 3, 4, 0) + (tempz - Pz) * VY( 3, 4, 1) +  ABtemp * (VY( 0, 4,0)-2 * CD*ABCD*VY( 0, 4,1));
                                VY( 9, 5, 0) = (Pz-RAz) * VY( 3, 5, 0) + (tempz - Pz) * VY( 3, 5, 1) +  ABtemp * (VY( 0, 5,0)-2 * CD*ABCD*VY( 0, 5,1)) + ABCD * VY( 3, 2, 1);
                                VY( 9, 6, 0) = (Pz-RAz) * VY( 3, 6, 0) + (tempz - Pz) * VY( 3, 6, 1) +  ABtemp * (VY( 0, 6,0)-2 * CD*ABCD*VY( 0, 6,1)) + ABCD * VY( 3, 1, 1);
                                VY( 9, 7, 0) = (Pz-RAz) * VY( 3, 7, 0) + (tempz - Pz) * VY( 3, 7, 1) +  ABtemp * (VY( 0, 7,0)-2 * CD*ABCD*VY( 0, 7,1));
                                VY( 9, 8, 0) = (Pz-RAz) * VY( 3, 8, 0) + (tempz - Pz) * VY( 3, 8, 1) +  ABtemp * (VY( 0, 8,0)-2 * CD*ABCD*VY( 0, 8,1));
                                VY( 9, 9, 0) = (Pz-RAz) * VY( 3, 9, 0) + (tempz - Pz) * VY( 3, 9, 1) +  ABtemp * (VY( 0, 9,0)-2 * CD*ABCD*VY( 0, 9,1)) + 2*ABCD * VY( 3, 3, 1);
                                
                            }
                        }*/
                        
                        LOC2(store, 0, 0, STOREDIM, STOREDIM) += VY( 0, 0, 0);
                        
                        //------------- FOR CASE I+J+K+L > 0 -----------------------------------
                        // includes NABCDTYPE = 1,10,2,20,11,12,21,22 situation
                        //----------------------------------------------------------------------
                        if (I+J+K+L != 0) {
                            int NABCDTYPE = (I+J)*10+K+L;
                            QUICKDouble tempx = (Px*AB+Qx*CD)*ABCD;
                            QUICKDouble tempy = (Py*AB+Qy*CD)*ABCD;
                            QUICKDouble tempz = (Pz*AB+Qz*CD)*ABCD;
                            
                            if (NABCDTYPE == 10) { //NABCDTYPE: 10
                                //PSSS(0, YVerticalTemp, Px-RAx, Py-RAy, Pz-RAz, (Px*AB+Qx*CD)*ABCD - Px, (Py*AB+Qy*CD)*ABCD - Py, tempz - Pz);
                                LOC2(store, 1, 0, STOREDIM, STOREDIM) += (Px-RAx) * VY( 0, 0, 0) + (tempx - Px) * VY( 0, 0, 1);
                                LOC2(store, 2, 0, STOREDIM, STOREDIM) += (Py-RAy) * VY( 0, 0, 0) + (tempy - Py) * VY( 0, 0, 1);
                                LOC2(store, 3, 0, STOREDIM, STOREDIM) += (Pz-RAz) * VY( 0, 0, 0) + (tempz - Pz) * VY( 0, 0, 1);
                                 
                            }
                            
                            if (NABCDTYPE == 1){ //NABCDTYPE: 10
                                //SSPS(0, YVerticalTemp, Qx-RCx, Qy-RCy, Qz-RCz, tempx - Qx, tempy - Qy, tempz - Qz);
                                LOC2(store, 0, 1, STOREDIM, STOREDIM) += (Qx-RCx) * VY( 0, 0, 0) + (tempx - Qx) * VY( 0, 0, 1);
                                LOC2(store, 0, 2, STOREDIM, STOREDIM) += (Qy-RCy) * VY( 0, 0, 0) + (tempy - Qy) * VY( 0, 0, 1);
                                LOC2(store, 0, 3, STOREDIM, STOREDIM) += (Qz-RCz) * VY( 0, 0, 0) + (tempz - Qz) * VY( 0, 0, 1);
                                 
                            }
                            
                            //------------- FOR CASE I+J+K+L > 1 -----------------------------------
                            // includes NABCDTYPE = 2,20,11,12,21,22 situation
                            //----------------------------------------------------------------------
                            if (I+J+K+L>1) { 
                                
                                if (K+L>0) {
                                    //SSPS(0, YVerticalTemp, Qx-RCx, Qy-RCy, Qz-RCz, tempx - Qx, tempy - Qy, tempz - Qz);
                                    VY( 0, 1, 0) = (Qx-RCx) * VY( 0, 0, 0) + (tempx - Qx) * VY( 0, 0, 1);
                                    VY( 0, 2, 0) = (Qy-RCy) * VY( 0, 0, 0) + (tempy - Qy) * VY( 0, 0, 1);
                                    VY( 0, 3, 0) = (Qz-RCz) * VY( 0, 0, 0) + (tempz - Qz) * VY( 0, 0, 1);
                                    if (I==0) {
                                    LOC2(store, 0, 1, STOREDIM, STOREDIM)+= VY( 0, 1, 0);
                                    LOC2(store, 0, 2, STOREDIM, STOREDIM)+= VY( 0, 2, 0);
                                    LOC2(store, 0, 3, STOREDIM, STOREDIM)+= VY( 0, 3, 0);
                                    }
                                }
                                
                                if (I+J>0) {
                                    //PSSS(0, YVerticalTemp, Px-RAx, Py-RAy, Pz-RAz, tempx - Px, tempy - Py, tempz - Pz);
                                    VY( 1, 0, 0) = (Px-RAx) * VY( 0, 0, 0) + (tempx - Px) * VY( 0, 0, 1);
                                    VY( 2, 0, 0) = (Py-RAy) * VY( 0, 0, 0) + (tempy - Py) * VY( 0, 0, 1);
                                    VY( 3, 0, 0) = (Pz-RAz) * VY( 0, 0, 0) + (tempz - Pz) * VY( 0, 0, 1);
                                    if(K==0) {
                                        LOC2(store, 1, 0, STOREDIM, STOREDIM)+= (Px-RAx) * VY( 0, 0, 0) + (tempx - Px) * VY( 0, 0, 1);
                                        LOC2(store, 2, 0, STOREDIM, STOREDIM)+= (Py-RAy) * VY( 0, 0, 0) + (tempy - Py) * VY( 0, 0, 1);
                                        LOC2(store, 3, 0, STOREDIM, STOREDIM)+= (Pz-RAz) * VY( 0, 0, 0) + (tempz - Pz) * VY( 0, 0, 1);
                                    }
                                }
                                
                                if ((I+J)>0 && (K+L)>0 || (K+L)>1) {
                                    //SSPS(1, YVerticalTemp, Qx-RCx, Qy-RCy, Qz-RCz, (Px*AB+Qx*CD)*ABCD - Qx, (Py*AB+Qy*CD)*ABCD - Qy, (Pz*AB+Qz*CD)*ABCD - Qz);
                                    VY( 0, 1, 1) = (Qx-RCx) * VY( 0, 0, 1) + (tempx - Qx) * VY( 0, 0, 2);
                                    VY( 0, 2, 1) = (Qy-RCy) * VY( 0, 0, 1) + (tempy - Qy) * VY( 0, 0, 2);
                                    VY( 0, 3, 1) = (Qz-RCz) * VY( 0, 0, 1) + (tempz - Qz) * VY( 0, 0, 2);
                                }
                                
                                if(NABCDTYPE == 2) {
                                    LOC2(store, 0, 4, STOREDIM, STOREDIM)+= (Qx-RCx) * VY( 0, 2, 0) + (tempx - Qx) * VY( 0, 2, 1);
                                    LOC2(store, 0, 5, STOREDIM, STOREDIM)+= (Qy-RCy) * VY( 0, 3, 0) + (tempy - Qy) * VY( 0, 3, 1);
                                    LOC2(store, 0, 6, STOREDIM, STOREDIM)+= (Qx-RCx) * VY( 0, 3, 0) + (tempx - Qx) * VY( 0, 3, 1);
                                    
                                    QUICKDouble tmp = 0.5 / CD *(VY( 0, 0, 0) - AB*ABCD * VY( 0, 0, 1));
                                    LOC2(store, 0, 7, STOREDIM, STOREDIM)+= (Qx-RCx) * VY( 0, 1, 0) + (tempx - Qx) * VY( 0, 1, 1) + tmp ;
                                    LOC2(store, 0, 8, STOREDIM, STOREDIM)+= (Qy-RCy) * VY( 0, 2, 0) + (tempy - Qy) * VY( 0, 2, 1) + tmp ;
                                    LOC2(store, 0, 9, STOREDIM, STOREDIM)+= (Qz-RCz) * VY( 0, 3, 0) + (tempz - Qz) * VY( 0, 3, 1) + tmp ;
                                     
                                }
                                
                                if ((I+J)>1) {
                                    //PSSS(1, YVerticalTemp, Px-RAx, Py-RAy, Pz-RAz, tempx - Px, tempy - Py, tempz - Pz);
                                    VY( 1, 0, 1) = (Px-RAx) * VY( 0, 0, 1) + (tempx - Px) * VY( 0, 0, 2);
                                    VY( 2, 0, 1) = (Py-RAy) * VY( 0, 0, 1) + (tempy - Py) * VY( 0, 0, 2);
                                    VY( 3, 0, 1) = (Pz-RAz) * VY( 0, 0, 1) + (tempz - Pz) * VY( 0, 0, 2);
                                }
                                
                                if(NABCDTYPE == 20) {
                                    LOC2(store, 4, 0, STOREDIM, STOREDIM)+= (Px-RAx) * VY( 2, 0, 0) + (tempx - Px) * VY( 2, 0, 1);
                                    LOC2(store, 5, 0, STOREDIM, STOREDIM)+= (Py-RAy) * VY( 3, 0, 0) + (tempy - Py) * VY( 3, 0, 1);
                                    LOC2(store, 6, 0, STOREDIM, STOREDIM)+= (Px-RAx) * VY( 3, 0, 0) + (tempx - Px) * VY( 3, 0, 1);
                                    
                                    QUICKDouble tmp = 0.5 / AB *(VY( 0, 0, 0) - CD * ABCD * VY( 0, 0, 1));
                                    LOC2(store, 7, 0, STOREDIM, STOREDIM)+= (Px-RAx) * VY( 1, 0, 0) + (tempx - Px) * VY( 1, 0, 1) + tmp;
                                    LOC2(store, 8, 0, STOREDIM, STOREDIM)+= (Py-RAy) * VY( 2, 0, 0) + (tempy - Py) * VY( 2, 0, 1) + tmp;
                                    LOC2(store, 9, 0, STOREDIM, STOREDIM)+= (Pz-RAz) * VY( 3, 0, 0) + (tempz - Pz) * VY( 3, 0, 1) + tmp;
                                     
                                }
                                
                                if((I+J) != 0 && (K+L) != 0) {// for case 11, 21, 12
                                
                                    QUICKDouble tmp = 0.5*ABCD * VY( 0, 0, 1);
                                    LOC2(store, 1, 1, STOREDIM, STOREDIM) += (Px-RAx) * VY( 0, 1, 0) + (tempx - Px) * VY( 0, 1, 1) + tmp;
                                    LOC2(store, 2, 1, STOREDIM, STOREDIM) += (Py-RAy) * VY( 0, 1, 0) + (tempy - Py) * VY( 0, 1, 1);
                                    LOC2(store, 3, 1, STOREDIM, STOREDIM) += (Pz-RAz) * VY( 0, 1, 0) + (tempz - Pz) * VY( 0, 1, 1);
                                    
                                    LOC2(store, 1, 2, STOREDIM, STOREDIM) += (Px-RAx) * VY( 0, 2, 0) + (tempx - Px) * VY( 0, 2, 1);
                                    LOC2(store, 2, 2, STOREDIM, STOREDIM) += (Py-RAy) * VY( 0, 2, 0) + (tempy - Py) * VY( 0, 2, 1) + tmp;
                                    LOC2(store, 3, 2, STOREDIM, STOREDIM) += (Pz-RAz) * VY( 0, 2, 0) + (tempz - Pz) * VY( 0, 2, 1);
                                    
                                    LOC2(store, 1, 3, STOREDIM, STOREDIM) += (Px-RAx) * VY( 0, 3, 0) + (tempx - Px) * VY( 0, 3, 1);
                                    LOC2(store, 2, 3, STOREDIM, STOREDIM) += (Py-RAy) * VY( 0, 3, 0) + (tempy - Py) * VY( 0, 3, 1);
                                    LOC2(store, 3, 3, STOREDIM, STOREDIM) += (Pz-RAz) * VY( 0, 3, 0) + (tempz - Pz) * VY( 0, 3, 1) + tmp;
                                }
                                
                                //------------- FOR CASE I+J+K+L > 2 -----------------------------------
                                // includes NABCDTYPE = 12,21,22 situation
                                //----------------------------------------------------------------------
                                if ((I+J+K+L)>2) {//NABCDTYPE: 21, 12, 22
                                    //PSPS(0, YVerticalTemp, Px-RAx, Py-RAy, Pz-RAz, (Px*AB+Qx*CD)*ABCD - Px, tempy - Py, tempz - Pz, 0.5*ABCD);
                                    
                                    QUICKDouble tmp = 0.5*ABCD * VY( 0, 0, 1);
                                    VY( 1, 1, 0) = (Px-RAx) * VY( 0, 1, 0) + (tempx - Px) * VY( 0, 1, 1) + tmp;
                                    VY( 2, 1, 0) = (Py-RAy) * VY( 0, 1, 0) + (tempy - Py) * VY( 0, 1, 1);
                                    VY( 3, 1, 0) = (Pz-RAz) * VY( 0, 1, 0) + (tempz - Pz) * VY( 0, 1, 1);
                                    
                                    VY( 1, 2, 0) = (Px-RAx) * VY( 0, 2, 0) + (tempx - Px) * VY( 0, 2, 1);
                                    VY( 2, 2, 0) = (Py-RAy) * VY( 0, 2, 0) + (tempy - Py) * VY( 0, 2, 1) + tmp;
                                    VY( 3, 2, 0) = (Pz-RAz) * VY( 0, 2, 0) + (tempz - Pz) * VY( 0, 2, 1);
                                    
                                    VY( 1, 3, 0) = (Px-RAx) * VY( 0, 3, 0) + (tempx - Px) * VY( 0, 3, 1);
                                    VY( 2, 3, 0) = (Py-RAy) * VY( 0, 3, 0) + (tempy - Py) * VY( 0, 3, 1);
                                    VY( 3, 3, 0) = (Pz-RAz) * VY( 0, 3, 0) + (tempz - Pz) * VY( 0, 3, 1) + tmp;
                                    
                                    if ((I+J)>1 && (K+L)>0) {
                                        //DSSS(0, YVerticalTemp, Px-RAx, Py-RAy, Pz-RAz, tempx - Px, tempy - Py, tempz - Pz, 0.5/AB, CD*ABCD);
                                        VY( 4, 0, 0) = (Px-RAx) * VY( 2, 0, 0) + (tempx - Px) * VY( 2, 0, 1);
                                        VY( 5, 0, 0) = (Py-RAy) * VY( 3, 0, 0) + (tempy - Py) * VY( 3, 0, 1);
                                        VY( 6, 0, 0) = (Px-RAx) * VY( 3, 0, 0) + (tempx - Px) * VY( 3, 0, 1);
                                        
                                        QUICKDouble tmp = 0.5 / AB *(VY( 0, 0, 0) - CD * ABCD * VY( 0, 0, 1));
                                        VY( 7, 0, 0) = (Px-RAx) * VY( 1, 0, 0) + (tempx - Px) * VY( 1, 0, 1) + tmp;
                                        VY( 8, 0, 0) = (Py-RAy) * VY( 2, 0, 0) + (tempy - Py) * VY( 2, 0, 1) + tmp;
                                        VY( 9, 0, 0) = (Pz-RAz) * VY( 3, 0, 0) + (tempz - Pz) * VY( 3, 0, 1) + tmp;
                                        
                                        
                                        //PSSS(2, YVerticalTemp, Px-RAx, Py-RAy, Pz-RAz, tempx - Px, tempy - Py, tempz - Pz);
                                        VY( 1, 0, 2) = (Px-RAx) * VY( 0, 0, 2) + (tempx - Px) * VY( 0, 0, 3);
                                        VY( 2, 0, 2) = (Py-RAy) * VY( 0, 0, 2) + (tempy - Py) * VY( 0, 0, 3);
                                        VY( 3, 0, 2) = (Pz-RAz) * VY( 0, 0, 2) + (tempz - Pz) * VY( 0, 0, 3);
                                        
                                        //DSSS(1, YVerticalTemp, Px-RAx, Py-RAy, Pz-RAz, tempx - Px, tempy - Py, tempz - Pz, 0.5/AB, CD*ABCD);
                                        VY( 4, 0, 1) = (Px-RAx) * VY( 2, 0, 1) + (tempx - Px) * VY( 2, 0, 2);
                                        VY( 5, 0, 1) = (Py-RAy) * VY( 3, 0, 1) + (tempy - Py) * VY( 3, 0, 2);
                                        VY( 6, 0, 1) = (Px-RAx) * VY( 3, 0, 1) + (tempx - Px) * VY( 3, 0, 2);
                                        
                                        tmp = 0.5 / AB *(VY( 0, 0, 1) - CD * ABCD * VY( 0, 0, 2));
                                        VY( 7, 0, 1) = (Px-RAx) * VY( 1, 0, 1) + (tempx - Px) * VY( 1, 0, 2) + tmp;
                                        VY( 8, 0, 1) = (Py-RAy) * VY( 2, 0, 1) + (tempy - Py) * VY( 2, 0, 2) + tmp;
                                        VY( 9, 0, 1) = (Pz-RAz) * VY( 3, 0, 1) + (tempz - Pz) * VY( 3, 0, 2) + tmp;
                                        
                                        LOC2(store, 4, 0, STOREDIM, STOREDIM) += VY( 4, 0, 0);
                                        LOC2(store, 5, 0, STOREDIM, STOREDIM) += VY( 5, 0, 0);
                                        LOC2(store, 6, 0, STOREDIM, STOREDIM) += VY( 6, 0, 0);
                                        LOC2(store, 7, 0, STOREDIM, STOREDIM) += VY( 7, 0, 0);
                                        LOC2(store, 8, 0, STOREDIM, STOREDIM) += VY( 8, 0, 0);
                                        LOC2(store, 9, 0, STOREDIM, STOREDIM) += VY( 9, 0, 0);
                                        
                                        tmp = 0.5 * ABCD * VY( 2, 0, 1);
                                        LOC2(store, 4, 1, STOREDIM, STOREDIM) += (Qx-RCx) * VY( 4, 0, 0) + (tempx - Qx) * VY( 4, 0, 1) + tmp;
                                        LOC2(store, 5, 3, STOREDIM, STOREDIM) += (Qz-RCz) * VY( 5, 0, 0) + (tempz - Qz) * VY( 5, 0, 1) + tmp;
                                        LOC2(store, 8, 2, STOREDIM, STOREDIM) += (Qy-RCy) * VY( 8, 0, 0) + (tempy - Qy) * VY( 8, 0, 1) + 2*tmp;
                                        
                                        tmp = 0.5 * ABCD * VY( 1, 0, 1);
                                        LOC2(store, 7, 1, STOREDIM, STOREDIM) += (Qx-RCx) * VY( 7, 0, 0) + (tempx - Qx) * VY( 7, 0, 1) + 2*tmp;
                                        LOC2(store, 4, 2, STOREDIM, STOREDIM) += (Qy-RCy) * VY( 4, 0, 0) + (tempy - Qy) * VY( 4, 0, 1) + tmp;
                                        LOC2(store, 6, 3, STOREDIM, STOREDIM) += (Qz-RCz) * VY( 6, 0, 0) + (tempz - Qz) * VY( 6, 0, 1) + tmp;
                                        
                                        tmp = 0.5 * ABCD * VY( 3, 0, 1);
                                        LOC2(store, 6, 1, STOREDIM, STOREDIM) += (Qx-RCx) * VY( 6, 0, 0) + (tempx - Qx) * VY( 6, 0, 1) + tmp;
                                        LOC2(store, 5, 2, STOREDIM, STOREDIM) += (Qy-RCy) * VY( 5, 0, 0) + (tempy - Qy) * VY( 5, 0, 1) + tmp;
                                        LOC2(store, 9, 3, STOREDIM, STOREDIM) += (Qz-RCz) * VY( 9, 0, 0) + (tempz - Qz) * VY( 9, 0, 1) + 2*tmp;
                                        
                                        LOC2(store, 5, 1, STOREDIM, STOREDIM) += (Qx-RCx) * VY( 5, 0, 0) + (tempx - Qx) * VY( 5, 0, 1);
                                        LOC2(store, 8, 1, STOREDIM, STOREDIM) += (Qx-RCx) * VY( 8, 0, 0) + (tempx - Qx) * VY( 8, 0, 1);
                                        LOC2(store, 9, 1, STOREDIM, STOREDIM) += (Qx-RCx) * VY( 9, 0, 0) + (tempx - Qx) * VY( 9, 0, 1);
                                        
                                        LOC2(store, 6, 2, STOREDIM, STOREDIM) += (Qy-RCy) * VY( 6, 0, 0) + (tempy - Qy) * VY( 6, 0, 1);
                                        LOC2(store, 7, 2, STOREDIM, STOREDIM) += (Qy-RCy) * VY( 7, 0, 0) + (tempy - Qy) * VY( 7, 0, 1);
                                        LOC2(store, 9, 2, STOREDIM, STOREDIM) += (Qy-RCy) * VY( 9, 0, 0) + (tempy - Qy) * VY( 9, 0, 1);
                                        
                                        LOC2(store, 4, 3, STOREDIM, STOREDIM) += (Qz-RCz) * VY( 4, 0, 0) + (tempz - Qz) * VY( 4, 0, 1);
                                        LOC2(store, 7, 3, STOREDIM, STOREDIM) += (Qz-RCz) * VY( 7, 0, 0) + (tempz - Qz) * VY( 7, 0, 1);
                                        LOC2(store, 8, 3, STOREDIM, STOREDIM) += (Qz-RCz) * VY( 8, 0, 0) + (tempz - Qz) * VY( 8, 0, 1);
                                           
                                    }
                                    
                                    if ((K+L)>1 && (I+J)>0) {
                                    
                                        //SSDS(0, YVerticalTemp, Qx-RCx, Qy-RCy, Qz-RCz, tempx - Qx, tempy - Qy, tempz - Qz, 0.5/CD, AB*ABCD);
                                        VY( 0, 4, 0) = (Qx-RCx) * VY( 0, 2, 0) + (tempx - Qx) * VY( 0, 2, 1);
                                        VY( 0, 5, 0) = (Qy-RCy) * VY( 0, 3, 0) + (tempy - Qy) * VY( 0, 3, 1);
                                        VY( 0, 6, 0) = (Qx-RCx) * VY( 0, 3, 0) + (tempx - Qx) * VY( 0, 3, 1);
                                        
                                        QUICKDouble tmp =  0.5 / CD *(VY( 0, 0, 0) - AB*ABCD * VY( 0, 0, 1));
                                        VY( 0, 7, 0) = (Qx-RCx) * VY( 0, 1, 0) + (tempx - Qx) * VY( 0, 1, 1) + tmp;
                                        VY( 0, 8, 0) = (Qy-RCy) * VY( 0, 2, 0) + (tempy - Qy) * VY( 0, 2, 1) + tmp;
                                        VY( 0, 9, 0) = (Qz-RCz) * VY( 0, 3, 0) + (tempz - Qz) * VY( 0, 3, 1) + tmp;                                
                                        
                                        //SSPS(2, YVerticalTemp, Qx-RCx, Qy-RCy, Qz-RCz, tempx - Qx, tempy - Qy, tempz - Qz);
                                        VY( 0, 1, 2) = (Qx-RCx) * VY( 0, 0, 2) + (tempx - Qx) * VY( 0, 0, 3);
                                        VY( 0, 2, 2) = (Qy-RCy) * VY( 0, 0, 2) + (tempy - Qy) * VY( 0, 0, 3);
                                        VY( 0, 3, 2) = (Qz-RCz) * VY( 0, 0, 2) + (tempz - Qz) * VY( 0, 0, 3);
                                        
                                        //SSDS(1, YVerticalTemp, Qx-RCx, Qy-RCy, Qz-RCz, tempx - Qx, tempy - Qy, tempz - Qz, 0.5/CD, AB*ABCD);
                                        VY( 0, 4, 1) = (Qx-RCx) * VY( 0, 2, 1) + (tempx - Qx) * VY( 0, 2, 2);
                                        VY( 0, 5, 1) = (Qy-RCy) * VY( 0, 3, 1) + (tempy - Qy) * VY( 0, 3, 2);
                                        VY( 0, 6, 1) = (Qx-RCx) * VY( 0, 3, 1) + (tempx - Qx) * VY( 0, 3, 2);
                                        
                                        tmp = 0.5 / CD *(VY( 0, 0, 1) - AB*ABCD * VY( 0, 0, 2)) ;
                                        VY( 0, 7, 1) = (Qx-RCx) * VY( 0, 1, 1) + (tempx - Qx) * VY( 0, 1, 2) + tmp;
                                        VY( 0, 8, 1) = (Qy-RCy) * VY( 0, 2, 1) + (tempy - Qy) * VY( 0, 2, 2) + tmp;
                                        VY( 0, 9, 1) = (Qz-RCz) * VY( 0, 3, 1) + (tempz - Qz) * VY( 0, 3, 2) + tmp;
                                        
                                        LOC2(store, 0, 4, STOREDIM, STOREDIM) += VY( 0, 4, 0);
                                        LOC2(store, 0, 5, STOREDIM, STOREDIM) += VY( 0, 5, 0);
                                        LOC2(store, 0, 6, STOREDIM, STOREDIM) += VY( 0, 6, 0);
                                        LOC2(store, 0, 7, STOREDIM, STOREDIM) += VY( 0, 7, 0);
                                        LOC2(store, 0, 8, STOREDIM, STOREDIM) += VY( 0, 8, 0);
                                        LOC2(store, 0, 9, STOREDIM, STOREDIM) += VY( 0, 9, 0);
                                        
                                        tmp = 0.5 * ABCD * VY( 0, 1, 1);
                                        LOC2(store, 2, 4, STOREDIM, STOREDIM) += (Py-RAy) * VY( 0, 4, 0) + (tempy - Py) * VY( 0, 4, 1) + tmp;
                                        LOC2(store, 3, 6, STOREDIM, STOREDIM) += (Pz-RAz) * VY( 0, 6, 0) + (tempz - Pz) * VY( 0, 6, 1) + tmp;
                                        LOC2(store, 1, 7, STOREDIM, STOREDIM) += (Px-RAx) * VY( 0, 7, 0) + (tempx - Px) * VY( 0, 7, 1) + 2*tmp;
                                        
                                        tmp = 0.5 * ABCD * VY( 0, 2, 1);
                                        LOC2(store, 1, 4, STOREDIM, STOREDIM) += (Px-RAx) * VY( 0, 4, 0) + (tempx - Px) * VY( 0, 4, 1) + tmp;
                                        LOC2(store, 2, 8, STOREDIM, STOREDIM) += (Py-RAy) * VY( 0, 8, 0) + (tempy - Py) * VY( 0, 8, 1) + 2*tmp;
                                        LOC2(store, 3, 5, STOREDIM, STOREDIM) += (Pz-RAz) * VY( 0, 5, 0) + (tempz - Pz) * VY( 0, 5, 1) + tmp;
                                        
                                        tmp = 0.5 * ABCD * VY( 0, 3, 1);
                                        LOC2(store, 1, 6, STOREDIM, STOREDIM) += (Px-RAx) * VY( 0, 6, 0) + (tempx - Px) * VY( 0, 6, 1) + tmp;
                                        LOC2(store, 2, 5, STOREDIM, STOREDIM) += (Py-RAy) * VY( 0, 5, 0) + (tempy - Py) * VY( 0, 5, 1) + tmp;
                                        LOC2(store, 3, 9, STOREDIM, STOREDIM) += (Pz-RAz) * VY( 0, 9, 0) + (tempz - Pz) * VY( 0, 9, 1) + 2*tmp;
                                    
                                        LOC2(store, 1, 5, STOREDIM, STOREDIM) += (Px-RAx) * VY( 0, 5, 0) + (tempx - Px) * VY( 0, 5, 1);
                                        LOC2(store, 1, 8, STOREDIM, STOREDIM) += (Px-RAx) * VY( 0, 8, 0) + (tempx - Px) * VY( 0, 8, 1);
                                        LOC2(store, 1, 9, STOREDIM, STOREDIM) += (Px-RAx) * VY( 0, 9, 0) + (tempx - Px) * VY( 0, 9, 1);
                                        
                                        LOC2(store, 2, 6, STOREDIM, STOREDIM) += (Py-RAy) * VY( 0, 6, 0) + (tempy - Py) * VY( 0, 6, 1);
                                        LOC2(store, 2, 7, STOREDIM, STOREDIM) += (Py-RAy) * VY( 0, 7, 0) + (tempy - Py) * VY( 0, 7, 1);
                                        LOC2(store, 2, 9, STOREDIM, STOREDIM) += (Py-RAy) * VY( 0, 9, 0) + (tempy - Py) * VY( 0, 9, 1);
                                        
                                        LOC2(store, 3, 4, STOREDIM, STOREDIM) += (Pz-RAz) * VY( 0, 4, 0) + (tempz - Pz) * VY( 0, 4, 1);
                                        LOC2(store, 3, 7, STOREDIM, STOREDIM) += (Pz-RAz) * VY( 0, 7, 0) + (tempz - Pz) * VY( 0, 7, 1);
                                        LOC2(store, 3, 8, STOREDIM, STOREDIM) += (Pz-RAz) * VY( 0, 8, 0) + (tempz - Pz) * VY( 0, 8, 1);
                                        
                                    }
                                    
                                    if(NABCDTYPE == 22) {
                                        
                                        //PSDS(0, YVerticalTemp, Px-RAx, Py-RAy, Pz-RAz, tempx - Px, tempy - Py, tempz - Pz, 0.5*ABCD);
                                        QUICKDouble tmp = 0.5 * ABCD * VY( 0, 1, 1);
                                        VY( 2, 4, 0) = (Py-RAy) * VY( 0, 4, 0) + (tempy - Py) * VY( 0, 4, 1) + tmp;
                                        VY( 3, 6, 0) = (Pz-RAz) * VY( 0, 6, 0) + (tempz - Pz) * VY( 0, 6, 1) + tmp;
                                        VY( 1, 7, 0) = (Px-RAx) * VY( 0, 7, 0) + (tempx - Px) * VY( 0, 7, 1) + 2*tmp;
                                        tmp = 0.5 * ABCD * VY( 0, 2, 1);
                                        VY( 1, 4, 0) = (Px-RAx) * VY( 0, 4, 0) + (tempx - Px) * VY( 0, 4, 1) + tmp;
                                        VY( 3, 5, 0) = (Pz-RAz) * VY( 0, 5, 0) + (tempz - Pz) * VY( 0, 5, 1) + tmp;
                                        VY( 2, 8, 0) = (Py-RAy) * VY( 0, 8, 0) + (tempy - Py) * VY( 0, 8, 1) + 2*tmp;
                                        tmp = 0.5 * ABCD * VY( 0, 3, 1);
                                        VY( 2, 5, 0) = (Py-RAy) * VY( 0, 5, 0) + (tempy - Py) * VY( 0, 5, 1) + tmp;
                                        VY( 1, 6, 0) = (Px-RAx) * VY( 0, 6, 0) + (tempx - Px) * VY( 0, 6, 1) + tmp;
                                        VY( 3, 9, 0) = (Pz-RAz) * VY( 0, 9, 0) + (tempz - Pz) * VY( 0, 9, 1) + 2*tmp;
                                        
                                        VY( 3, 4, 0) = (Pz-RAz) * VY( 0, 4, 0) + (tempz - Pz) * VY( 0, 4, 1) ;
                                        VY( 1, 5, 0) = (Px-RAx) * VY( 0, 5, 0) + (tempx - Px) * VY( 0, 5, 1);
                                        VY( 2, 6, 0) = (Py-RAy) * VY( 0, 6, 0) + (tempy - Py) * VY( 0, 6, 1);
                                        VY( 2, 7, 0) = (Py-RAy) * VY( 0, 7, 0) + (tempy - Py) * VY( 0, 7, 1);
                                        VY( 3, 7, 0) = (Pz-RAz) * VY( 0, 7, 0) + (tempz - Pz) * VY( 0, 7, 1);
                                        VY( 1, 8, 0) = (Px-RAx) * VY( 0, 8, 0) + (tempx - Px) * VY( 0, 8, 1);
                                        VY( 3, 8, 0) = (Pz-RAz) * VY( 0, 8, 0) + (tempz - Pz) * VY( 0, 8, 1);
                                        VY( 1, 9, 0) = (Px-RAx) * VY( 0, 9, 0) + (tempx - Px) * VY( 0, 9, 1);
                                        VY( 2, 9, 0) = (Py-RAy) * VY( 0, 9, 0) + (tempy - Py) * VY( 0, 9, 1);
                                        
                                        //PSSS(1, YVerticalTemp, Px-RAx, Py-RAy, Pz-RAz, tempx - Px, tempy - Py, tempz - Pz);
                                        VY( 1, 0, 1) = (Px-RAx) * VY( 0, 0, 1) + (tempx - Px) * VY( 0, 0, 2);
                                        VY( 2, 0, 1) = (Py-RAy) * VY( 0, 0, 1) + (tempy - Py) * VY( 0, 0, 2);
                                        VY( 3, 0, 1) = (Pz-RAz) * VY( 0, 0, 1) + (tempz - Pz) * VY( 0, 0, 2);
                                        
                                        //SSPS(3, YVerticalTemp, Qx-RCx, Qy-RCy, Qz-RCz, tempx - Qx, tempy - Qy, tempz - Qz);
                                        VY( 0, 1, 3) = (Qx-RCx) * VY( 0, 0, 3) + (tempx - Qx) * VY( 0, 0, 4);
                                        VY( 0, 2, 3) = (Qy-RCy) * VY( 0, 0, 3) + (tempy - Qy) * VY( 0, 0, 4);
                                        VY( 0, 3, 3) = (Qz-RCz) * VY( 0, 0, 3) + (tempz - Qz) * VY( 0, 0, 4);
                                        
                                        //SSDS(2, YVerticalTemp, Qx-RCx, Qy-RCy, Qz-RCz, tempx - Qx, tempy - Qy, tempz - Qz, 0.5/CD, AB*ABCD);
                                        VY( 0, 4, 2) = (Qx-RCx) * VY( 0, 2, 2) + (tempx - Qx) * VY( 0, 2, 3);
                                        VY( 0, 5, 2) = (Qy-RCy) * VY( 0, 3, 2) + (tempy - Qy) * VY( 0, 3, 3);
                                        VY( 0, 6, 2) = (Qx-RCx) * VY( 0, 3, 2) + (tempx - Qx) * VY( 0, 3, 3);
                                        
                                        tmp = 0.5 / CD *(VY( 0, 0, 2) - AB*ABCD * VY( 0, 0, 3)) ;
                                        VY( 0, 7, 2) = (Qx-RCx) * VY( 0, 1, 2) + (tempx - Qx) * VY( 0, 1, 3) + tmp;
                                        VY( 0, 8, 2) = (Qy-RCy) * VY( 0, 2, 2) + (tempy - Qy) * VY( 0, 2, 3) + tmp;
                                        VY( 0, 9, 2) = (Qz-RCz) * VY( 0, 3, 2) + (tempz - Qz) * VY( 0, 3, 3) + tmp;
                                        
                                        //PSDS(1, YVerticalTemp, Px-RAx, Py-RAy, Pz-RAz, tempx - Px, tempy - Py, tempz - Pz, 0.5*ABCD);
                                        tmp = 0.5 * ABCD * VY( 0, 1, 2);
                                        VY( 2, 4, 1) = (Py-RAy) * VY( 0, 4, 1) + (tempy - Py) * VY( 0, 4, 2) + tmp;
                                        VY( 3, 6, 1) = (Pz-RAz) * VY( 0, 6, 1) + (tempz - Pz) * VY( 0, 6, 2) + tmp;
                                        VY( 1, 7, 1) = (Px-RAx) * VY( 0, 7, 1) + (tempx - Px) * VY( 0, 7, 2) + 2 * tmp;
                                        
                                        tmp = 0.5 * ABCD * VY( 0, 2, 2);
                                        VY( 2, 8, 1) = (Py-RAy) * VY( 0, 8, 1) + (tempy - Py) * VY( 0, 8, 2) + 2 * tmp;
                                        VY( 1, 4, 1) = (Px-RAx) * VY( 0, 4, 1) + (tempx - Px) * VY( 0, 4, 2) + tmp;
                                        VY( 3, 5, 1) = (Pz-RAz) * VY( 0, 5, 1) + (tempz - Pz) * VY( 0, 5, 2) + tmp;
                                        
                                        tmp = 0.5 * ABCD * VY( 0, 3, 2);
                                        VY( 2, 5, 1) = (Py-RAy) * VY( 0, 5, 1) + (tempy - Py) * VY( 0, 5, 2) + tmp;
                                        VY( 1, 6, 1) = (Px-RAx) * VY( 0, 6, 1) + (tempx - Px) * VY( 0, 6, 2) + tmp;
                                        VY( 3, 9, 1) = (Pz-RAz) * VY( 0, 9, 1) + (tempz - Pz) * VY( 0, 9, 2) + 2 * tmp;
                                        VY( 3, 4, 1) = (Pz-RAz) * VY( 0, 4, 1) + (tempz - Pz) * VY( 0, 4, 2) ;
                                        
                                        VY( 1, 5, 1) = (Px-RAx) * VY( 0, 5, 1) + (tempx - Px) * VY( 0, 5, 2);
                                        VY( 2, 6, 1) = (Py-RAy) * VY( 0, 6, 1) + (tempy - Py) * VY( 0, 6, 2);
                                        VY( 2, 7, 1) = (Py-RAy) * VY( 0, 7, 1) + (tempy - Py) * VY( 0, 7, 2);
                                        VY( 3, 7, 1) = (Pz-RAz) * VY( 0, 7, 1) + (tempz - Pz) * VY( 0, 7, 2);
                                        VY( 1, 8, 1) = (Px-RAx) * VY( 0, 8, 1) + (tempx - Px) * VY( 0, 8, 2);
                                        VY( 3, 8, 1) = (Pz-RAz) * VY( 0, 8, 1) + (tempz - Pz) * VY( 0, 8, 2);
                                        
                                        VY( 1, 9, 1) = (Px-RAx) * VY( 0, 9, 1) + (tempx - Px) * VY( 0, 9, 2);
                                        VY( 2, 9, 1) = (Py-RAy) * VY( 0, 9, 1) + (tempy - Py) * VY( 0, 9, 2);
                                        
                                        //PSPS(1, YVerticalTemp, Px-RAx, Py-RAy, Pz-RAz, tempx - Px, tempy - Py, tempz - Pz, 0.5*ABCD);
                                        VY( 1, 1, 1) = (Px-RAx) * VY( 0, 1, 1) + (tempx - Px) * VY( 0, 1, 2) + 0.5*ABCD * VY( 0, 0, 2);
                                        VY( 2, 1, 1) = (Py-RAy) * VY( 0, 1, 1) + (tempy - Py) * VY( 0, 1, 2);
                                        VY( 3, 1, 1) = (Pz-RAz) * VY( 0, 1, 1) + (tempz - Pz) * VY( 0, 1, 2);
                                        
                                        VY( 1, 2, 1) = (Px-RAx) * VY( 0, 2, 1) + (tempx - Px) * VY( 0, 2, 2);
                                        VY( 2, 2, 1) = (Py-RAy) * VY( 0, 2, 1) + (tempy - Py) * VY( 0, 2, 2) + 0.5*ABCD * VY( 0, 0, 2);
                                        VY( 3, 2, 1) = (Pz-RAz) * VY( 0, 2, 1) + (tempz - Pz) * VY( 0, 2, 2);
                                        
                                        VY( 1, 3, 1) = (Px-RAx) * VY( 0, 3, 1) + (tempx - Px) * VY( 0, 3, 2);
                                        VY( 2, 3, 1) = (Py-RAy) * VY( 0, 3, 1) + (tempy - Py) * VY( 0, 3, 2);
                                        VY( 3, 3, 1) = (Pz-RAz) * VY( 0, 3, 1) + (tempz - Pz) * VY( 0, 3, 2) + 0.5*ABCD * VY( 0, 0, 2);
                                        
                                        LOC2(store, 4, 4, STOREDIM, STOREDIM) += (Px-RAx) * VY( 2, 4, 0) + (tempx - Px) * VY( 2, 4, 1) + 0.5 * ABCD * VY( 2, 2, 1);
                                        LOC2(store, 4, 5, STOREDIM, STOREDIM) += (Px-RAx) * VY( 2, 5, 0) + (tempx - Px) * VY( 2, 5, 1);
                                        LOC2(store, 4, 6, STOREDIM, STOREDIM) += (Px-RAx) * VY( 2, 6, 0) + (tempx - Px) * VY( 2, 6, 1) + 0.5 * ABCD * VY( 2, 3, 1);
                                        LOC2(store, 4, 7, STOREDIM, STOREDIM) += (Px-RAx) * VY( 2, 7, 0) + (tempx - Px) * VY( 2, 7, 1) + ABCD * VY( 2, 1, 1);
                                        LOC2(store, 4, 8, STOREDIM, STOREDIM) += (Px-RAx) * VY( 2, 8, 0) + (tempx - Px) * VY( 2, 8, 1);
                                        LOC2(store, 4, 9, STOREDIM, STOREDIM) += (Px-RAx) * VY( 2, 9, 0) + (tempx - Px) * VY( 2, 9, 1);
                                        
                                        LOC2(store, 5, 4, STOREDIM, STOREDIM) += (Py-RAy) * VY( 3, 4, 0) + (tempy - Py) * VY( 3, 4, 1) + 0.5 * ABCD * VY( 3, 1, 1);
                                        LOC2(store, 5, 5, STOREDIM, STOREDIM) += (Py-RAy) * VY( 3, 5, 0) + (tempy - Py) * VY( 3, 5, 1) + 0.5 * ABCD * VY( 3, 3, 1);
                                        LOC2(store, 5, 6, STOREDIM, STOREDIM) += (Py-RAy) * VY( 3, 6, 0) + (tempy - Py) * VY( 3, 6, 1);
                                        LOC2(store, 5, 7, STOREDIM, STOREDIM) += (Py-RAy) * VY( 3, 7, 0) + (tempy - Py) * VY( 3, 7, 1);
                                        LOC2(store, 5, 8, STOREDIM, STOREDIM) += (Py-RAy) * VY( 3, 8, 0) + (tempy - Py) * VY( 3, 8, 1) + ABCD * VY( 3, 2, 1);
                                        LOC2(store, 5, 9, STOREDIM, STOREDIM) += (Py-RAy) * VY( 3, 9, 0) + (tempy - Py) * VY( 3, 9, 1);
                                        
                                        LOC2(store, 6, 4, STOREDIM, STOREDIM) += (Px-RAx) * VY( 3, 4, 0) + (tempx - Px) * VY( 3, 4, 1) + 0.5 * ABCD * VY( 3, 2, 1);
                                        LOC2(store, 6, 5, STOREDIM, STOREDIM) += (Px-RAx) * VY( 3, 5, 0) + (tempx - Px) * VY( 3, 5, 1);
                                        LOC2(store, 6, 6, STOREDIM, STOREDIM) += (Px-RAx) * VY( 3, 6, 0) + (tempx - Px) * VY( 3, 6, 1) + 0.5 * ABCD * VY( 3, 3, 1);
                                        LOC2(store, 6, 7, STOREDIM, STOREDIM) += (Px-RAx) * VY( 3, 7, 0) + (tempx - Px) * VY( 3, 7, 1) + ABCD * VY( 3, 1, 1);
                                        LOC2(store, 6, 8, STOREDIM, STOREDIM) += (Px-RAx) * VY( 3, 8, 0) + (tempx - Px) * VY( 3, 8, 1);
                                        LOC2(store, 6, 9, STOREDIM, STOREDIM) += (Px-RAx) * VY( 3, 9, 0) + (tempx - Px) * VY( 3, 9, 1);
                                        
                                        tmp = 0.5/AB;
                                        VY( 0, 4,1) = tmp * (VY( 0, 4,0)-CD * ABCD*VY( 0, 4, 1));
                                        VY( 0, 5,1) = tmp * (VY( 0, 5,0)-CD * ABCD*VY( 0, 5, 1));
                                        VY( 0, 6,1) = tmp * (VY( 0, 6,0)-CD * ABCD*VY( 0, 6, 1));
                                        VY( 0, 7,1) = tmp * (VY( 0, 7,0)-CD * ABCD*VY( 0, 7, 1));
                                        VY( 0, 8,1) = tmp * (VY( 0, 8,0)-CD * ABCD*VY( 0, 8, 1));
                                        VY( 0, 9,1) = tmp * (VY( 0, 9,0)-CD * ABCD*VY( 0, 9, 1));
                                        
                                        LOC2(store, 7, 4, STOREDIM, STOREDIM) += (Px-RAx) * VY( 1, 4, 0) + (tempx - Px) * VY( 1, 4, 1) +  VY( 0, 4, 1) + 0.5 * ABCD * VY( 1, 2, 1);
                                        LOC2(store, 7, 5, STOREDIM, STOREDIM) += (Px-RAx) * VY( 1, 5, 0) + (tempx - Px) * VY( 1, 5, 1) +  VY( 0, 5, 1) ;
                                        LOC2(store, 7, 6, STOREDIM, STOREDIM) += (Px-RAx) * VY( 1, 6, 0) + (tempx - Px) * VY( 1, 6, 1) +  VY( 0, 6, 1) + 0.5 * ABCD * VY( 1, 3, 1);
                                        LOC2(store, 7, 7, STOREDIM, STOREDIM) += (Px-RAx) * VY( 1, 7, 0) + (tempx - Px) * VY( 1, 7, 1) +  VY( 0, 7, 1) + ABCD * VY( 1, 1, 1);
                                        LOC2(store, 7, 8, STOREDIM, STOREDIM) += (Px-RAx) * VY( 1, 8, 0) + (tempx - Px) * VY( 1, 8, 1) +  VY( 0, 8, 1) ;
                                        LOC2(store, 7, 9, STOREDIM, STOREDIM) += (Px-RAx) * VY( 1, 9, 0) + (tempx - Px) * VY( 1, 9, 1) +  VY( 0, 9, 1) ;
                                        
                                        LOC2(store, 8, 4, STOREDIM, STOREDIM) += (Py-RAy) * VY( 2, 4, 0) + (tempy - Py) * VY( 2, 4, 1) +  VY( 0, 4, 1) + 0.5 * ABCD * VY( 2, 1, 1);
                                        LOC2(store, 8, 5, STOREDIM, STOREDIM) += (Py-RAy) * VY( 2, 5, 0) + (tempy - Py) * VY( 2, 5, 1) +  VY( 0, 5, 1)  + 0.5 * ABCD * VY( 2, 3, 1);
                                        LOC2(store, 8, 6, STOREDIM, STOREDIM) += (Py-RAy) * VY( 2, 6, 0) + (tempy - Py) * VY( 2, 6, 1) +  VY( 0, 6, 1) ;
                                        LOC2(store, 8, 7, STOREDIM, STOREDIM) += (Py-RAy) * VY( 2, 7, 0) + (tempy - Py) * VY( 2, 7, 1) +  VY( 0, 7, 1) ;
                                        LOC2(store, 8, 8, STOREDIM, STOREDIM) += (Py-RAy) * VY( 2, 8, 0) + (tempy - Py) * VY( 2, 8, 1) +  VY( 0, 8, 1) + ABCD * VY( 2, 2, 1);
                                        LOC2(store, 8, 9, STOREDIM, STOREDIM) += (Py-RAy) * VY( 2, 9, 0) + (tempy - Py) * VY( 2, 9, 1) +  VY( 0, 9, 1) ;
                                        
                                        
                                        LOC2(store, 9, 4, STOREDIM, STOREDIM) += (Pz-RAz) * VY( 3, 4, 0) + (tempz - Pz) * VY( 3, 4, 1) +  VY( 0, 4, 1);
                                        LOC2(store, 9, 5, STOREDIM, STOREDIM) += (Pz-RAz) * VY( 3, 5, 0) + (tempz - Pz) * VY( 3, 5, 1) +  VY( 0, 5, 1) + 0.5 * ABCD * VY( 3, 2, 1);
                                        LOC2(store, 9, 6, STOREDIM, STOREDIM) += (Pz-RAz) * VY( 3, 6, 0) + (tempz - Pz) * VY( 3, 6, 1) +  VY( 0, 6, 1) + 0.5 * ABCD * VY( 3, 1, 1);
                                        LOC2(store, 9, 7, STOREDIM, STOREDIM) += (Pz-RAz) * VY( 3, 7, 0) + (tempz - Pz) * VY( 3, 7, 1) +  VY( 0, 7, 1) ;
                                        LOC2(store, 9, 8, STOREDIM, STOREDIM) += (Pz-RAz) * VY( 3, 8, 0) + (tempz - Pz) * VY( 3, 8, 1) +  VY( 0, 8, 1) ;
                                        LOC2(store, 9, 9, STOREDIM, STOREDIM) += (Pz-RAz) * VY( 3, 9, 0) + (tempz - Pz) * VY( 3, 9, 1) +  VY( 0, 9, 1) + ABCD * VY( 3, 3, 1);
                                    }
                                    
                                }
                            }
                        }
                    }   
                }
            }
        }
    }

    
    // IJKLTYPE is the I, J, K,L type
    int IJKLTYPE = (int) (1000 * I + 100 *J + 10 * K + L);
    
    RBx = LOC2(devSim.xyz, 0 , devSim.katom[JJ]-1, 3, devSim.natom);
    RBy = LOC2(devSim.xyz, 1 , devSim.katom[JJ]-1, 3, devSim.natom);
    RBz = LOC2(devSim.xyz, 2 , devSim.katom[JJ]-1, 3, devSim.natom);
    
    
    RDx = LOC2(devSim.xyz, 0 , devSim.katom[LL]-1, 3, devSim.natom);
    RDy = LOC2(devSim.xyz, 1 , devSim.katom[LL]-1, 3, devSim.natom);
    RDz = LOC2(devSim.xyz, 2 , devSim.katom[LL]-1, 3, devSim.natom);
        
    int III1 = LOC2(devSim.Qsbasis, II, I, nshell, 4);
    int III2 = LOC2(devSim.Qfbasis, II, I, nshell, 4);
    int JJJ1 = LOC2(devSim.Qsbasis, JJ, J, nshell, 4);
    int JJJ2 = LOC2(devSim.Qfbasis, JJ, J, nshell, 4);
    int KKK1 = LOC2(devSim.Qsbasis, KK, K, nshell, 4);
    int KKK2 = LOC2(devSim.Qfbasis, KK, K, nshell, 4);
    int LLL1 = LOC2(devSim.Qsbasis, LL, L, nshell, 4);
    int LLL2 = LOC2(devSim.Qfbasis, LL, L, nshell, 4);
    
    
    // maxIJKL is the max of I,J,K,L
    int maxIJKL = (int)MAX(MAX(I,J),MAX(K,L));
    
    if (((maxIJKL == 2)&&(J != 0 || L!=0)) || (maxIJKL >= 3)) {
        IJKLTYPE = 999;
    }
    
    if ((II < JJ) && (II < KK) && (KK < LL)) {
        for (int III = III1; III <= III2; III++) {
            for (int JJJ = JJJ1; JJJ <= JJJ2; JJJ++) {
                for (int KKK = KKK1; KKK <= KKK2 ; KKK++) {
                    for (int LLL = LLL1; LLL <= LLL2; LLL++) {
                        QUICKDouble Y = (QUICKDouble) hrrwhole(III, JJJ, KKK, LLL, IJKLTYPE, store, \
                                                               RAx, RAy, RAz, RBx, RBy, RBz, \
                                                               RCx, RCy, RCz, RDx, RDy, RDz);
                        QUICKDouble DENSEKI = (QUICKDouble) LOC2(devSim.dense, KKK-1, III-1, devSim.nbasis, devSim.nbasis);
                        QUICKDouble DENSEKJ = (QUICKDouble) LOC2(devSim.dense, KKK-1, JJJ-1, devSim.nbasis, devSim.nbasis);
                        QUICKDouble DENSELJ = (QUICKDouble) LOC2(devSim.dense, LLL-1, JJJ-1, devSim.nbasis, devSim.nbasis);
                        QUICKDouble DENSELI = (QUICKDouble) LOC2(devSim.dense, LLL-1, III-1, devSim.nbasis, devSim.nbasis);
                        QUICKDouble DENSELK = (QUICKDouble) LOC2(devSim.dense, LLL-1, KKK-1, devSim.nbasis, devSim.nbasis);
                        QUICKDouble DENSEJI = (QUICKDouble) LOC2(devSim.dense, JJJ-1, III-1, devSim.nbasis, devSim.nbasis);

                        // Find the (ij|kl) integrals where j>i, k>i, l>k, and k and j are equal.
                        QUICKULL val1 = (QUICKULL) (fabs(2.0*DENSELK*Y*OSCALE) + (QUICKDouble)0.5);
                        if ( DENSELK*Y < (QUICKDouble)0.0)
                        val1 = 0ull - val1;
                        
                        QUICKULL val2 = (QUICKULL) (fabs(2.0*DENSEJI*Y*OSCALE) + (QUICKDouble)0.5);
                        if ( DENSEJI*Y < (QUICKDouble)0.0)
                        val2 = 0ull - val2;
                        
                        QUICKULL val3 = (QUICKULL) (fabs(0.5*DENSELJ*Y*OSCALE) + (QUICKDouble)0.5);
                        if ( DENSELJ*Y < (QUICKDouble)0.0)
                        val3 = 0ull - val3;
                        
                        QUICKULL val4 = (QUICKULL) (fabs(0.5*DENSEKJ*Y*OSCALE) + (QUICKDouble)0.5);
                        if ( DENSEKJ*Y < (QUICKDouble)0.0)
                        val4 = 0ull - val4;
                        
                        QUICKULL val5 = (QUICKULL) (fabs(0.5*DENSELI*Y*OSCALE) + (QUICKDouble)0.5);
                        if ( DENSELI*Y < (QUICKDouble)0.0)
                        val5 = 0ull - val5;
                        
                        QUICKULL val6 = (QUICKULL) (fabs(0.5*DENSEKI*Y*OSCALE) + (QUICKDouble)0.5);
                        if ( DENSEKI*Y < (QUICKDouble)0.0)
                        val6 = 0ull - val6;
                        
                        QUICKADD(LOC2(devSim.oULL, JJJ-1, III-1, devSim.nbasis, devSim.nbasis), val1);
                        QUICKADD(LOC2(devSim.oULL, LLL-1, KKK-1, devSim.nbasis, devSim.nbasis), val2);
                        QUICKADD(LOC2(devSim.oULL, KKK-1, III-1, devSim.nbasis, devSim.nbasis), 0ull-val3);
                        QUICKADD(LOC2(devSim.oULL, LLL-1, III-1, devSim.nbasis, devSim.nbasis), 0ull-val4);
                        QUICKADD(LOC2(devSim.oULL, JJJ-1, KKK-1, devSim.nbasis, devSim.nbasis), 0ull-val5);
                        QUICKADD(LOC2(devSim.oULL, KKK-1, JJJ-1, devSim.nbasis, devSim.nbasis), 0ull-val5);
                        QUICKADD(LOC2(devSim.oULL, JJJ-1, LLL-1, devSim.nbasis, devSim.nbasis), 0ull-val6);
                        QUICKADD(LOC2(devSim.oULL, LLL-1, JJJ-1, devSim.nbasis, devSim.nbasis), 0ull-val6);               
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
                        QUICKDouble Y = (QUICKDouble) hrrwhole(III, JJJ, KKK, LLL, IJKLTYPE, store, \
                                                               RAx, RAy, RAz, RBx, RBy, RBz, \
                                                               RCx, RCy, RCz, RDx, RDy, RDz);

                            if ((III < JJJ)&&(KKK < LLL)) {
                                QUICKDouble DENSEKI = (QUICKDouble) LOC2(devSim.dense, KKK-1, III-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEKJ = (QUICKDouble) LOC2(devSim.dense, KKK-1, JJJ-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSELJ = (QUICKDouble) LOC2(devSim.dense, LLL-1, JJJ-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSELI = (QUICKDouble) LOC2(devSim.dense, LLL-1, III-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSELK = (QUICKDouble) LOC2(devSim.dense, LLL-1, KKK-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEJI = (QUICKDouble) LOC2(devSim.dense, JJJ-1, III-1, devSim.nbasis, devSim.nbasis);
                                // Find the (ij|kl) integrals where j>i, k>i, l>k, and k and j are equal.
                                QUICKULL val1 = (QUICKULL) (fabs(2.0*DENSELK*Y*OSCALE) + (QUICKDouble)0.5);
                                if ( DENSELK*Y < (QUICKDouble)0.0)
                                val1 = 0ull - val1;
                                
                                QUICKULL val2 = (QUICKULL) (fabs(2.0*DENSEJI*Y*OSCALE) + (QUICKDouble)0.5);
                                if ( DENSEJI*Y < (QUICKDouble)0.0)
                                val2 = 0ull - val2;

                                QUICKULL val3 = (QUICKULL) (fabs(0.5*DENSELJ*Y*OSCALE) + (QUICKDouble)0.5);
                                if ( DENSELJ*Y < (QUICKDouble)0.0)
                                val3 = 0ull - val3;
                                
                                QUICKULL val4 = (QUICKULL) (fabs(0.5*DENSEKJ*Y*OSCALE) + (QUICKDouble)0.5);
                                if ( DENSEKJ*Y < (QUICKDouble)0.0)
                                val4 = 0ull - val4;
                                
                                QUICKULL val5 = (QUICKULL) (fabs(0.5*DENSELI*Y*OSCALE) + (QUICKDouble)0.5);
                                if ( DENSELI*Y < (QUICKDouble)0.0)
                                val5 = 0ull - val5;
                                
                                QUICKULL val6 = (QUICKULL) (fabs(0.5*DENSEKI*Y*OSCALE) + (QUICKDouble)0.5);
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
                                
                            }else if ((III == JJJ)&&(KKK == LLL)) {
                                
                                // Find  all the (ii|jj) integrals.
                                QUICKDouble DENSEJI = LOC2(devSim.dense, KKK-1, III-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEJJ = LOC2(devSim.dense, KKK-1, KKK-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEII = LOC2(devSim.dense, III-1, III-1, devSim.nbasis, devSim.nbasis);

                                QUICKULL val1 = (QUICKULL) (fabs(DENSEJJ*Y*OSCALE) + (QUICKDouble)0.5);
                                if ( DENSEJJ*Y < (QUICKDouble)0.0)
                                val1 = 0ull - val1;                               
                                
                                QUICKULL val2 = (QUICKULL) (fabs(DENSEII*Y*OSCALE) + (QUICKDouble)0.5);
                                if ( DENSEII*Y < (QUICKDouble)0.0)
                                val2 = 0ull - val2;
                                
                                QUICKULL val3 = (QUICKULL) (fabs(0.5*DENSEJI*Y*OSCALE) + (QUICKDouble)0.5);
                                if ( DENSEJI*Y < (QUICKDouble)0.0)
                                val3 = 0ull - val3;
                                
                                QUICKADD(LOC2(devSim.oULL, III-1, III-1, devSim.nbasis, devSim.nbasis), val1);
                                QUICKADD(LOC2(devSim.oULL, KKK-1, KKK-1, devSim.nbasis, devSim.nbasis), val2);
                                QUICKADD(LOC2(devSim.oULL, KKK-1, III-1, devSim.nbasis, devSim.nbasis), 0ull-val3);

                            }else if ((JJJ == KKK)&&(JJJ==LLL)) {

                                // Find all the (ij|jj) integrals.
                                QUICKDouble DENSEJI = (QUICKDouble) LOC2(devSim.dense, JJJ-1, III-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEJJ = (QUICKDouble) LOC2(devSim.dense, JJJ-1, JJJ-1, devSim.nbasis, devSim.nbasis);

                                QUICKULL val1 = (QUICKULL) (fabs(0.5*DENSEJJ*Y*OSCALE) + (QUICKDouble)0.5);
                                if ( DENSEJJ*Y < (QUICKDouble)0.0)
                                val1 = 0ull - val1;                               
                                
                                QUICKULL val2 = (QUICKULL) (fabs(DENSEJI*Y*OSCALE) + (QUICKDouble)0.5);
                                if ( DENSEJI*Y < (QUICKDouble)0.0)
                                val2 = 0ull - val2;
                                
                                QUICKADD(LOC2(devSim.oULL, JJJ-1, III-1, devSim.nbasis, devSim.nbasis), val1);
                                QUICKADD(LOC2(devSim.oULL, JJJ-1, JJJ-1, devSim.nbasis, devSim.nbasis), val2);
                            
                            }else if ((KKK == LLL)&&(III<JJJ)&&(JJJ!=KKK)) {
                                
                                //Find all the (ij|kk) integrals where j>i, k>j.
                                QUICKDouble DENSEKI = (QUICKDouble) LOC2(devSim.dense, KKK-1, III-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEKJ = (QUICKDouble) LOC2(devSim.dense, KKK-1, JJJ-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEKK = (QUICKDouble) LOC2(devSim.dense, KKK-1, KKK-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEJI = (QUICKDouble) LOC2(devSim.dense, JJJ-1, III-1, devSim.nbasis, devSim.nbasis);
 
                                QUICKULL val1 = (QUICKULL) (fabs(DENSEKK*Y*OSCALE) + (QUICKDouble)0.5);
                                if ( DENSEKK*Y < (QUICKDouble)0.0)
                                val1 = 0ull - val1;                               
                                
                                QUICKULL val2 = (QUICKULL) (fabs(2.0*DENSEJI*Y*OSCALE) + (QUICKDouble)0.5);
                                if ( DENSEJI*Y < (QUICKDouble)0.0)
                                val2 = 0ull - val2;
                                
                                QUICKULL val3 = (QUICKULL) (fabs(0.5*DENSEKJ*Y*OSCALE) + (QUICKDouble)0.5);
                                if ( DENSEKJ*Y < (QUICKDouble)0.0)
                                val3 = 0ull - val3;
                                
                                QUICKULL val4 = (QUICKULL) (fabs(0.5*DENSEKI*Y*OSCALE) + (QUICKDouble)0.5);
                                if ( DENSEKI*Y < (QUICKDouble)0.0)
                                val4 = 0ull - val4;

                                QUICKADD(LOC2(devSim.oULL, JJJ-1, III-1, devSim.nbasis, devSim.nbasis), val1);
                                QUICKADD(LOC2(devSim.oULL, KKK-1, KKK-1, devSim.nbasis, devSim.nbasis), val2);
                                QUICKADD(LOC2(devSim.oULL, KKK-1, III-1, devSim.nbasis, devSim.nbasis), 0ull-val3);
                                QUICKADD(LOC2(devSim.oULL, KKK-1, JJJ-1, devSim.nbasis, devSim.nbasis), 0ull-val4);
                                QUICKADD(LOC2(devSim.oULL, JJJ-1, KKK-1, devSim.nbasis, devSim.nbasis), 0ull-val4);

                            }else if ((III==JJJ)&&(KKK<LLL)) {
                                
                                //Find all the (ii|jk) integrals where j>i, k>j.
                                QUICKDouble DENSEII = (QUICKDouble) LOC2(devSim.dense, III-1, III-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEJI = (QUICKDouble) LOC2(devSim.dense, KKK-1, III-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEKI = (QUICKDouble) LOC2(devSim.dense, LLL-1, III-1, devSim.nbasis, devSim.nbasis);
                                QUICKDouble DENSEKJ = (QUICKDouble) LOC2(devSim.dense, LLL-1, KKK-1, devSim.nbasis, devSim.nbasis);
                                
                                QUICKULL val1 = (QUICKULL) (fabs(DENSEII*Y*OSCALE) + (QUICKDouble)0.5);
                                if ( DENSEII*Y < (QUICKDouble)0.0)
                                val1 = 0ull - val1;                               
                                
                                QUICKULL val2 = (QUICKULL) (fabs(2.0*DENSEKJ*Y*OSCALE) + (QUICKDouble)0.5);
                                if ( DENSEKJ*Y < (QUICKDouble)0.0)
                                val2 = 0ull - val2;
                                
                                QUICKULL val3 = (QUICKULL) (fabs(0.5*DENSEKI*Y*OSCALE) + (QUICKDouble)0.5);
                                if ( DENSEKI*Y < (QUICKDouble)0.0)
                                val3 = 0ull - val3;
                                
                                QUICKULL val4 = (QUICKULL) (fabs(0.5*DENSEJI*Y*OSCALE) + (QUICKDouble)0.5);
                                if ( DENSEJI*Y < (QUICKDouble)0.0)
                                val4 = 0ull - val4;

                                QUICKADD(LOC2(devSim.oULL, LLL-1, KKK-1, devSim.nbasis, devSim.nbasis), val1);
                                QUICKADD(LOC2(devSim.oULL, III-1, III-1, devSim.nbasis, devSim.nbasis), val2);
                                QUICKADD(LOC2(devSim.oULL, KKK-1, III-1, devSim.nbasis, devSim.nbasis), 0ull-val3);
                                QUICKADD(LOC2(devSim.oULL, LLL-1, III-1, devSim.nbasis, devSim.nbasis), 0ull-val4);
                            } 
                            
                        }else {
                        
                            if (JJJ <= LLL) {
                                QUICKDouble Y = (QUICKDouble) hrrwhole(III, JJJ, KKK, LLL, IJKLTYPE, store, \
                                                               RAx, RAy, RAz, RBx, RBy, RBz, \
                                                               RCx, RCy, RCz, RDx, RDy, RDz);

                                if((III==JJJ)&&(III==KKK)&&(III==LLL)){
                                    // do all the (ii|ii) integrals
                                    QUICKDouble DENSEII = (QUICKDouble) LOC2(devSim.dense, III-1, III-1, devSim.nbasis, devSim.nbasis);
                                    QUICKULL val1 = (QUICKULL) (fabs(0.5*DENSEII*Y*OSCALE) + (QUICKDouble)0.5);
                                    if ( DENSEII*Y < (QUICKDouble)0.0)
                                    val1 = 0ull - val1;                               
                                    QUICKADD(LOC2(devSim.oULL, III-1, III-1, devSim.nbasis, devSim.nbasis), val1);
                                }else if ((III == JJJ) && (III == KKK) && (III < LLL)){
                                    
                                    QUICKDouble DENSEJI = (QUICKDouble) LOC2(devSim.dense, LLL-1, III-1, devSim.nbasis, devSim.nbasis);
                                    QUICKDouble DENSEII = (QUICKDouble) LOC2(devSim.dense, III-1, III-1, devSim.nbasis, devSim.nbasis);
                                                                      

                                    QUICKULL val1 = (QUICKULL) (fabs(0.5*DENSEII*Y*OSCALE) + (QUICKDouble)0.5);
                                    if ( DENSEII*Y < (QUICKDouble)0.0)
                                    val1 = 0ull - val1;                               
                                    
                                    QUICKULL val2 = (QUICKULL) (fabs(DENSEJI*Y*OSCALE) + (QUICKDouble)0.5);
                                    if ( DENSEJI*Y < (QUICKDouble)0.0)
                                    val2 = 0ull - val2;                               
                                    
                                    QUICKADD(LOC2(devSim.oULL, LLL-1, III-1, devSim.nbasis, devSim.nbasis), val1);
                                    QUICKADD(LOC2(devSim.oULL, III-1, III-1, devSim.nbasis, devSim.nbasis), val2);

                                }else if ((III == KKK) && (JJJ == LLL) && (III < JJJ)){
                                    
                                    QUICKDouble DENSEJI = (QUICKDouble) LOC2(devSim.dense, JJJ-1, III-1, devSim.nbasis, devSim.nbasis);
                                    QUICKDouble DENSEJJ = (QUICKDouble) LOC2(devSim.dense, JJJ-1, JJJ-1, devSim.nbasis, devSim.nbasis);
                                    QUICKDouble DENSEII = (QUICKDouble) LOC2(devSim.dense, III-1, III-1, devSim.nbasis, devSim.nbasis);
                                    QUICKULL val1 = (QUICKULL) (fabs(1.5*DENSEJI*Y*OSCALE) + (QUICKDouble)0.5);
                                    if ( DENSEJI*Y < (QUICKDouble)0.0)
                                    val1 = 0ull - val1;                               
                                    
                                    QUICKULL val2 = (QUICKULL) (fabs(0.5*DENSEII*Y*OSCALE) + (QUICKDouble)0.5);
                                    if ( DENSEII*Y < (QUICKDouble)0.0)
                                    val2 = 0ull - val2;                               
                                    
                                    QUICKULL val3 = (QUICKULL) (fabs(0.5*DENSEJJ*Y*OSCALE) + (QUICKDouble)0.5);
                                    if ( DENSEJJ*Y < (QUICKDouble)0.0)
                                    val3 = 0ull - val3;                               
                                    
                                    
                                    QUICKADD(LOC2(devSim.oULL, JJJ-1, III-1, devSim.nbasis, devSim.nbasis), val1);
                                    QUICKADD(LOC2(devSim.oULL, JJJ-1, JJJ-1, devSim.nbasis, devSim.nbasis), 0ull-val2);
                                    QUICKADD(LOC2(devSim.oULL, III-1, III-1, devSim.nbasis, devSim.nbasis), 0ull-val3);
                                }else if ((III == KKK) && (III <  JJJ) && (JJJ < LLL)){
                                    
                                    QUICKDouble DENSEKI = (QUICKDouble) LOC2(devSim.dense, LLL-1, III-1, devSim.nbasis, devSim.nbasis);
                                    QUICKDouble DENSEKJ = (QUICKDouble) LOC2(devSim.dense, LLL-1, JJJ-1, devSim.nbasis, devSim.nbasis);
                                    QUICKDouble DENSEII = (QUICKDouble) LOC2(devSim.dense, III-1, III-1, devSim.nbasis, devSim.nbasis);
                                    QUICKDouble DENSEJI = (QUICKDouble) LOC2(devSim.dense, JJJ-1, III-1, devSim.nbasis, devSim.nbasis);
                                   
                                    QUICKULL val1 = (QUICKULL) (fabs(1.5*DENSEKI*Y*OSCALE) + (QUICKDouble)0.5);
                                    if ( DENSEKI*Y < (QUICKDouble)0.0)
                                    val1 = 0ull - val1;                               
                                    
                                    QUICKULL val2 = (QUICKULL) (fabs(1.5*DENSEJI*Y*OSCALE) + (QUICKDouble)0.5);
                                    if ( DENSEJI*Y < (QUICKDouble)0.0)
                                    val2 = 0ull - val2;                               
                                    
                                    QUICKULL val3 = (QUICKULL) (fabs(1.0*DENSEKJ*Y*OSCALE) + (QUICKDouble)0.5);
                                    if ( DENSEKJ*Y < (QUICKDouble)0.0)
                                    val3 = 0ull - val3;                               
                                    
                                    QUICKULL val4 = (QUICKULL) (fabs(0.5*DENSEII*Y*OSCALE) + (QUICKDouble)0.5);
                                    if ( DENSEII*Y < (QUICKDouble)0.0)
                                    val4 = 0ull - val4;
                                    QUICKADD(LOC2(devSim.oULL, JJJ-1, III-1, devSim.nbasis, devSim.nbasis), val1);
                                    QUICKADD(LOC2(devSim.oULL, LLL-1, III-1, devSim.nbasis, devSim.nbasis), val2);
                                    QUICKADD(LOC2(devSim.oULL, III-1, III-1, devSim.nbasis, devSim.nbasis), 0ull-val3);
                                    QUICKADD(LOC2(devSim.oULL, LLL-1, JJJ-1, devSim.nbasis, devSim.nbasis), 0ull-val4);
								}
							}
                        }
                    }
                }
            }
        }
    }
	return;
}


#ifndef TEST
__device__
#endif
  void FmT(int MaxM, QUICKDouble X, QUICKDouble* YVerticalTemp)
{
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
        LOC3(YVerticalTemp, 0, 0, 0, VDIM1, VDIM2, VDIM3) = WW1;
        for (int m = 1; m<= MaxM; m++) {
            LOC3(YVerticalTemp, 0, 0, m, VDIM1, VDIM2, VDIM3) = (((2*m-1)*LOC3(YVerticalTemp, 0, 0, m-1, VDIM1, VDIM2, VDIM3))- E)*0.5*XINV;
        }
    }else {
        LOC3(YVerticalTemp, 0, 0, MaxM, VDIM1, VDIM2, VDIM3) = WW1;
        for (int m = MaxM-1; m >=0; m--) {
            LOC3(YVerticalTemp, 0, 0, m, VDIM1, VDIM2, VDIM3) = (2.0 * X * LOC3(YVerticalTemp, 0, 0, m+1, VDIM1, VDIM2, VDIM3) + E) / (QUICKDouble)(m*2+1);
        }
    }
	return;
}


#ifndef TEST
__device__
#endif
void PSSS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Ptempx, QUICKDouble Ptempy,QUICKDouble Ptempz, \
          QUICKDouble WPtempx, QUICKDouble WPtempy, QUICKDouble WPtempz)
{
    VY( 1, 0, mtemp) = Ptempx * VY( 0, 0, mtemp) + WPtempx * VY( 0, 0, mtemp+1);
    VY( 2, 0, mtemp) = Ptempy * VY( 0, 0, mtemp) + WPtempy * VY( 0, 0, mtemp+1);
    VY( 3, 0, mtemp) = Ptempz * VY( 0, 0, mtemp) + WPtempz * VY( 0, 0, mtemp+1);
	return;
}

#ifndef TEST
__device__
#endif
void SSPS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Qtempx, QUICKDouble Qtempy,QUICKDouble Qtempz, \
          QUICKDouble WQtempx, QUICKDouble WQtempy, QUICKDouble WQtempz)
{
    
    VY( 0, 1, mtemp) = Qtempx * VY( 0, 0, mtemp) + WQtempx * VY( 0, 0, mtemp+1);
    VY( 0, 2, mtemp) = Qtempy * VY( 0, 0, mtemp) + WQtempy * VY( 0, 0, mtemp+1);
    VY( 0, 3, mtemp) = Qtempz * VY( 0, 0, mtemp) + WQtempz * VY( 0, 0, mtemp+1);
    
	return;
}

#ifndef TEST
__device__
#endif
void PSPS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Ptempx, QUICKDouble Ptempy,QUICKDouble Ptempz, \
          QUICKDouble WPtempx, QUICKDouble WPtempy, QUICKDouble WPtempz, QUICKDouble ABCDtemp)
{
    
    VY( 1, 1, mtemp) = Ptempx * VY( 0, 1, mtemp) + WPtempx * VY( 0, 1, mtemp+1) + ABCDtemp * VY( 0, 0, mtemp+1);
    VY( 2, 1, mtemp) = Ptempy * VY( 0, 1, mtemp) + WPtempy * VY( 0, 1, mtemp+1);
    VY( 3, 1, mtemp) = Ptempz * VY( 0, 1, mtemp) + WPtempz * VY( 0, 1, mtemp+1);
    
    VY( 1, 2, mtemp) = Ptempx * VY( 0, 2, mtemp) + WPtempx * VY( 0, 2, mtemp+1);
    VY( 2, 2, mtemp) = Ptempy * VY( 0, 2, mtemp) + WPtempy * VY( 0, 2, mtemp+1) + ABCDtemp * VY( 0, 0, mtemp+1);
    VY( 3, 2, mtemp) = Ptempz * VY( 0, 2, mtemp) + WPtempz * VY( 0, 2, mtemp+1);
    
    VY( 1, 3, mtemp) = Ptempx * VY( 0, 3, mtemp) + WPtempx * VY( 0, 3, mtemp+1);
    VY( 2, 3, mtemp) = Ptempy * VY( 0, 3, mtemp) + WPtempy * VY( 0, 3, mtemp+1);
    VY( 3, 3, mtemp) = Ptempz * VY( 0, 3, mtemp) + WPtempz * VY( 0, 3, mtemp+1) + ABCDtemp * VY( 0, 0, mtemp+1);
    
	return;
}

#ifndef TEST
__device__
#endif
void DSSS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Ptempx, QUICKDouble Ptempy,QUICKDouble Ptempz, \
          QUICKDouble WPtempx, QUICKDouble WPtempy, QUICKDouble WPtempz, QUICKDouble ABtemp, QUICKDouble CDcom)
{
    VY( 4, 0, mtemp) = Ptempx * VY( 2, 0, mtemp) + WPtempx * VY( 2, 0, mtemp+1);
    VY( 5, 0, mtemp) = Ptempy * VY( 3, 0, mtemp) + WPtempy * VY( 3, 0, mtemp+1);
    VY( 6, 0, mtemp) = Ptempx * VY( 3, 0, mtemp) + WPtempx * VY( 3, 0, mtemp+1);
    
    VY( 7, 0, mtemp) = Ptempx * VY( 1, 0, mtemp) + WPtempx * VY( 1, 0, mtemp+1)+ ABtemp*(VY( 0, 0, mtemp) - CDcom * VY( 0, 0, mtemp+1));
    VY( 8, 0, mtemp) = Ptempy * VY( 2, 0, mtemp) + WPtempy * VY( 2, 0, mtemp+1)+ ABtemp*(VY( 0, 0, mtemp) - CDcom * VY( 0, 0, mtemp+1));
    VY( 9, 0, mtemp) = Ptempz * VY( 3, 0, mtemp) + WPtempz * VY( 3, 0, mtemp+1)+ ABtemp*(VY( 0, 0, mtemp) - CDcom * VY( 0, 0, mtemp+1));
    
	return;
}

#ifndef TEST
__device__
#endif
void SSDS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Qtempx, QUICKDouble Qtempy,QUICKDouble Qtempz, \
          QUICKDouble WQtempx, QUICKDouble WQtempy, QUICKDouble WQtempz, QUICKDouble CDtemp, QUICKDouble ABcom)
{
    VY( 0, 4, mtemp) = Qtempx * VY( 0, 2, mtemp) + WQtempx * VY( 0, 2, mtemp+1);
    VY( 0, 5, mtemp) = Qtempy * VY( 0, 3, mtemp) + WQtempy * VY( 0, 3, mtemp+1);
    VY( 0, 6, mtemp) = Qtempx * VY( 0, 3, mtemp) + WQtempx * VY( 0, 3, mtemp+1);
    
    VY( 0, 7, mtemp) = Qtempx * VY( 0, 1, mtemp) + WQtempx * VY( 0, 1, mtemp+1)+ CDtemp*(VY( 0, 0, mtemp) - ABcom * VY( 0, 0, mtemp+1));
    VY( 0, 8, mtemp) = Qtempy * VY( 0, 2, mtemp) + WQtempy * VY( 0, 2, mtemp+1)+ CDtemp*(VY( 0, 0, mtemp) - ABcom * VY( 0, 0, mtemp+1));
    VY( 0, 9, mtemp) = Qtempz * VY( 0, 3, mtemp) + WQtempz * VY( 0, 3, mtemp+1)+ CDtemp*(VY( 0, 0, mtemp) - ABcom * VY( 0, 0, mtemp+1));
    
	return;
}



#ifndef TEST
__device__
#endif
void DSPS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Qtempx, QUICKDouble Qtempy,QUICKDouble Qtempz, \
          QUICKDouble WQtempx, QUICKDouble WQtempy, QUICKDouble WQtempz, QUICKDouble ABCDtemp)
{
    
    VY( 4, 1, mtemp) = Qtempx * VY( 4, 0, mtemp) + WQtempx * VY( 4, 0, mtemp + 1) + ABCDtemp * VY( 2, 0, mtemp + 1);
    VY( 4, 2, mtemp) = Qtempy * VY( 4, 0, mtemp) + WQtempy * VY( 4, 0, mtemp + 1) + ABCDtemp * VY( 1, 0, mtemp + 1);
    VY( 4, 3, mtemp) = Qtempz * VY( 4, 0, mtemp) + WQtempz * VY( 4, 0, mtemp + 1);
    
    VY( 5, 1, mtemp) = Qtempx * VY( 5, 0, mtemp) + WQtempx * VY( 5, 0, mtemp + 1);
    VY( 5, 2, mtemp) = Qtempy * VY( 5, 0, mtemp) + WQtempy * VY( 5, 0, mtemp + 1) + ABCDtemp * VY( 3, 0, mtemp + 1);
    VY( 5, 3, mtemp) = Qtempz * VY( 5, 0, mtemp) + WQtempz * VY( 5, 0, mtemp + 1) + ABCDtemp * VY( 2, 0, mtemp + 1);
    
    VY( 6, 1, mtemp) = Qtempx * VY( 6, 0, mtemp) + WQtempx * VY( 6, 0, mtemp + 1) + ABCDtemp * VY( 3, 0, mtemp + 1);
    VY( 6, 2, mtemp) = Qtempy * VY( 6, 0, mtemp) + WQtempy * VY( 6, 0, mtemp + 1);
    VY( 6, 3, mtemp) = Qtempz * VY( 6, 0, mtemp) + WQtempz * VY( 6, 0, mtemp + 1) + ABCDtemp * VY( 1, 0, mtemp + 1);
    
    VY( 7, 1, mtemp) = Qtempx * VY( 7, 0, mtemp) + WQtempx * VY( 7, 0, mtemp + 1) + ABCDtemp * VY( 1, 0, mtemp + 1) * 2;
    VY( 7, 2, mtemp) = Qtempy * VY( 7, 0, mtemp) + WQtempy * VY( 7, 0, mtemp + 1);
    VY( 7, 3, mtemp) = Qtempz * VY( 7, 0, mtemp) + WQtempz * VY( 7, 0, mtemp + 1);
    
    VY( 8, 1, mtemp) = Qtempx * VY( 8, 0, mtemp) + WQtempx * VY( 8, 0, mtemp + 1);
    VY( 8, 2, mtemp) = Qtempy * VY( 8, 0, mtemp) + WQtempy * VY( 8, 0, mtemp + 1) + ABCDtemp * VY( 2, 0, mtemp + 1) * 2;
    VY( 8, 3, mtemp) = Qtempz * VY( 8, 0, mtemp) + WQtempz * VY( 8, 0, mtemp + 1);
    
    VY( 9, 1, mtemp) = Qtempx * VY( 9, 0, mtemp) + WQtempx * VY( 9, 0, mtemp + 1);
    VY( 9, 2, mtemp) = Qtempy * VY( 9, 0, mtemp) + WQtempy * VY( 9, 0, mtemp + 1);
    VY( 9, 3, mtemp) = Qtempz * VY( 9, 0, mtemp) + WQtempz * VY( 9, 0, mtemp + 1) + ABCDtemp * VY( 3, 0, mtemp + 1) * 2;            
    
	return;
}

#ifndef TEST
__device__
#endif
void PSDS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Ptempx, QUICKDouble Ptempy,QUICKDouble Ptempz, \
          QUICKDouble WPtempx, QUICKDouble WPtempy, QUICKDouble WPtempz, QUICKDouble ABCDtemp)
{
    
    
    VY( 1, 4, mtemp) = Ptempx * VY( 0, 4, mtemp) + WPtempx * VY( 0, 4, mtemp + 1) + ABCDtemp * VY( 0, 2, mtemp + 1);
    VY( 2, 4, mtemp) = Ptempy * VY( 0, 4, mtemp) + WPtempy * VY( 0, 4, mtemp + 1) + ABCDtemp * VY( 0, 1, mtemp + 1);
    VY( 3, 4, mtemp) = Ptempz * VY( 0, 4, mtemp) + WPtempz * VY( 0, 4, mtemp + 1);
    
    VY( 1, 5, mtemp) = Ptempx * VY( 0, 5, mtemp) + WPtempx * VY( 0, 5, mtemp + 1);
    VY( 2, 5, mtemp) = Ptempy * VY( 0, 5, mtemp) + WPtempy * VY( 0, 5, mtemp + 1) + ABCDtemp * VY( 0, 3, mtemp + 1);
    VY( 3, 5, mtemp) = Ptempz * VY( 0, 5, mtemp) + WPtempz * VY( 0, 5, mtemp + 1) + ABCDtemp * VY( 0, 2, mtemp + 1);
    
    VY( 1, 6, mtemp) = Ptempx * VY( 0, 6, mtemp) + WPtempx * VY( 0, 6, mtemp + 1) + ABCDtemp * VY( 0, 3, mtemp + 1);
    VY( 2, 6, mtemp) = Ptempy * VY( 0, 6, mtemp) + WPtempy * VY( 0, 6, mtemp + 1);
    VY( 3, 6, mtemp) = Ptempz * VY( 0, 6, mtemp) + WPtempz * VY( 0, 6, mtemp + 1) + ABCDtemp * VY( 0, 1, mtemp + 1);
    
    VY( 1, 7, mtemp) = Ptempx * VY( 0, 7, mtemp) + WPtempx * VY( 0, 7, mtemp + 1) + ABCDtemp * VY( 0, 1, mtemp + 1) * 2;
    VY( 2, 7, mtemp) = Ptempy * VY( 0, 7, mtemp) + WPtempy * VY( 0, 7, mtemp + 1);
    VY( 3, 7, mtemp) = Ptempz * VY( 0, 7, mtemp) + WPtempz * VY( 0, 7, mtemp + 1);
    
    VY( 1, 8, mtemp) = Ptempx * VY( 0, 8, mtemp) + WPtempx * VY( 0, 8, mtemp + 1);
    VY( 2, 8, mtemp) = Ptempy * VY( 0, 8, mtemp) + WPtempy * VY( 0, 8, mtemp + 1) + ABCDtemp * VY( 0, 2, mtemp + 1) * 2;
    VY( 3, 8, mtemp) = Ptempz * VY( 0, 8, mtemp) + WPtempz * VY( 0, 8, mtemp + 1);
    
    VY( 1, 9, mtemp) = Ptempx * VY( 0, 9, mtemp) + WPtempx * VY( 0, 9, mtemp + 1);
    VY( 2, 9, mtemp) = Ptempy * VY( 0, 9, mtemp) + WPtempy * VY( 0, 9, mtemp + 1);
    VY( 3, 9, mtemp) = Ptempz * VY( 0, 9, mtemp) + WPtempz * VY( 0, 9, mtemp + 1) + ABCDtemp * VY( 0, 3, mtemp + 1) * 2;    
    
	return;
}

#ifndef TEST
__device__
#endif
void DSDS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Ptempx, QUICKDouble Ptempy,QUICKDouble Ptempz, \
          QUICKDouble WPtempx, QUICKDouble WPtempy, QUICKDouble WPtempz, QUICKDouble ABCDtemp, QUICKDouble ABtemp, QUICKDouble CDcom)
{
    
    VY( 4, 4, mtemp) = Ptempx * VY( 2, 4, mtemp) + WPtempx * VY( 2, 4, mtemp+1) + ABCDtemp * VY( 2, 2, mtemp+1);
    VY( 4, 5, mtemp) = Ptempx * VY( 2, 5, mtemp) + WPtempx * VY( 2, 5, mtemp+1);
    VY( 4, 6, mtemp) = Ptempx * VY( 2, 6, mtemp) + WPtempx * VY( 2, 6, mtemp+1) + ABCDtemp * VY( 2, 3, mtemp+1);
    VY( 4, 7, mtemp) = Ptempx * VY( 2, 7, mtemp) + WPtempx * VY( 2, 7, mtemp+1) + 2 * ABCDtemp * VY( 2, 1, mtemp+1);
    VY( 4, 8, mtemp) = Ptempx * VY( 2, 8, mtemp) + WPtempx * VY( 2, 8, mtemp+1);
    VY( 4, 9, mtemp) = Ptempx * VY( 2, 9, mtemp) + WPtempx * VY( 2, 9, mtemp+1);
    
    VY( 5, 4, mtemp) = Ptempy * VY( 3, 4, mtemp) + WPtempy * VY( 3, 4, mtemp+1) + ABCDtemp * VY( 3, 1, mtemp+1);
    VY( 5, 5, mtemp) = Ptempy * VY( 3, 5, mtemp) + WPtempy * VY( 3, 5, mtemp+1) + ABCDtemp * VY( 3, 3, mtemp+1);
    VY( 5, 6, mtemp) = Ptempy * VY( 3, 6, mtemp) + WPtempy * VY( 3, 6, mtemp+1);
    VY( 5, 7, mtemp) = Ptempy * VY( 3, 7, mtemp) + WPtempy * VY( 3, 7, mtemp+1);
    VY( 5, 8, mtemp) = Ptempy * VY( 3, 8, mtemp) + WPtempy * VY( 3, 8, mtemp+1) + 2 * ABCDtemp * VY( 3, 2, mtemp+1);
    VY( 5, 9, mtemp) = Ptempy * VY( 3, 9, mtemp) + WPtempy * VY( 3, 9, mtemp+1);
    
    VY( 6, 4, mtemp) = Ptempx * VY( 3, 4, mtemp) + WPtempx * VY( 3, 4, mtemp+1) + ABCDtemp * VY( 3, 2, mtemp+1);
    VY( 6, 5, mtemp) = Ptempx * VY( 3, 5, mtemp) + WPtempx * VY( 3, 5, mtemp+1);
    VY( 6, 6, mtemp) = Ptempx * VY( 3, 6, mtemp) + WPtempx * VY( 3, 6, mtemp+1) + ABCDtemp * VY( 3, 3, mtemp+1);
    VY( 6, 7, mtemp) = Ptempx * VY( 3, 7, mtemp) + WPtempx * VY( 3, 7, mtemp+1) + 2 * ABCDtemp * VY( 3, 1, mtemp+1);
    VY( 6, 8, mtemp) = Ptempx * VY( 3, 8, mtemp) + WPtempx * VY( 3, 8, mtemp+1);
    VY( 6, 9, mtemp) = Ptempx * VY( 3, 9, mtemp) + WPtempx * VY( 3, 9, mtemp+1);
    
    VY( 7, 4, mtemp) = Ptempx * VY( 1, 4, mtemp) + WPtempx * VY( 1, 4, mtemp+1) +  ABtemp * (VY( 0, 4,mtemp)-CDcom*VY( 0, 4,mtemp+1)) + ABCDtemp * VY( 1, 2, mtemp+1);
    VY( 7, 5, mtemp) = Ptempx * VY( 1, 5, mtemp) + WPtempx * VY( 1, 5, mtemp+1) +  ABtemp * (VY( 0, 5,mtemp)-CDcom*VY( 0, 5,mtemp+1));
    VY( 7, 6, mtemp) = Ptempx * VY( 1, 6, mtemp) + WPtempx * VY( 1, 6, mtemp+1) +  ABtemp * (VY( 0, 6,mtemp)-CDcom*VY( 0, 6,mtemp+1)) + ABCDtemp * VY( 1, 3, mtemp+1);
    VY( 7, 7, mtemp) = Ptempx * VY( 1, 7, mtemp) + WPtempx * VY( 1, 7, mtemp+1) +  ABtemp * (VY( 0, 7,mtemp)-CDcom*VY( 0, 7,mtemp+1)) + 2 * ABCDtemp * VY( 1, 1, mtemp+1);
    VY( 7, 8, mtemp) = Ptempx * VY( 1, 8, mtemp) + WPtempx * VY( 1, 8, mtemp+1) +  ABtemp * (VY( 0, 8,mtemp)-CDcom*VY( 0, 8,mtemp+1));
    VY( 7, 9, mtemp) = Ptempx * VY( 1, 9, mtemp) + WPtempx * VY( 1, 9, mtemp+1) +  ABtemp * (VY( 0, 9,mtemp)-CDcom*VY( 0, 9,mtemp+1));
    
    
    VY( 8, 4, mtemp) = Ptempy * VY( 2, 4, mtemp) + WPtempy * VY( 2, 4, mtemp+1) +  ABtemp * (VY( 0, 4,mtemp)-CDcom*VY( 0, 4,mtemp+1)) + ABCDtemp * VY( 2, 1, mtemp+1);
    VY( 8, 5, mtemp) = Ptempy * VY( 2, 5, mtemp) + WPtempy * VY( 2, 5, mtemp+1) +  ABtemp * (VY( 0, 5,mtemp)-CDcom*VY( 0, 5,mtemp+1)) + ABCDtemp * VY( 2, 3, mtemp+1);
    VY( 8, 6, mtemp) = Ptempy * VY( 2, 6, mtemp) + WPtempy * VY( 2, 6, mtemp+1) +  ABtemp * (VY( 0, 6,mtemp)-CDcom*VY( 0, 6,mtemp+1));
    VY( 8, 7, mtemp) = Ptempy * VY( 2, 7, mtemp) + WPtempy * VY( 2, 7, mtemp+1) +  ABtemp * (VY( 0, 7,mtemp)-CDcom*VY( 0, 7,mtemp+1));
    VY( 8, 8, mtemp) = Ptempy * VY( 2, 8, mtemp) + WPtempy * VY( 2, 8, mtemp+1) +  ABtemp * (VY( 0, 8,mtemp)-CDcom*VY( 0, 8,mtemp+1)) + 2 * ABCDtemp * VY( 2, 2, mtemp+1);
    VY( 8, 9, mtemp) = Ptempy * VY( 2, 9, mtemp) + WPtempy * VY( 2, 9, mtemp+1) +  ABtemp * (VY( 0, 9,mtemp)-CDcom*VY( 0, 9,mtemp+1));
    
    VY( 9, 4, mtemp) = Ptempz * VY( 3, 4, mtemp) + WPtempz * VY( 3, 4, mtemp+1) +  ABtemp * (VY( 0, 4,mtemp)-CDcom*VY( 0, 4,mtemp+1));
    VY( 9, 5, mtemp) = Ptempz * VY( 3, 5, mtemp) + WPtempz * VY( 3, 5, mtemp+1) +  ABtemp * (VY( 0, 5,mtemp)-CDcom*VY( 0, 5,mtemp+1)) + ABCDtemp * VY( 3, 2, mtemp+1);
    VY( 9, 6, mtemp) = Ptempz * VY( 3, 6, mtemp) + WPtempz * VY( 3, 6, mtemp+1) +  ABtemp * (VY( 0, 6,mtemp)-CDcom*VY( 0, 6,mtemp+1)) + ABCDtemp * VY( 3, 1, mtemp+1);
    VY( 9, 7, mtemp) = Ptempz * VY( 3, 7, mtemp) + WPtempz * VY( 3, 7, mtemp+1) +  ABtemp * (VY( 0, 7,mtemp)-CDcom*VY( 0, 7,mtemp+1));
    VY( 9, 8, mtemp) = Ptempz * VY( 3, 8, mtemp) + WPtempz * VY( 3, 8, mtemp+1) +  ABtemp * (VY( 0, 8,mtemp)-CDcom*VY( 0, 8,mtemp+1));
    VY( 9, 9, mtemp) = Ptempz * VY( 3, 9, mtemp) + WPtempz * VY( 3, 9, mtemp+1) +  ABtemp * (VY( 0, 9,mtemp)-CDcom*VY( 0, 9,mtemp+1)) + 2 * ABCDtemp * VY( 3, 3, mtemp+1);
    
	return;
}


#ifndef TEST
__device__
#endif
QUICKDouble hrrwhole(int III, int JJJ, int KKK, int LLL, int IJKLTYPE, QUICKDouble* store, \
                     QUICKDouble RAx,QUICKDouble RAy,QUICKDouble RAz, \
                     QUICKDouble RBx,QUICKDouble RBy,QUICKDouble RBz, \
                     QUICKDouble RCx,QUICKDouble RCy,QUICKDouble RCz, \
                     QUICKDouble RDx,QUICKDouble RDy,QUICKDouble RDz)
{
    QUICKDouble Y;
    
    int NAx = LOC2(devSim.KLMN,0,III-1,3,devSim.nbasis);
    int NAy = LOC2(devSim.KLMN,1,III-1,3,devSim.nbasis);
    int NAz = LOC2(devSim.KLMN,2,III-1,3,devSim.nbasis);
    
    int NBx = LOC2(devSim.KLMN,0,JJJ-1,3,devSim.nbasis);
    int NBy = LOC2(devSim.KLMN,1,JJJ-1,3,devSim.nbasis);
    int NBz = LOC2(devSim.KLMN,2,JJJ-1,3,devSim.nbasis);
    
    int NCx = LOC2(devSim.KLMN,0,KKK-1,3,devSim.nbasis);
    int NCy = LOC2(devSim.KLMN,1,KKK-1,3,devSim.nbasis);
    int NCz = LOC2(devSim.KLMN,2,KKK-1,3,devSim.nbasis);
    
    int NDx = LOC2(devSim.KLMN,0,LLL-1,3,devSim.nbasis);
    int NDy = LOC2(devSim.KLMN,1,LLL-1,3,devSim.nbasis);
    int NDz = LOC2(devSim.KLMN,2,LLL-1,3,devSim.nbasis);
    
    
    int MA = LOC3(devTrans, NAx, NAy, NAz, TRANSDIM, TRANSDIM, TRANSDIM);
    int MB = LOC3(devTrans, NBx, NBy, NBz, TRANSDIM, TRANSDIM, TRANSDIM);
    int MC = LOC3(devTrans, NCx, NCy, NCz, TRANSDIM, TRANSDIM, TRANSDIM);
    int MD = LOC3(devTrans, NDx, NDy, NDz, TRANSDIM, TRANSDIM, TRANSDIM);

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
            if (NBx != 0) {
                Y = (QUICKDouble) LOC2(store, MB-1, 0, STOREDIM, STOREDIM) + (RAx-RBx)*LOC2(store, 0, 0, STOREDIM, STOREDIM);
            }else if (NBy != 0) {
                Y = (QUICKDouble) LOC2(store, MB-1, 0, STOREDIM, STOREDIM) + (RAy-RBy)*LOC2(store, 0, 0, STOREDIM, STOREDIM);
            }else if (NBz != 0) {
                Y = (QUICKDouble) LOC2(store, MB-1, 0, STOREDIM, STOREDIM) + (RAz-RBz)*LOC2(store, 0, 0, STOREDIM, STOREDIM);
            }
            break;
        }
        case 110:
        {
        
            if (NBx != 0) {
                Y = (QUICKDouble) LOC2(store, MB-1, MC-1, STOREDIM, STOREDIM) + (RAx-RBx)*LOC2(store, 0, MC-1, STOREDIM, STOREDIM);
            }else if (NBy != 0) {
                Y = (QUICKDouble) LOC2(store, MB-1, MC-1, STOREDIM, STOREDIM) + (RAy-RBy)*LOC2(store, 0, MC-1, STOREDIM, STOREDIM);
            }else if (NBz != 0) {
                Y = (QUICKDouble) LOC2(store, MB-1, MC-1, STOREDIM, STOREDIM) + (RAz-RBz)*LOC2(store, 0, MC-1, STOREDIM, STOREDIM);
            }
            break;
        }
        case 101:
        {
            QUICKDouble Y1,Y2;
            if (NDx != 0) {
                QUICKDouble c = (QUICKDouble) (RCx - RDx);
                Y1 = (QUICKDouble) LOC2(store, MB-1, MD-1 , STOREDIM, STOREDIM) + c * LOC2(store, MB-1,  0, STOREDIM, STOREDIM);
                Y2 = (QUICKDouble) LOC2(store,    0, MD-1 , STOREDIM, STOREDIM) + c * LOC2(store,    0,  0, STOREDIM, STOREDIM);
            }else if (NDy != 0) {
                QUICKDouble c = (QUICKDouble) (RCy - RDy);
                Y1 = (QUICKDouble) LOC2(store, MB-1, MD-1 , STOREDIM, STOREDIM) + c * LOC2(store, MB-1,  0, STOREDIM, STOREDIM);
                Y2 = (QUICKDouble) LOC2(store,    0, MD-1 , STOREDIM, STOREDIM) + c * LOC2(store,    0,  0, STOREDIM, STOREDIM);
            }else if (NDz != 0) {
                QUICKDouble c = (QUICKDouble) (RCz - RDz);
                Y1 = (QUICKDouble) LOC2(store, MB-1, MD-1 , STOREDIM, STOREDIM) + c * LOC2(store, MB-1,  0, STOREDIM, STOREDIM);
                Y2 = (QUICKDouble) LOC2(store,    0, MD-1 , STOREDIM, STOREDIM) + c * LOC2(store,    0,  0, STOREDIM, STOREDIM);
            }
            
            if (NBx != 0) {
                Y = Y1 + (RAx-RBx)*Y2;
            }else if (NBy != 0) {
                Y = Y1 + (RAy-RBy)*Y2;
            }else if (NBz != 0) {
                Y = Y1 + (RAz-RBz)*Y2;
            }
            break;
        }
        case 111:
        {
            QUICKDouble Y1,Y2;
            int MCD = (int) LOC3(devTrans, NCx+NDx, NCy+NDy, NCz+NDz, TRANSDIM, TRANSDIM, TRANSDIM);
            if (NDx != 0) {
                QUICKDouble c = (QUICKDouble) (RCx - RDx);
                Y1 = (QUICKDouble) LOC2(store, MB-1, MCD-1 , STOREDIM, STOREDIM) + c * LOC2(store, MB-1,  MC-1, STOREDIM, STOREDIM);
                Y2 = (QUICKDouble) LOC2(store,    0, MCD-1 , STOREDIM, STOREDIM) + c * LOC2(store,    0,  MC-1, STOREDIM, STOREDIM);
            }else if (NDy != 0) {
                QUICKDouble c = (QUICKDouble) (RCy - RDy);
                Y1 = (QUICKDouble) LOC2(store, MB-1, MCD-1 , STOREDIM, STOREDIM) + c * LOC2(store, MB-1,  MC-1, STOREDIM, STOREDIM);
                Y2 = (QUICKDouble) LOC2(store,    0, MCD-1 , STOREDIM, STOREDIM) + c * LOC2(store,    0,  MC-1, STOREDIM, STOREDIM);
            }else if (NDz != 0) {
                QUICKDouble c = (QUICKDouble) (RCz - RDz);
                Y1 = (QUICKDouble) LOC2(store, MB-1, MCD-1 , STOREDIM, STOREDIM) + c * LOC2(store, MB-1,  MC-1, STOREDIM, STOREDIM);
                Y2 = (QUICKDouble) LOC2(store,    0, MCD-1 , STOREDIM, STOREDIM) + c * LOC2(store,    0,  MC-1, STOREDIM, STOREDIM);
            }
            
            if (NBx != 0) {
                Y = Y1 + (RAx-RBx)*Y2;
            }else if (NBy != 0) {
                Y = Y1 + (RAy-RBy)*Y2;
            }else if (NBz != 0) {
                Y = Y1 + (RAz-RBz)*Y2;
            }
            break;
        }
        case 1100:
        {
            int MAB = (int) LOC3(devTrans, NAx+NBx, NAy+NBy, NAz+NBz, TRANSDIM, TRANSDIM, TRANSDIM);
            if (NBx != 0) {
                Y = (QUICKDouble) LOC2(store, MAB-1, 0 , STOREDIM, STOREDIM) + (RAx-RBx)*LOC2(store, MA-1, 0, STOREDIM, STOREDIM);
            }else if (NBy != 0) {
                Y = (QUICKDouble) LOC2(store, MAB-1, 0 , STOREDIM, STOREDIM) + (RAy-RBy)*LOC2(store, MA-1, 0, STOREDIM, STOREDIM);
            }else if (NBz != 0) {
                Y = (QUICKDouble) LOC2(store, MAB-1, 0 , STOREDIM, STOREDIM) + (RAz-RBz)*LOC2(store, MA-1, 0, STOREDIM, STOREDIM);
            }
            break;
        }
        case 1110:
        {   
            int MAB = (int) LOC3(devTrans, NAx+NBx, NAy+NBy, NAz+NBz, TRANSDIM, TRANSDIM, TRANSDIM);
            
            if (NBx != 0) {
                Y = (QUICKDouble) LOC2(store, MAB-1, MC-1 , STOREDIM, STOREDIM) + (RAx-RBx)*LOC2(store, MA-1, MC-1, STOREDIM, STOREDIM);
            }else if (NBy != 0) {
                Y = (QUICKDouble) LOC2(store, MAB-1, MC-1 , STOREDIM, STOREDIM) + (RAy-RBy)*LOC2(store, MA-1, MC-1, STOREDIM, STOREDIM);
            }else if (NBz != 0) {
                Y = (QUICKDouble) LOC2(store, MAB-1, MC-1 , STOREDIM, STOREDIM) + (RAz-RBz)*LOC2(store, MA-1, MC-1, STOREDIM, STOREDIM);
            }
            break;
        }
        case 1101:
        {
            QUICKDouble Y1,Y2;
            int MAB = (int) LOC3(devTrans, NAx+NBx, NAy+NBy, NAz+NBz, TRANSDIM, TRANSDIM, TRANSDIM);
            if (NDx != 0) {
                QUICKDouble c = (QUICKDouble) (RCx - RDx);
                Y1 = (QUICKDouble) LOC2(store, MAB-1, MD-1 , STOREDIM, STOREDIM) + c * LOC2(store, MAB-1,  0, STOREDIM, STOREDIM);
                Y2 = (QUICKDouble) LOC2(store,  MA-1, MD-1 , STOREDIM, STOREDIM) + c * LOC2(store,  MA-1,  0, STOREDIM, STOREDIM);
            }else if (NDy != 0) {
                QUICKDouble c = (QUICKDouble) (RCy - RDy);
                Y1 = (QUICKDouble) LOC2(store, MAB-1, MD-1 , STOREDIM, STOREDIM) + c * LOC2(store, MAB-1,  0, STOREDIM, STOREDIM);
                Y2 = (QUICKDouble) LOC2(store,  MA-1, MD-1 , STOREDIM, STOREDIM) + c * LOC2(store,  MA-1,  0, STOREDIM, STOREDIM);
            }else if (NDz != 0) {
                QUICKDouble c = (QUICKDouble) (RCz - RDz);
                Y1 = (QUICKDouble) LOC2(store, MAB-1, MD-1 , STOREDIM, STOREDIM) + c * LOC2(store, MAB-1,  0, STOREDIM, STOREDIM);
                Y2 = (QUICKDouble) LOC2(store,  MA-1, MD-1 , STOREDIM, STOREDIM) + c * LOC2(store,  MA-1,  0, STOREDIM, STOREDIM);
            }
            
            if (NBx != 0) {
                Y = Y1 + (RAx-RBx)*Y2;
            }else if (NBy != 0) {
                Y = Y1 + (RAy-RBy)*Y2;
            }else if (NBz != 0) {
                Y = Y1 + (RAz-RBz)*Y2;
            }
            break;
        }
        case 1111:
        {
            QUICKDouble Y1,Y2;
            int MAB = (int) LOC3(devTrans, NAx+NBx, NAy+NBy, NAz+NBz, TRANSDIM, TRANSDIM, TRANSDIM);
            int MCD = (int) LOC3(devTrans, NCx+NDx, NCy+NDy, NCz+NDz, TRANSDIM, TRANSDIM, TRANSDIM);

            if (NDx != 0) {
                QUICKDouble c = (QUICKDouble) (RCx - RDx);
                Y1 = (QUICKDouble) LOC2(store, MAB-1, MCD-1 , STOREDIM, STOREDIM) + c * LOC2(store, MAB-1, MC-1, STOREDIM, STOREDIM);
                Y2 = (QUICKDouble) LOC2(store,  MA-1, MCD-1 , STOREDIM, STOREDIM) + c * LOC2(store,  MA-1, MC-1, STOREDIM, STOREDIM);
            }else if (NDy != 0) {
                QUICKDouble c = (QUICKDouble) (RCy - RDy);
                Y1 = (QUICKDouble) LOC2(store, MAB-1, MCD-1 , STOREDIM, STOREDIM) + c * LOC2(store, MAB-1, MC-1, STOREDIM, STOREDIM);
                Y2 = (QUICKDouble) LOC2(store,  MA-1, MCD-1 , STOREDIM, STOREDIM) + c * LOC2(store,  MA-1, MC-1, STOREDIM, STOREDIM);
            }else if (NDz != 0) {
                QUICKDouble c = (QUICKDouble) (RCz - RDz);
                Y1 = (QUICKDouble) LOC2(store, MAB-1, MCD-1 , STOREDIM, STOREDIM) + c * LOC2(store, MAB-1, MC-1, STOREDIM, STOREDIM);
                Y2 = (QUICKDouble) LOC2(store,  MA-1, MCD-1 , STOREDIM, STOREDIM) + c * LOC2(store,  MA-1, MC-1, STOREDIM, STOREDIM);
            }
            
            if (NBx != 0) {
                Y = Y1 + (RAx-RBx)*Y2;
            }else if (NBy != 0) {
                Y = Y1 + (RAy-RBy)*Y2;
            }else if (NBz != 0) {
                Y = Y1 + (RAz-RBz)*Y2;
            }
            
            break;
        }
        case 1:
        {
            if (NDx != 0) {
                Y = (QUICKDouble) LOC2(store, 0, MD-1, STOREDIM, STOREDIM) + (RCx-RDx)*LOC2(store, 0, 0, STOREDIM, STOREDIM);
            }else if (NDy != 0) {
                Y = (QUICKDouble) LOC2(store, 0, MD-1, STOREDIM, STOREDIM) + (RCy-RDy)*LOC2(store, 0, 0, STOREDIM, STOREDIM);
            }else if (NDz != 0) {
                Y = (QUICKDouble) LOC2(store, 0, MD-1, STOREDIM, STOREDIM) + (RCz-RDz)*LOC2(store, 0, 0, STOREDIM, STOREDIM);
            }
            break;
        }
        case 11:
        {
            int MCD = (int) LOC3(devTrans, NCx+NDx, NCy+NDy, NCz+NDz, TRANSDIM, TRANSDIM, TRANSDIM);
            if (NDx != 0) {
                Y = (QUICKDouble) LOC2(store, 0, MCD-1, STOREDIM, STOREDIM) + (RCx-RDx)*LOC2(store, 0, MC-1, STOREDIM, STOREDIM);
            }else if (NDy != 0) {
                Y = (QUICKDouble) LOC2(store, 0, MCD-1, STOREDIM, STOREDIM) + (RCy-RDy)*LOC2(store, 0, MC-1, STOREDIM, STOREDIM);
            }else if (NDz != 0) {
                Y = (QUICKDouble) LOC2(store, 0, MCD-1, STOREDIM, STOREDIM) + (RCz-RDz)*LOC2(store, 0, MC-1, STOREDIM, STOREDIM);
            }
            break;
        }
        case 1001:
        {   
            if (NDx != 0) {
                Y = (QUICKDouble) LOC2(store, MA-1, MD-1, STOREDIM, STOREDIM) + (RCx-RDx)*LOC2(store, MA-1, 0, STOREDIM, STOREDIM);
            }else if (NDy != 0) {
                Y = (QUICKDouble) LOC2(store, MA-1, MD-1, STOREDIM, STOREDIM) + (RCy-RDy)*LOC2(store, MA-1, 0, STOREDIM, STOREDIM);
            }else if (NDz != 0) {
                Y = (QUICKDouble) LOC2(store, MA-1, MD-1, STOREDIM, STOREDIM) + (RCz-RDz)*LOC2(store, MA-1, 0, STOREDIM, STOREDIM);
            }
        }
        case 1011:
        {
            int MCD = (int) LOC3(devTrans, NCx+NDx, NCy+NDy, NCz+NDz, TRANSDIM, TRANSDIM, TRANSDIM);
            if (NDx != 0) {
                Y = (QUICKDouble) LOC2(store, MA-1, MCD-1, STOREDIM, STOREDIM) + (RCx-RDx)*LOC2(store, MA-1, MC-1, STOREDIM, STOREDIM);
            }else if (NDy != 0) {
                Y = (QUICKDouble) LOC2(store, MA-1, MCD-1, STOREDIM, STOREDIM) + (RCy-RDy)*LOC2(store, MA-1, MC-1, STOREDIM, STOREDIM);
            }else if (NDz != 0) {
                Y = (QUICKDouble) LOC2(store, MA-1, MCD-1, STOREDIM, STOREDIM) + (RCz-RDz)*LOC2(store, MA-1, MC-1, STOREDIM, STOREDIM);
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