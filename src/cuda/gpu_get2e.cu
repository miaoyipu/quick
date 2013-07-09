/*
 *  gpu_get2e.cpp
 *  new_quick
 *
 *  Created by Yipu Miao on 6/17/11.
 *  Copyright 2011 University of Florida.All rights reserved.
 *  
 *  Yipu Miao 9/15/11:  the first draft is released. And the GPUGP QM compuation can 
 *                      achieve as much as 15x faster at double precision level compared with CPU.
 */

#include "gpu.h"
#include <cuda.h>

#ifdef CUDA_SPDF
#include "int.h"
#endif

#include "int2.h"

/*
 Constant Memory in GPU is fast but quite limited and hard to operate, usually not allocatable and 
 readonly. So we put the following variables into constant memory:
 devSim: a gpu simluation type variable. which is to store to location of basic information about molecule and basis
 set. Note it only store the location, so it's mostly a set of pointer to GPU memory. and with some non-pointer
 value like the number of basis set. See gpu_type.h for details.
 devTrans : arrays to save the mapping index, will be elimited by hand writing unrolling code.
 Sumindex: a array to store refect how many temp variable needed in VRR. can be elimited by hand writing code.
 */
static __constant__ gpu_simulation_type devSim;
static __constant__ int devTrans[TRANSDIM*TRANSDIM*TRANSDIM];
static __constant__ int Sumindex[10]={0,0,1,4,10,20,35,56,84,120};


/*
 upload gpu simulation type to constant memory
 */
void upload_sim_to_constant(_gpu_type gpu){
    cudaError_t status;
	status = cudaMemcpyToSymbol(devSim, &gpu->gpu_sim, sizeof(gpu_simulation_type));
	PRINTERROR(status, " cudaMemcpyToSymbol, sim copy to constants failed")
}


#define int_spd
#undef int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#include "gpu_get2e_subs.h"

#ifdef CUDA_SPDF
#undef int_spd
#define int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#include "gpu_get2e_subs.h"

#undef int_spd
#undef int_spdf
#define int_spdf2
#undef int_spdf3
#undef int_spdf4
#include "gpu_get2e_subs.h"


#undef int_spd
#undef int_spdf
#undef int_spdf2
#define int_spdf3
#undef int_spdf4
#include "gpu_get2e_subs.h"


#undef int_spd
#undef int_spdf
#undef int_spdf2
#undef int_spdf3
#define int_spdf4
#include "gpu_get2e_subs.h"
#endif

#undef int_spd
#undef int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4


// totTime is the timer for GPU 2e time. Only on under debug mode
#ifdef DEBUG
static float totTime;
#endif

// =======   INTERFACE SECTION ===========================
// interface to call Kernel subroutine
void getAOInt(_gpu_type gpu, QUICKULL intStart, QUICKULL intEnd, cudaStream_t streamI, int streamID,  ERI_entry* aoint_buffer)
{
    QUICK_SAFE_CALL((getAOInt_kernel<<<gpu->blocks, gpu->twoEThreadsPerBlock, 0, streamI>>>(intStart, intEnd, aoint_buffer, streamID)));
#ifdef CUDA_SPDF
    // Part f-1
    QUICK_SAFE_CALL((getAOInt_kernel_spdf<<<gpu->blocks, gpu->twoEThreadsPerBlock, 0, streamI>>>( intStart, intEnd, aoint_buffer, streamID)));
    // Part f-2
    QUICK_SAFE_CALL((getAOInt_kernel_spdf2<<<gpu->blocks, gpu->twoEThreadsPerBlock, 0, streamI>>>( intStart, intEnd, aoint_buffer, streamID)));
    // Part f-3
    QUICK_SAFE_CALL((getAOInt_kernel_spdf3<<<gpu->blocks, gpu->twoEThreadsPerBlock, 0, streamI>>>( intStart, intEnd, aoint_buffer, streamID)));
    // Part f-4
    QUICK_SAFE_CALL((getAOInt_kernel_spdf4<<<gpu->blocks, gpu->twoEThreadsPerBlock, 0, streamI>>>( intStart, intEnd, aoint_buffer, streamID)));
#endif 
}

// interface to call Kernel subroutine
void get2e(_gpu_type gpu)
{
    // Part spd
    QUICK_SAFE_CALL((get2e_kernel<<<gpu->blocks, gpu->twoEThreadsPerBlock>>>()));
#ifdef CUDA_SPDF
    // Part f-1
    QUICK_SAFE_CALL((get2e_kernel_spdf<<<gpu->blocks, gpu->twoEThreadsPerBlock>>>()));
    // Part f-2
    QUICK_SAFE_CALL((get2e_kernel_spdf2<<<gpu->blocks, gpu->twoEThreadsPerBlock>>>()));
    // Part f-3
    QUICK_SAFE_CALL((get2e_kernel_spdf3<<<gpu->blocks, gpu->twoEThreadsPerBlock>>>()));
    // Part f-4
    QUICK_SAFE_CALL((get2e_kernel_spdf4<<<gpu->blocks, gpu->twoEThreadsPerBlock>>>()));
#endif 
}


// interface to call Kernel subroutine
void getAddInt(_gpu_type gpu, int bufferSize, ERI_entry* aoint_buffer)
{
    QUICK_SAFE_CALL((getAddInt_kernel<<<gpu->blocks, gpu->twoEThreadsPerBlock>>>(bufferSize, aoint_buffer)));
}



// =======   KERNEL SECTION ===========================
__global__ void getAddInt_kernel(int bufferSize, ERI_entry* aoint_buffer){
    unsigned int offside = blockIdx.x*blockDim.x+threadIdx.x;
    int totalThreads = blockDim.x*gridDim.x;
    int const batchSize = 20;
    ERI_entry a[batchSize];
    int j = 0;
    
    QUICKULL myInt = (QUICKULL) (bufferSize) / totalThreads;
    if ((bufferSize - myInt*totalThreads)> offside) myInt++;
    
    for (QUICKULL i = 1; i<=myInt; i++) {
        
        QUICKULL currentInt = totalThreads * (i-1) + offside;
        a[j] = aoint_buffer[currentInt];
        j++;
        if (j == batchSize || i == myInt) {
            
            for (int k = 0; k<j; k++) {
                int III = a[k].IJ / devSim.nbasis + 1;
                int JJJ = a[k].IJ % devSim.nbasis + 1;
                int KKK = a[k].KL / devSim.nbasis + 1;
                int LLL = a[k].KL % devSim.nbasis + 1;
                
                if (III <= devSim.nbasis && III >= 1 && JJJ <= devSim.nbasis && JJJ >= 1 && KKK <= devSim.nbasis && KKK >= 1 && LLL <= devSim.nbasis && LLL >= 1){
                    QUICKDouble hybrid_coeff = 0.0;
                    if (devSim.method == HF){
                        hybrid_coeff = 1.0;
                    }else if (devSim.method == B3LYP){
                        hybrid_coeff = 0.2;
                    }else if (devSim.method == DFT){
                        hybrid_coeff = 0.0;
                    }
                    
                    addint(devSim.oULL, a[k].value, III, JJJ, KKK, LLL, hybrid_coeff, devSim.dense, devSim.nbasis);
                }
            }
            j = 0;
        }
        
    }
    
}


__device__ __forceinline__ void addint(QUICKULL* oULL, QUICKDouble Y, int III, int JJJ, int KKK, int LLL,QUICKDouble hybrid_coeff,  QUICKDouble* dense, int nbasis){
    
    QUICKDouble DENSEKI = (QUICKDouble) LOC2(dense, KKK-1, III-1, nbasis, nbasis);
    QUICKDouble DENSEKJ = (QUICKDouble) LOC2(dense, KKK-1, JJJ-1, nbasis, nbasis);
    QUICKDouble DENSELJ = (QUICKDouble) LOC2(dense, LLL-1, JJJ-1, nbasis, nbasis);
    QUICKDouble DENSELI = (QUICKDouble) LOC2(dense, LLL-1, III-1, nbasis, nbasis);
    QUICKDouble DENSELK = (QUICKDouble) LOC2(dense, LLL-1, KKK-1, nbasis, nbasis);
    QUICKDouble DENSEJI = (QUICKDouble) LOC2(dense, JJJ-1, III-1, nbasis, nbasis);
    
    
    // ATOMIC ADD VALUE 1
    QUICKDouble _tmp = 2.0;
    if (KKK==LLL) {
        _tmp = 1.0;
    }
    
    QUICKDouble val1d = _tmp*DENSELK*Y;
    QUICKULL val1 = (QUICKULL) (fabs(val1d*OSCALE) + (QUICKDouble)0.5);
    if ( val1d < (QUICKDouble)0.0) val1 = 0ull - val1;
    QUICKADD(LOC2(oULL, JJJ-1, III-1, nbasis, nbasis), val1);
    
    
    // ATOMIC ADD VALUE 2
    if ((LLL != JJJ) || (III!=KKK)) {
        _tmp = 2.0;
        if (III==JJJ) {
            _tmp = 1.0;
        }
        
        QUICKDouble val2d = _tmp*DENSEJI*Y;
        QUICKULL val2 = (QUICKULL) (fabs(val2d*OSCALE) + (QUICKDouble)0.5);
        if ( val2d < (QUICKDouble)0.0) val2 = 0ull - val2;
        QUICKADD(LOC2(oULL, LLL-1, KKK-1, nbasis, nbasis), val2);
    }
    
    
    // ATOMIC ADD VALUE 3
    QUICKDouble val3d = hybrid_coeff*0.5*DENSELJ*Y;
    
    QUICKULL val3 = (QUICKULL) (fabs(val3d*OSCALE) + (QUICKDouble)0.5);
    if (((III == KKK) && (III <  JJJ) && (JJJ < LLL))) {
        val3 = (QUICKULL) (fabs(2*val3d*OSCALE) + (QUICKDouble)0.5);
    }
    if ( DENSELJ*Y < (QUICKDouble)0.0) val3 = 0ull - val3;
    QUICKADD(LOC2(oULL, KKK-1, III-1, nbasis, nbasis), 0ull-val3);
    
    // ATOMIC ADD VALUE 4
    if (KKK != LLL) {
        QUICKDouble val4d = hybrid_coeff*0.5*DENSEKJ*Y;
        
        QUICKULL val4 = (QUICKULL) (fabs(val4d*OSCALE) + (QUICKDouble)0.5);
        if ( val4d < (QUICKDouble)0.0) val4 = 0ull - val4;
        QUICKADD(LOC2(oULL, LLL-1, III-1, nbasis, nbasis), 0ull-val4);
    }
    
    
    
    // ATOMIC ADD VALUE 5
    QUICKDouble val5d = hybrid_coeff*0.5*DENSELI*Y;
    
    QUICKULL val5 = (QUICKULL) (fabs(val5d*OSCALE) + (QUICKDouble)0.5);
    if ( val5d < (QUICKDouble)0.0) val5 = 0ull - val5;
    
    if ((III != JJJ && III<KKK) || ((III == JJJ) && (III == KKK) && (III < LLL)) || ((III == KKK) && (III <  JJJ) && (JJJ < LLL))) {
        QUICKADD(LOC2(oULL, MAX(JJJ,KKK)-1, MIN(JJJ,KKK)-1, nbasis, nbasis), 0ull-val5);
    }
    
    
    // ATOMIC ADD VALUE 5 - 2
    if ( III != JJJ && JJJ == KKK) {
        QUICKADD(LOC2(oULL, JJJ-1, KKK-1, nbasis, nbasis), 0ull-val5);
    }
    
    // ATOMIC ADD VALUE 6
    if (III != JJJ) {
        if (KKK != LLL) {
            QUICKDouble val6d = hybrid_coeff*0.5*DENSEKI*Y;
            QUICKULL val6 = (QUICKULL) (fabs(val6d*OSCALE) + (QUICKDouble)0.5);
            if ( val6d < (QUICKDouble)0.0) val6 = 0ull - val6;
            
            QUICKADD(LOC2(oULL, MAX(JJJ,LLL)-1, MIN(JJJ,LLL)-1, devSim.nbasis, devSim.nbasis), 0ull-val6);
            
            // ATOMIC ADD VALUE 6 - 2
            if (JJJ == LLL && III!= KKK) {
                QUICKADD(LOC2(oULL, LLL-1, JJJ-1, nbasis, nbasis), 0ull-val6);
            }
        }
    }
}

__device__ __forceinline__ void FmT(int MaxM, QUICKDouble X, QUICKDouble* YVerticalTemp)
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

/*
 sqr for double precision. there no internal function to do that in fast-math-lib of CUDA
 */
__device__ __forceinline__ QUICKDouble quick_dsqr(QUICKDouble a)
{
    return a*a;
}


__device__ __forceinline__ QUICKDouble hrrwhole(int I, int J, int K, int L, \
                     int III, int JJJ, int KKK, int LLL, int IJKLTYPE, QUICKDouble* store, \
                     QUICKDouble RAx,QUICKDouble RAy,QUICKDouble RAz, \
                     QUICKDouble RBx,QUICKDouble RBy,QUICKDouble RBz, \
                     QUICKDouble RCx,QUICKDouble RCy,QUICKDouble RCz, \
                     QUICKDouble RDx,QUICKDouble RDy,QUICKDouble RDz)
{
    QUICKDouble Y;
#ifdef CUDA_SP
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
            break;
    }
#else
    
    int angularL[12], angularR[12];
    QUICKDouble coefAngularL[12], coefAngularR[12];
    Y = (QUICKDouble) 0.0;
    
    int numAngularL = lefthrr(RAx, RAy, RAz, RBx, RBy, RBz, 
                              LOC2(devSim.KLMN,0,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,III-1,3,devSim.nbasis),
                              LOC2(devSim.KLMN,0,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,JJJ-1,3,devSim.nbasis),
                              J, coefAngularL, angularL);
    int numAngularR = lefthrr(RCx, RCy, RCz, RDx, RDy, RDz,
                              LOC2(devSim.KLMN,0,KKK-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,KKK-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,KKK-1,3,devSim.nbasis),
                              LOC2(devSim.KLMN,0,LLL-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,LLL-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,LLL-1,3,devSim.nbasis),
                              L, coefAngularR, angularR);
    
    for (int i = 0; i<numAngularL; i++) {
        for (int j = 0; j<numAngularR; j++) {
            Y += coefAngularL[i] * coefAngularR[j] * LOC2(store, angularL[i]-1, angularR[j]-1 , STOREDIM, STOREDIM);
        }
    }
    //if (K == 2 && L == 3) Y = coefAngularL[0] * coefAngularR[2] * LOC2(store, angularL[0]-1, angularR[2]-1 , STOREDIM, STOREDIM);
    
    Y = Y * devSim.cons[III-1] * devSim.cons[JJJ-1] * devSim.cons[KKK-1] * devSim.cons[LLL-1];
#endif
    return Y;
}  


#ifndef CUDA_SP
__device__ __forceinline__ int lefthrr(QUICKDouble RAx, QUICKDouble RAy, QUICKDouble RAz, 
            QUICKDouble RBx, QUICKDouble RBy, QUICKDouble RBz,
            int KLMNAx, int KLMNAy, int KLMNAz,
            int KLMNBx, int KLMNBy, int KLMNBz,
            int IJTYPE,QUICKDouble* coefAngularL, int* angularL)
{           
    int numAngularL;
    switch (IJTYPE) {
        
        case 0:
        {
            numAngularL = 1;
            coefAngularL[0] = 1.0;
            angularL[0] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            break;
        }
        case 1:
        {
            coefAngularL[0] = 1.0;
            numAngularL = 2;
            angularL[0] = (int) LOC3(devTrans, KLMNAx + KLMNBx, KLMNAy + KLMNBy, KLMNAz + KLMNBz, TRANSDIM, TRANSDIM, TRANSDIM);
            angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            
            if (KLMNBx != 0) {
                coefAngularL[1] = RAx-RBx;
            }else if(KLMNBy !=0 ){
                coefAngularL[1] = RAy-RBy;
            }else if (KLMNBz != 0) {
                coefAngularL[1] = RAz-RBz;
            }
            break;
        }
        case 2:
        {
            coefAngularL[0] = 1.0;
            angularL[0] = (int) LOC3(devTrans, KLMNAx + KLMNBx, KLMNAy + KLMNBy, KLMNAz + KLMNBz, TRANSDIM, TRANSDIM, TRANSDIM);
            
            if (KLMNBx == 2) {
                numAngularL = 3;
                QUICKDouble tmp = RAx - RBx;
                coefAngularL[1] = 2 * tmp;
                angularL[1] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                coefAngularL[2]= tmp * tmp;
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if(KLMNBy == 2) {
                numAngularL = 3;
                QUICKDouble tmp = RAy - RBy;
                coefAngularL[1] = 2 * tmp;
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                coefAngularL[2]= tmp * tmp;
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBz == 2 ){
                numAngularL = 3;
                QUICKDouble tmp = RAz - RBz;
                coefAngularL[1] = 2 * tmp;
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                coefAngularL[2]= tmp * tmp;
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBx == 1 && KLMNBy == 1){
                numAngularL = 4;
                coefAngularL[1] = RAx - RBx;
                coefAngularL[2] = RAy - RBy;
                coefAngularL[3] = (RAx - RBx) * (RAy - RBy);
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                
            }else if (KLMNBx == 1 && KLMNBz == 1) {
                numAngularL = 4;
                coefAngularL[1] = RAx - RBx;
                coefAngularL[2] = RAz - RBz;
                coefAngularL[3] = (RAx - RBx) * (RAz - RBz);
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBy == 1 && KLMNBz == 1) {
                numAngularL = 4;
                coefAngularL[1] = RAy - RBy;
                coefAngularL[2] = RAz - RBz;
                coefAngularL[3] = (RAy - RBy) * (RAz - RBz);
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }
            break;
        }
        case 3:
        {
            coefAngularL[0] = 1.0;
            angularL[0] = (int) LOC3(devTrans, KLMNAx + KLMNBx, KLMNAy + KLMNBy, KLMNAz + KLMNBz, TRANSDIM, TRANSDIM, TRANSDIM);
            if (KLMNBx == 3) {
                numAngularL = 4;
                QUICKDouble tmp = RAx - RBx;
                
                coefAngularL[1] = 3 * tmp;
                coefAngularL[2] = 3 * tmp * tmp;
                coefAngularL[3] = tmp * tmp * tmp;
                
                angularL[1] = (int) LOC3(devTrans, KLMNAx+2, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBy == 3) {
                numAngularL = 4;
                QUICKDouble tmp = RAy - RBy;
                coefAngularL[1] = 3 * tmp;
                coefAngularL[2] = 3 * tmp * tmp;
                coefAngularL[3] = tmp * tmp * tmp;
                
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy+2, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBz == 3) {
                numAngularL = 4;
                
                QUICKDouble tmp = RAz - RBz;
                coefAngularL[1] = 3 * tmp;
                coefAngularL[2] = 3 * tmp * tmp;
                coefAngularL[3] = tmp * tmp * tmp;
                
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBx == 1 && KLMNBy ==2) { // case 120
                numAngularL = 6;
                QUICKDouble tmp = RAx - RBx;
                QUICKDouble tmp2 = RAy - RBy;
                
                coefAngularL[1] = tmp;
                coefAngularL[2] = 2 * tmp2;
                coefAngularL[3] = 2 * tmp * tmp2;
                coefAngularL[4] = tmp2 * tmp2;
                coefAngularL[5] = tmp * tmp2 * tmp2;
                
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy+2, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);   
                angularL[5] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBx == 1 && KLMNBz ==2) { // case 102
                numAngularL = 6;
                QUICKDouble tmp = RAx - RBx;
                QUICKDouble tmp2 = RAz - RBz;
                coefAngularL[1] = tmp;
                coefAngularL[2] = 2 * tmp2;
                coefAngularL[3] = 2 * tmp * tmp2;
                coefAngularL[4] = tmp2 * tmp2;
                coefAngularL[5] = tmp * tmp2 * tmp2;
                
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);   
                angularL[5] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBy == 1 && KLMNBz ==2) { // case 012
                numAngularL = 6;
                QUICKDouble tmp = RAy - RBy;
                QUICKDouble tmp2 = RAz - RBz;
                coefAngularL[1] = tmp;
                coefAngularL[2] = 2 * tmp2;
                coefAngularL[3] = 2 * tmp * tmp2;
                coefAngularL[4] = tmp2 * tmp2;
                coefAngularL[5] = tmp * tmp2 * tmp2;
                
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);   
                angularL[5] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBy== 2 && KLMNBz == 1) { // case 021
                numAngularL = 6;
                QUICKDouble tmp = RAz - RBz;
                QUICKDouble tmp2 = RAy - RBy;
                coefAngularL[1] = tmp;
                coefAngularL[2] = 2 * tmp2;
                coefAngularL[3] = 2 * tmp * tmp2;
                coefAngularL[4] = tmp2 * tmp2;
                coefAngularL[5] = tmp * tmp2 * tmp2;
                
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy+2, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBx == 2 && KLMNBy == 1) { // case 210
                numAngularL = 6;
                QUICKDouble tmp = RAy - RBy;
                QUICKDouble tmp2 = RAx - RBx;
                coefAngularL[1] = tmp;
                coefAngularL[2] = 2 * tmp2;
                coefAngularL[3] = 2 * tmp * tmp2;
                coefAngularL[4] = tmp2 * tmp2;
                coefAngularL[5] = tmp * tmp2 * tmp2;
                
                angularL[1] = (int) LOC3(devTrans, KLMNAx+2, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBx == 2 && KLMNBz ==1) { // case 201
                numAngularL = 6;
                QUICKDouble tmp = RAz - RBz;
                QUICKDouble tmp2 = RAx - RBx;
                coefAngularL[1] = tmp;
                coefAngularL[2] = 2 * tmp2;
                coefAngularL[3] = 2 * tmp * tmp2;
                coefAngularL[4] = tmp2 * tmp2;
                coefAngularL[5] = tmp * tmp2 * tmp2;
                
                angularL[1] = (int) LOC3(devTrans, KLMNAx+2, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBx == 1 && KLMNBy == 1) {
                numAngularL = 8;
                QUICKDouble tmp = RAx - RBx;
                QUICKDouble tmp2 = RAy - RBy;
                QUICKDouble tmp3 = RAz - RBz;
                
                coefAngularL[1] = tmp;
                coefAngularL[2] = tmp2;
                coefAngularL[3] = tmp3;
                coefAngularL[4] = tmp * tmp2;
                coefAngularL[5] = tmp * tmp3;                
                coefAngularL[6] = tmp2 * tmp3;
                coefAngularL[7] = tmp * tmp2 * tmp3;

                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[6] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[7] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }
            break;
            
        }
        case 4:
        {
            coefAngularL[0] = 1.0;
            angularL[0] = (int) LOC3(devTrans, KLMNAx + KLMNBx, KLMNAy + KLMNBy, KLMNAz + KLMNBz, TRANSDIM, TRANSDIM, TRANSDIM);
            if (KLMNBx == 4) {
                numAngularL = 5;
                QUICKDouble tmp = RAx - RBx;
                
                coefAngularL[1] = 4 * tmp;
                coefAngularL[2] = 6 * tmp * tmp;
                coefAngularL[3] = 4 * tmp * tmp * tmp;
                coefAngularL[4] = tmp * tmp * tmp * tmp;
                
                angularL[1] = (int) LOC3(devTrans, KLMNAx+3, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+2, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBy == 4) {
                numAngularL = 5;
                QUICKDouble tmp = RAy - RBy;
                coefAngularL[1] = 4 * tmp;
                coefAngularL[2] = 6 * tmp * tmp;
                coefAngularL[3] = 4 * tmp * tmp * tmp;
                coefAngularL[4] = tmp * tmp * tmp * tmp;
                
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy+3, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy+2, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                
            }else if (KLMNBz == 4) {
                numAngularL = 5;
                
                QUICKDouble tmp = RAz - RBz;
                coefAngularL[1] = 4 * tmp;
                coefAngularL[2] = 6 * tmp * tmp;
                coefAngularL[3] = 4 * tmp * tmp * tmp;
                coefAngularL[4] = tmp * tmp * tmp * tmp;
                
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+3, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                
            }else if (KLMNBx == 1 && KLMNBy ==3) {
                numAngularL = 8;
                QUICKDouble tmp = RAx - RBx;
                QUICKDouble tmp2 = RAy - RBy;
                
                coefAngularL[1] = tmp;
                coefAngularL[2] = 3 * tmp2;
                coefAngularL[3] = 3 * tmp * tmp2;
                coefAngularL[4] = 3 * tmp2 * tmp2;
                coefAngularL[5] = 3 * tmp * tmp2 * tmp2;
                coefAngularL[6] = tmp2 * tmp2 * tmp2;
                coefAngularL[7] = tmp * tmp2 * tmp2 * tmp2;
                
                angularL[1] = (int) LOC3(devTrans, KLMNAx,   KLMNAy+3, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy+2, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx,   KLMNAy+2, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int) LOC3(devTrans, KLMNAx,   KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[6] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[7] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBx == 1 && KLMNBz ==3) {
                numAngularL = 8;
                QUICKDouble tmp = RAx - RBx;
                QUICKDouble tmp2 = RAz - RBz;
                
                coefAngularL[1] = tmp;
                coefAngularL[2] = 3 * tmp2;
                coefAngularL[3] = 3 * tmp * tmp2;
                coefAngularL[4] = 3 * tmp2 * tmp2;
                coefAngularL[5] = 3 * tmp * tmp2 * tmp2;
                coefAngularL[6] = tmp2 * tmp2 * tmp2;
                coefAngularL[7] = tmp * tmp2 * tmp2 * tmp2;
                
                angularL[1] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz+3, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[6] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[7] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                
            }else if (KLMNBy == 1 && KLMNBz ==3) {
                numAngularL = 8;
                QUICKDouble tmp = RAy - RBy;
                QUICKDouble tmp2 = RAz - RBz;
                
                coefAngularL[1] = tmp;
                coefAngularL[2] = 3 * tmp2;
                coefAngularL[3] = 3 * tmp * tmp2;
                coefAngularL[4] = 3 * tmp2 * tmp2;
                coefAngularL[5] = 3 * tmp * tmp2 * tmp2;
                coefAngularL[6] = tmp2 * tmp2 * tmp2;
                coefAngularL[7] = tmp * tmp2 * tmp2 * tmp2;
                
                angularL[1] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz+3, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[6] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[7] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                
            }else if (KLMNBx == 2 && KLMNBy == 2) {
                numAngularL = 9;
                QUICKDouble tmp = RAx - RBx;
                QUICKDouble tmp2 = RAy - RBy;
                
                coefAngularL[1] = 2 * tmp;
                coefAngularL[2] = 2 * tmp2;
                coefAngularL[3] = 4 * tmp * tmp2;
                coefAngularL[4] = tmp * tmp;
                coefAngularL[5] = tmp2 * tmp2;
                coefAngularL[6] = 2 * tmp * tmp2 * tmp2;
                coefAngularL[7] = 2 * tmp * tmp * tmp2;
                coefAngularL[8] = tmp * tmp * tmp2 * tmp2;
                
                angularL[1] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy+2, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+2, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx,   KLMNAy+2, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int) LOC3(devTrans, KLMNAx+2, KLMNAy,   KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[6] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy,   KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[7] = (int) LOC3(devTrans, KLMNAx,   KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[8] = (int) LOC3(devTrans, KLMNAx,   KLMNAy,   KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBx == 2 && KLMNBz == 2) {
                numAngularL = 9;
                QUICKDouble tmp = RAx - RBx;
                QUICKDouble tmp2 = RAz - RBz;
                
                coefAngularL[1] = 2 * tmp;
                coefAngularL[2] = 2 * tmp2;
                coefAngularL[3] = 4 * tmp * tmp2;
                coefAngularL[4] = tmp * tmp;
                coefAngularL[5] = tmp2 * tmp2;
                coefAngularL[6] = 2 * tmp * tmp2 * tmp2;
                coefAngularL[7] = 2 * tmp * tmp * tmp2;
                coefAngularL[8] = tmp * tmp * tmp2 * tmp2;
                
                angularL[1] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+2, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int) LOC3(devTrans, KLMNAx+2, KLMNAy,   KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[6] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy,   KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[7] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[8] = (int) LOC3(devTrans, KLMNAx,   KLMNAy,   KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBy == 2 && KLMNBz == 2) {
                numAngularL = 9;
                QUICKDouble tmp = RAy - RBy;
                QUICKDouble tmp2 = RAz - RBz;
                
                coefAngularL[1] = 2 * tmp;
                coefAngularL[2] = 2 * tmp2;
                coefAngularL[3] = 4 * tmp * tmp2;
                coefAngularL[4] = tmp * tmp;
                coefAngularL[5] = tmp2 * tmp2;
                coefAngularL[6] = 2 * tmp * tmp2 * tmp2;
                coefAngularL[7] = 2 * tmp * tmp * tmp2;
                coefAngularL[8] = tmp * tmp * tmp2 * tmp2;
                
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy+2, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int) LOC3(devTrans, KLMNAx, KLMNAy+2,   KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[6] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1,   KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[7] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[8] = (int) LOC3(devTrans, KLMNAx,   KLMNAy,   KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBx == 1 && KLMNBy == 1) {
                numAngularL = 12;
                QUICKDouble tmp = RAx - RBx;
                QUICKDouble tmp2 = RAy - RBy;
                QUICKDouble tmp3 = RAz - RBz;
                
                coefAngularL[1] = tmp;
                coefAngularL[2] = tmp2;
                coefAngularL[3] = 2 * tmp3;
                coefAngularL[4] = tmp * tmp2;
                coefAngularL[5] = 2 * tmp * tmp3;
                coefAngularL[6] = 2 * tmp2 * tmp3;
                coefAngularL[7] = tmp3 * tmp3;
                coefAngularL[8] = 2 * tmp * tmp2 * tmp3;
                coefAngularL[9] = tmp * tmp3 * tmp3;
                coefAngularL[10] = tmp * tmp3 * tmp3;
                coefAngularL[11] = tmp * tmp2 * tmp3 * tmp3;
                
                angularL[1] = (int)  LOC3(devTrans, KLMNAx,   KLMNAy+1, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int)  LOC3(devTrans, KLMNAx+1, KLMNAy,   KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int)  LOC3(devTrans, KLMNAx+1, KLMNAy+1, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int)  LOC3(devTrans, KLMNAx,   KLMNAy,   KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int)  LOC3(devTrans, KLMNAx,   KLMNAy+1, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[6] = (int)  LOC3(devTrans, KLMNAx+1, KLMNAy,   KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[7] = (int)  LOC3(devTrans, KLMNAx+1, KLMNAy+1, KLMNAz,   TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[8] = (int)  LOC3(devTrans, KLMNAx,   KLMNAy,   KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[9] = (int)  LOC3(devTrans, KLMNAx,   KLMNAy+1, KLMNAz,   TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[10] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy,   KLMNAz,   TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[11] = (int) LOC3(devTrans, KLMNAx,   KLMNAy,   KLMNAz,   TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBx == 1 && KLMNBz == 1) {
                numAngularL = 12;
                QUICKDouble tmp = RAx - RBx;
                QUICKDouble tmp2 = RAz - RBz;
                QUICKDouble tmp3 = RAy - RBy;
                
                coefAngularL[1] = tmp;
                coefAngularL[2] = tmp2;
                coefAngularL[3] = 2 * tmp3;
                coefAngularL[4] = tmp * tmp2;
                coefAngularL[5] = 2 * tmp * tmp3;
                coefAngularL[6] = 2 * tmp2 * tmp3;
                coefAngularL[7] = tmp3 * tmp3;
                coefAngularL[8] = 2 * tmp * tmp2 * tmp3;
                coefAngularL[9] = tmp * tmp3 * tmp3;
                coefAngularL[10] = tmp * tmp3 * tmp3;
                coefAngularL[11] = tmp * tmp2 * tmp3 * tmp3;
                
                angularL[1] = (int)  LOC3(devTrans, KLMNAx,   KLMNAy+2, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int)  LOC3(devTrans, KLMNAx+1, KLMNAy+2, KLMNAz,   TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int)  LOC3(devTrans, KLMNAx+1, KLMNAy+1, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int)  LOC3(devTrans, KLMNAx,   KLMNAy+2, KLMNAz,   TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int)  LOC3(devTrans, KLMNAx,   KLMNAy+1, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[6] = (int)  LOC3(devTrans, KLMNAx+1, KLMNAy+1, KLMNAz,   TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[7] = (int)  LOC3(devTrans, KLMNAx+1, KLMNAy,   KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[8] = (int)  LOC3(devTrans, KLMNAx,   KLMNAy+1, KLMNAz,   TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[9] = (int)  LOC3(devTrans, KLMNAx,   KLMNAy,   KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[10] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy,   KLMNAz,   TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[11] = (int) LOC3(devTrans, KLMNAx,   KLMNAy,   KLMNAz,   TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBy == 1 && KLMNBz == 1) {
                numAngularL = 12;
                QUICKDouble tmp = RAy - RBy;
                QUICKDouble tmp2 = RAz - RBz;
                QUICKDouble tmp3 = RAx - RBx;
                
                coefAngularL[1] = tmp;
                coefAngularL[2] = tmp2;
                coefAngularL[3] = 2 * tmp3;
                coefAngularL[4] = tmp * tmp2;
                coefAngularL[5] = 2 * tmp * tmp3;
                coefAngularL[6] = 2 * tmp2 * tmp3;
                coefAngularL[7] = tmp3 * tmp3;
                coefAngularL[8] = 2 * tmp * tmp2 * tmp3;
                coefAngularL[9] = tmp * tmp3 * tmp3;
                coefAngularL[10] = tmp * tmp3 * tmp3;
                coefAngularL[11] = tmp * tmp2 * tmp3 * tmp3;
                
                angularL[1] = (int)  LOC3(devTrans, KLMNAx+2, KLMNAy,   KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int)  LOC3(devTrans, KLMNAx+2, KLMNAy+1, KLMNAz,   TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int)  LOC3(devTrans, KLMNAx+1, KLMNAy+1, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int)  LOC3(devTrans, KLMNAx+2, KLMNAy,   KLMNAz,   TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int)  LOC3(devTrans, KLMNAx+1, KLMNAy,   KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[6] = (int)  LOC3(devTrans, KLMNAx+1, KLMNAy+1, KLMNAz,   TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[7] = (int)  LOC3(devTrans, KLMNAx,   KLMNAy+1, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[8] = (int)  LOC3(devTrans, KLMNAx+1, KLMNAy,   KLMNAz,   TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[9] = (int)  LOC3(devTrans, KLMNAx,   KLMNAy+1, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[10] = (int) LOC3(devTrans, KLMNAx,   KLMNAy,   KLMNAz,   TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[11] = (int) LOC3(devTrans, KLMNAx,   KLMNAy,   KLMNAz,   TRANSDIM, TRANSDIM, TRANSDIM);
            }
            
            
            break;

            
        }
    }
    return numAngularL;
}

#endif

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
    // upload to trans device location
    cudaError_t status;

    status = cudaMemcpyToSymbol(devTrans, trans, sizeof(int)*TRANSDIM*TRANSDIM*TRANSDIM);
    PRINTERROR(status, " cudaMemcpyToSymbol, Trans copy to constants failed")

}

