//
//  gpu_get2e_subs.h
//  new_quick 2
//
//  Created by Yipu Miao on 6/18/13.
//
//
#include "gpu_common.h"


/*
 In the following kernel, we treat f orbital into 5 parts.
 
 type:   ss sp ps sd ds pp dd sf pf | df ff |
 ss                                 |       |
 sp                                 |       |
 ps                                 | zone  |
 sd                                 |  2    |
 ds         zone 0                  |       |
 pp                                 |       |
 dd                                 |       |
 sf                                 |       |
 pf                                 |       |
 -------------------------------------------
 df         zone 1                  | z | z |
 ff                                 | 3 | 4 |
 -------------------------------------------
 
 
 because the single f orbital kernel is impossible to compile completely, we treat VRR as:
 
 
 I+J  0 1 2 3 4 | 5 | 6 |
 0 ----------------------
 1|             |       |
 2|   Kernel    |  K2   |
 3|     0       |       |
 4|             |       |
 -----------------------|
 5|   Kernel    | K | K |
 6|     1       | 3 | 4 |
 ------------------------
 
 Their responses for
              I+J          K+L
 Kernel 0:   0-4           0-4
 Kernel 1:   0-4           5,6
 Kernel 2:   5,6           0-4
 Kernel 3:   5             5,6
 Kernel 4:   6             5,6
 
 Integrals in zone need kernel:
 zone 0: kernel 0
 zone 1: kernel 0,1
 zone 2: kernel 0,2
 zone 3: kernel 0,1,2,3
 zone 4: kernel 0,1,2,3,4
 
 so first, kernel 0: zone 0,1,2,3,4 (get2e_kernel()), if no f, then that's it.
 second,   kernel 1: zone 1,3,4(get2e_kernel_spdf())
 then,     kernel 2: zone 2,3,4(get2e_kernel_spdf2())
 then,     kernel 3: zone 3,4(get2e_kernel_spdf3())
 finally,  kernel 4: zone 4(get2e_kernel_spdf4())
 
 */
#ifdef int_spd
__global__ void get2e_kernel()
#elif defined int_spdf
__global__ void get2e_kernel_spdf()
#elif defined int_spdf2
__global__ void get2e_kernel_spdf2()
#elif defined int_spdf3
__global__ void get2e_kernel_spdf3()
#elif defined int_spdf4
__global__ void get2e_kernel_spdf4()
#endif
{
    unsigned int offside = blockIdx.x*blockDim.x+threadIdx.x;
    int totalThreads = blockDim.x*gridDim.x;
    
    
    
#ifdef int_spd
    QUICKULL jshell = (QUICKULL) devSim.sqrQshell;
    QUICKULL jshell2 = (QUICKULL) devSim.sqrQshell;
    
#elif defined int_spdf
    
    QUICKULL jshell = (QUICKULL) devSim.sqrQshell;
    QUICKULL jshell2 = (QUICKULL) devSim.sqrQshell - devSim.fStart;
    
#elif defined int_spdf2
    
    QUICKULL jshell = (QUICKULL) devSim.sqrQshell;
    QUICKULL jshell2 = (QUICKULL) devSim.sqrQshell - devSim.fStart;
    
#elif defined int_spdf3
    
    QUICKULL jshell0 = (QUICKULL) devSim.fStart;
    QUICKULL jshell = (QUICKULL) devSim.sqrQshell - jshell0;
    QUICKULL jshell2 = (QUICKULL) devSim.sqrQshell - jshell0;
    
#elif defined int_spdf4
    
    QUICKULL jshell0 = (QUICKULL) devSim.fStart;
    QUICKULL jshell00 = (QUICKULL) devSim.ffStart;
    QUICKULL jshell = (QUICKULL) devSim.sqrQshell - jshell00;
    QUICKULL jshell2 = (QUICKULL) devSim.sqrQshell - jshell0;
    
#endif
    
    for (QUICKULL i = offside; i<jshell2*jshell; i+= totalThreads) {
        
#ifdef int_spd
        
        // Zone 0
        QUICKULL a = (QUICKULL) i/jshell;
        QUICKULL b = (QUICKULL) (i - a*jshell);
        
#elif defined int_spdf
        
        
        // Zone 1
        QUICKULL b = (QUICKULL) i/jshell;
        QUICKULL a = (QUICKULL) (i - b*jshell);
        b = b + devSim.fStart;
        
#elif defined int_spdf2
        
        // Zone 2
        QUICKULL a = (QUICKULL) i/jshell;
        QUICKULL b = (QUICKULL) (i - a*jshell);
        a = a + devSim.fStart;
        
#elif defined int_spdf3
        
        // Zone 3
        QUICKULL a, b;
        if (jshell != 0 ) {
            a = (QUICKULL) i/jshell;
            b = (QUICKULL) (i - a*jshell);
            a = a + jshell0;
            b = b + jshell0;
        }else{
            a = 0;
            b = 0;
        }
#elif defined int_spdf4
        
        // Zone 4
        QUICKULL a, b;
        if (jshell2 != 0 ) {
            a = (QUICKULL) i/jshell2;
            b = (QUICKULL) (i - a*jshell2);
            a = a + jshell00;
            b = b + jshell0;
        }else{
            a = 0;
            b = 0;
        }
#endif
        
        int II = devSim.sorted_YCutoffIJ[a].x;
        int KK = devSim.sorted_YCutoffIJ[b].x;
        
        int ii = devSim.sorted_Q[II];
        int kk = devSim.sorted_Q[KK];
        
        if (ii<=kk){
            
            int JJ = devSim.sorted_YCutoffIJ[a].y;
            int LL = devSim.sorted_YCutoffIJ[b].y;
            
            int jj = devSim.sorted_Q[JJ];
            int ll = devSim.sorted_Q[LL];
            
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
#ifdef int_spd
                
                iclass(iii, jjj, kkk, lll, ii, jj, kk, ll, DNMax);
                
#elif defined int_spdf
                iclass_spdf(iii, jjj, kkk, lll, ii, jj, kk, ll, DNMax);
                
#elif defined int_spdf2
                
                iclass_spdf2(iii, jjj, kkk, lll, ii, jj, kk, ll, DNMax);
                
#elif defined int_spdf3
                
                iclass_spdf3(iii, jjj, kkk, lll, ii, jj, kk, ll, DNMax);
                
#elif defined int_spdf4
                
                iclass_spdf4(iii, jjj, kkk, lll, ii, jj, kk, ll, DNMax);
#endif
                
            }
        }
    }
}

/*
 iclass subroutine is to generate 2-electron intergral using HRR and VRR method, which is the most
 performance algrithem for electron intergral evaluation. See description below for details
 */
#ifdef int_spd
__device__ __forceinline__ void iclass
#elif defined int_spdf
__device__ __forceinline__ void iclass_spdf
#elif defined int_spdf2
__device__ __forceinline__ void iclass_spdf2
#elif defined int_spdf3
__device__ __forceinline__ void iclass_spdf3
#elif defined int_spdf4
__device__ __forceinline__ void iclass_spdf4
#endif
                                      (int I, int J, int K, int L, unsigned int II, unsigned int JJ, unsigned int KK, unsigned int LL, QUICKDouble DNMax)
{
    
    /*
     kAtom A, B, C ,D is the coresponding atom for shell ii, jj, kk, ll
     and be careful with the index difference between Fortran and C++,
     Fortran starts array index with 1 and C++ starts 0.
     
     
     RA, RB, RC, and RD are the coordinates for atom katomA, katomB, katomC and katomD,
     which means they are corrosponding coorinates for shell II, JJ, KK, and LL.
     And we don't need the coordinates now, so we will not retrieve the data now.
     */
    QUICKDouble RAx = LOC2(devSim.xyz, 0 , devSim.katom[II]-1, 3, devSim.natom);
    QUICKDouble RAy = LOC2(devSim.xyz, 1 , devSim.katom[II]-1, 3, devSim.natom);
    QUICKDouble RAz = LOC2(devSim.xyz, 2 , devSim.katom[II]-1, 3, devSim.natom);
    
    QUICKDouble RCx = LOC2(devSim.xyz, 0 , devSim.katom[KK]-1, 3, devSim.natom);
    QUICKDouble RCy = LOC2(devSim.xyz, 1 , devSim.katom[KK]-1, 3, devSim.natom);
    QUICKDouble RCz = LOC2(devSim.xyz, 2 , devSim.katom[KK]-1, 3, devSim.natom);
    
    /*
     kPrimI, J, K and L indicates the primtive gaussian function number
     kStartI, J, K, and L indicates the starting guassian function for shell I, J, K, and L.
     We retrieve from global memory and save them to register to avoid multiple retrieve.
     */
    int kPrimI = devSim.kprim[II];
    int kPrimJ = devSim.kprim[JJ];
    int kPrimK = devSim.kprim[KK];
    int kPrimL = devSim.kprim[LL];
    
    int kStartI = devSim.kstart[II]-1;
    int kStartJ = devSim.kstart[JJ]-1;
    int kStartK = devSim.kstart[KK]-1;
    int kStartL = devSim.kstart[LL]-1;
    
    
    /*
     store saves temp contracted integral as [as|bs] type. the dimension should be allocatable but because
     of cuda limitation, we can not do that now.
     
     See M.Head-Gordon and J.A.Pople, Jchem.Phys., 89, No.9 (1988) for VRR algrithem details.
     */
    QUICKDouble store[STOREDIM*STOREDIM];
    
    /*
     Initial the neccessary element for
     */
    for (int i = Sumindex[K+1]+1; i<= Sumindex[K+L+2]; i++) {
        for (int j = Sumindex[I+1]+1; j<= Sumindex[I+J+2]; j++) {
            LOC2(store, j-1, i-1, STOREDIM, STOREDIM) = 0;
        }
    }
    
    for (int i = 0; i<kPrimI*kPrimJ;i++){
        int JJJ = (int) i/kPrimI;
        int III = (int) i-kPrimI*JJJ;
        /*
         In the following comments, we have I, J, K, L denote the primitive gaussian function we use, and
         for example, expo(III, ksumtype(II)) stands for the expo for the IIIth primitive guassian function for II shell,
         we use I to express the corresponding index.
         AB = expo(I)+expo(J)
         --->                --->
         ->     expo(I) * xyz (I) + expo(J) * xyz(J)
         P  = ---------------------------------------
         expo(I) + expo(J)
         Those two are pre-calculated in CPU stage.
         
         */
        int ii_start = devSim.prim_start[II];
        int jj_start = devSim.prim_start[JJ];
        
        QUICKDouble AB = LOC2(devSim.expoSum, ii_start+III, jj_start+JJJ, devSim.prim_total, devSim.prim_total);
        QUICKDouble Px = LOC2(devSim.weightedCenterX, ii_start+III, jj_start+JJJ, devSim.prim_total, devSim.prim_total);
        QUICKDouble Py = LOC2(devSim.weightedCenterY, ii_start+III, jj_start+JJJ, devSim.prim_total, devSim.prim_total);
        QUICKDouble Pz = LOC2(devSim.weightedCenterZ, ii_start+III, jj_start+JJJ, devSim.prim_total, devSim.prim_total);
        
        /*
         X1 is the contracted coeffecient, which is pre-calcuated in CPU stage as well.
         cutoffprim is used to cut too small prim gaussian function when bring density matrix into consideration.
         */
        QUICKDouble cutoffPrim = DNMax * LOC2(devSim.cutPrim, kStartI+III, kStartJ+JJJ, devSim.jbasis, devSim.jbasis);
        QUICKDouble X1 = LOC4(devSim.Xcoeff, kStartI+III, kStartJ+JJJ, I - devSim.Qstart[II], J - devSim.Qstart[JJ], devSim.jbasis, devSim.jbasis, 2, 2);
        
        for (int j = 0; j<kPrimK*kPrimL; j++){
            int LLL = (int)j/kPrimK;
            int KKK = (int) j-kPrimK*LLL;
            
            if (cutoffPrim * LOC2(devSim.cutPrim, kStartK+KKK, kStartL+LLL, devSim.jbasis, devSim.jbasis) > devSim.primLimit) {
                /*
                 CD = expo(L)+expo(K)
                 ABCD = 1/ (AB + CD) = 1 / (expo(I)+expo(J)+expo(K)+expo(L))
                 AB * CD      (expo(I)+expo(J))*(expo(K)+expo(L))
                 Rou(Greek Letter) =   ----------- = ------------------------------------
                 AB + CD         expo(I)+expo(J)+expo(K)+expo(L)
                 
                 expo(I)+expo(J)                        expo(K)+expo(L)
                 ABcom = --------------------------------  CDcom = --------------------------------
                 expo(I)+expo(J)+expo(K)+expo(L)           expo(I)+expo(J)+expo(K)+expo(L)
                 
                 ABCDtemp = 1/2(expo(I)+expo(J)+expo(K)+expo(L))
                 */
                
                int kk_start = devSim.prim_start[KK];
                int ll_start = devSim.prim_start[LL];
                
                QUICKDouble CD = LOC2(devSim.expoSum, kk_start+KKK, ll_start+LLL, devSim.prim_total, devSim.prim_total);
                
                QUICKDouble ABCD = 1/(AB+CD);
                
                /*
                 X2 is the multiplication of four indices normalized coeffecient
                 */
                QUICKDouble X2 = sqrt(ABCD) * X1 * LOC4(devSim.Xcoeff, kStartK+KKK, kStartL+LLL, K - devSim.Qstart[KK], L - devSim.Qstart[LL], devSim.jbasis, devSim.jbasis, 2, 2);
                
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
                 
                 ->  -> 2
                 T = ROU * | P - Q|
                 */
                
                QUICKDouble Qx = LOC2(devSim.weightedCenterX, kk_start+KKK, ll_start+LLL, devSim.prim_total, devSim.prim_total);
                QUICKDouble Qy = LOC2(devSim.weightedCenterY, kk_start+KKK, ll_start+LLL, devSim.prim_total, devSim.prim_total);
                QUICKDouble Qz = LOC2(devSim.weightedCenterZ, kk_start+KKK, ll_start+LLL, devSim.prim_total, devSim.prim_total);
                
                //QUICKDouble T = AB * CD * ABCD * ( quick_dsqr(Px-Qx) + quick_dsqr(Py-Qy) + quick_dsqr(Pz-Qz));
                
                QUICKDouble YVerticalTemp[VDIM1*VDIM2*VDIM3];
                FmT(I+J+K+L, AB * CD * ABCD * ( quick_dsqr(Px-Qx) + quick_dsqr(Py-Qy) + quick_dsqr(Pz-Qz)), YVerticalTemp);
                for (int i = 0; i<=I+J+K+L; i++) {
                    VY(0, 0, i) = VY(0, 0, i) * X2;
                }
#ifdef int_spd
                vertical(I, J, K, L, YVerticalTemp, store, \
                         Px - RAx, Py - RAy, Pz - RAz, (Px*AB+Qx*CD)*ABCD - Px, (Py*AB+Qy*CD)*ABCD - Py, (Pz*AB+Qz*CD)*ABCD - Pz, \
                         Qx - RCx, Qy - RCy, Qz - RCz, (Px*AB+Qx*CD)*ABCD - Qx, (Py*AB+Qy*CD)*ABCD - Qy, (Pz*AB+Qz*CD)*ABCD - Qz, \
                         0.5 * ABCD, 0.5 / AB, 0.5 / CD, AB * ABCD, CD * ABCD);
#elif defined int_spdf
                
                vertical_spdf(I, J, K, L, YVerticalTemp, store, \
                         Px - RAx, Py - RAy, Pz - RAz, (Px*AB+Qx*CD)*ABCD - Px, (Py*AB+Qy*CD)*ABCD - Py, (Pz*AB+Qz*CD)*ABCD - Pz, \
                         Qx - RCx, Qy - RCy, Qz - RCz, (Px*AB+Qx*CD)*ABCD - Qx, (Py*AB+Qy*CD)*ABCD - Qy, (Pz*AB+Qz*CD)*ABCD - Qz, \
                         0.5 * ABCD, 0.5 / AB, 0.5 / CD, AB * ABCD, CD * ABCD);
#elif defined int_spdf2
                
                vertical_spdf2(I, J, K, L, YVerticalTemp, store, \
                              Px - RAx, Py - RAy, Pz - RAz, (Px*AB+Qx*CD)*ABCD - Px, (Py*AB+Qy*CD)*ABCD - Py, (Pz*AB+Qz*CD)*ABCD - Pz, \
                              Qx - RCx, Qy - RCy, Qz - RCz, (Px*AB+Qx*CD)*ABCD - Qx, (Py*AB+Qy*CD)*ABCD - Qy, (Pz*AB+Qz*CD)*ABCD - Qz, \
                              0.5 * ABCD, 0.5 / AB, 0.5 / CD, AB * ABCD, CD * ABCD);
#elif defined int_spdf3
                
                vertical_spdf3(I, J, K, L, YVerticalTemp, store, \
                              Px - RAx, Py - RAy, Pz - RAz, (Px*AB+Qx*CD)*ABCD - Px, (Py*AB+Qy*CD)*ABCD - Py, (Pz*AB+Qz*CD)*ABCD - Pz, \
                              Qx - RCx, Qy - RCy, Qz - RCz, (Px*AB+Qx*CD)*ABCD - Qx, (Py*AB+Qy*CD)*ABCD - Qy, (Pz*AB+Qz*CD)*ABCD - Qz, \
                              0.5 * ABCD, 0.5 / AB, 0.5 / CD, AB * ABCD, CD * ABCD);
#elif defined int_spdf4
                
                vertical_spdf4(I, J, K, L, YVerticalTemp, store, \
                               Px - RAx, Py - RAy, Pz - RAz, (Px*AB+Qx*CD)*ABCD - Px, (Py*AB+Qy*CD)*ABCD - Py, (Pz*AB+Qz*CD)*ABCD - Pz, \
                               Qx - RCx, Qy - RCy, Qz - RCz, (Px*AB+Qx*CD)*ABCD - Qx, (Py*AB+Qy*CD)*ABCD - Qy, (Pz*AB+Qz*CD)*ABCD - Qz, \
                               0.5 * ABCD, 0.5 / AB, 0.5 / CD, AB * ABCD, CD * ABCD);
#endif
                
            }
        }
    }
    
    
    // IJKLTYPE is the I, J, K,L type
    int IJKLTYPE = (int) (1000 * I + 100 *J + 10 * K + L);
    
    QUICKDouble RBx, RBy, RBz;
    QUICKDouble RDx, RDy, RDz;
    
    RBx = LOC2(devSim.xyz, 0 , devSim.katom[JJ]-1, 3, devSim.natom);
    RBy = LOC2(devSim.xyz, 1 , devSim.katom[JJ]-1, 3, devSim.natom);
    RBz = LOC2(devSim.xyz, 2 , devSim.katom[JJ]-1, 3, devSim.natom);
    
    
    RDx = LOC2(devSim.xyz, 0 , devSim.katom[LL]-1, 3, devSim.natom);
    RDy = LOC2(devSim.xyz, 1 , devSim.katom[LL]-1, 3, devSim.natom);
    RDz = LOC2(devSim.xyz, 2 , devSim.katom[LL]-1, 3, devSim.natom);
    
    int III1 = LOC2(devSim.Qsbasis, II, I, devSim.nshell, 4);
    int III2 = LOC2(devSim.Qfbasis, II, I, devSim.nshell, 4);
    int JJJ1 = LOC2(devSim.Qsbasis, JJ, J, devSim.nshell, 4);
    int JJJ2 = LOC2(devSim.Qfbasis, JJ, J, devSim.nshell, 4);
    int KKK1 = LOC2(devSim.Qsbasis, KK, K, devSim.nshell, 4);
    int KKK2 = LOC2(devSim.Qfbasis, KK, K, devSim.nshell, 4);
    int LLL1 = LOC2(devSim.Qsbasis, LL, L, devSim.nshell, 4);
    int LLL2 = LOC2(devSim.Qfbasis, LL, L, devSim.nshell, 4);
    
    
    // maxIJKL is the max of I,J,K,L
    int maxIJKL = (int)MAX(MAX(I,J),MAX(K,L));
    
    if (((maxIJKL == 2)&&(J != 0 || L!=0)) || (maxIJKL >= 3)) {
        IJKLTYPE = 999;
    }
    
    QUICKDouble hybrid_coeff = 0.0;
    if (devSim.method == HF){
        hybrid_coeff = 1.0;
    }else if (devSim.method == B3LYP){
        hybrid_coeff = 0.2;
    }else if (devSim.method == DFT){
        hybrid_coeff = 0.0;
    }
    
    for (int III = III1; III <= III2; III++) {
        for (int JJJ = MAX(III,JJJ1); JJJ <= JJJ2; JJJ++) {
            for (int KKK = MAX(III,KKK1); KKK <= KKK2; KKK++) {
                for (int LLL = MAX(KKK,LLL1); LLL <= LLL2; LLL++) {
                    
                    if (III < KKK ||
                        ((III == JJJ) && (III == LLL)) ||
                        ((III == JJJ) && (III  < LLL)) ||
                        ((JJJ == LLL) && (III  < JJJ)) ||
                        ((III == KKK) && (III  < JJJ)  && (JJJ < LLL))) {
                        
                        
                        QUICKDouble Y = (QUICKDouble) hrrwhole( I, J, K, L,\
                                                               III, JJJ, KKK, LLL, IJKLTYPE, store, \
                                                               RAx, RAy, RAz, RBx, RBy, RBz, \
                                                               RCx, RCy, RCz, RDx, RDy, RDz);
                        if (abs(Y) > devSim.integralCutoff) {
                            addint(devSim.oULL, Y, III, JJJ, KKK, LLL, hybrid_coeff, devSim.dense, devSim.nbasis);
                        }
                        
                    }
                }
            }
        }
    }
    return;
}



#ifdef int_spd
__global__ void getAOInt_kernel(QUICKULL intStart, QUICKULL intEnd, ERI_entry* aoint_buffer, int streamID)
#elif defined int_spdf
__global__ void getAOInt_kernel_spdf(QUICKULL intStart, QUICKULL intEnd, ERI_entry* aoint_buffer, int streamID)
#elif defined int_spdf2
__global__ void getAOInt_kernel_spdf2(QUICKULL intStart, QUICKULL intEnd, ERI_entry* aoint_buffer, int streamID)
#elif defined int_spdf3
__global__ void getAOInt_kernel_spdf3(QUICKULL intStart, QUICKULL intEnd, ERI_entry* aoint_buffer, int streamID)
#elif defined int_spdf4
__global__ void getAOInt_kernel_spdf4(QUICKULL intStart, QUICKULL intEnd, ERI_entry* aoint_buffer, int streamID)
#endif
{
    
    unsigned int offside = blockIdx.x*blockDim.x+threadIdx.x;
    int totalThreads = blockDim.x*gridDim.x;
    
    QUICKULL jshell         = (QUICKULL) devSim.sqrQshell;
    QUICKULL myInt          = (QUICKULL) (intEnd - intStart + 1) / totalThreads;
    
    
    
    if ((intEnd - intStart + 1 - myInt*totalThreads)> offside) myInt++;
    
    for (QUICKULL i = 1; i<=myInt; i++) {
        QUICKULL currentInt = totalThreads * (i-1) + offside + intStart;
        QUICKULL a = (QUICKULL) currentInt/jshell;
        QUICKULL b = (QUICKULL) (currentInt - a*jshell);
        
        int II = devSim.sorted_YCutoffIJ[a].x;
        int JJ = devSim.sorted_YCutoffIJ[a].y;
        int KK = devSim.sorted_YCutoffIJ[b].x;
        int LL = devSim.sorted_YCutoffIJ[b].y;
        
        int ii = devSim.sorted_Q[II];
        int jj = devSim.sorted_Q[JJ];
        int kk = devSim.sorted_Q[KK];
        int ll = devSim.sorted_Q[LL];
        
        if (ii<=kk) {
            int nshell = devSim.nshell;
            
            if ((LOC2(devSim.YCutoff, kk, ll, nshell, nshell) * LOC2(devSim.YCutoff, ii, jj, nshell, nshell))> devSim.leastIntegralCutoff) {
                
                int iii = devSim.sorted_Qnumber[II];
                int jjj = devSim.sorted_Qnumber[JJ];
                int kkk = devSim.sorted_Qnumber[KK];
                int lll = devSim.sorted_Qnumber[LL];
#ifdef int_spd
        //        if (!((iii + jjj) > 4 || (kkk + lll) > 4)) {
                    iclass_AOInt(iii, jjj, kkk, lll, ii, jj, kk, ll, 1.0, aoint_buffer, streamID);
        //        }
#elif defined int_spdf
                if ((iii + jjj) > 4 || (kkk + lll) > 4) {
                    iclass_AOInt_spdf(iii, jjj, kkk, lll, ii, jj, kk, ll, 1.0, aoint_buffer, streamID);
                }
#elif defined int_spdf2
                if ((iii + jjj) > 4 || (kkk + lll) > 4) {
                    iclass_AOInt_spdf2(iii, jjj, kkk, lll, ii, jj, kk, ll, 1.0, aoint_buffer, streamID);
                }
#elif defined int_spdf3
                if ((iii + jjj) > 4 || (kkk + lll) > 4) {
                    iclass_AOInt_spdf3(iii, jjj, kkk, lll, ii, jj, kk, ll, 1.0, aoint_buffer, streamID);
                }
#elif defined int_spdf4
                if ((iii + jjj) > 4 || (kkk + lll) > 4) {
                    iclass_AOInt_spdf4(iii, jjj, kkk, lll, ii, jj, kk, ll, 1.0, aoint_buffer, streamID);
                }
#endif
            }
        }
    }
}






/*
 iclass subroutine is to generate 2-electron intergral using HRR and VRR method, which is the most
 performance algrithem for electron intergral evaluation. See description below for details
 */
#ifdef int_spd
__device__ __forceinline__ void iclass_AOInt
#elif defined int_spdf
__device__ __forceinline__ void iclass_AOInt_spdf
#elif defined int_spdf2
__device__ __forceinline__ void iclass_AOInt_spdf2
#elif defined int_spdf3
__device__ __forceinline__ void iclass_AOInt_spdf3
#elif defined int_spdf4
__device__ __forceinline__ void iclass_AOInt_spdf4
#endif
                            (int I, int J, int K, int L, unsigned int II, unsigned int JJ, unsigned int KK, unsigned int LL, QUICKDouble DNMax, ERI_entry* aoint_buffer, int streamID)
{
    
    /*
     kAtom A, B, C ,D is the coresponding atom for shell ii, jj, kk, ll
     and be careful with the index difference between Fortran and C++,
     Fortran starts array index with 1 and C++ starts 0.
     
     
     RA, RB, RC, and RD are the coordinates for atom katomA, katomB, katomC and katomD,
     which means they are corrosponding coorinates for shell II, JJ, KK, and LL.
     And we don't need the coordinates now, so we will not retrieve the data now.
     */
    QUICKDouble RAx = LOC2(devSim.xyz, 0 , devSim.katom[II]-1, 3, devSim.natom);
    QUICKDouble RAy = LOC2(devSim.xyz, 1 , devSim.katom[II]-1, 3, devSim.natom);
    QUICKDouble RAz = LOC2(devSim.xyz, 2 , devSim.katom[II]-1, 3, devSim.natom);
    
    QUICKDouble RCx = LOC2(devSim.xyz, 0 , devSim.katom[KK]-1, 3, devSim.natom);
    QUICKDouble RCy = LOC2(devSim.xyz, 1 , devSim.katom[KK]-1, 3, devSim.natom);
    QUICKDouble RCz = LOC2(devSim.xyz, 2 , devSim.katom[KK]-1, 3, devSim.natom);
    
    /*
     kPrimI, J, K and L indicates the primtive gaussian function number
     kStartI, J, K, and L indicates the starting guassian function for shell I, J, K, and L.
     We retrieve from global memory and save them to register to avoid multiple retrieve.
     */
    int kPrimI = devSim.kprim[II];
    int kPrimJ = devSim.kprim[JJ];
    int kPrimK = devSim.kprim[KK];
    int kPrimL = devSim.kprim[LL];
    
    int kStartI = devSim.kstart[II]-1;
    int kStartJ = devSim.kstart[JJ]-1;
    int kStartK = devSim.kstart[KK]-1;
    int kStartL = devSim.kstart[LL]-1;
    
    
    /*
     store saves temp contracted integral as [as|bs] type. the dimension should be allocatable but because
     of cuda limitation, we can not do that now.
     
     See M.Head-Gordon and J.A.Pople, Jchem.Phys., 89, No.9 (1988) for VRR algrithem details.
     */
    QUICKDouble store[STOREDIM*STOREDIM];
    
    /*
     Initial the neccessary element for
     */
    for (int i = Sumindex[K+1]+1; i<= Sumindex[K+L+2]; i++) {
        for (int j = Sumindex[I+1]+1; j<= Sumindex[I+J+2]; j++) {
            LOC2(store, j-1, i-1, STOREDIM, STOREDIM) = 0;
        }
    }
    
    for (int i = 0; i<kPrimI*kPrimJ;i++){
        int JJJ = (int) i/kPrimI;
        int III = (int) i-kPrimI*JJJ;
        /*
         In the following comments, we have I, J, K, L denote the primitive gaussian function we use, and
         for example, expo(III, ksumtype(II)) stands for the expo for the IIIth primitive guassian function for II shell,
         we use I to express the corresponding index.
         AB = expo(I)+expo(J)
         --->                --->
         ->     expo(I) * xyz (I) + expo(J) * xyz(J)
         P  = ---------------------------------------
         expo(I) + expo(J)
         Those two are pre-calculated in CPU stage.
         
         */
        int ii_start = devSim.prim_start[II];
        int jj_start = devSim.prim_start[JJ];
        
        QUICKDouble AB = LOC2(devSim.expoSum, ii_start+III, jj_start+JJJ, devSim.prim_total, devSim.prim_total);
        QUICKDouble Px = LOC2(devSim.weightedCenterX, ii_start+III, jj_start+JJJ, devSim.prim_total, devSim.prim_total);
        QUICKDouble Py = LOC2(devSim.weightedCenterY, ii_start+III, jj_start+JJJ, devSim.prim_total, devSim.prim_total);
        QUICKDouble Pz = LOC2(devSim.weightedCenterZ, ii_start+III, jj_start+JJJ, devSim.prim_total, devSim.prim_total);
        
        /*
         X1 is the contracted coeffecient, which is pre-calcuated in CPU stage as well.
         cutoffprim is used to cut too small prim gaussian function when bring density matrix into consideration.
         */
        QUICKDouble cutoffPrim = DNMax * LOC2(devSim.cutPrim, kStartI+III, kStartJ+JJJ, devSim.jbasis, devSim.jbasis);
        QUICKDouble X1 = LOC4(devSim.Xcoeff, kStartI+III, kStartJ+JJJ, I - devSim.Qstart[II], J - devSim.Qstart[JJ], devSim.jbasis, devSim.jbasis, 2, 2);
        
        for (int j = 0; j<kPrimK*kPrimL; j++){
            int LLL = (int) j/kPrimK;
            int KKK = (int) j-kPrimK*LLL;
            
            if (cutoffPrim * LOC2(devSim.cutPrim, kStartK+KKK, kStartL+LLL, devSim.jbasis, devSim.jbasis) > devSim.primLimit) {
                /*
                 CD = expo(L)+expo(K)
                 ABCD = 1/ (AB + CD) = 1 / (expo(I)+expo(J)+expo(K)+expo(L))
                 
                 `````````````````````````AB * CD      (expo(I)+expo(J))*(expo(K)+expo(L))
                 Rou(Greek Letter) =   ----------- = ------------------------------------
                 `````````````````````````AB + CD         expo(I)+expo(J)+expo(K)+expo(L)
                 
                 ```````````````````expo(I)+expo(J)                        expo(K)+expo(L)
                 ABcom = --------------------------------  CDcom = --------------------------------
                 `````````expo(I)+expo(J)+expo(K)+expo(L)           expo(I)+expo(J)+expo(K)+expo(L)
                 
                 ABCDtemp = 1/2(expo(I)+expo(J)+expo(K)+expo(L))
                 */
                
                int kk_start = devSim.prim_start[KK];
                int ll_start = devSim.prim_start[LL];
                
                QUICKDouble CD = LOC2(devSim.expoSum, kk_start+KKK, ll_start+LLL, devSim.prim_total, devSim.prim_total);
                
                QUICKDouble ABCD = 1/(AB+CD);
                
                /*
                 X2 is the multiplication of four indices normalized coeffecient
                 */
                QUICKDouble X2 = sqrt(ABCD) * X1 * LOC4(devSim.Xcoeff, kStartK+KKK, kStartL+LLL, K - devSim.Qstart[KK], L - devSim.Qstart[LL], devSim.jbasis, devSim.jbasis, 2, 2);
                
                /*
                 Q' is the weighting center of K and L
                 ```````````````````````````--->           --->
                 ->  ------>       expo(K)*xyz(K)+expo(L)*xyz(L)
                 Q = P'(K,L)  = ------------------------------
                 `````````````````````````expo(K) + expo(L)
                 
                 W' is the weight center for I, J, K, L
                 
                 ```````````````--->             --->             --->            --->
                 ->     expo(I)*xyz(I) + expo(J)*xyz(J) + expo(K)*xyz(K) +expo(L)*xyz(L)
                 W = -------------------------------------------------------------------
                 `````````````````````````expo(I) + expo(J) + expo(K) + expo(L)
                 ``````->  ->  2
                 RPQ =| P - Q |
                 
                 ```````````->  -> 2
                 T = ROU * | P - Q|
                 */
                
                QUICKDouble Qx = LOC2(devSim.weightedCenterX, kk_start+KKK, ll_start+LLL, devSim.prim_total, devSim.prim_total);
                QUICKDouble Qy = LOC2(devSim.weightedCenterY, kk_start+KKK, ll_start+LLL, devSim.prim_total, devSim.prim_total);
                QUICKDouble Qz = LOC2(devSim.weightedCenterZ, kk_start+KKK, ll_start+LLL, devSim.prim_total, devSim.prim_total);
                
                QUICKDouble T = AB * CD * ABCD * ( quick_dsqr(Px-Qx) + quick_dsqr(Py-Qy) + quick_dsqr(Pz-Qz));
                
                QUICKDouble YVerticalTemp[VDIM1*VDIM2*VDIM3];
                FmT(I+J+K+L, T, YVerticalTemp);
                for (int i = 0; i<=I+J+K+L; i++) {
                    VY(0, 0, i) = VY(0, 0, i) * X2;
                }
                
#ifdef int_spd
                vertical(I, J, K, L, YVerticalTemp, store, \
                         Px - RAx, Py - RAy, Pz - RAz, (Px*AB+Qx*CD)*ABCD - Px, (Py*AB+Qy*CD)*ABCD - Py, (Pz*AB+Qz*CD)*ABCD - Pz, \
                         Qx - RCx, Qy - RCy, Qz - RCz, (Px*AB+Qx*CD)*ABCD - Qx, (Py*AB+Qy*CD)*ABCD - Qy, (Pz*AB+Qz*CD)*ABCD - Qz, \
                         0.5 * ABCD, 0.5 / AB, 0.5 / CD, AB * ABCD, CD * ABCD);
#elif defined int_spdf
                
                vertical_spdf(I, J, K, L, YVerticalTemp, store, \
                              Px - RAx, Py - RAy, Pz - RAz, (Px*AB+Qx*CD)*ABCD - Px, (Py*AB+Qy*CD)*ABCD - Py, (Pz*AB+Qz*CD)*ABCD - Pz, \
                              Qx - RCx, Qy - RCy, Qz - RCz, (Px*AB+Qx*CD)*ABCD - Qx, (Py*AB+Qy*CD)*ABCD - Qy, (Pz*AB+Qz*CD)*ABCD - Qz, \
                              0.5 * ABCD, 0.5 / AB, 0.5 / CD, AB * ABCD, CD * ABCD);
#elif defined int_spdf2
                
                vertical_spdf2(I, J, K, L, YVerticalTemp, store, \
                               Px - RAx, Py - RAy, Pz - RAz, (Px*AB+Qx*CD)*ABCD - Px, (Py*AB+Qy*CD)*ABCD - Py, (Pz*AB+Qz*CD)*ABCD - Pz, \
                               Qx - RCx, Qy - RCy, Qz - RCz, (Px*AB+Qx*CD)*ABCD - Qx, (Py*AB+Qy*CD)*ABCD - Qy, (Pz*AB+Qz*CD)*ABCD - Qz, \
                               0.5 * ABCD, 0.5 / AB, 0.5 / CD, AB * ABCD, CD * ABCD);
#elif defined int_spdf3
                
                vertical_spdf3(I, J, K, L, YVerticalTemp, store, \
                               Px - RAx, Py - RAy, Pz - RAz, (Px*AB+Qx*CD)*ABCD - Px, (Py*AB+Qy*CD)*ABCD - Py, (Pz*AB+Qz*CD)*ABCD - Pz, \
                               Qx - RCx, Qy - RCy, Qz - RCz, (Px*AB+Qx*CD)*ABCD - Qx, (Py*AB+Qy*CD)*ABCD - Qy, (Pz*AB+Qz*CD)*ABCD - Qz, \
                               0.5 * ABCD, 0.5 / AB, 0.5 / CD, AB * ABCD, CD * ABCD);
#elif defined int_spdf4
                
                vertical_spdf4(I, J, K, L, YVerticalTemp, store, \
                               Px - RAx, Py - RAy, Pz - RAz, (Px*AB+Qx*CD)*ABCD - Px, (Py*AB+Qy*CD)*ABCD - Py, (Pz*AB+Qz*CD)*ABCD - Pz, \
                               Qx - RCx, Qy - RCy, Qz - RCz, (Px*AB+Qx*CD)*ABCD - Qx, (Py*AB+Qy*CD)*ABCD - Qy, (Pz*AB+Qz*CD)*ABCD - Qz, \
                               0.5 * ABCD, 0.5 / AB, 0.5 / CD, AB * ABCD, CD * ABCD);
#endif
            }
        }
    }
    
    
    // IJKLTYPE is the I, J, K,L type
    int IJKLTYPE = (int) (1000 * I + 100 *J + 10 * K + L);
    
    QUICKDouble RBx, RBy, RBz;
    QUICKDouble RDx, RDy, RDz;
    
    RBx = LOC2(devSim.xyz, 0 , devSim.katom[JJ]-1, 3, devSim.natom);
    RBy = LOC2(devSim.xyz, 1 , devSim.katom[JJ]-1, 3, devSim.natom);
    RBz = LOC2(devSim.xyz, 2 , devSim.katom[JJ]-1, 3, devSim.natom);
    
    
    RDx = LOC2(devSim.xyz, 0 , devSim.katom[LL]-1, 3, devSim.natom);
    RDy = LOC2(devSim.xyz, 1 , devSim.katom[LL]-1, 3, devSim.natom);
    RDz = LOC2(devSim.xyz, 2 , devSim.katom[LL]-1, 3, devSim.natom);
    
    int III1 = LOC2(devSim.Qsbasis, II, I, devSim.nshell, 4);
    int III2 = LOC2(devSim.Qfbasis, II, I, devSim.nshell, 4);
    int JJJ1 = LOC2(devSim.Qsbasis, JJ, J, devSim.nshell, 4);
    int JJJ2 = LOC2(devSim.Qfbasis, JJ, J, devSim.nshell, 4);
    int KKK1 = LOC2(devSim.Qsbasis, KK, K, devSim.nshell, 4);
    int KKK2 = LOC2(devSim.Qfbasis, KK, K, devSim.nshell, 4);
    int LLL1 = LOC2(devSim.Qsbasis, LL, L, devSim.nshell, 4);
    int LLL2 = LOC2(devSim.Qfbasis, LL, L, devSim.nshell, 4);
    
    
    // maxIJKL is the max of I,J,K,L
    int maxIJKL = (int)MAX(MAX(I,J),MAX(K,L));
    
    if (((maxIJKL == 2)&&(J != 0 || L!=0)) || (maxIJKL >= 3)) {
        IJKLTYPE = 999;
    }
    
    
    // Store generated ERI to buffer
    for (int III = III1; III <= III2; III++) {
        for (int JJJ = MAX(III,JJJ1); JJJ <= JJJ2; JJJ++) {
            for (int KKK = MAX(III,KKK1); KKK <= KKK2; KKK++) {
                for (int LLL = MAX(KKK,LLL1); LLL <= LLL2; LLL++) {
                    if( (III < JJJ && III < KKK && KKK < LLL) ||
                       (III < KKK || JJJ <= LLL)){
                        QUICKDouble Y = (QUICKDouble) hrrwhole( I, J, K, L,\
                                                               III, JJJ, KKK, LLL, IJKLTYPE, store, \
                                                               RAx, RAy, RAz, RBx, RBy, RBz, \
                                                               RCx, RCy, RCz, RDx, RDy, RDz);
                        
                        if (abs(Y) > devSim.maxIntegralCutoff){
                            ERI_entry a;
                            a.value = Y;
                            a.IJ = (III - 1) * devSim.nbasis + JJJ - 1;
                            a.KL = (KKK - 1) * devSim.nbasis + LLL - 1;
                            
                            aoint_buffer[QUICKADD(devSim.intCount[streamID], 1)] = a;
                        }
                    }
                }
            }
        }
    }
}

/*
 iclass subroutine is to generate 2-electron intergral using HRR and VRR method, which is the most
 performance algrithem for electron intergral evaluation. See description below for details
 */
__device__ __forceinline__ void iclass_grad
(int I, int J, int K, int L, unsigned int II, unsigned int JJ, unsigned int KK, unsigned int LL, QUICKDouble DNMax)
{
    
    
    
    /*
     kAtom A, B, C ,D is the coresponding atom for shell ii, jj, kk, ll
     and be careful with the index difference between Fortran and C++,
     Fortran starts array index with 1 and C++ starts 0.
     
     
     RA, RB, RC, and RD are the coordinates for atom katomA, katomB, katomC and katomD,
     which means they are corrosponding coorinates for shell II, JJ, KK, and LL.
     And we don't need the coordinates now, so we will not retrieve the data now.
     */
    QUICKDouble RAx = LOC2(devSim.xyz, 0 , devSim.katom[II]-1, 3, devSim.natom);
    QUICKDouble RAy = LOC2(devSim.xyz, 1 , devSim.katom[II]-1, 3, devSim.natom);
    QUICKDouble RAz = LOC2(devSim.xyz, 2 , devSim.katom[II]-1, 3, devSim.natom);
    
    QUICKDouble RCx = LOC2(devSim.xyz, 0 , devSim.katom[KK]-1, 3, devSim.natom);
    QUICKDouble RCy = LOC2(devSim.xyz, 1 , devSim.katom[KK]-1, 3, devSim.natom);
    QUICKDouble RCz = LOC2(devSim.xyz, 2 , devSim.katom[KK]-1, 3, devSim.natom);
    
    /*
     kPrimI, J, K and L indicates the primtive gaussian function number
     kStartI, J, K, and L indicates the starting guassian function for shell I, J, K, and L.
     We retrieve from global memory and save them to register to avoid multiple retrieve.
     */
    int kPrimI = devSim.kprim[II];
    int kPrimJ = devSim.kprim[JJ];
    int kPrimK = devSim.kprim[KK];
    int kPrimL = devSim.kprim[LL];
    
    int kStartI = devSim.kstart[II]-1;
    int kStartJ = devSim.kstart[JJ]-1;
    int kStartK = devSim.kstart[KK]-1;
    int kStartL = devSim.kstart[LL]-1;
    
    
    /*
     store saves temp contracted integral as [as|bs] type. the dimension should be allocatable but because
     of cuda limitation, we can not do that now.
     
     See M.Head-Gordon and J.A.Pople, Jchem.Phys., 89, No.9 (1988) for VRR algrithem details.
     */
    QUICKDouble store[STOREDIM*STOREDIM];
    QUICKDouble storeAA[STOREDIM*STOREDIM];
    QUICKDouble storeBB[STOREDIM*STOREDIM];
    QUICKDouble storeCC[STOREDIM*STOREDIM];
    
    
    bool debut = true;
    
    for (int i = 0; i<kPrimI*kPrimJ;i++){
        int JJJ = (int) i/kPrimI;
        int III = (int) i-kPrimI*JJJ;
        /*
         In the following comments, we have I, J, K, L denote the primitive gaussian function we use, and
         for example, expo(III, ksumtype(II)) stands for the expo for the IIIth primitive guassian function for II shell,
         we use I to express the corresponding index.
         AB = expo(I)+expo(J)
         --->                --->
         ->     expo(I) * xyz (I) + expo(J) * xyz(J)
         P  = ---------------------------------------
         expo(I) + expo(J)
         Those two are pre-calculated in CPU stage.
         
         */
        int ii_start = devSim.prim_start[II];
        int jj_start = devSim.prim_start[JJ];
        
        QUICKDouble AA = LOC2(devSim.gcexpo, III , devSim.Ksumtype[II] - 1, MAXPRIM, devSim.nbasis);
        QUICKDouble BB = LOC2(devSim.gcexpo, JJJ , devSim.Ksumtype[JJ] - 1, MAXPRIM, devSim.nbasis);
        
        QUICKDouble AB = LOC2(devSim.expoSum, ii_start+III, jj_start+JJJ, devSim.prim_total, devSim.prim_total);
        QUICKDouble Px = LOC2(devSim.weightedCenterX, ii_start+III, jj_start+JJJ, devSim.prim_total, devSim.prim_total);
        QUICKDouble Py = LOC2(devSim.weightedCenterY, ii_start+III, jj_start+JJJ, devSim.prim_total, devSim.prim_total);
        QUICKDouble Pz = LOC2(devSim.weightedCenterZ, ii_start+III, jj_start+JJJ, devSim.prim_total, devSim.prim_total);
        
        /*
         X1 is the contracted coeffecient, which is pre-calcuated in CPU stage as well.
         cutoffprim is used to cut too small prim gaussian function when bring density matrix into consideration.
         */
        QUICKDouble cutoffPrim = DNMax * LOC2(devSim.cutPrim, kStartI+III, kStartJ+JJJ, devSim.jbasis, devSim.jbasis);
        QUICKDouble X1 = LOC4(devSim.Xcoeff, kStartI+III, kStartJ+JJJ, I - devSim.Qstart[II], J - devSim.Qstart[JJ], devSim.jbasis, devSim.jbasis, 2, 2);
        
        for (int j = 0; j<kPrimK*kPrimL; j++){
            int LLL = (int) j/kPrimK;
            int KKK = (int) j-kPrimK*LLL;
            
            if (cutoffPrim * LOC2(devSim.cutPrim, kStartK+KKK, kStartL+LLL, devSim.jbasis, devSim.jbasis) > devSim.gradCutoff) {
                
                QUICKDouble CC = LOC2(devSim.gcexpo, KKK , devSim.Ksumtype[KK] - 1, MAXPRIM, devSim.nbasis);
                /*
                 CD = expo(L)+expo(K)
                 ABCD = 1/ (AB + CD) = 1 / (expo(I)+expo(J)+expo(K)+expo(L))
                 AB * CD      (expo(I)+expo(J))*(expo(K)+expo(L))
                 Rou(Greek Letter) =   ----------- = ------------------------------------
                 AB + CD         expo(I)+expo(J)+expo(K)+expo(L)
                 
                 expo(I)+expo(J)                        expo(K)+expo(L)
                 ABcom = --------------------------------  CDcom = --------------------------------
                 expo(I)+expo(J)+expo(K)+expo(L)           expo(I)+expo(J)+expo(K)+expo(L)
                 
                 ABCDtemp = 1/2(expo(I)+expo(J)+expo(K)+expo(L))
                 */
                
                int kk_start = devSim.prim_start[KK];
                int ll_start = devSim.prim_start[LL];
                
                QUICKDouble CD = LOC2(devSim.expoSum, kk_start+KKK, ll_start+LLL, devSim.prim_total, devSim.prim_total);
                
                QUICKDouble ABCD = 1/(AB+CD);
                
                /*
                 X2 is the multiplication of four indices normalized coeffecient
                 */
                QUICKDouble X2 = sqrt(ABCD) * X1 * LOC4(devSim.Xcoeff, kStartK+KKK, kStartL+LLL, K - devSim.Qstart[KK], L - devSim.Qstart[LL], devSim.jbasis, devSim.jbasis, 2, 2);
                
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
                 
                 ->  -> 2
                 T = ROU * | P - Q|
                 */
                
                QUICKDouble Qx = LOC2(devSim.weightedCenterX, kk_start+KKK, ll_start+LLL, devSim.prim_total, devSim.prim_total);
                QUICKDouble Qy = LOC2(devSim.weightedCenterY, kk_start+KKK, ll_start+LLL, devSim.prim_total, devSim.prim_total);
                QUICKDouble Qz = LOC2(devSim.weightedCenterZ, kk_start+KKK, ll_start+LLL, devSim.prim_total, devSim.prim_total);
                
                //QUICKDouble T = AB * CD * ABCD * ( quick_dsqr(Px-Qx) + quick_dsqr(Py-Qy) + quick_dsqr(Pz-Qz));
                
                QUICKDouble YVerticalTemp[VDIM1*VDIM2*VDIM3];
                FmT(I+J+K+L+2, AB * CD * ABCD * ( quick_dsqr(Px-Qx) + quick_dsqr(Py-Qy) + quick_dsqr(Pz-Qz)), YVerticalTemp);
                
                for (int i = 0; i<=I+J+K+L+2; i++) {
                    VY(0, 0, i) = VY(0, 0, i) * X2;
                }
                
                
                QUICKDouble store2[STOREDIM*STOREDIM];
                
                vertical2(I, J + 1, K, L + 1, YVerticalTemp, store2, \
                         Px - RAx, Py - RAy, Pz - RAz, (Px*AB+Qx*CD)*ABCD - Px, (Py*AB+Qy*CD)*ABCD - Py, (Pz*AB+Qz*CD)*ABCD - Pz, \
                         Qx - RCx, Qy - RCy, Qz - RCz, (Px*AB+Qx*CD)*ABCD - Qx, (Py*AB+Qy*CD)*ABCD - Qy, (Pz*AB+Qz*CD)*ABCD - Qz, \
                         0.5 * ABCD, 0.5 / AB, 0.5 / CD, AB * ABCD, CD * ABCD);
                
                if (debut) {
                    
                    for (int i = Sumindex[K]; i< Sumindex[K+L+2]; i++) {
                        for (int j = Sumindex[I]; j< Sumindex[I+J+2]; j++) {
                            LOC2(store, j - Sumindex[I], i  - Sumindex[K], STOREDIM, STOREDIM) = LOC2(store2, j, i, STOREDIM, STOREDIM);
                        }
                    }
                    
                    
                    for (int i = Sumindex[K]; i< Sumindex[K+L+2]; i++) {
                        for (int j = Sumindex[I+1]; j< Sumindex[I+J+3]; j++) {
                            LOC2(storeAA, j - Sumindex[I], i  - Sumindex[K], STOREDIM, STOREDIM) = LOC2(store2, j, i, STOREDIM, STOREDIM) * AA * 2 ;
                            LOC2(storeBB, j - Sumindex[I], i  - Sumindex[K], STOREDIM, STOREDIM) = LOC2(store2, j, i, STOREDIM, STOREDIM) * BB * 2 ;
                        }
                    }
                    
                    for (int i = Sumindex[K+1]; i< Sumindex[K+L+3]; i++) {
                        for (int j = Sumindex[I]; j< Sumindex[I+J+2]; j++) {
                            LOC2(storeCC, j - Sumindex[I], i  - Sumindex[K], STOREDIM, STOREDIM) = LOC2(store2, j, i, STOREDIM, STOREDIM) * CC * 2 ;
                        }
                    }
                }else{
                    
                    for (int i = Sumindex[K]; i< Sumindex[K+L+2]; i++) {
                        for (int j = Sumindex[I]; j< Sumindex[I+J+2]; j++) {
                            LOC2(store, j - Sumindex[I], i  - Sumindex[K], STOREDIM, STOREDIM) += LOC2(store2, j, i, STOREDIM, STOREDIM);
                        }
                    }
                    
                    for (int i = Sumindex[K]; i<  Sumindex[K+L+2]; i++) {
                        for (int j = Sumindex[I+1]; j< Sumindex[I+J+3]; j++) {
                            LOC2(storeAA, j - Sumindex[I], i  - Sumindex[K], STOREDIM, STOREDIM) += LOC2(store2, j, i, STOREDIM, STOREDIM) * AA * 2 ;
                            LOC2(storeBB, j - Sumindex[I], i  - Sumindex[K], STOREDIM, STOREDIM) += LOC2(store2, j, i, STOREDIM, STOREDIM) * BB * 2 ;
                        }
                    }
                    for (int i = Sumindex[K+1]; i<  Sumindex[K+L+3]; i++) {
                        for (int j = Sumindex[I]; j< Sumindex[I+J+2]; j++) {
                            LOC2(storeCC, j - Sumindex[I], i  - Sumindex[K], STOREDIM, STOREDIM) += LOC2(store2, j, i, STOREDIM, STOREDIM) * CC * 2 ;
                        }
                    }
                }
                debut = false;
                
            }
        }
    }
    
    
    
    QUICKDouble AGradx = 0.0;
    QUICKDouble AGrady = 0.0;
    QUICKDouble AGradz = 0.0;
    QUICKDouble BGradx = 0.0;
    QUICKDouble BGrady = 0.0;
    QUICKDouble BGradz = 0.0;
    QUICKDouble CGradx = 0.0;
    QUICKDouble CGrady = 0.0;
    QUICKDouble CGradz = 0.0;
    
    int         AStart = (devSim.katom[II]-1) * 3;
    int         BStart = (devSim.katom[JJ]-1) * 3;
    int         CStart = (devSim.katom[KK]-1) * 3;
    int         DStart = (devSim.katom[LL]-1) * 3;
    
    
    QUICKDouble RBx, RBy, RBz;
    QUICKDouble RDx, RDy, RDz;
    
    RBx = LOC2(devSim.xyz, 0 , devSim.katom[JJ]-1, 3, devSim.natom);
    RBy = LOC2(devSim.xyz, 1 , devSim.katom[JJ]-1, 3, devSim.natom);
    RBz = LOC2(devSim.xyz, 2 , devSim.katom[JJ]-1, 3, devSim.natom);
    
    
    RDx = LOC2(devSim.xyz, 0 , devSim.katom[LL]-1, 3, devSim.natom);
    RDy = LOC2(devSim.xyz, 1 , devSim.katom[LL]-1, 3, devSim.natom);
    RDz = LOC2(devSim.xyz, 2 , devSim.katom[LL]-1, 3, devSim.natom);
    
    int III1 = LOC2(devSim.Qsbasis, II, I, devSim.nshell, 4);
    int III2 = LOC2(devSim.Qfbasis, II, I, devSim.nshell, 4);
    int JJJ1 = LOC2(devSim.Qsbasis, JJ, J, devSim.nshell, 4);
    int JJJ2 = LOC2(devSim.Qfbasis, JJ, J, devSim.nshell, 4);
    int KKK1 = LOC2(devSim.Qsbasis, KK, K, devSim.nshell, 4);
    int KKK2 = LOC2(devSim.Qfbasis, KK, K, devSim.nshell, 4);
    int LLL1 = LOC2(devSim.Qsbasis, LL, L, devSim.nshell, 4);
    int LLL2 = LOC2(devSim.Qfbasis, LL, L, devSim.nshell, 4);
    
    
    int    IJKLTYPE = 999;
    
    int         nbasis = devSim.nbasis;
    
    for (int III = III1; III <= III2; III++) {
        for (int JJJ = MAX(III,JJJ1); JJJ <= JJJ2; JJJ++) {
            for (int KKK = MAX(III,KKK1); KKK <= KKK2; KKK++) {
                for (int LLL = MAX(KKK,LLL1); LLL <= LLL2; LLL++) {
                    
                    if (III < KKK ||
                        ((III == JJJ) && (III == LLL)) ||
                        ((III == JJJ) && (III  < LLL)) ||
                        ((JJJ == LLL) && (III  < JJJ)) ||
                        ((III == KKK) && (III  < JJJ)  && (JJJ < LLL))) {
                        
                        QUICKDouble Yaax, Yaay, Yaaz;
                        QUICKDouble Ybbx, Ybby, Ybbz;
                        QUICKDouble Yccx, Yccy, Yccz;
                        
                        hrrwholegrad(&Yaax, &Yaay, &Yaaz, \
                                     &Ybbx, &Ybby, &Ybbz, \
                                     &Yccx, &Yccy, &Yccz, \
                                     I, J, K, L,\
                                     III, JJJ, KKK, LLL, IJKLTYPE, \
                                     store, storeAA, storeBB, storeCC, \
                                     RAx, RAy, RAz, RBx, RBy, RBz, \
                                     RCx, RCy, RCz, RDx, RDy, RDz);
                        
                        QUICKDouble constant = 0.0 ;
                        
                        QUICKDouble DENSEKI = (QUICKDouble) LOC2(devSim.dense, KKK-1, III-1, nbasis, nbasis);
                        QUICKDouble DENSEKJ = (QUICKDouble) LOC2(devSim.dense, KKK-1, JJJ-1, nbasis, nbasis);
                        QUICKDouble DENSELJ = (QUICKDouble) LOC2(devSim.dense, LLL-1, JJJ-1, nbasis, nbasis);
                        QUICKDouble DENSELI = (QUICKDouble) LOC2(devSim.dense, LLL-1, III-1, nbasis, nbasis);
                        QUICKDouble DENSELK = (QUICKDouble) LOC2(devSim.dense, LLL-1, KKK-1, nbasis, nbasis);
                        QUICKDouble DENSEJI = (QUICKDouble) LOC2(devSim.dense, JJJ-1, III-1, nbasis, nbasis);
                        
                        if (II < JJ && II < KK && KK < LL ||
                            ( III < KKK && III < JJJ && KKK < LLL)) {
                            constant = ( 4.0 * DENSEJI * DENSELK - DENSEKI * DENSELJ - DENSELI * DENSEKJ);
                        }else{
                            if (III < KKK) {
                                if( III == JJJ && KKK == LLL){
                                    constant = (DENSEJI * DENSELK - 0.5 * DENSEKI * DENSEKI);
                                }else if (JJJ == KKK && JJJ == LLL){
                                    constant = DENSELJ * DENSEJI;
                                }else if (KKK == LLL && III < JJJ && JJJ != KKK){
                                    constant = (2.0* DENSEJI * DENSELK - DENSEKI * DENSEKJ);
                                }else if ( III == JJJ && KKK < LLL){
                                    constant = (2.0* DENSELK * DENSEJI - DENSEKI * DENSELI);
                                }
                            }
                            else{
                                if (JJJ <= LLL) {
                                    if (III == JJJ && III == KKK && III == LLL) {
                                        // Do nothing
                                    }else if (III==JJJ && III==KKK && III < LLL){
                                        constant = DENSELI * DENSEJI;
                                    }else if (III==KKK && JJJ==LLL && III < JJJ){
                                        constant = (1.5 * DENSEJI * DENSEJI - 0.5 * DENSELJ * DENSEKI);
                                    }else if (III== KKK && III < JJJ && JJJ < LLL){
                                        constant = (3.0 * DENSEJI * DENSELI - DENSELJ * DENSEKI);
                                    }
                                }
                            }
                        }
                        
                        
                        AGradx += constant * Yaax;
                        AGrady += constant * Yaay;
                        AGradz += constant * Yaaz;
                        
                        BGradx += constant * Ybbx;
                        BGrady += constant * Ybby;
                        BGradz += constant * Ybbz;
                        
                        CGradx += constant * Yccx;
                        CGrady += constant * Yccy;
                        CGradz += constant * Yccz;
                    }
                }
            }
        }
    }
    
    
    GRADADD(devSim.gradULL[AStart], AGradx);
    GRADADD(devSim.gradULL[AStart + 1], AGrady);
    GRADADD(devSim.gradULL[AStart + 2], AGradz);
    
    
    GRADADD(devSim.gradULL[BStart], BGradx);
    GRADADD(devSim.gradULL[BStart + 1], BGrady);
    GRADADD(devSim.gradULL[BStart + 2], BGradz);
    
    
    GRADADD(devSim.gradULL[CStart], CGradx);
    GRADADD(devSim.gradULL[CStart + 1], CGrady);
    GRADADD(devSim.gradULL[CStart + 2], CGradz);
    
    
    GRADADD(devSim.gradULL[DStart], (-AGradx-BGradx-CGradx));
    GRADADD(devSim.gradULL[DStart + 1], (-AGrady-BGrady-CGrady));
    GRADADD(devSim.gradULL[DStart + 2], (-AGradz-BGradz-CGradz));
    
    //printf("II JJ KK LL= %i %i %i %i\n IJKL= %i %i %i %i\n %f %f %f \n %f %f %f \n %f %f %f \n %i %i %i %i\n", II, JJ, KK, LL, \
           I, J, K, L, AGradx, AGrady, AGradz, BGradx, BGrady, BGradz, CGradx, CGrady, CGradz, AStart, BStart, CStart, DStart);
    
    return;
}



#ifndef new_quick_2_gpu_get2e_subs_h
#define new_quick_2_gpu_get2e_subs_h



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
    
    Y = Y * devSim.cons[III-1] * devSim.cons[JJJ-1] * devSim.cons[KKK-1] * devSim.cons[LLL-1];
#endif
    return Y;
}



__device__ __forceinline__ void hrrwholegrad(QUICKDouble* Yaax, QUICKDouble* Yaay, QUICKDouble* Yaaz, \
                                             QUICKDouble* Ybbx, QUICKDouble* Ybby, QUICKDouble* Ybbz, \
                                             QUICKDouble* Yccx, QUICKDouble* Yccy, QUICKDouble* Yccz, \
                                             int I, int J, int K, int L, \
                                             int III, int JJJ, int KKK, int LLL, int IJKLTYPE,
                                             QUICKDouble* store, QUICKDouble* storeAA, QUICKDouble* storeBB, QUICKDouble* storeCC, \
                                             QUICKDouble RAx,QUICKDouble RAy,QUICKDouble RAz, \
                                             QUICKDouble RBx,QUICKDouble RBy,QUICKDouble RBz, \
                                             QUICKDouble RCx,QUICKDouble RCy,QUICKDouble RCz, \
                                             QUICKDouble RDx,QUICKDouble RDy,QUICKDouble RDz)
{
    int angularL[12], angularR[12];
    QUICKDouble coefAngularL[12], coefAngularR[12];
    
    *Yaax = 0.0;
    *Yaay = 0.0;
    *Yaaz = 0.0;
    *Ybbx = 0.0;
    *Ybby = 0.0;
    *Ybbz = 0.0;
    *Yccx = 0.0;
    *Yccy = 0.0;
    *Yccz = 0.0;
    
    QUICKDouble constant = devSim.cons[III-1] * devSim.cons[JJJ-1] * devSim.cons[KKK-1] * devSim.cons[LLL-1];
    int numAngularL, numAngularR;
    
    numAngularR = lefthrr(RCx, RCy, RCz, RDx, RDy, RDz, \
                          LOC2(devSim.KLMN,0,KKK-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,KKK-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,KKK-1,3,devSim.nbasis), \
                          LOC2(devSim.KLMN,0,LLL-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,LLL-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,LLL-1,3,devSim.nbasis), \
                          L, coefAngularR, angularR);
    
    
    //  Part A - x
    
    numAngularL = lefthrr(RAx, RAy, RAz, RBx, RBy, RBz, \
                          LOC2(devSim.KLMN,0,III-1,3,devSim.nbasis) + 1, LOC2(devSim.KLMN,1,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,III-1,3,devSim.nbasis), \
                          LOC2(devSim.KLMN,0,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,JJJ-1,3,devSim.nbasis), \
                          J, coefAngularL, angularL);
    for (int i = 0; i<numAngularL; i++) {
        for (int j = 0; j<numAngularR; j++) {
            *Yaax = *Yaax + coefAngularL[i] * coefAngularR[j] * LOC2(storeAA, angularL[i]-1  - Sumindex[I], angularR[j]-1  - Sumindex[K] , STOREDIM, STOREDIM) ;
        }
    }
    
    numAngularL = lefthrr(RAx, RAy, RAz, RBx, RBy, RBz, \
                          LOC2(devSim.KLMN,0,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,III-1,3,devSim.nbasis) + 1, LOC2(devSim.KLMN,2,III-1,3,devSim.nbasis), \
                          LOC2(devSim.KLMN,0,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,JJJ-1,3,devSim.nbasis), \
                          J, coefAngularL, angularL);
    for (int i = 0; i<numAngularL; i++) {
        for (int j = 0; j<numAngularR; j++) {
            *Yaay = *Yaay + coefAngularL[i] * coefAngularR[j] * LOC2(storeAA, angularL[i]-1-Sumindex[I], angularR[j]-1-Sumindex[K] , STOREDIM, STOREDIM) ;
        }
    }
    
    numAngularL = lefthrr(RAx, RAy, RAz, RBx, RBy, RBz, \
                          LOC2(devSim.KLMN,0,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,III-1,3,devSim.nbasis) + 1, \
                          LOC2(devSim.KLMN,0,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,JJJ-1,3,devSim.nbasis), \
                          J, coefAngularL, angularL);
    for (int i = 0; i<numAngularL; i++) {
        for (int j = 0; j<numAngularR; j++) {
            *Yaaz = *Yaaz + coefAngularL[i] * coefAngularR[j] * LOC2(storeAA, angularL[i]-1-Sumindex[I], angularR[j]-1-Sumindex[K] , STOREDIM, STOREDIM) ;
        }
    }
    
    numAngularL = lefthrr(RAx, RAy, RAz, RBx, RBy, RBz, \
                          LOC2(devSim.KLMN,0,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,III-1,3,devSim.nbasis), \
                          LOC2(devSim.KLMN,0,JJJ-1,3,devSim.nbasis) + 1, LOC2(devSim.KLMN,1,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,JJJ-1,3,devSim.nbasis), \
                          J + 1, coefAngularL, angularL);
    for (int i = 0; i<numAngularL; i++) {
        for (int j = 0; j<numAngularR; j++) {
            *Ybbx = *Ybbx + coefAngularL[i] * coefAngularR[j] * LOC2(storeBB, angularL[i]-1-Sumindex[I], angularR[j]-1-Sumindex[K] , STOREDIM, STOREDIM);
        }
    }
    
    numAngularL = lefthrr(RAx, RAy, RAz, RBx, RBy, RBz, \
                          LOC2(devSim.KLMN,0,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,III-1,3,devSim.nbasis), \
                          LOC2(devSim.KLMN,0,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,JJJ-1,3,devSim.nbasis) + 1, LOC2(devSim.KLMN,2,JJJ-1,3,devSim.nbasis), \
                          J + 1, coefAngularL, angularL);
    for (int i = 0; i<numAngularL; i++) {
        for (int j = 0; j<numAngularR; j++) {
            *Ybby = *Ybby + coefAngularL[i] * coefAngularR[j] * LOC2(storeBB, angularL[i]-1-Sumindex[I], angularR[j]-1-Sumindex[K] , STOREDIM, STOREDIM);
        }
    }
    
    numAngularL = lefthrr(RAx, RAy, RAz, RBx, RBy, RBz, \
                          LOC2(devSim.KLMN,0,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,III-1,3,devSim.nbasis), \
                          LOC2(devSim.KLMN,0,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,JJJ-1,3,devSim.nbasis) + 1, \
                          J + 1, coefAngularL, angularL);
    for (int i = 0; i<numAngularL; i++) {
        for (int j = 0; j<numAngularR; j++) {
            *Ybbz = *Ybbz + coefAngularL[i] * coefAngularR[j] * LOC2(storeBB, angularL[i]-1-Sumindex[I], angularR[j]-1-Sumindex[K] , STOREDIM, STOREDIM);
        }
    }
    
    
    
    
    if (LOC2(devSim.KLMN,0,III-1,3,devSim.nbasis) >= 1) {
        
        numAngularL = lefthrr(RAx, RAy, RAz, RBx, RBy, RBz,
                              LOC2(devSim.KLMN,0,III-1,3,devSim.nbasis) - 1, LOC2(devSim.KLMN,1,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,III-1,3,devSim.nbasis),
                              LOC2(devSim.KLMN,0,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,JJJ-1,3,devSim.nbasis),
                              J, coefAngularL, angularL);
        for (int i = 0; i<numAngularL; i++) {
            for (int j = 0; j<numAngularR; j++) {
                *Yaax = *Yaax - LOC2(devSim.KLMN,0,III-1,3,devSim.nbasis) * coefAngularL[i] * coefAngularR[j] * LOC2(store, angularL[i]-1-Sumindex[I], angularR[j]-1-Sumindex[K] , STOREDIM, STOREDIM);
            }
        }
        
    }
    
    if (LOC2(devSim.KLMN,1,III-1,3,devSim.nbasis) >= 1) {
        
        numAngularL = lefthrr(RAx, RAy, RAz, RBx, RBy, RBz,
                              LOC2(devSim.KLMN,0,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,III-1,3,devSim.nbasis) - 1, LOC2(devSim.KLMN,2,III-1,3,devSim.nbasis),
                              LOC2(devSim.KLMN,0,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,JJJ-1,3,devSim.nbasis),
                              J, coefAngularL, angularL);
        for (int i = 0; i<numAngularL; i++) {
            for (int j = 0; j<numAngularR; j++) {
                *Yaay = *Yaay - LOC2(devSim.KLMN,1,III-1,3,devSim.nbasis) * coefAngularL[i] * coefAngularR[j] * LOC2(store, angularL[i]-1-Sumindex[I], angularR[j]-1-Sumindex[K] , STOREDIM, STOREDIM);
            }
        }
    }
    
    if (LOC2(devSim.KLMN,2,III-1,3,devSim.nbasis) >= 1) {
        
        numAngularL = lefthrr(RAx, RAy, RAz, RBx, RBy, RBz,
                              LOC2(devSim.KLMN,0,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,III-1,3,devSim.nbasis) - 1,
                              LOC2(devSim.KLMN,0,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,JJJ-1,3,devSim.nbasis),
                              J, coefAngularL, angularL);
        for (int i = 0; i<numAngularL; i++) {
            for (int j = 0; j<numAngularR; j++) {
                *Yaaz = *Yaaz - LOC2(devSim.KLMN,2,III-1,3,devSim.nbasis) * coefAngularL[i] * coefAngularR[j] * LOC2(store, angularL[i]-1-Sumindex[I], angularR[j]-1-Sumindex[K] , STOREDIM, STOREDIM);
            }
        }
    }
    
    
    
    
    if (LOC2(devSim.KLMN,0,JJJ-1,3,devSim.nbasis) >= 1) {
        
        numAngularL = lefthrr(RAx, RAy, RAz, RBx, RBy, RBz,
                              LOC2(devSim.KLMN,0,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,III-1,3,devSim.nbasis),
                              LOC2(devSim.KLMN,0,JJJ-1,3,devSim.nbasis) - 1, LOC2(devSim.KLMN,1,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,JJJ-1,3,devSim.nbasis),
                              J - 1, coefAngularL, angularL);
        for (int i = 0; i<numAngularL; i++) {
            for (int j = 0; j<numAngularR; j++) {
                *Ybbx = *Ybbx - LOC2(devSim.KLMN,0,JJJ-1,3,devSim.nbasis) * coefAngularL[i] * coefAngularR[j] * LOC2(store, angularL[i]-1-Sumindex[I], angularR[j]-1-Sumindex[K] , STOREDIM, STOREDIM);
            }
        }
    }
    
    
    
    if (LOC2(devSim.KLMN,1,JJJ-1,3,devSim.nbasis) >= 1) {
        
        numAngularL = lefthrr(RAx, RAy, RAz, RBx, RBy, RBz,
                              LOC2(devSim.KLMN,0,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,III-1,3,devSim.nbasis),
                              LOC2(devSim.KLMN,0,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,JJJ-1,3,devSim.nbasis) - 1, LOC2(devSim.KLMN,2,JJJ-1,3,devSim.nbasis),
                              J - 1, coefAngularL, angularL);
        for (int i = 0; i<numAngularL; i++) {
            for (int j = 0; j<numAngularR; j++) {
                *Ybby = *Ybby - LOC2(devSim.KLMN,1,JJJ-1,3,devSim.nbasis) * coefAngularL[i] * coefAngularR[j] * LOC2(store, angularL[i]-1-Sumindex[I], angularR[j]-1-Sumindex[K] , STOREDIM, STOREDIM);
            }
        }
    }
    
    
    if (LOC2(devSim.KLMN,2,JJJ-1,3,devSim.nbasis) >= 1) {
        
        numAngularL = lefthrr(RAx, RAy, RAz, RBx, RBy, RBz,
                              LOC2(devSim.KLMN,0,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,III-1,3,devSim.nbasis),
                              LOC2(devSim.KLMN,0,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,JJJ-1,3,devSim.nbasis) - 1,
                              J - 1, coefAngularL, angularL);
        for (int i = 0; i<numAngularL; i++) {
            for (int j = 0; j<numAngularR; j++) {
                *Ybbz = *Ybbz - LOC2(devSim.KLMN,2,JJJ-1,3,devSim.nbasis) * coefAngularL[i] * coefAngularR[j] * LOC2(store, angularL[i]-1-Sumindex[I], angularR[j]-1-Sumindex[K] , STOREDIM, STOREDIM);
            }
        }
    }
    
    
    
    // KET PART =====================================
    
    // Part C - x
    
    
    numAngularL = lefthrr(RAx, RAy, RAz, RBx, RBy, RBz,
                          LOC2(devSim.KLMN,0,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,III-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,III-1,3,devSim.nbasis),
                          LOC2(devSim.KLMN,0,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,JJJ-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,JJJ-1,3,devSim.nbasis),
                          J, coefAngularL, angularL);
    
    numAngularR = lefthrr(RCx, RCy, RCz, RDx, RDy, RDz,
                          LOC2(devSim.KLMN,0,KKK-1,3,devSim.nbasis) + 1, LOC2(devSim.KLMN,1,KKK-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,KKK-1,3,devSim.nbasis),
                          LOC2(devSim.KLMN,0,LLL-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,LLL-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,LLL-1,3,devSim.nbasis),
                          L, coefAngularR, angularR);
    
    
    for (int i = 0; i<numAngularL; i++) {
        for (int j = 0; j<numAngularR; j++) {
            *Yccx = *Yccx + coefAngularL[i] * coefAngularR[j] * LOC2(storeCC, angularL[i]-1-Sumindex[I], angularR[j]-1-Sumindex[K] , STOREDIM, STOREDIM) ;
        }
    }
    
    if (LOC2(devSim.KLMN,0,KKK-1,3,devSim.nbasis) >= 1) {
        
        numAngularR = lefthrr(RCx, RCy, RCz, RDx, RDy, RDz,
                              LOC2(devSim.KLMN,0,KKK-1,3,devSim.nbasis) - 1, LOC2(devSim.KLMN,1,KKK-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,KKK-1,3,devSim.nbasis),
                              LOC2(devSim.KLMN,0,LLL-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,LLL-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,LLL-1,3,devSim.nbasis),
                              L, coefAngularR, angularR);
        
        for (int i = 0; i<numAngularL; i++) {
            for (int j = 0; j<numAngularR; j++) {
                *Yccx = *Yccx - LOC2(devSim.KLMN,0,KKK-1,3,devSim.nbasis) * coefAngularL[i] * coefAngularR[j] * LOC2(store, angularL[i]-1-Sumindex[I], angularR[j]-1-Sumindex[K] , STOREDIM, STOREDIM);
            }
        }
    }
    
    
    // Part C - y
    
    numAngularR = lefthrr(RCx, RCy, RCz, RDx, RDy, RDz,
                          LOC2(devSim.KLMN,0,KKK-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,KKK-1,3,devSim.nbasis) + 1, LOC2(devSim.KLMN,2,KKK-1,3,devSim.nbasis),
                          LOC2(devSim.KLMN,0,LLL-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,LLL-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,LLL-1,3,devSim.nbasis),
                          L, coefAngularR, angularR);
    
    
    for (int i = 0; i<numAngularL; i++) {
        for (int j = 0; j<numAngularR; j++) {
            *Yccy = *Yccy + coefAngularL[i] * coefAngularR[j] * LOC2(storeCC, angularL[i]-1-Sumindex[I], angularR[j]-1-Sumindex[K] , STOREDIM, STOREDIM) ;
        }
    }
    
    if (LOC2(devSim.KLMN,1,KKK-1,3,devSim.nbasis) >= 1) {
        
        numAngularR = lefthrr(RCx, RCy, RCz, RDx, RDy, RDz,
                              LOC2(devSim.KLMN,0,KKK-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,KKK-1,3,devSim.nbasis) - 1, LOC2(devSim.KLMN,2,KKK-1,3,devSim.nbasis),
                              LOC2(devSim.KLMN,0,LLL-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,LLL-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,LLL-1,3,devSim.nbasis),
                              L, coefAngularR, angularR);
        
        for (int i = 0; i<numAngularL; i++) {
            for (int j = 0; j<numAngularR; j++) {
                *Yccy = *Yccy - LOC2(devSim.KLMN,1,KKK-1,3,devSim.nbasis) * coefAngularL[i] * coefAngularR[j] * LOC2(store, angularL[i]-1-Sumindex[I], angularR[j]-1-Sumindex[K] , STOREDIM, STOREDIM);
            }
        }
    }
    
    // Part C - z
    
    numAngularR = lefthrr(RCx, RCy, RCz, RDx, RDy, RDz,
                          LOC2(devSim.KLMN,0,KKK-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,KKK-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,KKK-1,3,devSim.nbasis) + 1,
                          LOC2(devSim.KLMN,0,LLL-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,LLL-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,LLL-1,3,devSim.nbasis),
                          L, coefAngularR, angularR);
    
    
    for (int i = 0; i<numAngularL; i++) {
        for (int j = 0; j<numAngularR; j++) {
            *Yccz = *Yccz + coefAngularL[i] * coefAngularR[j] * LOC2(storeCC, angularL[i]-1-Sumindex[I], angularR[j]-1-Sumindex[K] , STOREDIM, STOREDIM) ;
        }
    }
    
    if (LOC2(devSim.KLMN,2,KKK-1,3,devSim.nbasis) >= 1) {
        
        numAngularR = lefthrr(RCx, RCy, RCz, RDx, RDy, RDz,
                              LOC2(devSim.KLMN,0,KKK-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,KKK-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,KKK-1,3,devSim.nbasis) - 1,
                              LOC2(devSim.KLMN,0,LLL-1,3,devSim.nbasis), LOC2(devSim.KLMN,1,LLL-1,3,devSim.nbasis), LOC2(devSim.KLMN,2,LLL-1,3,devSim.nbasis),
                              L, coefAngularR, angularR);
        
        for (int i = 0; i<numAngularL; i++) {
            for (int j = 0; j<numAngularR; j++) {
                *Yccz = *Yccz - LOC2(devSim.KLMN,2,KKK-1,3,devSim.nbasis) * coefAngularL[i] * coefAngularR[j] * LOC2(store, angularL[i]-1-Sumindex[I], angularR[j]-1-Sumindex[K] , STOREDIM, STOREDIM);
            }
        }
    }
    
    
    *Yaax = *Yaax * constant;
    *Yaay = *Yaay * constant;
    *Yaaz = *Yaaz * constant;
    
    
    *Ybbx = *Ybbx * constant;
    *Ybby = *Ybby * constant;
    *Ybbz = *Ybbz * constant;
    
    
    *Yccx = *Yccx * constant;
    *Yccy = *Yccy * constant;
    *Yccz = *Yccz * constant;
    
    
    return;
    
}


#ifndef CUDA_SP

__device__ __forceinline__ int lefthrr(QUICKDouble RAx, QUICKDouble RAy, QUICKDouble RAz,
                                       QUICKDouble RBx, QUICKDouble RBy, QUICKDouble RBz,
                                       int KLMNAx, int KLMNAy, int KLMNAz,
                                       int KLMNBx, int KLMNBy, int KLMNBz,
                                       int IJTYPE,QUICKDouble* coefAngularL, int* angularL)
{
    int numAngularL;
    
    coefAngularL[0] = 1.0;
    angularL[0] = (int) LOC3(devTrans, KLMNAx + KLMNBx, KLMNAy + KLMNBy, KLMNAz + KLMNBz, TRANSDIM, TRANSDIM, TRANSDIM);
    
    if (IJTYPE == 0) {
        numAngularL = 1;
        return numAngularL;
    }
    else if (IJTYPE == 1)
    {
        numAngularL = 2;
        
        if (KLMNBx != 0) {
            coefAngularL[1] = RAx-RBx;
        }else if(KLMNBy !=0 ){
            coefAngularL[1] = RAy-RBy;
        }else if (KLMNBz != 0) {
            coefAngularL[1] = RAz-RBz;
        }
    }
    else if (IJTYPE == 2)
    {
        if (KLMNBx == 2) {
            numAngularL = 3;
            QUICKDouble tmp = RAx - RBx;
            coefAngularL[1] = 2 * tmp;
            angularL[1] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            coefAngularL[2]= tmp * tmp;
        }else if(KLMNBy == 2) {
            numAngularL = 3;
            QUICKDouble tmp = RAy - RBy;
            coefAngularL[1] = 2 * tmp;
            angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            coefAngularL[2]= tmp * tmp;
        }else if (KLMNBz == 2 ){
            numAngularL = 3;
            QUICKDouble tmp = RAz - RBz;
            coefAngularL[1] = 2 * tmp;
            angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
            coefAngularL[2]= tmp * tmp;
        }else if (KLMNBx == 1 && KLMNBy == 1){
            numAngularL = 4;
            coefAngularL[1] = RAx - RBx;
            coefAngularL[2] = RAy - RBy;
            coefAngularL[3] = (RAx - RBx) * (RAy - RBy);
            angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            
        }else if (KLMNBx == 1 && KLMNBz == 1) {
            numAngularL = 4;
            coefAngularL[1] = RAx - RBx;
            coefAngularL[2] = RAz - RBz;
            coefAngularL[3] = (RAx - RBx) * (RAz - RBz);
            angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
            angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
        }else if (KLMNBy == 1 && KLMNBz == 1) {
            numAngularL = 4;
            coefAngularL[1] = RAy - RBy;
            coefAngularL[2] = RAz - RBz;
            coefAngularL[3] = (RAy - RBy) * (RAz - RBz);
            angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
            angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
        }
    }
    else if (IJTYPE == 3)
    {
        if (KLMNBx == 3) {
            numAngularL = 4;
            QUICKDouble tmp = RAx - RBx;
            
            coefAngularL[1] = 3 * tmp;
            coefAngularL[2] = 3 * tmp * tmp;
            coefAngularL[3] = tmp * tmp * tmp;
            
            angularL[1] = (int) LOC3(devTrans, KLMNAx+2, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
        }else if (KLMNBy == 3) {
            numAngularL = 4;
            QUICKDouble tmp = RAy - RBy;
            coefAngularL[1] = 3 * tmp;
            coefAngularL[2] = 3 * tmp * tmp;
            coefAngularL[3] = tmp * tmp * tmp;
            
            angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy+2, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
        }else if (KLMNBz == 3) {
            numAngularL = 4;
            
            QUICKDouble tmp = RAz - RBz;
            coefAngularL[1] = 3 * tmp;
            coefAngularL[2] = 3 * tmp * tmp;
            coefAngularL[3] = tmp * tmp * tmp;
            
            angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
            angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
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
        }
        
    }
    else if (IJTYPE == 4)
    {
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
        }
    }


    angularL[numAngularL - 1] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);

    return numAngularL;
}


__device__ __forceinline__ int lefthrr2(QUICKDouble RAx, QUICKDouble RAy, QUICKDouble RAz,
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
            angularL[0] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            break;
        }
        case 1:
        {
            numAngularL = 2;
            angularL[0] = (int) LOC3(devTrans, KLMNAx + KLMNBx, KLMNAy + KLMNBy, KLMNAz + KLMNBz, TRANSDIM, TRANSDIM, TRANSDIM);
            angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            break;
        }
        case 2:
        {
            angularL[0] = (int) LOC3(devTrans, KLMNAx + KLMNBx, KLMNAy + KLMNBy, KLMNAz + KLMNBz, TRANSDIM, TRANSDIM, TRANSDIM);
            
            if (KLMNBx == 2) {
                numAngularL = 3;
                angularL[1] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if(KLMNBy == 2) {
                numAngularL = 3;
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBz == 2 ){
                numAngularL = 3;
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBx == 1 && KLMNBy == 1){
                numAngularL = 4;
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                
            }else if (KLMNBx == 1 && KLMNBz == 1) {
                numAngularL = 4;
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBy == 1 && KLMNBz == 1) {
                numAngularL = 4;
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }
            break;
        }
        case 3:
        {
            angularL[0] = (int) LOC3(devTrans, KLMNAx + KLMNBx, KLMNAy + KLMNBy, KLMNAz + KLMNBz, TRANSDIM, TRANSDIM, TRANSDIM);
            if (KLMNBx == 3) {
                numAngularL = 4;
                angularL[1] = (int) LOC3(devTrans, KLMNAx+2, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBy == 3) {
                numAngularL = 4;
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy+2, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBz == 3) {
                numAngularL = 4;
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBx == 1 && KLMNBy ==2) { // case 120
                numAngularL = 6;
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy+2, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBx == 1 && KLMNBz ==2) { // case 102
                numAngularL = 6;
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBy == 1 && KLMNBz ==2) { // case 012
                numAngularL = 6;
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBy== 2 && KLMNBz == 1) { // case 021
                numAngularL = 6;
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy+2, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBx == 2 && KLMNBy == 1) { // case 210
                numAngularL = 6;
                angularL[1] = (int) LOC3(devTrans, KLMNAx+2, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBx == 2 && KLMNBz ==1) { // case 201
                numAngularL = 6;
                angularL[1] = (int) LOC3(devTrans, KLMNAx+2, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBx == 1 && KLMNBy == 1) {
                numAngularL = 8;
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
                angularL[1] = (int) LOC3(devTrans, KLMNAx+3, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+2, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBy == 4) {
                numAngularL = 5;
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy+3, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy+2, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                
            }else if (KLMNBz == 4) {
                numAngularL = 5;
                angularL[1] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+3, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                
            }else if (KLMNBx == 1 && KLMNBy ==3) {
                numAngularL = 8;
                
                angularL[1] = (int) LOC3(devTrans, KLMNAx,   KLMNAy+3, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy+2, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx,   KLMNAy+2, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int) LOC3(devTrans, KLMNAx,   KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[6] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[7] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
            }else if (KLMNBx == 1 && KLMNBz ==3) {
                numAngularL = 8;
                
                angularL[1] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz+3, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[6] = (int) LOC3(devTrans, KLMNAx+1, KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[7] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                
            }else if (KLMNBy == 1 && KLMNBz ==3) {
                numAngularL = 8;
                
                angularL[1] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz+3, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[2] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[3] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz+2, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[4] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[5] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz+1, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[6] = (int) LOC3(devTrans, KLMNAx, KLMNAy+1, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                angularL[7] = (int) LOC3(devTrans, KLMNAx,   KLMNAy, KLMNAz, TRANSDIM, TRANSDIM, TRANSDIM);
                
            }else if (KLMNBx == 2 && KLMNBy == 2) {
                numAngularL = 9;
                
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
#endif
