/*
 *  gpu_startup.h
 *  new_quick
 *
 *  Created by Yipu Miao on 4/20/11.
 *  Copyright 2011 University of Florida. All rights reserved.
 *
 */

#include "../config.h"
#include "gpu_type.h"
//#include "cuda.h"

// device initial and shutdown operation
extern "C" void gpu_set_device_(int* gpu_dev_id);
extern "C" void gpu_startup_(void);
extern "C" void gpu_init_(void);
extern "C" void gpu_shutdown_(void);


// molecule, basis sets, and some other information
extern "C" void gpu_upload_atom_and_chg_(int* atom, QUICKDouble* atom_chg);
extern "C" void gpu_upload_cutoff_(QUICKDouble* cutMatrix, QUICKDouble* integralCutoff,QUICKDouble* primLimit);
extern "C" void gpu_upload_cutoff_matrix_(QUICKDouble* YCutoff,QUICKDouble* cutPrim);

extern "C" void gpu_upload_calculated_(QUICKDouble* o, QUICKDouble* co, QUICKDouble* vec, QUICKDouble* dense);
extern "C" void gpu_upload_basis_(int* nshell, int* nprim, int* jshell, int* jbasis, int* maxcontract, \
                                  int* ncontract, int* itype,     QUICKDouble* aexp,      QUICKDouble* dcoeff,\
                                  int* first_basis_function, int* last_basis_function, int* first_shell_basis_function, int* last_shell_basis_function, \
                                  int* ncenter,   int* kstart,    int* katom,     int* ktype,     int* kprim,  int* kshell, int* Ksumtype, \
                                  int* Qnumber,   int* Qstart,    int* Qfinal,    int* Qsbasis,   int* Qfbasis,\
                                  QUICKDouble* gccoeff,           QUICKDouble* cons,      QUICKDouble* gcexpo, int* KLMN);


extern "C" void gpu_get2e_(QUICKDouble* o);
void get2e(_gpu_type gpu);
void upload_sim_to_constant(_gpu_type gpu);
void upload_para_to_const();
__global__ void get2e_kernel();

__device__ void gpu_shell(unsigned int II, unsigned int JJ, unsigned int KK, unsigned int LL);

__device__ void FmT(int MaxM, QUICKDouble X, QUICKDouble* vals, QUICKDouble sqrtABCD);
__device__ void vertical(int NABCDTYPE, QUICKDouble* YVerticalTemp, QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,\
              QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz,\
              QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,\
              QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
              QUICKDouble ABCDtemp,QUICKDouble ABtemp,QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);

__device__ void iclass(int I, int J, int K, int L, unsigned int II, unsigned int JJ, unsigned int KK, unsigned int LL, QUICKDouble DNMax);


__device__ void PSSS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Ptempx, QUICKDouble Ptempy,QUICKDouble Ptempz, \
           QUICKDouble WPtempx, QUICKDouble WPtempy, QUICKDouble WPtempz);
__device__ void SSPS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Qtempx, QUICKDouble Qtempy,QUICKDouble Qtempz, \
           QUICKDouble WQtempx, QUICKDouble WQtempy, QUICKDouble WQtempz);
__device__ void PSPS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Ptempx, QUICKDouble Ptempy,QUICKDouble Ptempz, \
           QUICKDouble WPtempx, QUICKDouble WPtempy, QUICKDouble WPtempz, QUICKDouble ABCDtemp);
__device__ void DSSS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Ptempx, QUICKDouble Ptempy,QUICKDouble Ptempz, \
           QUICKDouble WPtempx, QUICKDouble WPtempy, QUICKDouble WPtempz, QUICKDouble ABtemp, QUICKDouble CDcom);
__device__ void SSDS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Qtempx, QUICKDouble Qtempy,QUICKDouble Qtempz, \
           QUICKDouble WQtempx, QUICKDouble WQtempy, QUICKDouble WQtempz, QUICKDouble CDtemp, QUICKDouble ABcom);
__device__ void DSPS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Qtempx, QUICKDouble Qtempy,QUICKDouble Qtempz, \
           QUICKDouble WQtempx, QUICKDouble WQtempy, QUICKDouble WQtempz, QUICKDouble ABCDtemp);
__device__ void PSDS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Ptempx, QUICKDouble Ptempy,QUICKDouble Ptempz, \
           QUICKDouble WPtempx, QUICKDouble WPtempy, QUICKDouble WPtempz, QUICKDouble ABCDtemp);
__device__ void DSDS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble Ptempx, QUICKDouble Ptempy,QUICKDouble Ptempz, \
           QUICKDouble WPtempx, QUICKDouble WPtempy, QUICKDouble WPtempz, QUICKDouble ABCDtemp, QUICKDouble ABtemp, QUICKDouble CDcom);
__device__ QUICKDouble hrrwhole(int III, int JJJ, int KKK, int LLL, int IJKLTYPE, QUICKDouble* store, \
                     QUICKDouble RAx,QUICKDouble RAy,QUICKDouble RAz, \
                     QUICKDouble RBx,QUICKDouble RBy,QUICKDouble RBz, \
                     QUICKDouble RCx,QUICKDouble RCy,QUICKDouble RCz, \
                     QUICKDouble RDx,QUICKDouble RDy,QUICKDouble RDz);

