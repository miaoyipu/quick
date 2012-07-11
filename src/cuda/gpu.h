/*
 *  gpu_startup.h
 *  new_quick
 *
 *  Created by Yipu Miao on 4/20/11.
 *  Copyright 2011 University of Florida. All rights reserved.
 *
 */
#ifndef QUICK_GPU_H
#define QUICK_GPU_H

#include "../config.h"
#include "gpu_type.h"
//#include "cuda.h"

// device initial and shutdown operation
extern "C" void gpu_set_device_(int* gpu_dev_id);
extern "C" void gpu_startup_(void);
extern "C" void gpu_init_(void);
extern "C" void gpu_shutdown_(void);

extern "C" void gpu_get_device_info_(int* gpu_dev_count, int* gpu_dev_id,int* gpu_dev_mem,
                                     int* gpu_num_proc,double* gpu_core_freq,char* gpu_dev_name,int* name_len, int* majorv, int* minorv);


// molecule, basis sets, and some other information
extern "C" void gpu_upload_method_(int* quick_method);
extern "C" void gpu_upload_atom_and_chg_(int* atom, QUICKDouble* atom_chg);
extern "C" void gpu_upload_cutoff_(QUICKDouble* cutMatrix, QUICKDouble* integralCutoff,QUICKDouble* primLimit, QUICKDouble* DMCutoff);
extern "C" void gpu_upload_cutoff_matrix_(QUICKDouble* YCutoff,QUICKDouble* cutPrim);
extern "C" void gpu_upload_energy_(QUICKDouble* E);
extern "C" void gpu_upload_calculated_(QUICKDouble* o, QUICKDouble* co, QUICKDouble* vec, QUICKDouble* dense);
extern "C" void gpu_upload_basis_(int* nshell, int* nprim, int* jshell, int* jbasis, int* maxcontract, \
                                  int* ncontract, int* itype,     QUICKDouble* aexp,      QUICKDouble* dcoeff,\
                                  int* first_basis_function, int* last_basis_function, int* first_shell_basis_function, int* last_shell_basis_function, \
                                  int* ncenter,   int* kstart,    int* katom,     int* ktype,     int* kprim,  int* kshell, int* Ksumtype, \
                                  int* Qnumber,   int* Qstart,    int* Qfinal,    int* Qsbasis,   int* Qfbasis,\
                                  QUICKDouble* gccoeff,           QUICKDouble* cons,      QUICKDouble* gcexpo, int* KLMN);



void get2e(_gpu_type gpu);
void getxc(_gpu_type gpu);
void get2e_MP2(_gpu_type gpu);

extern "C" void gpu_get2e_(QUICKDouble* o);
extern "C" void gpu_getxc_(int* isg, QUICKDouble* sigrad2, QUICKDouble* Eelxc, QUICKDouble* aelec, QUICKDouble* belec, QUICKDouble *o);
extern "C" void gpu_mp2_(QUICKDouble* EMP2);

__global__ void get2e_kernel();
__global__ void getxc_kernel();
__global__ void get2e_MP2_kernel();


__device__ void iclass(int I, int J, int K, int L, unsigned int II, unsigned int JJ, unsigned int KK, unsigned int LL, QUICKDouble DNMax);
__device__ void iclass_MP2(int I, int J, int K, int L, unsigned int II, unsigned int JJ, unsigned int KK, unsigned int LL, QUICKDouble DNMax);

void upload_sim_to_constant(_gpu_type gpu);
void upload_sim_to_constant_dft(_gpu_type gpu);
void upload_sim_to_constant_MP2(_gpu_type gpu);

void upload_para_to_const();
void upload_para_to_const_MP2();


//__device__ void gpu_shell(unsigned int II, unsigned int JJ, unsigned int KK, unsigned int LL);

__device__ void FmT(int MaxM, QUICKDouble X, QUICKDouble* vals);
__device__ void vertical(int NABCDTYPE, QUICKDouble* YVerticalTemp, QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,\
              QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz,\
              QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,\
              QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
              QUICKDouble ABCDtemp,QUICKDouble ABtemp,QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);

__device__ void FmT_MP2(int MaxM, QUICKDouble X, QUICKDouble* vals);
__device__ void vertical_MP2(int NABCDTYPE, QUICKDouble* YVerticalTemp, QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,\
                         QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz,\
                         QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,\
                         QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                         QUICKDouble ABCDtemp,QUICKDouble ABtemp,QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);



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

__device__ QUICKDouble hrrwhole(int I, int J, int K, int L, \
                                int III, int JJJ, int KKK, int LLL, int IJKLTYPE, QUICKDouble* store, \
                                QUICKDouble RAx,QUICKDouble RAy,QUICKDouble RAz, \
                                QUICKDouble RBx,QUICKDouble RBy,QUICKDouble RBz, \
                                QUICKDouble RCx,QUICKDouble RCy,QUICKDouble RCz, \
                                QUICKDouble RDx,QUICKDouble RDy,QUICKDouble RDz);                                                         

__device__ void vertical(int I, int J, int K, int L, QUICKDouble* YVerticalTemp, QUICKDouble* store, \
              QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
              QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
              QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
              QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
              QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
              QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);

__device__ int lefthrr(QUICKDouble RAx, QUICKDouble RAy, QUICKDouble RAz, 
                       QUICKDouble RBx, QUICKDouble RBy, QUICKDouble RBz,
                       int KLMNAx, int KLMNAy, int KLMNAz,
                       int KLMNBx, int KLMNBy, int KLMNBz,
                       int IJTYPE,QUICKDouble* coefAngularL, int* angularL);

__device__ QUICKDouble hrrwhole_MP2(int I, int J, int K, int L, \
                                int III, int JJJ, int KKK, int LLL, int IJKLTYPE, QUICKDouble* store, \
                                QUICKDouble RAx,QUICKDouble RAy,QUICKDouble RAz, \
                                QUICKDouble RBx,QUICKDouble RBy,QUICKDouble RBz, \
                                QUICKDouble RCx,QUICKDouble RCy,QUICKDouble RCz, \
                                QUICKDouble RDx,QUICKDouble RDy,QUICKDouble RDz);                                                         

__device__ void vertical_MP2(int I, int J, int K, int L, QUICKDouble* YVerticalTemp, QUICKDouble* store, \
                         QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
                         QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                         QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
                         QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                         QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                         QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);

__device__ int lefthrr_MP2(QUICKDouble RAx, QUICKDouble RAy, QUICKDouble RAz, 
                       QUICKDouble RBx, QUICKDouble RBy, QUICKDouble RBz,
                       int KLMNAx, int KLMNAy, int KLMNAz,
                       int KLMNBx, int KLMNBy, int KLMNBz,
                       int IJTYPE,QUICKDouble* coefAngularL, int* angularL);
/*

__device__ void vertical_case_12(int I, int J, int K, int L, QUICKDouble* YVerticalTemp, QUICKDouble* store,
                                 QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz, \
                                 QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                                 QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz, \
                                 QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                                 QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                                 QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);

__device__ void vertical_case_21(int I, int J, int K, int L, QUICKDouble* YVerticalTemp, QUICKDouble* store,
                                 QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
                                 QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                                 QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
                                 QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                                 QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                                 QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);


__device__ void vertical_case_22(int I, int J, int K, int L, QUICKDouble* YVerticalTemp, QUICKDouble* store,
                                 QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
                                 QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                                 QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
                                 QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                                 QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                                 QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);

__device__ void vertical_case_30(int I, int J, int K, int L, QUICKDouble* YVerticalTemp, QUICKDouble* store,
                                 QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
                                 QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                                 QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
                                 QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                                 QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                                 QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);

__device__ void vertical_case_3(int I, int J, int K, int L, QUICKDouble* YVerticalTemp, QUICKDouble* store,
                                QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
                                QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                                QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
                                QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                                QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                                QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);

__device__ void vertical_case_40(int I, int J, int K, int L, QUICKDouble* YVerticalTemp, QUICKDouble* store,
                                 QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
                                 QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                                 QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
                                 QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                                 QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                                 QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);

__device__ void vertical_case_4(int I, int J, int K, int L, QUICKDouble* YVerticalTemp, QUICKDouble* store,
                                QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
                                QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                                QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
                                QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                                QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                                QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);                                                         

__device__ void vertical_case_31(int I, int J, int K, int L, QUICKDouble* YVerticalTemp, QUICKDouble* store,
                                 QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
                                 QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                                 QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
                                 QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                                 QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                                 QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);                                                         

__device__ void vertical_case_13(int I, int J, int K, int L, QUICKDouble* YVerticalTemp, QUICKDouble* store,
                                 QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
                                 QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                                 QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
                                 QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                                 QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                                 QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);                                                        


__device__ void vertical_case_41(int I, int J, int K, int L, QUICKDouble* YVerticalTemp, QUICKDouble* store,
                                 QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
                                 QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                                 QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
                                 QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                                 QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                                 QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);                                                         


__device__ void vertical_case_14(int I, int J, int K, int L, QUICKDouble* YVerticalTemp, QUICKDouble* store,
                                 QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
                                 QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                                 QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
                                 QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                                 QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                                 QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);                                                         

__device__ void vertical_case_32(int I, int J, int K, int L, QUICKDouble* YVerticalTemp, QUICKDouble* store,
                                 QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
                                 QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                                 QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
                                 QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                                 QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                                 QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);                                                         

__device__ void vertical_case_23(int I, int J, int K, int L, QUICKDouble* YVerticalTemp, QUICKDouble* store,
                                 QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
                                 QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                                 QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
                                 QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                                 QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                                 QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);                                                         

__device__ void vertical_case_42(int I, int J, int K, int L, QUICKDouble* YVerticalTemp, QUICKDouble* store,
                                 QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
                                 QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                                 QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
                                 QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                                 QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                                 QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);                                                         

__device__ void vertical_case_24(int I, int J, int K, int L, QUICKDouble* YVerticalTemp, QUICKDouble* store,
                                 QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
                                 QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                                 QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
                                 QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                                 QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                                 QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);                                                         

__device__ void vertical_case_33(int I, int J, int K, int L, QUICKDouble* YVerticalTemp, QUICKDouble* store,
                                 QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
                                 QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                                 QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
                                 QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                                 QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                                 QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);                                                         

__device__ void vertical_case_43(int I, int J, int K, int L, QUICKDouble* YVerticalTemp, QUICKDouble* store,
                                 QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz,  \
                                 QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                                 QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  \
                                 QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                                 QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                                 QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);                                                         

__device__ void vertical_case_34(int I, int J, int K, int L, QUICKDouble* YVerticalTemp, QUICKDouble* store,
                                 QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz, \
                                 QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                                 QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz, \
                                 QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                                 QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                                 QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);                                                         

__device__ void vertical_case_44(int I, int J, int K, int L, QUICKDouble* YVerticalTemp, QUICKDouble* store, 
                                 QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz, \
                                 QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, \
                                 QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz, \
                                 QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, \
                                 QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
                                 QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);
*/

//__device__ void gpu_grid_b3lyp(int irad, int iradtemp, int iatm);
__device__ void gpu_grid_xc(int irad, int iradtemp, int iatm, QUICKDouble XAng, QUICKDouble YAng, QUICKDouble ZAng, QUICKDouble WAng);

__device__ int gridFormSG1(int iitype, QUICKDouble distance, \
                           QUICKDouble* XAng, QUICKDouble* YAng, QUICKDouble* ZAng, QUICKDouble* WAng);

__device__ QUICKDouble SSW( QUICKDouble gridx, QUICKDouble gridy, QUICKDouble gridz, int atm);
__device__ void pteval(QUICKDouble gridx, QUICKDouble gridy, QUICKDouble gridz, 
                       QUICKDouble* phi, QUICKDouble* dphidx, QUICKDouble* dphidy,  QUICKDouble* dphidz, 
                       int ibas);

__device__ void denspt(QUICKDouble gridx, QUICKDouble gridy, QUICKDouble gridz, QUICKDouble* density, QUICKDouble* densityb, 
                       QUICKDouble* gax,   QUICKDouble* gay,   QUICKDouble* gaz,   QUICKDouble* gbx,     QUICKDouble* gby,     QUICKDouble* gbz);
__device__ QUICKDouble b3lyp_e(QUICKDouble rho, QUICKDouble sigma);

__device__ QUICKDouble b3lypf(QUICKDouble rho, QUICKDouble sigma, QUICKDouble* dfdr);
__device__ int gen_oh(int code, int num, QUICKDouble* x, QUICKDouble* y, QUICKDouble* z, QUICKDouble* w, QUICKDouble a, QUICKDouble b, QUICKDouble v);

__device__ void LD0006(QUICKDouble* x, QUICKDouble* y, QUICKDouble* z, QUICKDouble* w, int N);
__device__ void LD0014(QUICKDouble* x, QUICKDouble* y, QUICKDouble* z, QUICKDouble* w, int N);
__device__ void LD0026(QUICKDouble* x, QUICKDouble* y, QUICKDouble* z, QUICKDouble* w, int N);
__device__ void LD0038(QUICKDouble* x, QUICKDouble* y, QUICKDouble* z, QUICKDouble* w, int N);
__device__ void LD0050(QUICKDouble* x, QUICKDouble* y, QUICKDouble* z, QUICKDouble* w, int N);
__device__ void LD0074(QUICKDouble* x, QUICKDouble* y, QUICKDouble* z, QUICKDouble* w, int N);
__device__ void LD0086(QUICKDouble* x, QUICKDouble* y, QUICKDouble* z, QUICKDouble* w, int N);
__device__ void LD0110(QUICKDouble* x, QUICKDouble* y, QUICKDouble* z, QUICKDouble* w, int N);
__device__ void LD0146(QUICKDouble* x, QUICKDouble* y, QUICKDouble* z, QUICKDouble* w, int N);
__device__ void LD0170(QUICKDouble* x, QUICKDouble* y, QUICKDouble* z, QUICKDouble* w, int N);
__device__ void LD0194(QUICKDouble* x, QUICKDouble* y, QUICKDouble* z, QUICKDouble* w, int N);


__device__ void lyp(QUICKDouble pa, QUICKDouble pb, QUICKDouble gax, QUICKDouble gay, QUICKDouble gaz, QUICKDouble gbx, QUICKDouble gby, QUICKDouble gbz,
                    QUICKDouble* dfdr, QUICKDouble* dfdgg, QUICKDouble* dfdggo);
__device__ void becke(QUICKDouble density, QUICKDouble gx, QUICKDouble gy, QUICKDouble gz, QUICKDouble gotherx, QUICKDouble gothery, QUICKDouble gotherz,
                      QUICKDouble* dfdr, QUICKDouble* dfdgg, QUICKDouble* dfdggo);
__device__ QUICKDouble lyp_e(QUICKDouble pa, QUICKDouble pb, QUICKDouble gax, QUICKDouble gay, QUICKDouble gaz,
                             QUICKDouble gbx,     QUICKDouble gby,      QUICKDouble gbz);

__device__ QUICKDouble becke_e(QUICKDouble density, QUICKDouble densityb, QUICKDouble gax, QUICKDouble gay, QUICKDouble gaz,
                               QUICKDouble gbx,     QUICKDouble gby,      QUICKDouble gbz);
#endif