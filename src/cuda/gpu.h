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
#include "cuda.h"

extern "C" void gpu_set_device_(int* gpu_dev_id);
extern "C" void gpu_startup_(void);
extern "C" void gpu_init_(void);
extern "C" void gpu_shutdown_(void);

extern "C" void gpu_upload_atom_and_chg_(int* atom, QUICKDouble* atom_chg);
extern "C" void gpu_upload_cutoff_(QUICKDouble* cutMatrix, QUICKDouble* YCutoff, QUICKDouble* integralCutoff,\
                                   QUICKDouble* cutPrim, QUICKDouble* primLimit);
extern "C" void gpu_upload_calculated_(QUICKDouble* o, QUICKDouble* co, QUICKDouble* vec, QUICKDouble* dense);
extern "C" void gpu_upload_basis_(int* nshell, int* nprim, int* jshell, int* jbasis, int* maxcontract, \
                                  int* ncontract, int* itype,     QUICKDouble* aexp,      QUICKDouble* dcoeff,\
                                  int* first_basis_function, int* last_basis_function, int* first_shell_basis_function, int* last_shell_basis_function, \
                                  int* ncenter,   int* kstart,    int* katom,     int* ktype,     int* kprim,  int* kshell, int* Ksumtype, \
                                  int* Qnumber,   int* Qstart,    int* Qfinal,    int* Qsbasis,   int* Qfbasis,\
                                  QUICKDouble* gccoeff,           QUICKDouble* cons,      QUICKDouble* gcexpo, int* KLMN);

extern "C" void gpu_get2e_(QUICKDouble* o);
void upload_sim_to_constant(_gpu_type gpu);
void upload_para_to_const();
#ifndef TEST
__global__
#endif
void get2e_kernel();

#ifndef TEST
__device__
#endif
void gpu_shell(unsigned int II, unsigned int JJ, unsigned int KK, unsigned int LL);
#ifndef TEST
__device__
#endif
void FmT(int MaxM, QUICKDouble X, QUICKDouble* vals);
#ifndef TEST
__device__
#endif
void vertical(int NABCDTYPE, QUICKDouble* YVerticalTemp, QUICKDouble* Ptemp, QUICKDouble* WPtemp, \
              QUICKDouble* Qtemp, QUICKDouble* WQtemp,QUICKDouble ABCDtemp,QUICKDouble ABtemp, \
              QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom);
#ifndef TEST
__device__
#endif
void iclass(int I, int J, int K, int L, unsigned int II, unsigned int JJ, unsigned int KK, unsigned int LL, QUICKDouble DNMax);
#ifndef TEST
__device__
#endif
void PSSS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble* Ptemp, QUICKDouble* WPtemp);
#ifndef TEST
__device__
#endif
void SSPS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble* Qtemp, QUICKDouble* WQtemp);
#ifndef TEST
__device__
#endif
void PSPS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble* Ptemp, QUICKDouble* WPtemp, QUICKDouble ABCDtemp);
#ifndef TEST
__device__
#endif
void DSSS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble* Ptemp, QUICKDouble* WPtemp, QUICKDouble ABtemp, QUICKDouble CDcom);

#ifndef TEST
__device__
#endif
void SSDS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble* Qtemp, QUICKDouble* WQtemp, QUICKDouble CDtemp, QUICKDouble ABcom);


#ifndef TEST
__device__
#endif
void DSPS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble* Qtemp, QUICKDouble* WQtemp, QUICKDouble ABCDtemp);

#ifndef TEST
__device__
#endif
void PSDS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble* Ptemp, QUICKDouble* WPtemp, QUICKDouble ABCDtemp);

#ifndef TEST
__device__
#endif
void DSDS(int mtemp, QUICKDouble* YVerticalTemp, QUICKDouble* Ptemp, QUICKDouble* WPtemp, QUICKDouble ABCDtemp, QUICKDouble ABtemp, QUICKDouble CDcom);

#ifndef TEST
__device__
#endif
QUICKDouble hrrwhole(int III, int JJJ, int KKK, int LLL, int IJKLTYPE, QUICKDouble* store, QUICKDouble* RA, QUICKDouble* RB, QUICKDouble* RC, QUICKDouble* RD);



