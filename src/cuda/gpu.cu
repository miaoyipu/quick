/*
 *  gpu_startup.cu
 *  new_quick
 *
 *  Created by Yipu Miao on 4/20/11.
 *  Copyright 2011 University of Florida. All rights reserved.
 *
 */

#include <stdio.h>

#include "gpu.h"

//-----------------------------------------------
// Set up specified device and be ready to ignite
//-----------------------------------------------
extern "C" void gpu_set_device_(int* gpu_dev_id)
{
    gpu->gpu_dev_id = *gpu_dev_id;
}

//-----------------------------------------------
// create gpu class
//-----------------------------------------------
extern "C" void gpu_startup_(void)
{
	PRINTDEBUG("BEGIN TO WARM UP")
#ifdef DEBUG
    debugFile = fopen("DEBUG", "w+");
#endif
    gpu = new gpu_type;
	PRINTDEBUG("CREATE NEW GPU")
}


//-----------------------------------------------
// Initialize the device
//-----------------------------------------------
extern "C" void gpu_init_(void)
{

    PRINTDEBUG("BEGIN TO INIT")

    int device = -1;
    int gpuCount = 0;
    cudaError_t status;
    cudaDeviceProp deviceProp;
    status = cudaGetDeviceCount(&gpuCount);
    PRINTERROR(status,"cudaGetDeviceCount gpu_init failed!");
    if (gpuCount == 0)
    {
        printf("NO CUDA-Enabled GPU FOUND.\n");
        cudaThreadExit();
        exit(-1);
    }
    
    cudaGetDeviceProperties(&deviceProp, gpu->gpu_dev_id);
//    if ( (deviceProp.major >=2) || ((deviceProp.major == 1) && (deviceProp.minor == 3))) 
        device = gpu->gpu_dev_id;
/*    else {
        printf("SELECT GPU HAS CUDA SUPPORTING VERSION UNDER 1.3. EXITING. \n");
        cudaThreadExit;
        exit(-1);
    }
  */ 
    if (device == -1) {
        printf("NO CUDA 1.3 SUPPORTED GPU IS FOUND\n");
        gpu_shutdown_();
        exit(-1);
    }

    status = cudaSetDevice(device);
    PRINTERROR(status, "cudaSetDevice gpu_init failed!");
    cudaThreadSynchronize();
    
    gpu->blocks = deviceProp.multiProcessorCount;
    if (deviceProp.major ==1) {
        switch (deviceProp.minor) {
/*
            case 0:
            case 1:
            case 2:
            case 5:
                printf("GPU SM VERSION SHOULD BE HIGHER THAN 1.3\n");
                gpu_shutdown_();
                exit(-1);
                break;
*/
            default:
                gpu -> sm_version         =   SM_13;
                gpu -> threadsPerBlock    =   SM_13_THREADS_PER_BLOCK;
                break;
        }
    }else {
        gpu -> sm_version               = SM_2X;
        gpu -> threadsPerBlock          = SM_2X_THREADS_PER_BLOCK;
    }

    PRINTDEBUG("FINISH INIT")

    return;
}

//-----------------------------------------------
// shutdonw gpu and terminate gpu calculation part
//-----------------------------------------------
extern "C" void gpu_shutdown_(void)
{
	PRINTDEBUG("BEGIN TO SHUTDOWN")
#ifdef DEBUG
    fclose(debugFile);
#endif
    delete gpu;
    cudaThreadExit();
	PRINTDEBUG("SHUTDOWN NORMALLY")
    return;
}

//-----------------------------------------------
//  Setup up basic infomation of the system
//-----------------------------------------------
extern "C" void gpu_setup_(int* natom, int* nbasis, int* nElec, int* imult, int* molchg, int* iAtomType)
{

    PRINTDEBUG("BEGIN TO SETUP")

    gpu -> natom                    =   *natom;
    gpu -> nbasis                   =   *nbasis;
    gpu -> nElec                    =   *nElec;
    gpu -> imult                    =   *imult;
    gpu -> molchg                   =   *molchg;
    gpu -> iAtomType                =   *iAtomType;
    gpu -> gpu_calculated           =   new gpu_calculated_type;
    gpu -> gpu_basis                =   new gpu_basis_type;
    gpu -> gpu_cutoff               =   new gpu_cutoff_type;
    gpu -> gpu_calculated -> natom  =   *natom;
    gpu -> gpu_basis -> natom       =   *natom;
    gpu -> gpu_calculated -> nbasis =   *nbasis;
    gpu -> gpu_basis -> nbasis      =   *nbasis;
    
    gpu -> gpu_sim.natom            =   *natom;
    gpu -> gpu_sim.nbasis           =   *nbasis;
    gpu -> gpu_sim.nElec            =   *nElec;
    gpu -> gpu_sim.imult            =   *imult;
    gpu -> gpu_sim.molchg           =   *molchg;
    gpu -> gpu_sim.iAtomType        =   *iAtomType;
        
    PRINTDEBUG("FINISH SETUP")
    upload_para_to_const();
}

//-----------------------------------------------
//  upload coordinates
//-----------------------------------------------
extern "C" void gpu_upload_xyz_(QUICKDouble* atom_xyz)
{
    PRINTDEBUG("BEGIN TO UPLOAD COORDINATES")
    gpu -> xyz = new cuda_buffer_type<QUICKDouble>(atom_xyz, 3, gpu->natom);
    gpu -> gpu_basis -> xyz = new cuda_buffer_type<QUICKDouble>(atom_xyz, 3, gpu->natom);
    gpu -> gpu_calculated -> distance = new cuda_buffer_type<QUICKDouble>(gpu->natom, gpu->natom);
    gpu -> xyz -> Upload();

	gpu -> gpu_basis -> xyz ->Upload();
    gpu -> gpu_sim.xyz =  gpu -> gpu_basis -> xyz -> _devData;
/*
    for (int i = 0; i< gpu->natom; i++) {
        for (int j = i; j<gpu->natom; j++) {
            QUICKDouble distance = 0;
            for (int k = 0; k<3; k++) {
                    distance += pow(LOC2(gpu->xyz->_hostData, k, i, gpu->natom, gpu->natom)
                                   -LOC2(gpu->xyz->_hostData, k, j, gpu->natom, gpu->natom),2);
            }
            LOC2(gpu->gpu_calculated->distance->_hostData, i, j, gpu->natom, gpu->natom) = distance;
            LOC2(gpu->gpu_calculated->distance->_hostData, j, i, gpu->natom, gpu->natom) = distance;
        }
    }
    gpu->gpu_calculated->distance->Upload();
    gpu->gpu_sim.distance =  gpu->gpu_calculated->distance->_devData;
*/
    PRINTDEBUG("COMPLETE UPLOADING COORDINATES")

}


//-----------------------------------------------
//  upload molecule infomation
//-----------------------------------------------
extern "C" void gpu_upload_atom_and_chg_(int* atom, QUICKDouble* atom_chg)
{

    PRINTDEBUG("BEGIN TO UPLOAD ATOM AND CHARGE")
/*    
    gpu -> iattype = new cuda_buffer_type<int>(atom, gpu->natom);
    gpu -> chg     = new cuda_buffer_type<QUICKDouble>(atom_chg, gpu->natom);
    gpu -> iattype -> Upload();
    gpu -> chg     -> Upload();
*/
    PRINTDEBUG("COMPLETE UPLOADING ATOM AND CHARGE")
}

extern "C" void gpu_upload_cutoff_(QUICKDouble* cutMatrix, QUICKDouble* YCutoff, QUICKDouble* integralCutoff,\
                                   QUICKDouble* cutPrim, QUICKDouble* primLimit)
{
    PRINTDEBUG("BEGIN TO UPLOAD CUTOFF")
    
    gpu -> gpu_cutoff -> natom      = gpu -> natom;
    gpu -> gpu_cutoff -> cutMatrix  = new cuda_buffer_type<QUICKDouble>(cutMatrix, gpu->nshell, gpu->nshell);
    gpu -> gpu_cutoff -> YCutoff    = new cuda_buffer_type<QUICKDouble>(YCutoff, gpu->nshell, gpu->nshell);
    gpu -> gpu_cutoff -> integralCutoff = *integralCutoff;
    gpu -> gpu_cutoff -> cutPrim    = new cuda_buffer_type<QUICKDouble>(cutPrim, gpu->jbasis, gpu->jbasis);
    gpu -> gpu_cutoff -> primLimit  = *primLimit;
    
    gpu -> gpu_cutoff -> cutMatrix  -> Upload();
    gpu -> gpu_cutoff -> YCutoff    -> Upload();
    gpu -> gpu_cutoff -> cutPrim    -> Upload();
    
    gpu -> gpu_sim.cutMatrix        = gpu -> gpu_cutoff -> cutMatrix -> _devData;
    gpu -> gpu_sim.YCutoff          = gpu -> gpu_cutoff -> YCutoff -> _devData;
    gpu -> gpu_sim.cutPrim          = gpu -> gpu_cutoff -> cutPrim -> _devData;
    gpu -> gpu_sim.integralCutoff   = gpu -> gpu_cutoff -> integralCutoff;
    gpu -> gpu_sim.primLimit        = gpu -> gpu_cutoff -> primLimit;
    PRINTDEBUG("COMPLETE UPLOADING CUTOFF")
}

//-----------------------------------------------
//  upload calculated information,
//  o is the operator matrix,
//  co is the coeffecient matrix,
//  vec is the eigenvector matrix
//  dense it the density matrix.
//-----------------------------------------------
extern "C" void gpu_upload_calculated_(QUICKDouble* o, QUICKDouble* co, QUICKDouble* vec, QUICKDouble* dense)
{
    PRINTDEBUG("BEGIN TO UPLOAD O MATRIX")
    
    gpu -> gpu_calculated -> o        =   new cuda_buffer_type<QUICKDouble>(o,      gpu->nbasis, gpu->nbasis);
//    gpu -> gpu_calculated -> co       =   new cuda_buffer_type<QUICKDouble>(co,     gpu->nbasis, gpu->nbasis);
//    gpu -> gpu_calculated -> vec      =   new cuda_buffer_type<QUICKDouble>(vec,    gpu->nbasis, gpu->nbasis);
    gpu -> gpu_calculated -> dense    =   new cuda_buffer_type<QUICKDouble>(dense,  gpu->nbasis, gpu->nbasis);
    gpu -> gpu_calculated -> oULL     =   new cuda_buffer_type<QUICKULL>(gpu->nbasis, gpu->nbasis);
    gpu -> gpu_calculated -> o        -> Upload();
//    gpu -> gpu_calculated -> co       -> Upload();
//    gpu -> gpu_calculated -> vec      -> Upload();
    gpu -> gpu_calculated -> dense    -> Upload();
    
    for (int i = 0; i<gpu->nbasis; i++) {
        for (int j = 0; j<gpu->nbasis; j++) {
            QUICKULL valUII = (QUICKULL) (fabs ( LOC2( gpu->gpu_calculated->o->_hostData, i, j, gpu->nbasis, gpu->nbasis)*OSCALE + \
                                                 (QUICKDouble)0.5));
            if (LOC2( gpu->gpu_calculated->o->_hostData, i, j, gpu->nbasis, gpu->nbasis)<(QUICKDouble)0.0)
            valUII = 0ull - valUII;
            LOC2( gpu->gpu_calculated->oULL->_hostData, i, j, gpu->nbasis, gpu->nbasis) = valUII;
            
        }
    }
    gpu -> gpu_calculated -> oULL     -> Upload();
    
    gpu -> gpu_sim.o                 =  gpu -> gpu_calculated -> o -> _devData;
//    gpu -> gpu_sim.co                =  gpu -> gpu_calculated -> co -> _devData;
//    gpu -> gpu_sim.vec               =  gpu -> gpu_calculated -> vec -> _devData;
    gpu -> gpu_sim.dense             =  gpu -> gpu_calculated -> dense -> _devData;
    gpu -> gpu_sim.oULL              =  gpu -> gpu_calculated -> oULL -> _devData;
    
    PRINTDEBUG("COMPLETE UPLOADING O MATRIX")
}

/*
//-----------------------------------------------
//  upload calculated information,
//  o is the operator matrix,
//  co is the coeffecient matrix,
//  vec is the eigenvector matrix
//  dense it the density matrix.
//-----------------------------------------------
extern "C" void gpu_download_o_matrix_(QUICKDouble* o)
{
    PRINTDEBUG("BEGIN TO UPLOAD O MATRIX")
    
    gpu -> gpu_calculated -> o        -> Download(o);
    
    PRINTDEBUG("COMPLETE UPLOADING O MATRIX")
}
*/


//-----------------------------------------------
//  upload basis set information
//-----------------------------------------------
extern "C" void gpu_upload_basis_(int* nshell, int* nprim, int* jshell, int* jbasis, int* maxcontract, \
int* ncontract, int* itype,     QUICKDouble* aexp,      QUICKDouble* dcoeff,\
int* first_basis_function, int* last_basis_function, int* first_shell_basis_function, int* last_shell_basis_function, \
int* ncenter,   int* kstart,    int* katom,     int* ktype,     int* kprim,  int* kshell, int* Ksumtype, \
int* Qnumber,   int* Qstart,    int* Qfinal,    int* Qsbasis,   int* Qfbasis,\
QUICKDouble* gccoeff,           QUICKDouble* cons,      QUICKDouble* gcexpo, int* KLMN)
{
    PRINTDEBUG("BEGIN TO UPLOAD BASIS")
    
    gpu -> gpu_basis -> nshell          =   *nshell;
    gpu -> gpu_basis -> nprim           =   *nprim;
    gpu -> gpu_basis -> jshell          =   *jshell;
    gpu -> gpu_basis -> jbasis          =   *jbasis;
    gpu -> gpu_basis -> maxcontract     =   *maxcontract;
    
    gpu -> nshell                       =   *nshell;
    gpu -> nprim                        =   *nprim;
    gpu -> jshell                       =   *jshell;
    gpu -> jbasis                       =   *jbasis;

    gpu -> gpu_sim.nshell                   =   *nshell;
    gpu -> gpu_sim.nprim                    =   *nprim;
    gpu -> gpu_sim.jshell                   =   *jshell;
    gpu -> gpu_sim.jbasis                   =   *jbasis;
    gpu -> gpu_sim.maxcontract              =   *maxcontract;

    gpu -> gpu_basis -> ncontract       =   new cuda_buffer_type<int>(ncontract, gpu->nbasis);
    gpu -> gpu_basis -> itype           =   new cuda_buffer_type<int>(itype, 3, gpu->nbasis);
    gpu -> gpu_basis -> aexp            =   new cuda_buffer_type<QUICKDouble>(aexp, gpu->gpu_basis->maxcontract, gpu->nbasis);
    gpu -> gpu_basis -> dcoeff          =   new cuda_buffer_type<QUICKDouble>(dcoeff, gpu->gpu_basis->maxcontract, gpu->nbasis);

    gpu -> gpu_basis -> first_basis_function        =   new cuda_buffer_type<int>(first_basis_function, gpu->natom);
    gpu -> gpu_basis -> last_basis_function         =   new cuda_buffer_type<int>(last_basis_function,  gpu->natom);

    gpu -> gpu_basis -> first_shell_basis_function  =   new cuda_buffer_type<int>(first_shell_basis_function, gpu->gpu_basis->nshell);
    gpu -> gpu_basis -> last_shell_basis_function   =   new cuda_buffer_type<int>(last_shell_basis_function,  gpu->gpu_basis->nshell);
    gpu -> gpu_basis -> ncenter                     =   new cuda_buffer_type<int>(ncenter,                    gpu->gpu_basis->nbasis);

    gpu -> gpu_basis -> kstart                      =   new cuda_buffer_type<int>(kstart,   gpu->gpu_basis->nshell);
    gpu -> gpu_basis -> katom                       =   new cuda_buffer_type<int>(katom,    gpu->gpu_basis->nshell);
    gpu -> gpu_basis -> ktype                       =   new cuda_buffer_type<int>(ktype,    gpu->gpu_basis->nshell);
    gpu -> gpu_basis -> kprim                       =   new cuda_buffer_type<int>(kprim,    gpu->gpu_basis->nshell);
    gpu -> gpu_basis -> kshell                      =   new cuda_buffer_type<int>(kshell,   93);
    gpu -> gpu_basis -> Ksumtype                    =   new cuda_buffer_type<int>(Ksumtype, gpu->gpu_basis->nshell+1);

    gpu -> gpu_basis -> Qnumber                     =   new cuda_buffer_type<int>(Qnumber,  gpu->gpu_basis->nshell);
    gpu -> gpu_basis -> Qstart                      =   new cuda_buffer_type<int>(Qstart,   gpu->gpu_basis->nshell);
    gpu -> gpu_basis -> Qfinal                      =   new cuda_buffer_type<int>(Qfinal,   gpu->gpu_basis->nshell);
    gpu -> gpu_basis -> Qsbasis                     =   new cuda_buffer_type<int>(Qsbasis,  gpu->gpu_basis->nshell, 4);
    gpu -> gpu_basis -> Qfbasis                     =   new cuda_buffer_type<int>(Qfbasis,  gpu->gpu_basis->nshell, 4);
    gpu -> gpu_basis -> gccoeff                     =   new cuda_buffer_type<QUICKDouble>(gccoeff, 6, gpu->nbasis);

    gpu -> gpu_basis -> cons                        =   new cuda_buffer_type<QUICKDouble>(cons, gpu->nbasis);
    gpu -> gpu_basis -> gcexpo                      =   new cuda_buffer_type<QUICKDouble>(gcexpo, 6, gpu->nbasis);
    gpu -> gpu_basis -> KLMN                        =   new cuda_buffer_type<int>(KLMN, 3, gpu->nbasis);
    
    gpu -> gpu_basis -> Xcoeff                      =   new cuda_buffer_type<QUICKDouble>(4*gpu->jbasis,4*gpu->jbasis);
    
    
    for (int i = 0; i<gpu->jshell; i++) {
        for (int j = 0; j<gpu->jshell; j++) {
            int kAtomI = gpu->gpu_basis->katom->_hostData[i];
            int kAtomJ = gpu->gpu_basis->katom->_hostData[j];
            int KsumtypeI = gpu->gpu_basis->Ksumtype->_hostData[i];
            int KsumtypeJ = gpu->gpu_basis->Ksumtype->_hostData[j];
            int kstartI = gpu->gpu_basis->kstart->_hostData[i];
            int kstartJ = gpu->gpu_basis->kstart->_hostData[j];
            
            QUICKDouble distance = 0;
            for (int k = 0; k<3; k++) {
                    distance += pow(LOC2(gpu->xyz->_hostData, k, kAtomI-1, gpu->natom, gpu->natom)
                                   -LOC2(gpu->xyz->_hostData, k, kAtomJ-1, gpu->natom, gpu->natom),2);
            }
            
            QUICKDouble DIJ = distance;
            
            for (int ii = 0; ii<gpu->gpu_basis->kprim->_hostData[i]; ii++) {
                for (int jj = 0; jj<gpu->gpu_basis->kprim->_hostData[j]; jj++) {
                    
                    QUICKDouble II = LOC2(gpu->gpu_basis->gcexpo->_hostData, ii , KsumtypeI-1, 6, gpu->nbasis);
                    QUICKDouble JJ = LOC2(gpu->gpu_basis->gcexpo->_hostData, jj , KsumtypeJ-1, 6, gpu->nbasis);
                    
                    
                    QUICKDouble X = exp(-II*JJ/(II+JJ)*DIJ)/(II+JJ);
                    for (int itemp = gpu->gpu_basis->Qstart->_hostData[i]; itemp <= gpu->gpu_basis->Qfinal->_hostData[i]; itemp++) {
                        for (int itemp2 = gpu->gpu_basis->Qstart->_hostData[j]; itemp2 <= gpu->gpu_basis->Qfinal->_hostData[j]; itemp2++) {
                            LOC4(gpu->gpu_basis->Xcoeff->_hostData, kstartI+ii-1, kstartJ+jj-1, itemp, itemp2, gpu->jbasis, gpu->jbasis, 4, 4)
                            = X0 * X * LOC2(gpu->gpu_basis->gccoeff->_hostData, ii, KsumtypeI+itemp-1, 6, gpu->nbasis) \
                                     * LOC2(gpu->gpu_basis->gccoeff->_hostData, jj, KsumtypeJ+itemp2-1, 6, gpu->nbasis);
                        }
                    }
                }
            }
        }
    }
    
    
    gpu -> gpu_basis -> Xcoeff   -> Upload();
    
    gpu -> gpu_basis -> upload_all();

    gpu -> gpu_sim.Xcoeff                       =   gpu -> gpu_basis -> Xcoeff -> _devData;
    gpu -> gpu_sim.ncontract                    =   gpu -> gpu_basis -> ncontract -> _devData;
    gpu -> gpu_sim.first_basis_function         =   gpu -> gpu_basis -> first_basis_function -> _devData;
    gpu -> gpu_sim.last_basis_function          =   gpu -> gpu_basis -> last_basis_function -> _devData;
    gpu -> gpu_sim.first_shell_basis_function   =   gpu -> gpu_basis -> first_shell_basis_function -> _devData;
    gpu -> gpu_sim.last_shell_basis_function    =   gpu -> gpu_basis -> last_shell_basis_function -> _devData;
    gpu -> gpu_sim.ncenter                      =   gpu -> gpu_basis -> ncenter -> _devData;
    gpu -> gpu_sim.kstart                       =   gpu -> gpu_basis -> kstart -> _devData;    
    gpu -> gpu_sim.katom                        =   gpu -> gpu_basis -> katom -> _devData;
    gpu -> gpu_sim.ktype                        =   gpu -> gpu_basis -> ktype -> _devData;
    gpu -> gpu_sim.kprim                        =   gpu -> gpu_basis -> kprim -> _devData;
    gpu -> gpu_sim.kshell                       =   gpu -> gpu_basis -> kshell -> _devData;    
    gpu -> gpu_sim.Ksumtype                     =   gpu -> gpu_basis -> Ksumtype -> _devData;
    gpu -> gpu_sim.Qnumber                      =   gpu -> gpu_basis -> Qnumber -> _devData;
    gpu -> gpu_sim.Qstart                       =   gpu -> gpu_basis -> Qstart -> _devData;
    gpu -> gpu_sim.Qfinal                       =   gpu -> gpu_basis -> Qfinal -> _devData;    
    gpu -> gpu_sim.Qsbasis                      =   gpu -> gpu_basis -> Qsbasis -> _devData;
    gpu -> gpu_sim.Qfbasis                      =   gpu -> gpu_basis -> Qfbasis -> _devData;
    gpu -> gpu_sim.gccoeff                      =   gpu -> gpu_basis -> gccoeff -> _devData;
    gpu -> gpu_sim.cons                         =   gpu -> gpu_basis -> cons -> _devData;
    gpu -> gpu_sim.gcexpo                       =   gpu -> gpu_basis -> gcexpo -> _devData;
    gpu -> gpu_sim.KLMN                         =   gpu -> gpu_basis -> KLMN -> _devData;    

    PRINTDEBUG("COMPLETE UPLOADING BASIS")
}


__global__ void test(QUICKDouble* d);
extern "C" void gpu_get2e_(QUICKDouble* o)
{
    PRINTDEBUG("BEGIN TO RUN GET2E")
    upload_sim_to_constant(gpu);
    
    
    
    PRINTDEBUG("BEGIN TO RUN KERNEL") 
#ifndef TEST
    get2e(gpu);
#else
    for (int i = 1; i<= gpu->gpu_basis->jshell; i++) {
        for (int j = i; j<= gpu->gpu_basis->jshell; j++) {
            for (int k = i; k<= gpu->gpu_basis->jshell; k++) {
                for (int l = k; l<= gpu->gpu_basis->jshell; l++) {
                    gpu_shell(i,j,k,l);
                }
            }
        }
    }
#endif

    PRINTDEBUG("COMPLETE KERNEL")
    gpu -> gpu_calculated -> oULL -> Download();
    
    for (int i = 0; i< gpu->nbasis; i++) {
        for (int j = i; j< gpu->nbasis; j++) {
            QUICKULL valULL = LOC2(gpu->gpu_calculated->oULL->_hostData, j, i, gpu->nbasis, gpu->nbasis);
            QUICKDouble valDB;
            
            if (valULL >= 0x8000000000000000ull) {
                valDB  = -(QUICKDouble)(valULL ^ 0xffffffffffffffffull);
            }
            else
            {
                valDB  = (QUICKDouble) valULL;
            }
            LOC2(gpu->gpu_calculated->o->_hostData,i,j,gpu->nbasis, gpu->nbasis) = (QUICKDouble)valDB*ONEOVEROSCALE;
            LOC2(gpu->gpu_calculated->o->_hostData,j,i,gpu->nbasis, gpu->nbasis) = (QUICKDouble)valDB*ONEOVEROSCALE;
        }
    }
    
    gpu -> gpu_calculated -> o    -> Download(o);

    PRINTDEBUG("DELETE TEMP VARIABLES")
	delete gpu->gpu_cutoff->cutMatrix;
	delete gpu->gpu_cutoff->YCutoff;
	delete gpu->gpu_cutoff->cutPrim;
	
	delete gpu->gpu_calculated->o;
//	delete gpu->gpu_calculated->co;
//	delete gpu->gpu_calculated->vec;
	delete gpu->gpu_calculated->dense;
	delete gpu->gpu_calculated->oULL;


    PRINTDEBUG("COMPLETE RUNNING GET2E")
}
