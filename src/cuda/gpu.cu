/*
 *  gpu_startup.cu
 *  new_quick
 *
 *  Created by Yipu Miao on 4/20/11.
 *  Copyright 2011 University of Florida. All rights reserved.
 *
 */

#include <stdio.h>
#include <string>
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
    
    if (gpu->gpu_dev_id == -1){
        device = 0;
        // if gpu count is greater than 1(multi-gpu) select one with bigger free memory, or available. 
        if (gpuCount > 1) {
            size_t maxMem = 0;
            for (int i = gpuCount-1; i>=0; i--) {
                /*status = cudaSetDevice(i);
                 size_t free_mem = 0;
                 size_t tot_mem  = 0;
                 
                 status = cudaMemGetInfo(&free_mem, &tot_mem); // If error returns, that is to say this device is unavailable.
                 // Else, use one with larger memory.
                 if (free_mem >= maxMem) {
                 maxMem = free_mem;
                 device = i;
                 }
                 cudaThreadExit();*/
                cudaGetDeviceProperties(&deviceProp, i);
                
                if (((deviceProp.major >= 2) || ((deviceProp.major == 1) && (deviceProp.minor == 3))) &&
                    (deviceProp.totalGlobalMem >= maxMem))
                {
                    maxMem                          = deviceProp.totalGlobalMem;
                    device                          = i;
                }
                
            }
    	    gpu->gpu_dev_id = device;
        }       
        
    }else{
        if (gpu->gpu_dev_id >= gpuCount)
        {
            printf("GPU ID IS ILLEGAL, PLEASE SELECT FROM 0 TO %i.\n", gpuCount-1);
            cudaThreadExit();
            exit(-1);
        }
        
    	cudaGetDeviceProperties(&deviceProp, gpu->gpu_dev_id);
    	if ( (deviceProp.major >=2) || ((deviceProp.major == 1) && (deviceProp.minor == 3)))
        	device = gpu->gpu_dev_id;
    	else {
        	printf("SELECT GPU HAS CUDA SUPPORTING VERSION UNDER 1.3. EXITING. \n");
        	cudaThreadExit();
        	exit(-1);
    	}
        device = gpu->gpu_dev_id;
    }
    
    if (device == -1) {
        printf("NO CUDA 1.3 (OR ABOVE) SUPPORTED GPU IS FOUND\n");
        gpu_shutdown_();
        exit(-1);
    }
   
    status = cudaSetDevice(device);
    PRINTERROR(status, "cudaSetDevice gpu_init failed!");
    cudaThreadSynchronize();
    
    gpu->blocks = deviceProp.multiProcessorCount;
    if (deviceProp.major ==1) {
        switch (deviceProp.minor) {
            case 0:
            case 1:
            case 2:
            case 5:
                printf("GPU SM VERSION SHOULD BE HIGHER THAN 1.3\n");
                gpu_shutdown_();
                exit(-1);
                break;
            default:
                gpu -> sm_version           =   SM_13;
                gpu -> threadsPerBlock      =   SM_13_THREADS_PER_BLOCK;
                gpu -> twoEThreadsPerBlock  =   SM_13_2E_THREADS_PER_BLOCK;
                gpu -> XCThreadsPerBlock    =   SM_13_XC_THREADS_PER_BLOCK;
                break;
        }
    }else {
        gpu -> sm_version               = SM_2X;
        gpu -> threadsPerBlock          = SM_2X_THREADS_PER_BLOCK;
        gpu -> twoEThreadsPerBlock      = SM_2X_2E_THREADS_PER_BLOCK;
        gpu -> XCThreadsPerBlock        = SM_2X_XC_THREADS_PER_BLOCK;
    }

    PRINTDEBUG("FINISH INIT")

    return;
}

extern "C" void gpu_get_device_info_(int* gpu_dev_count, int* gpu_dev_id,int* gpu_dev_mem,
                                     int* gpu_num_proc,double* gpu_core_freq,char* gpu_dev_name,int* name_len, int* majorv, int* minorv)
{
    cudaError_t cuda_error;
    cudaDeviceProp prop;
    size_t device_mem;
    
    *gpu_dev_id = gpu->gpu_dev_id;  // currently one single GPU is supported
    cuda_error = cudaGetDeviceCount(gpu_dev_count);
    PRINTERROR(cuda_error,"cudaGetDeviceCount gpu_get_device_info failed!");
    if (*gpu_dev_count == 0) 
    {
        printf("NO CUDA DEVICE FOUNDED \n");
        cudaThreadExit();
        exit(-1);
    }
    cudaGetDeviceProperties(&prop,*gpu_dev_id);
    device_mem = (prop.totalGlobalMem/(1024*1024));
    *gpu_dev_mem = (int) device_mem;
    *gpu_num_proc = (int) (prop.multiProcessorCount);
    *gpu_core_freq = (double) (prop.clockRate * 1e-6f);
    strcpy(gpu_dev_name,prop.name);
    *name_len = strlen(gpu_dev_name);
    *majorv = prop.major;
    *minorv = prop.minor;
    
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
;
//-----------------------------------------------
//  Setup up basic infomation of the system
//-----------------------------------------------
extern "C" void gpu_setup_(int* natom, int* nbasis, int* nElec, int* imult, int* molchg, int* iAtomType)
{

#ifdef DEBUG
    cudaEvent_t start,end;
    cudaEventCreate(&start);
    cudaEventCreate(&end);
    cudaEventRecord(start, 0);
#endif

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

	upload_para_to_const();

#ifdef DEBUG
    cudaEventRecord(end, 0);
    cudaEventSynchronize(end);
    float time;
    cudaEventElapsedTime(&time, start, end);
    PRINTUSINGTIME("UPLOAD PARA TO CONST",time);
    cudaEventDestroy(start);
    cudaEventDestroy(end);
#endif

    PRINTDEBUG("FINISH SETUP")    
}

extern "C" void gpu_upload_method_(int* quick_method)
{
    if (*quick_method == 0) {
        gpu -> gpu_sim.method = HF;
    }else if (*quick_method == 1) {
        gpu -> gpu_sim.method = B3LYP;
    }else if (*quick_method == 2) {
        gpu -> gpu_sim.method = DFT;
    }
}

//-----------------------------------------------
//  upload coordinates
//-----------------------------------------------
extern "C" void gpu_upload_xyz_(QUICKDouble* atom_xyz)
{
#ifdef DEBUG
    cudaEvent_t start,end;
    cudaEventCreate(&start);
    cudaEventCreate(&end);
    cudaEventRecord(start, 0);
#endif

    PRINTDEBUG("BEGIN TO UPLOAD COORDINATES")
//    gpu -> gpu_basis -> xyz = new cuda_buffer_type<QUICKDouble>(atom_xyz, 3, gpu->natom);
//	gpu -> gpu_basis -> xyz ->Upload();
    gpu -> gpu_calculated -> distance = new cuda_buffer_type<QUICKDouble>(gpu->natom, gpu->natom);

    gpu -> xyz = new cuda_buffer_type<QUICKDouble>(atom_xyz, 3, gpu->natom);
    
    for (int i = 0; i < gpu->natom; i++) {
        for (int j = 0; j < gpu->natom; j++) {
            QUICKDouble distance = 0;
            for (int k = 0; k<3; k++) {
                distance += pow(LOC2(gpu->xyz->_hostData, k, i, 3, gpu->natom)
                                -LOC2(gpu->xyz->_hostData, k, j, 3, gpu->natom),2);
            }
            LOC2(gpu->gpu_calculated->distance->_hostData, i, j, gpu->natom, gpu->natom) = sqrt(distance);
        }
    }
    
    gpu -> xyz -> Upload();
    gpu -> gpu_calculated -> distance -> Upload();

    gpu -> gpu_sim.xyz =  gpu -> xyz -> _devData;
    gpu -> gpu_sim.distance = gpu -> gpu_calculated -> distance -> _devData;

#ifdef DEBUG
    cudaEventRecord(end, 0);
    cudaEventSynchronize(end);
    float time;
    cudaEventElapsedTime(&time, start, end);
    PRINTUSINGTIME("UPLOAD XYZ",time);
    cudaEventDestroy(start);
    cudaEventDestroy(end);
#endif

    PRINTDEBUG("COMPLETE UPLOADING COORDINATES")

}


//-----------------------------------------------
//  upload molecule infomation
//-----------------------------------------------
extern "C" void gpu_upload_atom_and_chg_(int* atom, QUICKDouble* atom_chg)
{

    PRINTDEBUG("BEGIN TO UPLOAD ATOM AND CHARGE")
    
    gpu -> iattype = new cuda_buffer_type<int>(atom, gpu->natom);
    gpu -> chg     = new cuda_buffer_type<QUICKDouble>(atom_chg, gpu->natom);
    gpu -> iattype -> Upload();
    gpu -> chg     -> Upload();
    
    
    gpu -> gpu_sim.chg              = gpu -> chg -> _devData;
    gpu -> gpu_sim.iattype          = gpu -> iattype -> _devData;
    
    PRINTDEBUG("COMPLETE UPLOADING ATOM AND CHARGE")
}


//-----------------------------------------------
//  upload cutoff criteria, will update every 
//  interation
//-----------------------------------------------
extern "C" void gpu_upload_cutoff_(QUICKDouble* cutMatrix, QUICKDouble* integralCutoff,QUICKDouble* primLimit, QUICKDouble* DMCutoff)
{

#ifdef DEBUG
    cudaEvent_t start,end;
    cudaEventCreate(&start);
    cudaEventCreate(&end);
    cudaEventRecord(start, 0);
#endif

    PRINTDEBUG("BEGIN TO UPLOAD CUTOFF")
    
    gpu -> gpu_cutoff -> integralCutoff = *integralCutoff;
    gpu -> gpu_cutoff -> primLimit      = *primLimit;
    gpu -> gpu_cutoff -> DMCutoff       = *DMCutoff;
    
    gpu -> gpu_cutoff -> cutMatrix  = new cuda_buffer_type<QUICKDouble>(cutMatrix, gpu->nshell, gpu->nshell);
    
    gpu -> gpu_cutoff -> cutMatrix  -> Upload();

    gpu -> gpu_cutoff -> cutMatrix  -> DeleteCPU();

    gpu -> gpu_sim.cutMatrix        = gpu -> gpu_cutoff -> cutMatrix -> _devData;
    gpu -> gpu_sim.integralCutoff   = gpu -> gpu_cutoff -> integralCutoff;
    gpu -> gpu_sim.primLimit        = gpu -> gpu_cutoff -> primLimit;
    gpu -> gpu_sim.DMCutoff         = gpu -> gpu_cutoff -> DMCutoff;

#ifdef DEBUG
    cudaEventRecord(end, 0);
    cudaEventSynchronize(end);
    float time;
    cudaEventElapsedTime(&time, start, end);
    PRINTUSINGTIME("UPLOAD CUTOFF",time);
    cudaEventDestroy(start);
    cudaEventDestroy(end);
#endif

    PRINTDEBUG("COMPLETE UPLOADING CUTOFF")
}


//-----------------------------------------------
//  upload cutoff matrix, only update at first
//  interation
//-----------------------------------------------
extern "C" void gpu_upload_cutoff_matrix_(QUICKDouble* YCutoff,QUICKDouble* cutPrim)
{

#ifdef DEBUG
    cudaEvent_t start,end;
    cudaEventCreate(&start);
    cudaEventCreate(&end);
    cudaEventRecord(start, 0);
#endif

    PRINTDEBUG("BEGIN TO UPLOAD CUTOFF")
    
    gpu -> gpu_cutoff -> natom      = gpu -> natom;
    gpu -> gpu_cutoff -> YCutoff    = new cuda_buffer_type<QUICKDouble>(YCutoff, gpu->nshell, gpu->nshell);
    gpu -> gpu_cutoff -> cutPrim    = new cuda_buffer_type<QUICKDouble>(cutPrim, gpu->jbasis, gpu->jbasis);
    
    gpu -> gpu_cutoff -> YCutoff    -> Upload();
    gpu -> gpu_cutoff -> cutPrim    -> Upload();
    
    gpu -> gpu_cutoff -> sqrQshell  = (gpu -> gpu_basis -> Qshell) * (gpu -> gpu_basis -> Qshell);
    gpu -> gpu_cutoff -> sorted_YCutoffIJ           = new cuda_buffer_type<int2>(gpu->gpu_cutoff->sqrQshell);
    
    
    int a = 0;
    bool flag = true;
    int2 temp; 
    
    for (int q = 0; q <= 2; q++) {
        for (int p = 0; p <= 2; p++) {
            
            // First to order ERI type
            // Second to order primitive Gaussian function number
            // Third to order Schwartz cutoff upbound
            
            int b=0;
            for (int i = 0; i < gpu->gpu_basis->Qshell; i++) {
                for (int j = 0; j<gpu->gpu_basis->Qshell; j++) {
                    if (gpu->gpu_basis->sorted_Qnumber->_hostData[i] == q && gpu->gpu_basis->sorted_Qnumber->_hostData[j] == p) {
                        if (LOC2(YCutoff, gpu->gpu_basis->sorted_Q->_hostData[i], gpu->gpu_basis->sorted_Q->_hostData[j], gpu->nshell, gpu->nshell) > 1E-9 && 
                            gpu->gpu_basis->sorted_Q->_hostData[i] <= gpu->gpu_basis->sorted_Q->_hostData[j]) {
                            gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[a].x = i;
                            gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[a].y = j;
                            a++;
                            b++;
                        }
                    }
                }
            }    

            PRINTDEBUG("FINISH STEP 1")  
            printf("a=%i b=%i\n", a, b); 
            for (int i = 0; i < b - 1; i ++)
            {
                flag = true;
                for (int j = 0; j < b - i - 1; j ++)
                {
                    if ((LOC2(YCutoff, gpu->gpu_basis->sorted_Q->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+a-b].x], \
                              gpu->gpu_basis->sorted_Q->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+a-b].y], gpu->nshell, gpu->nshell) < \
                         LOC2(YCutoff, gpu->gpu_basis->sorted_Q->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+1+a-b].x], \
                              gpu->gpu_basis->sorted_Q->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+1+a-b].y], gpu->nshell, gpu->nshell)))
                        //&&
                           //gpu->gpu_basis->sorted_Qnumber->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+1].x] == q &&  \
                             //gpu->gpu_basis->sorted_Qnumber->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+1].y]== p &&  \
                             //gpu->gpu_basis->sorted_Qnumber->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j].x] == q && \
                             //gpu->gpu_basis->sorted_Qnumber->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j].y] == p )
                    {
                        temp = gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+a-b];
                        gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+a-b] = gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j + 1+a-b];
                        gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j + 1+a-b] = temp;
                        flag = false;
                    }
                } 
                
                if (flag == true)
                    break;
            }
            
            PRINTDEBUG("FINISH STEP 2")
            flag = true;
             
            for (int i = 0; i < b - 1; i ++)
            {
                flag = true;
                for (int j = 0; j < b - i - 1; j ++)
                { 
                    if (gpu->gpu_basis->kprim->_hostData[gpu->gpu_basis->sorted_Q->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+a-b].x]] *
                        gpu->gpu_basis->kprim->_hostData[gpu->gpu_basis->sorted_Q->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+a-b].y]] <
                        gpu->gpu_basis->kprim->_hostData[gpu->gpu_basis->sorted_Q->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+1+a-b].x]] *
                        gpu->gpu_basis->kprim->_hostData[gpu->gpu_basis->sorted_Q->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+1+a-b].y]])
                    {
                        temp = gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+a-b];
                        gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+a-b] = gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+a-b + 1];
                        gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+a-b + 1] = temp;
                        flag = false;
                    }
                    else if (gpu->gpu_basis->kprim->_hostData[gpu->gpu_basis->sorted_Q->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+a-b].x]] *
                              gpu->gpu_basis->kprim->_hostData[gpu->gpu_basis->sorted_Q->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+a-b].y]] ==
                              gpu->gpu_basis->kprim->_hostData[gpu->gpu_basis->sorted_Q->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+a-b+1].x]] *
                              gpu->gpu_basis->kprim->_hostData[gpu->gpu_basis->sorted_Q->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+a-b+1].y]])
                    {
                        if (gpu->gpu_basis->kprim->_hostData[gpu->gpu_basis->sorted_Q->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+a-b].x]]<
                            gpu->gpu_basis->kprim->_hostData[gpu->gpu_basis->sorted_Q->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+1+a-b].x]]) {
                            temp = gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+a-b];
                            gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j+a-b] = gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j + 1+a-b];
                            gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[j + 1+a-b] = temp;
                            flag = false;
                        }
                    }
                } 
                
                if (flag == true)
                    break;
            }
            
            flag = true;
            PRINTDEBUG("FINISH STEP 3") 
        }
    }
    
    printf("a = %i, total = %i, pect= %f\n", a, gpu->gpu_basis->Qshell * (gpu->gpu_basis->Qshell+1)/2, (float) 2*a/(gpu->gpu_basis->Qshell*(gpu->gpu_basis->Qshell)));
        
    gpu->gpu_cutoff->sqrQshell  = a;
   /* 
    printf("SS = %i\n",a);
    for (int i = 0; i<a; i++) {
        printf("%8i %4i %4i %18.13f Q=%4i %4i %4i %4i prim = %4i %4i\n", i, \
        gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[i].x, \
        gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[i].y, \
        LOC2(YCutoff, gpu->gpu_basis->sorted_Q->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[i].x], gpu->gpu_basis->sorted_Q->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[i].y], gpu->nshell, gpu->nshell),\
        gpu->gpu_basis->sorted_Q->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[i].x], \
        gpu->gpu_basis->sorted_Q->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[i].y], \
        gpu->gpu_basis->sorted_Qnumber->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[i].x], \
        gpu->gpu_basis->sorted_Qnumber->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[i].y], \
        gpu->gpu_basis->kprim->_hostData[gpu->gpu_basis->sorted_Q->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[i].x]],
        gpu->gpu_basis->kprim->_hostData[gpu->gpu_basis->sorted_Q->_hostData[gpu->gpu_cutoff->sorted_YCutoffIJ ->_hostData[i].y]]);
    }
    */
    gpu -> gpu_cutoff -> sorted_YCutoffIJ  -> Upload();
    gpu -> gpu_sim.sqrQshell        = gpu -> gpu_cutoff -> sqrQshell;
    gpu -> gpu_sim.YCutoff          = gpu -> gpu_cutoff -> YCutoff -> _devData;
    gpu -> gpu_sim.cutPrim          = gpu -> gpu_cutoff -> cutPrim -> _devData;
    gpu -> gpu_sim.sorted_YCutoffIJ = gpu -> gpu_cutoff -> sorted_YCutoffIJ  -> _devData;
    
    
    gpu -> gpu_cutoff -> YCutoff -> DeleteCPU();
    gpu -> gpu_cutoff -> cutPrim -> DeleteCPU();
    gpu -> gpu_cutoff -> sorted_YCutoffIJ -> DeleteCPU();
 
#ifdef DEBUG
    cudaEventRecord(end, 0);
    cudaEventSynchronize(end);
    float time;
    cudaEventElapsedTime(&time, start, end);
    PRINTUSINGTIME("UPLOAD CUTOFF",time);
    cudaEventDestroy(start);
    cudaEventDestroy(end);
#endif

    PRINTDEBUG("COMPLETE UPLOADING CUTOFF")
}

//-----------------------------------------------
//  upload calculated information
//-----------------------------------------------
extern "C" void gpu_upload_calculated_(QUICKDouble* o, QUICKDouble* co, QUICKDouble* vec, QUICKDouble* dense)
{

#ifdef DEBUG
    cudaEvent_t start,end;
    cudaEventCreate(&start);
    cudaEventCreate(&end);
    cudaEventRecord(start, 0);
#endif

    PRINTDEBUG("BEGIN TO UPLOAD O MATRIX")
    
    gpu -> gpu_calculated -> o        =   new cuda_buffer_type<QUICKDouble>(o,      gpu->nbasis, gpu->nbasis);
    gpu -> gpu_calculated -> o        ->  DeleteGPU();
    gpu -> gpu_calculated -> dense    =   new cuda_buffer_type<QUICKDouble>(dense,  gpu->nbasis, gpu->nbasis);
    gpu -> gpu_calculated -> oULL     =   new cuda_buffer_type<QUICKULL>(gpu->nbasis, gpu->nbasis);
    
    
    /*
        oULL is the unsigned long long int type of O matrix. The reason to do so is because 
        Atomic Operator for CUDA 2.0 is only available for integer. So for double precision type, 
        an comprimise way is to multiple a very large number (OSCALE), first and divided it
        after atomic operator.
     */
    for (int i = 0; i<gpu->nbasis; i++) {
        for (int j = 0; j<gpu->nbasis; j++) {
            QUICKULL valUII = (QUICKULL) (fabs ( LOC2( gpu->gpu_calculated->o->_hostData, i, j, gpu->nbasis, gpu->nbasis)*OSCALE + (QUICKDouble)0.5));

            if (LOC2( gpu->gpu_calculated->o->_hostData, i, j, gpu->nbasis, gpu->nbasis)<(QUICKDouble)0.0)
            {
                valUII = 0ull - valUII;
            }
            
            LOC2( gpu->gpu_calculated->oULL->_hostData, i, j, gpu->nbasis, gpu->nbasis) = valUII;
        }
    }
    
//    gpu -> gpu_calculated -> o        -> Upload();
    gpu -> gpu_calculated -> dense    -> Upload();
    gpu -> gpu_calculated -> oULL     -> Upload();
    
//    gpu -> gpu_sim.o                 =  gpu -> gpu_calculated -> o -> _devData;
    gpu -> gpu_sim.dense             =  gpu -> gpu_calculated -> dense -> _devData;
    gpu -> gpu_sim.oULL              =  gpu -> gpu_calculated -> oULL -> _devData;
    
    
#ifdef DEBUG
    cudaEventRecord(end, 0);
    cudaEventSynchronize(end);
    float time;
    cudaEventElapsedTime(&time, start, end);
    PRINTUSINGTIME("UPLOAD CALCULATE",time);
    cudaEventDestroy(start);
    cudaEventDestroy(end);
#endif

    PRINTDEBUG("COMPLETE UPLOADING O MATRIX")
}

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

#ifdef DEBUG
    cudaEvent_t start,end;
    cudaEventCreate(&start);
    cudaEventCreate(&end);
    cudaEventRecord(start, 0);
#endif

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


    gpu -> gpu_basis -> ncontract                   =   new cuda_buffer_type<int>(ncontract, gpu->nbasis);//gpu->nbasis);
    gpu -> gpu_basis -> itype                       =   new cuda_buffer_type<int>(itype, 3,  gpu->nbasis);//3, gpu->nbasis);
    gpu -> gpu_basis -> aexp                        =   new cuda_buffer_type<QUICKDouble>(aexp, gpu->gpu_basis->maxcontract, gpu->nbasis);//gpu->gpu_basis->maxcontract, gpu->nbasis);
    gpu -> gpu_basis -> dcoeff                      =   new cuda_buffer_type<QUICKDouble>(dcoeff, gpu->gpu_basis->maxcontract, gpu->nbasis);//gpu->gpu_basis->maxcontract, gpu->nbasis);
/*
    gpu -> gpu_basis -> first_basis_function        =   new cuda_buffer_type<int>(first_basis_function, 1);//gpu->natom);
    gpu -> gpu_basis -> last_basis_function         =   new cuda_buffer_type<int>(last_basis_function,  1);//gpu->natom);

    gpu -> gpu_basis -> first_shell_basis_function  =   new cuda_buffer_type<int>(first_shell_basis_function, 1);//gpu->gpu_basis->nshell);
    gpu -> gpu_basis -> last_shell_basis_function   =   new cuda_buffer_type<int>(last_shell_basis_function,  1);//gpu->gpu_basis->nshell);
 
    gpu -> gpu_basis -> ktype                       =   new cuda_buffer_type<int>(ktype,    gpu->gpu_basis->nshell);
    gpu -> gpu_basis -> kshell                      =   new cuda_buffer_type<int>(kshell,   93);
*/
    gpu -> gpu_basis -> ncenter                     =   new cuda_buffer_type<int>(ncenter,  gpu->gpu_basis->nbasis);

    gpu -> gpu_basis -> kstart                      =   new cuda_buffer_type<int>(kstart,   gpu->gpu_basis->nshell);
    gpu -> gpu_basis -> katom                       =   new cuda_buffer_type<int>(katom,    gpu->gpu_basis->nshell);
    gpu -> gpu_basis -> kprim                       =   new cuda_buffer_type<int>(kprim,    gpu->gpu_basis->nshell);
    gpu -> gpu_basis -> Ksumtype                    =   new cuda_buffer_type<int>(Ksumtype, gpu->gpu_basis->nshell+1);

    gpu -> gpu_basis -> Qnumber                     =   new cuda_buffer_type<int>(Qnumber,  gpu->gpu_basis->nshell);
    gpu -> gpu_basis -> Qstart                      =   new cuda_buffer_type<int>(Qstart,   gpu->gpu_basis->nshell);
    gpu -> gpu_basis -> Qfinal                      =   new cuda_buffer_type<int>(Qfinal,   gpu->gpu_basis->nshell);
    gpu -> gpu_basis -> Qsbasis                     =   new cuda_buffer_type<int>(Qsbasis,  gpu->gpu_basis->nshell, 4);
    gpu -> gpu_basis -> Qfbasis                     =   new cuda_buffer_type<int>(Qfbasis,  gpu->gpu_basis->nshell, 4);
    gpu -> gpu_basis -> gccoeff                     =   new cuda_buffer_type<QUICKDouble>(gccoeff, MAXPRIM, gpu->nbasis);

    gpu -> gpu_basis -> cons                        =   new cuda_buffer_type<QUICKDouble>(cons, gpu->nbasis);
    gpu -> gpu_basis -> gcexpo                      =   new cuda_buffer_type<QUICKDouble>(gcexpo, MAXPRIM, gpu->nbasis);
    gpu -> gpu_basis -> KLMN                        =   new cuda_buffer_type<int>(KLMN, 3, gpu->nbasis);
    
    gpu -> gpu_basis -> prim_start                  =   new cuda_buffer_type<int>(gpu->gpu_basis->nshell);
    gpu -> gpu_basis -> prim_total = 0;
    
    for (int i = 0 ; i < gpu->gpu_basis->nshell; i++) {
        gpu -> gpu_basis -> prim_start -> _hostData[i] = gpu -> gpu_basis -> prim_total;
        gpu -> gpu_basis -> prim_total += gpu -> gpu_basis -> kprim -> _hostData[i];
    }
    
    for (int i = 0; i<gpu->gpu_basis->nshell; i++) {
        printf("for %i prim= %i, start= %i\n", i, gpu -> gpu_basis -> kprim -> _hostData[i], gpu -> gpu_basis -> prim_start -> _hostData[i]);
    }
    printf("total=%i\n", gpu -> gpu_basis -> prim_total);
    int prim_total = gpu -> gpu_basis -> prim_total;
    gpu -> gpu_sim.prim_total = gpu -> gpu_basis -> prim_total;
    
    gpu -> gpu_basis -> Xcoeff                      =   new cuda_buffer_type<QUICKDouble>(2*gpu->jbasis, 2*gpu->jbasis);
    gpu -> gpu_basis -> expoSum                     =   new cuda_buffer_type<QUICKDouble>(prim_total, prim_total);
    gpu -> gpu_basis -> weightedCenterX             =   new cuda_buffer_type<QUICKDouble>(prim_total, prim_total);
    gpu -> gpu_basis -> weightedCenterY             =   new cuda_buffer_type<QUICKDouble>(prim_total, prim_total);
    gpu -> gpu_basis -> weightedCenterZ             =   new cuda_buffer_type<QUICKDouble>(prim_total, prim_total);
    
    
    /*
        After uploading basis set information, we want to do some more things on CPU so that will accelarate GPU.
        The very first is to sort orbital type. In this case, we will calculate s orbitals then p, d, and etc.
        Here Qshell is the number of shell orbtials, for example, sp orbitals account for 2 shell orbitals, and s orbital accounts
        1 shell orbital.
     */
    gpu->gpu_basis->Qshell = 0;
    for (int i = 0; i<gpu->nshell; i++) {
        gpu->gpu_basis->Qshell += gpu->gpu_basis->Qfinal->_hostData[i] - gpu->gpu_basis->Qstart->_hostData[i] + 1;
    }
    
    for (int i = 0; i<gpu->gpu_basis->nshell; i++) {
        for (int j = 0; j<4; j++) {
            LOC2(gpu->gpu_basis->Qsbasis->_hostData, i, j, gpu->gpu_basis->nshell, 4) += gpu->gpu_basis->Ksumtype->_hostData[i];
            LOC2(gpu->gpu_basis->Qfbasis->_hostData, i, j, gpu->gpu_basis->nshell, 4) += gpu->gpu_basis->Ksumtype->_hostData[i];
        }
    }
    
    gpu -> gpu_sim.Qshell = gpu->gpu_basis->Qshell;
    
    gpu -> gpu_basis -> sorted_Q                    =   new cuda_buffer_type<int>( gpu->gpu_basis->Qshell);
    gpu -> gpu_basis -> sorted_Qnumber              =   new cuda_buffer_type<int>( gpu->gpu_basis->Qshell);
    
    /*
        Now because to sort, sorted_Q stands for the shell no, and sorted_Qnumber is the shell orbital type (or angular momentum).
        For instance:
        
        original: s sp s s s sp s s
        sorteed : s s  s s s s  s s p p
        
        move p orbital to the end of the sequence. so the Qshell stands for the length of sequence after sorting.
     */
    int a = 0;
    for (int i = 0; i<gpu->gpu_basis->nshell; i++) {
        for (int j = gpu->gpu_basis->Qstart->_hostData[i]; j<= gpu->gpu_basis->Qfinal->_hostData[i]; j++) {

            if (a == 0) {
                gpu->gpu_basis->sorted_Q->_hostData[0] = i;
                gpu->gpu_basis->sorted_Qnumber->_hostData[0] = j;
            }else {
                for (int k = 0; k<a; k++) {
                    if (j<gpu->gpu_basis->sorted_Qnumber->_hostData[k]) {
                    
                        int kk = k;
                        for (int l = a; l> kk; l--) {
                            gpu->gpu_basis->sorted_Q->_hostData[l] = gpu->gpu_basis->sorted_Q->_hostData[l-1];
                            gpu->gpu_basis->sorted_Qnumber->_hostData[l] = gpu->gpu_basis->sorted_Qnumber->_hostData[l-1];
                        }
                        
                        gpu->gpu_basis->sorted_Q->_hostData[kk] = i;
                        gpu->gpu_basis->sorted_Qnumber->_hostData[kk] = j;
                        break;
                    }
                    gpu->gpu_basis->sorted_Q->_hostData[a] = i;
                    gpu->gpu_basis->sorted_Qnumber->_hostData[a] = j;
                }
            }
            a++;
        }
    }
    
    
    
    /*
    for (int i = 0; i<gpu->gpu_basis->Qshell; i++) {
        for (int j = i; j<gpu->gpu_basis->Qshell; j++) {
            if (gpu->gpu_basis->sorted_Qnumber->_hostData[i] == gpu->gpu_basis->sorted_Qnumber->_hostData[j]) {
                if (gpu->gpu_basis->kprim->_hostData[gpu->gpu_basis->sorted_Q->_hostData[i]] < gpu->gpu_basis->kprim->_hostData[gpu->gpu_basis->sorted_Q->_hostData[j]]) {
                    int temp = gpu->gpu_basis->sorted_Q->_hostData[j];
                    gpu->gpu_basis->sorted_Q->_hostData[j] = gpu->gpu_basis->sorted_Q->_hostData[i];
                    gpu->gpu_basis->sorted_Q->_hostData[i] = temp;
                }
            }
        }
    }*/
    
    printf("Pre-Sorted orbitals:\n");
    printf("Qshell = %i\n", gpu->gpu_basis->Qshell);
    for (int i = 0; i<gpu->gpu_basis->Qshell; i++) {
        printf("i= %i, Q=%i, Qnumber= %i, nprim = %i \n", i, gpu->gpu_basis->sorted_Q->_hostData[i], gpu->gpu_basis->sorted_Qnumber->_hostData[i],
                                                             gpu->gpu_basis->kprim->_hostData[gpu->gpu_basis->sorted_Q->_hostData[i]]);
    }
    
    
    /*
        some pre-calculated variables includes
        
        expoSum(i,j) = expo(i)+expo(j)
        ------------->                 ->          ->
        weightedCenter(i,j) = (expo(i)*i + expo(j)*j)/(expo(i)+expo(j))
     */
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
                    distance += pow(LOC2(gpu->xyz->_hostData, k, kAtomI-1, 3, gpu->natom)
                                   -LOC2(gpu->xyz->_hostData, k, kAtomJ-1, 3, gpu->natom),2);
            }
            
            QUICKDouble DIJ = distance;
            
            for (int ii = 0; ii<gpu->gpu_basis->kprim->_hostData[i]; ii++) {
                for (int jj = 0; jj<gpu->gpu_basis->kprim->_hostData[j]; jj++) {
                    
                    QUICKDouble II = LOC2(gpu->gpu_basis->gcexpo->_hostData, ii , KsumtypeI-1, MAXPRIM, gpu->nbasis);
                    QUICKDouble JJ = LOC2(gpu->gpu_basis->gcexpo->_hostData, jj , KsumtypeJ-1, MAXPRIM, gpu->nbasis);
                    
                    int ii_start = gpu->gpu_basis->prim_start->_hostData[i];
                    int jj_start = gpu->gpu_basis->prim_start->_hostData[j];
                    
                    //expoSum(i,j) = expo(i)+expo(j)
                    LOC2(gpu->gpu_basis->expoSum->_hostData, ii_start+ii, jj_start+jj, prim_total, prim_total) = II + JJ;
                    
                    
                    //        ------------->                 ->          ->
                    //        weightedCenter(i,j) = (expo(i)*i + expo(j)*j)/(expo(i)+expo(j))
                    LOC2(gpu->gpu_basis->weightedCenterX->_hostData, ii_start+ii, jj_start+jj, prim_total, prim_total) = \
                        (LOC2(gpu->xyz->_hostData, 0, kAtomI-1, 3, gpu->natom) * II + LOC2(gpu->xyz->_hostData, 0, kAtomJ-1, 3, gpu->natom)*JJ)/(II+JJ);
                    LOC2(gpu->gpu_basis->weightedCenterY->_hostData, ii_start+ii, jj_start+jj, prim_total, prim_total) = \
                        (LOC2(gpu->xyz->_hostData, 1, kAtomI-1, 3, gpu->natom) * II + LOC2(gpu->xyz->_hostData, 1, kAtomJ-1, 3, gpu->natom)*JJ)/(II+JJ);
                    LOC2(gpu->gpu_basis->weightedCenterZ->_hostData, ii_start+ii, jj_start+jj, prim_total, prim_total) = \
                        (LOC2(gpu->xyz->_hostData, 2, kAtomI-1, 3, gpu->natom) * II + LOC2(gpu->xyz->_hostData, 2, kAtomJ-1, 3, gpu->natom)*JJ)/(II+JJ);
                    
                    
                    // Xcoeff = exp(-II*JJ/(II+JJ) * DIJ) / (II+JJ) * coeff(i) * coeff(j) * X0
                    QUICKDouble X = exp(-II*JJ/(II+JJ)*DIJ)/(II+JJ);
                    
                    for (int itemp = gpu->gpu_basis->Qstart->_hostData[i]; itemp <= gpu->gpu_basis->Qfinal->_hostData[i]; itemp++) {
                        for (int itemp2 = gpu->gpu_basis->Qstart->_hostData[j]; itemp2 <= gpu->gpu_basis->Qfinal->_hostData[j]; itemp2++) {
                            LOC4(gpu->gpu_basis->Xcoeff->_hostData, kstartI+ii-1, kstartJ+jj-1, \
                                 itemp-gpu->gpu_basis->Qstart->_hostData[i], itemp2-gpu->gpu_basis->Qstart->_hostData[j], gpu->jbasis, gpu->jbasis, 2, 2)
                            = X0 * X * LOC2(gpu->gpu_basis->gccoeff->_hostData, ii, KsumtypeI+itemp-1, MAXPRIM, gpu->nbasis) \
                                     * LOC2(gpu->gpu_basis->gccoeff->_hostData, jj, KsumtypeJ+itemp2-1, MAXPRIM, gpu->nbasis);
                        }
                    }
                }
            }
        }
    }
    
    gpu -> gpu_basis -> upload_all();
    
    gpu -> gpu_sim.expoSum                      =   gpu -> gpu_basis -> expoSum -> _devData;
    gpu -> gpu_sim.weightedCenterX              =   gpu -> gpu_basis -> weightedCenterX -> _devData;
    gpu -> gpu_sim.weightedCenterY              =   gpu -> gpu_basis -> weightedCenterY -> _devData;
    gpu -> gpu_sim.weightedCenterZ              =   gpu -> gpu_basis -> weightedCenterZ -> _devData;
    gpu -> gpu_sim.sorted_Q                     =   gpu -> gpu_basis -> sorted_Q -> _devData;
    gpu -> gpu_sim.sorted_Qnumber               =   gpu -> gpu_basis -> sorted_Qnumber -> _devData;
     
    gpu -> gpu_sim.Xcoeff                       =   gpu -> gpu_basis -> Xcoeff -> _devData;

    gpu -> gpu_sim.ncontract                    =   gpu -> gpu_basis -> ncontract -> _devData;
    gpu -> gpu_sim.dcoeff                       =   gpu -> gpu_basis -> dcoeff -> _devData;
    gpu -> gpu_sim.aexp                         =   gpu -> gpu_basis -> aexp -> _devData;
    gpu -> gpu_sim.ncenter                      =   gpu -> gpu_basis -> ncenter -> _devData;
    gpu -> gpu_sim.itype                        =   gpu -> gpu_basis -> itype -> _devData;
    gpu -> gpu_sim.prim_start                   =   gpu -> gpu_basis -> prim_start -> _devData;
/*
    gpu -> gpu_sim.first_basis_function         =   gpu -> gpu_basis -> first_basis_function -> _devData;
    gpu -> gpu_sim.last_basis_function          =   gpu -> gpu_basis -> last_basis_function -> _devData;
    gpu -> gpu_sim.first_shell_basis_function   =   gpu -> gpu_basis -> first_shell_basis_function -> _devData;
    gpu -> gpu_sim.last_shell_basis_function    =   gpu -> gpu_basis -> last_shell_basis_function -> _devData;
    gpu -> gpu_sim.ktype                        =   gpu -> gpu_basis -> ktype -> _devData;
    gpu -> gpu_sim.kshell                       =   gpu -> gpu_basis -> kshell -> _devData;
  */
    gpu -> gpu_sim.kstart                       =   gpu -> gpu_basis -> kstart -> _devData;    
    gpu -> gpu_sim.katom                        =   gpu -> gpu_basis -> katom -> _devData;
    gpu -> gpu_sim.kprim                        =   gpu -> gpu_basis -> kprim -> _devData;    
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


    gpu -> gpu_basis -> expoSum -> DeleteCPU();
    gpu -> gpu_basis -> weightedCenterX -> DeleteCPU();
    gpu -> gpu_basis -> weightedCenterY -> DeleteCPU();
    gpu -> gpu_basis -> weightedCenterZ -> DeleteCPU();
    gpu -> gpu_basis -> Xcoeff -> DeleteCPU();
    
    gpu -> gpu_basis -> ncontract -> DeleteCPU();
    gpu -> gpu_basis -> dcoeff -> DeleteCPU();
    gpu -> gpu_basis -> aexp -> DeleteCPU();
    gpu -> gpu_basis -> ncenter -> DeleteCPU();
    gpu -> gpu_basis -> itype -> DeleteCPU();
    
    gpu -> gpu_basis -> kstart -> DeleteCPU();
    gpu -> gpu_basis -> katom -> DeleteCPU();
    //kprim can not be deleted since it will be used later
    //gpu -> gpu_basis -> kprim -> DeleteCPU();
    gpu -> gpu_basis -> Ksumtype -> DeleteCPU();
    gpu -> gpu_basis -> prim_start -> DeleteCPU();
    
    gpu -> gpu_basis -> Qnumber -> DeleteCPU();
    gpu -> gpu_basis -> Qstart -> DeleteCPU();
    gpu -> gpu_basis -> Qfinal -> DeleteCPU();

    gpu -> gpu_basis -> Qsbasis -> DeleteCPU();
    gpu -> gpu_basis -> Qfbasis -> DeleteCPU();
    gpu -> gpu_basis -> gccoeff -> DeleteCPU();
    gpu -> gpu_basis -> cons -> DeleteCPU();
    gpu -> gpu_basis -> gcexpo -> DeleteCPU();
    gpu -> gpu_basis -> KLMN -> DeleteCPU();


#ifdef DEBUG
    cudaEventRecord(end, 0);
    cudaEventSynchronize(end);
    float time;
    cudaEventElapsedTime(&time, start, end);
    PRINTUSINGTIME("UPLOAD BASIS",time);
    cudaEventDestroy(start);
    cudaEventDestroy(end);
#endif

    PRINTDEBUG("COMPLETE UPLOADING BASIS")
}


//-----------------------------------------------
//  core part, compute 2-e integrals
//-----------------------------------------------
extern "C" void gpu_get2e_(QUICKDouble* o)
{
    PRINTDEBUG("BEGIN TO RUN GET2E")

    upload_sim_to_constant(gpu);

    PRINTDEBUG("BEGIN TO RUN KERNEL") 

    get2e(gpu);

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
   
#ifdef DEBUG
    cudaEvent_t start,end;
    cudaEventCreate(&start);
    cudaEventCreate(&end);
    cudaEventRecord(start, 0);
#endif
 
    gpu -> gpu_calculated -> o    -> Download(o);

#ifdef DEBUG
    cudaEventRecord(end, 0);    
    cudaEventSynchronize(end);
    float time;
    cudaEventElapsedTime(&time, start, end);
    PRINTUSINGTIME("DOWNLOAD O",time);
    cudaEventDestroy(start);
    cudaEventDestroy(end);
#endif 

    PRINTDEBUG("DELETE TEMP VARIABLES")
    
    delete gpu->gpu_calculated->o;
    delete gpu->gpu_calculated->dense;
    delete gpu->gpu_calculated->oULL;

    delete gpu->gpu_cutoff->cutMatrix;

    PRINTDEBUG("COMPLETE RUNNING GET2E")
}

extern "C" void gpu_getxc_(int* isg, QUICKDouble* sigrad2, QUICKDouble* Eelxc, QUICKDouble* aelec, QUICKDouble* belec, QUICKDouble *o)
{
    PRINTDEBUG("BEGIN TO RUN GETXC")

    
    gpu -> gpu_sim.isg = *isg;
    gpu -> gpu_basis -> sigrad2 = new cuda_buffer_type<QUICKDouble>(sigrad2, gpu->nbasis);
    gpu -> gpu_basis -> sigrad2 -> Upload();
    gpu -> gpu_sim.sigrad2      = gpu->gpu_basis->sigrad2->_devData;
    
    gpu -> DFT_calculated       = new cuda_buffer_type<DFT_calculated_type>(1, 1);
    
    
    QUICKULL valUII = (QUICKULL) (fabs ( *Eelxc * OSCALE + (QUICKDouble)0.5));

    if (*Eelxc<(QUICKDouble)0.0)
    {
                valUII = 0ull - valUII;
    }
            
    gpu -> DFT_calculated -> _hostData[0].Eelxc = valUII;
    
    valUII = (QUICKULL) (fabs ( *aelec * OSCALE + (QUICKDouble)0.5));

    if (*aelec<(QUICKDouble)0.0)
    {
                valUII = 0ull - valUII;
    }
    gpu -> DFT_calculated -> _hostData[0].aelec = valUII;
    
    valUII = (QUICKULL) (fabs ( *belec * OSCALE + (QUICKDouble)0.5));

    if (*belec<(QUICKDouble)0.0)
    {
                valUII = 0ull - valUII;
    }
    
    gpu -> DFT_calculated -> _hostData[0].belec = valUII;
    
    gpu -> DFT_calculated -> Upload();
    gpu -> gpu_sim.DFT_calculated= gpu -> DFT_calculated->_devData;
    
    upload_sim_to_constant_dft(gpu);
    PRINTDEBUG("BEGIN TO RUN KERNEL")
    
    getxc(gpu);
    gpu -> gpu_calculated -> oULL -> Download();
    gpu -> DFT_calculated -> Download();
    
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

    
    QUICKULL valULL = gpu->DFT_calculated -> _hostData[0].Eelxc;
    QUICKDouble valDB;
    
    if (valULL >= 0x8000000000000000ull) {
        valDB  = -(QUICKDouble)(valULL ^ 0xffffffffffffffffull);
    }
    else
    {
        valDB  = (QUICKDouble) valULL;
    }
    *Eelxc = (QUICKDouble)valDB*ONEOVEROSCALE;
    
    valULL = gpu->DFT_calculated -> _hostData[0].aelec;
    
    if (valULL >= 0x8000000000000000ull) {
        valDB  = -(QUICKDouble)(valULL ^ 0xffffffffffffffffull);
    }
    else
    {
        valDB  = (QUICKDouble) valULL;
    }
    *aelec = (QUICKDouble)valDB*ONEOVEROSCALE;
    
    valULL = gpu->DFT_calculated -> _hostData[0].belec;
    
    if (valULL >= 0x8000000000000000ull) {
        valDB  = -(QUICKDouble)(valULL ^ 0xffffffffffffffffull);
    }
    else
    {
        valDB  = (QUICKDouble) valULL;
    }
    *belec = (QUICKDouble)valDB*ONEOVEROSCALE;
    
    
    PRINTDEBUG("DELETE TEMP VARIABLES")
    
	delete gpu->gpu_calculated->o;
	delete gpu->gpu_calculated->dense;
	delete gpu->gpu_calculated->oULL;
}



char *trim(char *s) {
    char *ptr;
    if (!s)
        return NULL;   // handle NULL string
    if (!*s)
        return s;      // handle empty string
    for (ptr = s + strlen(s) - 1; (ptr >= s) && isspace(*ptr); --ptr);
    ptr[1] = '\0';
    return s;
}


extern "C" void gpu_aoint_(QUICKDouble* leastIntegralCutoff, QUICKDouble* maxIntegralCutoff, int* intNum, char* intFileName)
{
    PRINTDEBUG("BEGIN TO RUN AOINT")
    
    ERI_entry a;
    FILE *intFile;
    intFile = fopen(trim(intFileName), "wb");
    if (! intFile) {
        printf("UNABLE TO OPEN INT FILE\n");
    }
 	
    int iBatchCount = 0;
    int const availableMem = 400000000;
    int const availableERI = availableMem/sizeof(ERI_entry);
    int nBatchStart[1000], nBatchEnd[1000], nBatchSize[1000];
    int maxIntCount = 0;
    int currentCount = 0;

    
    nBatchStart[0] = 0;
    
    /* 
        fill up the GPU memory and if it is full, run another batch
     */
    
    
    for (int i = 0; i < gpu -> gpu_cutoff -> sqrQshell; i++) {
       
        int intCount = 0;
        
        intCount = (gpu -> gpu_cutoff -> sqrQshell ) * 5;
        
        if (currentCount + intCount < availableERI) {
            currentCount = currentCount + intCount;
        }else{
            nBatchStart[iBatchCount + 1] = i + 1;
            nBatchEnd  [iBatchCount]     = i ;
            nBatchSize [iBatchCount]     = currentCount;
            iBatchCount++;
            currentCount = intCount;
        }
        
    }
    
    
    nBatchEnd[iBatchCount] = gpu -> gpu_cutoff -> sqrQshell - 1  ;
    nBatchSize[iBatchCount]= currentCount;
    iBatchCount++;
    
    for (int i = 0; i < iBatchCount; i++) {
        if (maxIntCount < nBatchSize[i]) {
            maxIntCount = nBatchSize[i];
        }
    }
    
    
    printf("batch count = %i\n", iBatchCount);
    printf("max int count = %i\n", maxIntCount * sizeof(ERI_entry));
    for (int i = 0; i<iBatchCount; i++) {
        printf(" %i from %i to %i %i\n", i, nBatchStart[i], nBatchEnd[i], nBatchSize[i] * sizeof(ERI_entry));
    }
    
    int nBatchERICount = maxIntCount; 
    
    gpu -> aoint_buffer                 = new cuda_buffer_type<ERI_entry>( nBatchERICount, true );
    gpu -> gpu_sim.aoint_buffer         = gpu -> aoint_buffer -> _devData;
    gpu -> gpu_sim.leastIntegralCutoff  = *leastIntegralCutoff;
    gpu -> gpu_sim.maxIntegralCutoff    = *maxIntegralCutoff;
    gpu -> gpu_sim.iBatchSize           = nBatchERICount;
    gpu -> intCount                     = new cuda_buffer_type<QUICKULL>(1);
    gpu -> intCount -> _hostData[0]     = 0;
    gpu -> intCount -> Upload();
    gpu -> gpu_sim.intCount = gpu->intCount->_devData;
    
    upload_sim_to_constant(gpu);
    
#ifdef DEBUG
    float time_downloadERI, time_kernel, time_io;
    time_downloadERI = 0;
    time_io = 0;
    time_kernel = 0;
    cudaEvent_t start_tot,end_tot;
    cudaEventCreate(&start_tot);
    cudaEventCreate(&end_tot);
    cudaEventRecord(start_tot, 0);
#endif

    
    for (int iBatch = 0; iBatch < iBatchCount; iBatch++) {
        
        printf("batch %i start %i end %i\n", iBatch, nBatchStart[iBatch], nBatchEnd[iBatch]);

#ifdef DEBUG
        cudaEvent_t start,end;
        cudaEventCreate(&start);
        cudaEventCreate(&end);
        cudaEventRecord(start, 0);
#endif
        
        gpu -> intCount -> _hostData[0] = 0;
        gpu -> intCount -> Upload();

        // calculate ERI, kernel part
        getAOInt(gpu, nBatchStart[iBatch], nBatchEnd[iBatch]);

#ifdef DEBUG
        cudaEventRecord(end, 0);
        cudaEventSynchronize(end);
        float time;
        cudaEventElapsedTime(&time, start, end);
        PRINTUSINGTIME("KERNEL",time);
        time_kernel += time;
        cudaEventDestroy(start);
        cudaEventDestroy(end);
#endif

        
#ifdef DEBUG
        cudaEventCreate(&start);
        cudaEventCreate(&end);
        cudaEventRecord(start, 0);
#endif
        
        gpu -> intCount -> Download();
        // download ERI from GPU, this is time-consuming part, that need to be reduced
        cudaMemcpy(gpu->aoint_buffer->_hostData, gpu->aoint_buffer->_devData, gpu->intCount->_hostData[0]*sizeof(ERI_entry), cudaMemcpyDeviceToHost);
        
#ifdef DEBUG
        cudaEventRecord(end, 0);
        cudaEventSynchronize(end);
        cudaEventElapsedTime(&time, start, end);
        PRINTUSINGTIME("DOWNLOAD ERI",time);
        time_downloadERI += time;
        cudaEventDestroy(start);
        cudaEventDestroy(end);
#endif
        
        
#ifdef DEBUG
        cudaEventCreate(&start);
        cudaEventCreate(&end);
        cudaEventRecord(start, 0);
#endif
        
        
        printf("intCount = %i\n", gpu->intCount->_hostData[0]);
        // write to disk. there is no way to avoid this part currently.
        for (int i = 0; i < gpu->intCount->_hostData[0]  ; i++) {
            
            a = gpu -> aoint_buffer -> _hostData[i];
            
            if (abs(a.value) > *maxIntegralCutoff) {
                fwrite(&a, sizeof(ERI_entry), 1, intFile);
                *intNum = *intNum + 1;
            }
            
        }
        
        
#ifdef DEBUG
        cudaEventRecord(end, 0);
        cudaEventSynchronize(end);
        cudaEventElapsedTime(&time, start, end);
        PRINTUSINGTIME("IO",time);
        time_io += time;
        cudaEventDestroy(start);
        cudaEventDestroy(end);
#endif
        
        
    }
    
    
    delete gpu->aoint_buffer;

#ifdef DEBUG
    
    cudaEvent_t start,end;
    cudaEventCreate(&start);
    cudaEventCreate(&end);
    cudaEventRecord(start, 0);
#endif
    
    fclose(intFile);

#ifdef DEBUG
    float time;
    cudaEventRecord(end, 0);
    cudaEventSynchronize(end);
    cudaEventElapsedTime(&time, start, end);
    PRINTUSINGTIME("IO FLUSHING",time);
    time_io += time;
    cudaEventDestroy(start);
    cudaEventDestroy(end);
#endif

    
    printf(" TOTAL INT = %i \n", *intNum);
    PRINTDEBUG("END TO RUN AOINT KERNEL")
    
#ifdef DEBUG
    cudaEventRecord(end_tot, 0);
    cudaEventSynchronize(end_tot);
    float time_tot = 0;
    cudaEventElapsedTime(&time_tot, start_tot, end_tot);
    PRINTUSINGTIME("KERNEL",time_kernel);
    PRINTUSINGTIME("DOWNLOAD ERI", time_downloadERI);
    PRINTUSINGTIME("IO", time_io);
    PRINTUSINGTIME("TOTAL",time_tot);
    cudaEventDestroy(start_tot);
    cudaEventDestroy(end_tot);
#endif

}
