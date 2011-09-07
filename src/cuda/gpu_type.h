/*
 *  gpu_type.h
 *  new_quick
 *
 *  Created by Yipu Miao on 6/1/11.
 *  Copyright 2011 University of Florida. All rights reserved.
 *
 */

/*
 * this head file includes following type define
    a. common variables             : see details below
    b. gpu type and buffer type     : for communication between GPU and CPU
 */

#include <stdio.h>
#include "../config.h"
#include "gpu_common.h"

// CUDA-C includes
#include <cuda.h>
//#include <cuda_runtime_api.h>

/*
 ****************************************************************
 *  gpu type and buffer type
 ****************************************************************
 */
template <typename T> struct cuda_buffer_type;
struct gpu_calculated_type {
    int                             natom;
    int                             nbasis;
    cuda_buffer_type<QUICKDouble>*  o;      // O matrix
    cuda_buffer_type<QUICKDouble>*  dense;  // Density Matrix
    cuda_buffer_type<QUICKULL>*     oULL;   // Unsigned long long int type O matrix
    
    cuda_buffer_type<QUICKDouble>*  distance;
};

struct gpu_cutoff_type {
    int                             natom;
    int                             nbasis;
    int                             nshell;
    
    int                             sqrQshell;
    cuda_buffer_type<int2>*         sorted_YCutoffIJ;
    cuda_buffer_type<QUICKDouble>*  cutMatrix;
    cuda_buffer_type<QUICKDouble>*  YCutoff;
    cuda_buffer_type<QUICKDouble>*  cutPrim;
    QUICKDouble                     integralCutoff;
    QUICKDouble                     primLimit;
    
};

struct gpu_simulation_type {
    int                             natom;
    int                             nbasis;
    int                             nshell;
    int                             nprim;
    int                             jshell;
    int                             jbasis;
    int                             nElec;
    int                             imult;
    int                             molchg;
    int                             iAtomType;
    int                             maxcontract;
    int                             Qshell;
    // Gaussian Type function
    /*
    int*                            ncontract;
    int*                            itype;
    QUICKDouble*                    aexp;
    QUICKDouble*                    dcoeff;
    */
    // Some more infos about basis function
    QUICKDouble*                    xyz;
/*
    int*                            first_basis_function;
    int*                            last_basis_function;
    int*                            first_shell_basis_function;
    int*                            last_shell_basis_function;
    int*                            ncenter;
  */
    int*                            kstart;
    int*                            katom;
//    int*                            ktype;
    int*                            kprim;
//    int*                            kshell;
    int*                            Ksumtype;
    int*                            Qnumber;
    int*                            Qstart;
    int*                            Qfinal;
    int*                            Qsbasis;
    int*                            Qfbasis;
    int*                            sorted_Qnumber;
    int*                            sorted_Q;
    QUICKDouble*                    gccoeff;
    QUICKDouble*                    cons;
    QUICKDouble*                    gcexpo;
    int*                            KLMN;
    
    // Some more infos about pre-calculated values
    QUICKDouble*                    o;
    QUICKULL*                       oULL;
    QUICKDouble*                    dense;
    
    QUICKDouble*                    distance;
    QUICKDouble*                    Xcoeff;
    QUICKDouble*                    expoSum;
    QUICKDouble*                    weightedCenterX;
    QUICKDouble*                    weightedCenterY;
    QUICKDouble*                    weightedCenterZ;
    
    // cutoff
    int                             sqrQshell;
    int2*                           sorted_YCutoffIJ;
    QUICKDouble*                    cutMatrix;
    QUICKDouble*                    YCutoff;
    QUICKDouble*                    cutPrim;
    QUICKDouble                     integralCutoff;
    QUICKDouble                     primLimit;
};

struct gpu_basis_type {
    int                             natom;
    int                             nbasis;
    int                             nshell;
    int                             nprim;
    int                             jshell;
    int                             jbasis;
    int                             Qshell;
    int                             maxcontract;
    // Gaussian Type function
/*
    cuda_buffer_type<int>*          ncontract;
    cuda_buffer_type<int>*          itype;
    cuda_buffer_type<QUICKDouble>*  aexp;
    cuda_buffer_type<QUICKDouble>*  dcoeff;
  */  
    // Some more infos about basis function
/*
    cuda_buffer_type<QUICKDouble>*  xyz;
    cuda_buffer_type<int>*          first_basis_function;
    cuda_buffer_type<int>*          last_basis_function;
    cuda_buffer_type<int>*          first_shell_basis_function;
    cuda_buffer_type<int>*          last_shell_basis_function;
    cuda_buffer_type<int>*          ncenter;
    cuda_buffer_type<int>*          ktype;
    cuda_buffer_type<int>*          kshell;
  */
    cuda_buffer_type<int>*          kstart;
    cuda_buffer_type<int>*          katom;  
    cuda_buffer_type<int>*          kprim;
    cuda_buffer_type<int>*          Ksumtype;
    cuda_buffer_type<int>*          Qnumber;
    cuda_buffer_type<int>*          Qstart;
    cuda_buffer_type<int>*          Qfinal;
    cuda_buffer_type<int>*          Qsbasis;
    cuda_buffer_type<int>*          Qfbasis;
    cuda_buffer_type<int>*          sorted_Qnumber;
    cuda_buffer_type<int>*          sorted_Q;
    cuda_buffer_type<QUICKDouble>*  gccoeff;
    cuda_buffer_type<QUICKDouble>*  Xcoeff;                     // 4-dimension one
    cuda_buffer_type<QUICKDouble>*  expoSum;                    // 4-dimension one
    cuda_buffer_type<QUICKDouble>*  weightedCenterX;            // 4-dimension one
    cuda_buffer_type<QUICKDouble>*  weightedCenterY;            // 4-dimension one
    cuda_buffer_type<QUICKDouble>*  weightedCenterZ;            // 4-dimension one
    cuda_buffer_type<QUICKDouble>*  cons;
    cuda_buffer_type<QUICKDouble>*  gcexpo;
    cuda_buffer_type<int>*          KLMN;
    cuda_buffer_type<QUICKDouble>*  Apri;
    cuda_buffer_type<QUICKDouble>*  Kpri;
    cuda_buffer_type<QUICKDouble>*  PpriX;
    cuda_buffer_type<QUICKDouble>*  PpriY;
    cuda_buffer_type<QUICKDouble>*  PpriZ;
    
    void upload_all();
    
};



// a type to define a graphic card
struct gpu_type {
    SM_VERSION                      sm_version;
    
    // Memory parameters
    long long int                   totalCPUMemory; // total CPU memory allocated by CUDA part
    long long int                   totalGPUMemory; // total GPU memory allocated by CUDA part
    
    // Launch parameters
    int                             gpu_dev_id;  // set 0 for master GPU
    unsigned int                    blocks;
    unsigned int                    threadsPerBlock;
    
    // Molecule specification part
    int                             natom;
    int                             nbasis;
    int                             nElec;
    int                             imult;
    int                             molchg;
    int                             iAtomType;
    
    int                             nshell;
    int                             nprim;
    int                             jshell;
    int                             jbasis;
    
    cuda_buffer_type<int>*          iattype;
    cuda_buffer_type<QUICKDouble>*  xyz;
    cuda_buffer_type<QUICKDouble>*  chg;

    gpu_calculated_type*            gpu_calculated;
    gpu_basis_type*                 gpu_basis;
    gpu_cutoff_type*                gpu_cutoff;
    gpu_simulation_type             gpu_sim;
    
/*    
    // Method
    cuda_gpu_type();
    ~cuda_gpu_type();
 */
};

typedef struct gpu_type *_gpu_type;
static _gpu_type gpu = NULL;

// template to pack buffered data for GPU-CPU communication
template <typename T>
struct cuda_buffer_type {
    unsigned int    _length;     // length of the data
    unsigned int    _length2;    // length 2 is the row, and if it's not equals to 0, the data is a matrix
    T*              _hostData;   // data on host (CPU)
    T*              _devData;    // data on device (GPU)
    T*              _f90Data;    // if constructed from f90 array, it is the pointer
    
    // constructor
    cuda_buffer_type(int length);   
    cuda_buffer_type(unsigned int length);
    cuda_buffer_type(int length, int length2);
    cuda_buffer_type(unsigned int length, unsigned int length2);
    cuda_buffer_type(T* f90data, int length, int length2);
    cuda_buffer_type(T* f90Data, unsigned int length, unsigned int length2);
    cuda_buffer_type(T* f90data, int length);
    cuda_buffer_type(T* f90Data, unsigned int length);
    
    // destructor
    virtual ~cuda_buffer_type();    
    
    // allocate and deallocate data
    void Allocate();
    void Deallocate();
    
    // use pinned data communication method. Upload and Download from host to device
    void Upload();      
    void Download();
    void Download(T* f90Data);
};


template <typename T>
cuda_buffer_type<T> :: cuda_buffer_type(int length) :
_length(length), _length2(1), _hostData(NULL), _devData(NULL), _f90Data(NULL)
{
    Allocate();
}

template <typename T>
cuda_buffer_type<T> :: cuda_buffer_type(unsigned int length) :
_length(length), _length2(1), _hostData(NULL), _devData(NULL), _f90Data(NULL)
{
    Allocate();
}

template <typename T>
cuda_buffer_type<T> :: cuda_buffer_type(int length, int length2) :
_length(length), _length2(length2), _hostData(NULL), _devData(NULL), _f90Data(NULL)
{
    Allocate();
}

template <typename T>
cuda_buffer_type<T> :: cuda_buffer_type(unsigned int length, unsigned int length2) :
_length(length), _length2(length2), _hostData(NULL), _devData(NULL), _f90Data(NULL)
{
    Allocate();
}

template <typename T>
cuda_buffer_type<T> :: cuda_buffer_type(T* f90data, unsigned int length, unsigned int length2) :
_length(length), _length2(length2), _hostData(NULL), _devData(NULL), _f90Data(f90data)
{
    Allocate();
}

template <typename T>
cuda_buffer_type<T> :: cuda_buffer_type(T* f90data, int length, int length2) :
_length(length), _length2(length2), _hostData(NULL), _devData(NULL), _f90Data(f90data)
{
    Allocate();
}

template <typename T>
cuda_buffer_type<T> :: cuda_buffer_type(T* f90data, unsigned int length) :
_length(length), _length2(1), _hostData(NULL), _devData(NULL), _f90Data(f90data)
{
    Allocate();
}

template <typename T>
cuda_buffer_type<T> :: cuda_buffer_type(T* f90data, int length) :
_length(length), _length2(1), _hostData(NULL), _devData(NULL), _f90Data(f90data)
{
    Allocate();
}


template <typename T>
cuda_buffer_type<T> :: ~cuda_buffer_type()
{
    Deallocate();
}

template <typename T>
void cuda_buffer_type<T> :: Allocate()
{
    
    PRINTDEBUG(">>BEGIN TO ALLOCATE TEMPLATE")

    if (! _f90Data) // if not constructed from f90 array
    {
        cudaError_t status;
        
        //Allocate GPU memeory
        status = cudaMalloc((void**)&_devData,_length*_length2*sizeof(T));
//		status = cudaHostAlloc((void**)&_hostData, _length*_length2*sizeof(T),cudaHostAllocMapped);
        PRINTERROR(status, " cudaMalloc cuda_buffer_type :: Allocate failed!");
        gpu->totalGPUMemory   += _length*_length2*sizeof(T);
        gpu->totalCPUMemory   += _length*_length2*sizeof(T);        
        
//		status = cudaHostGetDevicePointer((void **)&_devData, (void *)_hostData, 0);
//		PRINTERROR(status, " cudaGetDevicePointer cuda_buffer_type :: Allocate failed!");
//		memset(_hostData, 0, _length*_length2*sizeof(T));
        
        //Allocate CPU emembory
  		_hostData = new T[_length*_length2];      
		memset(_hostData, 0, _length*_length2*sizeof(T));
    }else {
        
        // if constructed from f90 array
        cudaError_t status;

        //Allocate GPU memeory
        status = cudaMalloc((void**)&_devData,_length*_length2*sizeof(T));
        PRINTERROR(status, "cudaMalloc cuda_buffertype :: Allocate failed!");
        //Allocate CPU emembory
//        _hostData = (T*) malloc(_length*_length2*sizeof(T));
        _hostData = new T[_length*_length2];
        gpu->totalGPUMemory   += _length*_length2*sizeof(T);
        gpu->totalCPUMemory   += _length*_length2*sizeof(T);        
//     status = cudaHostAlloc((void**)&_hostData, _length*_length2*sizeof(T),cudaHostAllocMapped);
//	   PRINTERROR(status, " cudaMalloc cuda_buffer_type :: Allocate failed!");

//	   status = cudaHostGetDevicePointer((void **)&_devData, (void *)_hostData, 0);
//     PRINTERROR(status, " cudaGetDevicePointer cuda_buffer_type :: Allocate failed!");

        // copy f90 data to _hostData
        size_t index_c = 0;
        size_t index_f = 0;
        for (size_t j=0; j<_length2; j++) {
            for (size_t i=0; i<_length; i++) {
                index_c = i*_length2+j;
                _hostData[index_c] = _f90Data[index_f++];
            }
        }
        
    }
    PRINTMEM("ALLOCATE GPU MEMORY",(unsigned long long int)_length*_length2*sizeof(T))
    PRINTMEM("GPU++",gpu->totalGPUMemory);
    PRINTMEM("CPU++",gpu->totalCPUMemory);
    PRINTDEBUG("<<FINISH ALLOCATION TEMPLATE")
}

template <typename T>
void cuda_buffer_type<T> :: Deallocate()
{

    PRINTDEBUG(">>BEGIN TO DEALLOCATE TEMPLATE")

    cudaError_t status;
    status = cudaFree(_devData);
//	status = cudaFreeHost(_hostData);
    PRINTERROR(status, " cudaFree cuda_buffer_type :: Deallocate failed!");
    
    free(_hostData);
    _hostData = NULL;
    _devData = NULL;
    _f90Data = NULL;
#ifdef DEBUG
    gpu->totalGPUMemory -= _length*_length2*sizeof(T);
    gpu->totalCPUMemory -= _length*_length2*sizeof(T);
    
    PRINTMEM("GPU--",gpu->totalGPUMemory);
    PRINTMEM("CPU--",gpu->totalCPUMemory);
    
#endif
    PRINTDEBUG("<<FINSH DEALLOCATION TEMPLATE")

}

template <typename T>
void cuda_buffer_type<T> :: Upload()
{

    PRINTDEBUG(">>BEGIN TO UPLOAD TEMPLATE")

    cudaError_t status;
    status = cudaMemcpy(_devData,_hostData,_length*_length2*sizeof(T),cudaMemcpyHostToDevice);
    PRINTERROR(status, " cudaMemcpy cuda_buffer_type :: Upload failed!");
    PRINTDEBUG("<<FINISH UPLOADING TEMPLATE")

}

template <typename T>
void cuda_buffer_type<T> :: Download()
{

    PRINTDEBUG(">>BEGIN TO DOWNLOAD TEMPLATE")
    cudaError_t status;
    status = cudaMemcpy(_hostData,_devData,_length*_length2*sizeof(T),cudaMemcpyDeviceToHost);
    PRINTERROR(status, " cudaMemcpy cuda_buffer_type :: Download failed!");
    PRINTDEBUG("<<FINISH DOWNLOADING TEMPLATE")
}

template <typename T>
void cuda_buffer_type<T> :: Download(T* f90Data)
{
    PRINTDEBUG(">>BEGIN TO DOWNLOAD TEMPLATE TO FORTRAN ARRAY")
    size_t index_c = 0;
    size_t index_f;
    for (size_t i = 0; i < _length; i++) {
        for (size_t j = 0; j <_length2; j++) {
            index_f = j*_length+i;
            f90Data[index_f] = _hostData[index_c++];
        }
    }
    PRINTDEBUG("<<FINISH DOWNLOADING TEMPLATE TO FORTRAN ARRAY")
}








