template < int ignoreAdetection > __global__ void MCd(MemStruct);
__device__ float rand_MWC_oc(unsigned long long * , unsigned int * );
__device__ float rand_MWC_co(unsigned long long * , unsigned int * );
__device__ void LaunchPhoton(PhotonStruct * , unsigned long long * , unsigned int * );
__global__ void LaunchPhoton_Global(MemStruct);
__device__ void Spin(PhotonStruct * , float, unsigned long long * , unsigned int * );
__device__ unsigned int Reflect(PhotonStruct * , int, unsigned long long * , unsigned int * );
__device__ unsigned int PhotonSurvive(PhotonStruct * , unsigned long long * , unsigned int * );
__device__ void AtomicAddULL(unsigned long long * address, unsigned int add);