template < int ignoreAdetection > __global__ void MCd(MemStruct);
__device__ float rand_MWC_oc(unsigned long long * , unsigned int * );
__device__ float rand_MWC_co(unsigned long long * , unsigned int * );
__device__ void LaunchPhoton(PhotonStruct * , unsigned long long * , unsigned int * );
__global__ void LaunchPhoton_Global(MemStruct);
__device__ void Spin(PhotonStruct * , float, unsigned long long * , unsigned int * );
__device__ unsigned int Reflect(PhotonStruct * , int, unsigned long long * , unsigned int * );
__device__ unsigned int PhotonSurvive(PhotonStruct * , unsigned long long * , unsigned int * );
__device__ void AtomicAddULL(unsigned long long * address, unsigned int add);

__device__ void LaunchPhoton(PhotonStruct * p, unsigned long long * x, unsigned int * a) {
  p -> x = 0.0f;
  p -> y = 0.0f;
  p -> z = 0.0f;
  p -> dx = 0.0f;
  p -> dy = 0.0f;
  p -> dz = 1.0f;
  p -> layer = 1;
  p -> weight = * start_weight_dc;
}