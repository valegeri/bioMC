template < int ignoreAdetection > __global__ void MCd(MemStruct);
__device__ float rand_MWC_oc(unsigned long long * , unsigned int * );
__device__ float rand_MWC_co(unsigned long long * , unsigned int * );
__device__ void LaunchPhoton(PhotonStruct * , unsigned long long * , unsigned int * );
__global__ void LaunchPhoton_Global(MemStruct);
__device__ void Spin(PhotonStruct * , float, unsigned long long * , unsigned int * );
__device__ unsigned int Reflect(PhotonStruct * , int, unsigned long long * , unsigned int * );
__device__ unsigned int PhotonSurvive(PhotonStruct * , unsigned long long * , unsigned int * );
__device__ void AtomicAddULL(unsigned long long * address, unsigned int add);

template < int ignoreAdetection > __global__ void MCd(MemStruct DeviceMem) {
  int bx = blockIdx.x;
  int tx = threadIdx.x;
  int init = NUM_THREADS_PER_BLOCK * bx;

  unsigned long long int x = DeviceMem.x[init + tx]; 
  unsigned int a = DeviceMem.a[init + tx]; 

  float s;

  unsigned int index, w, index_old;
  unsigned int w_temp;
  index_old = 0;
  w = 0;

  PhotonStruct p = DeviceMem.p[init + tx];

  int new_layer;
  unsigned int ii = 0;

  if (!DeviceMem.thread_active[init + tx]) ii = NUMSTEPS_GPU;

  for (; ii < NUMSTEPS_GPU; ii++) {
      
    if (layers_dc[p.layer].mu_total != FLT_MAX)
      s = -__logf(rand_MWC_oc( & x, & a)) * layers_dc[p.layer].mu_total;
    else
      s = 100.0f;

    new_layer = p.layer;

    if (p.z + s * p.dz < layers_dc[p.layer].z_min) {
      new_layer--;
      s = __fdividef(layers_dc[p.layer].z_min - p.z, p.dz);
    }

    if (p.z + s * p.dz > layers_dc[p.layer].z_max) {
      new_layer++;
      s = __fdividef(layers_dc[p.layer].z_max - p.z, p.dz);
    }

    p.x += p.dx * s;
    p.y += p.dy * s;
    p.z += p.dz * s;

    if (p.z > layers_dc[p.layer].z_max) p.z = layers_dc[p.layer].z_max; 
    if (p.z < layers_dc[p.layer].z_min) p.z = layers_dc[p.layer].z_min;

    if (new_layer != p.layer) {
      s = 0.0f;

      if (Reflect( & p, new_layer, & x, & a) == 0u) { 
        if (new_layer == 0) { 
          index = __float2int_rz(acosf(-p.dz) * 2.0f * RPI * det_dc[0].na) * det_dc[0].nr + min(__float2int_rz(__fdividef(sqrtf(p.x * p.x + p.y * p.y), det_dc[0].dr)), (int) det_dc[0].nr - 1);
          AtomicAddULL( & DeviceMem.Rd_ra[index], p.weight);
          p.weight = 0; 
        }
        if (new_layer > * n_layers_dc) { 
          index = __float2int_rz(acosf(p.dz) * 2.0f * RPI * det_dc[0].na) * det_dc[0].nr + min(__float2int_rz(__fdividef(sqrtf(p.x * p.x + p.y * p.y), det_dc[0].dr)), (int) det_dc[0].nr - 1);
          AtomicAddULL( & DeviceMem.Tt_ra[index], p.weight);
          p.weight = 0;
        }
      }
    }

    if (s > 0.0f) {
      w_temp = __float2uint_rn(layers_dc[p.layer].mu * layers_dc[p.layer].mu_total * __uint2float_rn(p.weight));
      p.weight -= w_temp;

      if (ignoreAdetection == 0) {
        index = (min(__float2int_rz(__fdividef(p.z, det_dc[0].dz)), (int) det_dc[0].nz - 1) * det_dc[0].nr + min(__float2int_rz(__fdividef(sqrtf(p.x * p.x + p.y * p.y), det_dc[0].dr)), (int) det_dc[0].nr - 1));
        if (index == index_old) {
          w += w_temp;
        } else {
          AtomicAddULL( & DeviceMem.A_rz[index_old], w);
          index_old = index;
          w = w_temp;
        }
      }
      Spin( & p, layers_dc[p.layer].a_factor, & x, & a);
    }

    if (!PhotonSurvive( & p, & x, & a)) {
      if (atomicAdd(DeviceMem.num_terminated_photons, 1u) < ( * num_photons_dc - NUM_THREADS)) { 
        LaunchPhoton( & p, & x, & a); 
      } else { 
        DeviceMem.thread_active[init + tx] = 0u;
        ii = NUMSTEPS_GPU;
      }
    }
  }

  if (ignoreAdetection == 1 && w != 0)
    AtomicAddULL( & DeviceMem.A_rz[index_old], w);

  __syncthreads();

  DeviceMem.p[init + tx] = p;
  DeviceMem.x[init + tx] = x;
}

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

__global__ void LaunchPhoton_Global(MemStruct DeviceMem) {
  int bx = blockIdx.x;
  int tx = threadIdx.x;
  int init = NUM_THREADS_PER_BLOCK * bx;

  PhotonStruct p;
  unsigned long long int x = DeviceMem.x[init + tx]; 
  unsigned int a = DeviceMem.a[init + tx]; 

  LaunchPhoton( & p, & x, & a);

  DeviceMem.p[init + tx] = p; 
}

__device__ void Spin(PhotonStruct * p, float g, unsigned long long * x, unsigned int * a) {
  float cost, sint; // cosine and sine of the polar deflection angle theta
  float cosp, sinp; // cosine and sine of the azimuthal angle psi. 
  float temp;

  float tempdir = p -> dx;

  temp = __fdividef((1.0f - (g) * (g)), (1.0f - (g) + 2.0f * (g) * rand_MWC_co(x, a)));
  cost = __fdividef((1.0f + (g) * (g) - temp * temp), (2.0f * (g)));
  if (g == 0.0f)
    cost = 2.0f * rand_MWC_co(x, a) - 1.0f;

  sint = sqrtf(1.0f - cost * cost);

  __sincosf(2.0f * PI * rand_MWC_co(x, a), & sinp, & cosp);

  temp = sqrtf(1.0f - p -> dz * p -> dz);

  if (temp == 0.0f) {
    //normal incident.
    p -> dx = sint * cosp;
    p -> dy = sint * sinp;
    p -> dz = copysignf(cost, p -> dz * cost);
  } else {
    // regular incident.
    p -> dx = __fdividef(sint * (p -> dx * p -> dz * cosp - p -> dy * sinp), temp) + p -> dx * cost;
    p -> dy = __fdividef(sint * (p -> dy * p -> dz * cosp + tempdir * sinp), temp) + p -> dy * cost;
    p -> dz = -sint * cosp * temp + p -> dz * cost;
  }

  temp = rsqrtf(p -> dx * p -> dx + p -> dy * p -> dy + p -> dz * p -> dz);
  p -> dx = p -> dx * temp;
  p -> dy = p -> dy * temp;
  p -> dz = p -> dz * temp;
}

__device__ unsigned int Reflect(PhotonStruct * p, int new_layer, unsigned long long * x, unsigned int * a) {
  float n1 = layers_dc[p -> layer].n;
  float n2 = layers_dc[new_layer].n;
  float r;
  float cos_angle_i = fabsf(p -> dz);

  if (n1 == n2) {
    p -> layer = new_layer;
    return 0u;
  }

  if (n1 > n2 && n2 * n2 < n1 * n1 * (1 - cos_angle_i * cos_angle_i)) {
    p -> dz *= -1.0f;
    return 1u;
  }

  if (cos_angle_i == 1.0f) {
    r = __fdividef((n1 - n2), (n1 + n2));
    if (rand_MWC_co(x, a) <= r * r) {
      p -> dz *= -1.0f;
      return 1u;
    } else {
      p -> layer = new_layer;
      return 0u;
    }
  }

  float e = __fdividef(n1 * n1, n2 * n2) * (1.0f - cos_angle_i * cos_angle_i); 
  r = 2 * sqrtf((1.0f - cos_angle_i * cos_angle_i) * (1.0f - e) * e * cos_angle_i * cos_angle_i);
  e = e + (cos_angle_i * cos_angle_i) * (1.0f - 2.0f * e); 
  r = e * __fdividef((1.0f - e - r), ((1.0f - e + r) * (e + r)));

  if (rand_MWC_co(x, a) <= r) {
    p -> dz *= -1.0f;
    return 1u;
  } else {
    r = __fdividef(n1, n2);
    e = r * r * (1.0f - cos_angle_i * cos_angle_i);
    p -> dx *= r;
    p -> dy *= r;
    p -> dz = copysignf(sqrtf(1 - e), p -> dz);
    p -> layer = new_layer;
    return 0u;
  }
}

__device__ unsigned int PhotonSurvive(PhotonStruct * p, unsigned long long * x, unsigned int * a) { 

  if (p -> weight > WEIGHTI) return 1u; 
  if (p -> weight == 0u) return 0u;

  if (rand_MWC_co(x, a) < CHANCE) {
    p -> weight = __float2uint_rn(__fdividef((float) p -> weight, CHANCE));
    return 1u;
  }
  return 0u;
}

__device__ void AtomicAddULL(unsigned long long * address, unsigned int add) {
  if (atomicAdd((unsigned int * ) address, add) + add < add)
    atomicAdd(((unsigned int * ) address) + 1, 1u);
}