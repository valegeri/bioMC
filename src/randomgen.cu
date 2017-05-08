int init_RNG(unsigned long long * x, unsigned int * a,
  const unsigned int n_rng,
  const char * safeprimes_file, unsigned long long xinit) {
  FILE * fp;
  unsigned int begin = 0u;
  unsigned int fora, tmp1, tmp2;

  if (strlen(safeprimes_file) == 0) {
    safeprimes_file = "safeprimes.txt";
  }

  fp = fopen(safeprimes_file, "r");

  if (fp == NULL) {
    printf("Could not find the file of safeprimes (%s)! Terminating!\n", safeprimes_file);
    return 1;
  }

  fscanf(fp, "%u %u %u", & begin, & tmp1, & tmp2);

  if ((xinit == 0ull) | (((unsigned int)(xinit >> 32)) >= (begin - 1)) | (((unsigned int) xinit) >= 0xfffffffful)) {
    printf("%llu not a valid seed! Terminating!\n", xinit);
    return 1;
  }

  for (unsigned int i = 0; i < n_rng; i++) {
    fscanf(fp, "%u %u %u", & fora, & tmp1, & tmp2);
    a[i] = fora;
    x[i] = 0;
    while ((x[i] == 0) | (((unsigned int)(x[i] >> 32)) >= (fora - 1)) | (((unsigned int) x[i]) >= 0xfffffffful)) {
      xinit = (xinit & 0xffffffffull) * (begin) + (xinit >> 32);
      x[i] = (unsigned int) floor((((double)((unsigned int) xinit)) / (double) 0x100000000) * fora);
      x[i] = x[i] << 32;
      xinit = (xinit & 0xffffffffull) * (begin) + (xinit >> 32); 
      x[i] += (unsigned int) xinit;
    }
  }

  fclose(fp);

  return 0;
}

__device__ float rand_MWC_co(unsigned long long * x, unsigned int * a) {
  // Generate a random number [0,1)
  * x = ( * x & 0xffffffffull) * ( * a) + ( * x >> 32);
  return __fdividef(__uint2float_rz((unsigned int)( * x)), (float) 0x100000000);
}

__device__ float rand_MWC_oc(unsigned long long * x, unsigned int * a) {
  // Generate a random number (0,1]
  return 1.0f - rand_MWC_co(x, a);
}
