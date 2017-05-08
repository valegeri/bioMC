int InitMemStructs(MemStruct * HostMem, MemStruct * DeviceMem, SimulationStruct * sim) {
  int rz_size, ra_size;
  rz_size = sim -> det.nr * sim -> det.nz;
  ra_size = sim -> det.nr * sim -> det.na;

  cudaMalloc((void ** ) & DeviceMem -> p, NUM_THREADS * sizeof(PhotonStruct));

  HostMem -> A_rz = (unsigned long long * ) malloc(rz_size * sizeof(unsigned long long));

  if (HostMem -> A_rz == NULL) {
    printf("Error allocating HostMem->A_rz");
    exit(1);
  }

  cudaMalloc((void ** ) & DeviceMem -> A_rz, rz_size * sizeof(unsigned long long));
  cudaMemset(DeviceMem -> A_rz, 0, rz_size * sizeof(unsigned long long));

  HostMem -> Rd_ra = (unsigned long long * ) malloc(ra_size * sizeof(unsigned long long));

  if (HostMem -> Rd_ra == NULL) {
    printf("Error allocating HostMem->Rd_ra");
    exit(1);
  }

  cudaMalloc((void ** ) & DeviceMem -> Rd_ra, ra_size * sizeof(unsigned long long));
  cudaMemset(DeviceMem -> Rd_ra, 0, ra_size * sizeof(unsigned long long));

  HostMem -> Tt_ra = (unsigned long long * ) malloc(ra_size * sizeof(unsigned long long));

  if (HostMem -> Tt_ra == NULL) {
    printf("Error allocating HostMem->Tt_ra");
    exit(1);
  }

  cudaMalloc((void ** ) & DeviceMem -> Tt_ra, ra_size * sizeof(unsigned long long));
  cudaMemset(DeviceMem -> Tt_ra, 0, ra_size * sizeof(unsigned long long));

  cudaMalloc((void ** ) & DeviceMem -> x, NUM_THREADS * sizeof(unsigned long long));
  cudaMemcpy(DeviceMem -> x, HostMem -> x, NUM_THREADS * sizeof(unsigned long long), cudaMemcpyHostToDevice);
  cudaMalloc((void ** ) & DeviceMem -> a, NUM_THREADS * sizeof(unsigned int));
  cudaMemcpy(DeviceMem -> a, HostMem -> a, NUM_THREADS * sizeof(unsigned int), cudaMemcpyHostToDevice);

  HostMem -> thread_active = (unsigned int * ) malloc(NUM_THREADS * sizeof(unsigned int));

  if (HostMem -> thread_active == NULL) {
    printf("Error allocating HostMem->thread_active");
    exit(1);
  }

  for (int i = 0; i < NUM_THREADS; i++) HostMem -> thread_active[i] = 1u;

  cudaMalloc((void ** ) & DeviceMem -> thread_active, NUM_THREADS * sizeof(unsigned int));
  cudaMemcpy(DeviceMem -> thread_active, HostMem -> thread_active, NUM_THREADS * sizeof(unsigned int), cudaMemcpyHostToDevice);

  HostMem -> num_terminated_photons = (unsigned int * ) malloc(sizeof(unsigned int));

  if (HostMem -> num_terminated_photons == NULL) {
    printf("Error allocating HostMem->num_terminated_photons");
    exit(1);
  }

  * HostMem -> num_terminated_photons = 0;

  cudaMalloc((void ** ) & DeviceMem -> num_terminated_photons, sizeof(unsigned int));
  cudaMemcpy(DeviceMem -> num_terminated_photons, HostMem -> num_terminated_photons, sizeof(unsigned int), cudaMemcpyHostToDevice);

  return 1;
}
