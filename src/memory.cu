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

int CopyDeviceToHostMem(MemStruct * HostMem, MemStruct * DeviceMem, SimulationStruct * sim) {
  int rz_size = sim -> det.nr * sim -> det.nz;
  int ra_size = sim -> det.nr * sim -> det.na;

  cudaMemcpy(HostMem -> A_rz, DeviceMem -> A_rz, rz_size * sizeof(unsigned long long), cudaMemcpyDeviceToHost);
  cudaMemcpy(HostMem -> Rd_ra, DeviceMem -> Rd_ra, ra_size * sizeof(unsigned long long), cudaMemcpyDeviceToHost);
  cudaMemcpy(HostMem -> Tt_ra, DeviceMem -> Tt_ra, ra_size * sizeof(unsigned long long), cudaMemcpyDeviceToHost);
  cudaMemcpy(HostMem -> x, DeviceMem -> x, NUM_THREADS * sizeof(unsigned long long), cudaMemcpyDeviceToHost);

  return 0;
}

int InitDCMem(SimulationStruct * sim) {
  cudaMemcpyToSymbol(det_dc, & (sim -> det), sizeof(DetStruct));
  cudaMemcpyToSymbol(n_layers_dc, & (sim -> n_layers), sizeof(unsigned int));
  cudaMemcpyToSymbol(start_weight_dc, & (sim -> start_weight), sizeof(unsigned int));
  cudaMemcpyToSymbol(layers_dc, sim -> layers, (sim -> n_layers + 2) * sizeof(LayerStruct));
  cudaMemcpyToSymbol(num_photons_dc, & (sim -> number_of_photons), sizeof(unsigned int));

  return 0;
}

void FreeMemStructs(MemStruct * HostMem, MemStruct * DeviceMem) {
  free(HostMem -> A_rz);
  free(HostMem -> Rd_ra);
  free(HostMem -> Tt_ra);
  free(HostMem -> thread_active);
  free(HostMem -> num_terminated_photons);

  cudaFree(DeviceMem -> A_rz);
  cudaFree(DeviceMem -> Rd_ra);
  cudaFree(DeviceMem -> Tt_ra);
  cudaFree(DeviceMem -> x);
  cudaFree(DeviceMem -> a);
  cudaFree(DeviceMem -> thread_active);
  cudaFree(DeviceMem -> num_terminated_photons);
}

void FreeSimulationStruct(SimulationStruct * sim, int n_simulations) {
  for (int i = 0; i < n_simulations; i++) free(sim[i].layers);
  free(sim);
}