#define NFLOATS 5
#define NINTS 5

#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>

int isnumeric(char a) {
  if (a >= (char) 48 && a <= (char) 57) return 1;
  else return 0;
}

int readfloats(int n_floats, float * temp, FILE * pFile) {
  int ii = 0;
  char mystring[200];

  if (n_floats > NFLOATS) return 0; 

  while (ii <= 0) {
    if (feof(pFile)) return 0; 
    fgets(mystring, 200, pFile);
    memset(temp, 0, NFLOATS * sizeof(float));
    ii = sscanf(mystring, "%f %f %f %f %f", & temp[0], & temp[1], & temp[2], & temp[3], & temp[4]);
    if (ii > n_floats) return 0; 
  }
  return 1;
}

int readints(int n_ints, int * temp, FILE * pFile) {
  int ii = 0;
  char mystring[STR_LEN];

  if (n_ints > NINTS) return 0;

  while (ii <= 0) {
    if (feof(pFile)) return 0; 
    fgets(mystring, STR_LEN, pFile);
    memset(temp, 0, NINTS * sizeof(int));
    ii = sscanf(mystring, "%d %d %d %d %d", & temp[0], & temp[1], & temp[2], & temp[3], & temp[4]);
    if (ii > n_ints) return 0; 
  }
  return 1; 
}

int ischar(char a) {
  if ((a >= (char) 65 && a <= (char) 90) || (a >= (char) 97 && a <= (char) 122)) return 1;
  else return 0;
}

int read_simulation_data(char * filename, SimulationStruct ** simulations, int ignoreAdetection) {
  int i = 0;
  int ii = 0;
  unsigned long number_of_photons;
  unsigned int start_weight;
  int n_simulations = 0;
  int n_layers = 0;
  FILE * pFile;
  char mystring[STR_LEN];
  char str[STR_LEN];
  char AorB;
  float dtot = 0;

  float ftemp[NFLOATS]; 
  int itemp[NINTS];

  pFile = fopen(filename, "r");

  if (pFile == NULL) {
    perror("Error opening file");
    return 0;
  }

  if (!readints(1, itemp, pFile)) {
    perror("Error reading number of runs");
    return 0;
  }

  n_simulations = itemp[0];

  * simulations = (SimulationStruct * ) malloc(sizeof(SimulationStruct) * n_simulations);
  if ( * simulations == NULL) {
    perror("Failed to malloc simulations.\n");
    return 0;
  } 

  for (i = 0; i < n_simulations; i++) {
    strcpy(( * simulations)[i].inp_filename, filename);
    printf("Input filename: %s\n",filename);

    ( * simulations)[i].ignoreAdetection = ignoreAdetection;

    ii = 0;
    while (ii <= 0) {
      ( * simulations)[i].begin = ftell(pFile);
      fgets(mystring, STR_LEN, pFile);
      ii = sscanf(mystring, "%s %c", str, & AorB);
      if (feof(pFile) || ii > 2) {
        perror("Error reading output filename");
        return 0;
      }
      if (ii > 0) ii = ischar(str[0]);
    }

    strcpy(( * simulations)[i].outp_filename, str);
    ( * simulations)[i].AorB = AorB;

    ii = 0;
    while (ii <= 0) {
      fgets(mystring, STR_LEN, pFile);
      number_of_photons = 0;
      ii = sscanf(mystring, "%lu", & number_of_photons);
      if (feof(pFile) || ii > 1) {
        perror("Error reading number of photons");
        return 0;
      }
    }

    ( * simulations)[i].number_of_photons = number_of_photons;

    if (!readfloats(2, ftemp, pFile)) {
      perror("Error reading dr and dz");
      return 0;
    }

    ( * simulations)[i].det.dz = ftemp[0];
    ( * simulations)[i].det.dr = ftemp[1];

    if (!readints(3, itemp, pFile)) {
      perror("Error reading No. of dz, dr and da");
      return 0;
    }

    ( * simulations)[i].det.nz = itemp[0];
    ( * simulations)[i].det.nr = itemp[1];
    ( * simulations)[i].det.na = itemp[2];

    if (!readints(1, itemp, pFile)) {
      perror("Error reading No. of layers");
      return 0;
    }

    printf("No. of layers=%d\n", itemp[0]);
    n_layers = itemp[0];
    ( * simulations)[i].n_layers = itemp[0];

    ( * simulations)[i].layers = (LayerStruct * ) malloc(sizeof(LayerStruct) * (n_layers + 2));
    if (( * simulations)[i].layers == NULL) {
      perror("Failed to malloc layers.\n");
      return 0;
    } 

    if (!readfloats(1, ftemp, pFile)) {
      perror("Error reading upper refractive index");
      return 0;
    }
    printf("Upper refractive index=%f\n", ftemp[0]);
    ( * simulations)[i].layers[0].n = ftemp[0];

    dtot = 0;
    for (ii = 1; ii <= n_layers; ii++) {
      if (!readfloats(5, ftemp, pFile)) {
        perror("Error reading layer data");
        return 0;
      }
      printf("n=%f, mu=%f, mu_total=%f,a_factor=%f, d=%f\n", ftemp[0], ftemp[1], ftemp[2], ftemp[3], ftemp[4]);
      ( * simulations)[i].layers[ii].n = ftemp[0];
      ( * simulations)[i].layers[ii].mu = ftemp[1];
      ( * simulations)[i].layers[ii].a_factor = ftemp[3];
      ( * simulations)[i].layers[ii].z_min = dtot;
      dtot += ftemp[4];
      ( * simulations)[i].layers[ii].z_max = dtot;
      if (ftemp[2] == 0.0f)( * simulations)[i].layers[ii].mu_total = FLT_MAX;
      else( * simulations)[i].layers[ii].mu_total = 1.0f / (ftemp[1] + ftemp[2]);
    } 

    if (!readfloats(1, ftemp, pFile)) {
      perror("Error reading lower refractive index");
      return 0;
    }
    printf("Lower refractive index=%f\n", ftemp[0]);
    ( * simulations)[i].layers[n_layers + 1].n = ftemp[0];
    ( * simulations)[i].end = ftell(pFile);

    double n1 = ( * simulations)[i].layers[0].n;
    double n2 = ( * simulations)[i].layers[1].n;
    double r = (n1 - n2) / (n1 + n2);
    r = r * r;
    start_weight = (unsigned int)((double) 0xffffffff * (1 - r));
    ( * simulations)[i].start_weight = start_weight;
  } 
  return n_simulations;
}