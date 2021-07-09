# Monte Carlo simulation of photon propagation in biological tissue

This is CUDA implementation of Monte Carlo simulation of photon propagation in complex multilayerd media (i.e. biological tissue).
The main objective is to obtain the map of the internal photon absorption probability distribution (based on three main matrices with 
reflected, absorbed and transmitted photon weights on each step of propagation).

## General structure

Project consists of 6 CUDA-files:

+ io.cu

    This file contains read-function for input file with initial parameters of medium and number of photons and write-function for output simulation data.

+ structs.h

    SimulationStruct contains all necessary fields that define external properties of one simulation.
    PhotonStruct describes the state of a single photon, i.e. position (x, y, z), direction (dx, dy, dz), weight and current layer.
    DetStruct describes photon detection grid.
    LayerStruct contains fields that describe a single layer.
    MemStruct is the main structure that will be eventuelly stored in global device memory (contains pointer to PhotonStruct).

+ main.cu

    First of all, main.cu file is the entry point of the program.
    The core function in this case is RunSimulation function that does all the tasks: allocates device memory for MemStruct for each thread,
    performs the whole simulation, frees memory and writes output data.
    The main function also reads an input file (via functions from io.cu).

+ memory.cu

    InitMemStructs initializes MemStruct fields in global device memory for each running thread.
    InitDCMem initializes SimulationStruct fields in constant device memory (those fields won't be changed that's why "constant").

+ randomgen.cu

    This is an implementation of random number generator (Multiply-With-Carry version) for seeding all the instances independently.
    All the suitable multipliers for MWC-version are in trueprimes.txt. Each sequence of multipliers will be stored in 2 32-bit SM registers.

+ transport.cu 

    In this file the most prominent function is MCd-function. This kernel function copies data from MemStruct to SM registers and performs
    1000 steps of photon migration for each thread. After that, MCd-function copies ... state back to global memory.

    While LaunchPhoton is just a device function that initializes photons states, LaunchPhoton_Global stores these data in global device memory.

## Running program

First of all, NVIDIA drivers and nvcc compiler is the first and central requirements.

Secondly, in structs.h files `NUM_BLOCKS` and `NUM_THREADS_PER_BLOCK` can be changed, but this depends solely on GPU configuration.
Just keep in mind that current implementation uses on average 64 SM's registers (on Ubuntu). For example, in the case of 64k total number of registers per SM,
1024 thread can run simultaneously on one streaming multiprocessor. 
Number of blocks (aka SMs) is purely based on GPU physical configuration.

If a path to nvcc binary is included in PATH environment variable, then:

```bash
export out="your-name-of-executable"
nvcc src/main.cu -o $out
```

Input parameters of media and photon number can be changed in input.txt. Current version contains homogenous media just for simplicity.

Also, some arguments can be included, like `-S` argument that's related to the seed of the initial state of the RNG. 
If the internal photon distribution (a-matrix) is not on the agenda, then `-A` can be used.
For example:

```bash
$out input.txt -S8 -A
```

Finally, output file with necessary data has been generated.