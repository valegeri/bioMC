#define NUM_BLOCKS 200 					// Depending on the GPU this constant can be changed
#define NUM_THREADS_PER_BLOCK 1024		
#define NUM_THREADS 204800				// NUM_THREADS = NUM_BLOCKS * NUM_THREADS_PER_BLOCK

#define NUMSTEPS_GPU 1000				// Number of steps per one photon
#define PI 3.141592654f
#define RPI 0.318309886f
#define MAX_LAYERS 100
#define STR_LEN 200
#define WEIGHTI 429497u 
#define CHANCE 0.1f


typedef struct __align__(16) {
	float z_min;		// Layer z_min [cm]
	float z_max;		// Layer z_max [cm]
	float mu_total;		// Reciprocal mu_total [cm]
	float mu;			// Absorption coefficient [1/cm]
	float a_factor;		// Anisotropy factor [-]
	float n;			// Refractive index [-]
} LayerStruct;

typedef struct __align__(16) {
	float x;		// Global x coordinate [cm]
	float y;		// Global y coordinate [cm]
	float z;		// Global z coordinate [cm]
	float dx;		// Delta offset, x-direction
	float dy;		// Delta offset, y-direction
	float dz;		// Delta offset, z-direction

	unsigned int weight;	// Photon weight

	int layer;				// Current layer
} PhotonStruct;

typedef struct __align__(16) {
	float dr;		// Detection grid resolution, r-direction [cm]
	float dz;		// Detection grid resolution, z-direction [cm]
	
	int na;			// Number of grid elements in angular-direction [-]
	int nr;			// Number of grid elements in r-direction
	int nz;			// Number of grid elements in z-direction
} DetStruct;


typedef struct {
	unsigned long number_of_photons;
	int ignoreAdetection;
	unsigned int n_layers;
	unsigned int start_weight;
	char outp_filename[STR_LEN];
	char inp_filename[STR_LEN];
	long begin,end;
	char AorB;
	DetStruct det;
	LayerStruct* layers;
} SimulationStruct;


typedef struct {
	PhotonStruct* p;						// Pointer to structure array containing all the photon data
	unsigned long long* x;					// Pointer to the array containing all the WMC x's
	unsigned int* a;						// Pointer to the array containing all the WMC a's
	unsigned int* thread_active;			// Pointer to the array containing the thread active status
	unsigned int* num_terminated_photons;	// Pointer to a scalar keeping track of the number of terminated photons

	unsigned long long* Rd_ra;				// Pointer to the 2D detection matrix used to store reflected photon weights
	unsigned long long* A_rz;				// Pointer to the 2D detection matrix used to store absorbed photon weights
	unsigned long long* Tt_ra;				// Pointer to the 2D detection matrix used to store transmitted photon weights
} MemStruct;
