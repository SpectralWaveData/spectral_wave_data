#ifndef SWD_API_H_INCLUDED
#define SWD_API_H_INCLUDED
/*
The C interface for the spectral-wave-data ocean wave model.

For all functions we apply the naming convention: swd_api_*

Ver-3.0: 
Coded by: Jens B. Helmers DNVGL,  2019.07.31

*/

#ifdef __cplusplus
extern "C" {
#endif

#ifndef __cplusplus
#include <stdbool.h>
#endif    

// Floating point model for interfacing with spectral_wave_data
#ifdef SWD_API_FLOAT
   typedef float real_swd; 
#else
   typedef double real_swd; 
#endif

typedef struct {
    real_swd x; 
    real_swd y; 
    real_swd z;
} vector_swd;

typedef struct {
    real_swd xx;   // Direction of derivatives
    real_swd xy; 
    real_swd xz; 
    real_swd yy; 
    real_swd yz; 
    real_swd zz;
} vector_2nd_phi_swd;

typedef struct {
    real_swd xx;   // Direction of derivatives
    real_swd xy;
    real_swd yy;
} vector_2nd_elev_swd;


// constructor
void *swd_api_allocate(const char *file_swd, real_swd x0, real_swd y0, 
                       real_swd t0, real_swd beta, real_swd rho, 
                       int nsumx, int nsumy, int impl, int ipol,
                       int norder, bool dc_bias);
    /*
    file_swd:          Name of actual swd file

    x0, y0, t0, beta:  Relation between SWD and application coordinates.
                       beta in degree.

    rho:               Density of water(applied for pressure calculations)

    nsumx, nsumy       Number of spectral components to apply (<0: apply all)

    impl               Index to determine actual derived class to apply
                       0 = Automatic selection of best class (default)
                      <0 = In-house and experimental implementations
                      >0 = Validated implementations available open software

    ipol               Index to request actual temporal interpolation scheme
                       0 = C^2 continous scheme (default)
                       1 = C^1 continous
                       2 = C^3 continous

    norder             Expansion order to apply in kinematics for z>0
                       0 = Apply expansion order specified in swd file (default)
                      <0 = Apply exp(kj z)
                      >0 = Apply expansion order = norder 

    dc_bias            Control application of zero-frequency bias present in SWD file  
                       false = Suppress contribution from zero frequency amplitudes (default)
                       true  = Apply zero frequency amplitudes from SWD file.
                       
    Error signals has been raised if the function 'swd_api_error_raised()' return true.
    That function should be called after this constructor to check if successful creation.
    The actual error codes is returned using the function 'swd_api_error_get_id()':

    1001: SWD file not able to open
    1002: SWD file has wrong binary convention (not float/little endian)
    1003: SWD file data read error
    1004: invalid input parameters
    1005: Not able to allocate memory

   */

// destructor
void swd_api_close(void *swd);


// ===================================================================
//  Field calculations at current application time:
// ===================================================================

// apply current application time
// Eventual error signals (ref constructor): 1003, 1004
void swd_api_update_time(void *swd, real_swd time);

// wave potential
// Eventual error signals (ref constructor): None
real_swd swd_api_phi(void *swd, real_swd x, real_swd y, real_swd z);

// stream function
// Eventual error signals (ref constructor): None
real_swd swd_api_stream(void *swd, real_swd x, real_swd y, real_swd z);

// time derivative of wave potential (Euler derivative, earth fixed)
// Eventual error signals (ref constructor): None
real_swd swd_api_phi_t(void *swd, real_swd x, real_swd y, real_swd z);

// particle velocity
// Eventual error signals (ref constructor): None
vector_swd swd_api_grad_phi(void *swd, real_swd x, real_swd y, real_swd z);

// 2nd order gradients of potential
// Eventual error signals (ref constructor): None
vector_2nd_phi_swd swd_api_grad_phi_2nd(void *swd, real_swd x, real_swd y, real_swd z);

// Local (Euler) acceleration
// Eventual error signals (ref constructor): None
vector_swd swd_api_acc_euler(void *swd, real_swd x, real_swd y, real_swd z);

// Particle acceleration
// Eventual error signals (ref constructor): None
vector_swd swd_api_acc_particle(void *swd, real_swd x, real_swd y, real_swd z);

// wave elevation
// Eventual error signals (ref constructor): None
real_swd swd_api_elev(void *swd, real_swd x, real_swd y);

// Local time derivative of wave elevation
// Eventual error signals (ref constructor): None
real_swd swd_api_elev_t(void *swd, real_swd x, real_swd y);

// Gradient of wave elevation (slopes)
// Eventual error signals (ref constructor): None
vector_swd swd_api_grad_elev(void *swd, real_swd x, real_swd y);

//  2nd order gradients of elevation
// Eventual error signals (ref constructor): None
vector_2nd_elev_swd swd_api_grad_elev_2nd(void *swd, real_swd x, real_swd y);

// Complete Bernoulli pressure
// Eventual error signals (ref constructor): None
real_swd swd_api_pressure(void *swd, real_swd x, real_swd y, real_swd z);

// Vertical distance from z=0 to sea floor (<0 if infinite)
// Eventual error signals (ref constructor): None
real_swd swd_api_bathymetry(void *swd, real_swd x, real_swd y);

// Unit normal vector of sea floor into the ocean at(x, y)
// Eventual error signals (ref constructor): None
vector_swd swd_api_bathymetry_nvec(void *swd, real_swd x, real_swd y);

// For a specific location return a csv-file on how velocity, elevation
// and pressure converge as a function of number of spectral components.
// Eventual error signals (ref constructor): 1001
void swd_api_convergence(void *swd, real_swd x, real_swd y, real_swd z, const char *csv);

// To save storage for an interesting event you may create a new SWD file
// containing only the time steps within the time window [tmin, tmax].
// The name of the new SWD file is defined by file_swd_new
// Eventual error signals (ref constructor): 1001, 1003
void swd_api_strip(void *swd, real_swd tmin, real_swd tmax, const char *file_swd_new);

// ===================================================================
//  Provide parameters from the swd-file:

// Extract the character parameter 'name'
// Eventual error signals (ref constructor): 1004
// **NOTE**: This function is not thread-safe. Make a copy of the
// result if you need to store it before calling this function again.
char *swd_api_get_chr(void *swd, const char *name);

// Extract the int parameter 'name'
// Eventual error signals (ref constructor): 1004
int swd_api_get_int(void *swd, const char *name);

// Extract the bool parameter 'name'
// Eventual error signals (ref constructor): 1004
bool swd_api_get_bool(void *swd, const char *name);

// Extract the real parameter 'name'
// Eventual error signals (ref constructor): 1004
real_swd swd_api_get_real(void *swd, const char *name);

// ===================================================================
//  Provide error handling:

// Check if SWD object has signaled an error
bool swd_api_error_raised(void *swd);

// Extract error id from SWD obhect
int swd_api_error_get_id(void *swd);

// Extract error message from SWD object
char *swd_api_error_get_msg(void *swd);

// Clear the error flag from SWD object
void *swd_api_error_clear(void *swd);

// ===================================================================

#ifdef __cplusplus
}
#endif

#endif