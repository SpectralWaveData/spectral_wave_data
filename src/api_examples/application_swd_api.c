#include <stdbool.h> /* bool */
#include <stdio.h> /* for fprintf and stderr */
#include <stdlib.h> /* for exit */

#include "spectral_wave_data.h"    // Namespace convention: swd_api_*

int main() {

// First we define some constructor parameters...

// Name of actual SWD-file
char *file_swd = "StokesTest_shallowWater.swd";

// Relate the application and SWD coordinate systems
double x0 = 0.0, y0 = 0.0, t0 = 0.0, beta = 0.0;

// impl=0: apply recommended implementation based on header of swd file.
int impl = 0;

// Number of spectral components to apply in calculations. 
int nsumx = -1, nsumy = -1;   // Negative: apply all components from the swd-file.

// Density of water needed for pressure calculations
double rho = 1025.0;

// A pointer to the object containing the SWD class
void *swd;

// Select temporal interpolation scheme
int ipol = 0;   // C^2 continous

// Select expansion order for kinematics above z=0
int norder = 0;   // Apply same order as specified in SWD file

// Control handling of zero-frequency components in SWD file
bool dc_bias = false;   // Suppress contributions from zero-frequency

// Constructor
swd = swd_api_allocate(file_swd, x0, y0, t0, beta, rho, nsumx, nsumy, 
                       impl, ipol, norder, dc_bias);
if (swd_api_error_raised(swd)) {
    fprintf(stderr, "%s", swd_api_error_get_msg(swd));
    exit(EXIT_FAILURE); /* indicate failure.*/
}

// Time simulation loop...
double t = 0.0, dt = 0.1, tmax = 2.0;
while (t < tmax) {
    // Tell the swd object current application time...
    swd_api_update_time(swd, t);
    if (swd_api_error_raised(swd)) {   // t is too big...
        fprintf(stderr, "%s", swd_api_error_get_msg(swd));
        exit(EXIT_FAILURE); 
    }

    // Application coordinate where we need kinematics...
    double x = 4.3, y = 5.4, z = -1.7;

    // Current set of provided kinematic quantities...

    // Velocity potential
    double phi = swd_api_phi(swd, x, y, z);

    // Stream function
    double varphi =  swd_api_stream(swd, x, y, z);

    // partial(phi) / partial t
    double phi_t = swd_api_phi_t(swd, x, y, z);

    // Particle velocity  (vector_swd is a (x,y,z) struct)
    vector_swd  grad_phi = swd_api_grad_phi(swd, x, y, z);

    // 2nd order gradients of potential  (vector_2nd_phi_swd is a (xx,xy,xz,yy,yz,zz) struct)
    vector_2nd_phi_swd  grad_phi_2nd = swd_api_grad_phi_2nd(swd, x, y, z);

    // Local (Euler) acceleration
    vector_swd  acc_euler = swd_api_acc_euler(swd, x, y, z);

    // Particle acceleration
    vector_swd  acc_particle = swd_api_acc_particle(swd, x, y, z);

    // Wave elevation
    double elev = swd_api_elev(swd, x, y);

    // Local time derivative of wave elevation
    double elev_t = swd_api_elev_t(swd, x, y);

    // Spatial gradient of wave elevation
    vector_swd  grad_elev = swd_api_grad_elev(swd, x, y);

    // 2nd order gradients of elevation  (vector_2nd_elev_swd is a (xx,xy,yy) struct)
    vector_2nd_elev_swd  grad_elev_2nd = swd_api_grad_elev_2nd(swd, x, y);

    // Total pressure
    double wave_pressure = swd_api_pressure(swd, x, y, z);
    
    // Vertical distance from z=0 to sea floor (<0 if infinite)
    double local_depth = swd_api_bathymetry(swd, x, y);

    // Unit normal vector of sea floor into the ocean at(x, y)
    vector_swd  nvec_sf = swd_api_bathymetry_nvec(swd, x, y);

    // For a specific(x, y, z) at this t, return a CSV file on how particle velocity,
    // elevation and pressure converge as a function of number of spectral components
    swd_api_convergence(swd, x, y, z, "dump.csv");
    if (swd_api_error_raised(swd)) { /* in rare cases you could not open the CSV file...*/
        fprintf(stderr, "%s", swd_api_error_get_msg(swd));
        exit(EXIT_FAILURE);
    }

    t += dt;
}

{

// To save storage for an interesting event you may create a new SWD file
// containing only the time steps within a specified time window.
// In this case we only care about the time interval [100.0, 200.0]
swd_api_strip(swd, 100.0, 200.0, "stripped.swd");
if (swd_api_error_raised(swd)) { /* if t=[100.0, 200.0] is not a proper interval...*/
    fprintf(stderr, "%s", swd_api_error_get_msg(swd));
    exit(EXIT_FAILURE);
}

/*
============================================================
There are 4 functions returning a value from the swd file...
   swd_api_get_chr(swd, name)
   swd_api_get_int(swd, name)
   swd_api_get_bool(swd, name)
   swd_api_get_real(swd, name)
where name is the requested parameter in the swd file.
Five examples are given below...
============================================================
*/

// Extract the cid string from the SWD file
// (contains typical the content of the input file applied in the wave generator)
char *ciddata = swd_api_get_chr(swd, "cid");

// Name of actual implementation class
char *cls_name = swd_api_get_chr(swd, "class");

// The nsteps parameter from the swd file
int nsteps = swd_api_get_int(swd, "nsteps");

// The nsteps parameter from the swd file
bool dc_bias = swd_api_get_bool(swd, "dc_bias");

// Time step in SWD file
double dt_swd_file = swd_api_get_real(swd, "dt");

}

// Destructor
swd_api_close(swd);

return 0;
}
