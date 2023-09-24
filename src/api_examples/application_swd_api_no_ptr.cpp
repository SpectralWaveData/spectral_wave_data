#include <string>
#include <iostream>
#include <stdlib.h> /* for exit */

#include "SpectralWaveData.h"    // Version for C++

int main()
{

// Actual swd file
std::string file_swd = "test_nx3.swd";

// Constructor parameters to relate the application and SWD coordinate systems
double x0 = 0.0, y0 = 0.0, t0 = 0.0, beta = 0.0;

// impl=0: apply recommended implementation based on header of swd file.
int impl = 0;  

// Density of water needed for pressure calculations
double rho = 1025.0;

// Number of spectral components to apply in calculations. 
int nsumx = -1, nsumy = -1;   // Negative: apply all components from swd-file.

// Optional: Select temporal interpolation scheme
int ipol = 0;   // C^2 continous

// Optional: Select expansion order for kinematics above z=0
int norder = 0;   // Apply same order as specified in SWD file

// Optional: Control handling of zero-frequency components in SWD file
bool dc_bias = false;   // Suppress contributions from zero-frequency

//============
// Constructor
//============
SpectralWaveData swd = SpectralWaveData(file_swd, x0, y0, t0, beta,
                           rho, nsumx, nsumy, impl, ipol,
                           norder, dc_bias);

//======================================
// Time domain simulation...
//======================================
double t = 0.0, dt = 0.1, tmax = 2.0;
while (t < tmax) {

    // Tell the swd object current application time...
    try {
        swd.UpdateTime(t);
    } catch (SwdException& e) {  //Could be t > tmax from file.
        std::cerr << typeid(e).name() << std::endl << e.what() << std::endl;
        // If we will try again with a new value of t
        // we first need to call: swd.ExceptionClear()
        exit(EXIT_FAILURE);  // In this case we just abort.
    }

    // Application coordinate where we need kinematics...
    double x = 4.3, y = 5.4, z = -1.7;

    // Current set of provided kinematic quantities...

    // Velocity potential
    double wave_potential = swd.Phi(x, y, z);

    // Stream function
    double stream_function = swd.Stream(x, y, z);

    // partial(phi) / partial t
    double wave_ddt_potential = swd.DdtPhi(x, y, z);

    // Particle velocity  (vector_swd is a (x,y,z) struct)
    vector_swd wave_velocity = swd.GradPhi(x, y, z);

    // 2nd order gradients of potential  (vector_2nd_phi_swd is a (xx,xy,xz,yy,yz,zz) struct)
    vector_2nd_phi_swd  grad_phi_2nd = swd.GradPhi2nd(x, y, z);

    // Local acceleration
    vector_swd euler_acceleration = swd.AccEuler(x, y, z);

    // Particle acceleration
    vector_swd particle_acceleration = swd.AccParticle(x, y, z);

    // Wave elevation
    double wave_elevation = swd.Elev(x, y);

    // Local time derivative of wave elevation
    double wave_ddt_elevation = swd.DdtElev(x, y);

    // Spatial gradient of wave elevation
    vector_swd wave_grad_elevation = swd.GradElev(x, y);

    // 2nd order gradients of elevation  (vector_2nd_elev_swd is a (xx,xy,yy) struct)
    vector_2nd_elev_swd grad_elev_2nd = swd.GradElev2nd(x, y);

    // Total pressure
    double wave_pressure = swd.Pressure(x, y, z);

    // Distance from z=0 to sea floor (<0 if infinite)
    double local_depth = swd.Bathymetry(x, y);

    // Unit normal vector of sea floor into the ocean at(x, y)
    vector_swd  nvec_sf = swd.BathymetryNvec(x, y);

    // For a specific(x, y, z) at this t, return a CSV file on how particle velocity,
    // elevation and pressure converge as a function of number of spectral components
    swd.Convergence(x, y, z, "dump.csv");

    t += dt;
}

// To save storage for an interesting event you may create a new SWD file
// containing only the time steps within a specified time window.
// In this case we only care about the time interval [0.2, 0.9]
swd.Strip(0.2, 0.9, "my_new.swd");

/*
===========================================================
There are 4 methods returning metadata from object...
  swd.GetChr(swd, name)
  swd.GetInt(swd, name)
  swd.GetBool(swd, name)
  swd.GetReal(swd, name)
where name is the requested parameter.
Five examples are given below...
===========================================================
*/

// Extract the cid string from the SWD file
// (contains typical the content of the input file applied in the wave generator)
std::string cid = swd.GetChr("cid");

// Name of actual implementation class
std::string cls_name = swd.GetChr("class");

// Extract shp parameter from object
int shp = swd.GetInt("shp");

// Extract dc_bias parameter from object
bool dc_bias2 = swd.GetBool("dc_bias");

// Extract time step as applied in SWD file
double dt_swd_file = swd.GetReal("dt");

// Automatic destructor of swd object

return 0;
}