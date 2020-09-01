#ifndef SWD_CPP_H_INCLUDED
#define SWD_CPP_H_INCLUDED

/*
The C++ interface for the spectral-wave-data ocean wave model.

Ver-3.0: 
Coded by: Jens B. Helmers DNVGL,  2019.08.11

*/
#include <string>
#include <stdexcept>

#include "spectral_wave_data.h"  // The core C interface

// Specific Exception Classes...

class SwdException : public std::runtime_error
{
public:
    SwdException(const char* msg) : std::runtime_error(msg) { }
};

class SwdFileCantOpenException : public SwdException
{
public:
    SwdFileCantOpenException(const char* msg) : SwdException(msg) { }
};

class SwdFileBinaryException : public SwdException
{
public:
    SwdFileBinaryException(const char* msg) : SwdException(msg) { }
};

class SwdFileDataException : public SwdException
{
public:
    SwdFileDataException(const char* msg) : SwdException(msg) { }
};

class SwdInputValueException : public SwdException
{
public:
    SwdInputValueException(const char* msg) : SwdException(msg) { }
};

class SwdAllocateException : public SwdException
{
public:
    SwdAllocateException(const char* msg) : SwdException(msg) { }
};


// The main class...

class SpectralWaveData
{
    void* obj; // Wrapper to the C-object of the ocean wave model

public:
    SpectralWaveData(std::string file_swd, real_swd x0, real_swd y0, 
                     real_swd t0, real_swd beta, real_swd rho=1025.0, 
                     int nsumx=-1, int nsumy=-1, int impl=0, int ipol=0,
                     int norder=0, bool dc_bias=false);
    /*
    file_swd:          Name of actual swd file

    x0, y0, t0, beta:  Relation between SWD and application coordinates.
                       beta in degree.

    rho:               Density of water(applied for pressure calculations)
    nsumx, nsumy       Number of spectral components to apply (<0: apply all)

    impl               Index to determine actual derived class
                       0 = Default
                      <0 = In-house and experimental implementations
                      >0 = Validated implementations available open software

    ipol               Index to request actual temporal interpolation scheme
                       0 = Default (C^2 continous scheme)
                       1 = C^1 continous
                       2 = C^3 continous

    norder             Expansion order to apply in kinematics for z>0
                       0 = Apply expansion order specified in swd file (default)
                       <0 = Apply exp(kj z)
                       >0 = Apply expansion order = norder

    dc_bias            Control application of zero-frequency bias present in SWD file
                       false = Suppress contribution from zero frequency amplitudes (default)
                       true  = Apply zero frequency amplitudes from SWD file.

    C++ exceptions the constructor may throw:  (should be catched in application)

    SwdException:               Base class for the following exceptions
    SwdFileCantOpenException:   Typical if SWD file is not an existing file
    SwdFileBinaryException:     SWD file does not apply float/little-endian
    SwdFileDataException:       Error during reading and checking data
    SwdInputValueException:     Input arguments for class methods are not sound
    SwdAllocateException:       Not able to allocate internal SWD storage

    */

    ~SpectralWaveData();

    // =================================================================
    // Methods for evaluating kinematics...
    // =================================================================

    // Apply current application time
    // Possible exceptions thrown: SwdFileDataException, SwdInputValueException
    void UpdateTime(real_swd time);

    // Calculate velocity potential
    // Possible exceptions thrown: None
    real_swd Phi(real_swd x, real_swd y, real_swd z);

    // Calculate stream function
    // Possible exceptions thrown: None
    real_swd Stream(real_swd x, real_swd y, real_swd z);
    
    // Calculate time derivative of wave potential (earth fixed observer)
    // Possible exceptions thrown: None
    real_swd DdtPhi(real_swd x, real_swd y, real_swd z);

    // Calculate particle velocity
    // Possible exceptions thrown: None
    vector_swd GradPhi(real_swd x, real_swd y, real_swd z);

    // 2nd order gradients of potential
    // Possible exceptions thrown: None
    vector_2nd_phi_swd GradPhi2nd(real_swd x, real_swd y, real_swd z);

    // Calculate Euler acceleration: grad(phi_t)
    // Possible exceptions thrown: None
    vector_swd AccEuler(real_swd x, real_swd y, real_swd z);

    // Calculate particle acceleration
    // Possible exceptions thrown: None
    vector_swd AccParticle(real_swd x, real_swd y, real_swd z);
    
    // Calculate wave elevation
    // Possible exceptions thrown: None
    real_swd Elev(real_swd x, real_swd y);

    // Calculate time derivative of wave elevation (Euler derivative, earth fixed)
    // Possible exceptions thrown: None
    real_swd DdtElev(real_swd x, real_swd y);

    // Calculate gradient of wave elevation (slopes)
    // Possible exceptions thrown: None
    vector_swd GradElev(real_swd x, real_swd y);

    //  2nd order gradients of elevation
    // Possible exceptions thrown: None
    vector_2nd_elev_swd GradElev2nd(real_swd x, real_swd y);
    
    // Complete Bernoulli pressure
    // Possible exceptions thrown: None
    real_swd Pressure(real_swd x, real_swd y, real_swd z);

    // Vertical distance from z=0 to sea floor (<0 if infinite)
    // Possible exceptions thrown: None
    real_swd Bathymetry(real_swd x, real_swd y);

    // Unit normal vector of sea floor into the ocean at(x, y)
    // Possible exceptions thrown: None
    vector_swd BathymetryNvec(real_swd x, real_swd y);

    // For a specific location return a csv-file on how velocity, elevation
    // and pressure converge as a function of number of spectral components
    // Possible exceptions thrown: SwdFileCantOpenException
    void Convergence(real_swd x, real_swd y, real_swd z, std::string csv);

    // To save storage for an interesting event you may create a new SWD file
    // containing only the time steps within the time window [tmin, tmax].
    // The name of the new file is defined by file_swd_new
    // Possible exceptions thrown: SwdFileCantOpenException
    //                             SwdInputValueException
    void Strip(real_swd tmin, real_swd tmax, std::string file_swd_new);

    // ===================================================================
    //  Provide parameters from the swd-file:

    // Extract the character parameter 'name' from object
    // Possible exceptions thrown: SwdInputValueException
    std::string GetChr(std::string const& name);

    // Extract the int parameter 'name' from object
    // Possible exceptions thrown: SwdInputValueException
    int GetInt(std::string const& name);

	// Extract the int parameter 'name' from object
	// Possible exceptions thrown: SwdInputValueException
	bool GetBool(std::string const& name);

	// Extract the real parameter 'name' from object
    // Possible exceptions thrown: SwdInputValueException
    real_swd GetReal(std::string const& name);
    // ===================================================================

    // Clear error signal in C/Fortran code in case of recovering from
    // advanced exception handling
    // Possible exceptions thrown: None
    void ExceptionClear();

};
#endif // SWD_CPP_H_INCLUDED
