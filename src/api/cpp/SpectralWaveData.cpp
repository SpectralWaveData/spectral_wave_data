#include <string>
#include "SpectralWaveData.h"

/*
SpectralWaveData is a C++ class for the spectral-wave-data ocean wave model.

This implementation is based on establishing wrappers to the C-implementation.

Ver-3.0:
Coded by: Jens B. Helmers DNVGL,  2019.08.11
*/

SpectralWaveData::SpectralWaveData(
    std::string file_swd,     // Name of actual SWD file
    real_swd x0, real_swd y0, // Parameters relating the SWD and application
    real_swd t0,              // coordinate system. beta is in degree.
    real_swd beta,
    real_swd rho,             // Density of water (only relevant for pressures)
    int nsumx, int nsumy,     // Number of applied spectral components. (negative=> apply all)
    int impl,                 // Index to determine actual derived class
                              //    0 = Default
                              //   <0 = In-house and experimental implementations
                              //   >0 = Validated implementations available open software
    int ipol,                 // Index to request actual temporal interpolation scheme
                              //    0 = Default (C^2 continous scheme)
                              //    1 = C^1 continous
                              //    2 = C^3 continous
    int norder,               // Expansion order to apply in kinematics for z>0
                              //    0 = Apply expansion order specified in swd file(default)
                              //   <0 = Apply exp(kj z)
                              //   >0 = Apply expansion order = norder
    bool dc_bias              // Control application of zero - frequency bias present in SWD file
                              //  false = Suppress contribution from zero frequency amplitudes(default)
                              //  true = Apply zero frequency amplitudes from SWD file.
) {
    obj = swd_api_allocate(file_swd.c_str(), x0, y0, t0, beta, rho, nsumx, nsumy, 
                           impl, ipol, norder, dc_bias);
	if (swd_api_error_raised(obj)) {
		char *msg = swd_api_error_get_msg(obj);
		switch (swd_api_error_get_id(obj)) {
		case 1001:
			throw SwdFileCantOpenException(msg);
			break;
		case 1002:
			throw SwdFileBinaryException(msg);
			break;
		case 1003:
			throw SwdFileDataException(msg);
			break;
		case 1004:
			throw SwdInputValueException(msg);
			break;
		case 1005:
			throw SwdAllocateException(msg);
			break;
		default:
			throw SwdException(msg);
			break;
		}
	}
}

SpectralWaveData::~SpectralWaveData() {
    swd_api_close(obj);
    //delete obj;
}

// =================================================================
// Methods for evaluating kinematics...
// =================================================================

// update current user time
void SpectralWaveData::UpdateTime(real_swd time)
{
    swd_api_update_time(obj, time);
	if (swd_api_error_raised(obj)) {
		throw SwdException(swd_api_error_get_msg(obj));
	}
}

// Calculate wave potential
real_swd SpectralWaveData::Phi(real_swd x, real_swd y, real_swd z)
{
    return swd_api_phi(obj, x, y, z);
}

// Calculate stream function
real_swd SpectralWaveData::Stream(real_swd x, real_swd y, real_swd z)
{
    return swd_api_stream(obj, x, y, z);
}

// Calculate time derivative of wave potential (earth fixed observer)
real_swd SpectralWaveData::DdtPhi(real_swd x, real_swd y, real_swd z)
{
    return swd_api_phi_t(obj, x, y, z);
}

// Calculate particle velocity
vector_swd SpectralWaveData::GradPhi(real_swd x, real_swd y, real_swd z)
{
    return swd_api_grad_phi(obj, x, y, z);
}

// 2nd order gradients of potential
vector_2nd_phi_swd SpectralWaveData::GradPhi2nd(real_swd x, real_swd y, real_swd z)
{
    return swd_api_grad_phi_2nd(obj, x, y, z);
}

// Calculate Euler acceleration: grad(phi_t)
vector_swd SpectralWaveData::AccEuler(real_swd x, real_swd y, real_swd z)
{
    return swd_api_acc_euler(obj, x, y, z);
}

// Calculate particle acceleration
vector_swd SpectralWaveData::AccParticle(real_swd x, real_swd y, real_swd z)
{
    return swd_api_acc_particle(obj, x, y, z);
}

// Calculate wave elevation
real_swd SpectralWaveData::Elev(real_swd x, real_swd y)
{
    return swd_api_elev(obj, x, y);
}

// Calculate time derivative of wave elevation (Euler derivative, earth fixed)
real_swd SpectralWaveData::DdtElev(real_swd x, real_swd y)
{
    return swd_api_elev_t(obj, x, y);
}

// Calculate gradient of wave elevation (slopes)
vector_swd SpectralWaveData::GradElev(real_swd x, real_swd y)
{
    return swd_api_grad_elev(obj, x, y);
}

// 2nd order gradients of elevation
vector_2nd_elev_swd SpectralWaveData::GradElev2nd(real_swd x, real_swd y)
{
    return swd_api_grad_elev_2nd(obj, x, y);
}

// Calculate complete Bernoulli pressure
real_swd SpectralWaveData::Pressure(real_swd x, real_swd y, real_swd z)
{
    return swd_api_pressure(obj, x, y, z);
}

// Vertical distance from z=0 to sea floor (<0 if infinite)
real_swd SpectralWaveData::Bathymetry(real_swd x, real_swd y)
{
    return swd_api_bathymetry(obj, x, y);
}

// Unit normal vector of sea floor into the ocean at(x, y)
vector_swd SpectralWaveData::BathymetryNvec(real_swd x, real_swd y)
{
    return swd_api_bathymetry_nvec(obj, x, y);
}

// For a specific location return a csv-file on how velocity, elevation
// and pressure converge as a function of number of spectral components
void SpectralWaveData::Convergence(real_swd x, real_swd y, real_swd z, std::string csv)
{
    swd_api_convergence(obj, x, y, z, csv.c_str());
	if (swd_api_error_raised(obj)) {
		throw SwdException(swd_api_error_get_msg(obj));
	}
}

// To save storage for an interesting event you may create a new SWD file
// containing only the time steps within the time window [tmin, tmax].
// The name of the new file is defined by file_swd_new
void SpectralWaveData::Strip(real_swd tmin, real_swd tmax, std::string file_swd_new)
{
    swd_api_strip(obj, tmin, tmax, file_swd_new.c_str());
	if (swd_api_error_raised(obj)) {
		throw SwdException(swd_api_error_get_msg(obj));
	}
}

// ===================================================================
//  Provide parameters from the swd-file:

// Extract the character parameter 'name' from object
std::string SpectralWaveData::GetChr(std::string const& name)
{
    return swd_api_get_chr(obj, name.c_str());
}

// Extract the int parameter 'name' from object
int SpectralWaveData::GetInt(std::string const& name)
{
    return swd_api_get_int(obj, name.c_str());
}

// Extract the int parameter 'name' from object
bool SpectralWaveData::GetBool(std::string const& name)
{
	return swd_api_get_bool(obj, name.c_str());
}

// Extract the real parameter 'name' from object
real_swd SpectralWaveData::GetReal(std::string const& name)
{
    return swd_api_get_real(obj, name.c_str());
}
// ===================================================================

// Clear error flag in C/Fortran implementation.
// Only applied in case of advanced exception handling. (recovery)
void SpectralWaveData::ExceptionClear()
{
	swd_api_error_clear(obj);
}
