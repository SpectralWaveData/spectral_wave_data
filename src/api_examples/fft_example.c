#include <stdbool.h> /* bool */
#include <stdio.h> /* for fprintf and stderr */
#include <stdlib.h> /* for exit */
#include <assert.h>
#include <string.h>
#include <math.h>

#include "ISO_Fortran_binding.h"
#include "spectral_wave_data.h"    // Namespace convention: swd_api_*

#define MAX(a, b) ((a) > (b) ? (a) : (b))

int main() {
void *swd;

// First we define constructor parameters (see application_swd_api.c)
char *file_swd = "n5.swd";
double x0 = 0.0, y0 = 0.0, t0 = 0.0, beta = 0.0;
int impl = 0;
int nsumx = -1, nsumy = -1;   // Negative: apply all components from the swd-file.
double rho = 1025.0;
int ipol = 0;   // C^2 continous
int norder = 0;   // Apply same order as specified in SWD file
bool dc_bias = false;   // Suppress contributions from zero-frequency

// Constructor
swd = swd_api_allocate(file_swd, x0, y0, t0, beta, rho, nsumx, nsumy, 
                       impl, ipol, norder, dc_bias);
if (swd_api_error_raised(swd)) {
    fprintf(stderr, "%s", swd_api_error_get_msg(swd));
    exit(EXIT_FAILURE); /* indicate failure.*/
}

swd_api_update_time(swd, 0.0);

double x, y, eta_fft, eta;

CFI_cdesc_t *elev_fft = swd_api_elev_fft(swd, -1, -1);
CFI_cdesc_t *x_fft = swd_api_x_fft(swd, -1);
CFI_cdesc_t *y_fft = swd_api_y_fft(swd, -1);

if (swd_api_error_raised(swd)) {
    fprintf(stderr, "%s", swd_api_error_get_msg(swd));
    exit(EXIT_FAILURE); /* indicate failure.*/
}

// get the array dimensions
int nx_fft = x_fft->dim[0].extent;
int ny_fft = y_fft->dim[0].extent;
assert(elev_fft->dim[0].extent == nx_fft);
assert(elev_fft->dim[1].extent == ny_fft);

// get the lower bounds (important for the indexing below)
// NOTE: Currently this seems to follow Fortran indexing with 1 as lower bounds, not sure if this
//       is intended or if it a "bug" in gfortran/gcc-10? Anyway, as long as the lower bounds are taken
//       into account in the loops below it should be ok.
int lbx = x_fft->dim[0].lower_bound;  // seems to always be 1
int lby = y_fft->dim[0].lower_bound;  // seems to always be 1
assert(elev_fft->dim[0].lower_bound == lbx);
assert(elev_fft->dim[1].lower_bound == lby);

// for indexing arrays
CFI_index_t sub_x[1], sub_y[1], sub_xy[2];

double err_max = 0.0;
double err;
for (int j = lby; j < ny_fft + lby; j++) {
    sub_y[0] = j;
    sub_xy[1] = j;
    y = *((double *) CFI_address( y_fft,  sub_y));
    // for (int i = lbx; i < nx_fft + lbx; i++) {
    for (int i = lbx; i < 4 + lbx; i++) {
        sub_x[0] = i;
        sub_xy[0] = i;
        x = *((double *) CFI_address( x_fft,  sub_x));
        eta_fft = *((double *) CFI_address( elev_fft, sub_xy));

        // make sure we get the same result from standard elev-routine'
        eta = swd_api_elev(swd, x, y);
        err = fabs(eta - eta_fft);
        err_max = MAX(err_max, err);
    }  
}
printf("Largest error between elev() and elev_fft(): err_max = %g \n", err_max);

// close and deallocate
swd_api_fft_deallocate(elev_fft);
swd_api_fft_deallocate(x_fft);
swd_api_fft_deallocate(y_fft);
swd_api_close(swd);

return 0;
}
