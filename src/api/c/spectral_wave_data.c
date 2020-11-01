#include <stdlib.h>
#include <assert.h>
#include "ISO_Fortran_binding.h"
#include "spectral_wave_data.h"


CFI_cdesc_t *swd_api_elev_fft_obj(void *swd, int nx_fft, int ny_fft) {
    CFI_cdesc_t *desc_elev_arr = (CFI_cdesc_t *) malloc(sizeof(CFI_CDESC_T(2)));

    int rc = CFI_establish(desc_elev_arr,
		    NULL,
            CFI_attribute_allocatable,
            CFI_type_double,
            sizeof(double),
            (CFI_rank_t)2,
            NULL);
    assert(CFI_SUCCESS == rc);

    swd_api_elev_fft_(swd, nx_fft, ny_fft, desc_elev_arr);

    return desc_elev_arr;

}

void swd_api_elev_fft_obj_deallocate(CFI_cdesc_t *desc_elev_arr) {
    int rc = CFI_deallocate(desc_elev_arr);
    assert(CFI_SUCCESS == rc);
}
