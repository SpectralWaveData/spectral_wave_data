#include <stdlib.h>
#include <assert.h>
#include "ISO_Fortran_binding.h"
#include "spectral_wave_data.h"


CFI_cdesc_t *swd_api_elev_fft(void *swd, int nx_fft, int ny_fft) {
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

void swd_api_elev_fft_deallocate(CFI_cdesc_t *desc_elev_arr) {
    int rc = CFI_deallocate(desc_elev_arr);
    assert(CFI_SUCCESS == rc);
}

CFI_cdesc_t *swd_api_x_fft(void *swd, int nx_fft) {
    CFI_cdesc_t *desc_x_arr = (CFI_cdesc_t *) malloc(sizeof(CFI_CDESC_T(1)));

    int rc = CFI_establish(desc_x_arr,
		                   NULL,
                           CFI_attribute_allocatable,
                           CFI_type_double,
                           sizeof(double),
                           (CFI_rank_t)1,
                           NULL);
    assert(CFI_SUCCESS == rc);

    swd_api_x_fft_(swd, nx_fft, desc_x_arr);

    return desc_x_arr;

}

void swd_api_x_fft_deallocate(CFI_cdesc_t *desc_x_arr) {
    int rc = CFI_deallocate(desc_x_arr);
    assert(CFI_SUCCESS == rc);
}

CFI_cdesc_t *swd_api_y_fft(void *swd, int ny_fft) {
    CFI_cdesc_t *desc_y_arr = (CFI_cdesc_t *) malloc(sizeof(CFI_CDESC_T(1)));

    int rc = CFI_establish(desc_y_arr,
		                   NULL,
                           CFI_attribute_allocatable,
                           CFI_type_double,
                           sizeof(double),
                           (CFI_rank_t)1,
                           NULL);
    assert(CFI_SUCCESS == rc);

    swd_api_y_fft_(swd, ny_fft, desc_y_arr);

    return desc_y_arr;

}

void swd_api_y_fft_deallocate(CFI_cdesc_t *desc_y_arr) {
    int rc = CFI_deallocate(desc_y_arr);
    assert(CFI_SUCCESS == rc);
}
