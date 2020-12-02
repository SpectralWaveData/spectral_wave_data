#include <stdlib.h>
#include <assert.h>
#include "ISO_Fortran_binding.h"
#include "spectral_wave_data.h"


CFI_cdesc_t *swd_api_elev_fft(void *swd, int nx_fft, int ny_fft) {
    CFI_cdesc_t *desc_elev_arr = (CFI_cdesc_t *) malloc(sizeof(CFI_CDESC_T(2)));

    CFI_establish(desc_elev_arr,
                  NULL,
                  CFI_attribute_allocatable,
                  CFI_type_double,
                  sizeof(double),
                  (CFI_rank_t)2,
                  NULL);
    
    swd_api_elev_fft_(swd, nx_fft, ny_fft, desc_elev_arr);

    return desc_elev_arr;

}

CFI_cdesc_t *swd_api_grad_phi_fft(void *swd, double z, int nx_fft, int ny_fft) {
    CFI_cdesc_t *desc_grad_phi_arr = (CFI_cdesc_t *) malloc(sizeof(CFI_CDESC_T(3)));

    CFI_establish(desc_grad_phi_arr,
                  NULL,
                  CFI_attribute_allocatable,
                  CFI_type_double,
                  sizeof(double),
                  (CFI_rank_t)3,
                  NULL);
    
    swd_api_grad_phi_fft_(swd, z, nx_fft, ny_fft, desc_grad_phi_arr);

    return desc_grad_phi_arr;

}

CFI_cdesc_t *swd_api_x_fft(void *swd, int nx_fft, int ny_fft) {
    CFI_CDESC_T(2) tmp_y;
    CFI_cdesc_t *desc_x_arr = (CFI_cdesc_t *) malloc(sizeof(CFI_CDESC_T(2)));
    CFI_cdesc_t *desc_y_arr = (CFI_cdesc_t *) &tmp_y;

    CFI_establish(desc_x_arr,
                  NULL,
                  CFI_attribute_allocatable,
                  CFI_type_double,
                  sizeof(double),
                  (CFI_rank_t)2,
                  NULL);
    CFI_establish(desc_y_arr,
                  NULL,
                  CFI_attribute_allocatable,
                  CFI_type_double,
                  sizeof(double),
                  (CFI_rank_t)2,
                  NULL);
    
    swd_api_xy_fft_(swd, desc_x_arr, desc_y_arr, nx_fft, ny_fft);
    CFI_deallocate(desc_y_arr);

    return desc_x_arr;

}

CFI_cdesc_t *swd_api_y_fft(void *swd, int nx_fft, int ny_fft) {
    CFI_CDESC_T(2) tmp_x;
    CFI_cdesc_t *desc_x_arr = (CFI_cdesc_t *) &tmp_x;
    CFI_cdesc_t *desc_y_arr = (CFI_cdesc_t *) malloc(sizeof(CFI_CDESC_T(2)));

    CFI_establish(desc_x_arr,
                  NULL,
                  CFI_attribute_allocatable,
                  CFI_type_double,
                  sizeof(double),
                  (CFI_rank_t)2,
                  NULL);
    CFI_establish(desc_y_arr,
                  NULL,
                  CFI_attribute_allocatable,
                  CFI_type_double,
                  sizeof(double),
                  (CFI_rank_t)2,
                  NULL);
    
    swd_api_xy_fft_(swd, desc_x_arr, desc_y_arr, nx_fft, ny_fft);
    CFI_deallocate(desc_x_arr);

    return desc_y_arr;

}

void swd_api_fft_deallocate(CFI_cdesc_t *desc_arr) {
    CFI_deallocate(desc_arr);
}
