from spectral_wave_data import SpectralWaveData
from spectral_wave_data.tools import airy

"""
This function was used for testing of proper garbage collection in case someone
forget to close the swd object before the object goes out of scope.

To do this check again, you need to temporarily do the following:

Include a print statement in the top of the "__del__" method
in spectral_wave_data.py to see if the destruction is executed.
"""

def test_check_out_of_scope_destructor_message():
    file_swd = "my.swd"
    depth = 17.3
    grav = 9.81
    amps_inp = [1.2, 0.4, 0.8]
    dirs_inp = [173.2, 25.0, -130.0]
    phases_inp = [0.0, 210.0, 70.0]
    Twaves_inp = [3.0, 11.0, 70.0]
    airy.write_swd(file_swd, amps_inp, dirs_inp, phases_inp, Twaves=Twaves_inp, depth=depth, grav=grav)
    swd = SpectralWaveData(file_swd)
    t = 3.0
    x = 1.4
    y = 4.5
    swd.update_time(t)
    zeta = swd.elev(x, y)
    print("Leaving the function...")


test_check_out_of_scope_destructor_message()
print("Leaving the script...")
