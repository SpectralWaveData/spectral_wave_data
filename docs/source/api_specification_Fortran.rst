************
Fortran-2008
************

Using standard ISO Fortran-2008 the name of the abstract base class is
:class:`spectral_wave_data` which is available from the module
:mod:`spectral_wave_data_def` defined in :doc:`spectral_wave_data.f90<fortran_swd_base>`.
The abstract interface block in this module defines the signature for all applied methods.

The constructor for allocating and initializing a proper
specialized class is defined in :doc:`spectral_wave_data_allocate.f90<fortran_swd_allocate>`.
All constructor arguments are explained in the header of this constructor.

.. code-block:: fortran

   program my_application

   use spectral_wave_data_def,   only: spectral_wave_data
   use spectral_wave_data_allocate_def,   only: spectral_wave_data_allocate
   ...
   class(spectral_wave_data), allocatable :: swd   ! Type is not known at this stage
   ...
   file_swd = 'my_wave.swd'    ! File based on any supported spectral formulation
   ...
   ! Allocate and initialize the swd using the default implementation for the
   ! actual spectral formulation as defined in the swd file...
   call spectral_wave_data_allocate(swd, file_swd, x0, y0, t0, beta)
   if (swd % error % raised()) then
       print*, swd % error % get_msg()
       stop
   end if

   t = 0.0_wp; dt = 0.1_wp; tmax = swd % get_real('tmax')
   do
      if (t > tmax) exit
      call swd % update_time(t)

     ! apply generic swd % methods() according to the interface defined in spectral_wave_data.f90

     t = t + dt
   end do
   call swd % close()

   end program my_application

Application developers may modify the :attr:`SELECT CASE` block in
:doc:`spectral_wave_data_allocate.f90<fortran_swd_allocate>` to include support for new implementations
derived from the base class.

.. toctree::
   :hidden:

   fortran_swd_base.rst
   fortran_swd_allocate.rst

