**********
SWD format
**********

In this section we describe the format of the SWD-file.
Symbols follow closely the definitions in the :doc:`theory<theory>` section.

The units of all quantities are required to be consistent with the
`SI-system <https://en.wikipedia.org/wiki/International_System_of_Units>`_.


The SWD-file is a C-stream binary file (see the :ref:`rationals<swd-c-stream-why>`).

.. note::

  In order to ensure portability between various computer systems SWD files are required
  to be stored using the 'little endian' byte order convention. This is the default
  convention for all Intel/AMD CPUs running Windows, Linux or Apple.
  PowerPC CPUs (older Macintosh machines and IBM clusters and some supercomputers) use
  big endian. Wave generators making SWD files on big endian systems are required to
  specify a flag to signal that the SWD file should be written in little endian. This
  is usually specifed in the file open statement, but can also be defined using environment
  variables.

The binary stream in the SWD file is outlined in the following pseudo code:

.. code-block:: fortran

   magic
   fmt
   shp, amp
   prog
   date
   nid
   cid
   grav, lscale
   nstrip, nsteps, dt
   order
   select case(shp)
   case(1)
      n, dk
      do i = 1, nsteps
          h(0), h(1), ..., h(n)
          ht(0), ht(1), ..., ht(n)
          if amp < 3 then
             c(0), c(1), ..., c(n)
             ct(0), ct(1), ..., ct(n)
          end if
      end do
   case(2)
      n, dk, d
      do i = 1, nsteps
          h(0), h(1), ..., h(n)
          ht(0), ht(1), ..., ht(n)
          if amp < 3 then
             c(0), c(1), ..., c(n)
             ct(0), ct(1), ..., ct(n)
          end if
      end do
   case(3)
      n, nh, dk, isf
      nsf
      xsf(1), xsf(2), ..., xsf(nsf)
      zsf(1), zsf(2), ..., zsf(nsf)
      do i = 1, nsteps
          h(0), h(1), ..., h(n)
          ht(0), ht(1), ..., ht(n)
          if amp < 3 then
              c(0), c(1), ..., c(n)
              ct(0), ct(1), ..., ct(n)
              if nsf > 1 then
                 ch(0), ch(1), ..., ch(nh)
                 cht(0), cht(1), ..., cht(nh)
              end if
          end if
      end do
   case(4)
      nx, ny, dkx, dky
      do i = 1, nsteps
          h(-ny,0), h(-ny+1,0), ..., h(ny-1,nx), h(ny,nx)        ! Fortran element order
          ht(-ny,0), ht(-ny+1,0), ..., ht(ny-1,nx), ht(ny,nx)
          if amp < 3 then
              c(-ny,0), c(-ny+1,0), ..., c(ny-1,nx), c(ny,nx)
              ct(-ny,0), ct(-ny+1,0), ..., ct(ny-1,nx), ct(ny,nx)
          end if
      end do
   case(5)
      nx, ny, dkx, dky, d
      do i = 1, nsteps
          h(-ny,0), h(-ny+1,0), ..., h(ny-1,nx), h(ny,nx)        ! Fortran element order
          ht(-ny,0), ht(-ny+1,0), ..., ht(ny-1,nx), ht(ny,nx)
          if amp < 3 then
              c(-ny,0), c(-ny+1,0), ..., c(ny-1,nx), c(ny,nx)
              ct(-ny,0), ct(-ny+1,0), ..., ct(ny-1,nx), ct(ny,nx)
          end if
      end do
   case(6)
      n, d
      do i = 1, n
          amp(i), kw(i), gam(i), phs(n)
      end do
   end select

where

.. list-table::
   :widths: 10 10 10 70
   :header-rows: 1

   * - name
     - C type
     - bytes
     - description
   * - magic
     - float
     - 4
     - | magic = 37.0221
       | Constant decimal number in all SWD files. (future proof)
   * - fmt
     - int
     - 4
     - | Integer to identify version of the file format.
       | In this version: fmt = 100
   * - shp
     - int
     - 4
     - | Actual shape functions as defined in the :doc:`theory <theory>` section.
       | (e.g. shp = 3 implies :doc:`Shape 3 <shape_3>` is applied)
   * - amp
     - int
     - 4
     - | Flag to indicate which temporal amplitudes are stored in the swd-file
       |  1: All temporal amplitudes related to the shape class is specified
       |     as defined in the :doc:`theory <theory>` section.
       |  2: All amplitudes related to the shape class is specified. However,
       |     the velocity potential is only interpreted on the free surface.
       |     The :math:`z`-dependencies are removed from all formulas.
       |     This option is introduced for reasearch on how to map the
       |     potential to other vertical locations. More accurate calculations
       |     are possible but at significantly higher computational cost.
       |     The current implementation of the API does not support this feature.
       |  3: Functions related to the velocity potential are not specified.
       |     It is only possible to evaluate surface elevation quantities
       |     with this option. For other calculations :math:`\phi\equiv 0` is applied.
   * - prog
     - 30 char
     - 30
     - The name of the program building
       this file. (Including version number)
   * - date
     - 20 char
     - 20
     - Text providing the date and time this file was build.
   * - nid
     - int
     - 4
     - Number of characters describing the next field. (nid > 0)
   * - cid
     - nid char
     - nid
     - | Text to describe the wave field. It is expected to be the complete
       | content of the input file applied in the wave generator.
   * - grav
     - float
     - 4
     - | Acceleration of gravity applied in the wave generator
       | (Not applied in current version)
   * - lscale
     - float
     - 4
     - | Number of length units per meter applied in the wave generator.
       | (Not applied in current version)
   * - nstrip
     - int
     - 4
     - | Number of initial time steps stripped off from the original simulation.
       | Default nstrip=0. The strip() method will remove some initial and trailing
       | time steps. nstrip can be used to deduce the original time reference.
   * - nsteps
     - int
     - 4
     - Number of time steps stored in this swd file.
   * - dt
     - float
     - 4
     - Constant time step for spectral amplitudes stored in this file
   * - order
     - int
     - 4
     - | Perturbation order applied in the :ref:`wave generator<wave-generator>`.
       | (<0 if fully non-linear)
   * - n, nx
     - int
     - 4
     - Number of spectral components, :math:`n` or :math:`n_x`, in the :math:`x`-direction.
   * - ny
     - int
     - 4
     - Number of spectral components :math:`n_y` in the :math:`y`-direction.
   * - nh
     - int
     - 4
     - | Number of auxiliary spectral components :math:`\hat{n}` in the
       | :math:`x`-direction in case of bathymetry.
   * - dk, dkx
     - float
     - 4
     - Spacing of wave numbers, :math:`\Delta k` or :math:`\Delta k_x`, in the :math:`x`-direction.
   * - dky
     - float
     - 4
     - Spacing of wave numbers :math:`\Delta k_y` in the :math:`y`-direction.
   * - d
     - float
     - 4
     - Constant or average water depth :math:`d`. (<0 if infinite)
   * - isf
     - int
     - 4
     - | Flag to indicate geometric description of the sea floor
       |  0: Piecewise linear sea floor
   * - nsf
     - int
     - 4
     - | Number of offset points defining the sea floor in :math:`x`-direction.
       |  0: Infinite water depth
       |  1: Constant water depth
       |  >1: Varying water depth (bathymetry)
   * - xsf()
     - float
     - 4
     - | :math:`x`-positions of offset points defining the sea floor.
       | xsf() should cover the range :math:`x\in[0, 2\pi/\Delta k]`.
   * - zsf()
     - float
     - 4
     - :math:`z`-positions of offset points defining the sea floor.
   * - h()
     - complex
     - 4+4
     - Spectral amplitudes :math:`h()`  (real and imaginary part)
   * - ht()
     - complex
     - 4+4
     - Spectral amplitudes :math:`\frac{dh}{dt}()`  (real and imaginary part)
   * - c()
     - complex
     - 4+4
     - Spectral amplitudes :math:`c()`  (real and imaginary part)
   * - ct()
     - complex
     - 4+4
     - Spectral amplitudes :math:`\frac{dc}{dt}()`  (real and imaginary part)
   * - ch()
     - complex
     - 4+4
     - Spectral amplitudes :math:`\hat{c}()`  (real and imaginary part)
   * - cht()
     - complex
     - 4+4
     - Spectral amplitudes :math:`\frac{d\hat{c}}{dt}()`  (real and imaginary part)
   * - amp()
     - float
     - 4
     - Single amplitudes :math:`A_j` for Airy model (shp=6).
   * - kw()
     - float
     - 4
     - Wave number :math:`k_j` for Airy model (shp=6).
   * - gam()
     - float
     - 4
     - Wave direction :math:`\gamma_j` (rad) for Airy model (shp=6).
   * - phs()
     - float
     - 4
     - Wave phase :math:`\delta_j` (rad) for Airy model (shp=6).

