===============
Wave generators
===============

There are no restrictions on programming styles applied in
:ref:`wave generators<wave-generator>` as long as it can produce a sound SWD-file.

:ref:`Wave generators<wave-generator>` are typical limited to provide wave fields
described by a small set of spectral formulations. Often only one or two formulations.
Consequently, only a small subset of the :doc:`SWD-specification <swd_format>` needs
to be considered when including support for writing SWD files.

-------
Fortran
-------

Fortran programmers may apply the following section for how to
write a sound SWD file from a :ref:`wave generator<wave-generator>`.
The handling of making a C binary stream is demonstrated.

For each spectral formulation this repository provides source code demonstrating
how to create proper SWD-files with minimal effort.
Implementation or adaption of this source
code to the actual :ref:`wave generator<wave-generator>` is
expected to be a very small job.

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Spectral model
     - Example source files for wave generators
   * - shape 1 or 2
     - :doc:`swd_write_shape_1_or_2.f90 <fortran_swd_write_shape_1_or_2>`
   * - shape 3
     - :doc:`swd_write_shape_3.f90 <fortran_swd_write_shape_3>`
   * - shape 4 or 5
     - :doc:`swd_write_shape_4_or_5.f90 <fortran_swd_write_shape_4_or_5>`
   * - shape 6
     - :doc:`swd_write_shape_6.f90 <fortran_swd_write_shape_6>`

.. toctree::
   :hidden:

   fortran_swd_write_shape_1_or_2
   fortran_swd_write_shape_3
   fortran_swd_write_shape_4_or_5
   fortran_swd_write_shape_6

.. note::
  If your :ref:`wave generator<wave-generator>` applies a different numerical precision than
  :attr:`double precision` you need to adjust the kind parameter :attr:`wp` defined in the beginning
  of the source code.

Assuming your wave generator will support the spectral shape_X formulation
where X is one of the shape models as defined in the :doc:`theory <theory>` section,
you may write something like this in your :ref:`wave generator<wave-generator>`.

.. code-block:: fortran

   program my_wave_generator
   ...
   use swd_write_shape_X_def, only: swd_write_shape_X   ! Replace X with the index of the actual model
   ...
   type(swd_write_shape_X) :: swd_X
   ...
   swd_X % init('my_wave.swd', nx, dk, ...)  ! All time independent swd-format parameters for shape X
   ...
   do (temporal integration loop)
       ...
       swd_X % update(h, ht, c, ct)   ! Relevant spectral amplitudes for this time step for shape X
   end do
   swd_X % close()
   ...
   program my_wave_generator

If your :ref:`wave generator<wave-generator>` also supports another spectral formulation :attr:`shp=y`,
just add a similar :attr:`swd_y`
variable of :attr:`type(swd_write_shape_y)` and provide the relevant parameters for
that model too.

Following this approach, Fortran source code dealing with the low level binary C-stream output
for the various spectral formulations are listed in this table


----------------
C / C++ / Python
----------------

Writing C-streams using these programming languages are straight forward.
The logical ordering is the same as for the fortran examples given above.
An example from open source is the `raschii <https://raschii.readthedocs.io/>`_ Python
package for making non-linear regular waves.

.. note::
   Any 2-dimensional arrays should be written using the Fortran element ordering described
   in the :doc:`SWD format<swd_format>` description.
