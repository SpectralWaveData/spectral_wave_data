**********
Rationales
**********

Why spectral waves?
-------------------

There are several advanced wave models in academic and industrial applications.
All with pros and cons. A particular challenge is that a wide range of length
and time scales are relevant for ocean wave kinematics. Some effects like modulation
instabilities may develop over a long range of time and space, while formation of
breaking waves may occur very locally and quickly. Consequently, there are no
model practical for all assessments.

This API is exclusively dedicated to spectral wave formulations for the following reasons

 - It provides a continuous and smooth description of the wave field.
   Hence, :ref:`application programs<application-program>` using different
   numerical models and discretizations  may evaluate the same wave field
   at arbitrary locations without field interpolation.
 - The wave field is completely described by a relative small set of data.
   This makes it possible to do large scale computations in time and space.
 - Spectral wave fields are typical applied as boundary conditions in
   advanced :ref:`application programs<application-program>`. Such programs may simulate
   overturning and breaking waves and/or responses of marine structures.
 - Spectral waves can describe all relevant classical wave models (Airy, Stokes, etc)
   and linear and non-linear sea states. Capturing several important non-linear physical effects
   are documented for all levels of water depths including bathymetry, for both long
   and short-crested waves.
 - The required set of kinematics to evaluate may differ between
   :ref:`application programs<application-program>`.
   E.g. perturbation based programs may need to evaluate higher-order gradients of the velocity potential
   while a CFD program may only need surface elevation and particle velocities.
   Without any such considerations when first constructing the field,
   the `spectral_wave_data <https://github.com/SpectralWaveData/spectral_wave_data>`_ API
   can provide all relevant kinematics without introducing further mathematical approximations.
 - :ref:`Wave generators<wave-generator>` are not limited to e.g.
   Higher-Order-Spectral-Method HOS(M) or Boussinesq solvers. Output from other numerical
   or experimental wave fields may be fitted to a spectral formulation. This
   is e.g. relevant when doing wave reconstruction based on experimental data from
   surface elevation probes. Such wave fields may then be analyzed using any
   :ref:`application program<application-program>`.

Why introduce the SWD file?
---------------------------

The main reasons for introducing the SWD file and not to couple programs directly at the source code
level are

 - The spectral wave field analyzed by :ref:`application programs<application-program>`
   may come from very different sources. E.g.

   - Different institutions (academic and industry)
   - Different theories (e.g. HOSM vs Boussinesq)
   - Fitting of wave elevation probes (different model basins)

 - It is technically complicated to couple :ref:`wave generators<wave-generator>` and
   :ref:`application programs<application-program>` at the source code level. Especially if you
   need wave fields from several sources as described above.

 - It takes much less effort to assess different sources using the file approach, which we think
   benefit both the research and industry communities.

 - There might be severe copyright issues when coupling software at the source code.

 - If you need to assess the same wave field several times (e.g. assess two different loading
   conditions of a ship) it is faster to decode the SWD file than to resolve the same ambient wave field
   for every case.

 - Reassessing short critical time windows from from a long simulation, using a more advanced solver
   (e.g. a CFD-solver instead of a BEM-solver applied for the long simultion) is complicated if not
   using the interface file approach.

 - No, the SWD file is usually not huge. Typical size for a 20-min simulation of irregular
   seas is 10-100 Mb.


.. _swd-c-stream-why:

Why use a binary C stream format for the SWD file?
--------------------------------------------------

The reasons for this decision are

 - Standard Fortran, C, C++, Python and Matlab support I/O of binary C streams.
 - No additional libraries are required to deal with the
   `spectral_wave_data <https://github.com/SpectralWaveData/spectral_wave_data>`_ API.
 - The file format is simple and I/O is efficiently handled with minimum storage requirements.
 - Nothing wrong with HDF5 and other advanced libraries, some people just don't like the idea of
   compiling and linking very large and complicated libraries for smaller tasks.
   This API should be simple to apply in all
   :ref:`wave generators<wave-generator>` and :ref:`application programs<application-program>`.
   Small is beautiful.

Why don't you implement vector operations?
------------------------------------------

For the following reasons most kinematics is evaluated for individual locations only.

 - The cost of function calls is small compared to kinematic calculations for irregular seas.
 - Where to vectorize or parallelize the code, is highly dependent on the actual application program.
 - Parallel computations are more complicated in multi-language software
 - Simplicity is important in the open source.

You may always optimize the implementation for your specific application as you like,
and still enjoy the benefits of the open API.

Will there be FFT based evaluations?
------------------------------------

Yes, for a set of specific features, like global grid evaluation of surface elevation, there will be
FFT based methods. Different API and implementations are considered.
The first version with such features is expected to be 1.1.0.

For other kinematics the implemented recursive schemes are superior with respect to speed
and accuracy when doing exact evaluation at arbitrary locations.

Why Fortran with C/C++/Python wrappers?
---------------------------------------

This API is comprehensive with many methods supporting many spectral formulations that may use
several numerical implementations.
That is why the initial release only provides one
native implementation and apply wrappers for the other languages.

Fortran was selected because
it is fast and has an ISO standardized interface to C which again has well established interfaces to
C++ and Python. It would be difficult to support Fortran based on other native implementations.


If somebody is willing to make and maintain a complete native implementation in other languages than
Fortran-2008 that is great, please contribute. Not only C/C++/Python but other languages like Matlab
may be relevant.

For your own implementations you may of course apply any language you like and still enjoy
the benefits of the open API.














