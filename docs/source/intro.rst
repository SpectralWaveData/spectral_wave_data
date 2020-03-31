************
Introduction
************

The open-source GitHub organization `SpectralWaveData <https://github.com/SpectralWaveData>`_ hosts several
repositories related to the `spectral_wave_data <https://github.com/SpectralWaveData/spectral_wave_data>`_ API.

The main purpose of the **spectral_wave_data** API is to establish an open interface for how to exchange
spectral ocean wave kinematics between computer programs.

Calculation of ocean wave kinematics is not only relevant for environmental assessments, but is
also highly important for predicting wave induced motions and loads on marine structures.

In lack of proper tools, knowledge and computer capacity, realistic critical wave trains have typical
been replaced by idealized regular waves or linear sea states, when assessing responses of marine structures.
However, extensive recent research documents
important physical effects not present in these simplistic wave models.

The simplistic wave models can be described by a small set of parameters.
Consequently, it has not been important to have standardized field interfaces between programs.
Ad-hoc input transformations, like adjustment of phase definitions, were sufficient prior to actual calculations.

The more realistic critical wave trains are complicated to describe.

 - What are the important characteristics of such wave trains?
 - What is the appropriate mathematical model(s)?
 - How should my software deal with such wave fields?

These questions are still challenging, and the answers may differ depending on the
considered phenomena.
However, the **spectral_wave_data** API is an attempt to push forward in this respect.
Not to provide the answers to all these questions, but to facilitate such investigations.

