^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Case 1: The HOS-Ocean formulation in relation to SWD
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The formulation applied in `HOS-Ocean <https://github.com/LHEEA/HOS-ocean/wiki>`_,
for long-crested waves propagating over a constant water depth  :math:`d`,
is explained in :cite:`HOSOceanOrg`. From equations (6) and (8) in that publication,
the velocity potential :math:`\phi` and surface elevation :math:`\zeta` are defined as

.. math::
   \phi(\mathbf{x},t)= \sum_{j=0}^n \mathcal{Re}
     \Bigl\{A_j(t)\,  \frac{\cosh k_j(z+d)}{\cosh k_j d} e^{i k_j x} \Bigr\}, \qquad
      k_j = j \Delta k_x

.. math::
   \zeta(x, y, t)= \sum_{j=0}^n \mathcal{Re} \Bigl\{B_j^{\zeta}(t)\, e^{i k_j x} \Bigr\}

We have adjusted the symbols for wave elevation, water depth and loop index
to simplify the comparison with the SWD documentation.
When `HOS-Ocean` has calculated the temporal spectral amplitudes :math:`A_j(t)`
and :math:`B_j^{\zeta}(t)` the kinematics can be calculated everywhere using
the above equations.

We compare the formulation above to :doc:`shape 2 <shape_2>` as defined in the theory
section of our documentation:

.. math::
   \phi(\mathbf{x}, t)= \sum_{j=0}^n \mathcal{Re} \Bigl\{c_j(t)\, e^{-i k_j x} \Bigr\}
           \frac{\cosh k_j(z+d)}{\cosh k_j d}

.. math::
   \zeta(x, y, t)= \sum_{j=0}^n \mathcal{Re} \Bigl\{h_j(t)\, e^{-i k_j x} \Bigr\}

It is noted that the sign in the exponential terms are different in theses formulations.
We can apply the SWD formulation if and only if the kinematics in the physical space is
identical to the `HOS-Ocean` formulation.
To satisfy this requirement it follows from pure algebra that:

.. math::
           c_j(t) = \bar{A}_j(t), \qquad h_j(t) = \bar{B}_j^{\zeta}(t)

:math:`\bar{A}_j(t)` and :math:`\bar{B}_j^{\zeta}(t)` denote the complex conjugate of
:math:`A_j(t)` and :math:`B_j^{\zeta}(t)` respectively.

Hence, when `HOS-Ocean` has established :math:`A_j(t)` and :math:`B_j^{\zeta}(t)`
at the end of each time step,
the program may in an SWD output routine calculate :math:`c_j(t)` and :math:`h_j(t)`
from above equations and write these SWD specific amplitudes to the
:doc:`SWD-file <swd_format>`.

In addition SWD requires output of :math:`\frac{dc_j(t)}{dt}` and :math:`\frac{dh_j(t)}{dt}`.
It follows from differentiation that

.. math::
           \frac{dc_j(t)}{dt} = \frac{d\bar{A}_j(t)}{dt}, \qquad
           \frac{dh_j(t)}{dt} = \frac{d\bar{B}_j^{\zeta}(t)}{dt}

It may be tempting to evaluate the right-hand side of these equations using some kind
of finite-difference schemes. However, as pointed
out in the SWD theory section, higher accuracy is obtained if these expressions are evaluated
from the dynamic and kinematic free surface conditions. This is because these expressions
can be evaluated using analytical spatial gradients as a replacement for the temporal gradient.
In practice these terms can often be extracted from the right-hand side of the ODE equation system.

Noting that the physical boundary conditions explicitly define the temporal slopes
of these amplitudes, this is also important information regarding temporal interpolation
to be applied in the application program. This is accounted for in the SWD API assuming
sound input of :math:`\frac{dc_j(t)}{dt}` and :math:`\frac{dh_j(t)}{dt}`.

In practice, a copy of the subroutine
:doc:`swd_write_shape_1_or_2.f90 <fortran_swd_write_shape_1_or_2>`
may be adapted for convenient SWD output. This code takes care of the low level C-stream output.

.. note::
  It was assumed that the wave train propagates in the positive x-direction
  in `HOS-Ocean`, the same convention as assumed in SWD. If the wave moves in the negative
  direction in `HOS-Ocean`, we must flip the sign of the x-axis to comply with the
  SWD convention. Hence, we then must apply
  :math:`c_j(t)=A_j(t)` and :math:`h_j(t)=B_j^{\zeta}(t)` as a result of conjugating twice.



.. bibliography:: swd_references.bib
   :filter: docname in docnames
