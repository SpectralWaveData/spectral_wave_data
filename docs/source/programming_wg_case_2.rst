^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Case 2: The variable depth formulation by Gouin et al. in relation to SWD
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The formulation applied by :cite:`HOSvardepthGouin`
considers long-crested potential flow waves propagating over a variable water depth.
The local water depth :math:`d(x)` is defined as

.. math::
   d(x) = d_0 - \beta(x), \qquad d_0 > 0, \qquad d(x) > 0

where :math:`d_0` is a constant and the spatial variation is prescribed by
the function :math:`\beta(x)`. From equations (8) and (9) in :cite:`HOSvardepthGouin`
the wave potential is split into two components

.. math::
   \phi(\mathbf{x}, t) = \phi_{d_0}(\mathbf{x}, t) + \phi_{\beta}(\mathbf{x}, t)

.. math::
   \phi_{d_0}(\mathbf{x}, t) = \sum_{j=0}^n \mathcal{Re}
     \Bigl\{A_j(t)\,  \frac{\cosh k_j(z+d_0)}{\cosh k_j d_0} e^{i k_j x} \Bigr\}, \qquad
      k_j = j \Delta k_x

.. math::
   \phi_{\beta}(\mathbf{x}, t) = \sum_{j=0}^n \mathcal{Re}
     \Bigl\{B_j(t)\,  \frac{\sinh k_jz}{\cosh k_j d_0} e^{i k_j x} \Bigr\}

.. math::
   \zeta(x, y, t)= \sum_{j=0}^n \mathcal{Re} \Bigl\{H_j(t)\, e^{i k_j x} \Bigr\}

The motivation for this split is that :math:`\phi_{d_0}` represents the classical formulation
if :math:`\beta(x)\equiv 0`. The potential :math:`\phi_{\beta}` has only minor contribution
at the free surface, in fact the potential vanish at :math:`z=0`. At the sea bed the fraction
:math:`\frac{\sinh k_jz}{\cosh k_j d_0}` is close to -1. Hence :math:`\phi_{\beta}` has similar
numerical properties at the sea bed as :math:`\phi_{d_0}` has close to the free surface and vice versa.
Attractive numerical properties using this split are demonstrated in :cite:`HOSvardepthGouin`.

If we compare the above formulation to :doc:`shape 3 <shape_3>` as defined in the theory section of our
documentation:

.. math::
   \phi(\mathbf{x}, t)= \sum_{j=0}^n \mathcal{Re} \Bigl\{c_j(t)\, e^{-i k_j x} \Bigr\} e^{k_j z} +
                  \sum_{j=0}^{\hat{n}} \mathcal{Re}\Bigl\{\hat{c}_j(t)\, e^{-i k_j x} \Bigr\} e^{-k_j z}

.. math::
   \zeta(x, t)= \sum_{j=0}^n \mathcal{Re} \Bigl\{h_j(t)\, e^{-i k_j x} \Bigr\}

we can apply the SWD formulation if and only if the kinematics in the physical space
is identical to the original formulation presented above by :cite:`HOSvardepthGouin`.
To satisfy this requirement it follows from pure algebra that:

.. math::
           c_j(t) = \gamma_j(t) + \hat{\gamma}_j(t), \qquad
           \hat{c}_j(t) = \gamma_j(t)e^{-2k_j d_0} - \hat{\gamma}_j(t), \qquad
           h_j(t) = \bar{H}_j(t)

.. math::
   \hat{n}=n, \qquad   \gamma_j(t) = \frac{1+\tanh k_j d_0}{2}\bar{A}_j(t), \qquad
      \hat{\gamma}_j(t) = \frac{e^{-k_j d_0}}{1 + e^{-2k_j d_0}}\bar{B}_j(t)

:math:`\bar{A}_j(t)`, :math:`\bar{B}_j(t)` and :math:`\bar{H}_j(t)` denote the
complex conjugate of :math:`A_j(t)`, :math:`B_j(t)` and :math:`H_j(t)` respectively.
The conjugate account for the change of sign present in the exponential terms related to :math:`x`.

Very similar numerical properties are observed for the SWD formulation. For large wave
numbers the main contribution to :math:`c_j(t)` derives from :math:`A_j(t)`, while the main
contribution to :math:`\hat{c}_j(t)` derives from :math:`B_j(t)`. At wave crests
the main contribution derives from :math:`c_j(t)` and at the sea floor from :math:`\hat{c}_j(t)`.

When the wave generator has established :math:`A_j(t)`, :math:`B_j(t)` and :math:`H_j(t)`
at the end of each time step, the program may in an SWD output routine calculate
:math:`c_j(t)`, :math:`\hat{c}_j(t)` and :math:`h_j(t)` from above equations.
Then write the SWD specific amplitudes to the :doc:`SWD file <swd_format>`.

In addition SWD requires output of :math:`\frac{dc_j(t)}{dt}`, :math:`\frac{d\hat{c}_j(t)}{dt}`
and :math:`\frac{dh_j(t)}{dt}`. It follows from differentiation that

.. math::
           \frac{dc_j(t)}{dt} = \frac{d\gamma_j(t)}{dt} + \frac{d\hat{\gamma}_j(t)}{dt}, \qquad
           \frac{d\hat{c}_j(t)}{dt} = e^{-2k_j d_0}\frac{d\gamma_j(t)}{dt} - \frac{d\hat{\gamma}_j(t)}{dt}

.. math::
           \frac{d\gamma_j(t)}{dt} = \frac{1+\tanh k_j d_0}{2}\frac{d\bar{A}_j(t)}{dt}, \qquad
           \frac{d\hat{\gamma}_j(t)}{dt} = \frac{e^{-k_j d_0}}{1 + e^{-2k_j d_0}}\frac{d\bar{B}_j(t)}{dt}

.. math::
           \frac{dh_j(t)}{dt} = \frac{d\bar{H}_j(t)}{dt}

It may be tempting to evaluate :math:`\frac{dA_j(t)}{dt}`, :math:`\frac{dB_j(t)}{dt}` and
:math:`\frac{dH_j(t)}{dt}` using some kind
of finite-difference schemes. However, as pointed
out in the SWD theory section, higher accuracy is obtained if these expressions are evaluated
from the dynamic and kinematic free surface conditions. This is because these expressions
can be evaluated using analytical spatial gradients as a replacement for the temporal gradient.

In practice these terms can often be extracted
from the right-hand side of the ODE equation system. Consequently, the reconstructed flow
field in SWD, satisfies the free surface boundary conditions to a higher accuracy.

A copy of the subroutine :doc:`swd_write_shape_3.f90 <fortran_swd_write_shape_3>`
may be adapted for convenient SWD output. This code takes care of the low level C-stream output.

.. note::
  It was assumed that the wave train propagates in the positive x-direction
  in the original formulation, the same convention as assumed in SWD.
  If the wave moves in the negative direction, it is in general **not possible**
  to just flip the wave propagation direction due to the varying sea floor, unless
  the sequence of sea floor offset points are reversed when writing the SWD file.

.. note::
  For most sea floor topologies, it is mainly the low frequency components of :math:`\hat{c}_j(t)`
  or (:math:`B_j(t)`) that contributes to the flow field. Consequently, the wave generator
  may specify :math:`\hat{n}<n` to reduce the computational cost.

.. bibliography:: swd_references.bib
   :filter: docname in docnames
