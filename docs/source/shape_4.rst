Shape class 4
-------------

This shape class describes general short crested spectral waves propagating in infinite water depth.

.. math::
   \phi(x,y,z,t) = \sum_{j_x=0}^{n_x}\sum_{j_y=-n_y}^{n_y}
                  \mathcal{Re} \Bigl\{c_{j_y,j_x}(t)\, X_{j_x}(x)\,Y_{j_y}(y)\Bigr\}\, Z_{j_y,j_x}(z)

.. math::
  \zeta(x,y,t) = \sum_{j_x=0}^{n_x}\sum_{j_y=-n_y}^{n_y}
                  \mathcal{Re} \Bigl\{h_{j_y,j_x}(t) \,X_{j_x}(x)\,Y_{j_y}(y)\Bigr\}

.. math::
  X_{j_x}(x) = e^{- i k_{j_x} x}, \quad
  Y_{j_y}(y) = e^{- i k_{j_y} y}, \quad
  Z_{j_y,j_x}(z) = e^{k_{j_y,j_x} z}

.. math::
  k_{j_x} =  j_x\cdot\Delta k_x, \quad k_{j_y} = j_y\cdot\Delta k_y, \quad
  k_{j_y,j_x} = \sqrt{k_{j_x}^2+k_{j_y}^2}, \quad i = \sqrt{-1}

The set of real constants :math:`k_{j_x}` and :math:`k_{j_y}` resemble wave numbers
in the :math:`x` and :math:`y` directions respectively.
It follows that the kinematics is periodic in space

.. math::
   \phi(x + L_x, y + L_y, z, t) = \phi(x, y, z, t), \qquad
   \zeta(x + L_x, y + L_y, t) = \zeta(x, y, t)

.. math::
   L_x = \frac{2\pi}{\Delta k_x}, \qquad L_y = \frac{2\pi}{\Delta k_y}, \qquad
   \lambda_{\min} = \frac{2\pi}{\sqrt{(n_x \Delta k_x)^2+(n_y \Delta k_y)^2}}

where :math:`\lambda_{\min}` is the shortest wave lengths resolved.
The actual set of shape functions is uniquely defined by the five input parameters
:math:`\Delta k_x`, :math:`\Delta k_y`, :math:`n_x`, :math:`n_y` and :math:`d`.

.. note::

  The fields related to :math:`j_x=j_y=0` are uniform in space (DC bias). Non-zero values of
  :math:`h_{0,0}(t)` violates mass conservation. The amplitude :math:`c_{0,0}(t)`
  adds a uniform time varying ambient pressure field not influencing the flow field.
  Consequently, these components will by default be suppressed in the kinematic
  calculations. However, there is an option in the API for including all DC values
  provided by the :ref:`wave generator<wave-generator>`.

  The fields related to :math:`j_x=n_x` and :math:`j_y=\pm n_y` are expected to correspond
  to the Nyquist frequencies of the physical resolution applied in the
  :ref:`wave generator<wave-generator>`.
  Hence, typical :math:`n_x=\lfloor n_{x,fft}/2 \rfloor` and
  :math:`n_y=\lfloor n_{y,fft}/2 \rfloor` where :math:`n_{x,fft}` and
  :math:`n_{y,fft}` are the physical spatial resolutions applied in the
  :ref:`wave generator<wave-generator>`, in the :math:`x` and :math:`y`
  directions respectively.


Evaluation of kinematics in short-crested seas is in general computational demanding.
Consequently, this API provides several alternative implementations in order
to exploit eventual symmetric properties or numerical approximations.

 - :doc:`Shape 4, impl 1: <shape_4_impl_1>` A general implementation.
 - :doc:`Shape 4, impl 2: <shape_4_impl_2>` Optimized for symmetric resolution
   :math:`\Delta k_x=\Delta k_y` and :math:`n_x=n_y`.

.. toctree::
   :hidden:
   :caption: Shapes 4, Implementations

   shape_4_impl_1
   shape_4_impl_2



