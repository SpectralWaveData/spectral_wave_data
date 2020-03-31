Shape class 5
-------------

This shape class describes general short crested spectral waves propagating in constant water depth :math:`d`.

.. math::
   \phi(x,y,z,t) = \sum_{j_x=0}^{n_x}\sum_{j_y=-n_y}^{n_y}
                  \mathcal{Re} \Bigl\{c_{j_y,j_x}(t)\, X_{j_x}(x)\,Y_{j_y}(y)\Bigr\}\, Z_{j_y,j_x}(z)

.. math::
  \zeta(x,y,t) = \sum_{j_x=0}^{n_x}\sum_{j_y=-n_y}^{n_y}
                  \mathcal{Re} \Bigl\{h_{j_y,j_x}(t) \,X_{j_x}(x)\,Y_{j_y}(y)\Bigr\}

.. math::
  X_{j_x}(x) = e^{- i k_{j_x} x}, \quad
  Y_{j_y}(y) = e^{- i k_{j_y} y}, \quad
  Z_{j_y,j_x}(z) = \frac{\cosh k_{j_y,j_x}(z+d)}{\cosh k_{j_y,j_x} d}

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

  A special implementation optimized for the case :math:`\Delta k_x=\Delta k_y` and
  :math:`n_x = n_y` should be implemented. Such a class will improve the computational speed.

Kinematics
^^^^^^^^^^

In order to exploit some symmetric properties in the :math:`y`-direction of this class
we apply the following equivalent formulation for more efficient evaluations.

.. math::
   \phi(x,y,z,t) = \sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
                  \mathcal{Re} \Bigl\{C_{1,j_y,j_x}(y, t)\, X_{j_x}(x)\Bigr\}\, Z_{j_y,j_x}(z)

.. math::
  \zeta(x,y,t) = \sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
                  \mathcal{Re} \Bigl\{H_{1, j_y,j_x}(y, t) \,X_{j_x}(x)\Bigr\}

.. math::
   C_{1,j_y,j_x}(y, t) = c_{1,j_y,j_x}(t)\,Y_{j_y}(y) +
                         c_{2,j_y,j_x}(t)\, \bar{Y}_{j_y}(y)

.. math::
   H_{1,j_y,j_x}(y, t) = h_{1,j_y,j_x}(t)\,Y_{j_y}(y) +
                         h_{2,j_y,j_x}(t)\, \bar{Y}_{j_y}(y)

.. math::
   c_{1,j_y,j_x}(t) = c_{j_y,j_x}(t), \qquad
   c_{2,j_y,j_x}(t) = \begin{cases}
                                      c_{-j_y,j_x}(t),  & \text{$j_y>0$}, \\
                                      0,                & \text{$j_y=0$}
                      \end{cases} \\
   h_{1,j_y,j_x}(t) = h_{j_y,j_x}(t), \qquad
   h_{2,j_y,j_x}(t) = \begin{cases}
                                      h_{-j_y,j_x}(t),  & \text{$j_y>0$}, \\
                                      0,                & \text{$j_y=0$}
                      \end{cases}

where :math:`\bar{Y}_{j_y}(y)` denotes the complex conjugate of :math:`Y_{j_y}(y)`.

Given the definitions above we obtain the following explicit kinematics:

.. math::
   \phi(\bar{x},\bar{y},\bar{z},\bar{t})= \sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
              \mathcal{Re} \Bigl\{C_{1,j_y,j_x}(y, t)\, X_{j_x}(x)\Bigr\} Z_{j_y,j_x}(z)

.. math::
   \varphi(\bar{x},\bar{y},\bar{z},\bar{t}) \equiv 0

.. math::
  \frac{\partial\phi}{\partial \bar{t}}(\bar{x},\bar{y},\bar{z},\bar{t}) = \sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
              \mathcal{Re} \Bigl\{\frac{d C_{1,j_y,j_x}(y, t)}{dt} \, X_{j_x}(x)\Bigr\} Z_{j_y,j_x}(z)

.. math::
   \zeta(\bar{x},\bar{y},\bar{t})= \sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
              \mathcal{Re}\Bigl\{H_{1,j_y,j_x}(y, t)\, X_{j_x}(x)\Bigr\}

.. math::
  \frac{\partial\zeta}{\partial \bar{t}}(\bar{x},\bar{y},\bar{t}) = \sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
              \mathcal{Re} \Bigl\{\frac{d H_{1,j_y,j_x}(y, t)}{dt} \, X_{j_x}(x)\Bigr\}

.. math::
   \frac{\partial\zeta}{\partial \bar{x}}(\bar{x},\bar{y},\bar{t}) = \zeta_x\cos\beta - \zeta_y\sin\beta, \qquad
   \frac{\partial\zeta}{\partial \bar{y}}(\bar{x},\bar{y},\bar{t}) = \zeta_x\sin\beta + \zeta_y\cos\beta

.. math::
   \zeta_x =\sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
            k_{j_x} \mathcal{Im} \Bigl\{H_{1,j_y,j_x}(y, t)\, X_{j_x}(x)\Bigr\}

.. math::
   \zeta_y = -\sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
            k_{j_y} \mathcal{Im} \Bigl\{H_{2,j_y,j_x}(y, t)\, X_{j_x}(x)\Bigr\}

.. math::
   \bar{\nabla}\phi(\bar{x},\bar{y},\bar{z},\bar{t}) =
              [\phi_x\cos\beta - \phi_y\sin\beta, \phi_x\sin\beta + \phi_y\cos\beta,\phi_z]^T

.. math::
   \phi_x = \sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
          k_{j_x}\mathcal{Im} \Bigl\{C_{1,j_y,j_x}(y, t)\, X_{j_x}(x)\Bigr\} \, Z_{j_y,j_x}(z)

.. math::
   \phi_y = \sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
          k_{j_y}\mathcal{Im} \Bigl\{C_{2,j_y,j_x}(y, t)\, X_{j_x}(x)\Bigr\} \, Z_{j_y,j_x}(z)

.. math::
   \phi_z = \sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
          \mathcal{Re} \Bigl\{C_{1,j_y,j_x}(y, t)\, X_{j_x}(x)\Bigr\} \, \frac{d Z_{j_y,j_x}(z)}{dz}

.. math::
  \frac{\partial\bar{\nabla}\phi}{\partial \bar{t}}(\bar{x},\bar{y},\bar{z},\bar{t}) =
              [\phi_{xt}\cos\beta - \phi_{yt}\sin\beta, \phi_{xt}\sin\beta + \phi_{yt}\cos\beta,\phi_z]^T

.. math::
   \phi_{xt} = \sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
          k_{j_x}\mathcal{Im} \Bigl\{\frac{d C_{1,j_y,j_x}(y, t)}{dt}\, X_{j_x}(x)\Bigr\} \, Z_{j_y,j_x}(z)

.. math::
   \phi_{yt} = \sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
          k_{j_y}\mathcal{Im} \Bigl\{\frac{d C_{2,j_y,j_x}(y, t)}{dt}\, X_{j_x}(x)\Bigr\} \, Z_{j_y,j_x}(z)

.. math::
   \phi_{zt} = \sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
          \mathcal{Re} \Bigl\{\frac{d C_{1,j_y,j_x}(y, t)}{dt}\, X_{j_x}(x)\Bigr\} \, \frac{d Z_{j_y,j_x}(z)}{dz}

.. math::
  \frac{d\bar{\nabla}\phi}{d\bar{t}}(\bar{x},\bar{y},\bar{z},\bar{t}) =
           \frac{\partial\bar{\nabla}\phi}{\partial \bar{t}} +
  \bar{\nabla}\phi \cdot \bar{\nabla}\bar{\nabla}\phi

.. math::
   \bar{\nabla}\bar{\nabla}\phi (\bar{x},\bar{y},\bar{z},\bar{t}) =
     \begin{bmatrix}
       \phi_{\bar{x},\bar{x}}  & \phi_{\bar{x},\bar{y}} & \phi_{\bar{x},\bar{z}} \\
       \phi_{\bar{x},\bar{y}}  & \phi_{\bar{y},\bar{y}} & \phi_{\bar{y},\bar{z}} \\
       \phi_{\bar{x},\bar{z}}  & \phi_{\bar{y},\bar{z}} & \phi_{\bar{z},\bar{z}}
     \end{bmatrix}

.. math::
   \phi_{\bar{x},\bar{x}} = \phi_{xx}\cos^2\beta - \phi_{xy}\sin(2\beta) + \phi_{yy}\sin^2\beta

.. math::
   \phi_{\bar{x},\bar{y}} = \phi_{xy}(\cos^2\beta - \sin^2\beta) + (\phi_{xx} - \phi_{yy})\sin\beta\cos\beta

.. math::
   \phi_{\bar{x},\bar{z}} = \phi_{xz}\cos\beta - \phi_{yz}\sin\beta

.. math::
   \phi_{\bar{y},\bar{y}} = \phi_{yy}\cos^2\beta + \phi_{xy}\sin(2\beta) + \phi_{xx}\sin^2\beta

.. math::
   \phi_{\bar{y},\bar{z}} = \phi_{yz}\cos\beta + \phi_{xz}\sin\beta

.. math::
   \phi_{\bar{z},\bar{z}} = \phi_{zz} = -\phi_{xx} -\phi_{yy}

.. math::
   \phi_{xx} = - \sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
       k_{j_x}^2 \mathcal{Re} \Bigl\{C_{1,j_y,j_x}(y, t) \, X_{j_x}(x)\Bigr\} Z_{j_y,j_x}(z)

.. math::
   \phi_{xy} = - \sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
       k_{j_x} k_{j_y} \mathcal{Re} \Bigl\{C_{2,j_y,j_x}(y, t) \, X_{j_x}(x)\Bigr\} Z_{j_y,j_x}(z)

.. math::
   \phi_{xz} = \sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
       k_{j_x} \mathcal{Im} \Bigl\{C_{1,j_y,j_x}(y, t) \, X_{j_x}(x)\Bigr\} \frac{d Z_{j_y,j_x}(z)}{dz}

.. math::
   \phi_{yy} = - \sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
       k_{j_y}^2 \mathcal{Re} \Bigl\{C_{1,j_y,j_x}(y, t) \, X_{j_x}(x)\Bigr\} Z_{j_y,j_x}(z)

.. math::
   \phi_{yz} = \sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
       k_{j_y} \mathcal{Im} \Bigl\{C_{2,j_y,j_x}(y, t) \, X_{j_x}(x)\Bigr\} \frac{d Z_{j_y,j_x}(z)}{dz}

.. math::
   \phi_{zz} = \sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
       k_{j_y,j_x}^2 \mathcal{Re} \Bigl\{C_{1,j_y,j_x}(y, t) \, X_{j_x}(x)\Bigr\} Z_{j_y,j_x}(z)
       = -\phi_{xx} - \phi_{yy}

.. math::
   \frac{\partial^2\zeta}{\partial \bar{x}^2}(\bar{x},\bar{y},\bar{t}) =
      \zeta_{xx}\cos^2\beta - \zeta_{xy}\sin(2\beta) + \zeta_{yy}\sin^2\beta

.. math::
   \frac{\partial^2\zeta}{\partial\bar{x}\partial\bar{y}}(\bar{x},\bar{y},\bar{t}) =
      \zeta_{xy}(\cos^2\beta - \sin^2\beta) + (\zeta_{xx} - \zeta_{yy})\sin\beta\cos\beta

.. math::
   \frac{\partial^2\zeta}{\partial\bar{y}^2}(\bar{x},\bar{y},\bar{t}) =
      \zeta_{yy}\cos^2\beta  + \zeta_{xy}\sin(2\beta) + \zeta_{xx}\sin^2\beta

.. math::
   \zeta_{xx} = -\sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
      k_{j_x}^2 \mathcal{Re} \Bigl\{H_{1,j_y,j_x}(y, t) \, X_{j_x}(x)\Bigr\}

.. math::
   \zeta_{xy} = -\sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
      k_{j_x} k_{j_y} \mathcal{Re} \Bigl\{H_{2,j_y,j_x}(y, t) \, X_{j_x}(x)\Bigr\}

.. math::
   \zeta_{yy} = -\sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
      k_{j_y}^2 \mathcal{Re} \Bigl\{H_{1,j_y,j_x}(y, t) \, X_{j_x}(x)\Bigr\}

.. math::
   p = -\rho\frac{\partial\phi}{\partial \bar{t}}
       -\frac{1}{2}\rho\bar{\nabla}\phi\cdot\bar{\nabla}\phi
       -\rho g \bar{z}

where :math:`\bar{\nabla}` denotes gradients with respect to
:math:`\bar{x}`, :math:`\bar{y}` and :math:`\bar{z}`.
:math:`\mathcal{Im}\{\alpha\}` denotes the imaginary part of a complex number :math:`\alpha` and

.. math::
  C_{2, j_y, j_x}(y,t) = c_{1, j_y, j_x}(t)\,Y_{j_y}(y) - c_{2, j_y, j_x}(t)\,\bar{Y}_{j_y}(y)

.. math::
  H_{2, j_y, j_x}(y,t) = h_{1, j_y, j_x}(t)\,Y_{j_y}(y) - h_{2, j_y, j_x}(t)\,\bar{Y}_{j_y}(y)




The particle acceleration is labeled :math:`\frac{d\bar{\nabla}\phi}{d\bar{t}}`.

The stream function :math:`\varphi` is not relevant for short crested seas.
Hence, we apply the dummy definition :math:`\varphi=0` for all locations.



Implementation notes
^^^^^^^^^^^^^^^^^^^^

Evaluation of costly transcendental functions (:math:`\cos`, :math:`\sin`, :math:`\exp`, :math:`\cosh`, ...)
is significantly reduced by exploiting the following recursive relations

.. math::
   X_{j_x}(x) = X_1(x)\, X_{j_x-1}(x), \qquad
   Y_{j_y}(y) = Y_1(y)\, Y_{j_y-1}(y)

.. math::

   Z_{j_y,j_x}(z) = U_{j_y,j_x}S_{j_y,j_x} +  V_{j_y,j_x}T_{j_y,j_x}, \qquad
   \frac{dZ_{j_y,j_x}(z)}{dz} = k_{j_y,j_x}(U_{j_y,j_x}S_{j_y,j_x} - V_{j_y,j_x}T_{j_y,j_x})

.. math::
  U_{j_y,j_x} = \frac{1+R_{j_y,j_x}}{2}, \qquad  V_{j_y,j_x}  = 1- U_{j_y,j_x}

.. math::
  S_{j_y,j_x} \equiv e^{k_{j_y,j_x} z}, \quad
  T_{j_y,j_x} \equiv e^{-k_{j_y,j_x} z} = 1/S_{j_y,j_x}, \quad
  R_{j_y,j_x} \equiv \tanh k_{j_y,j_x}d

It should be noted that contrary to long crested seas,
there are no trivial recursive relations for the :math:`z`-dependent terms
:math:`S_{j_y,j_x}`, :math:`T_{j_y,j_x}` and :math:`R_{j_y,j_x}`.
This makes calculations of surface elevations significantly faster than calculations
of other kinematics for short crested seas.

In case the :ref:`wave generator<wave-generator>`
applies a perturbation theory of
order :math:`q` we apply the following Taylor expansion above the calm free surface.


.. math::
   S_{j_y, j_x}(z) = 1 + \sum_{p=1}^{q-1}\frac{(k_{j_y, j_x} z)^p}{p!}, \qquad z > 0
