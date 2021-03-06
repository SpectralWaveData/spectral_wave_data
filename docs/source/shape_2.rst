Shape class 2
-------------

This shape class describes long crested waves propagating in constant water depth :math:`d`.

.. math::
   \phi(x, z, t)= \sum_{j=0}^n \mathcal{Re} \Bigl\{c_j(t)\, X_j(x) \Bigr\} Z_j(z)

.. math::
   \zeta(x, t)= \sum_{j=0}^n \mathcal{Re} \Bigl\{h_j(t)\, X_j(x) \Bigr\}

.. math::
   X_j(x) = e^{-i k_j x}, \quad Z_j(z) = \frac{\cosh k_j(z+d)}{\cosh k_j d}, \quad k_j = j\cdot\Delta k, \quad i=\sqrt{-1}

The set of real constants :math:`k_j` resemble wave numbers. It follows that the
kinematics is periodic in space

.. math::
   \phi(x + \lambda_{\max}, z, t) = \phi(x, z, t), \qquad
   \zeta(x + \lambda_{\max}, t) = \zeta(x, t)

.. math::
   \lambda_{\max} = \frac{2\pi}{\Delta k}, \qquad \lambda_{\min} = \frac{\lambda_{\max}}{n}

where :math:`\lambda_{\min}` and :math:`\lambda_{\max}` are the shortest and longest
wave lengths resolved respectively.

The actual set of shape functions is uniquely defined by the three input parameters
:math:`\Delta k` and :math:`n` and :math:`d`.

.. note::

  The fields related to :math:`j=0` are uniform in space (DC bias). Non-zero values of
  :math:`h_0(t)` violates mass conservation. The amplitude :math:`c_0(t)`
  adds a uniform time varying ambient pressure field not influencing the flow field.
  Consequently, these components will by default be suppressed in the kinematic
  calculations. However, there is an option in the API for including all DC values
  provided by the :ref:`wave generator<wave-generator>`.

  The fields related to :math:`j=n` are expected to correspond to the Nyquist
  frequency of the physical resolution applied in the
  :ref:`wave generator<wave-generator>`.
  Hence, typical :math:`n=\lfloor n_{fft}/2 \rfloor` where :math:`n_{fft}` is the physical
  spatial resolution applied in the :ref:`wave generator<wave-generator>`.

Kinematics
^^^^^^^^^^

Given the definitions above we obtain the following explicit kinematics:

.. math::
   \phi(\bar{x},\bar{y},\bar{z},\bar{t})= \sum_{j=0}^n \mathcal{Re} \Bigl\{c_j(t)\, X_j(x)\Bigr\} Z_j(z)

.. math::
   \varphi(\bar{x},\bar{y},\bar{z},\bar{t})= \sum_{j=0}^n \mathcal{Im} \Bigl\{c_j(t)\, X_j(x)\Bigr\} \hat{Z}_j(z)

.. math::
  \frac{\partial\phi}{\partial \bar{t}}(\bar{x},\bar{y},\bar{z},\bar{t}) = \sum_{j=0}^n \mathcal{Re}
               \Bigl\{\frac{d c_j(t)}{dt} \, X_j(x)\Bigr\} Z_j(z)

.. math::
   \zeta(\bar{x},\bar{y},\bar{t})= \sum_{j=0}^n \mathcal{Re}
               \Bigl\{h_j(t)\, X_j(x)\Bigr\}

.. math::
  \frac{\partial\zeta}{\partial \bar{t}}(\bar{x},\bar{y},\bar{t}) = \sum_{j=0}^n \mathcal{Re}
               \Bigl\{\frac{d h_j(t)}{dt} \, X_j(x)\Bigr\}

.. math::
   \frac{\partial\zeta}{\partial \bar{x}}(\bar{x},\bar{y},\bar{t}) = \zeta_x\cos\beta, \qquad
   \frac{\partial\zeta}{\partial \bar{y}}(\bar{x},\bar{y},\bar{t}) = \zeta_x\sin\beta

.. math::
   \zeta_x = \sum_{j=0}^n k_j\mathcal{Im} \Bigl\{h_j(t)\, X_j(x)\Bigr\}

.. math::
   \bar{\nabla}\phi(\bar{x},\bar{y},\bar{z},\bar{t}) = [\phi_x\cos\beta,\phi_x\sin\beta,\phi_z]^T

.. math::
   \phi_x = \sum_{j=0}^n k_j\mathcal{Im} \Bigl\{c_j(t)\, X_j(x)\Bigr\} \, Z_j(z)

.. math::
   \phi_z = \sum_{j=0}^n \mathcal{Re} \Bigl\{c_j(t)\, X_j(x)\Bigr\} \, \frac{dZ_j(z)}{dz}

.. math::
  \frac{\partial\bar{\nabla}\phi}{\partial \bar{t}}(\bar{x},\bar{y},\bar{z},\bar{t}) =
           [\phi_{xt}\cos\beta,\phi_{xt}\sin\beta,\phi_{zt}]^T

.. math::
   \phi_{xt} = \sum_{j=0}^n k_j \mathcal{Im} \Bigl\{\frac{d c_j(t)}{dt} \, X_j(x)\Bigr\} Z_j(z)

.. math::
   \phi_{zt} = \sum_{j=0}^n \mathcal{Re} \Bigl\{\frac{d c_j(t)}{dt} \, X_j(x)\Bigr\} \frac{dZ_j(z)}{dz}

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
   \phi_{\bar{x},\bar{x}} = \phi_{xx}\cos^2\beta, \qquad
   \phi_{\bar{x},\bar{y}} = \phi_{xx}\sin\beta\cos\beta, \qquad
   \phi_{\bar{x},\bar{z}} = \phi_{xz}\cos\beta

.. math::
   \phi_{\bar{y},\bar{y}} = \phi_{xx}\sin^2\beta, \qquad
   \phi_{\bar{y},\bar{z}} = \phi_{xz}\sin\beta, \qquad
   \phi_{\bar{z},\bar{z}} = \phi_{zz} = -\phi_{xx}

.. math::
   \phi_{xx} = -\sum_{j=0}^n k_j^2 \mathcal{Re} \Bigl\{c_j(t) \, X_j(x)\Bigr\} Z_j(z)

.. math::
   \phi_{zz} = \sum_{j=0}^n k_j^2 \mathcal{Re} \Bigl\{c_j(t) \, X_j(x)\Bigr\} Z_j(z)
   = - \phi_{xx}

.. math::
   \phi_{xz} = \sum_{j=0}^n k_j \mathcal{Im} \Bigl\{c_j(t) \, X_j(x)\Bigr\} \frac{dZ_j(z)}{dz}

.. math::
   \frac{\partial^2\zeta}{\partial \bar{x}^2}(\bar{x},\bar{y},\bar{t}) = \zeta_{xx}\cos^2\beta
   \qquad
   \frac{\partial^2\zeta}{\partial \bar{y}^2}(\bar{x},\bar{y},\bar{t}) = \zeta_{xx}\sin^2\beta

.. math::
   \frac{\partial^2\zeta}{\partial\bar{x}\partial\bar{y}}(\bar{x},\bar{y},\bar{t}) =
                \zeta_{xx}\sin\beta\cos\beta

.. math::
   \zeta_{xx} = -\sum_{j=0}^n k_j^2 \mathcal{Re} \Bigl\{h_j(t) \, X_j(x)\Bigr\}

.. math::
   p = -\rho\frac{\partial\phi}{\partial \bar{t}}
       -\frac{1}{2}\rho\bar{\nabla}\phi\cdot\bar{\nabla}\phi
       -\rho g \bar{z}

where :math:`\bar{\nabla}` denotes gradients with respect to
:math:`\bar{x}`, :math:`\bar{y}` and :math:`\bar{z}`. The particle acceleration
is labeled :math:`\frac{d\bar{\nabla}\phi}{d\bar{t}}`.

The stream function :math:`\varphi` is related to the velocity potential  :math:`\phi`.
Hence :math:`\partial \phi/\partial x = \partial \varphi/\partial z`
and  :math:`\partial \phi/\partial z = -\partial \varphi/\partial x`.
Note that for the stream function evaluation we apply the function

.. math::
    \hat{Z}_j(z) = \frac{\sinh k_j(z+d)}{\cosh k_j d}

Implementation notes
^^^^^^^^^^^^^^^^^^^^

Evaluation of costly transcendental functions (:math:`\cos`, :math:`\sin`, :math:`\exp`, ...)
are almost eliminated by exploiting the following recursive relations

.. math::
   X_j(x) = X_1(x)\, X_{j-1}(x), \qquad
   Z_j(z) = U_jS_j +  V_jT_j, \qquad
   \frac{dZ_j(z)}{dz} = k_j(U_jS_j - V_jT_j)

.. math::
   U_j = \frac{1+R_j}{2}, \qquad  V_j  = 1 - U_j, \qquad
   R_j \equiv \tanh k_jd = \frac{R_1 + R_{j-1}}{1 + R_1\,R_{j-1}}

.. math::
  S_j \equiv e^{k_j z} = S_1 S_{j-1}, \qquad   T_j \equiv e^{-k_j z} = T_1 T_{j-1} = 1/S_j, \qquad
  \hat{Z}_j(z) = U_jS_j - V_jT_j

When :math:`1 - R_j < 100 \epsilon_m`, where :math:`\epsilon_m` is the machine precision, we
apply the deep water approximations :math:`Z_j=S_j` and :math:`\frac{dZ_j(z)}{dz} = k_j S_j`.
Consequently, only the first :math:`\hat{j}` spectral components should be treated as shallow
water components

.. math::
  \hat{j} = \operatorname{int}\large(\frac{\tanh^{-1}(1-100\epsilon_m)}{\Delta k\cdot d}\large)

This number can be precomputed prior to any kinematical calculations to avoid testing and numerical issues.

In case the :ref:`wave generator<wave-generator>` applies a perturbation theory of
order :math:`q` we apply the following Taylor expansion above the calm free surface.


.. math::
   S_j(z) = 1 + \sum_{p=1}^{q-1}\frac{(k_j z)^p}{p!}, \qquad z > 0
