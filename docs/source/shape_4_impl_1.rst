Shape class 4, impl=1
---------------------

This implementation of :doc:`Shape 4 <shape_4>` evaluates the more general
non-symmetric spatial resolution where :math:`\Delta k_x` may differ from
:math:`\Delta k_y` or :math:`n_x` from :math:`n_y`.

However, since :math:`Y_{-j_y}(y) = \bar{Y}_{j_y}(y)`,
where :math:`\bar{Y}_{j_y}(y)` denotes the complex conjugate of :math:`Y_{j_y}(y)`,
and :math:`Z_{-j_y,j_x}(z) = Z_{j_y,j_x}(z)` we may always apply the following
equivalent wave field formulation.

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


Kinematics
^^^^^^^^^^

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
          k_{j_y,j_x} \mathcal{Re} \Bigl\{C_{1,j_y,j_x}(y, t)\, X_{j_x}(x)\Bigr\} \, Z_{j_y,j_x}(z)

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
          k_{j_y,j_x} \mathcal{Re} \Bigl\{\frac{d C_{1,j_y,j_x}(y, t)}{dt}\, X_{j_x}(x)\Bigr\} \, Z_{j_y,j_x}(z)

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
       k_{j_x} k_{j_y,j_x} \mathcal{Im} \Bigl\{C_{1,j_y,j_x}(y, t) \, X_{j_x}(x)\Bigr\} Z_{j_y,j_x}(z)

.. math::
   \phi_{yy} = - \sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
       k_{j_y}^2 \mathcal{Re} \Bigl\{C_{1,j_y,j_x}(y, t) \, X_{j_x}(x)\Bigr\} Z_{j_y,j_x}(z)

.. math::
   \phi_{yz} = \sum_{j_x=0}^{n_x}\sum_{j_y=0}^{n_y}
       k_{j_y} k_{j_y,j_x} \mathcal{Im} \Bigl\{C_{2,j_y,j_x}(y, t) \, X_{j_x}(x)\Bigr\} Z_{j_y,j_x}(z)

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
We also apply

.. math::
  C_{2, j_y, j_x}(y,t) = c_{1, j_y, j_x}(t)\,Y_{j_y}(y) - c_{2, j_y, j_x}(t)\,\bar{Y}_{j_y}(y)

.. math::
  H_{2, j_y, j_x}(y,t) = h_{1, j_y, j_x}(t)\,Y_{j_y}(y) - h_{2, j_y, j_x}(t)\,\bar{Y}_{j_y}(y)




The particle acceleration is labeled :math:`\frac{d\bar{\nabla}\phi}{d\bar{t}}`.

The stream function :math:`\varphi` is not relevant for short crested seas.
Hence, we apply the dummy definition :math:`\varphi=0` for all locations.



Implementation notes
^^^^^^^^^^^^^^^^^^^^

Evaluation of costly transcendental functions (:math:`\cos`, :math:`\sin`, :math:`\exp`, ...)
is significantly reduced by exploiting the following recursive relations

.. math::
   X_{j_x}(x) = X_1(x)\, X_{j_x-1}(x), \qquad
   Y_{j_y}(y) = Y_1(y)\, Y_{j_y-1}(y)

It should be noted that contrary to long crested seas,
there are no trivial recursive relations for the :math:`z`-dependent term
:math:`Z_{j_y,j_x}(z)`.
This makes calculations of surface elevations significantly faster than calculations
of other kinematics for short crested seas.

In case the :ref:`wave generator<wave-generator>`
applies a perturbation theory of
order :math:`q` we apply the following Taylor expansion above the calm free surface.


.. math::
   Z_{j_y, j_x}(z) = 1 + \sum_{p=1}^{q-1}\frac{(k_{j_y, j_x} z)^p}{p!}, \qquad z > 0
