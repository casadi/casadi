Inverted Pendulum
=================

In this example, a mode predictive controller (MPC) is used to stabilize an
inverted pendulum mounted on a moving cart.

.. raw:: html
   :file: ../../sphinxstatic/inverted-pendulum.html

.. raw:: html

    <script>
    window.addEventListener('load',
        function() {
            document.querySelectorAll(".anim-buttons > button:nth-child(6)")
                .forEach(btn => { btn.click(); });
        }, false);
    </script>

.. figure:: ../../sphinxstatic/inverted-pendulum-graphs.svg
   :scale: 8%
   :align: center
   :alt: Figure of the MPC solution to the inverted pendulum problem.

   Plot of the states and the control signal of the MPC solution.

The state vector consist of the angle of the pendulum :math:`\theta`, its 
angular velocity :math:`\omega`, the position of the cart :math:`x`, and its
velocity :math:`v`. We use the following model:

.. math::

    \newcommand{\pend}[1]{#1_\mathrm{pend}}
    \newcommand{\cart}[1]{#1_\mathrm{cart}}
    \begin{aligned}
        \pend{a} &= \pend{l} \sin(\theta)\, \omega^2 - g \cos(\theta) \sin(\theta) \\
        a &= \frac{F - \cart{b}\, v + \pend{m}\, \pend{a}}{\cart{m} + \pend{m}} \\
        \alpha &= \frac{g \sin(\theta) - a \cos(\theta)}{\pend{l}} \\
        \begin{pmatrix} \dot \theta \\ \dot \omega \\ \dot x \\ \dot v \end{pmatrix}
        &= \begin{pmatrix} \omega \\ \alpha \\ v \\ a \end{pmatrix}
    \end{aligned}

A quadratic cost function is used, applied to the sine of half the pendulum angle,
to make it periodic with period :math:`2\pi`:

.. math::

    \ell(x, v, \theta, \omega, F) = q_\theta \sin^2(\theta/2) + q_\omega\, \omega^2 + q_x\, x^2 + q_v\, v^2 + r_F\, F^2

The force :math:`F` exerted on the cart is limited to :math:`\pm 2 \mathrm{N}`.

.. literalinclude:: ../../../../../examples/Python/mpc/inverted-pendulum-mpc.py
    :language: python
    :linenos:
