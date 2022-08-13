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

.. literalinclude:: ../../../../../examples/Python/mpc/inverted-pendulum-mpc.py
    :language: python
    :linenos:
