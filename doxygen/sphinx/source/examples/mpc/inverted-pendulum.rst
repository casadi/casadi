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

.. literalinclude:: ../../../../../examples/Python/mpc/inverted-pendulum-mpc.py
    :language: python
    :linenos:
