Hanging Chain
=============

In this example, a mode predictive controller (MPC) is used to stabilize a 
system of weights connected by springs. The rightmost weight is fixed in place 
at the origin, whereas the velocity of the leftmost weight can be controlled by
an actuator. The six weights in the middle move under the influence of gravity 
and the forces of the springs between them.

The goal of the controller is to stabilize the system (i.e. drive the velocity
of all weights to zero) with the rightmost weight at position :math:`(1, 0)`. 
Additionally, a non-convex cubic constraint on the weights' position is imposed,
shown in green on the figure below.

.. raw:: html
   :file: ../../sphinxstatic/hanging-chain.html

.. raw:: html

    <script>
    window.addEventListener('load',
        function() {
            document.querySelectorAll(".anim-buttons > button:nth-child(6)")
                .forEach(btn => { btn.click(); });
        }, false);
    </script>

.. literalinclude:: ../../../../examples/mpc/python/hanging-chain/hanging-chain-mpc.py
    :language: python
    :linenos:

.. literalinclude:: ../../../../examples/mpc/python/hanging-chain/hanging_chain_dynamics.py
    :language: python
    :linenos:
