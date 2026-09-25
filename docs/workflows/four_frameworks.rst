Four frameworks
===============

Part of :doc:`/workflows/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.cli``
* Callable: ``qdex.cli.main``
* CLI: ``--config``
* YAML: ``system.*, physics.*, namd.*``

.. code-block:: python

   main()


Tutorial 4: Comparing the Four Excited-State Frameworks
-------------------------------------------------------

In this tutorial, we solve for the optical excitations using all four frameworks on the same system to observe the physical differences.


Compare via CLI
~~~~~~~~~~~~~~~

1. **Bare DFT Transitions**:
   
   .. code-block:: bash

      qdex --config tutorial3_gw.yaml --excitation-mode independent_dft --nroots 10

2. **Quasiparticle Transitions (Independent QP)**:

   .. code-block:: bash

      qdex --config tutorial3_gw.yaml --excitation-mode independent_qp --nroots 10

3. **Diagonal BSE**:

   .. code-block:: bash

      qdex --config tutorial3_gw.yaml --excitation-mode diagonal_bse --kernel resta --nroots 10

4. **Full BSE (sTDA with Davidson Solver)**:

   .. code-block:: bash

      qdex --config tutorial3_gw.yaml --excitation-mode bse --kernel resta --nroots 10


Comparison Summary
~~~~~~~~~~~~~~~~~~

* **Gap Ordering**: :math:`E_{\mathrm{opt}}(\text{DFT}) \ll E_{\mathrm{opt}}(\text{sTDA}) \approx E_{\mathrm{opt}}(\text{Diag-BSE}) < E_{\mathrm{opt}}(\text{QP})`.
* **Binding Energy**: Subtracting :math:`E_1(\text{sTDA})` from :math:`E_1(\text{QP})` gives a model-dependent estimate of binding relative to the chosen QP reference, not a universal 280 meV value.
* **Superradiance**: Notice that in full BSE, oscillator strength concentrates strongly into Root 1 due to coherent transition dipole redistribution, whereas independent transitions remain diffuse.

