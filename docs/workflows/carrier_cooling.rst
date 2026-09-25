Carrier cooling
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


Tutorial 6: Carrier Cooling Dynamics & Photoluminescence (NAMD)
---------------------------------------------------------------

In this final tutorial, we simulate non-adiabatic hot electron and hole cooling across an *ab initio* molecular dynamics trajectory of 500 frames, compute cooling rates, and calculate the Photoluminescence Quantum Yield.


Directory Layout
~~~~~~~~~~~~~~~~

.. code-block:: text

   trajectory_dir/
   ├── frame_0001/ (contains CsPbBr3.xyz, MOs.mbse)
   ├── frame_0002/
   ...
   └── frame_0500/


Configuration File (``tutorial6_namd.yaml``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   system:
     basis_txt: "BASIS_MOLOPT"
     basis_name: "DZVP-MOLOPT-PBE-GTH"

   physics:
     excitation_mode: "diagonal_bse"
     kernel: "resta"
     qp_gap: "gw"
     material: "CSPBBR3"
     nhomos: 25
     nlumos: 25

   namd:
     trajectory_dir: "trajectory_dir"
     dt_fs: 1.0
     temperature_k: 300.0
     engine: "master_equation"
     tau_dec_fs: "cumulant"
     
     recombination:
       include_ground_state: true
       radiative: true
       tau_nr_ns: 20.0
     
     storage:
       precompute_dir: "namd_precomputed"
       output_dir: "namd_results"


Step 1: Trajectory Precomputation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   qdex --config tutorial6_namd.yaml --namd-precompute

This step computes cross-frame overlaps :math:`S_{IJ}(t, t+\Delta t)`, applies gauge phase corrections, and tracks trivial state crossings via the Hungarian algorithm.


Step 2: Carrier Cooling Simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   qdex --config tutorial6_namd.yaml --namd-run


Results & Visualizations
~~~~~~~~~~~~~~~~~~~~~~~~

In ``namd_results/``:

* ``carrier_cooling_populations.csv``: Time-dependent populations of electron and hole levels.
* ``cooling_curves.png``: Transient excess energy :math:`\Delta E(t)` showing hot hole cooling (:math:`\tau_{\mathrm{hole}} \approx 120\text{ fs}`) and hot electron cooling (:math:`\tau_{\mathrm{elec}} \approx 350\text{ fs}`).
* ``photoluminescence_yield.txt``: Reports the calculated Einstein radiative rate :math:`k_{\mathrm{rad}} \approx 4.5 \times 10^7\text{ s}^{-1}` and the total **PLQY** :math:`\approx 47.4\%`.
