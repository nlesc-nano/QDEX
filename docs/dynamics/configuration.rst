Configuration
=============

Part of :doc:`/dynamics/index`.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.namd.precompute``
* Callable: ``qdex.namd.precompute.precompute_namd_data``
* CLI: ``--namd``
* YAML: ``namd.engine, namd.dt_fs``

.. code-block:: python

   precompute_namd_data(config)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part6_namd/index.rst:1361-1363``

11. CLI Flags & YAML Configuration Reference
--------------------------------------------


.. rubric:: From ``docs/part6_namd/index.rst:1364-1416``

Command-Line Arguments
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - CLI Flag
     - Default
     - Description
   * - ``--namd-precompute``
     - ``False``
     - Execute Stage 1 trajectory precomputation (cross-overlaps, NACs, phase tracking, caching).
   * - ``--namd-run``
     - ``False``
     - Execute Stage 2 NAMD carrier cooling simulation from precomputed data.
   * - ``--namd-decoherence [dir]``
     - ``None``
     - Compute and cache state-pair pure-dephasing matrices (:math:`\tau_{ij}`) into ``decoherence_times.npz``.
   * - ``--namd-compact [dir]``
     - ``None``
     - Compress precomputed directory, eliminating redundant duplicate arrays.
   * - ``--namd-soc``
     - ``False``
     - Enable relativistic Spin-Orbit Coupling across NAMD precomputation.
   * - ``--namd-ta``
     - ``False``
     - Compute ultrafast pump-probe transient absorption (TA) spectra from NAMD dynamics.
   * - ``--namd-ta-sigma <float>``
     - ``0.03``
     - Gaussian line broadening in eV for transient absorption probe spectra.
   * - ``--namd-ta-plot``
     - ``False``
     - Generate 2D false-color TA map and 1S bleach rise kinetics plot.
   * - ``--namd-ecsh-auger``
     - ``False``
     - Enable Energy-Conserving Surface Hopping (ECSH) for two-body Auger processes during NAMD.
   * - ``--namd-biexciton``
     - ``False``
     - Initialize NAMD from a biexciton state (XX) to simulate Auger annihilation dynamics.
   * - ``--namd-initial-conditions <mode>``
     - ``"single"``
     - Set initial condition sampling: ``"single"`` (start from :math:`t_0 = 0`) or ``"multiple"`` (automated ensemble sampling across uncorrelated trajectory origins).
   * - ``--namd-multi-init``
     - ``False``
     - Convenience shortcut for ``--namd-initial-conditions multiple``.
   * - ``--namd-origins <int>``
     - ``None`` (auto)
     - Explicit number of ensemble origins (when unset, calibrated automatically from :math:`\Delta t_0 \ge 2\tau_{\mathrm{corr}}`).
   * - ``--namd-window-fs <float>``
     - ``None`` (auto)
     - Simulation window duration in fs per origin (when unset, calibrated automatically from pilot cooling :math:`3\tau_{\mathrm{cool}}`).


.. rubric:: From ``docs/part6_namd/index.rst:1417-1448``

YAML Configuration Example
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   namd:
     trajectory:
       dir: "./trajectory"               # Directory containing frame_* subdirectories
       dt_nuc_fs: 2.0                    # Nuclear MD time step in femtoseconds
       start_frame: 1
       end_frame: 500

     dynamics:
       method: "dish"                    # "master_equation" (PME), "dish" (DISH), or "cpa_fssh" (FSSH-EDC)
       initial_conditions: "multiple"    # "single" (t0 = 0) or "multiple" (auto-calibrated ensemble)
       temperature_k: 300.0              # Lattice temperature for detailed balance
       tau_dec_fs: "cumulant"            # "cumulant" (ab initio), "edc", or fixed float in fs
       decoherence: "edc"                # Decoherence scheme for FSSH (continuous EDC)
       n_trajectories: 1000              # Trajectory count (split evenly across origins in multi-mode)
       detailed_balance: true            # Enforce Boltzmann detailed balance factor

     integration:
       integrator: "strang"              # Unitary Strang operator splitting
       n_substeps: 2                     # Electronic sub-steps per nuclear interval

     transient_absorption:
       run: true                         # Enable pump-probe transient absorption calculation
       sigma: 0.03                       # Probe spectral broadening in eV
       plot: true                        # Generate 2D TA map and kinetics figure
       plot_file: "transient_absorption_map.png"
       csv_file: "ta_bleach_kinetics.csv"
