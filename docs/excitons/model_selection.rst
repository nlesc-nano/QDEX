Model selection
===============

Part of :doc:`/excitons/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.solver``
* Callable: ``qdex.solver.ExcitonSolver.solve``
* CLI: ``--excitation-mode, --include-direct-eh, --include-exchange``
* YAML: ``physics.excitation_mode, physics.include_direct_eh, physics.include_exchange``

.. code-block:: python

   solve(self, nroots=10, full_diag=False, tol=1e-05, excitation_mode='bse')


8. Comprehensive Exciton Calculation Matrix
-------------------------------------------

.. list-table::
   :widths: 18 15 15 14 38
   :header-rows: 1

   * - ``excitation_mode``
     - ``2e-integrals``
     - ``kernel``
     - Complexity
     - Recommended Purpose
   * - ``bse``
     - ``mnok``
     - ``resta``
     - :math:`O(N_{\mathrm{atoms}}^2)`
     - **Production Standard**: Absorption spectra for large colloidal quantum dots (up to 10,000 atoms).
   * - ``bse``
     - ``mnok``
     - ``dim``
     - :math:`O(N_{\mathrm{atoms}}^2)`
     - Core/shell or anisotropic nanocrystals with polarizable dielectric boundaries.
   * - ``diagonal_bse``
     - ``mnok``
     - ``resta``
     - :math:`O(N_{\mathrm{pairs}})`
     - **NAMD Production**: Ultrafast excited-state dynamics and non-adiabatic trajectories.
   * - ``bse``
     - ``xs``
     - ``rpa``
     - :math:`O(N_{\mathrm{ao}}^4)`
     - **Ab Initio Benchmark**: Analytical AO density-pair Coulomb integrals and model RPA screening; benchmark accuracy requires validation.
   * - ``bse``
     - ``xs``
     - ``dim``
     - :math:`O(N_{\mathrm{ao}}^4)`
     - Exact Gaussian integrals with atomistic polarizable dipole screening.
   * - ``independent_qp``
     - None
     - None
     - :math:`O(N_{\mathrm{pairs}})`
     - Single-particle joint density of states with quasiparticle gap correction.
   * - ``independent_dft``
     - None
     - None
     - :math:`O(1)`
     - Uncorrected baseline Kohn-Sham single-particle transitions.

