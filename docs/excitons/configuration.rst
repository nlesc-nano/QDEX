Configuration
=============

Part of :doc:`/excitons/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.solver``
* Callable: ``qdex.solver.ExcitonSolver.solve``
* CLI: ``--excitation-mode, --include-direct-eh, --include-exchange``
* YAML: ``physics.excitation_mode, physics.include_direct_eh, physics.include_exchange``

.. code-block:: python

   solve(self, nroots=10, full_diag=False, tol=1e-05, excitation_mode='bse')


10. CLI Flags & YAML Configuration Reference
--------------------------------------------


Command-Line Arguments
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - CLI Flag
     - Default
     - Description
   * - ``--excitation-mode <choice>``
     - ``bse``
     - Excitation framework: ``bse`` (full sTDA), ``diagonal_bse``, ``independent_qp``, or ``independent_dft``.
   * - ``--2e-integrals <choice>``
     - ``mnok``
     - Two-electron integral representation: ``mnok`` (semi-empirical atom-centered) or ``xs`` (analytical Gaussian AO density-pair integrals).
   * - ``--kernel <choice>``
     - ``bse``
     - Dielectric screening kernel: ``resta``, ``dim``, ``rpa``, ``sbse``, ``xs-resta``, ``xs-dim``, ``xs-rpa``, or ``bse``.
   * - ``--nhomos <int>``
     - ``25``
     - Number of occupied frontier orbitals to include in the BSE active space.
   * - ``--nlumos <int>``
     - ``25``
     - Number of virtual frontier orbitals to include in the BSE active space.
   * - ``--e_thresh <float>``
     - ``None``
     - Energy threshold (in eV) to automatically select active pairs with :math:`\varepsilon_a - \varepsilon_i \le E_{\mathrm{thresh}}`.
   * - ``--f_thresh <float>``
     - ``0.0``
     - Minimum oscillator strength threshold to print and log excited states.
   * - ``--nroots <int>``
     - ``10``
     - Number of lowest exciton roots to compute via Davidson diagonalization.
   * - ``--full-diag``
     - ``False``
     - Force full dense LAPACK diagonalization instead of the iterative Davidson solver.
   * - ``--tol <float>``
     - ``1e-5``
     - Convergence tolerance for the Davidson solver residual norm.
   * - ``--triplet``
     - ``False``
     - Perform triplet excited-state BSE calculation (omitting :math:`2K^x`).
   * - ``--charge_type <choice>``
     - ``mulliken``
     - Transition charge partitioning: ``mulliken`` or ``lowdin``.
