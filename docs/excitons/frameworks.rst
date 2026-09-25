Frameworks
==========

Part of :doc:`/excitons/index`.

.. important::

   ``diagonal_bse`` omits all off-diagonal transition mixing. It can miss a substantial part of binding; it cannot generally be assumed to reproduce the bulk Wannier exciton. Exchange and direct attraction have independent switches.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.exciton_hamiltonian``
* Callable: ``qdex.exciton_hamiltonian.independent_transition_energies``
* CLI: ``--excitation-mode, --include-direct-eh, --include-exchange``
* YAML: ``physics.excitation_mode, physics.include_direct_eh, physics.include_exchange``

.. code-block:: python

   independent_transition_energies(self, mode)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part4_excited_states/index.rst:314-338``

5. The Four Excitation Frameworks (``excitation_mode``)
-------------------------------------------------------

``QDEX`` provides four progressive levels of physical theory via ``--excitation-mode``:

.. list-table::
   :widths: 22 28 50
   :header-rows: 1

   * - Framework Mode
     - Energy Expression
     - Physical Characteristics
   * - **(A) `independent_dft`**
     - :math:`\Omega_{ia} = \varepsilon_a^{\mathrm{DFT}} - \varepsilon_i^{\mathrm{DFT}}`
     - Bare Kohn-Sham transitions. Completely ignores quasiparticle self-energy corrections and electron-hole Coulomb interactions.
   * - **(B) `independent_qp`**
     - :math:`\Omega_{ia} = \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}}`
     - Non-interacting quasiparticles. Applies the scaled GW scissor shift :math:`\Delta_{\mathrm{GW}}`, opening the gap to experimental values, but neglects electron-hole binding (:math:`E_b = 0`).
   * - **(C) `diagonal_bse`**
     - :math:`\Omega_{ia} = (\varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}}) + 2 K_{ia,ia}^x - K_{ia,ia}^d`
     - Diagonal BSE. Accounts for both quasiparticle self-energy and diagonal electron-hole Coulomb binding, but omits off-diagonal configuration mixing. Extremely fast.
   * - **(D) `bse` (sTDA)**
     - Full diagonalization of :math:`A_{ia, jb}`
     - Fully coupled configuration interaction. Solves the complete resonant matrix, capturing spatial exciton delocalization, state mixing, and oscillator strength redistribution.


.. rubric:: From ``docs/part4_excited_states/index.rst:339-352``

Why Diagonal BSE Works in Nanocrystals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In standard BSE, diagonalizing the full :math:`A_{ia, jb}` matrix requires :math:`O(N_{\mathrm{pairs}}^3)` operations. In a quantum dot with 100 occupied and 100 virtual orbitals, :math:`N_{\mathrm{pairs}} = 10,000`, requiring gigabytes of memory and long compute times.

The **Diagonal BSE** framework omits off-diagonal configuration interaction (:math:`ia \neq jb`), evaluating:

.. math::

   \Omega_{ia} = \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}} + 2 K_{ia, ia}^x - K_{ia, ia}^d.

* **Dominance of Diagonal Coulomb Attraction**: In quantum dots, the diagonal term :math:`K_{ia, ia}^d` represents the direct electrostatic attraction between the electron distribution :math:`|\phi_a|^2` and hole distribution :math:`|\phi_i|^2`, whose share of the coupled exciton binding energy is system-dependent; the supplied CdSe example gives 58.3% for diagonal BSE.
* **Essential for NAMD**: In non-adiabatic molecular dynamics simulations where excited states must be evaluated at every time step (e.g. 5,000 steps), full BSE diagonalization is computationally prohibitive. Diagonal BSE provides a less expensive approximate surface; its error must be checked against coupled BSE for the chosen active space.


.. rubric:: From ``docs/part4_excited_states/index.rst:353-368``

Implementation in QDEX (``--excitation-mode``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The excitation frameworks are implemented across :mod:`qdex.cli`, :mod:`qdex.exciton_hamiltonian`, and :mod:`qdex.solver`:

1. **Bare DFT Transitions (``independent_dft``)**:
   Returns raw differences :math:`\Omega_{ia} = \varepsilon_a^{\mathrm{DFT}} - \varepsilon_i^{\mathrm{DFT}}` directly from Kohn-Sham eigenvalues without building two-body Coulomb matrices.
2. **Quasiparticle Transitions (``independent_qp``)**:
   Applies the QP scissor shifts :math:`\varepsilon^{\mathrm{QP}}` from Part 3: :math:`\Omega_{ia} = \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}}`. Evaluates single-particle transition dipoles :math:`\boldsymbol{\mu}_{ia}` and oscillator strengths without electron-hole binding.
3. **Diagonal BSE (``diagonal_bse``)**:
   Implemented via :func:`qdex.solver.ExcitonSolver.solve`. Builds only the diagonal elements :math:`A_{ia, ia} = (\varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}}) + 2 K_{ia,ia}^x - K_{ia,ia}^d` in :math:`O(N_{\mathrm{pairs}})` time. Completely bypasses matrix diagonalization, making it the ideal engine for multi-thousand step non-adiabatic dynamics.
4. **Coupled Bethe-Salpeter Equation (``bse``)**:
   Implemented via :func:`qdex.davidson.davidson`. Constructs the active space transition basis, applies energy truncation thresholds (``--e_thresh``), and solves for the lowest :math:`N_{\mathrm{roots}}` exciton eigenvectors using a block-Davidson iterative subspace algorithm.

---
