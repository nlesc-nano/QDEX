Model selection
===============

Part of :doc:`/excitons/index`.

The QP model fixes the BSE kernel (:doc:`screened_kernel`). What remains to choose is the framework,
the representation and the active space.

.. list-table::
   :header-rows: 1
   :widths: 22 16 14 48

   * - ``qp_gap``
     - ``excitation_mode``
     - ``2e-integrals``
     - Purpose
   * - ``sgw-resta`` / ``sgw-dim``
     - ``bse``
     - ``mnok``
     - **Production:** absorption spectra and size series of colloidal dots, in the solvent. DIM for
       shape, ligand or shell effects.
   * - ``gw``
     - ``bse``
     - ``mnok``
     - Cheapest consistent estimate (kernel ``resta-sphere``).
   * - ``qsgw-*``
     - ``bse``
     - ``mnok``
     - Orbital relaxation; small and medium dots.
   * - any of the above
     - ``bse``
     - ``xs``
     - Small clusters and molecules, where the short range matters.
   * - any of the above
     - ``diagonal_bse``
     - ``mnok``
     - Non-adiabatic dynamics; check against ``bse``.
   * - ``brus``
     - ``bse``
     - ``mnok``
     - ΔW = 0 limit; use ``kernel: resta``.
   * - any
     - ``independent_qp`` / ``independent_dft``
     - —
     - Joint density of states; reference only.

**Active space.** S₁ converges slowly with the number of transitions at 2 nm (2.499 → 2.443 eV from
25 × 25 to 100 × 100). Keep the active space fixed within a comparison.
