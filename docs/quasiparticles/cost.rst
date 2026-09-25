Cost
====

Part of :doc:`/quasiparticles/index`.

Full G₀W₀ scales as N⁴ (localized basis) and needs the response of thousands of empty states. That
is prohibitive for dots of 10³–10⁴ atoms beyond single calculations (:doc:`theory`, section 3). The
QDEX models cost, with n_ao basis functions and N_at atoms:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Step
     - Cost
   * - Populations of all orbitals (Mulliken)
     - one n_ao³ product (S C)
   * - ΔW (Resta / DIM)
     - N_at²; DIM adds one 3N_at × 3N_at inversion
   * - One-shot ΔCOHSEX, all orbitals
     - one cached eigendecomposition of S + three n_ao³ products
   * - ``evgw-*``
     - the ΔW step per iteration (10–20 iterations), then ΔCOHSEX once
   * - ``qsgw-*``
     - one n_ao³ diagonalization per iteration
   * - xs representation
     - exact (μμ|νν) integrals and several n_ao² matrices (about 1.4 GB each at 13k basis functions)

For 10⁴ basis functions the QP step is minutes with mnok. The orthonormality check of the MO file
(one n_ao³ product) can be skipped with ``skip_orthonormality_check``.
