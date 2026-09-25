QP × BSE sweep
==============

Part of :doc:`/validation/index`.

``benchmarks/qp_bse_sweep.py`` runs every consistent combination of two-electron representation, QP
model and BSE kernel for one quantum dot. Copy it into a folder that holds ``config.yaml`` and its
input files:

.. code-block:: bash

   cp benchmarks/qp_bse_sweep.py tests/CdSe/2.0nm/
   cd tests/CdSe/2.0nm
   python qp_bse_sweep.py --list                  # 69 cases
   python qp_bse_sweep.py --dry-run               # write sweep/<case>/config.yaml only
   python qp_bse_sweep.py --jobs 2 --nthreads 4   # run (resumable)
   python qp_bse_sweep.py --collect               # sweep/summary.md, summary.csv

Each case folder contains the complete ``config.yaml`` of that combination, so it can be rerun with
``qdex --config config.yaml``.

One W for QP and BSE
--------------------

.. list-table::
   :header-rows: 1

   * - QP model
     - W of the QP correction
     - BSE kernel
   * - ``gw``, ``qp_polarization: sphere``
     - bulk Resta + dielectric-sphere reaction field
     - ``resta-sphere`` (the same W)
   * - ``gw``, ``qp_polarization: legacy``
     - none (κ/(R+ℓ) curve)
     - ``resta`` (bulk; not consistent, legacy reference)
   * - ``brus``
     - none (kinetic confinement)
     - ``resta`` (bulk)
   * - ``sgw-resta``, ``evgw-resta``, ``qsgw-resta``
     - Resta W with the Penn ε\ :sub:`in`\ (R) + solvent term
     - ``qp`` (the same W)
   * - ``sgw-dim``, ``evgw-dim``, ``qsgw-dim``
     - DIM/Thole W + solvent term
     - ``qp`` (the same W)

**Two-electron representation.** ``two_electron_integrals`` selects the representation of W for both
the QP correction and the kernel.

* ``mnok``: atom-pair Ohno–Klopman interaction γ\ :sub:`AB`; atomic populations.
* ``xs``: exact AO density-pair integrals :math:`(\mu\mu|\nu\nu)`. The model's screening ratio
  :math:`W_{AB}/\gamma_{AB}` multiplies each AO block, the additive classical terms (solvent, sphere
  image) are expanded to AO blocks, and the QP corrections use AO populations. qsGW builds its
  Δ-COHSEX matrix from the same AO ΔW.

**QP levels.** ``qp_levels: orbital`` (default) corrects every orbital p of a window of
max(100, nhomos, nlumos) + 10 occupied and virtual orbitals:

.. math::

   \varepsilon_i^{\mathrm{QP}} = \varepsilon_i - f_b\Delta_{\mathrm{bulk}} - Z_i\sigma_i,\qquad
   \varepsilon_a^{\mathrm{QP}} = \varepsilon_a + (1-f_b)\Delta_{\mathrm{bulk}} + Z_a\sigma_a,\qquad
   \sigma_p = \tfrac12\,\mathbf q_p^{\mathsf T}\Delta W\,\mathbf q_p .

This is the classical form (``qp_selfenergy: classical``). The default for the Delta-W models is the
one-shot static ΔCOHSEX diagonal, the same form qsGW uses, applied to every orbital without an orbital
update (``qp_selfenergy: cohsex``):

.. math::

   \Sigma_n = \tfrac12\sum_\mu c_{\mu n}^2\,\Delta W_{\mu\mu}
   \;-\;\tfrac12\sum_{\mu\nu} c_{\mu n}c_{\nu n}\,P_{\mu\nu}\,\Delta W_{\mu\nu},
   \qquad c = S^{1/2}C,\; P = 2c_{\mathrm{occ}}c_{\mathrm{occ}}^{\mathsf T}.

The classical term is its limit when the occupied states act as a complete set. The difference, the
non-classical screened exchange, does not cancel against the binding. At 1.2 nm it lifts the
``sgw-resta`` QP gap from 5.66 to 5.96 eV (evGW: 6.03 eV) and S₁ by about 0.4 eV, to the qsGW value.
The calibrated anchor residual (:doc:`/quasiparticles/anchor`) is then added.

* **qsGW:** already returns orbital energies and orbitals, which are used as they are.
* **gw sphere model:** the frontier shifts are pinned to the per-edge anchor curves. Every other
  orbital adds its image term relative to the frontier orbital: σ\ :sub:`p` − σ\ :sub:`H` for
  occupied orbitals and σ\ :sub:`p` − σ\ :sub:`L` for virtual ones.
* **Models without W** (``gw`` legacy, ``brus``): a rigid scissor on the virtual orbitals.

``qp_levels: rigid`` restores the scissor for every model.

Groups
------

.. list-table::
   :header-rows: 1

   * - group
     - question
     - cases
   * - A
     - core matrix: {mnok, xs} × 9 QP models × {vacuum, solvent}
     - 36
   * - B
     - Z = 1 and Z = 0.8 against the derived Z
     - 7
   * - C
     - rigid scissor against orbital-resolved levels
     - 6
   * - D
     - two-anchor sensitivity: residual power p, bulk-only kernel
     - 4
   * - E
     - SOC in the solvent, main models, mnok and xs
     - 8
   * - F
     - active space 50 × 50 and 100 × 100 (Davidson)
     - 4
   * - G
     - Löwdin transition charges
     - 2
   * - H
     - triplet (singlet–triplet splitting)
     - 2

For Cd₁₆Se₁₃Cl₆ the model's own QP gap matches evGW directly, since it is the anchor
(:doc:`anchor_evgw_benchmark`). For larger dots, compare the bright SOC state in the solvent with
experiment (``compare_models.py --exp-ref``).
