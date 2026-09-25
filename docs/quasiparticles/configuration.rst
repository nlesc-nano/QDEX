Configuration
=============

Part of :doc:`/quasiparticles/index`. All keys go in the ``physics`` section of the YAML file; the CLI
flag is given in brackets.

.. list-table::
   :header-rows: 1
   :widths: 26 16 58

   * - Key [flag]
     - Default
     - Meaning
   * - ``qp_gap`` [``--qp_gap``]
     - ``brus``
     - QP model: ``gw``, ``sgw-resta``, ``sgw-dim``, ``evgw-resta``, ``evgw-dim``, ``qsgw-resta``,
       ``qsgw-dim``; also ``sgw``, ``brus``, ``pbe`` or a gap in eV (:doc:`models`).
   * - ``kernel`` [``--kernel``]
     - model default
     - ``qp`` for the Delta-W models and ``resta-sphere`` for ``gw``, set automatically. Other kernels
       only for models that define no W.
   * - ``two_electron_integrals`` [``--two-electron-integrals``]
     - ``mnok``
     - Representation of W for the QP correction and the kernel: ``mnok`` (atom pairs) or ``xs``
       (exact AO density pairs).
   * - ``eps_out`` [``--eps-out``]
     - ``2.0``
     - Optical dielectric constant of the environment (vacuum 1, toluene 2.24).
   * - ``qp_selfenergy`` [``--qp-selfenergy``]
     - ``cohsex``
     - Delta-W levels: one-shot ΔCOHSEX (``cohsex``) or the classical ½ qᵀΔWq (``classical``).
   * - ``qp_solvent_term`` [``--qp-solvent-term``]
     - ``sphere``
     - Environment part of ΔW in the Delta-W models: dielectric-sphere reaction field (``sphere``) or
       the earlier softened Born term (``born``).
   * - ``qp_levels`` [``--qp-levels``]
     - ``orbital``
     - Correct every orbital (``orbital``) or apply one scissor (``rigid``).
   * - ``qp_z`` [``--qp-z``]
     - ``derived``
     - Quasiparticle weight: plasmon pole of the model's ε (``derived``) or a fixed number.
   * - ``charge_type`` [``--charge_type``]
     - ``mulliken``
     - Transition charges of the BSE; the QP populations follow the same partition (xs uses Löwdin).
   * - ``qp_anchor_residual`` [``--qp-anchor-residual``]
     - ``on``
     - Add the calibrated per-edge anchor residual to the Delta-W levels.
   * - ``qp_residual_scaling`` [``--qp-residual-scaling``]
     - ``econf``
     - Size scaling of the anchor residual: E_conf(R)/E_conf(R₀) (``econf``) or (R₀/R)^p (``power``).
   * - ``qp_residual_power`` [``--qp-residual-power``]
     - ``2.0``
     - p of the ``power`` scaling.
   * - ``qp_polarization`` [``--qp-polarization``]
     - ``sphere``
     - Finite-size term of ``gw``: dielectric sphere (``sphere``) or the older κ/(R + ℓ) (``legacy``).
   * - ``qp_edge_split`` [``--qp-edge-split``]
     - ``anchor``
     - HOMO/LUMO split for absolute IP/EA: anchor curves (``anchor``) or the model's own (``model``).
   * - [``--qp-anchor-calibrate``]
     - off
     - Run on the anchor cluster in vacuum to store this model's residual in
       ``qdex/data/dw_anchor_residuals.json`` (or ``--qp-anchor-table``).
   * - ``skip_orthonormality_check`` [``--skip-orthonormality-check``]
     - off
     - Skip the full Cᵀ S C test for MO files known to be orthonormal (saves one n_ao³ product).

Example for a size series in toluene:

.. code-block:: yaml

   physics:
     qp_gap: sgw-resta          # kernel qp, ΔCOHSEX levels, anchor residual: defaults
     eps_out: 2.24
     two_electron_integrals: mnok
     nhomos: 25
     nlumos: 25
     skip_orthonormality_check: true
