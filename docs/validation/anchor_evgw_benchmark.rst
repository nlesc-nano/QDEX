evGW anchor benchmark (Cd\ :sub:`16`\ Se\ :sub:`13`\ Cl\ :sub:`6`)
=====================================================================

Part of :doc:`/validation/index`.

The CdSe "monomer" of ``MATERIAL_DB`` is the 35-atom cluster Cd\ :sub:`16`\ Se\ :sub:`13`\ Cl\ :sub:`6`
(:math:`D\approx1.2` nm). Its PBE and evGW frontier energies, computed with CP2K, are the finite-size
anchor of the ``gw`` model. For this cluster the reference is therefore a real evGW calculation, and
every QP model can be compared with it directly. Inputs are in ``tests/CdSe/1.2nm``; the MO file is
stored as ``MOs_cleaned_12ang.txt.gz``.

Reference (``MATERIAL_DB["CDSE"]`` monomer entries, vacuum, eV):

.. list-table::
   :header-rows: 1

   * -
     - HOMO
     - LUMO
     - gap
   * - PBE
     - −6.4196
     - −3.7852
     - 2.634
   * - evGW
     - −7.8177
     - −1.7884
     - 6.029
   * - shift
     - −1.398
     - +1.997
     - +3.395 (41 % HOMO / 59 % LUMO)

The DFT gap of the supplied MO file is 2.639 eV, 5 meV from the stored PBE gap. The bulk GW–PBE
opening is 1.270 eV, so evGW requires a finite-size correction of 2.125 eV on top of it.

Results with the current QP layer
---------------------------------

This section reflects the September 2026 revision:

* sphere polarization with an exact anchor for ``gw``;
* Penn-gap :math:`\epsilon_{\mathrm{in}}(R)` for the Resta family;
* Z derived from the model's own plasmon pole (the default);
* HOMO/LUMO split from the per-edge anchor curves.

Vacuum, spin-free, 25 × 25. Δ = model − evGW in eV. S\ :sub:`1` uses each model's consistent
kernel: ``qp`` for the Delta-W models and ``resta-sphere`` for ``gw``.

.. list-table::
   :header-rows: 1

   * - qp_gap
     - Z
     - QP gap
     - Δ gap
     - Δ HOMO
     - Δ LUMO
     - S\ :sub:`1` (vac / toluene)
   * - ``gw``
     - —
     - 6.034
     - +0.005
     - −0.002
     - +0.002
     - 3.177 / 3.132
   * - ``sgw-resta``
     - derived 0.95
     - 5.664
     - −0.37
     - +0.15
     - −0.21
     - 3.325 / 3.400
   * - ``sgw-resta``
     - 1.0
     - 5.764
     - −0.26
     - +0.11
     - −0.16
     - 3.314 / 3.398
   * - ``sgw-resta-pure``
     - derived 0.95
     - 5.444
     - −0.58
     - +0.24
     - −0.34
     - 3.373
   * - ``sgw-dim``
     - derived 0.95
     - 5.576
     - −0.45
     - +0.19
     - −0.27
     - 3.369 / 3.439
   * - ``evgw-resta``
     - derived 0.94
     - 5.880
     - −0.15
     - +0.06
     - −0.09
     - 3.286
   * - ``evgw-dim``
     - derived 0.95
     - 5.565
     - −0.46
     - +0.19
     - −0.27
     - 3.366
   * - ``qsgw-resta``
     - derived 0.94
     - 6.251
     - +0.22
     - −0.09
     - +0.13
     - 3.677
   * - ``qsgw-dim``
     - derived 0.95
     - 5.868
     - −0.16
     - +0.07
     - −0.09
     - 3.686
   * - ``sgw`` (sBSE)
     - derived 0.98
     - 6.684
     - +0.66
     - −0.27
     - +0.38
     - 4.524 / 4.488

What changed:

* **``gw`` reproduces the anchor exactly.** :math:`R_0` was recomputed with the code's radius
  definition (5.3133 Å). The remaining +5 meV is the PBE-gap difference of the MO file. Both
  frontier levels match evGW to 2 meV.
* **The Resta interior term is now modest.** With the Penn gap (:math:`E_P` = 6.2 eV),
  :math:`\epsilon_{\mathrm{in}}` is 3.97 instead of 1.48. Resta and DIM now agree to 0.1 eV:
  −0.37 and −0.45 eV. Both underestimate the evGW gap by 0.3–0.5 eV.
* **The closest models are ``evgw-resta`` (−0.15 eV) and ``qsgw-dim`` (−0.16 eV).**
* **Derived Z.** Z ≈ 0.95 and changes the gap by 0.1 eV relative to Z = 1. It is no longer the
  dominant knob.
* **IP/EA.** With the anchor-calibrated split, the HOMO and LUMO errors are each about half of the
  gap error, with the right signs. The old 50/50 split put up to 0.8 eV error on one edge.
* **S₁ consistency.** All consistent Delta-W routes give S₁ = 3.29–3.37 eV and ``gw`` gives
  3.18 eV. The qsGW models and sBSE are higher. With the bulk-only Resta kernel, ``gw`` would give
  5.58 eV, because the kernel lacks the surface polarization that the QP gap contains.

The tables below are the results before this revision and are kept for reference.

Results before the revision
---------------------------

Vacuum, spin-free. Each Delta-W model is run with its shared-:math:`W` BSE kernel (``kernel: qp``) and
the gap-only models with ``resta``. The QP gap and edges do not depend on the kernel. Δ is the deviation
from evGW in eV.

.. list-table::
   :header-rows: 1

   * - ``qp_gap``
     - Z
     - QP gap
     - Δ gap
     - Δ HOMO
     - Δ LUMO
   * - ``gw`` (anchor)
     - —
     - 6.009
     - −0.02
     - +0.01
     - −0.01
   * - ``brus``
     - —
     - 6.220
     - +0.19
     - −0.08
     - +0.11
   * - ``sgw-resta``
     - 0.8
     - 6.215
     - +0.19
     - −0.45
     - −0.26
   * - ``sgw-resta``
     - 1.0
     - 6.792
     - +0.76
     - −0.74
     - +0.02
   * - ``sgw-resta-pure``
     - 0.8
     - 5.201
     - −0.83
     - +0.07
     - −0.76
   * - ``sgw-resta-pure``
     - 1.0
     - 5.524
     - −0.51
     - −0.10
     - −0.61
   * - ``sgw-dim``
     - 0.8
     - 5.316
     - −0.71
     - −0.14
     - −0.86
   * - ``sgw-dim``
     - 1.0
     - 5.667
     - −0.36
     - −0.34
     - −0.71
   * - ``evgw-resta``
     - derived (0.93)
     - 6.338
     - +0.31
     - −0.51
     - −0.20
   * - ``evgw-dim``
     - derived (0.92)
     - 5.514
     - −0.52
     - −0.21
     - −0.72
   * - ``qsgw-resta``
     - derived (0.93)
     - 6.726
     - +0.70
     - −0.83
     - −0.13
   * - ``qsgw-dim``
     - derived (0.92)
     - 5.810
     - −0.22
     - −0.23
     - −0.44
   * - ``sgw`` (sBSE)
     - 0.8
     - 6.182
     - +0.15
     - −0.43
     - −0.27
   * - ``sgw`` (sBSE)
     - 1.0
     - 6.750
     - +0.72
     - −0.72
     - +0.00

Interpretation
--------------

* **``gw`` is not tested here.** It is interpolated through this point. Before the revision, the
  remaining −20 meV came from the radius: QDEX gives :math:`R_{\mathrm{eff}}` = 5.313 Å for this
  geometry, while ``MATERIAL_DB`` stored :math:`R_0` = 5.258 Å. It is now 5.3133 Å.
* **The interior term separates the models.** The solvent/surface term is nearly the same for Resta and
  DIM (1.29 eV with Z = 0.8, 1.61 eV with Z = 1). The interior contrast differs by a factor of nine:
  1.01 eV for Resta with its size-scaled :math:`\epsilon_{\mathrm{eff}}`, and 0.11 eV for DIM. Only the
  Resta scaling brings the total near the required 2.125 eV.
* **Z is as large an effect as the choice of model.** ``sgw-resta`` goes from +0.19 eV (Z = 0.8) to
  +0.76 eV (Z = 1). Its agreement at Z = 0.8 relies on the factor 0.8 compensating an interior term that
  is too large. With Z = 1 the closest Delta-W model is ``sgw-dim`` at −0.36 eV; ``sgw-resta`` is at
  +0.76 eV.
* **Iterating does not improve the gap systematically.** ``evgw-*`` adds 0.1–0.2 eV. ``qsgw-resta``
  overshoots by 0.7 eV.
* **The HOMO/LUMO split is wrong in all Delta-W models.** evGW places 41 % of the correction on the HOMO;
  the charging term :math:`\tfrac12 q^T\Delta W q` places about 50 %, because it is nearly symmetric for
  electrons and holes. The asymmetry of the true self-energy is non-classical exchange–correlation,
  which these models do not contain. Absolute IP/EA from the Delta-W models are therefore
  unreliable by several tenths of an eV even when the gap is right.
* **Brus lands near evGW by coincidence.** It is only the effective-mass kinetic term added to the
  experimental bulk gap.

Scope
-----

This is one cluster at the size where the ``gw`` anchor is defined. It tests the finite-size term of
each model, not its transferability to other sizes. evGW@PBE is itself an approximation to the true
quasiparticle gap. Rerun it with:

.. code-block:: bash

   python benchmarks/compare_models.py --system tests/CdSe/1.2nm --profile full --only A
