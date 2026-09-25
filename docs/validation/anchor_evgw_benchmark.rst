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

Results
-------

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

* **``gw`` is not tested here.** It is interpolated through this point. The remaining −20 meV comes from
  the radius: QDEX gives :math:`R_{\mathrm{eff}}` = 5.313 Å for this geometry, while ``MATERIAL_DB``
  stores :math:`R_0` = 5.258 Å.
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
