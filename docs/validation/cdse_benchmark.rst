CdSe benchmark provenance
=========================


.. figure:: /_static/figures/cdse_energy_ladder.svg
   :width: 100%
   :alt: cdse energy ladder

   Energy decomposition for the CdSe test run (gw QP model, Resta direct kernel, 25×25 active space, spin-free).

``MATERIAL_DB["CDSE"]`` currently stores a bulk PBE gap of 0.64 eV, a bulk GW
gap of 1.91 eV, and an optical dielectric constant of 6.2. The bulk gap
opening used by the corresponding calibrated models is therefore 1.27 eV.
Those values form one named parameter set; numbers from a different CdSe
calculation must not be mixed into it without a separate source and convention.

The original audit request supplied alternative bulk PBE gaps of 1.14 and
1.25 eV and a GW gap of 2.76 eV. These imply openings of 1.62 and 1.51 eV,
respectively, and cannot both define the same bulk correction. Its proposed
55:45 frontier split is a benchmark hypothesis, not a CdSe material constant.
The split depends on orbitals, surface, screening model, and charge
partition. Likewise, a diagonal BSE binding energy depends on the chosen
transition basis and cannot be substituted for a converged coupled BSE result.

A reproducible 2 nm CdSe comparison needs the exact geometry and diameter
definition, surface/passivation, pseudopotential and basis, DFT functional,
spin-orbit convention, frontier orbitals and energy reference, dielectric
frequency convention, QP model, BSE active-space convergence, and independent
reference data. A fixed 10 Å radius in an interpolation formula is only an
arithmetic illustration; QDEX derives its effective radius from geometry.

The audit's numerical comparison and remaining limitations are recorded in
``audit/SCIENTIFIC_AUDIT.md``. The implementation sources are
``qdex.hardness.MATERIAL_DB`` and ``qdex.hardness.estimate_gw_qp_gap``;
use CLI ``--material CDSE --qp_gap sgw-anchor`` or YAML ``system.material``
and ``physics.qp_gap``.

Available 2.0 nm example
------------------------

The example in ``tests/CdSe/2.0nm`` (MO file stored as ``MOs_cleaned_20ang.txt.gz``)
contains a 149-atom, chloride-passivated geometry (68 Cd, 55 Se, 26 Cl), a
CP2K MO file, 2,753 AO basis functions, configurations and saved output.
Its QDEX provenance file reports a geometry-derived effective radius of
9.221 Å and a DFT frontier gap of 1.4568 eV. The DIM result at
``eps_out=2.4`` reports a 1.6721 eV QP opening, a 3.1289 eV resulting
gap, and HOMO/LUMO allocation of 55.16%/44.84%. These are model outputs
for this passivated structure, not independent GW reference values.

The saved audit result file records a different DIM fraction (61.13%/38.87%)
and a Resta fraction (52.81%/47.19%) from its own calls. This difference
illustrates why the fraction needs its full configuration, run version and
screening convention. The same file reports 0.2774 eV diagonal and
0.4756 eV coupled BSE binding relative to its 1.7255 eV QP transition,
so the diagonal result captures 58.3% in that run. This is a calculation
record, not a converged physical binding-energy benchmark.

The example's ``audit_tests/run_audit_tests.py`` records observations rather
than making eight pass/fail assertions. In particular, its stored result for
the solvent pathology is a stale fixed string, and its “discarded in CLI”
field is hard-coded ``true`` despite the current CLI wiring repair. Use the
actual printed values and current source when assessing those two fixes.


Results of the 25 September 2026 run
------------------------------------

Spin-free, 25 × 25 active space, dense diagonalization, ``eps_out = 1`` unless
stated. Full tables and discussion: ``audit/AUDIT_2026-09-25_QP_EXCITED_STATES.md``.
Reference values are enforced by the opt-in test
``QDEX_RUN_CDSE=1 pytest tests/test_cdse_integration.py``.

.. list-table:: Quasiparticle model (Resta direct kernel)
   :header-rows: 1

   * - ``qp_gap``
     - QP gap (eV)
     - S\ :sub:`1` (eV)
     - binding (eV)
   * - ``gw``
     - 3.861
     - 3.634
     - 0.227
   * - ``gw``, ``eps_out = 2.4``
     - 3.203
     - 2.977
     - 0.226
   * - ``sgw-resta``
     - 3.762
     - 3.535
     - 0.227
   * - ``sgw-dim``
     - 3.593
     - 3.366
     - 0.227
   * - ``qsgw-dim``
     - 3.955
     - 3.728
     - 0.227
   * - ``qsgw-resta``
     - 4.234
     - 4.005
     - 0.229

.. list-table:: Direct kernel at fixed QP gap (``gw``, 3.861 eV)
   :header-rows: 1

   * - kernel
     - S\ :sub:`1` (eV)
     - binding (eV)
   * - Resta
     - 3.634
     - 0.227
   * - DIM
     - 3.540
     - 0.321
   * - ``xs-resta``
     - 3.594
     - 0.267
   * - sBSE (atom)
     - 2.482
     - 1.379
   * - bare MNOK
     - 2.184
     - 1.677

For comparison, the first-exciton absorption of CdSe dots with
D ≈ 1.6–2.0 nm lies near 2.7–3.0 eV (Yu, Qu, Guo, Peng, Chem. Mater. 15, 2854
(2003) sizing curve; to be cross-checked with Aubert et al., Nano Lett. 22, 1778
(2022)). These runs predate the shared W; with it the solvent dependence of S₁ largely cancels
(:doc:`/excitons/cancellation`). Current results: :doc:`qp_bse_sweep` and
:doc:`anchor_evgw_benchmark`.
