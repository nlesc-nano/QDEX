CdSe: comparison with experiment
================================

Part of :doc:`/validation/index`.

This page compares the current QP and excitation models with the measured first-exciton energy of
colloidal CdSe dots. All numbers were produced with the code as it is, with the default settings:

* bulk QSGW correction (Δ_bulk = 1.57 eV for CdSe, valence share f_b = 41.2 %);
* shared W in the QP correction and in K\ :sup:`d`;
* one-shot ΔCOHSEX for all orbitals, plasmon-pole Z;
* sphere reaction field for the environment, with the SAXS radius;
* MNOK integrals with γ_AA = IP − EA and exponent 2;
* no term fitted to a cluster calculation.

The raw results are in ``benchmarks/results/cdse_validation_2026-09-27.csv``.

Systems and settings
--------------------

.. list-table::
   :header-rows: 1

   * - Folder
     - Cluster
     - SAXS diameter d (formula-unit, hull)
     - KS (PBE) gap
   * - ``tests/CdSe/1.2nm``
     - Cd₁₆Se₁₃Cl₆ (35 atoms)
     - 1.12 nm (1.15, 1.06)
     - 2.639 eV
   * - ``tests/CdSe/2.0nm``
     - Cd₆₈Se₅₅Cl₂₆ (149 atoms, 2,753 basis functions)
     - 1.93 nm (1.87, 1.84)
     - 1.457 eV

d is the diameter a SAXS measurement would give: the Debye intensity of the inorganic atoms (Cd, Se
and the Cl surface) fitted with a homogeneous-sphere form factor (``qdex.cluster_size``). The sizing
curves use the SAXS core diameter. Organic ligands would be invisible to SAXS and are excluded in the
same way. Every QDEX run prints this size block.

MNOK integrals, Mulliken charges, 25 × 25 active space, dense diagonalization.
"Toluene" is ε_out = 2.24. With SOC, the first bright state is the lowest state with at least 10 % of
the largest oscillator strength.

Experimental references
-----------------------

``benchmarks/experimental_sizing.yaml``, evaluated at the SAXS diameter (for 2.4 and 3.2 nm at the
nominal diameter; these dots are not part of this page yet):

.. list-table::
   :header-rows: 1

   * - Reference
     - d = 1.93 nm
     - d = 2.4 nm
     - d = 3.2 nm
   * - Aubert, Hens et al., Nano Lett. 22, 1778 (2022), zinc blende
     - 3.31
     - 2.89
     - 2.46
   * - same, wurtzite
     - 3.04
     - 2.69
     - 2.35
   * - same, all CdSe
     - 3.22
     - 2.83
     - 2.43
   * - Yu, Qu, Guo, Peng, Chem. Mater. 15, 2854 (2003)
     - 2.77
     - 2.45
     - 2.22

At 1.93 nm both curves are extrapolated: the Hens fits start at 2.65–2.88 nm and the Yu fit at
450 nm (2.76 eV). The spread between them, 2.77–3.31 eV, comes mostly from the size calibration (TEM
for Yu, SAXS for Hens) and is the realistic experimental window. With the formula-unit diameter
(1.87 nm) the window would be 2.82–3.38 eV. There is no sizing reference at 1.12 nm; that cluster is
compared with evGW instead (:doc:`evgw_cluster`).

Results at 2.0 nm
-----------------

.. list-table::
   :header-rows: 1

   * - Model
     - QP gap vac / tol
     - S₁ vac / tol (spin-free)
     - bright S₁, tol, SOC
   * - ``bulk`` + ``sbse`` (Resta)
     - 3.027 / 3.027
     - 2.779 / 2.779
     - 2.652
   * - ``bulk`` + ``sbse`` (DIM)
     - 3.027 / 3.027
     - 2.675 / 2.675
     -
   * - ``bulk`` + ``stda``, a_x = 1/ε∞
     - 3.027 / 3.027
     - 2.699 / 2.699
     - 2.576
   * - ``none`` + ``stda``, a_x = 0 (PBE)
     - 1.457 / 1.457
     - 1.501 / 1.501
     - 1.391
   * - ``brus`` + Resta
     - 3.713 / 3.713
     - 3.465 / 3.465
     -
   * - ``sgw-resta``
     - 4.561 / 3.707
     - 3.038 / 2.964
     - 2.833
   * - ``sgw-dim``
     - 4.623 / 3.771
     - 3.070 / 2.997
     - 2.851
   * - ``evgw-resta``
     - 4.760 / 3.791
     - 3.114 / 2.998
     -
   * - ``evgw-dim``
     - 4.600 / 3.753
     - 3.059 / 2.989
     -
   * - ``qsgw-resta``
     - 4.826 / 3.830
     - 3.125 / 3.009
     -
   * - ``qsgw-dim``
     - 4.603 / 3.753
     - 3.053 / 2.986
     - 2.833

**Against experiment (toluene, bright state with SOC).** The Resta, DIM and qsGW models give
2.83–2.85 eV, inside the window 2.77–3.31 eV, 0.06–0.08 eV above the Yu curve and 0.19–0.48 eV below
the Hens curves. The sBSE gives 2.65 eV, 0.12 eV below the window, and the sTDA with a_x = 1/ε∞
2.58 eV, 0.19 eV below it. ``brus`` is far too high: the effective-mass kinetic term overestimates
confinement at this size.

**What the table shows.**

* **All ΔW models agree on S₁ within 0.05 eV**, although their QP gaps spread over 0.12 eV in toluene
  and 0.27 eV in vacuum. The extra QP opening of evGW and qsGW is compensated by a stronger K\ :sup:`d`.
* **Solvent.** The QP gap drops by 0.85–1.00 eV from vacuum to toluene; S₁ by 0.07–0.12 eV
  (:doc:`/excitons/cancellation`).
* **SOC** lowers the bright state by 0.13–0.15 eV. With SOC the lowest state is dark, 30–45 meV below
  the bright one (band-edge fine structure).
* **ΔW models vs sBSE.** The ΔW models are 0.19–0.23 eV above the sBSE. This is what does not cancel between
  the QP correction and K\ :sup:`d`: the image multipoles and the non-classical screened exchange. The
  sBSE, which drops both, misses it.
* **Resta vs DIM** differ by less than 0.03 eV in the ΔW models. In the sBSE, where only the bulk
  kernel changes, DIM binds 0.10 eV more.
* **sTDA** (:doc:`/excitons/stda`). With PBE orbitals the faithful setting is a_x = 0: no electron–hole
  attraction, so S₁ is the KS gap plus a small exchange term (1.50 eV), far below experiment. With the
  bulk GW gap and a_x = 1/ε∞, S₁ lies 0.08 eV below the sBSE. sTDA screens the short range more than
  Resta but leaves the 1/R tail unscreened; at 2 nm it binds 0.33 eV against 0.25 eV for the sBSE, and
  the gap to the sBSE will grow with size. sTDA is Grimme's method as published only with hybrid MO
  files, which were not used here. sTDA uses its own Grimme integrals and is not affected by the MNOK
  settings (:doc:`/integrals/representation`).

Results at 1.2 nm
-----------------

.. list-table::
   :header-rows: 1

   * - Model
     - QP gap vac / tol
     - S₁ vac / tol (spin-free)
     - bright vac / tol (spin-free)
     - bright S₁, tol, SOC
   * - ``bulk`` + ``sbse`` (Resta)
     - 4.209 / 4.209
     - 3.323 / 3.323
     - 3.617 / 3.617
     - 3.551
   * - ``bulk`` + ``sbse`` (DIM)
     - 4.209 / 4.209
     - 3.137 / 3.137
     - 3.442 / 3.442
     -
   * - ``bulk`` + ``stda``, a_x = 1/ε∞
     - 4.209 / 4.209
     - 3.808 / 3.808
     -
     - 3.668
   * - ``none`` + ``stda``, a_x = 0 (PBE)
     - 2.639 / 2.639
     - 2.748 / 2.748
     -
     - 2.626
   * - ``brus`` + Resta
     - 5.926 / 5.926
     - 5.040 / 5.040
     - 5.333 / 5.333
     -
   * - ``sgw-resta``
     - 6.767 / 5.407
     - 3.454 / 3.406
     - 3.819 / 3.746
     - 3.305
   * - ``sgw-dim``
     - 6.635 / 5.265
     - 3.476 / 3.426
     - 3.845 / 3.771
     - 3.347
   * - ``evgw-resta``
     - 7.124 / 5.578
     - 3.402 / 3.386
     - 3.785 / 3.733
     -
   * - ``evgw-dim``
     - 6.619 / 5.251
     - 3.476 / 3.427
     - 3.846 / 3.771
     -
   * - ``qsgw-resta``
     - 7.197 / 5.628
     - 3.420 / 3.397
     - 3.752 / 3.718
     -
   * - ``qsgw-dim``
     - 6.626 / 5.250
     - 3.437 / 3.408
     - 3.809 / 3.755
     - 3.362

At this size S₁ is dark in the spin-free calculation: the exchange splits the band-edge manifold and
the first bright state lies 0.30–0.37 eV above S₁. The rest of the pattern holds:

* the ΔW models agree on S₁ within 0.04 eV and on the bright state within 0.05 eV, while their QP
  gaps spread over 0.38 eV in toluene;
* they lie 0.06–0.10 eV (S₁) and 0.10–0.15 eV (bright) above the sBSE;
* the sTDA with a_x = 1/ε∞ lies 0.49 eV above the sBSE here (it binds 0.40 eV against 0.89 eV: at this
  size the on-site exchange and the strongly screened short range dominate, at 2 nm the unscreened tail);
* S₁ moves by 0.02–0.05 eV with the solvent against 1.4–1.6 eV for the QP gap.

How the results moved toward experiment
---------------------------------------

S₁ at 2 nm in toluene (spin-free, ``sgw-resta``) through the successive versions of the QP layer:

.. list-table::
   :header-rows: 1
   :widths: 46 12 12 30

   * - Version
     - QP gap
     - S₁
     - Change
   * - QP correction with its own W, bulk Resta kernel in the BSE (Δ_bulk = 1.27 eV)
     - 3.253
     - 3.026
     - agreement by error cancellation
   * - Shared W, classical self-energy (Δ_bulk = 1.27 eV)
     - 3.385
     - 2.515
     - −0.51: the missing mutual image restored in K\ :sup:`d`
   * - Shared W, ΔCOHSEX, sphere reaction field, anchor residual (Δ_bulk = 1.27 eV)
     - 3.211
     - 2.478
     - −0.04
   * - Bulk QSGW (Δ_bulk = 1.57 eV), no anchor residual
     - 3.734
     - 2.999
     - +0.52
   * - SAXS radius for the dielectric sphere
     - 3.681
     - 2.968
     - −0.03
   * - Current: MNOK on-site γ_AA = IP − EA (benchmarked against exact integrals)
     - 3.707
     - 2.964
     - −0.004

* **The first row matched experiment for the wrong reason.** The QP gap contained the surface
  polarization and the kernel did not contain the matching electron–hole image, so the whole
  self-image ended up in S₁, which also followed the solvent by 0.5 eV.
* **Sharing W removed that error.** S₁ became nearly solvent-independent, but fell 0.3–0.4 eV below
  experiment. With a consistent static W, S₁ ≈ KS gap + Δ_bulk − binding with bulk screening, so the
  remaining error had to be in Δ_bulk or in terms that do not cancel.
* **Bulk QSGW** raises Δ_bulk by 0.30 eV. Δ_bulk enters S₁ in full: +0.30 eV.
* **Dropping the anchor residual** removes −0.54 eV × s(R), with s = 0.41 at this size: +0.22 eV. The
  residual had been fitted to evGW\@PBE0 at 1.2 nm and was scaled to other sizes by assumption.
* **Together**, +0.52 eV moves the bright state with SOC from about 2.35 eV to 2.87 eV, into the
  experimental window.
* **SAXS radius.** The dielectric sphere now has the SAXS radius of the dot (9.67 Å instead of the hull
  value 9.22 Å): the image terms shrink, the QP gap drops by 0.05 eV (0.10 eV in vacuum) and S₁ by
  0.03 eV, to 2.84 eV for the bright state with SOC. The experimental references are evaluated at the
  same SAXS diameter.
* **MNOK on-site value.** The on-site integral is now IP − EA (twice the hardness η), the choice that
  reproduces exact (μμ|νν) integrals (:doc:`/integrals/representation`). At 2 nm S₁ barely moves
  (−0.004 eV spin-free, −0.01 eV for the bright state with SOC, now 2.83 eV), because the direct
  terms are dominated by the long range. At 1.2 nm the stronger on-site exchange and attraction lower
  S₁ by 0.5 eV and the bright state by 0.2 eV.

Limits
------

* **Two sizes only.** The 2.4 and 3.2 nm dots (6k and 13k basis functions) are the decisive test of
  the size dependence; their targets are in the table above (Hens 2.69–2.89 and 2.35–2.46 eV).
* **Active space.** S₁ decreases by about 0.06 eV from 25 × 25 to 100 × 100 at 2 nm.
* **Extrapolated references** at 1.93 nm, and a 0.5 eV spread between them. The size definition
  matters at the 0.05 eV level: the formula-unit diameter (1.87 nm) raises the references by 0.05–0.07 eV.
* **The QP gap is not validated by experiment here.** At 1.2 nm it is compared with evGW\@PBE0 in
  :doc:`evgw_cluster`; the ΔW models are 0.6–1.2 eV above it.

Rerun with ``benchmarks/qp_bse_sweep.py`` (:doc:`qp_bse_sweep`) or with the configs in
``tests/CdSe``.
