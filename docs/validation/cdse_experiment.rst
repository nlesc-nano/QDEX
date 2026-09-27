CdSe: comparison with experiment
================================

Part of :doc:`/validation/index`.

This page compares the current QP and excitation models with the measured first-exciton energy of
colloidal CdSe dots. All numbers were produced with the code as it is, with the default settings:

* bulk QSGW correction (Δ_bulk = 1.57 eV for CdSe, valence share f_b = 41.2 %);
* shared W in the QP correction and in K\ :sup:`d`;
* one-shot ΔCOHSEX for all orbitals, plasmon-pole Z;
* sphere reaction field for the environment;
* no term fitted to a cluster calculation.

The raw results are in ``benchmarks/results/cdse_validation_2026-09-27.csv``.

Systems and settings
--------------------

.. list-table::
   :header-rows: 1

   * - Folder
     - Cluster
     - Core diameter d
     - KS (PBE) gap
   * - ``tests/CdSe/1.2nm``
     - Cd₁₆Se₁₃Cl₆ (35 atoms)
     - 1.15 nm
     - 2.639 eV
   * - ``tests/CdSe/2.0nm``
     - Cd₆₈Se₅₅Cl₂₆ (149 atoms, 2,753 basis functions)
     - 1.87 nm
     - 1.457 eV

d is the volume-equivalent diameter of the inorganic core, the quantity SAXS measures and the one
the sizing curves use. MNOK integrals, Mulliken charges, 25 × 25 active space, dense diagonalization.
"Toluene" is ε_out = 2.24. With SOC, the first bright state is the lowest state with at least 10 % of
the largest oscillator strength.

Experimental references
-----------------------

``benchmarks/experimental_sizing.yaml``, evaluated at the core diameter:

.. list-table::
   :header-rows: 1

   * - Reference
     - d = 1.87 nm
     - d = 2.4 nm
     - d = 3.2 nm
   * - Aubert, Hens et al., Nano Lett. 22, 1778 (2022), zinc blende
     - 3.39
     - 2.89
     - 2.46
   * - same, wurtzite
     - 3.11
     - 2.69
     - 2.35
   * - same, all CdSe
     - 3.30
     - 2.83
     - 2.43
   * - Yu, Qu, Guo, Peng, Chem. Mater. 15, 2854 (2003)
     - 2.83
     - 2.45
     - 2.22

At 1.87 nm both curves are extrapolated: the Hens fits start at 2.65–2.88 nm and the Yu fit at
450 nm (2.76 eV). The spread between them, 2.83–3.39 eV, comes mostly from the size calibration (TEM
for Yu, SAXS for Hens) and is the realistic experimental window. There is no sizing reference at
1.15 nm; that cluster is compared with evGW instead (:doc:`evgw_cluster`).

Results at 2.0 nm
-----------------

.. list-table::
   :header-rows: 1

   * - Model
     - QP gap vac / tol
     - S₁ vac / tol (spin-free)
     - bright S₁, tol, SOC
   * - ``none`` + ``sbse`` (Resta)
     - 3.027 / 3.027
     - 2.800 / 2.800
     - 2.680
   * - ``none`` + ``sbse`` (DIM)
     - 3.027 / 3.027
     - 2.706 / 2.706
     -
   * - ``brus`` + Resta
     - 3.856 / 3.856
     - 3.629 / 3.629
     -
   * - ``sgw-resta``
     - 4.642 / 3.734
     - 3.091 / 2.999
     - 2.873
   * - ``sgw-dim``
     - 4.694 / 3.788
     - 3.115 / 3.025
     - 2.897
   * - ``evgw-resta``
     - 4.806 / 3.801
     - 3.138 / 3.021
     -
   * - ``evgw-dim``
     - 4.676 / 3.773
     - 3.109 / 3.019
     -
   * - ``qsgw-resta``
     - 4.854 / 3.830
     - 3.140 / 3.025
     -
   * - ``qsgw-dim``
     - 4.677 / 3.773
     - 3.101 / 3.016
     - 2.912

**Against experiment (toluene, bright state with SOC).** The Resta, DIM and qsGW models give
2.87–2.91 eV, inside the window 2.83–3.39 eV, 0.05–0.08 eV above the Yu curve and 0.2–0.5 eV below
the Hens curves. The sBSE gives 2.68 eV, 0.15 eV below the window. ``brus`` is far too high: the
effective-mass kinetic term overestimates confinement at this size.

**What the table shows.**

* **All ΔW models agree on S₁ within 0.04 eV**, although their QP gaps spread over 0.17 eV in toluene
  and 0.21 eV in vacuum. The extra QP opening of evGW and qsGW is compensated by a stronger K\ :sup:`d`.
* **Solvent.** The QP gap drops by 0.9–1.0 eV from vacuum to toluene; S₁ by 0.09–0.12 eV
  (:doc:`/excitons/cancellation`).
* **SOC** lowers the bright state by 0.11–0.13 eV.
* **ΔW models vs sBSE.** The ΔW models are 0.2 eV above the sBSE. This is what does not cancel between
  the QP correction and K\ :sup:`d`: the image multipoles and the non-classical screened exchange. The
  sBSE, which drops both, misses it.
* **Resta vs DIM** differ by less than 0.03 eV in the ΔW models. In the sBSE, where only the bulk
  kernel changes, DIM binds 0.09 eV more.

Results at 1.2 nm
-----------------

.. list-table::
   :header-rows: 1

   * - Model
     - QP gap vac / tol
     - S₁ vac / tol (spin-free)
     - bright S₁, tol, SOC
   * - ``none`` + ``sbse`` (Resta)
     - 4.209 / 4.209
     - 3.756 / 3.756
     - 3.645
   * - ``none`` + ``sbse`` (DIM)
     - 4.209 / 4.209
     - 3.620 / 3.620
     -
   * - ``brus`` + Resta
     - 6.220 / 6.220
     - 5.766 / 5.766
     -
   * - ``sgw-resta``
     - 6.875 / 5.424
     - 4.025 / 3.941
     - 3.799
   * - ``sgw-dim``
     - 6.776 / 5.317
     - 4.032 / 3.945
     - 3.823
   * - ``evgw-resta``
     - 7.165 / 5.557
     - 4.011 / 3.939
     -
   * - ``evgw-dim``
     - 6.763 / 5.305
     - 4.030 / 3.944
     -
   * - ``qsgw-resta``
     - 7.227 / 5.596
     - 4.018 / 3.943
     -
   * - ``qsgw-dim``
     - 6.770 / 5.305
     - 4.012 / 3.936
     - 3.826

The same pattern holds, larger:

* the ΔW models agree on S₁ within 0.02 eV while their QP gaps spread over 0.46 eV;
* they lie 0.27 eV above the sBSE;
* S₁ moves by 0.08 eV with the solvent against 1.5 eV for the QP gap.

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
   * - Current: bulk QSGW (Δ_bulk = 1.57 eV), no anchor residual
     - 3.734
     - 2.999
     - +0.52

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

Limits
------

* **Two sizes only.** The 2.4 and 3.2 nm dots (6k and 13k basis functions) are the decisive test of
  the size dependence; their targets are in the table above (Hens 2.69–2.89 and 2.35–2.46 eV).
* **Active space.** S₁ decreases by about 0.06 eV from 25 × 25 to 100 × 100 at 2 nm.
* **Extrapolated references** at 1.87 nm, and a 0.5 eV spread between them.
* **The QP gap is not validated by experiment here.** At 1.2 nm it is compared with evGW\@PBE0 in
  :doc:`evgw_cluster`; the ΔW models are 0.7–1.2 eV above it.

Rerun with ``benchmarks/qp_bse_sweep.py`` (:doc:`qp_bse_sweep`) or with the configs in
``tests/CdSe``.
