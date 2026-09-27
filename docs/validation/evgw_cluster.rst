Cd₁₆Se₁₃Cl₆ against evGW\@PBE0
==============================

Part of :doc:`/validation/index`.

For the 1.2 nm cluster there is an explicit GW calculation (CP2K: PBE → PBE0 → evGW, ``EV_GW_ITER 4``).
It is used here as a check of the QP models, not as a calibration: no model is fitted to it.

Reference (vacuum, eV; shifts relative to the PBE eigenvalues, which are the QDEX input):

.. list-table::
   :header-rows: 1

   * -
     - HOMO
     - LUMO
     - gap
   * - PBE
     - −6.420
     - −3.785
     - 2.634
   * - evGW\@PBE0
     - −7.818
     - −1.788
     - 6.029
   * - shift
     - −1.398
     - +1.997
     - +3.395

Models (vacuum, spin-free, 25 × 25; Δ = model − evGW):

.. list-table::
   :header-rows: 1

   * - Model
     - QP gap
     - Δ gap
     - HOMO shift
     - LUMO shift
     - S₁
   * - ``bulk``
     - 4.209
     - −1.82
     - −0.646
     - +0.924
     - 3.756 (sBSE)
   * - ``sgw-resta``
     - 6.875
     - +0.85
     - −1.988
     - +2.248
     - 4.025
   * - ``sgw-dim``
     - 6.776
     - +0.75
     - −1.913
     - +2.224
     - 4.032
   * - ``evgw-resta``
     - 7.165
     - +1.14
     - −2.173
     - +2.353
     - 4.011
   * - ``evgw-dim``
     - 6.763
     - +0.73
     - −1.902
     - +2.221
     - 4.030
   * - ``qsgw-resta``
     - 7.227
     - +1.20
     - —
     - —
     - 4.018
   * - ``qsgw-dim``
     - 6.770
     - +0.74
     - —
     - —
     - 4.012

**What it shows.**

* **The ΔW models open the gap 0.7–1.2 eV more than evGW\@PBE0.** Most of the excess is on the HOMO
  (0.5–0.8 eV too deep); the LUMO is 0.2–0.35 eV too high.
* **Where it comes from.**

  - Static COHSEX overestimates gaps (no dynamical screening beyond one plasmon pole).
  - The bulk QSGW Δ_bulk is 0.30 eV larger than the G₀W₀-type value used before.
  - The reference itself starts from PBE0 and is not converged with respect to self-consistency.
* **S₁ does not inherit the QP spread.** All ΔW models give 4.01–4.03 eV, because the extra opening
  is matched by a stronger K\ :sup:`d` (:doc:`/excitons/cancellation`).

The ΔW models are therefore not a substitute for GW when the absolute QP gap or IP/EA of a small
cluster is needed. For optical energies the error largely cancels.
