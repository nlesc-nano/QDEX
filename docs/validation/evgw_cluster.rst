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

Models (vacuum, spin-free, 25 × 25, MNOK integrals with γ_AA = IP − EA; Δ = model − evGW). At this
size S₁ is dark; the first bright state is given after it:

.. list-table::
   :header-rows: 1

   * - Model
     - QP gap
     - Δ gap
     - HOMO shift
     - LUMO shift
     - S₁ / bright
   * - ``bulk``
     - 4.209
     - −1.82
     - −0.646
     - +0.924
     - 3.323 / 3.617 (sBSE)
   * - ``sgw-resta``
     - 6.767
     - +0.74
     - −1.956
     - +2.171
     - 3.454 / 3.819
   * - ``sgw-dim``
     - 6.635
     - +0.61
     - −1.856
     - +2.140
     - 3.476 / 3.845
   * - ``evgw-resta``
     - 7.124
     - +1.10
     - −2.190
     - +2.294
     - 3.402 / 3.785
   * - ``evgw-dim``
     - 6.619
     - +0.59
     - −1.843
     - +2.137
     - 3.476 / 3.846
   * - ``qsgw-resta``
     - 7.197
     - +1.17
     - —
     - —
     - 3.420 / 3.752
   * - ``qsgw-dim``
     - 6.626
     - +0.60
     - —
     - —
     - 3.437 / 3.809

**What it shows.**

* **The ΔW models open the gap 0.59–1.17 eV more than evGW\@PBE0.** Most of the excess is on the HOMO
  (0.45–0.79 eV too deep); the LUMO is 0.14–0.30 eV too high.
* **Where it comes from.**

  - Static COHSEX overestimates gaps (no dynamical screening beyond one plasmon pole).
  - The bulk QSGW Δ_bulk is 0.30 eV larger than the G₀W₀-type value used before.
  - The reference itself starts from PBE0 and is not converged with respect to self-consistency.
* **S₁ does not inherit the QP spread.** All ΔW models give S₁ = 3.40–3.48 eV and a bright state at 3.75–3.85 eV, because the extra opening
  is matched by a stronger K\ :sup:`d` (:doc:`/excitons/cancellation`).

The ΔW models are therefore not a substitute for GW when the absolute QP gap or IP/EA of a small
cluster is needed. For optical energies the error largely cancels.
