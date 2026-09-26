Why ΔW in the QP gap and in K\ :sup:`d` does not cancel completely
==================================================================

Part of :doc:`/excitons/index`.

.. figure:: /_static/figures/dielectric_sphere.svg
   :width: 100%
   :alt: dielectric sphere

   Dielectric confinement. A single added carrier is repelled by its own surface polarization (the QP
   gap opens). In a neutral exciton this self-energy and the polarization-enhanced electron–hole
   attraction largely cancel.

For one dominant transition h → e, the lowest exciton is

.. math::

   S_1 \approx \underbrace{(\varepsilon_e - \varepsilon_h) + \Delta_{\mathrm{bulk}}}_{\text{KS + bulk}}
   + \underbrace{Z_e\Delta\Sigma_e - Z_h\Delta\Sigma_h}_{\text{QP: }\Delta W}
   + \underbrace{r_L s - r_H s}_{\text{residual}}
   - \underbrace{(hh|W^{\mathrm{bulk}}|ee)}_{\text{bulk binding}}
   - \underbrace{\bar Z\,(hh|\Delta W|ee)}_{\text{K}^d\text{: }\Delta W}
   + 2(he|v|he).

ΔW appears twice with opposite signs. It raises the QP gap through the self-energies and increases
the binding through the direct term. This page shows when the two cancel and what is left.

1. Where the cancellation is exact
----------------------------------

In the classical limit (:doc:`/quasiparticles/gw`, section 6):

.. math::

   \Delta\Sigma_e - \Delta\Sigma_h = \tfrac12\langle\Delta W(\mathbf r,\mathbf r)\rangle_e
   + \tfrac12\langle\Delta W(\mathbf r,\mathbf r)\rangle_h
   \qquad\text{(self-images)},

.. math::

   (hh|\Delta W|ee) = \iint |\psi_h(\mathbf r)|^2\,\Delta W(\mathbf r,\mathbf r')\,|\psi_e(\mathbf r')|^2
   \qquad\text{(mutual image)}.

For a **constant** ΔW = c, both equal c and S₁ does not change. The Born (l = 0) term of the sphere
reaction field is such a constant inside the dot. It shifts the QP gap by
:math:`(1/\epsilon_{\mathrm{out}} - 1/\epsilon_\infty)\,e^2/R`, which is the largest part of the
solvent dependence, and it drops out of S₁ exactly.

Everything below is what makes ΔW differ from a constant, or makes the two terms differ in some other
way.

2. Image multipoles (l ≥ 1)
---------------------------

The higher multipoles of the reaction field grow toward the surface as :math:`(rr'/R^2)^l
P_l(\cos\theta)`.

* **In the self-image** (r = r′), :math:`P_l(1) = 1` and the terms add up.
* **In the mutual image** the electron and hole sit at different points, and :math:`P_l(\cos\theta)`
  averages toward zero.

The self-image is therefore larger than the mutual image, and S₁ keeps a positive classical size
effect of about 0.1 e²/R. For CdSe with 1S envelopes in vacuum:

.. list-table::
   :header-rows: 1

   * -
     - 1.2 nm
     - 2.0 nm
   * - self-images added to the QP gap
     - 2.54 eV
     - 1.46 eV
   * - mutual image added to the binding
     - 2.28 eV
     - 1.29 eV
   * - net effect on S₁
     - +0.26 eV
     - +0.17 eV

This part is classical, and every shared-W model contains it: ``gw`` through ``resta-sphere``, and
Resta and DIM through W_add. The older softened Born form (``qp_solvent_term: born``) has no
multipoles, so it cancels completely.

3. Non-classical screened exchange
----------------------------------

ΔSEX (:doc:`/quasiparticles/gw`, section 5) is

.. math::

   \Delta\mathrm{SEX}_n = -\sum_{m\in\mathrm{occ}}(nm|\Delta W|mn).

It equals the classical −⟨ΔW(r, r)⟩ only if the occupied states are complete for ΔW. In a small dot
they are not. The difference

* depends on how each orbital overlaps the occupied manifold, which differs strongly between the
  HOMO and the LUMO;
* has **no counterpart in K**\ :sup:`d`, which contains only the densities of the electron and the
  hole, not the other occupied states.

It therefore goes straight into S₁. At 1.2 nm (``sgw-resta``), ΔSEX is −3.01 eV for the HOMO and
−0.22 eV for the LUMO. It raises the QP gap by 0.33 eV over the classical value, and S₁ by 0.43 eV
before calibration. It is the physical reason why a microscopic ΔW (Resta, DIM) changes S₁ and a
purely classical one does not. The ``gw`` model has no such term; its non-classical part is only in
the residual A.

4. Terms that are only in the QP gap
------------------------------------

* **Δ_bulk.** The bulk GW correction is not ΔW. Its binding counterpart is the bulk W in K\ :sup:`d`,
  which is small (0.2–0.3 eV at 2 nm).
* **The anchor residual** :math:`(r_L - r_H)\,s(R)`. By construction it holds what the static
  ΔCOHSEX misses: band stretching, dynamical effects beyond Z, and the PBE → PBE0 starting point of
  the reference. None of it is represented in the kernel. For CdSe it lowers S₁ by 0.54 eV at 1.2 nm
  and 0.22 eV at 2 nm (``sgw-resta``).

5. The exciton is not one transition
------------------------------------

The QP levels contain ΔΣ in full. K\ :sup:`d` acts only within the active space and is averaged over
the correlated exciton :math:`\sum X_{ia}|ia\rangle`.

* **Correlation.** Mixing of transitions lets the electron and hole correlate their positions, which
  changes the mutual image relative to the product of 1S densities.
* **Finite active space.** Missing transitions leave part of the binding out. S₁ moves by 0.06 eV
  from 25 × 25 to 100 × 100 at 2 nm.

6. Z
----

The QP shift of each orbital is scaled by its own Z_n, the kernel by the frontier average Z̄. The
difference is second order, because :math:`|\Delta\Sigma_H| \approx |\Delta\Sigma_L|` (1.46 and 1.45 eV at
1.2 nm) and hence Z_H ≈ Z_L.

7. What is left
---------------

At 2 nm (CdSe, calibrated models), going from vacuum to toluene:

.. list-table::
   :header-rows: 1

   * -
     - QP gap (eV)
     - S₁ (eV)
   * - ``gw`` + ``resta-sphere``
     - −0.93
     - −0.04
   * - ``sgw-resta``
     - −0.91
     - −0.09
   * - ``sgw-dim``
     - −0.91
     - −0.09
   * - ``evgw-resta``
     - −1.00
     - −0.12
   * - ``qsgw-dim``
     - −0.91
     - −0.09

* **About 90 % of the solvent shift of the QP gap cancels in S₁.**
* **The rest is physics of the model.** The multipoles (section 2) are in all models. The screened
  exchange of the image term (section 3) is only in the Resta and DIM models, which is why their
  solvent shift is larger than that of ``gw``.

8. When W is not shared
-----------------------

If the QP correction contains ΔW but the kernel is the bulk W, the whole self-image ends up in S₁.
For CdSe at 2 nm in vacuum that adds 1.3 eV, and S₁ then follows the solvent almost as strongly as
the QP gap. Improving only the classical W does not help either: a size-dependent ε(R), atomistic
dipoles or a smoother boundary change the self- and the mutual image almost equally. What changes S₁
beyond "KS gap + bulk correction" is sections 2–4, not the details of the classical screening.
