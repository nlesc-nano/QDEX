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
effect of about 0.1 e²/R. For CdSe with 1S envelopes in vacuum (SAXS radius 5.6 and 9.7 Å):

.. list-table::
   :header-rows: 1

   * -
     - 1.2 nm
     - 2.0 nm
   * - self-images added to the QP gap
     - 2.41 eV
     - 1.39 eV
   * - mutual image added to the binding
     - 2.16 eV
     - 1.23 eV
   * - net effect on S₁
     - +0.25 eV
     - +0.16 eV

This part is classical, and the Resta and DIM models contain it through W_add. The older softened Born form (``quasiparticles.solvent_term: born``) has no
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

It therefore goes straight into S₁. At 1.2 nm (``sgw-resta``), ΔSEX is −2.74 eV for the HOMO and
−0.13 eV for the LUMO. It raises the QP gap by 0.26 eV over the classical value, and S₁ by 0.36 eV. It is the physical reason why a microscopic ΔW (Resta, DIM) changes S₁ and a
purely classical one does not.

4. The bulk correction
----------------------

Δ_bulk is not ΔW and is only in the orbital energies. Its binding counterpart is the bulk W in
K\ :sup:`d`, which is small (0.2–0.3 eV at 2 nm).

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
difference is second order, because :math:`|\Delta\Sigma_H| \approx |\Delta\Sigma_L|` (1.37 and 1.33 eV at
1.2 nm) and hence Z_H ≈ Z_L.

7. What is left
---------------

At 2 nm (CdSe), going from vacuum to toluene:

.. list-table::
   :header-rows: 1

   * -
     - QP gap (eV)
     - S₁ (eV)
   * - ``sgw-resta``
     - −0.86
     - −0.08
   * - ``sgw-dim``
     - −0.85
     - −0.07
   * - ``evgw-resta``
     - −0.95
     - −0.10
   * - ``qsgw-dim``
     - −0.85
     - −0.07

* **About 90 % of the solvent shift of the QP gap cancels in S₁.**
* **The rest is physics of the model.** The multipoles (section 2) are in all models. The screened
  exchange of the image term (section 3) adds to it; a classical sphere alone gives about −0.04 eV.

8. When W is not shared
-----------------------

If the QP correction contains ΔW but the kernel is the bulk W, the whole self-image ends up in S₁.
For CdSe at 2 nm in vacuum that adds 1.3 eV, and S₁ then follows the solvent almost as strongly as
the QP gap. Improving only the classical W does not help either: a size-dependent ε(R), atomistic
dipoles or a smoother boundary change the self- and the mutual image almost equally. What changes S₁
beyond "KS gap + bulk correction" is sections 2 and 3, not the details of the classical screening.
