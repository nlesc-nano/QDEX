sTDA (``stda``, ``diagonal_stda``)
==================================

Part of :doc:`/excitons/index`.

.. code-block:: yaml

   quasiparticles:
     model: none               # DFT orbital energies as they are

   excitations:
     mode: stda                # or diagonal_stda
     functional: pbe0          # functional of the MO file -> a_x; or ax: 0.25

The simplified Tamm–Dancoff approximation of Grimme (J. Chem. Phys. 138, 244104 (2013)) with his
interaction parameters, taken from the ``std2`` source. The solvers are those of :doc:`bse` and
:doc:`diagonal_bse`.

**Matrix.**

.. math::

   A_{ia,jb} = (\varepsilon_a - \varepsilon_i)\,\delta_{ij}\delta_{ab} + 2\,(ia|jb)_K - (ij|ab)_J ,
   \qquad (\text{triplet: no } K \text{ term}),

with Löwdin transition charges q and two MNOK-type interactions:

.. math::

   \gamma^K_{AB} = \big(R_{AB}^{\alpha} + \eta_{AB}^{-\alpha}\big)^{-1/\alpha},\qquad
   \gamma^J_{AB} = \big(R_{AB}^{\beta} + (a_x\,\eta_{AB})^{-\beta}\big)^{-1/\beta},

.. math::

   \alpha = 1.42 + 0.48\,a_x,\qquad \beta = 0.20 + 1.83\,a_x,\qquad \eta_{AB} = \tfrac12(\eta_A + \eta_B).

* :math:`\eta_A` is twice the Ghosh–Islam hardness of ``HARDNESS_DICT``, i.e. (ii|ii) = IP − EA, in
  atomic units (R in bohr).
* Both interactions tend to 1/R at long range; γ\ :sup:`J` does so slowly (β ≈ 0.66 for a_x = 0.25).
* a_x enters only through γ\ :sup:`J`: the exponent β and the damping length a_x η. There is no
  further prefactor in TDA. For a_x = 0, γ\ :sup:`J` = 0 and there is no electron–hole attraction.
* sTDA is a monopole method; it runs only with ``integrals.representation: mnok``.

**a_x and the orbital energies.** a_x is the fraction of Fock exchange of the functional that produced
the MO file. It must match the orbital energies: a hybrid gap is already partly opened by exact
exchange, and a_x (ij|ab) is the attraction that belongs to it.

.. list-table::
   :header-rows: 1
   :widths: 28 22 50

   * - MO file
     - Setting
     - Meaning
   * - PBE0 (or another hybrid)
     - ``model: none``, ``functional: pbe0``
     - faithful sTDA (a_x = 0.25)
   * - PBE
     - ``model: none``, ``functional: pbe``
     - faithful sTDA with a_x = 0: K\ :sup:`x` only, no binding
   * - PBE
     - ``model: bulk``, ``ax: dielectric``
     - variant: bulk GW gap and a_x = 1/ε∞ (0.161 for CdSe) in the short range of γ\ :sup:`J`

``functional`` accepts pbe, blyp, bp86, tpss, r2scan (0), tpssh (0.10), b3lyp, b3pw91 (0.20), pbe0
(0.25), m06 (0.27), pw6b95 (0.28), bhlyp (0.50), m06-2x (0.54), hf (1). ``ax`` overrides it.

**Checks.** A QP model that defines its own W (``sgw-*``, ``evgw-*``, ``qsgw-*``) is rejected; ``bulk``
with a non-PBE functional and a_x = 0 are reported.

**What a_x does and does not do.** In TDA, a_x only sets the short range of γ\ :sup:`J`: the on-site
value is a_x η and the exponent β grows with a_x. At long range γ\ :sup:`J` → 1/R, the bare interaction,
for any a_x > 0. sTDA therefore has no long-range screening of the electron–hole attraction. That is
appropriate for molecules, but not for a dielectric: in a large dot the attraction should tend to
1/(ε∞ R). Setting a_x = 1/ε∞ weakens the short range only; it does not screen the tail.

**CdSe results** (PBE MO files, 25 × 25, spin-free S₁ in eV; QP model and a_x as indicated):

.. list-table::
   :header-rows: 1

   * - Setting
     - 1.2 nm
     - 2.0 nm
     - 2.0 nm, toluene, bright, SOC
   * - KS gap
     - 2.639
     - 1.457
     -
   * - ``none``, a_x = 0 (faithful for PBE)
     - 2.748
     - 1.501
     - 1.391
   * - ``none``, a_x = 0.25 (hybrid attraction on a PBE gap)
     - 1.838
     - 0.842
     -
   * - ``bulk``, a_x = 1/ε∞ = 0.161
     - 3.808
     - 2.699
     - 2.576
   * - sBSE: ``bulk``, bulk Resta W
     - 3.756
     - 2.800
     - 2.680
   * - ΔW models (``sgw-*``, ``qsgw-dim``), toluene
     - 3.90–3.91
     - 2.97–2.99
     - 2.84–2.87

* **a_x = 0 on PBE orbitals** gives the KS gap plus a small exchange term (+0.11 eV at 1.2 nm,
  +0.04 eV at 2 nm): no binding.
* **a_x = 0.25 on PBE orbitals** gives S₁ below the KS gap. It is the inconsistent combination: a
  hybrid-strength attraction on a gap that exact exchange has not opened.
* **Dielectric variant vs sBSE.** The binding relative to the bulk-corrected gap is 0.40 eV (sTDA) vs
  0.45 eV (sBSE) at 1.2 nm, but 0.33 vs 0.23 eV at 2 nm. sTDA screens the short range more (on-site
  a_x η instead of the unscreened Resta value) and the long range not at all. The long range grows in
  weight with size, so sTDA binds increasingly more than the sBSE in larger dots and does not reach the
  bulk exciton limit.
* **Solvent.** sTDA contains no environment term; S₁ is independent of ε_out.
* **Against experiment** (2 nm, window 2.77–3.31 eV): the dielectric sTDA gives 2.58 eV, 0.1 eV below
  the sBSE and 0.3 eV below the ΔW models.
* **diagonal_stda** recovers most of the coupled result here (2.72 vs 2.70 eV, dielectric, 2 nm).

For dots, sTDA with PBE orbitals is therefore a reference, not a replacement for the sBSE or the ΔW
models. With hybrid (PBE0) MO files it is Grimme's method as published.

**Transition selection.** Grimme's perturbative selection is available for all coupled solvers
(``excitations.selection: perturbative``; :doc:`bse`).
