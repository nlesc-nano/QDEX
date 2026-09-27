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
     - dielectric variant: bulk GW gap and a_x = 1/ε∞ (0.161 for CdSe)

``functional`` accepts pbe, blyp, bp86, tpss, r2scan (0), tpssh (0.10), b3lyp, b3pw91 (0.20), pbe0
(0.25), m06 (0.27), pw6b95 (0.28), bhlyp (0.50), m06-2x (0.54), hf (1). ``ax`` overrides it.

**Checks.** A QP model that defines its own W (``sgw-*``, ``evgw-*``, ``qsgw-*``) is rejected; ``bulk``
with a non-PBE functional and a_x = 0 are reported.

**Cd₁₆Se₁₃Cl₆** (PBE MO file, vacuum, 25 × 25, spin-free S₁; KS gap 2.64 eV):

.. list-table::
   :header-rows: 1

   * - Setting
     - S₁ (eV)
   * - ``none``, a_x = 0 (faithful for PBE)
     - 2.75
   * - ``none``, a_x = 0.25 (inconsistent: hybrid attraction on a PBE gap)
     - 1.84
   * - ``bulk``, a_x = 1/ε∞ (dielectric variant)
     - 3.81
   * - for comparison: sBSE (``bulk``, bulk Resta W)
     - 3.76

The dielectric variant lands close to the sBSE: both put bulk screening into K\ :sup:`d`, uniformly
(sTDA) or with the Resta distance dependence (sBSE). sTDA becomes genuinely different with hybrid MO
files.

**Not implemented yet:** Grimme's perturbative selection of high-energy transitions. The active space
is set by ``nhomos``, ``nlumos`` and ``e_thresh`` as for the other solvers.
