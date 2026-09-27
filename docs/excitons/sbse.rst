sBSE without finite-size QP terms (``sbse``, ``diagonal_sbse``)
===============================================================

Part of :doc:`/excitons/index`.

.. code-block:: yaml

   quasiparticles:
     model: bulk               # PBE energies + bulk GW correction, no ΔW (PBE orbitals only)

   excitations:
     mode: sbse                # or diagonal_sbse
     kernel: resta             # or dim: bulk W

**Energies.** The solvers are those of :doc:`bse` (``sbse``) and :doc:`diagonal_bse`
(``diagonal_sbse``). What changes is the input: ΔW is removed from both the orbital energies and the
kernel,

.. math::

   A_{ia,jb} = \big[\varepsilon_a - \varepsilon_i + \Delta_{\mathrm{bulk}}\big]\,\delta_{ij}\delta_{ab}
   + 2K^x_{ia,jb} - (ij|W^{\mathrm{bulk}}|ab).

**Why it works.** By the Delerue–Lannoo–Allan result (PRL 84, 2457 (2000); PRL 90, 076803 (2003)),
the finite-size self-energy that ΔW adds to the QP gap and the extra attraction that ΔW adds to
K\ :sup:`d` cancel in a neutral excitation to leading order (:doc:`cancellation`, section 1). Removing
ΔW from both sides keeps S₁ and avoids computing it. The confinement is carried by the KS orbitals
themselves.

**What it misses.** The parts of ΔW that do not cancel (:doc:`cancellation`, sections 2, 3 and 5):

* the image multipoles, about +0.17 eV in S₁ at 2 nm for CdSe in vacuum;
* the non-classical screened exchange of a small dot;
* the (weak) solvent dependence of S₁, which is zero by construction.

**What it is for.** Optical spectra of large dots, where these terms are small, at the cost of a bulk
kernel. The QP energies of the run are only KS + bulk correction; they are not the QP gap of the dot.
