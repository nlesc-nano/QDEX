Independent quasiparticles (``independent_qp``)
===============================================

Part of :doc:`/excitons/index`.

.. code-block:: yaml

   excitations:
     mode: independent_qp

**Energies.** Differences of QP energies, still with no electron–hole interaction:

.. math::

   \Omega_{ia} = \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}} ,

with every orbital corrected by the QP model (:doc:`/quasiparticles/models`). Intensities as for
:doc:`independent_dft`.

**What it describes.** An electron added to the dot and a hole created in it, far apart: the lowest
Ω is the fundamental gap, IP − EA. It is the right quantity for charging, transport and
photoemission, not for absorption.

**Why it is not the optical gap.** The QP energies contain the full self-energy of ΔW, including the
surface polarization (the self-image of each carrier). In a neutral excitation that polarization is
compensated by the electron–hole attraction, which this framework does not have. The energies are
therefore too high, by about 1.5 eV at 2 nm in vacuum for CdSe, and they follow the solvent as
strongly as the QP gap (about 0.9 eV between vacuum and toluene).

**What it is for.** Checking the QP correction itself (gap, IP/EA, level ordering) before adding the
interaction.

**Cost.** The QP step only.
