Quasiparticle models
====================

QDEX corrects Kohn–Sham energies to quasiparticle (QP) energies with models derived from GW. The
screened interaction W built here is also the screened interaction of the BSE
(:doc:`/excitons/index`). Reading order:

1. :doc:`gw`: the full GW self-energy, how it is computed, and how the quantum-dot problem reduces
   to the finite-size correction ΔW on top of bulk GW (ΔCOHSEX).
2. :doc:`models`: the approximations for ΔW, from no finite-size term (``brus``) to the two-anchor
   ``gw`` and the atomistic Resta and DIM models, and why each is needed.
3. :doc:`anchor`: the evGW anchor, the dielectric-sphere polarization and the residual calibration.
4. :doc:`dynamic_z`: the quasiparticle weight Z.
5. :doc:`representation`: how the integrals are represented, MNOK or ZDO xs, and the limits of each.
6. :doc:`configuration` and :doc:`cost`: keys, defaults and scaling.

.. toctree::
   :maxdepth: 1

   gw
   models
   anchor
   dynamic_z
   representation
   configuration
   cost
