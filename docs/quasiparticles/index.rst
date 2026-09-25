Quasiparticle models
====================

QDEX corrects Kohn–Sham energies to quasiparticle (QP) energies with models derived from GW, and
passes the same screened interaction W to the BSE. Reading order:

1. :doc:`theory`: from full GW–BSE to quantum dots, why W must be shared between QP and BSE, and
   why each approximation is made.
2. :doc:`models`: the formulas of every QP model as implemented (``sgw-resta``, ``sgw-dim``, their
   ``evgw``/``qsgw`` variants and the two-anchor ``gw``).
3. :doc:`anchor`: the evGW anchor, the dielectric-sphere polarization and the residual calibration.
4. :doc:`dynamic_z`: the quasiparticle weight Z and its use in the BSE kernel.
5. :doc:`configuration` and :doc:`cost`: keys, defaults and scaling.

.. toctree::
   :maxdepth: 1

   theory
   models
   anchor
   dynamic_z
   configuration
   cost
