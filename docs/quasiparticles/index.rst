Quasiparticles
==============

YAML sections: ``quasiparticles`` and ``environment``.

QDEX corrects the Kohn–Sham energies to quasiparticle (QP) energies with models derived from GW. The
screened interaction W built here is also the screened interaction of the excited states
(:doc:`/excitons/index`). The integrals are represented as set in :doc:`/integrals/index`.

1. :doc:`gw`: the full GW self-energy, how it is computed, and how the quantum-dot problem reduces
   to the finite-size correction ΔW on top of bulk GW (ΔCOHSEX).
2. :doc:`models`: the approximations for ΔW, from no correction (``none``) and the bulk correction (``bulk``, ``brus``) to the
   classical sphere and the atomistic Resta and DIM models, and why each is needed.
3. :doc:`dynamic_z`: the quasiparticle weight Z.
4. :doc:`configuration` and :doc:`cost`: keys and scaling.

.. toctree::
   :maxdepth: 1

   gw
   models
   dynamic_z
   configuration
   cost
