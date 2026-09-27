Configuration
=============

Part of :doc:`/quasiparticles/index`. The QP correction is set in the ``quasiparticles`` section, the
dielectric environment in ``environment`` and the integral representation in ``integrals``. All keys
are listed in :doc:`/getting_started/configuration`.

Example for a size series in toluene:

.. code-block:: yaml

   environment:
     eps_out: 2.24

   quasiparticles:
     model: sgw-resta          # or sgw-dim, evgw-*, qsgw-*
     # z: derived              # defaults: plasmon-pole Z, one-shot ΔCOHSEX,
     # selfenergy: cohsex      # all orbitals, sphere solvent term,
     # levels: orbital         # anchor residual scaled by E_conf
     # solvent_term: sphere
     # anchor_residual: on

   integrals:
     representation: mnok      # or xs

Calibrating a model on the anchor cluster (vacuum) stores its residual:

.. code-block:: bash

   qdex --config config.yaml --eps-out 1.0 --qp-anchor-calibrate
