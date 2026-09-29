Configuration
=============

Part of :doc:`/quasiparticles/index`. All keys are listed in :doc:`/getting_started/configuration`.

Example for a size series in toluene:

.. code-block:: yaml

   environment:
     eps_out: 2.24

   quasiparticles:
     model: sgw-resta          # none, brus, sgw-resta, sgw-dim, evgw-*, qsgw-*
     # z: derived              # defaults: plasmon-pole Z of the same dielectric model,
     # selfenergy: cohsex      # one-shot ΔCOHSEX for all orbitals,
     # levels: orbital         # sphere reaction field for the environment
     # solvent_term: sphere
