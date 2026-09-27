QP presets
==========

Part of :doc:`/workflows/index`. Each preset shows only the sections that differ from the example in
:doc:`/getting_started/configuration`. The BSE kernel is set by the QP model and is not repeated.

Size series of colloidal dots (production)
------------------------------------------

.. code-block:: yaml

   environment:
     eps_out: 2.24               # toluene

   quasiparticles:
     model: sgw-resta            # spherical dots; sgw-dim for shape, ligand or shell effects

   integrals:
     representation: mnok

Atom-resolved surface screening
-------------------------------

.. code-block:: yaml

   quasiparticles:
     model: sgw-dim

Eigenvalue self-consistent ΔW
-----------------------------

.. code-block:: yaml

   quasiparticles:
     model: evgw-resta           # or evgw-dim: ΔW iterated against the QP gap

Orbital relaxation
------------------

.. code-block:: yaml

   quasiparticles:
     model: qsgw-dim             # ΔCOHSEX as a matrix, diagonalized self-consistently

No finite-size correction (sBSE)
--------------------------------

KS energies plus the bulk GW correction, with the bulk W in the BSE:

.. code-block:: yaml

   quasiparticles:
     model: none

   excitations:
     mode: sbse
     kernel: resta               # or dim

Small clusters
--------------

.. code-block:: yaml

   integrals:
     representation: xs          # ZDO, exact (μμ|νν); more binding at short range
