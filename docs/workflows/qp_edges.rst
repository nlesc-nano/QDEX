Qp edges
========

Part of :doc:`/workflows/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.cli``
* Callable: ``qdex.cli.main``
* CLI: ``--config``
* YAML: ``system``, ``environment``, ``quasiparticles``, ``integrals``, ``excitations``, ``namd``

.. code-block:: python

   main()


Tutorial 3: Scaled GW Band Gap Correction & Absolute Edges
----------------------------------------------------------

Semi-local DFT underestimates the band gap. Here we activate the Scaled GW model with dielectric solvent screening (:math:`\epsilon_{\mathrm{out}} = 2.25` for toluene) and predict absolute IP and EA levels.


Configuration File (``tutorial3_gw.yaml``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   system:
     mo_file: "CsPbBr3_MOs.mbse"
     xyz: "CsPbBr3_QD.xyz"
     basis_txt: "BASIS_MOLOPT"
     basis_name: "DZVP-MOLOPT-PBE-GTH"
     cif: "CsPbBr3_bulk.cif"

     material: "CSPBBR3"

   environment:
     eps_out: 2.25

   quasiparticles:
     model: sgw-resta
     energy_reference: vacuum    # absolute IP/EA

   analysis:
     run_fuzzy: true
     dashboard_energy_mode: both
     pdos_atoms: ["Pb", "Br"]


Execution
~~~~~~~~~

.. code-block:: bash

   qdex --config tutorial3_gw.yaml


Key Output
~~~~~~~~~~

The terminal prints the absolute band edges of the dot:

.. code-block:: text

   [Absolute Band Edges (IP & EA)]
     PBE HOMO / LUMO  :  e_H / e_L eV (CP2K eigenvalues; on the vacuum scale only for PERIODIC NONE ...)
     QP shifts        : HOMO d_H eV, LUMO d_L eV (orbital-resolved QP levels)
     QP HOMO (IP)     :  e_H + d_H eV   -> IP = -(e_H + d_H) eV
     QP LUMO (EA)     :  e_L + d_L eV   -> EA = -(e_L + d_L) eV

The edges are the dot's own PBE eigenvalues plus the QP shift of each edge:

.. math::

   \mathrm{IP} = -\left(\varepsilon_\mathrm{HOMO}^\mathrm{PBE} + \Delta_\mathrm{HOMO}\right), \qquad
   \mathrm{EA} = -\left(\varepsilon_\mathrm{LUMO}^\mathrm{PBE} + \Delta_\mathrm{LUMO}\right).

A finite cluster computed with ``PERIODIC NONE`` and an isolated Poisson solver has its eigenvalues on
the vacuum scale (``MULTIPOLE`` is within about 0.04 eV of ``WAVELET`` and ``MT`` for the 1.2 nm CdSe
cluster in a 28 Å box), so dots of different sizes compare directly; keep a similar amount of vacuum
around each dot. A run with a periodic Poisson solver has an arbitrary eigenvalue zero, and its IP/EA are
meaningless (the 1.2 nm MO file in ``tests/CdSe/1.2nm`` lies rigidly 0.66 eV above the same cluster
computed with ``PERIODIC NONE`` and ``MULTIPOLE``).
With orbital-resolved QP levels (the Delta-W models, ``quasiparticles.levels: orbital``)
:math:`\Delta_\mathrm{HOMO}` and :math:`\Delta_\mathrm{LUMO}` are the model's own shifts of the two
orbitals; with rigid levels (``bulk``, ``levels: rigid``) the gap correction is split between the edges
with the model's HOMO fraction. Periodic runs have no vacuum level, and no IP/EA is assigned.

Both ``fuzzy_dashboard_dft.html`` and ``fuzzy_dashboard_qp.html`` are created, allowing direct side-by-side comparison of DFT and quasiparticle band structures.

