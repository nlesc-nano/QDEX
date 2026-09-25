Qp edges
========

Part of :doc:`/workflows/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.cli``
* Callable: ``qdex.cli.main``
* CLI: ``--config``
* YAML: ``system.*, physics.*, namd.*``

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

   physics:
     qp_gap: "gw"
     material: "CSPBBR3"
     eps_out: 2.25
     qp_regularization_length: 1.0
     qp_residual_power: 2.0

   fuzzy:
     run: true
     dashboard_energy_mode: "both"
     qp_energy_reference: "vacuum"
     pdos_atoms: ["Pb", "Br"]


Execution
~~~~~~~~~

.. code-block:: bash

   qdex --config tutorial3_gw.yaml


Key Output
~~~~~~~~~~

The terminal prints the reconstructed absolute band edges:

.. code-block:: text

   [Absolute Band Edges (IP & EA)]
     Raw CP2K HOMO    :  -3.2140 eV (Floating Vacuum)
     Modeled PBE HOMO :  -5.4210 eV (anchor-reconstructed)
     -> Shift Split   : HOMO takes 43.1%, LUMO takes 56.9%
     QP HOMO (IP)     :  -5.8920 eV   -> IP = 5.8920 eV
     QP LUMO (EA)     :  -3.4110 eV   -> EA = 3.4110 eV

Both ``fuzzy_dashboard_dft.html`` and ``fuzzy_dashboard_qp.html`` are created, allowing direct side-by-side comparison of DFT and quasiparticle band structures.

