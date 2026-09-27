Exciton analysis
================

Part of :doc:`/workflows/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.cli``
* Callable: ``qdex.cli.main``
* CLI: ``--config``
* YAML: ``system``, ``environment``, ``quasiparticles``, ``integrals``, ``excitations``, ``namd``

.. code-block:: python

   main()


Tutorial 5: Exciton Wavefunction Descriptors & NTO Analysis
------------------------------------------------------------

Here we perform deep spatial wavefunction analysis on the lowest exciton states using the Plasser-Dreuw framework and Natural Transition Orbitals.


Configuration File (``tutorial5_analysis.yaml``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   system:
     mo_file: "CsPbBr3_MOs.mbse"
     xyz: "CsPbBr3_QD.xyz"
     basis_txt: "BASIS_MOLOPT"
     basis_name: "DZVP-MOLOPT-PBE-GTH"

     material: "CSPBBR3"

   quasiparticles:
     model: sgw-resta

   excitations:
     mode: bse
     nhomos: 30
     nlumos: 30
     nroots: 10

   output:
     nto: true
     nto_states: [1, 2, 3]
     nto_top: 2
     nto_csv: true
     plot: true
     cube: true
     bse_states: [1]
     cube_spacing: 0.5


Execution
~~~~~~~~~

.. code-block:: bash

   qdex --config tutorial5_analysis.yaml


Analysis Results
~~~~~~~~~~~~~~~~

* Open ``exciton_analysis.html`` in your browser to explore the 6-panel interactive dashboard.
* Inspect ``nto_results.csv``: Root 1 typically shows :math:`\lambda_1^2 > 0.90`, confirming that a single NTO pair captures the optical transition.
* Render ``exciton_S1_hole.cube`` (red) and ``exciton_S1_elec.cube`` (green) in VMD to visualize the spatial extent of the exciton Bohr radius.

