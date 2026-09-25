Soc
===

Part of :doc:`/workflows/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.cli``
* Callable: ``qdex.cli.main``
* CLI: ``--config``
* YAML: ``system.*, physics.*, namd.*``

.. code-block:: python

   main()


Tutorial 2: Incorporating Relativistic Spin-Orbit Coupling
----------------------------------------------------------

Lead (:math:`Z = 82`) exhibits intense relativistic spin-orbit coupling. In this tutorial, we compute the two-component spinor electronic structure.


Configuration File (``tutorial2_soc.yaml``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   system:
     mo_file: "CsPbBr3_MOs.mbse"
     xyz: "CsPbBr3_QD.xyz"
     basis_txt: "BASIS_MOLOPT"
     basis_name: "DZVP-MOLOPT-SR-GTH"
     cif: "CsPbBr3_bulk.cif"

   physics:
     soc: true
     soc_window_ev: 8.0

   fuzzy:
     run: true
     pdos_atoms: ["Pb", "Br"]
     ewin: [-4.0, 4.0]

   cube:
     export: true
     nhomos: 2
     nlumos: 2


Execution
~~~~~~~~~

.. code-block:: bash

   qdex --config tutorial2_soc.yaml


What to Observe
~~~~~~~~~~~~~~~

1. Notice the **conduction band splitting**: the Pb :math:`6p` conduction band splits into a lower :math:`j = 1/2` doublet and an upper :math:`j = 3/2` quartet.
2. The fundamental band gap contracts by :math:`\approx 0.65\text{ eV}` compared to Tutorial 1.
3. ``spinor_sp_HOMO_density.cube`` and ``spinor_sp_LUMO_density.cube`` are exported, containing the relativistic spinor densities :math:`\rho(\mathbf{r}) = |\psi^\alpha|^2 + |\psi^\beta|^2`.

