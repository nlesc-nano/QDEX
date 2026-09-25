Exciton presets
===============

Part of :doc:`/workflows/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.cli``
* Callable: ``qdex.cli.main``
* CLI: ``--config``
* YAML: ``system.*, physics.*, namd.*``

.. code-block:: python

   main()


9. Recommended Workflow Presets
-------------------------------


Preset 1: Standard Colloidal QD Absorption Spectrum
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The recommended default for production calculations on colloidal nanocrystals:

.. code-block:: yaml

   physics:
     excitation_mode: "bse"
     2e-integrals: "mnok"              # Fast semi-empirical atom-centered representation
     kernel: "resta"                   # Thomas-Fermi electronic screening
     charge_type: "mulliken"           # "mulliken" or "lowdin"
     qp_gap: "sgw-dim"                 # Microscopic polarizable dipole QP gap
     dynamic_z: true                   # Empirical state-dependent Z_p damping
     nhomos: 50
     nlumos: 50

   bse:
     nroots: 15
     full_diag: false
     tol: 1.0e-5


Preset 2: Ultrafast Non-Adiabatic Molecular Dynamics (NAMD)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Optimized for high-speed trajectory propagation across thousands of frames:

.. code-block:: yaml

   physics:
     excitation_mode: "diagonal_bse"   # Diagonal e-h attraction without CI mixing
     2e-integrals: "mnok"
     kernel: "resta"
     qp_gap: "sgw-anchor"              # Sub-millisecond two-anchor scaled GW
     nhomos: 30
     nlumos: 30


Preset 3: Benchmark First-Principles Calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Exact analytical Gaussian two-electron integrals with microscopic RPA screening:

.. code-block:: yaml

   physics:
     excitation_mode: "bse"
     2e-integrals: "xs"                # Exact 4-center Gaussian integrals via Libint2
     kernel: "xs-rpa"                  # Parameter-free microscopic ZDO-RPA screening
     qp_gap: "qsgw-dim"                # Full AO orbital relaxation (qsGW)
     update_orbitals: true
     nhomos: 25
     nlumos: 25

   bse:
     nroots: 10
     full_diag: false


Preset 4: Spin-Orbit Coupling & Dark Excitons
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Relativistic 2-component spinor Bethe-Salpeter calculation:

.. code-block:: yaml

   physics:
     excitation_mode: "bse"
     2e-integrals: "mnok"
     kernel: "resta"
     soc: true                         # Enable 2-component spinor Hamiltonian
     soc_window_ev: 8.0                # Active window around Fermi level in eV
     qp_gap: "sgw-anchor"
     nhomos: 40
     nlumos: 40

   bse:
     nroots: 20
     full_diag: false

