Qp presets
==========

Part of :doc:`/workflows/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.cli``
* Callable: ``qdex.cli.main``
* CLI: ``--config``
* YAML: ``system.*, physics.*, namd.*``

.. code-block:: python

   main()


13. Recommended Workflow Presets
--------------------------------


Preset 1: Ultra-Fast NAMD Trajectory Dynamics (``sgw-anchor``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For non-adiabatic dynamics across 1,000–10,000 time steps, ``sgw-anchor`` provides sub-millisecond evaluation while reproducing experimental band gaps and solvatochromic shifts:

.. code-block:: yaml

   physics:
     qp_gap: "sgw-anchor"      # Two-anchor scaled GW (backward-compatible alias: "gw")
     material: "CSPBBR3"
     eps_out: 2.25
     excitation_mode: "diagonal_bse"
     kernel: "resta"
     2e-integrals: "mnok"


Preset 2: First-Principles Nanocrystal Screening (``sgw-dim``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For accurate single-point spectra without relying on calibrated monomer anchors, use atomistic polarizable dipoles with dynamic :math:`Z_p` from the plasmon-pole :math:`f`-sum rule:

.. code-block:: yaml

   physics:
     qp_gap: "sgw-dim"         # Microscopic Delta-W with Atomistic Polarizable Dipoles
     dynamic_z: true           # Empirical state-dependent Z_p damping
     material: "CDSE"
     eps_out: 2.40
     excitation_mode: "bse"
     kernel: "qp"              # same W as the QP model (required)
     2e-integrals: "mnok"


Preset 3: Eigenvalue Self-Consistency (``evgw-dim``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To remove starting-point DFT eigenvalue bias without the computational overhead of full orbital updates:

.. code-block:: yaml

   physics:
     qp_gap: "evgw-dim"        # Effective-gap self-consistent screening model with DIM screening
     material: "CDSE"
     eps_out: 2.40
     excitation_mode: "bse"
     kernel: "qp"              # same W as the QP model (required)
     2e-integrals: "mnok"


Preset 4: Full Quasiparticle Self-Consistency (``qsgw-dim``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For core/shell quantum dots, type-II heterojunctions, or surface-passivated dots where orbital polarization is critical:

.. code-block:: yaml

   physics:
     qp_gap: "qsgw-dim"        # Full AO-basis static orbital relaxation model
     update_orbitals: true     # Diagonalizes H_eff to relax wavefunctions
     dynamic_z: true
     material: "CDSE"
     eps_out: 2.40
     excitation_mode: "bse"
     kernel: "qp"              # same W as the QP model (required)
     2e-integrals: "mnok"

