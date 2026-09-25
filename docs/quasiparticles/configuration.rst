Configuration
=============

Part of :doc:`/quasiparticles/index`.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.estimate_gw_qp_gap``
* CLI: ``--qp_gap, --dynamic_z``
* YAML: ``physics.qp_gap, physics.dynamic_z``

.. code-block:: python

   estimate_gw_qp_gap(coords, atom_symbols, material_name, eps_out, return_details=False, regularization_length_ang=1.0, residual_power=2.0, strict=False)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part3_gw_scissor/index.rst:1147-1149``

15. CLI Flags & YAML Configuration Reference
--------------------------------------------


.. rubric:: From ``docs/part3_gw_scissor/index.rst:1150-1187``

Command-Line Arguments
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - CLI Flag
     - Default
     - Description
   * - ``--qp_gap <choice>``
     - ``brus``
     - Quasiparticle model: ``sgw-anchor`` (alias ``gw``), ``sgw-dim``, ``sgw-resta``, ``evgw-dim``, ``evgw-resta``, ``qsgw-dim``, ``qsgw-resta``, ``brus``, ``pbe``, or numeric gap in eV.
   * - ``--dynamic_z``
     - ``False``
     - Compute state-dependent dynamic renormalization :math:`Z_p` from the plasmon-pole :math:`f`-sum rule.
   * - ``--update_orbitals``
     - ``False``
     - Run the static AO-basis COHSEX-like orbital-relaxation model.
   * - ``--material <name>``
     - ``DEFAULT``
     - Material key in ``MATERIAL_DB`` (e.g. ``CSPBBR3``, ``CDSE``, ``INAS``).
   * - ``--eps-out <float>``
     - ``2.0``
     - Optical dielectric constant :math:`\epsilon_{\mathrm{out}}` of surrounding solvent or matrix.
   * - ``--qp-regularization-length <float>``
     - ``1.0``
     - Short-range regularization length :math:`\ell` (in Å) in the anchor confinement formula.
   * - ``--qp-residual-power <float>``
     - ``2.0``
     - Exponent :math:`p` for the quantum confinement power-law decay.
   * - ``--dashboard_energy_mode <choice>``
     - ``dft``
     - Energy axis for Fuzzy Band and PDOS dashboards: ``dft``, ``qp``, or ``both``.
   * - ``--qp_energy_reference <choice>``
     - ``vacuum``
     - Reference zero for QP spectra: ``vacuum`` (absolute IP/EA) or ``fermi`` (:math:`E_F = 0`).


.. rubric:: From ``docs/part3_gw_scissor/index.rst:1188-1205``

YAML Configuration Reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   physics:
     qp_gap: "sgw-anchor"      # "sgw-anchor", "sgw-dim", "evgw-dim", "qsgw-dim", "brus", "pbe"
     dynamic_z: true           # Apply empirical state-dependent Z_p damping
     update_orbitals: false    # Full AO orbital relaxation (qsGW)
     material: "CSPBBR3"
     eps_out: 2.25
     qp_regularization_length: 1.0
     qp_residual_power: 2.0

   fuzzy:
     run: true
     dashboard_energy_mode: "both"
     qp_energy_reference: "vacuum"
