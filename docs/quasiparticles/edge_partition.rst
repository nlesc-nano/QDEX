Edge partition
==============

Part of :doc:`/quasiparticles/index`.

.. important::

   Frontier fractions are model allocations of a gap correction. Signed solvent contributions are evaluated separately; the fractions are bounded and default to 50:50 when the relevant contrast vanishes. Absolute IP/EA still require a vacuum reference or a declared monomer-anchor reconstruction.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.estimate_gw_qp_gap``
* CLI: ``--qp_gap, --dynamic_z``
* YAML: ``physics.qp_gap, physics.dynamic_z``

.. code-block:: python

   estimate_gw_qp_gap(coords, atom_symbols, material_name, eps_out, return_details=False, regularization_length_ang=1.0, residual_power=2.0, strict=False)


Approach B: Microscopic Wavefunction Asymmetry (:math:`f_H^{\mathrm{micro}}, f_L^{\mathrm{micro}}`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A major advantage of the microscopic :math:`\Delta W` models (``sgw-dim``, ``sgw-resta``, ``evgw``, ``qsgw``) over the two-anchor model is that they **completely skip the need for calibrated monomer anchor clusters**.

Instead of inheriting empirical fractions from a database, the state-specific self-energies are computed directly from the actual frontier wavefunctions:

.. math::

   \Delta \varepsilon_{\mathrm{HOMO}}^{\mathrm{conf}} = Z_H \, \sigma_H^{\mathrm{stat}} = \frac{1}{2} Z_H \sum_{A, B=1}^{N_{\mathrm{atoms}}} q_A^{\mathrm{HOMO}} \, \Delta W_{AB} \, q_B^{\mathrm{HOMO}}

.. math::

   \Delta \varepsilon_{\mathrm{LUMO}}^{\mathrm{conf}} = Z_L \, \sigma_L^{\mathrm{stat}} = \frac{1}{2} Z_L \sum_{A, B=1}^{N_{\mathrm{atoms}}} q_A^{\mathrm{LUMO}} \, \Delta W_{AB} \, q_B^{\mathrm{LUMO}}.

Because the HOMO (typically anion :math:`p` orbitals) and LUMO (typically cation :math:`s` orbitals) possess distinct spatial delocalizations and chemical hardnesses, their dielectric self-energies are naturally asymmetric. ``QDEX`` defines the **microscopic asymmetry fractions**:

.. math::

   f_H^{\mathrm{micro}} = \frac{\Delta \varepsilon_{\mathrm{HOMO}}^{\mathrm{conf}}}{\Delta \varepsilon_{\mathrm{HOMO}}^{\mathrm{conf}} + \Delta \varepsilon_{\mathrm{LUMO}}^{\mathrm{conf}}}, \quad
   f_L^{\mathrm{micro}} = \frac{\Delta \varepsilon_{\mathrm{LUMO}}^{\mathrm{conf}}}{\Delta \varepsilon_{\mathrm{HOMO}}^{\mathrm{conf}} + \Delta \varepsilon_{\mathrm{LUMO}}^{\mathrm{conf}}}

subject to :math:`f_H^{\mathrm{micro}} + f_L^{\mathrm{micro}} = 1.0`.

The bulk shift :math:`\Delta_{\mathrm{bulk}} = E_g^{\mathrm{GW, bulk}} - E_g^{\mathrm{PBE, bulk}}` is then partitioned using the same model-dependent microscopic fractions:

.. math::

   \varepsilon_i^{\mathrm{QP}} = \varepsilon_i^{\mathrm{DFT}} - f_H^{\mathrm{micro}} \Delta_{\mathrm{total}} \quad (i \in \mathrm{occ})

.. math::

   \varepsilon_a^{\mathrm{QP}} = \varepsilon_a^{\mathrm{DFT}} + f_L^{\mathrm{micro}} \Delta_{\mathrm{total}} \quad (a \in \mathrm{virt})

where :math:`\Delta_{\mathrm{total}} = \Delta_{\mathrm{bulk}} + \Delta \varepsilon_{\mathrm{HOMO}}^{\mathrm{conf}} + \Delta \varepsilon_{\mathrm{LUMO}}^{\mathrm{conf}}`.

This provides:

1. **Anchor-Free Evaluation**: Can be applied to any material, core/shell geometry, or surface ligand shell without requiring pre-computed vacuum Wulff clusters.
2. **Consistent Level Alignment**: Accurately predicts absolute Ionization Potentials (IP) and Electron Affinities (EA) tailored to the specific nanocrystal shape and surface termination.
3. **Model bulk projector allocation**: Partitions the bulk reference Hamiltonian :math:`\mathbf{H}_{\mathrm{bulk}} = -f_H^{\mathrm{micro}} \Delta_{\mathrm{bulk}} (0.5 \mathbf{P}_{\mathrm{occ}}) + f_L^{\mathrm{micro}} \Delta_{\mathrm{bulk}} \mathbf{Q}_{\mathrm{virt}}`.
