Alignment
=========

Part of :doc:`/quasiparticles/index`.

.. important::

   The current CLI may reconstruct absolute edges from monomer entries in ``MATERIAL_DB`` even for an anchor-free gap model. Distinguish a computed gap from calibrated vacuum-level alignment.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.estimate_gw_qp_gap``
* CLI: ``--qp_gap, --dynamic_z``
* YAML: ``physics.qp_gap, physics.dynamic_z``

.. code-block:: python

   estimate_gw_qp_gap(coords, atom_symbols, material_name, eps_out, return_details=False, regularization_length_ang=1.0, residual_power=2.0, strict=False)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part3_gw_scissor/index.rst:601-618``

10. Rigid Scissor Operator & Frontier Splitting Strategies (Approach A vs. Approach B)
--------------------------------------------------------------------------------------

To predict absolute valence and conduction band edge alignments, the total quasiparticle gap scissor shift :math:`\Delta_{\mathrm{total}}` must be partitioned between occupied and virtual manifolds:

.. math::

   \varepsilon_i^{\mathrm{QP}} = \varepsilon_i^{\mathrm{DFT}} - f_{\mathrm{homo}} \, \Delta_{\mathrm{total}} \quad (i \in \mathrm{occ}), \qquad
   \varepsilon_a^{\mathrm{QP}} = \varepsilon_a^{\mathrm{DFT}} + f_{\mathrm{lumo}} \, \Delta_{\mathrm{total}} \quad (a \in \mathrm{virt})

where :math:`f_{\mathrm{homo}} + f_{\mathrm{lumo}} = 1.0`. By referencing eigenvalues to the vacuum zero, ``QDEX`` yields the true absolute **Ionization Potential (IP)** and **Electron Affinity (EA)**:

.. math::

   \mathrm{IP} = -\varepsilon_{\mathrm{HOMO}}^{\mathrm{QP}}, \qquad \mathrm{EA} = -\varepsilon_{\mathrm{LUMO}}^{\mathrm{QP}}.

``QDEX`` provides two distinct physical strategies for determining :math:`f_{\mathrm{homo}}` and :math:`f_{\mathrm{lumo}}`:


.. rubric:: From ``docs/part3_gw_scissor/index.rst:619-637``

Approach A: Database Monomer Anchor Frontier Splitting (Used by ``sgw-anchor``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the two-anchor model (``sgw-anchor``), individual self-energy corrections for the HOMO and LUMO levels are extracted from high-level vacuum GW calculations on the monomer anchor cluster stored in ``MATERIAL_DB``:

.. math::

   \delta_h = \varepsilon_{\mathrm{HOMO}}^{\mathrm{GW, mono}} - \varepsilon_{\mathrm{HOMO}}^{\mathrm{PBE, mono}} \le 0, \quad
   \delta_l = \varepsilon_{\mathrm{LUMO}}^{\mathrm{GW, mono}} - \varepsilon_{\mathrm{LUMO}}^{\mathrm{PBE, mono}} \ge 0.

The asymmetry fractions are:

.. math::

   f_{\mathrm{homo}}^{\mathrm{anchor}} = \frac{-\delta_h}{\delta_l - \delta_h}, \quad
   f_{\mathrm{lumo}}^{\mathrm{anchor}} = \frac{\delta_l}{\delta_l - \delta_h}.

*Applicability*: Fast and reliable for pristine, stoichiometric nanocrystals with tabulated material data.


.. rubric:: From ``docs/part3_gw_scissor/index.rst:638-668``

Approach B: Microscopic Wavefunction Asymmetry (Used by ``sgw-dim``, ``sgw-resta``, ``evgw``, ``qsgw``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the microscopic :math:`\Delta W` models, ``QDEX`` **completely skips the need for monomer anchor clusters** by calculating the self-energy shifts directly on the actual frontier wavefunctions of the quantum dot:

.. math::

   \Delta \varepsilon_{\mathrm{HOMO}}^{\mathrm{conf}} = \frac{1}{2} Z_H \sum_{A, B=1}^{N_{\mathrm{atoms}}} q_A^{\mathrm{HOMO}} \, \Delta W_{AB} \, q_B^{\mathrm{HOMO}}, \quad
   \Delta \varepsilon_{\mathrm{LUMO}}^{\mathrm{conf}} = \frac{1}{2} Z_L \sum_{A, B=1}^{N_{\mathrm{atoms}}} q_A^{\mathrm{LUMO}} \, \Delta W_{AB} \, q_B^{\mathrm{LUMO}}.

The state-specific asymmetry fractions are computed as:

.. math::

   f_H^{\mathrm{micro}} = \frac{\Delta \varepsilon_{\mathrm{HOMO}}^{\mathrm{conf}}}{\Delta \varepsilon_{\mathrm{HOMO}}^{\mathrm{conf}} + \Delta \varepsilon_{\mathrm{LUMO}}^{\mathrm{conf}}}, \quad
   f_L^{\mathrm{micro}} = \frac{\Delta \varepsilon_{\mathrm{LUMO}}^{\mathrm{conf}}}{\Delta \varepsilon_{\mathrm{HOMO}}^{\mathrm{conf}} + \Delta \varepsilon_{\mathrm{LUMO}}^{\mathrm{conf}}}.

Because the HOMO (typically anion :math:`p` orbitals) and LUMO (typically cation :math:`s` orbitals) possess distinct spatial delocalizations and chemical hardnesses, :math:`f_H^{\mathrm{micro}}` and :math:`f_L^{\mathrm{micro}}` capture the true physical asymmetry of the dot.

Furthermore, in the static AO-basis orbital-relaxation model (``qsgw-*``), this same model allocation is used to partition the bulk reference Hamiltonian:

.. math::

   \mathbf{H}_{\mathrm{bulk}} = -f_H^{\mathrm{micro}} \Delta_{\mathrm{bulk}} (0.5 \mathbf{P}_{\mathrm{occ}}) + f_L^{\mathrm{micro}} \Delta_{\mathrm{bulk}} \mathbf{Q}_{\mathrm{virt}}.

*Why Approach B is the Recommended Default*:

1. **No Anchor Database Required**: Works for any chemical composition, core/shell hetero-interface, facet termination, or organic ligand shell without requiring pre-computed monomer cluster data.
2. **First-Principles Consistency**: Dynamically evolves as nanocrystal size, shape, and aspect ratio change.
3. **Physical Soundness**: Reflects the actual orbital localization of the frontier states rather than an idealized small-molecule surrogate.


.. rubric:: From ``docs/part3_gw_scissor/index.rst:669-699``

Implementation in QDEX (Frontier Alignment & CLI)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Frontier level alignment is orchestrated between :mod:`qdex.hardness` and the main command-line entry point :mod:`qdex.cli`:

1. **Provenance Inspection**:
   When any microscopic model (``sgw-dim``, ``sgw-resta``, ``evgw-dim``, ``evgw-resta``, ``qsgw-dim``, ``qsgw-resta``) is selected, :mod:`qdex.hardness` returns ``f_homo_micro`` and ``f_lumo_micro`` in the ``provenance`` dictionary.
   
2. **Automatic Dispatching in CLI**:
   - If ``f_homo_micro`` is present, :mod:`qdex.cli` automatically activates **Approach B** and logs:

     .. code-block:: text

        [Absolute Band Edges (IP & EA - Microscopic Wavefunction Asymmetry)]
          -> Shift Split: HOMO takes 55.2%, LUMO takes 44.8% (SGW_DIM)

   - If ``sgw-anchor`` (or ``gw``) is selected, the CLI falls back to **Approach A**, reading :math:`f_{\mathrm{homo}}, f_{\mathrm{lumo}}` from ``MATERIAL_DB``.

3. **Rigid Scissor Application**:
   Occupied Kohn-Sham orbital energies are lowered by :math:`-f_{\mathrm{homo}} \Delta_{\mathrm{total}}` and virtual energies are raised by :math:`+f_{\mathrm{lumo}} \Delta_{\mathrm{total}}`:

   .. math::

      \varepsilon_i^{\mathrm{QP}} = \varepsilon_i^{\mathrm{DFT}} - f_{\mathrm{homo}} \Delta_{\mathrm{total}} \quad (i \le \mathrm{HOMO}), \qquad
      \varepsilon_a^{\mathrm{QP}} = \varepsilon_a^{\mathrm{DFT}} + f_{\mathrm{lumo}} \Delta_{\mathrm{total}} \quad (a \ge \mathrm{LUMO}).

4. **Absolute IP & EA Output**:
   When ``--qp_energy_reference vacuum`` is set (the default), the CLI reports the absolute Ionization Potential (:math:`\mathrm{IP} = -\varepsilon_{\mathrm{HOMO}}^{\mathrm{QP}}`) and Electron Affinity (:math:`\mathrm{EA} = -\varepsilon_{\mathrm{LUMO}}^{\mathrm{QP}}`), establishing direct contact with ultraviolet photoelectron spectroscopy (UPS) and cyclic voltammetry experiments.

---
