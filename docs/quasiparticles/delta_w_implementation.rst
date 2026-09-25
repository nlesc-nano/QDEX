Delta w implementation
======================

Part of :doc:`/quasiparticles/index`.

.. important::

   DIM maps a three-field polarizability response to pairwise screening with additional normalization. The mapping is heuristic and must be benchmarked against a microscopic dielectric response.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.estimate_sgw_dim_qp_gap``
* CLI: ``--qp_gap, --dynamic_z``
* YAML: ``physics.qp_gap, physics.dynamic_z``

.. code-block:: python

   estimate_sgw_dim_qp_gap(coords, atom_symbols, material_name=None, eps_out=2.4, C_occ_low=None, C_virt_low=None, eps_occ=None, eps_virt=None, atom_ao_ranges=None, alpha=1.0, Z=0.8, dynamic_z=False, self_consistent=False, max_iter=25, tol=0.0001, damping=0.5, return_details=False)


Implementation in QDEX (``sgw-dim``, ``sgw-resta``, ``sgw``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The microscopic :math:`\Delta W` algorithms are implemented in :mod:`qdex.hardness`:

1. **Discrete Dipole Interaction Model (``sgw-dim``)**:
   Implemented in :func:`qdex.hardness.estimate_sgw_dim_qp_gap`:
   - Builds the :math:`3N \times 3N` Thole dipole interaction matrix :math:`\mathbf{M} = \boldsymbol{\alpha}^{-1} + \mathbf{T}` from the polarizability table ``POLARIZABILITY_TABLE_AU``.
   - Inverts :math:`\mathbf{M}` to solve for induced dipoles under 3D probe fields, obtaining atom-specific effective polarizabilities :math:`p_A^{\mathrm{eff}}`.
   - Computes pairwise screening factors :math:`S_{AB}` and assembles the screened potential :math:`W_{AB}^{\mathrm{DIM}} = S_{AB} \gamma_{AB}^{\mathrm{Ohno}}`.
   - Evaluates the confinement contrast matrix :math:`\Delta W_{AB} = \max(0, W_{AB}^{\mathrm{DIM}} - W_{AB}^{\mathrm{bulk}}) + W_{AB}^{\mathrm{solv}}`.
   - Contracts with Löwdin atomic populations :math:`\mathbf{q}_H, \mathbf{q}_L` to evaluate static shifts :math:`\sigma_H^{\mathrm{stat}}, \sigma_L^{\mathrm{stat}}`.
   - Evaluates :math:`f_H^{\mathrm{micro}}, f_L^{\mathrm{micro}}` (Approach B) and exports them in the ``provenance`` dictionary.

2. **Resta Screened Dielectric Model (``sgw-resta``)**:
   Implemented in :func:`qdex.hardness.estimate_sgw_resta_qp_gap`:
   - Determines the median core nearest-neighbor bond distance :math:`d_{\mathrm{NN}}`.
   - Evaluates the Penn-scaled core dielectric constant :math:`\epsilon_{\mathrm{eff}}(R)`.
   - Computes the Thomas-Fermi screening wavevector :math:`k_s = \sqrt{\epsilon_{\mathrm{eff}} - 1} / d_{\mathrm{NN}}`.
   - Assembles the Resta screened potential :math:`W_{AB}^{\mathrm{Resta}} = [c_\infty + (1 - c_\infty) e^{-k_s R_{AB}}] \gamma_{AB}^{\mathrm{Ohno}}`.
   - Computes :math:`\Delta W` and evaluates state shifts and Approach B splitting fractions.

3. **Simplified BSE Polarizability Model (``sgw``)**:
   Implemented in :func:`qdex.hardness.estimate_sgw_qp_gap`:
   - Evaluates the non-interacting transition polarizability matrix :math:`\boldsymbol{\Pi}^0` from active Kohn-Sham orbital pairs.
   - Evaluates the dynamic inversion :math:`\mathbf{W} = (\mathbf{I} + \mathbf{J}_{\mathrm{solv}} \boldsymbol{\Pi}^0)^{-1} \mathbf{J}_{\mathrm{solv}}`.

*CLI & YAML Invocation*:

.. code-block:: bash

   # DIM polarizable dipoles (recommended default)
   qdex --mos ground_state.mos --material CDSE --qp_gap sgw-dim --dynamic_z --eps-out 2.40

   # Resta Thomas-Fermi screening
   qdex --mos ground_state.mos --material CDSE --qp_gap sgw-resta --dynamic_z --eps-out 2.40

