Orbital iteration
=================

Part of :doc:`/quasiparticles/index`.

.. figure:: /_static/figures/qsgw_loop.svg
   :width: 100%
   :alt: qsgw loop

   The ``qsgw-*`` orbital-relaxation loop as implemented in ``estimate_qsgw_dim_qp_gap`` / ``estimate_qsgw_resta_qp_gap``.


.. important::

   The ``qsgw-*`` implementation relaxes orbitals with a static atom-block COHSEX-like operator. It is distinct from conventional QSGW. The CLI now passes its returned QP eigenvalues to the BSE solver with zero additional scissor, while retaining DFT eigenvalues separately for the DFT framework and filtering.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.estimate_qsgw_dim_qp_gap``
* CLI: ``--qp_gap, --dynamic_z``
* YAML: ``physics.qp_gap, physics.dynamic_z``

.. code-block:: python

   estimate_qsgw_dim_qp_gap(coords, atom_symbols, C, eps, S, atom_ao_ranges, homo_index, material_name=None, eps_out=2.4, alpha=1.0, dynamic_z=True, max_iter=25, tol=0.0001, damping=0.5, return_details=False)


Quasiparticle Self-Consistent GW (``qsgw-dim``, ``qsgw-resta``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When orbital wavefunctions undergo significant polarization (e.g. in core/shell nanocrystals, type-II heterostructures, or surface defect states), the single-particle wavefunctions must be relaxed.

``QDEX`` implements full AO-basis **static COHSEX-like orbital relaxation (``qsgw-*``)**:

1. Constructs the effective non-local Hamiltonian in the full :math:`N_{\mathrm{AO}} \times N_{\mathrm{AO}}` atomic orbital basis:

   .. math::

      \mathbf{H}_{\mathrm{eff}} = \mathbf{H}_{\mathrm{DFT}} + \Delta \mathbf{H}_{\mathrm{bulk}} + Z \left( \boldsymbol{\Sigma}^{\mathrm{SEX}}[\mathbf{P}] + \boldsymbol{\Sigma}^{\mathrm{COH}} \right)

   where:

   * :math:`\mathbf{P} = \mathbf{C}_{\mathrm{occ}} \mathbf{C}_{\mathrm{occ}}^T` is the single-particle density matrix.
   * :math:`\Sigma_{\mu \nu}^{\mathrm{SEX}} = -P_{\mu \nu} \, \Delta W_{AB}` for :math:`\mu \in A, \nu \in B`.
   * :math:`\Sigma_{\mu \nu}^{\mathrm{COH}} = \frac{1}{2} \, \Delta W_{AA} \, S_{\mu \nu}`.

2. Solves the generalized Hermitian eigenvalue problem:

   .. math::

      \mathbf{H}_{\mathrm{eff}} \, \mathbf{C}_{\mathrm{QP}} = \mathbf{S} \, \mathbf{C}_{\mathrm{QP}} \, \mathbf{E}_{\mathrm{QP}}.

3. Evaluates orbital fidelity to monitor wavefunction reconstruction:

   .. math::

      \mathcal{F}_p = |\langle \psi_p^{(0)} | \psi_p^{(n)} \rangle| = | (\mathbf{C}_p^{(0)})^\dagger \mathbf{S} \, \mathbf{C}_p^{(n)} |.


Implementation in QDEX (``qsgw-dim``, ``qsgw-resta``, ``--update_orbitals``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Full AO-basis orbital relaxation is implemented in :func:`qdex.hardness.estimate_qsgw_dim_qp_gap` and :func:`qdex.hardness.estimate_qsgw_resta_qp_gap`:

1. **Löwdin Symmetric Orthogonalization**: Constructs the transformation :math:`\mathbf{S}^{-1/2} = \mathbf{U}_S \boldsymbol{\Lambda}_S^{-1/2} \mathbf{U}_S^T` from the overlap matrix :math:`\mathbf{S}` and transforms initial DFT orbitals into the orthogonal Löwdin basis: :math:`\mathbf{C}_{\mathrm{Low}}^{(0)} = \mathbf{S}^{1/2} \mathbf{C}_{\mathrm{DFT}}`.
2. **Hamiltonian in Löwdin Basis**: Evaluates :math:`\mathbf{H}_{\mathrm{DFT}}^{\mathrm{Low}} = \mathbf{C}_{\mathrm{Low}}^{(0)} \operatorname{diag}(\boldsymbol{\varepsilon}^{\mathrm{DFT}}) (\mathbf{C}_{\mathrm{Low}}^{(0)})^T`.
3. **AO Screening Contrast**: Expands atom-by-atom screening contrast :math:`\Delta W_{AB}` to atomic orbital block pairs :math:`\Delta W_{\mu \nu} = \Delta W_{AB}` for :math:`\mu \in A, \nu \in B`.
4. **Iterative Self-Consistent Loop** (up to ``max_iter=30``):
   - Constructs single-particle density matrix :math:`\mathbf{P}^{\mathrm{Low}} = 2 \mathbf{C}_{\mathrm{occ}}^{\mathrm{Low}} (\mathbf{C}_{\mathrm{occ}}^{\mathrm{Low}})^T` and projector :math:`\mathbf{Q}^{\mathrm{Low}} = \mathbf{I} - 0.5 \mathbf{P}^{\mathrm{Low}}`.
   - Computes frontier charges :math:`\mathbf{q}_H, \mathbf{q}_L` and dynamic asymmetry fractions :math:`f_H^{\mathrm{micro}}, f_L^{\mathrm{micro}}` (Approach B).
   - Assembles bulk reference operator: :math:`\mathbf{H}_{\mathrm{bulk}}^{\mathrm{Low}} = -f_H^{\mathrm{micro}} \Delta_{\mathrm{bulk}} (0.5 \mathbf{P}^{\mathrm{Low}}) + f_L^{\mathrm{micro}} \Delta_{\mathrm{bulk}} \mathbf{Q}^{\mathrm{Low}}`.
   - Builds non-local self-energies: :math:`\boldsymbol{\Sigma}^{\mathrm{SEX}} = -0.5 \mathbf{P}^{\mathrm{Low}} \odot \Delta \mathbf{W}_{\mathrm{AO}}` and :math:`\boldsymbol{\Sigma}^{\mathrm{COH}} = 0.5 \operatorname{diag}(\Delta \mathbf{W}_{\mathrm{AO}})`.
   - Updates target effective Hamiltonian: :math:`\mathbf{H}_{\mathrm{eff}}^{\mathrm{target}} = \mathbf{H}_{\mathrm{DFT}}^{\mathrm{Low}} + \mathbf{H}_{\mathrm{bulk}}^{\mathrm{Low}} + Z (\boldsymbol{\Sigma}^{\mathrm{SEX}} + \boldsymbol{\Sigma}^{\mathrm{COH}})`.
   - Mixes with previous Hamiltonian using linear damping :math:`\mathbf{H}_{\mathrm{eff}} = (1 - \alpha) \mathbf{H}_{\mathrm{eff}}^{\mathrm{prev}} + \alpha \mathbf{H}_{\mathrm{eff}}^{\mathrm{target}}` (:math:`\alpha = 0.5`).
   - Diagonalizes :math:`\mathbf{H}_{\mathrm{eff}} \mathbf{C}^{\mathrm{Low}} = \mathbf{C}^{\mathrm{Low}} \boldsymbol{\varepsilon}^{\mathrm{QP}}` via LAPACK/NumPy ``eigh``.
   - Aligns eigenvector phases: :math:`C_{\mu p} \leftarrow C_{\mu p} \times \operatorname{sign}(\langle \psi_p^{\mathrm{prev}} | \psi_p \rangle)` to eliminate arbitrary sign flips.
   - Monitors wavefunction reconstruction fidelity: :math:`\mathcal{F}_p = |(\mathbf{C}_p^{(0)})^\dagger \mathbf{C}_p^{(n)}|^2`.
5. **Back-Transformation to AO Basis**: Converged eigenvectors are transformed back to the non-orthogonal AO basis: :math:`\mathbf{C}_{\mathrm{QP}} = \mathbf{S}^{-1/2} \mathbf{C}^{\mathrm{Low}}`.
6. **Integration with BSE**: In :mod:`qdex.cli`, the relaxed wavefunctions :math:`\mathbf{C}_{\mathrm{QP}}` and quasiparticle eigenvalues :math:`\boldsymbol{\varepsilon}^{\mathrm{QP}}` directly replace the DFT starting point for subsequent Bethe-Salpeter excited-state calculations!

*CLI & YAML Invocation*:

.. code-block:: bash

   # Full qsGW orbital update with DIM polarizable dipoles
   qdex --mos ground_state.mos --material CDSE --qp_gap qsgw-dim --update_orbitals --dynamic_z --eps-out 2.40

