Gap iteration
=============

Part of :doc:`/quasiparticles/index`.

.. important::

   The ``evgw-*`` implementation iterates an effective gap and screening model; it does not update a complete Green function or the state-resolved GW self-energy.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.estimate_evgw_dim_qp_gap``
* CLI: ``--qp_gap, --dynamic_z``
* YAML: ``physics.qp_gap, physics.dynamic_z``

.. code-block:: python

   estimate_evgw_dim_qp_gap(coords, atom_symbols, material_name=None, eps_out=2.4, C_occ_low=None, C_virt_low=None, eps_occ=None, eps_virt=None, atom_ao_ranges=None, alpha=1.0, max_iter=25, tol=0.0001, damping=0.5, return_details=False)


Effective-gap Self-Consistent Screening Model (``evgw-dim``, ``evgw-resta``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Standard one-shot :math:`G_0W_0` depends on the starting DFT functional. The **Effective-gap Self-Consistent Screening Model** (:math:`evGW`) model eliminates this bias by updating the quasiparticle eigenvalues in the Green's function iteratively:

.. math::

   \varepsilon_p^{(n+1)} = \varepsilon_p^{\mathrm{DFT}} + \Delta \varepsilon_p^{\mathrm{QP}}\big( \{\varepsilon_q^{(n)}\} \big)

until :math:`|\varepsilon_p^{(n+1)} - \varepsilon_p^{(n)}| < 10^{-4}\text{ eV}` (typically 3–5 iterations).


Implementation in QDEX (``evgw-dim``, ``evgw-resta``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Eigenvalue self-consistency is implemented in :func:`qdex.hardness.estimate_evgw_dim_qp_gap` and :func:`qdex.hardness.estimate_sgw_resta_qp_gap` (with ``self_consistent=True``):

- Runs an iterative Dyson loop up to ``max_iter=25`` with a default energy convergence tolerance of ``tol=1e-4`` eV.
- In each cycle :math:`n`, updates the Penn core dielectric permittivity :math:`\epsilon_{\mathrm{eff}}(E_g^{(n)})`, recomputes the screening contrast :math:`\Delta W^{(n)}`, and updates dynamic renormalization factors :math:`Z_H^{(n)}, Z_L^{(n)}`.
- Stabilizes iteration dynamics via linear damping: :math:`E_g^{(n+1)} = (1 - \alpha) E_g^{(n)} + \alpha E_g^{\mathrm{target}}` with mixing factor :math:`\alpha = 0.5`.
- Successfully converges in 3–5 iterations with negligible computational overhead (:math:`< 0.5\text{ s}`).
