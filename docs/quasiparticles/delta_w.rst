Delta w
=======

Part of :doc:`/quasiparticles/index`.

.. important::

   These routines use static atom-centered screened-interaction contrasts and charging energies. They do not evaluate a dynamical GW self-energy; ``v_xc`` cancellation is an approximation and is least reliable near reconstructed or trapped surfaces.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.estimate_sgw_resta_qp_gap``
* CLI: ``--qp_gap, --dynamic_z``
* YAML: ``physics.qp_gap, physics.dynamic_z``

.. code-block:: python

   estimate_sgw_resta_qp_gap(coords, atom_symbols, material_name=None, eps_out=2.4, dft_gap=None, C_occ_low=None, C_virt_low=None, eps_occ=None, eps_virt=None, atom_ao_ranges=None, alpha=1.0, Z=0.8, penn_scaling=True, dynamic_z=False, self_consistent=False, max_iter=25, tol=0.0001, damping=0.5, return_details=False)


7. Avenue 2: Microscopic Dielectric Shift Model (:math:`\Delta W`)
------------------------------------------------------------------

When finite anchor clusters are unavailable or when state-specific orbital self-energies are required, the Microscopic :math:`\Delta W` Model computes the quasiparticle shift directly from the difference between the confined nanocrystal screened interaction :math:`W^{\mathrm{QD}}` and the bulk crystal screened interaction :math:`W^{\mathrm{bulk}}`.


The Physical Rationale: Cancellation of :math:`v_{xc}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In Hedin's :math:`G_0W_0` quasiparticle equation, the energy shift relative to DFT is:

.. math::

   \Delta \varepsilon_p^{\mathrm{QP}} = Z_p \, \langle \psi_p | \Sigma - v_{xc}^{\mathrm{DFT}} | \psi_p \rangle.

In semiconductors, the short-range exchange-correlation potential :math:`v_{xc}^{\mathrm{DFT}}(\mathbf{r})` is primarily determined by local atomic density and core-valence overlap, which is nearly identical between the quantum dot interior and the bulk crystal. By referencing the quantum dot self-energy to the known bulk quasiparticle correction :math:`\Delta \varepsilon_p^{\mathrm{bulk}} = Z^{\mathrm{bulk}} \langle \psi_p | \Sigma^{\mathrm{bulk}} - v_{xc}^{\mathrm{DFT}} | \psi_p \rangle`, the local potential :math:`v_{xc}^{\mathrm{DFT}}` cancels out identically:

.. math::

   \Delta \Sigma(\mathbf{r}, \mathbf{r}') = \Sigma^{\mathrm{QD}}(\mathbf{r}, \mathbf{r}') - \Sigma^{\mathrm{bulk}}(\mathbf{r}, \mathbf{r}') \approx \Delta \Sigma^{\mathrm{COH}}(\mathbf{r}, \mathbf{r}') + \Delta \Sigma^{\mathrm{SEX}}(\mathbf{r}, \mathbf{r}').


Screened COHSEX Operator
~~~~~~~~~~~~~~~~~~~~~~~~

In the static Coulomb-Hole plus Screened-Exchange (COHSEX) approximation:

.. math::

   \Delta W(\mathbf{r}, \mathbf{r}') = W^{\mathrm{QD}}(\mathbf{r}, \mathbf{r}') - W^{\mathrm{bulk}}(\mathbf{r}, \mathbf{r}') + W^{\mathrm{solv}}(\mathbf{r}, \mathbf{r}')

where:

1. **Screened Exchange** (:math:`\Delta \Sigma^{\mathrm{SEX}}`):

   .. math::

      \Delta \Sigma^{\mathrm{SEX}}(\mathbf{r}, \mathbf{r}') = -\rho(\mathbf{r}, \mathbf{r}') \, \Delta W(\mathbf{r}, \mathbf{r}')

2. **Coulomb Hole** (:math:`\Delta \Sigma^{\mathrm{COH}}`):

   .. math::

      \Delta \Sigma^{\mathrm{COH}}(\mathbf{r}, \mathbf{r}') = \frac{1}{2} \, \delta(\mathbf{r} - \mathbf{r}') \, \Delta W(\mathbf{r}, \mathbf{r}).

Projected onto state :math:`p`:

.. math::

   \Delta \varepsilon_p^{\mathrm{conf}} = \frac{Z_p}{2} \sum_{A=1}^{N_{\mathrm{atoms}}} q_A^p \, \Delta W_{AA} - Z_p \sum_{A,B} q_A^p \, \Delta W_{AB} \, q_B^{\mathrm{occ}}

where :math:`q_A^p = \sum_{\mu \in A, \nu} C_{\mu p} S_{\mu \nu} C_{\nu p}` is the Mulliken or Lowdin atomic orbital population of state :math:`p`.

The total quasiparticle scissor is:

.. math::

   \Delta_{\mathrm{sGW}} = \Delta_{\mathrm{bulk}} + \Delta \varepsilon_{\mathrm{LUMO}}^{\mathrm{conf}} - \Delta \varepsilon_{\mathrm{HOMO}}^{\mathrm{conf}}.
