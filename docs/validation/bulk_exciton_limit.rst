Bulk exciton limit
==================

Part of :doc:`/validation/index`.

.. important::

   The Wannier–Mott limit requires converged coupled BSE transitions, effective masses and a macroscopic screened tail. It is not guaranteed by a pair potential alone, and diagonal BSE does not generally recover finite bulk binding.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.solver``
* Callable: ``qdex.solver.solve``
* CLI: ``--nroots, --tol``
* YAML: ``bse.nroots, bse.tol``

.. code-block:: python

   solve(self, nroots=10, full_diag=False, tol=1e-05, excitation_mode='bse')

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part4_excited_states/index.rst:291-313``

Spatial Asymptotics & Wannier-Mott Bulk Limit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As established in the quasiparticle theory (:doc:`../part3_gw_scissor/index`), the screened interaction :math:`W(r)` connects short-range atomic scales to the macroscopic crystal:

1. **Short-Range Limit** (:math:`r \to 0`): :math:`W(r) \to v(r) = 1/r` (:math:`\epsilon \to 1`). Electronic screening ceases at sub-nanometer distances because the valence electrons cannot instantaneously compress inside an atomic core. This prevents the catastrophic underestimation of singlet-triplet exchange splitting and on-site Coulomb repulsion.
2. **Nanocrystal Boundary** (:math:`r \sim R_{\mathrm{QD}}`): Dielectric mismatch between the dot (:math:`\epsilon_\infty`) and the solvent (:math:`\epsilon_{\mathrm{out}}`) generates an image-charge reaction field :math:`W^{\mathrm{solv}}`. In a consistent charged/neutral treatment this boundary can change optical energies; QDEX currently does not include ``eps_out`` in its Resta/DIM BSE direct kernel.
3. Asymptotic Bulk Limit (:math:`r \to \infty`, :math:`R_{\mathrm{QD}} \to \infty`): :math:`W(r) \to \frac{1}{\epsilon_\infty r}`. In this limit, the Bethe-Salpeter equation continuously reduces to the hydrogenic Wannier-Mott exciton equation:

   .. math::

      \left( -\frac{\hbar^2 \nabla_{\mathbf{r}}^2}{2\mu} - \frac{e^2}{\epsilon_\infty r} \right) \phi_{\mathrm{exc}}(\mathbf{r}) = -E_b^{\mathrm{bulk}} \phi_{\mathrm{exc}}(\mathbf{r})

   recovering the bulk Rydberg binding energy:

   .. math::

      E_b^{\mathrm{bulk}} = \frac{\mu e^4}{2 \hbar^2 \epsilon_\infty^2} = R_y^*

   with bulk exciton Bohr radius :math:`a_{\mathrm{exc}} = a_0 \epsilon_\infty (m_0 / \mu)`. Consequently, as the nanocrystal diameter surpasses the Bohr radius (:math:`R_{\mathrm{QD}} \gg a_{\mathrm{exc}}`), a converged bulk BSE optical transition should approach the bulk band edge minus the Wannier-Mott binding energy: :math:`\Omega_1 \to E_g^{\mathrm{bulk}} - E_b^{\mathrm{bulk}}`.

---
