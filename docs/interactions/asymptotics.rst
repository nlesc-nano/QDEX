Asymptotics
===========

Part of :doc:`/interactions/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.build_gamma``
* CLI: ``--2e-integrals, --kernel, --eps-out``
* YAML: ``physics.2e-integrals, physics.kernel, physics.eps_out``

.. code-block:: python

   build_gamma(atom_symbols, coords, alpha, beta=0.0, eta_dict=HARDNESS_DICT)


Spatial Asymptotics of the Dielectric Kernel: Why Screening Fits Nanocrystals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The interaction kernel :math:`W(\mathbf{r}, \mathbf{r}')` exhibits three essential spatial regimes that make it uniquely suited for semiconductor nanocrystals:

1. Short-Range Limit (:math:`r \to 0`, intra-atomic / on-site):
   At sub-nanometer distances, the valence electron gas cannot displace fast enough to screen charge fluctuations without violating quantum kinetic energy constraints (Pauli exclusion / Thomas-Fermi cutoff). Therefore, the screening function satisfies:

   .. math::

      \lim_{r \to 0} W(r) = v(r) = \frac{1}{r}, \quad (\epsilon \to 1).

   Standard macroscopic dielectric approximations (:math:`W = v / \epsilon_\infty`) artificially divide on-site Coulomb repulsion by :math:`\epsilon_\infty \approx 6 - 10`, severely underestimating atomic charging energies and exciton binding energies. The Resta and DIM kernels preserve the true unscreened Coulomb interaction at short range.

2. **Intermediate Confinement & Solvent Boundary** (:math:`r \sim R_{\mathrm{QD}}`):
   Within the nanocrystal core, quantum confinement suppresses electronic polarizability, reducing the effective internal permittivity :math:`\epsilon_{\mathrm{in}}(R) < \epsilon_\infty`. At the nanocrystal-solvent interface, the dielectric mismatch with the surrounding medium (:math:`\epsilon_{\mathrm{out}}`) induces image charges that screen carriers through the reaction field :math:`W^{\mathrm{solv}}`.

3. Asymptotic Bulk Limit (:math:`r \to \infty` or :math:`R_{\mathrm{QD}} \to \infty`):
   Across large distances inside the crystal, dielectric polarization fully develops:

   .. math::

      \lim_{r \to \infty} W(r) = \frac{1}{\epsilon_\infty r}.

   Crucially, as the nanocrystal size grows (:math:`R_{\mathrm{QD}} \to \infty`), the cluster screening profile converges to the periodic crystal limit:

   .. math::

      \lim_{R_{\mathrm{QD}} \to \infty} W^{\mathrm{QD}}(\mathbf{r}, \mathbf{r}') = W^{\mathrm{bulk}}(\mathbf{r}, \mathbf{r}').

   Consequently, the net confinement shift :math:`\Delta W = W^{\mathrm{QD}} - W^{\mathrm{bulk}} + W^{\mathrm{solv}} \to 0` vanishes identically for large crystals, naturally recovering the bulk quasiparticle band gap and bulk exciton spectrum.
