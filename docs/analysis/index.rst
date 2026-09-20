Electronic Structure Analysis
=============================

`miniBSE` includes advanced analysis tools for characterizing excited-state wavefunctions, chemical bonding, and electronic band structures.

1. Exciton Spatial Descriptors
------------------------------

Based on the Plasser-Dreuw wavefunction analysis framework, `miniBSE` decomposes the two-particle exciton wavefunction $|\Psi_S\rangle$ into quantitative spatial metrics:

* **Electron-Hole Separation Distance ($d_{eh}$)**:
  The root-mean-square distance between the electron and the hole:

  .. math::

     d_{eh} = \sqrt{ \langle (\mathbf{r}_e - \mathbf{r}_h)^2 \rangle_S }

* **Pearson Correlation Coefficient ($R_{eh}$)**:
  Measures whether the electron and hole move together coherently ($R > 0$, bound Wannier-Mott exciton) or independently ($R \approx 0$ or $R < 0$, charge-transfer state).

* **Charge-Transfer Distance ($d_{\mathrm{CT}}$)**:
  The distance between the centroid of the hole density and the centroid of the electron density:

  .. math::

     d_{\mathrm{CT}} = |\langle \mathbf{r}_e \rangle_S - \langle \mathbf{r}_h \rangle_S|

* **Volumetric Transition Densities (.cube)**:
  Export 3D volumetric `.cube` files of the hole density $\rho_h(\mathbf{r})$ and electron density $\rho_e(\mathbf{r})$ for visualization in VMD, PyMOL, or ChimeraX:

  .. code-block:: bash

     minibse --config config.yaml --cube

2. Projected DOS & COOP Analysis
--------------------------------

`miniBSE` computes atom- and orbital-projected densities of states (PDOS) and Crystal Orbital Overlap Population (COOP) curves:

* **PDOS**: Decomposes the electronic density of states into specific atomic species (e.g. Pb $6s$, Br $4p$, Cs $5p$).
* **COOP**: Quantifies whether orbital interactions between specific pairs of atoms (e.g. Pb–Br) are bonding (COOP > 0) or antibonding (COOP < 0) as a function of energy.

Run via CLI:

.. code-block:: bash

   minibse --pdos --mo_file MOs.mbse --xyz structure.xyz --basis_txt BASIS_MOLOPT --basis_name DZVP-MOLOPT-PBE-GTH

3. Fuzzy Band Structure (Supercell Unfolding)
---------------------------------------------

For large nanocrystals or disordered clusters that lack periodic boundary conditions, `miniBSE` computes **Fuzzy Bands** by projecting cluster molecular orbitals onto plane-wave basis states:

.. math::

   P_{\mathbf{k}}(\varepsilon) = \sum_m |\langle e^{i \mathbf{k} \cdot \mathbf{r}} | \phi_m \rangle|^2 \, \delta(\varepsilon - \varepsilon_m)

This unfolds the dense, discrete cluster spectrum into an effective bulk-like dispersion $E(\mathbf{k})$ along high-symmetry paths in the Brillouin zone.
