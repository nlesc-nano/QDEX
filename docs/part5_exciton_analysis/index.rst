Part 5: Excited-State Wavefunction Analysis (Plasser-Dreuw)
===========================================================

Solving the Bethe-Salpeter Equation yields excitation energies :math:`\Omega_S` and eigenvectors :math:`\mathbf{X}_S = (X_{ia}^S)`. However, understanding the physical nature of an exciton—whether it is a tightly bound Wannier-Mott exciton, a localized Frenkel exciton, a surface-trap state, or a spatial charge-transfer (CT) excitation—requires quantitative real-space wavefunction analysis.

``miniBSE`` incorporates the rigorous **Plasser-Dreuw exciton analysis framework** along with **Natural Transition Orbitals (NTOs)** to decompose complex multi-configurational exciton wavefunctions into intuitive, publication-ready physical descriptors.

---

1. The Two-Particle Transition Density Matrix
---------------------------------------------

The two-particle exciton wavefunction is expanded in electron-hole configurations:

.. math::

   |\Psi_S\rangle = \sum_{i \in \mathrm{occ}} \sum_{a \in \mathrm{virt}} X_{ia}^S \, a_a^\dagger a_i |0\rangle

In real space, the transition density :math:`\gamma_0^S(\mathbf{r}_h, \mathbf{r}_e)` describes the joint probability amplitude of finding the hole at :math:`\mathbf{r}_h` and the electron at :math:`\mathbf{r}_e`:

.. math::

   \gamma_0^S(\mathbf{r}_h, \mathbf{r}_e) = \sum_{i \in \mathrm{occ}} \sum_{a \in \mathrm{virt}} X_{ia}^S \, \phi_i(\mathbf{r}_h) \, \phi_a(\mathbf{r}_e)

Reduced Hole and Electron Densities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Integrating out one of the quasi-particles yields the reduced spatial densities of the hole (:math:`\rho_h`) and electron (:math:`\rho_e`):

.. math::

   \rho_h(\mathbf{r}) = \int |\gamma_0^S(\mathbf{r}, \mathbf{r}_e)|^2 \, d\mathbf{r}_e = \sum_{i, j \in \mathrm{occ}} D_{ij}^h \, \phi_i(\mathbf{r}) \phi_j^*(\mathbf{r})

.. math::

   \rho_e(\mathbf{r}) = \int |\gamma_0^S(\mathbf{r}_h, \mathbf{r})|^2 \, d\mathbf{r}_h = \sum_{a, b \in \mathrm{virt}} D_{ab}^e \, \phi_a^*(\mathbf{r}) \phi_b(\mathbf{r})

where the single-particle reduced density matrices are:

.. math::

   \mathbf{D}^h = \mathbf{X} \mathbf{X}^\dagger, \quad \mathbf{D}^e = \mathbf{X}^\dagger \mathbf{X}

Vectorized Low-Memory Mulliken Population
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In standard implementations, computing :math:`\rho_h` and :math:`\rho_e` requires forming full :math:`N_{\mathrm{ao}} \times N_{\mathrm{ao}}` AO density matrices. In ``miniBSE``, atomic populations are evaluated directly in the active molecular orbital subspace:

.. math::

   q_A^h = \sum_{\mu \in A} \operatorname{Re}\left[ \left( \mathbf{C}_{\mathrm{occ}} \mathbf{D}^h \right)_{\mu, :} (\mathbf{S} \mathbf{C}_{\mathrm{occ}})_{\mu, :}^* \right]

.. math::

   q_A^e = \sum_{\mu \in A} \operatorname{Re}\left[ \left( \mathbf{C}_{\mathrm{virt}} \mathbf{D}^e \right)_{\mu, :} (\mathbf{S} \mathbf{C}_{\mathrm{virt}})_{\mu, :}^* \right]

This vectorized reduction avoids large intermediate arrays and runs in milliseconds even for large nanoclusters.

---

2. Rigorous Plasser-Dreuw Spatial Descriptors
---------------------------------------------

From the normalized atomic hole populations :math:`q_A^h` and electron populations :math:`q_A^e`, ``miniBSE`` computes eight rigorous physical metrics:

1. Spatial Centroids (:math:`\langle \mathbf{r}_h \rangle, \langle \mathbf{r}_e \rangle`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The center-of-mass positions of the hole and electron probability distributions:

.. math::

   \langle \mathbf{r}_h \rangle = \sum_{A=1}^{N_{\mathrm{atoms}}} q_A^h \mathbf{R}_A, \quad
   \langle \mathbf{r}_e \rangle = \sum_{A=1}^{N_{\mathrm{atoms}}} q_A^e \mathbf{R}_A

2. Charge-Transfer Distance (:math:`d_{\mathrm{CT}}`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The absolute vector distance between the electron centroid and the hole centroid:

.. math::

   d_{\mathrm{CT}} = \left| \langle \mathbf{r}_e \rangle - \langle \mathbf{r}_h \rangle \right|

* :math:`d_{\mathrm{CT}} \approx 0\text{ Å}`: Localized or spherically symmetric excitation.
* :math:`d_{\mathrm{CT}} > 3\text{ Å}`: Macroscopic dipolar charge separation (e.g. core-to-surface or type-II heterojunction state).

3. Root-Mean-Square Particle Sizes (:math:`\sigma_h, \sigma_e`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The spatial spread (standard deviation) of the hole and electron clouds around their respective centers:

.. math::

   \sigma_h = \sqrt{ \sum_{A=1}^{N_{\mathrm{atoms}}} q_A^h \left| \mathbf{R}_A - \langle \mathbf{r}_h \rangle \right|^2 }

.. math::

   \sigma_e = \sqrt{ \sum_{A=1}^{N_{\mathrm{atoms}}} q_A^e \left| \mathbf{R}_A - \langle \mathbf{r}_e \rangle \right|^2 }

In quantum dots, comparing :math:`\sigma_e` with :math:`\sigma_h` reveals whether the electron is more delocalized than the heavier hole, identifying quantum confinement asymmetry.

4. Electron-Hole Spatial Covariance (:math:`\mathrm{Cov}_{eh}`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Measures the extent to which the electron and hole coordinate fluctuations move together:

.. math::

   \mathrm{Cov}_{eh} = \sum_{i \in \mathrm{occ}} \sum_{a \in \mathrm{virt}} |X_{ia}^S|^2 \left( \mathbf{r}_a - \langle \mathbf{r}_e \rangle \right) \cdot \left( \mathbf{r}_i - \langle \mathbf{r}_h \rangle \right)

where :math:`\mathbf{r}_i = \sum_A P_{A, i} \mathbf{R}_A` is the spatial center of orbital :math:`\phi_i`.

5. Pearson Correlation Coefficient (:math:`R_{eh}`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Normalizes the covariance by the individual particle spreads:

.. math::

   R_{eh} = \frac{\mathrm{Cov}_{eh}}{\sigma_h \, \sigma_e}

* :math:`R_{eh} > 0` (typically :math:`+0.3` to :math:`+0.9`): **Bound Exciton**. The electron and hole are spatially correlated and track each other's motion (Wannier-Mott exciton).
* :math:`R_{eh} \approx 0`: **Uncorrelated Transitions**. The electron and hole distributions are independent.
* :math:`R_{eh} < 0`: **Anti-correlated / Charge Transfer**. When the hole localizes on one facet, the electron is displaced to the opposing facet.

6. True Root-Mean-Square Exciton Size (:math:`d_{eh}`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The true RMS separation between the electron and hole coordinates:

.. math::

   d_{eh} = \sqrt{ \langle |\mathbf{r}_e - \mathbf{r}_h|^2 \rangle } = \sqrt{ \sigma_h^2 + \sigma_e^2 + d_{\mathrm{CT}}^2 - 2 \mathrm{Cov}_{eh} }

Note that :math:`d_{eh} \neq d_{\mathrm{CT}}`. For a centrosymmetric Wannier exciton, :math:`d_{\mathrm{CT}} = 0`, but :math:`d_{eh}` correctly reports the finite Bohr radius of the electron-hole pair.

7. Charge-Transfer Ratio (:math:`\mathrm{CT}`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The fraction of the exciton size attributable to net dipolar displacement:

.. math::

   \mathrm{CT} = \frac{d_{\mathrm{CT}}}{d_{eh}} \in [0, 1]

8. Exciton Participation Ratio (:math:`\mathrm{PR}`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Quantifies configuration interaction delocalization across electron-hole basis pairs:

.. math::

   \mathrm{PR} = \frac{1}{\sum_{ia} |X_{ia}^S|^4}

* :math:`\mathrm{PR} = 1.0`: Pure single-configuration transition (:math:`|i \to a\rangle`).
* :math:`\mathrm{PR} \gg 1`: Strong multi-configurational mixing, characteristic of delocalized Wannier excitons in semiconductor nanocrystals.

---

3. Natural Transition Orbitals (NTOs)
-------------------------------------

In large systems, the BSE expansion vector :math:`X_{ia}^S` often contains hundreds of small non-zero coefficients, obscuring the primary orbital character of the excitation.

**Natural Transition Orbitals (NTOs)** solve this problem by applying a **Singular Value Decomposition (SVD)** to the transition density matrix :math:`\mathbf{X}^S`:

.. math::

   \mathbf{X}^S = \mathbf{U} \, \boldsymbol{\Lambda} \, \mathbf{V}^\dagger = \sum_{k=1}^{K_{\max}} \lambda_k \, \mathbf{u}_k \, \mathbf{v}_k^\dagger

where:
* :math:`\lambda_k \ge 0` are the singular values, normalized such that :math:`\sum_k \lambda_k^2 = 1`.
* :math:`\mathbf{u}_k` defines the :math:`k`-th **NTO hole orbital**: :math:`\tilde{\phi}_k^h(\mathbf{r}) = \sum_i U_{ik} \phi_i(\mathbf{r})`.
* :math:`\mathbf{v}_k` defines the :math:`k`-th **NTO electron orbital**: :math:`\tilde{\phi}_k^e(\mathbf{r}) = \sum_a V_{ak} \phi_a(\mathbf{r})`.

The complex BSE excitation is transformed into a compact sum of uncoupled transitions:

.. math::

   |\Psi_S\rangle = \sum_{k} \lambda_k \, |\tilde{\phi}_k^h \to \tilde{\phi}_k^e\rangle

In most semiconductor quantum dots, **a single NTO pair accounts for :math:`> 85\%` of the total transition weight** (:math:`\lambda_1^2 > 0.85`), providing an immediate visual interpretation of the excitation.

NTO Compactness Metrics
~~~~~~~~~~~~~~~~~~~~~~~

``miniBSE`` outputs comprehensive NTO statistics:
* **Lead Weight (:math:`\lambda_1^2`)**: Fraction of the transition described by the dominant NTO pair.
* **NTO Participation Ratio**: :math:`\mathrm{PR}_{\mathrm{NTO}} = \frac{1}{\sum_k \lambda_k^4}`.
* **NTO Shannon Entropy**: :math:`S_{\mathrm{NTO}} = -\sum_k \lambda_k^2 \ln(\lambda_k^2)`.
* **:math:`N_{90}` and :math:`N_{99}`**: Number of NTO pairs required to recover 90% and 99% of the transition weight.

---

4. 3D Volumetric Visualization (.cube Files)
--------------------------------------------

``miniBSE`` generates 3D volumetric Gaussian ``.cube`` files for direct visualization in VMD, PyMOL, or ChimeraX:

* **Hole Density**: ``exciton_S1_hole.cube``: :math:`\rho_h(\mathbf{r}) = \sum_{ia} |X_{ia}^1|^2 |\phi_i(\mathbf{r})|^2`
* **Electron Density**: ``exciton_S1_elec.cube``: :math:`\rho_e(\mathbf{r}) = \sum_{ia} |X_{ia}^1|^2 |\phi_a(\mathbf{r})|^2`
* **Difference Density**: ``exciton_S1_diff.cube``: :math:`\Delta \rho(\mathbf{r}) = \rho_e(\mathbf{r}) - \rho_h(\mathbf{r})`

Positive values in the difference cube indicate regions of net electron accumulation (photo-induced negative charge), while negative values indicate regions of hole accumulation (photo-induced positive charge).

---

5. Interactive Plotly 6-Panel Dashboard
---------------------------------------

When running with ``--plot``, ``miniBSE`` exports a self-contained, interactive HTML dashboard (``exciton_analysis.html``) featuring six coordinated subplots:

1. **Participation Ratio (PR)** vs. Energy (identifying multiconfigurational states).
2. **True Exciton Size (:math:`d_{eh}`)** vs. Energy (characterizing Bohr radius scaling).
3. **Particle Spread (:math:`\sigma_h, \sigma_e`)** (comparing electron vs. hole delocalization).
4. **Charge-Transfer Distance (:math:`d_{\mathrm{CT}}`)** color-coded by CT ratio.
5. **Spatial Pearson Correlation (:math:`R_{eh}`)** (distinguishing bound from dissociated excitons).
6. **Simulated UV-Vis Absorption Spectrum** with oscillator strength stick spectra and Gaussian broadening.

---

6. CLI Flags & YAML Configuration Reference
-------------------------------------------

Command-Line Arguments
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - CLI Flag
     - Default
     - Description
   * - ``--nto``
     - ``False``
     - Run Natural Transition Orbital (NTO) analysis after solving BSE.
   * - ``--nto-states <list>``
     - ``1 2 3``
     - Specific exciton states to analyze with NTOs (1-indexed, e.g. ``--nto-states 1 5 10``).
   * - ``--nto-top <int>``
     - ``3``
     - Number of dominant NTO pairs to report per state.
   * - ``--nto-csv``
     - ``False``
     - Export detailed NTO descriptors and weights to ``nto_results.csv``.
   * - ``--bse_states <list>``
     - ``1 2 3``
     - Specific exciton roots to export as 3D volumetric ``.cube`` files.
   * - ``--nbse <int>``
     - ``3``
     - Number of lowest exciton states to export as ``.cube`` files if ``--bse_states`` is omitted.
   * - ``--plot``
     - ``False``
     - Generate publication-ready figures and interactive Plotly HTML dashboards.
   * - ``--show``
     - ``False``
     - Display interactive plots in the web browser upon calculation completion.

YAML Configuration Example
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   analysis:
     nto: true
     nto_states: [1, 2, 3]
     nto_top: 3
     nto_csv: true
     plot: true

   cube:
     export: true
     bse_states: [1, 2]
     spacing_ang: 0.5
