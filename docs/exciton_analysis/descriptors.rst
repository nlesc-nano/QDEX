Descriptors
===========

Part of :doc:`/exciton_analysis/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.exciton_analysis``
* Callable: ``qdex.exciton_analysis.ExcitonAnalyzer``
* CLI: ``--nto, --nto-states``
* YAML: ``analysis.nto, analysis.nto_states``


2. Rigorous Plasser-Dreuw Spatial Descriptors
---------------------------------------------

From the normalized atomic hole populations :math:`q_A^h` and electron populations :math:`q_A^e`, ``QDEX`` computes eight rigorous physical metrics:


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

