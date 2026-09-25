Trajectory parameters
=====================

Part of :doc:`/recombination/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.auger``
* Callable: ``qdex.auger.calculate_auger_rates``
* CLI: ``--auger, --auger-channel``
* YAML: ``auger.run, auger.channel``

.. code-block:: python

   calculate_auger_rates(C: np.ndarray, eps: np.ndarray, S: np.ndarray, atom_ao_ranges: List[Tuple[int, int]], coords: np.ndarray, atom_symbols: List[str], homo_idx: int, W_resta: Optional[np.ndarray]=None, material_name: Optional[str]='DEFAULT', eps_out: float=2.0, sigma_ev: float=0.05, broadening_mode: str='gaussian', channel: str='all', n_initial_elec: int=1, n_initial_hole: int=1, e_search_sigma_factor: float=4.0, spinor: bool=False, U_spinor_alpha: Optional[np.ndarray]=None, U_spinor_beta: Optional[np.ndarray]=None, lambda_reorg_ev: Optional[float]=None, temperature_k: float=300.0, eps_eff: Optional[float]=None, verbose: bool=True)


4. Non-Empirical Extraction of Optical Phonon Energy from NAMD Spectral Density
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rather than relying on empirical phonon frequencies, ``QDEX`` extracts :math:`\hbar \omega_{\mathrm{LO}}` directly from the **Phonon Spectral Density** :math:`J(\omega)` of the NAMD trajectory:

1. The instantaneous energy gap fluctuation of the lowest transition along the MD trajectory is tracked:

   .. math::

      \delta E_g(t) = E_g(t) - \langle E_g \rangle

2. The normalized gap time-autocorrelation function is evaluated:

   .. math::

      C(t) = \frac{\langle \delta E_g(0) \delta E_g(t) \rangle}{\sigma_E^2}

3. A Hann-windowed Fourier transform stores :math:`\max(\mathrm{Re}\,\mathrm{FFT}[C], 0)` as :math:`J(\omega)` in wavenumbers (:math:`\text{cm}^{-1}`). The spacing of an independent frequency bin is :math:`1/T`:

   .. math::

      J(\omega) = \frac{1}{2\pi} \int_{-\infty}^{\infty} C(t) W(t) \, e^{i \omega t} \, dt

4. The dominant optical phonon mode is identified from the primary peak of :math:`J(\omega)` (excluding low-frequency acoustic noise :math:`< 30\text{ cm}^{-1}`):

   .. math::

      \tilde{\nu}_{\mathrm{LO}} = \operatorname{argmax}_{\tilde{\nu} \ge 30\text{ cm}^{-1}} J(\tilde{\nu})

5. Converting from wavenumber to energy gives the optical phonon quantum:

   .. math::

      \hbar \omega_{\mathrm{LO}} = h c \, \tilde{\nu}_{\mathrm{LO}} = (1.23984 \times 10^{-4}\text{ eV}\cdot\text{cm}) \times \tilde{\nu}_{\mathrm{LO}}

   For example, in lead halide perovskites (:math:`\text{CsPbBr}_3`), the dominant peak at :math:`\tilde{\nu} \approx 150\text{ cm}^{-1}` yields :math:`\hbar \omega_{\mathrm{LO}} = 18.6\text{ meV}`. In CdSe nanocrystals (:math:`\tilde{\nu} \approx 210\text{ cm}^{-1}`), it yields :math:`\hbar \omega_{\mathrm{LO}} = 26.0\text{ meV}`.


5. Derivation of Huang-Rhys Factor S and Reorganization Energy λ from Trajectory Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A central parameter in multi-phonon transitions is the **Huang-Rhys factor** :math:`S`, which quantifies the average number of phonons emitted during electronic transition. In ``QDEX``, :math:`S` and the nuclear reorganization energy :math:`\lambda` are determined non-empirically via the **Fluctuation-Dissipation Theorem / Marcus linear response theory**:

1. **Thermal Gap Variance**:
   Along the ab initio trajectory at temperature :math:`T`, the classical variance of the energy gap is computed directly:

   .. math::

      \sigma_E^2 = \langle (E_g(t) - \langle E_g \rangle)^2 \rangle = \frac{1}{N_{\mathrm{frames}}} \sum_{k=1}^{N_{\mathrm{frames}}} \delta E_g(t_k)^2

2. **Nuclear Reorganization Energy** (:math:`\lambda`):
   In linear response theory for a harmonic bath in the classical limit (:math:`k_B T \gg \hbar \omega / 2`), the energy gap variance is directly proportional to the reorganization energy:

   .. math::

      \sigma_E^2 = 2 \lambda \, k_B T \implies \lambda = \frac{\sigma_E^2}{2 \, k_B T}

   (Quantum mechanically, this corresponds to :math:`\sigma_E^2 = \int_0^\infty \frac{2}{\pi} \hbar \omega J(\omega) \coth\left(\frac{\hbar \omega}{2 k_B T}\right) d\omega`).

3. **Huang-Rhys Factor** (:math:`S`):
   Because the total reorganization energy partitioned into the dominant optical phonon mode of frequency :math:`\hbar \omega_{\mathrm{LO}}` is :math:`\lambda = S \, \hbar \omega_{\mathrm{LO}}`, we solve directly for :math:`S`:

   .. math::

      S = \frac{\lambda}{\hbar \omega_{\mathrm{LO}}} = \frac{\sigma_E^2}{2 \, k_B T \, \hbar \omega_{\mathrm{LO}}}

**Physical Significance**: Both :math:`\lambda` and :math:`S` are extracted directly from the NAMD trajectory without any adjustable parameters. Soft, polar perovskite lattices with large thermal gap fluctuations (:math:`\sigma_E \approx 60\text{ meV}`) yield :math:`\lambda \approx 70\text{ meV}` and :math:`S \approx 3.8`, reflecting significant electron-phonon coupling, whereas rigid covalent nanocrystals exhibit :math:`S \approx 0.5 - 1.5`.
