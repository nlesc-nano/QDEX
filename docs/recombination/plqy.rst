Plqy
====

Part of :doc:`/recombination/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.auger``
* Callable: ``qdex.auger.calculate_auger_rates``
* CLI: ``--auger, --auger-channel``
* YAML: ``auger.run, auger.channel``

.. code-block:: python

   calculate_auger_rates(C: np.ndarray, eps: np.ndarray, S: np.ndarray, atom_ao_ranges: List[Tuple[int, int]], coords: np.ndarray, atom_symbols: List[str], homo_idx: int, W_resta: Optional[np.ndarray]=None, material_name: Optional[str]='DEFAULT', eps_out: float=2.0, sigma_ev: float=0.05, broadening_mode: str='gaussian', channel: str='all', n_initial_elec: int=1, n_initial_hole: int=1, e_search_sigma_factor: float=4.0, spinor: bool=False, U_spinor_alpha: Optional[np.ndarray]=None, U_spinor_beta: Optional[np.ndarray]=None, lambda_reorg_ev: Optional[float]=None, temperature_k: float=300.0, eps_eff: Optional[float]=None, verbose: bool=True)


8. Photoluminescence Quantum Yield (PLQY)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The total Photoluminescence Quantum Yield (PLQY) represents the branching ratio between radiative photon emission and non-radiative dissipation:

.. math::

   \mathrm{PLQY} = \frac{\langle k_{\mathrm{rad}} \rangle}{\langle k_{\mathrm{rad}} \rangle + k_{\mathrm{nr}}} \times 100\%

* In pristine, defect-free nanocrystals where :math:`E_g \gg \hbar \omega_{\mathrm{LO}}`, the intrinsic non-radiative rate is negligible (:math:`k_{\mathrm{nr}} \ll 10^3\text{ s}^{-1}`), leading to near-unity PLQY (:math:`\sim 99\%`).
* In the presence of surface traps (:math:`\tau_{\mathrm{nr}} \sim 10 - 50\text{ ns}`), non-radiative decay competes directly with radiative emission (:math:`\tau_{\mathrm{rad}} \sim 5 - 20\text{ ns}`), yielding realistic PLQYs between :math:`20\%` and :math:`70\%`.

