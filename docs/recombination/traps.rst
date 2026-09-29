Traps
=====

Part of :doc:`/recombination/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.auger``
* Callable: ``qdex.auger.calculate_auger_rates``
* CLI: ``--auger, --auger-channel``
* YAML: ``auger.run, auger.channel``

.. code-block:: python

   calculate_auger_rates(C: np.ndarray, eps: np.ndarray, S: np.ndarray, atom_ao_ranges: List[Tuple[int, int]], coords: np.ndarray, atom_symbols: List[str], homo_idx: int, W_resta: Optional[np.ndarray]=None, material_name: Optional[str]='DEFAULT', eps_out: float=2.0, sigma_ev: float=0.05, broadening_mode: str='gaussian', channel: str='all', n_initial_elec: int=1, n_initial_hole: int=1, e_search_sigma_factor: float=4.0, spinor: bool=False, U_spinor_alpha: Optional[np.ndarray]=None, U_spinor_beta: Optional[np.ndarray]=None, lambda_reorg_ev: Optional[float]=None, temperature_k: float=300.0, eps_eff: Optional[float]=None, verbose: bool=True)


7. Defect Trap-Assisted Recombination (Shockley-Read-Hall)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In real quantum dots with unpassivated surfaces or vacancies, non-radiative recombination is overwhelmingly accelerated by **Shockley-Read-Hall (SRH) mid-gap traps**. Instead of bridging a single large gap of :math:`2.0\text{ eV}`, carriers drop into an intermediate trap state (:math:`\Delta E \approx 0.5 - 1.0\text{ eV}`), where multi-phonon tunneling is orders of magnitude faster.

In ``QDEX``, trap-assisted recombination can be configured via the YAML input:

.. code-block:: yaml

   namd:
     recombination:
       include_ground_state: true
       tau_nr_ns: 25.0  # Trap non-radiative lifetime in nanoseconds
