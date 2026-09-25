Auger configuration
===================

Part of :doc:`/recombination/index`.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.auger``
* Callable: ``qdex.auger.calculate_auger_rates``
* CLI: ``--auger, --auger-channel``
* YAML: ``auger.run, auger.channel``

.. code-block:: python

   calculate_auger_rates(C: np.ndarray, eps: np.ndarray, S: np.ndarray, atom_ao_ranges: List[Tuple[int, int]], coords: np.ndarray, atom_symbols: List[str], homo_idx: int, W_resta: Optional[np.ndarray]=None, material_name: Optional[str]='DEFAULT', eps_out: float=2.0, sigma_ev: float=0.05, broadening_mode: str='gaussian', channel: str='all', n_initial_elec: int=1, n_initial_hole: int=1, e_search_sigma_factor: float=4.0, spinor: bool=False, U_spinor_alpha: Optional[np.ndarray]=None, U_spinor_beta: Optional[np.ndarray]=None, lambda_reorg_ev: Optional[float]=None, temperature_k: float=300.0, eps_eff: Optional[float]=None, verbose: bool=True)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part8_auger/index.rst:583-585``

9. Configuration Reference (YAML & CLI)
---------------------------------------


.. rubric:: From ``docs/part8_auger/index.rst:586-607``

YAML Configuration Options
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   auger:
     run: true               # Enable static Auger recombination calculation
     sigma: 0.05             # Energy conservation broadening width in eV
     channel: "all"          # "all" (both eeh and hhe), "eeh", or "hhe"
     n_initial_states: 1     # Number of frontier band-edge carriers to consider
     lineshape: "gaussian"   # "gaussian" or "fcwd"
     # eps_eff is ignored. Resta already screens as eps(r). Do not stack a second factor.

   namd:
     dynamics:
       method: "cpa_fssh"
       ecsh_auger: true        # Enable Energy-Conserving Surface Hopping for Auger
       ecsh_window_ev: 0.026   # Resonance window in eV for Auger transitions (default: k_B*T)
       initial_state: "biexciton" # "biexciton" (XX -> X) or "ratio_eg" (single exciton)
       tau_auger_ps: 50.0      # Biexciton Auger lifetime in ps (or estimated from Resta)
       trajectory_loops: 10    # Number of times to loop MD trajectory to reach long timescales


.. rubric:: From ``docs/part8_auger/index.rst:608-623``

Command-Line Arguments
~~~~~~~~~~~~~~~~~~~~~~

Static Auger Calculation:
* ``--auger``: Enable Auger recombination rate calculation.
* ``--auger-sigma <float>``: Energy conservation broadening width in eV (default: ``0.05``).
* ``--auger-channel {all,eeh,hhe}``: Recombination channel to compute (default: ``all``).
* ``--auger-states <int>``: Number of band-edge frontier states to consider as initial carriers (default: ``1``).
* ``--auger-lineshape {gaussian,fcwd}``: Energy conservation line shape (default: ``gaussian``).
* ``--auger-eps-eff <float>``: Effective dielectric constant for dynamic screening at :math:`\hbar\omega = E_g`.

Dynamic ECSH Auger Simulation:
* ``--namd-ecsh-auger``: Enable Energy-Conserving Surface Hopping (ECSH) for Auger in NAMD.
* ``--namd-ecsh-window <float>``: Resonance energy window in eV for ECSH Auger transitions (default: :math:`k_B T`).
* ``--namd-trajectory-loops <int>``: Number of times to loop precomputed MD trajectory to reach long Auger timescales.
* ``--namd-biexciton``: Initialize NAMD from a biexciton state (:math:`XX`) to simulate Auger annihilation dynamics.
