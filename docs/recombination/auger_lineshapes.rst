Auger lineshapes
================

Part of :doc:`/recombination/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.auger``
* Callable: ``qdex.auger.calculate_auger_rates``
* CLI: ``--auger, --auger-channel``
* YAML: ``auger.run, auger.channel``

.. code-block:: python

   calculate_auger_rates(C: np.ndarray, eps: np.ndarray, S: np.ndarray, atom_ao_ranges: List[Tuple[int, int]], coords: np.ndarray, atom_symbols: List[str], homo_idx: int, W_resta: Optional[np.ndarray]=None, material_name: Optional[str]='DEFAULT', eps_out: float=2.0, sigma_ev: float=0.05, broadening_mode: str='gaussian', channel: str='all', n_initial_elec: int=1, n_initial_hole: int=1, e_search_sigma_factor: float=4.0, spinor: bool=False, U_spinor_alpha: Optional[np.ndarray]=None, U_spinor_beta: Optional[np.ndarray]=None, lambda_reorg_ev: Optional[float]=None, temperature_k: float=300.0, eps_eff: Optional[float]=None, verbose: bool=True)


5. Energy Conservation Line Shapes: Gaussian vs. FCWD
-----------------------------------------------------

In Fermi's Golden Rule, energy conservation is satisfied when the ejected carrier lands in a virtual state located approximately one bandgap above the conduction edge:


1. Gaussian Line Shape Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Thermal broadening is represented by a Gaussian envelope:

.. math::

   \rho(\Delta E) = \frac{1}{\sqrt{2\pi}\sigma} \exp\left( -\frac{\Delta E^2}{2\sigma^2} \right)

where :math:`\sigma` (configured via ``--auger-sigma``, default :math:`0.05\text{ eV}`) matches the thermal energy scale :math:`2 k_B T` at 300 K.


2. Multi-Phonon Marcus / Jortner Line Shape (FCWD)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In nanocrystals, Auger recombination is phonon-assisted. When the electronic energy mismatch :math:`\Delta E` is non-zero, lattice phonons supply or absorb the discrepancy:

.. math::

   \operatorname{FCWD}(\Delta E) = \frac{1}{\sqrt{4\pi \lambda k_B T}} \exp\left( -\frac{(\Delta E - \lambda)^2}{4\lambda k_B T} \right)

where the nuclear reorganization energy :math:`\lambda` and optical phonon energy :math:`\hbar \omega_{\mathrm{LO}}` are derived directly from the **Phonon Spectral Density** :math:`J(\omega)` of the NAMD trajectory (see Part 6).

