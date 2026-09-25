Auger workflows
===============

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

.. rubric:: From ``docs/part8_auger/index.rst:432-434``

7. Static Calculations vs. Trajectory Dynamics in NAMD
------------------------------------------------------


.. rubric:: From ``docs/part8_auger/index.rst:435-448``

Static Calculation
~~~~~~~~~~~~~~~~~~

To compute Auger recombination rates for a single quantum dot frame:

.. code-block:: bash

   qdex --config test_auger.yaml

Sample formatted output in nanoseconds:

.. code-block:: text

   =======================================================================================================

                                      QDEX AUGER RECOMBINATION REPORT                                     
   =======================================================================================================
     Fundamental Bandgap (E_g)       :   3.3050 eV
     Energy Conservation Line Shape  : GAUSSIAN (sigma = 50.0 meV)

     Evaluated Active Pathways       : 177 (eeh) | 48 (hhe)
   -------------------------------------------------------------------------------------------------------

     Channel                 Rate (s^-1)        Rate (ns^-1)        Lifetime (ns)          Lifetime (ps)  
   -------------------------------------------------------------------------------------------------------
     Negative Trion (eeh)       4.6613e+10        4.6613e+01            0.0215 ns           21.45 ps
     Positive Trion (hhe)       8.6182e+08        8.6182e-01            1.1603 ns         1160.34 ps

     Biexciton (XX)             1.8990e+11        1.8990e+02            0.0053 ns            5.27 ps
   -------------------------------------------------------------------------------------------------------


.. rubric:: From ``docs/part8_auger/index.rst:462-474``

Trajectory-Averaged Dynamic Auger Rates in NAMD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Along an *ab initio* molecular dynamics (AIMD) trajectory, thermal vibrations modulate orbital energies and continuum crossings. ``QDEX`` evaluates the trajectory average:

.. math::

   \langle \Gamma_{XX} \rangle = \frac{1}{N_{\mathrm{frames}}} \sum_{k=1}^{N_{\mathrm{frames}}} \Gamma_{XX}(t_k)

This trajectory average naturally samples the true vibronic density of states without depending on arbitrary Gaussian broadening widths.

---
