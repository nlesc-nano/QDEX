Auger implementation
====================

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

.. rubric:: From ``docs/part8_auger/index.rst:336-383``

4. How QDEX Data Are Used to Compute Every Formula Term
-------------------------------------------------------

Below is the exact mapping showing how each term in the mathematical formulas is computed from the underlying ``QDEX`` data structures:

.. list-table::
   :widths: 25 35 40
   :header-rows: 1

   * - Mathematical Term
     - QDEX Data Structure / Source
     - Computational Implementation
   * - **Molecular Orbitals** :math:`C_{\mu p}`
     - ``C`` from ``read_mos_mbse(mo_file, n_ao)``
     - :math:`(N_{\mathrm{ao}}, N_{\mathrm{mo}})` coefficient matrix from CP2K DFT.
   * - **AO Overlap** :math:`S_{\mu \nu}`
     - ``S = libint_cpp.overlap(shells, nthreads)``
     - Analytically computed using Libint2 primitive Gaussian overlaps.
   * - **Atomic AO Ranges**
     - ``atom_ao_ranges = build_atom_ao_ranges(shells)``
     - List of :math:`( \text{start\_ao}, \text{stop\_ao} )` for each atom :math:`A`.
   * - **Transition Charges** :math:`q_A^{ij}`
     - ``q_mat[A, :] = np.sum(C[ao_slice, i] * (S[ao_slice, :] @ C[:, j]), axis=0)``
     - Evaluated via vectorized matrix multiplications over atomic blocks :math:`A`.
   * - **Resta Kernel** :math:`W_{AB}`
     - ``W_resta = build_resta_mnok(syms, coords, ...)``
     - :math:`(N_{\mathrm{atoms}}, N_{\mathrm{atoms}})` matrix combining distance, hardness, and :math:`\epsilon_\infty`.
   * - **Direct Integral** :math:`V_{\mathrm{dir}}`
     - ``v_dir = q_recomb @ W_resta @ q_eject``
     - Vector-matrix-vector contraction in :math:`O(N_{\mathrm{atoms}}^2)` time.
   * - **Exchange Integral** :math:`V_{\mathrm{exch}}`
     - ``v_exch = q_cross @ W_resta @ q_eject_alt``
     - Permuted transition charges contracted with Resta kernel.
   * - **Matrix Element** :math:`|M_{if}|^2`
     - ``(v_dir - v_exch)**2 + v_dir**2`` (spatial) or ``|v_dir - v_exch|^2`` (spinor)
     - Spin-channel summation over all final states.
   * - **Energy Mismatch** :math:`\Delta E`
     - ``dE = (eps[e_prime] - eps[e1]) - (eps[e2] - eps[h])``
     - Instantaneous quasiparticle orbital energy differences.
   * - **Line Shape** :math:`\rho(\Delta E)`
     - ``inv_sqrt2pi_sigma * np.exp(- 0.5 * (dE / sigma)**2)``
     - Gaussian or Marcus/Jortner FCWD evaluated at mismatch :math:`\Delta E`.
   * - **Total Rates** :math:`\Gamma_{eeh}, \Gamma_{XX}`
     - ``rate_eeh_fs = 2*pi/HBAR * sum(M_sq * rho)``
     - Scaled to nanosecond units: :math:`k_{\mathrm{ns}} = \Gamma_{\mathrm{fs}} \times 10^6`.

---
