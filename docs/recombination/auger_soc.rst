Auger soc
=========

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

.. rubric:: From ``docs/part8_auger/index.rst:413-431``

6. Relativistic Spin-Orbit Coupling (SOC) in Auger Scattering
-------------------------------------------------------------

When heavy elements (Pb, Bi, I) are present, Spin-Orbit Coupling splits degenerate valence and conduction manifolds into two-component relativistic spinors:

.. math::

   \psi_k(\mathbf{r}) = \begin{pmatrix} \psi_k^\alpha(\mathbf{r}) \\ \psi_k^\beta(\mathbf{r}) \end{pmatrix}

``QDEX`` computes transition charges directly in the spinor basis:

.. math::

   q_A^{IJ} = \sum_{\mu \in A} \sum_{\nu=1}^{N_{\mathrm{ao}}} \left[ (U_{\mu I}^\alpha)^* S_{\mu \nu} U_{\nu J}^\alpha + (U_{\mu I}^\beta)^* S_{\mu \nu} U_{\nu J}^\beta \right]

Because the spinor coefficients mix spin-up and spin-down components, the anti-symmetrized amplitude :math:`M_{if} = V_{\mathrm{dir}} - V_{\mathrm{exch}}` **automatically captures all spin-flip Auger pathways**, eliminating empirical spin selection rules.

---
