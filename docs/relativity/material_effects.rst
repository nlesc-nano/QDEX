Material effects
================

Part of :doc:`/relativity/index`.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.soc_utils``
* Callable: ``qdex.soc_utils.compute_spinor_subspace``
* CLI: ``--soc_flag, --gth_file``
* YAML: ``physics.soc, physics.soc_window_ev``

.. code-block:: python

   compute_spinor_subspace(atom_symbols, coords_ang, shells, C_AO, eps_Ha, S_AO, active_indices, gth_file, nthreads=1, soc_cache=None, assume_orthonormal=False, SC_AO=None, device='numpy', verbose=True)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part2_soc/index.rst:40-44``

2. Impact of SOC on Semiconductor Nanocrystals
----------------------------------------------

In lead halide perovskites, the frontier conduction band arises from the hybridization of empty lead :math:`6p` orbitals with halogen :math:`np` states. 


.. rubric:: From ``docs/part2_soc/index.rst:45-62``

Band Inversion & Giant Gap Contraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the absence of SOC, the :math:`p`-like conduction band is triply degenerate (six-fold degenerate including spin). The operator :math:`\mathbf{L} \cdot \mathbf{S}` acts on total angular momentum :math:`\mathbf{J} = \mathbf{L} + \mathbf{S}`:

.. math::

   \mathbf{L} \cdot \mathbf{S} = \frac{1}{2} \left( \mathbf{J}^2 - \mathbf{L}^2 - \mathbf{S}^2 \right) = \frac{\hbar^2}{2} \left[ j(j+1) - l(l+1) - s(s+1) \right]

For :math:`p`-orbitals (:math:`l = 1, s = 1/2`):
* Lower :math:`j = 1/2` doublet: Eigenvalue :math:`-\hbar^2`. The conduction band minimum shifts downward in energy by :math:`\approx 0.65\text{ eV}` in :math:`\text{CsPbBr}_3`.
* Upper :math:`j = 3/2` quartet: Eigenvalue :math:`+\frac{1}{2}\hbar^2`.

As a direct consequence, standard non-relativistic DFT severely misidentifies the nature of the conduction band minimum. When SOC is activated:
1. The fundamental band gap contracts dramatically.
2. The effective mass of conduction band electrons decreases.
3. Carrier cooling dynamics within the conduction band manifold accelerate significantly due to dense non-adiabatic couplings between spinor levels.


.. rubric:: From ``docs/part2_soc/index.rst:63-75``

Rashba-Dresselhaus Splitting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At the surface of a quantum dot or in the presence of asymmetric ligand termination, the local inversion symmetry is broken. The resulting macroscopic gradient :math:`\boldsymbol{\nabla} V \neq \mathbf{0}` couples with electron momentum via the Rashba mechanism:

.. math::

   \hat{H}_{\mathrm{Rashba}} = \alpha_R \left( \boldsymbol{\sigma} \times \mathbf{k} \right) \cdot \hat{\mathbf{z}}

This splits spin-degenerate states in :math:`k`-space, protecting carriers against non-radiative electron-hole recombination and dramatically prolonging photoluminescence lifetimes.

---
