Transition density
==================

Part of :doc:`/exciton_analysis/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.nto``
* Callable: ``qdex.nto.analyze_nto_state``
* CLI: ``--nto, --nto-states``
* YAML: ``analysis.nto, analysis.nto_states``

.. code-block:: python

   analyze_nto_state(solver, vec, energy_ev, f_osc, state_index, coords, symbols, mu_ia=None, soc_U=None, top_n=3, context=None)


Solving the Bethe-Salpeter Equation yields excitation energies :math:`\Omega_S` and eigenvectors :math:`\mathbf{X}_S = (X_{ia}^S)`. However, understanding the physical nature of an exciton—whether it is a tightly bound Wannier-Mott exciton, a localized Frenkel exciton, a surface-trap state, or a spatial charge-transfer (CT) excitation—requires quantitative real-space wavefunction analysis.

``QDEX`` incorporates the rigorous **Plasser-Dreuw exciton analysis framework** along with **Natural Transition Orbitals (NTOs)** to decompose complex multi-configurational exciton wavefunctions into intuitive, publication-ready physical descriptors.


1. The Two-Particle Transition Density Matrix
---------------------------------------------

The two-particle exciton wavefunction is expanded in electron-hole configurations:

.. math::

   |\Psi_S\rangle = \sum_{i \in \mathrm{occ}} \sum_{a \in \mathrm{virt}} X_{ia}^S \, a_a^\dagger a_i |0\rangle

In real space, the transition density :math:`\gamma_0^S(\mathbf{r}_h, \mathbf{r}_e)` describes the joint probability amplitude of finding the hole at :math:`\mathbf{r}_h` and the electron at :math:`\mathbf{r}_e`:

.. math::

   \gamma_0^S(\mathbf{r}_h, \mathbf{r}_e) = \sum_{i \in \mathrm{occ}} \sum_{a \in \mathrm{virt}} X_{ia}^S \, \phi_i(\mathbf{r}_h) \, \phi_a(\mathbf{r}_e)


Reduced Hole and Electron Densities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Integrating out one of the quasi-particles yields the reduced spatial densities of the hole (:math:`\rho_h`) and electron (:math:`\rho_e`):

.. math::

   \rho_h(\mathbf{r}) = \int |\gamma_0^S(\mathbf{r}, \mathbf{r}_e)|^2 \, d\mathbf{r}_e = \sum_{i, j \in \mathrm{occ}} D_{ij}^h \, \phi_i(\mathbf{r}) \phi_j^*(\mathbf{r})

.. math::

   \rho_e(\mathbf{r}) = \int |\gamma_0^S(\mathbf{r}_h, \mathbf{r})|^2 \, d\mathbf{r}_h = \sum_{a, b \in \mathrm{virt}} D_{ab}^e \, \phi_a^*(\mathbf{r}) \phi_b(\mathbf{r})

where the single-particle reduced density matrices are:

.. math::

   \mathbf{D}^h = \mathbf{X} \mathbf{X}^\dagger, \quad \mathbf{D}^e = \mathbf{X}^\dagger \mathbf{X}


Vectorized Low-Memory Mulliken Population
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In standard implementations, computing :math:`\rho_h` and :math:`\rho_e` requires forming full :math:`N_{\mathrm{ao}} \times N_{\mathrm{ao}}` AO density matrices. In ``QDEX``, atomic populations are evaluated directly in the active molecular orbital subspace:

.. math::

   q_A^h = \sum_{\mu \in A} \operatorname{Re}\left[ \left( \mathbf{C}_{\mathrm{occ}} \mathbf{D}^h \right)_{\mu, :} (\mathbf{S} \mathbf{C}_{\mathrm{occ}})_{\mu, :}^* \right]

.. math::

   q_A^e = \sum_{\mu \in A} \operatorname{Re}\left[ \left( \mathbf{C}_{\mathrm{virt}} \mathbf{D}^e \right)_{\mu, :} (\mathbf{S} \mathbf{C}_{\mathrm{virt}})_{\mu, :}^* \right]

This vectorized reduction avoids large intermediate arrays and runs in milliseconds even for large nanoclusters.

