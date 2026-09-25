Dipoles oscillators
===================

Part of :doc:`/spectroscopy/index`.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.namd.transient_absorption``
* Callable: ``qdex.namd.transient_absorption.compute_transient_absorption``
* CLI: ``--sigma``
* YAML: ``namd.transient_absorption``

.. code-block:: python

   compute_transient_absorption(times_fs: np.ndarray, populations: np.ndarray, E_pairs: np.ndarray, f_pairs: np.ndarray, i_pairs: np.ndarray, a_pairs: np.ndarray, sigma_ev: float=0.03, e_range: Optional[Tuple[float, float]]=None, n_e_points: int=300, include_se: bool=False, all_energies: Optional[np.ndarray]=None, state_degeneracy: float=2.0, include_esa: bool=False, f_elec_esa: Optional[np.ndarray]=None, e_elec_esa: Optional[np.ndarray]=None)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part4_excited_states/index.rst:369-395``

6. Transition Dipoles, Oscillator Strengths & Superradiance
-----------------------------------------------------------

The single-particle transition dipole moment between occupied orbital :math:`i` and virtual orbital :math:`a` is:

.. math::

   \boldsymbol{\mu}_{ia} = \langle \phi_i | \mathbf{r} | \phi_a \rangle = \sum_{\mu \nu} C_{\mu i} \, \mathbf{D}_{\mu \nu} \, C_{\nu a}

where :math:`\mathbf{D}_{\mu \nu} = \langle \chi_\mu | \mathbf{r} | \chi_\nu \rangle` is the AO dipole matrix evaluated analytically via Libint2.

Under full configuration interaction, the exciton transition dipole moment is a coherent linear superposition:

.. math::

   \boldsymbol{\mu}_S = \sum_{ia} X_{ia}^S \, \boldsymbol{\mu}_{ia}.

The corresponding dimensionless oscillator strength is:

.. math::

   f_S = \frac{2}{3} \, \Omega_S \, |\boldsymbol{\mu}_S|^2.

This coherent summation describes **superradiance** and intensity borrowing, where optical strength from high-energy transitions is transferred into the lowest bright exciton.

---
