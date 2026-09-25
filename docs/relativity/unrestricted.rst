Unrestricted
============

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

.. rubric:: From ``docs/part2_soc/index.rst:205-226``

6. Spinor Representation & Unrestricted Kohn-Sham (UKS)
-------------------------------------------------------

Diagonalization of :math:`\mathbf{H}_{\mathrm{total}}` yields the relativistic spinor eigenvalues :math:`\varepsilon_k^{\mathrm{spinor}}` and the unitary expansion matrix :math:`\mathbf{U}`:

.. math::

   \mathbf{H}_{\mathrm{total}} \mathbf{U} = \mathbf{U} \operatorname{diag}\left( \varepsilon_1^{\mathrm{spinor}}, \dots, \varepsilon_{2N_{\mathrm{act}}}^{\mathrm{spinor}} \right)

The :math:`k`-th spinor wavefunction is partitioned into its :math:`\alpha` and :math:`\beta` spin components:

.. math::

   \psi_k^\alpha(\mathbf{r}) = \sum_m U_{mk}^\alpha \phi_m(\mathbf{r}), \quad
   \psi_k^\beta(\mathbf{r}) = \sum_m U_{mk}^\beta \phi_m(\mathbf{r})

The local spinor probability density is given by:

.. math::

   \rho_k^{\mathrm{spinor}}(\mathbf{r}) = |\psi_k^\alpha(\mathbf{r})|^2 + |\psi_k^\beta(\mathbf{r})|^2


.. rubric:: From ``docs/part2_soc/index.rst:227-242``

UKS Spin-Preserving Framework
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When starting from an Unrestricted Kohn-Sham (UKS) calculation with different molecular orbitals for alpha and beta spins (:math:`\mathbf{C}_\alpha \neq \mathbf{C}_\beta`), ``QDEX`` maps the active spaces independently:

.. math::

   \mathbf{H}_{\mathrm{total}}^{\mathrm{UKS}} = \begin{pmatrix}
     \mathbf{H}_{0, \alpha} - \frac{1}{2} \mathbf{H}_{z, \alpha \alpha} & -\frac{1}{2} (\mathbf{H}_{x, \alpha \beta} - i \mathbf{H}_{y, \alpha \beta}) \\
     -\frac{1}{2} (\mathbf{H}_{x, \beta \alpha} + i \mathbf{H}_{y, \beta \alpha}) & \mathbf{H}_{0, \beta} + \frac{1}{2} \mathbf{H}_{z, \beta \beta}
   \end{pmatrix}

This allows studying doped quantum dots, open-shell radicals, or spin-polarized nanocrystals without loss of relativistic accuracy.

---
