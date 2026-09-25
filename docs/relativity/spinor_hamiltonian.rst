Spinor hamiltonian
==================

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

.. rubric:: From ``docs/part2_soc/index.rst:111-131``

4. Mathematical Formulation: Separable GTH Pseudopotential
----------------------------------------------------------

In the GTH relativistic pseudopotential formalism, the spin-orbit potential :math:`\hat{V}_{\mathrm{SO}}` is expressed as a sum of separable semi-local projectors centered on each atom :math:`I`:

.. math::

   \hat{V}_{\mathrm{SO}} = \sum_I \sum_{l=1}^{l_{\max}} \sum_{m, m'=-l}^{l} \sum_{i=1}^{N_{\mathrm{proj}}} \sum_{j=1}^{N_{\mathrm{proj}}} |p_i^{Ilm}\rangle \, k_{ij}^{Il} \, \left( \mathbf{L} \cdot \mathbf{S} \right)_{mm'} \, \langle p_j^{Ilm'}|

where:
* :math:`|p_i^{Ilm}\rangle` are atom-centered Gaussian-type projector functions with angular momentum :math:`(l, m)` and radial index :math:`i`.
* :math:`k_{ij}^{Il}` are material-specific relativistic coupling coefficients tabulated in ``GTH_SOC_POTENTIALS.txt``.
* :math:`\mathbf{L} = (L_x, L_y, L_z)` are orbital angular momentum operators in the complex spherical harmonic basis.
* :math:`\mathbf{S} = \frac{1}{2} (\sigma_x, \sigma_y, \sigma_z)` are the Pauli spin matrices:

.. math::

   \sigma_x = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad
   \sigma_y = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}, \quad
   \sigma_z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}


.. rubric:: From ``docs/part2_soc/index.rst:132-150``

Angular Momentum Matrices in Spherical Harmonics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The matrix elements of :math:`L_z`, :math:`L_+`, and :math:`L_-` in the standard complex spherical harmonic basis :math:`|l, m\rangle` are:

.. math::

   \langle l, m | L_z | l, m' \rangle = m \, \delta_{m, m'}

.. math::

   \langle l, m \pm 1 | L_\pm | l, m \rangle = \sqrt{l(l+1) - m(m \pm 1)}

The Cartesian components are obtained via:

.. math::

   L_x = \frac{1}{2} (L_+ + L_-), \quad L_y = \frac{1}{2i} (L_+ - L_-)


.. rubric:: From ``docs/part2_soc/index.rst:151-186``

Two-Component Spinor Hamiltonian Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Expanding the total Hamiltonian in the two-component spinor basis :math:`(\alpha, \beta)`:

.. math::

   |\psi_k^{\mathrm{spinor}}\rangle = \sum_{m=1}^{N_{\mathrm{mo}}} \left[ U_{mk}^\alpha |\phi_m\rangle \otimes |\alpha\rangle + U_{mk}^\beta |\phi_m\rangle \otimes |\beta\rangle \right]

The complete Hamiltonian takes the :math:`2 N_{\mathrm{mo}} \times 2 N_{\mathrm{mo}}` block form:

.. math::

   \mathbf{H}_{\mathrm{total}} = \begin{pmatrix}
     \mathbf{H}_0 - \frac{1}{2} \mathbf{H}_z & -\frac{1}{2} (\mathbf{H}_x - i \mathbf{H}_y) \\
     -\frac{1}{2} (\mathbf{H}_x + i \mathbf{H}_y) & \mathbf{H}_0 + \frac{1}{2} \mathbf{H}_z
   \end{pmatrix}

where :math:`\mathbf{H}_0 = \operatorname{diag}(\varepsilon_1^{\mathrm{DFT}}, \dots, \varepsilon_{N_{\mathrm{mo}}}^{\mathrm{DFT}})` is the diagonal matrix of spin-free Kohn-Sham orbital energies, and the Cartesian spin-orbit blocks in the molecular orbital active space are:

.. math::

   \mathbf{H}_x = \frac{1}{2} \mathbf{B}_{\mathrm{mo}} \, \mathbf{K}_x \, \mathbf{B}_{\mathrm{mo}}^T

.. math::

   \mathbf{H}_y = \frac{i}{2} \mathbf{B}_{\mathrm{mo}} \, \tilde{\mathbf{K}}_y \, \mathbf{B}_{\mathrm{mo}}^T

.. math::

   \mathbf{H}_z = \frac{1}{2} \mathbf{B}_{\mathrm{mo}} \, \mathbf{K}_z \, \mathbf{B}_{\mathrm{mo}}^T

Here, :math:`\mathbf{B}_{\mathrm{mo}} = \mathbf{C}_{\mathrm{act}}^T \mathbf{B}_{\mathrm{raw}}` represents the projection of the atomic orbital-to-projector overlap matrix into the active space.

---
