Spinor hamiltonian
==================

Part of :doc:`/relativity/index`.

.. figure:: /_static/figures/soc_spinors.svg
   :width: 100%
   :alt: soc spinors

   Construction of the SOC spinor subspace: spatial MOs in the SOC window are doubled into spin-orbitals, H_SOC is added from the GTH projectors, and the resulting spinors define the spinor BSE.


.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.soc_utils``
* Callable: ``qdex.soc_utils.compute_spinor_subspace``
* CLI: ``--soc_flag, --gth_file``
* YAML: ``soc.enabled``, ``soc.window``

.. code-block:: python

   compute_spinor_subspace(atom_symbols, coords_ang, shells, C_AO, eps_Ha, S_AO, active_indices, gth_file, nthreads=1, soc_cache=None, assume_orthonormal=False, SC_AO=None, device='numpy', verbose=True)


4. Mathematical Formulation: Separable GTH Pseudopotential
----------------------------------------------------------

In the GTH relativistic pseudopotential formalism, the spin-orbit potential :math:`\hat{V}_{\mathrm{SO}}` is expressed as a sum of separable semi-local projectors centered on each atom :math:`I`:

.. math::

   \hat{V}_{\mathrm{SO}} = \sum_I \sum_{l=1}^{l_{\max}} \sum_{m, m'=-l}^{l} \sum_{i=1}^{N_{\mathrm{proj}}} \sum_{j=1}^{N_{\mathrm{proj}}} |p_i^{Ilm}\rangle \, k_{ij}^{Il} \, \left( \mathbf{L} \cdot \mathbf{S} \right)_{mm'} \, \langle p_j^{Ilm'}|

where:

* :math:`|p_i^{Ilm}\rangle = p_i^l(r) R_{lm}(\hat{r})` are the atom-centered HGH projectors (Hartwigsen, Goedecker, Hutter, PRB 58, 3641 (1998)), with the normalized radial functions

  .. math::

     p_i^l(r) = \frac{\sqrt{2}\, r^{l+2(i-1)} \, e^{-r^2/2r_l^2}}{r_l^{\,l+(4i-1)/2}\, \sqrt{\Gamma\!\left(l+\frac{4i-1}{2}\right)}}

  and :math:`R_{lm}` the real spherical harmonics in libint order (:math:`m=-l..l`; :math:`(y, z, x)` for :math:`l=1`). The overlaps :math:`\langle \chi_\mu | p_i^{Ilm}\rangle` are libint overlaps with a Gaussian of exponent :math:`1/2r_l^2`; the :math:`r^{2(i-1)}` factor of the higher projectors comes from exponent derivatives of the *unnormalized* Gaussian.
* :math:`k_{ij}^{Il}` are the SOC constants of the GTH pseudopotential, read from ``GTH_SOC_POTENTIALS.txt`` for the functional of the DFT calculation (entry ``GTH-<functional>-q<n>``, ``soc.gth_functional``, default PBE).
* :math:`\mathbf{L} = (L_x, L_y, L_z)` are the orbital angular momentum matrices in the same real-harmonic basis as the projectors. They are purely imaginary and antisymmetric (odd under time reversal), which makes every spinor level a Kramers pair.
* :math:`\mathbf{S} = \frac{1}{2} (\sigma_x, \sigma_y, \sigma_z)` are the Pauli spin matrices:

.. math::

   \sigma_x = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad
   \sigma_y = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}, \quad
   \sigma_z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}


Angular Momentum Matrices in Spherical Harmonics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The matrix elements of :math:`L_z`, :math:`L_+`, and :math:`L_-` in the complex spherical harmonic basis :math:`|l, m\rangle` (Condon-Shortley phases) are:

.. math::

   \langle l, m | L_z | l, m' \rangle = m \, \delta_{m, m'}

.. math::

   \langle l, m \pm 1 | L_\pm | l, m \rangle = \sqrt{l(l+1) - m(m \pm 1)}

The Cartesian components are obtained via:

.. math::

   L_x = \frac{1}{2} (L_+ + L_-), \quad L_y = \frac{1}{2i} (L_+ - L_-)

and transformed to the real harmonics :math:`R_{lm} = \sum_{m'} U_{mm'} Y_{lm'}`
(:math:`R_{l,m>0} = [Y_{l,-m} + (-1)^m Y_{lm}]/\sqrt{2}`, :math:`R_{l,m<0} = i[Y_{l,m} - (-1)^m Y_{l,-m}]/\sqrt{2}`):

.. math::

   \mathbf{L}^{\mathrm{real}}_\kappa = U^* \mathbf{L}_\kappa U^T


Two-Component Spinor Hamiltonian Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Expanding the total Hamiltonian in the two-component spinor basis :math:`(\alpha, \beta)`:

.. math::

   |\psi_k^{\mathrm{spinor}}\rangle = \sum_{m=1}^{N_{\mathrm{mo}}} \left[ U_{mk}^\alpha |\phi_m\rangle \otimes |\alpha\rangle + U_{mk}^\beta |\phi_m\rangle \otimes |\beta\rangle \right]

The complete Hamiltonian takes the :math:`2 N_{\mathrm{mo}} \times 2 N_{\mathrm{mo}}` block form:

.. math::

   \mathbf{H}_{\mathrm{total}} = \begin{pmatrix}
     \mathbf{H}_0 + \mathbf{H}_z & \mathbf{H}_x - i \mathbf{H}_y \\
     \mathbf{H}_x + i \mathbf{H}_y & \mathbf{H}_0 - \mathbf{H}_z
   \end{pmatrix}

that is :math:`\hat{V}_{\mathrm{SO}} = \sum_\kappa \sigma_\kappa \otimes \mathbf{H}_\kappa`, where :math:`\mathbf{H}_0 = \operatorname{diag}(\varepsilon_1^{\mathrm{DFT}}, \dots, \varepsilon_{N_{\mathrm{mo}}}^{\mathrm{DFT}})` is the diagonal matrix of spin-free Kohn-Sham orbital energies, and the Cartesian spin-orbit blocks in the molecular orbital active space are (the :math:`\frac{1}{2}` is :math:`\mathbf{S} = \boldsymbol{\sigma}/2`):

.. math::

   \mathbf{H}_\kappa = \frac{1}{2} \mathbf{B}_{\mathrm{mo}} \left( \mathbf{k} \otimes \mathbf{L}^{\mathrm{real}}_\kappa \right) \mathbf{B}_{\mathrm{mo}}^T
   = \frac{i}{2} \mathbf{B}_{\mathrm{mo}} \, \mathbf{A}_\kappa \, \mathbf{B}_{\mathrm{mo}}^T, \qquad \kappa = x, y, z

with :math:`\mathbf{A}_\kappa = \operatorname{Im}(\mathbf{k} \otimes \mathbf{L}^{\mathrm{real}}_\kappa)` real and antisymmetric, so each :math:`\mathbf{H}_\kappa` is purely imaginary and Hermitian.

Here, :math:`\mathbf{B}_{\mathrm{mo}} = \mathbf{C}_{\mathrm{act}}^T \mathbf{B}_{\mathrm{raw}}` represents the projection of the atomic orbital-to-projector overlap matrix into the active space.

For a single atom and a complete :math:`p` shell this gives the :math:`j = 1/2` doublet at :math:`-\lambda` and the :math:`j = 3/2` quartet at :math:`+\lambda/2`, with :math:`\lambda = \sum_{ij} b_i k_{ij} b_j` and :math:`b_i` the radial overlaps with the projectors (``tests/test_soc_integrals.py``).
