Foundations
===========

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

.. rubric:: From ``docs/part2_soc/index.rst:1-9``


In materials containing heavy chemical elements—such as lead halide perovskites (:math:`\text{CsPb}X_3`), bismuth double perovskites (:math:`\text{Cs}_2\text{AgBiBr}_6`), or heavy chalcogenide quantum dots (:math:`\text{InAs}`, :math:`\text{PbS}`, :math:`\text{CdSe}`)—non-relativistic Schrödinger quantum mechanics breaks down. 

Relativistic **Spin-Orbit Coupling (SOC)** fundamentally alters the electronic band structure, lifting orbital degeneracies, reshaping the density of states, and opening or contracting the optical band gap.

---


.. rubric:: From ``docs/part2_soc/index.rst:10-26``

1. Physical Origin of Spin-Orbit Coupling
-----------------------------------------

The spin-orbit interaction arises naturally from the relativistic **Dirac equation** for a spin-1/2 fermion in an electrostatic potential :math:`V(\mathbf{r})`. Under the low-velocity (Pauli) reduction to order :math:`(v/c)^2`, an electron moving with velocity :math:`\mathbf{v}` experiences the electric field of the nucleus :math:`\mathbf{E} = -\boldsymbol{\nabla} V(\mathbf{r})` transformed in its rest frame into an effective magnetic field :math:`\mathbf{B}_{\mathrm{eff}} = -\frac{1}{c} \mathbf{v} \times \mathbf{E}`.

The Zeeman coupling of the electron's intrinsic magnetic dipole moment :math:`\boldsymbol{\mu}_s = -g_s \frac{e}{2m_e} \mathbf{S}` to this relativistic field yields the spin-orbit Hamiltonian:

.. math::

   \hat{H}_{\mathrm{SO}} = -\boldsymbol{\mu}_s \cdot \mathbf{B}_{\mathrm{eff}} = \frac{\hbar}{4 m_e^2 c^2} \left( \boldsymbol{\nabla} V \times \mathbf{p} \right) \cdot \boldsymbol{\sigma}

For a central spherically symmetric potential :math:`V(r)` where :math:`\boldsymbol{\nabla} V(r) = \frac{\mathbf{r}}{r} \frac{dV}{dr}`, recognizing the orbital angular momentum operator :math:`\mathbf{L} = \mathbf{r} \times \mathbf{p}` and the spin operator :math:`\mathbf{S} = \frac{\hbar}{2} \boldsymbol{\sigma}` yields:

.. math::

   \hat{H}_{\mathrm{SO}} = \frac{1}{2 m_e^2 c^2} \frac{1}{r} \frac{dV(r)}{dr} \, \mathbf{L} \cdot \mathbf{S} = \xi(r) \, \mathbf{L} \cdot \mathbf{S}


.. rubric:: From ``docs/part2_soc/index.rst:27-39``

The :math:`Z^4` Scaling Law
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Near the atomic core, the Coulomb potential behaves as :math:`V(r) \approx -Z e^2 / r`, leading to :math:`\frac{1}{r} \frac{dV}{dr} \approx Z e^2 / r^3`. Evaluating the expectation value :math:`\langle r^{-3} \rangle` with hydrogenic wavefunctions produces a strong dependence on the atomic number :math:`Z`:

.. math::

   \langle \xi(r) \rangle \propto Z^4

This steep scaling explains why SOC is negligible for carbon (:math:`Z = 6`) or oxygen (:math:`Z = 8`), moderate for sulfur (:math:`Z = 16`), but overwhelmingly dominant in heavy elements such as iodine (:math:`Z = 53`), cesium (:math:`Z = 55`), lead (:math:`Z = 82`), and bismuth (:math:`Z = 83`).

---
