Orbital cubes
=============

Part of :doc:`/electronic_structure/index`.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.pdos_coop``
* Callable: ``qdex.pdos_coop.compute_pdos_and_coop``
* CLI: ``--charge_type, --run_fuzzy``
* YAML: ``physics.charge_type, fuzzy.run``

.. code-block:: python

   compute_pdos_and_coop(C, S, eps_eV, shells, pdos_atoms, coop_pairs, ewin, sigma=0.03, is_soc=False, prefix='sf', pops=None, population_bars=None, device='numpy')

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part1_ground_state/index.rst:219-223``

6. Visualizing Orbitals via 3D Gaussian .cube Files
---------------------------------------------------

To inspect spatial orbital distributions, HOMO/LUMO wavefunctions, and surface states, ``QDEX`` generates standard volumetric Gaussian ``.cube`` files.


.. rubric:: From ``docs/part1_ground_state/index.rst:224-234``

Grid Discretization
~~~~~~~~~~~~~~~~~~~

Given the Cartesian bounding box of the system :math:`[\mathbf{r}_{\min} - \Delta r_{\mathrm{margin}}, \mathbf{r}_{\max} + \Delta r_{\mathrm{margin}}]`, a uniform 3D grid of spacing :math:`h` (default 0.5 Å) is constructed:

.. math::

   N_x = \left\lceil \frac{x_{\max} - x_{\min}}{h} \right\rceil, \quad
   N_y = \left\lceil \frac{y_{\max} - y_{\min}}{h} \right\rceil, \quad
   N_z = \left\lceil \frac{z_{\max} - z_{\min}}{h} \right\rceil


.. rubric:: From ``docs/part1_ground_state/index.rst:235-251``

C++ Libint Acceleration
~~~~~~~~~~~~~~~~~~~~~~~

Evaluating :math:`\phi_m(\mathbf{r}_g) = \sum_\mu C_{\mu m} \chi_\mu(\mathbf{r}_g)` across :math:`10^6` grid points in Python is computationally slow. ``QDEX`` uses an optimized multithreaded C++ evaluator (``libint_cpp.evaluate_mos_on_grid``) that processes the grid in contiguous chunks:

.. math::

   \phi_m(\mathbf{r}_g) \quad \text{for } g = 1, \dots, N_{\mathrm{grid}}

The resulting volumetric grids are written to:
* ``spatial_MO_HOMO.cube``, ``spatial_MO_LUMO.cube``: Signed scalar wavefunctions :math:`\phi_m(\mathbf{r})`.
* ``spinor_sp_HOMO_density.cube``: Relativistic spinor electron density :math:`\rho_k(\mathbf{r}) = |\psi_k^\alpha(\mathbf{r})|^2 + |\psi_k^\beta(\mathbf{r})|^2`.

These files can be loaded directly into **VMD**, **PyMOL**, or **ChimeraX** for publication-quality rendering.

---
