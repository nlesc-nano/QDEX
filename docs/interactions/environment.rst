Environment
===========

Part of :doc:`/interactions/index`.

.. important::

   The solvent reaction term is a softened dielectric-boundary model. Its sign reverses when the exterior electronic permittivity exceeds the bulk value. It is not the complete Green function of a dielectric sphere.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.build_gamma``
* CLI: ``--2e-integrals, --kernel, --eps-out``
* YAML: ``physics.2e-integrals, physics.kernel, physics.eps_out``

.. code-block:: python

   build_gamma(atom_symbols, coords, alpha, beta=0.0, eta_dict=HARDNESS_DICT)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part3_gw_scissor/index.rst:581-585``

9. Environmental Dielectric Polarization
----------------------------------------

In practical applications, colloidal quantum dots are dispersed in liquid solvents (toluene :math:`\epsilon_{\mathrm{out}} = 2.38`, hexane :math:`\epsilon_{\mathrm{out}} = 1.88`) or embedded in polymer or dielectric matrices.


.. rubric:: From ``docs/part3_gw_scissor/index.rst:586-600``

Classical Image Charge Solvation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When a carrier charge :math:`q` is placed inside a sphere of core permittivity :math:`\epsilon_\infty` surrounded by solvent permittivity :math:`\epsilon_{\mathrm{out}}`, image charges are induced at the boundary:

.. math::

   \Delta_{\mathrm{pol}}(R, \epsilon_{\mathrm{out}}) = -\frac{e^2}{R} \left( \frac{\epsilon_\infty - \epsilon_{\mathrm{out}}}{\epsilon_\infty + \epsilon_{\mathrm{out}}} \right).

* **Vacuum** (:math:`\epsilon_{\mathrm{out}} = 1.0`): Image potential is repulsive; gap opens due to lack of external dielectric screening.
* **Dielectric Matching** (:math:`\epsilon_{\mathrm{out}} = \epsilon_\infty`): Boundary polarization vanishes identically.
* **High-Dielectric Matrix** (:math:`\epsilon_{\mathrm{out}} > \epsilon_\infty`): External screening reduces carrier charging energies, compressing the fundamental gap.

---
