Environment
===========

Part of :doc:`/interactions/index`.

.. figure:: /_static/figures/dielectric_sphere.svg
   :width: 100%
   :alt: dielectric sphere

   Dielectric confinement. A single added carrier is repelled by its own surface polarization (QP gap opens); for a neutral exciton the electron and hole self-energies and the polarization-enhanced e–h attraction largely cancel.


.. important::

   The solvent reaction term is a softened dielectric-boundary model. Its sign reverses when the exterior electronic permittivity exceeds the bulk value. It is not the complete Green function of a dielectric sphere.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.build_gamma``
* CLI: ``--2e-integrals, --kernel, --eps-out``
* YAML: ``physics.2e-integrals, physics.kernel, physics.eps_out``

.. code-block:: python

   build_gamma(atom_symbols, coords, alpha, beta=0.0, eta_dict=HARDNESS_DICT)


9. Environmental Dielectric Polarization
----------------------------------------

In practical applications, colloidal quantum dots are dispersed in liquid solvents (toluene :math:`\epsilon_{\mathrm{out}} = 2.38`, hexane :math:`\epsilon_{\mathrm{out}} = 1.88`) or embedded in polymer or dielectric matrices.


Classical Image Charge Solvation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When a carrier charge :math:`q` is placed inside a sphere of core permittivity :math:`\epsilon_\infty` surrounded by solvent permittivity :math:`\epsilon_{\mathrm{out}}`, image charges are induced at the boundary:

.. math::

   \Delta_{\mathrm{pol}}(R, \epsilon_{\mathrm{out}}) = -\frac{e^2}{R} \left( \frac{\epsilon_\infty - \epsilon_{\mathrm{out}}}{\epsilon_\infty + \epsilon_{\mathrm{out}}} \right).

* **Vacuum** (:math:`\epsilon_{\mathrm{out}} = 1.0`): Image potential is repulsive; gap opens due to lack of external dielectric screening.
* **Dielectric Matching** (:math:`\epsilon_{\mathrm{out}} = \epsilon_\infty`): Boundary polarization vanishes identically.
* **High-Dielectric Matrix** (:math:`\epsilon_{\mathrm{out}} > \epsilon_\infty`): External screening reduces carrier charging energies, compressing the fundamental gap.

For a model in which the QP gap and the BSE share one environment-aware interaction, see :doc:`/quasiparticles/environment_model`.
