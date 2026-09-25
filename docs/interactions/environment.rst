Environment
===========

Part of :doc:`/interactions/index`.

.. figure:: /_static/figures/dielectric_sphere.svg
   :width: 100%
   :alt: dielectric sphere

   Dielectric confinement. A single added carrier is repelled by its own surface polarization (QP gap opens); for a neutral exciton the electron and hole self-energies and the polarization-enhanced e–h attraction largely cancel.


.. important::

   Two environment terms are in use. The Delta-W models (``sgw-*``, ``evgw-*``, ``qsgw-*``) use a
   softened Born term, (1/ε_out − 1/ε∞) e²/√(r² + R²); it is not the full Green function of a
   dielectric sphere. The two-anchor ``gw`` model uses the full multipole Green function of the sphere
   (:doc:`/quasiparticles/anchor`). In both cases the same term enters the QP correction and the BSE
   kernel (:doc:`/quasiparticles/theory`, section 4). Both change sign when ε_out exceeds ε∞.

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

