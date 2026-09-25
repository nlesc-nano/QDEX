Dft reference
=============

Part of :doc:`/quasiparticles/index`.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.estimate_gw_qp_gap``
* CLI: ``--qp_gap, --dynamic_z``
* YAML: ``physics.qp_gap, physics.dynamic_z``

.. code-block:: python

   estimate_gw_qp_gap(coords, atom_symbols, material_name, eps_out, return_details=False, regularization_length_ang=1.0, residual_power=2.0, strict=False)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part3_gw_scissor/index.rst:135-153``

3. The DFT Band Gap Problem
---------------------------

The failure of Kohn-Sham DFT with semi-local functionals (LDA, PBE) to predict fundamental band gaps arises from two fundamental physical sources:

1. **Missing Derivative Discontinuity** (:math:`\Delta_{xc}`):
   The true fundamental band gap is:

   .. math::

      E_g = \mathrm{IP} - \mathrm{EA} = E_g^{\mathrm{KS}} + \Delta_{xc}

   where :math:`\Delta_{xc} = \left. \frac{\delta E_{xc}}{\delta \rho} \right|_{N+\delta} - \left. \frac{\delta E_{xc}}{\delta \rho} \right|_{N-\delta}` is the integer discontinuity of the exchange-correlation functional. In standard LDA/GGA functionals, :math:`E_{xc}[\rho]` is a continuous function of density, identically yielding :math:`\Delta_{xc} = 0`. Consequently, :math:`E_g^{\mathrm{PBE}}` underestimates the true quasiparticle gap by 30–50%.

2. **Self-Interaction Error (SIE)**:
   In semi-local functionals, an electron spuriously interacts with its own charge density through the Hartree term. This unphysical electrostatic repulsion artificially destabilizes occupied states and over-delocalizes frontier wavefunctions.

---
