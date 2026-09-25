Anchor
======

Part of :doc:`/quasiparticles/index`.

.. important::

   ``sgw-anchor`` interpolates a bulk gap correction and one finite vacuum anchor. Its exponent and regularization length are model parameters, not uniquely fixed by two endpoints.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.estimate_gw_qp_gap``
* CLI: ``--qp_gap, --dynamic_z``
* YAML: ``physics.qp_gap, physics.dynamic_z``

.. code-block:: python

   estimate_gw_qp_gap(coords, atom_symbols, material_name, eps_out, return_details=False, regularization_length_ang=1.0, residual_power=2.0, strict=False)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part3_gw_scissor/index.rst:189-207``

6. Avenue 1: Two-Anchor Scaled GW (``sgw-anchor``)
--------------------------------------------------

Instead of guessing empirical parameters, ``sgw-anchor`` (historically called ``gw``) anchors the quasiparticle gap between a finite-cluster calibration and a bulk reference:

.. list-table::
   :widths: 25 25 50
   :header-rows: 1

   * - Anchor Limit
     - Physical System
     - Theoretical Characterization
   * - **Anchor 1: Smallest Vacuum Anchor** (:math:`R_0`)
     - Monomer or smallest stoichiometric Wulff cluster
     - Relaxed cluster computed in vacuum with hybrid DFT (:math:`\text{PBE0}`) and eigenvalue-self-consistent GW (``EV_GW_ITER 4``).
   * - **Anchor 2: Bulk Limit** (:math:`R \to \infty`)
     - Periodic crystal
     - High-accuracy bulk :math:`G_0W_0` quasiparticle gap (:math:`E_g^{\mathrm{GW, bulk}}`), calibrated against experimental ARPES.


.. rubric:: From ``docs/part3_gw_scissor/index.rst:208-241``

The Confinement Interpolation Formula
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The radius :math:`R` is the equivalent-volume radius of the inorganic core, :math:`R = (3V/4\pi)^{1/3}`. The scissor added to a PBE eigenvalue difference is:

.. math::

   \Delta_{\mathrm{GW}}(R)
   = \Delta_{\mathrm{bulk}}
   + \frac{\kappa_{\mathrm{out}}}{R+\ell}
   + A\left(\frac{R_0}{R}\right)^{p}

with the bulk shift defined by:

.. math::

   \Delta_{\mathrm{bulk}} = E_g^{\mathrm{GW,bulk}} - E_g^{\mathrm{PBE,bulk}}.

The dielectric prefactor is:

.. math::

   \kappa_{\mathrm{out}} = 11.52\,\mathrm{eV\,\AA}\left(\frac{1}{\epsilon_{\mathrm{out}}} - \frac{1}{\epsilon_\infty}\right).

In vacuum (:math:`\epsilon_{\mathrm{out}}=1`), :math:`\kappa_{\mathrm{vac}} = 11.52\,\mathrm{eV\,\AA}\,(1 - 1/\epsilon_\infty)`. The standard defaults are regularization length :math:`\ell = 1.0\,\text{Å}` and power :math:`p = 2`.

The anchor amplitude :math:`A` is fixed by the finite cluster opening:

.. math::

   A = \Delta(R_0) - \Delta_{\mathrm{bulk}} - \frac{\kappa_{\mathrm{vac}}}{R_0+\ell}

where :math:`\Delta(R_0) = E_g^{\mathrm{GW,cluster}} - E_g^{\mathrm{PBE,cluster}}`.


.. rubric:: From ``docs/part3_gw_scissor/index.rst:242-260``

Implementation in QDEX (``sgw-anchor``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The two-anchor scaled GW model is implemented in :func:`qdex.hardness.estimate_gw_qp_gap`:

1. **Nanocrystal Core Sizing**: Evaluates the equivalent spherical core radius :math:`R = (3V/4\pi)^{1/3}` or gyro-radius from the 3D atomic coordinates :math:`\mathbf{R}_A` of the inorganic core.
2. **Database Query**: Queries tabulated bulk and monomer anchor parameters from ``MATERIAL_DB`` (:math:`E_g^{\mathrm{PBE, bulk}}, E_g^{\mathrm{GW, bulk}}, R_0, \Delta(R_0), f_{\mathrm{homo}}, f_{\mathrm{lumo}}`).
3. **Electrostatic Image Charges**: Computes dielectric prefactors :math:`\kappa_{\mathrm{vac}} = 11.52 (1 - 1/\epsilon_\infty)` and :math:`\kappa_{\mathrm{out}} = 11.52 (1/\epsilon_{\mathrm{out}} - 1/\epsilon_\infty)` in :math:`\mathrm{eV\cdot\mathring{A}}`.
4. **Analytic Scissor Calculation**: Evaluates the total opening :math:`\Delta_{\mathrm{GW}}(R) = \Delta_{\mathrm{bulk}} + \frac{\kappa_{\mathrm{out}}}{R+\ell} + A (R_0/R)^p`.
5. **Level Alignment & Provenance**: Partitions :math:`\Delta_{\mathrm{GW}}` across occupied and virtual manifolds using database fractions :math:`f_{\mathrm{homo}}, f_{\mathrm{lumo}}` (Approach A) and stores detailed diagnostics in the returned ``details`` dictionary.

*CLI & YAML Invocation*:

.. code-block:: bash

   qdex --mos ground_state.mos --material CSPBBR3 --qp_gap sgw-anchor --eps-out 2.25

---
