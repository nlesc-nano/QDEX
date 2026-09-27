Anchor calibration
==================

Part of :doc:`/quasiparticles/index`.

.. important::

   The anchor is an evGW\@PBE0 calculation of the smallest cluster (Cd₁₆Se₁₃Cl₆ for CdSe). It fixes
   the residuals r_H, r_L of the Resta and DIM models and the bulk split f_b. Its size scaling to larger
   dots, E_conf(R)/E_conf(R₀), is an assumption.

.. rubric:: QDEX implementation

* Module: ``qdex.qp_levels`` (residual table), ``qdex.hardness`` (``anchor_residual_scale``,
  ``anchor_bulk_homo_fraction``, ``anchor_edge_curves``)
* YAML: ``quasiparticles.anchor_residual``, ``quasiparticles.residual_scaling``,
  ``quasiparticles.anchor_calibrate``, ``quasiparticles.anchor_table``, ``quasiparticles.edge_split``

The two references
------------------

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Reference
     - Calculation
   * - Anchor cluster (R₀)
     - Relaxed Cd₁₆Se₁₃Cl₆ in vacuum: PBE single point → PBE0 single point → eigenvalue-self-consistent
       GW (``EV_GW_ITER 4``), i.e. evGW\@PBE0 (CP2K). The QP shifts are taken relative to the PBE
       eigenvalues, because the QDEX input orbitals are PBE: d_h0 = 1.398 eV (HOMO down) and
       d_l0 = 1.997 eV (LUMO up) for CdSe.
   * - Bulk (R → ∞)
     - Bulk QP gap in ``MATERIAL_DB``; Δ_bulk = E_g^GW(bulk) − E_g^PBE(bulk)
       (:doc:`/reference/materials`).

**Bulk split.** The valence share of Δ_bulk is taken from the anchor,
:math:`f_b = d_{h0}/(d_{h0} + d_{l0})` (41.2 % for CdSe).

Residual of a Resta or DIM model
--------------------------------

Running a model once on the anchor cluster in vacuum,

.. code-block:: bash

   qdex --config config.yaml --eps-out 1.0 --qp-anchor-calibrate

stores its HOMO and LUMO errors against evGW,

.. math::

   r_H = d_{h0} - \big(f_b\,\Delta_{\mathrm{bulk}} + Z_H\,\Delta\Sigma_H\big),\qquad
   r_L = d_{l0} - \big((1-f_b)\,\Delta_{\mathrm{bulk}} + Z_L\,\Delta\Sigma_L\big),

in ``qdex/data/dw_anchor_residuals.json`` (or ``anchor_table``). The key is material, model,
self-energy, representation, populations, Z and solvent term, so each combination has its own
residual. Later runs add :math:`r_H\,s(R)` to all occupied and :math:`r_L\,s(R)` to all virtual
orbitals; ``anchor_residual: off`` disables it. Without a calibrated entry the run prints a notice and
adds nothing.

The residual contains what the static ΔCOHSEX misses (dynamics beyond one plasmon pole, band
stretching) and the PBE → PBE0 starting point of the reference.

.. note::

   The residual depends on Δ_bulk and f_b through the formula above. After a change of the bulk data
   in ``MATERIAL_DB`` the table must be recalibrated, otherwise the anchor is no longer reproduced.

Size scaling
------------

.. math::

   s(R) = \frac{E_{\mathrm{conf}}(R)}{E_{\mathrm{conf}}(R_0)},\qquad
   E_{\mathrm{conf}} = E_g^{\mathrm{PBE}}(\mathrm{cluster}) - E_g^{\mathrm{PBE}}(\mathrm{bulk}),
   \qquad s \in [0, 1].

The non-classical part is taken to be mostly band stretching, an energy-dependent bulk GW correction.
It grows with how far the confined levels lie from the band edges, which the cluster's own PBE gap
measures. The form needs no radius definition and no free exponent. For CdSe, s = 0.41 at 2 nm and
0.26 at 3.2 nm. ``residual_scaling: power`` uses (R₀/R)^p instead.

Absolute levels (IP/EA)
-----------------------

The gap correction of a Resta or DIM model is split between HOMO and LUMO with per-edge curves that
reproduce the anchor exactly (``edge_split: anchor``, default):

.. math::

   \delta_H(R) = f_b\,\Delta_{\mathrm{bulk}} + \tfrac12 P(R) + A_H\,s(R),\qquad
   \delta_L(R) = (1-f_b)\,\Delta_{\mathrm{bulk}} + \tfrac12 P(R) + A_L\,s(R),

with P(R) the classical sphere polarization (:doc:`models`, section 2), symmetric in electron and
hole, and A_H, A_L fixed at R₀. The asymmetry sits in the residuals and fades toward the bulk split
f_b for large dots. ``edge_split: model`` uses the model's own split.
