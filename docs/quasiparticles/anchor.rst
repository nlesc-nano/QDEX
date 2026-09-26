Anchor calibration
==================

Part of :doc:`/quasiparticles/index`.

.. figure:: /_static/figures/anchor_model.svg
   :width: 100%
   :alt: anchor model

   Anchor-scaled model evaluated with the CdSe database entry (p = 2). The vacuum curve passes through the monomer GW anchor at R₀; below R₀ the radius is clamped. The square marks the 2 nm CdSe test cluster.


.. important::

   The anchor is an evGW\@PBE0 calculation of the smallest cluster (Cd₁₆Se₁₃Cl₆ for CdSe). It fixes
   the residual A of the two-anchor ``gw`` model and the residuals r_H, r_L of the Resta and DIM
   models. Its size scaling to larger dots, E_conf(R)/E_conf(R₀), is an assumption.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.estimate_gw_qp_gap``
* CLI: ``--qp_gap gw``, ``--qp-polarization sphere|legacy``, ``--qp-residual-power``, ``--kernel resta-sphere``
* YAML: ``physics.qp_gap``, ``physics.qp_polarization``, ``physics.qp_residual_power``

.. code-block:: python

   estimate_gw_qp_gap(coords, atom_symbols, material_name, eps_out, return_details=False, regularization_length_ang=1.0, residual_power=2.0, strict=False, polarization_model="sphere")


The two anchors
---------------

The two-anchor model ``gw`` (also ``sgw-anchor``) places the QP correction between a finite-cluster
calibration and a bulk reference:

.. list-table::
   :widths: 25 25 50
   :header-rows: 1

   * - Anchor Limit
     - Physical System
     - Theoretical Characterization
   * - **Anchor 1: Smallest Vacuum Anchor** (:math:`R_0`)
     - Monomer or smallest stoichiometric Wulff cluster
     - Relaxed cluster in vacuum: PBE single point → PBE0 single point → eigenvalue-self-consistent GW (``EV_GW_ITER 4``), i.e. evGW\@PBE0 (CP2K). The QP shifts are taken relative to the PBE eigenvalues, because the QDEX input orbitals are PBE.
   * - **Anchor 2: Bulk Limit** (:math:`R \to \infty`)
     - Periodic crystal
     - High-accuracy bulk :math:`G_0W_0` quasiparticle gap (:math:`E_g^{\mathrm{GW, bulk}}`), calibrated against experimental ARPES.


The Confinement Interpolation Formula
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The radius :math:`R` is ``R_eff_hull``: the equivalent-volume radius of the convex hull of the
inorganic core atoms, :math:`(3V/4\pi)^{1/3}`, plus 1.25 Å. The scissor added to a PBE eigenvalue
difference is:

.. math::

   \Delta_{\mathrm{GW}}(R)
   = \Delta_{\mathrm{bulk}}
   + P(R;\epsilon_\infty,\epsilon_{\mathrm{out}})
   + A\,s(R),
   \qquad \Delta_{\mathrm{bulk}} = E_g^{\mathrm{GW,bulk}} - E_g^{\mathrm{PBE,bulk}} .

**Polarization term.** :math:`P` is the classical self-polarization energy of the electron plus the
hole in a dielectric sphere. The sphere has :math:`\epsilon_\infty` inside and
:math:`\epsilon_{\mathrm{out}}` outside, and both carriers are averaged over the 1S envelope
:math:`|j_0(\pi r/R)|^2`:

.. math::

   P(R) = F(\epsilon_\infty,\epsilon_{\mathrm{out}})\,\frac{e^2}{R},\qquad
   F = \Big\langle \sum_{n\ge0}
   \frac{(\epsilon_\infty-\epsilon_{\mathrm{out}})(n+1)}{\epsilon_\infty\,[n\epsilon_\infty+(n+1)\epsilon_{\mathrm{out}}]}
   \Big(\frac rR\Big)^{2n}\Big\rangle_{1S}.

This is the image-charge series of Böttcher and Brus (J. Chem. Phys. 80, 4403 (1984)).

* The :math:`n=0` term is the Born term, :math:`1/\epsilon_{\mathrm{out}} - 1/\epsilon_\infty`.
* The higher multipoles add about 12 % for CdSe in vacuum: F = 0.937 against 0.839.
* F vanishes for :math:`\epsilon_{\mathrm{out}}=\epsilon_\infty`.

Tight-binding GW calculations of Si nanocrystals find that the finite-size self-energy
correction is dominated by this surface polarization term, evaluated with the *bulk*
:math:`\epsilon_\infty` inside (Delerue, Lannoo and Allan, PRL 84, 2457 (2000); PRL 90, 076803
(2003); PRB 68, 115411 (2003)). The term therefore uses :math:`\epsilon_\infty`, not a size-reduced
:math:`\epsilon(R)`. Its large-R limit, :math:`F e^2/R`, is exact classical electrostatics.

**Residual.** The anchor amplitude is fixed by the finite cluster opening in vacuum:

.. math::

   A = \Delta(R_0) - \Delta_{\mathrm{bulk}} - P(R_0;\epsilon_\infty,1),
   \qquad \Delta(R_0) = E_g^{\mathrm{GW,cluster}} - E_g^{\mathrm{PBE,cluster}} .

:math:`A` collects everything that is not classical polarization of a sharp sphere:

* reduced screening near the surface;
* non-locality of the polarization;
* exchange–correlation effects.

It is fitted in vacuum, kept unchanged in a solvent, and scaled with size by s(R) (below). For Cd₁₆Se₁₃Cl₆, :math:`A = -0.42` eV.

**Residual scaling** (``qp_residual_scaling``). By default the residual is scaled by the PBE confinement
energy of the cluster instead of a power of the radius:

.. math::

   A(R) = A\,\frac{E_{\mathrm{conf}}(R)}{E_{\mathrm{conf}}(R_0)},\qquad
   E_{\mathrm{conf}} = E_g^{\mathrm{PBE}}(\mathrm{cluster}) - E_g^{\mathrm{PBE}}(\mathrm{bulk}).

The non-classical part is taken to be mostly band stretching, i.e. an energy-dependent bulk GW
correction. That correction is proportional to how far the confined levels lie from the band edges,
and each cluster's own PBE gap measures that distance. The form needs no radius definition and no
power p; it is clipped to [0, 1]. ``qp_residual_scaling: power`` restores (R₀/R)^p. For CdSe the scale
is 0.41 at 2 nm (0.33 with p = 2) and 0.26 at 3.2 nm (0.10 with p = 2).

**Resta and DIM models.** The same idea applies to the Resta and DIM models. Running a model once on the anchor
cluster with ``--qp-anchor-calibrate`` stores its per-edge error against evGW in
``qdex/data/dw_anchor_residuals.json``. The key is material, model, self-energy, representation,
populations, Z and solvent term. Later runs add residual × E_conf(R)/E_conf(R₀) to all occupied
orbitals (HOMO residual) and all virtual orbitals (LUMO residual). ``qp_anchor_residual: off``
disables it.

With the default sphere solvent term, the calibrated gap residuals (LUMO minus HOMO residual) are:

* ΔCOHSEX: −0.54 eV (``sgw-resta``), −0.44 eV (``sgw-dim``), −0.83 eV (``evgw-resta``).
* Classical levels: −0.21 / −0.11 / −0.48 eV.

They have the same sign and size as the residual of the ``gw`` model (−0.42 eV). With the softened
Born term (``qp_solvent_term: born``), the ΔCOHSEX residuals are +0.08 / +0.18 / −0.15 eV. They are
smaller only because the missing l ≥ 1 image multipoles offset the non-classical correction.

**Exact anchor.** :math:`R_0` in ``MATERIAL_DB`` must be computed with the same radius definition as
the target cluster (``get_cluster_size_metrics``). The CdSe entry is 5.3133 Å, recomputed from the
anchor geometry ``tests/CdSe/1.2nm/geom.xyz`` (previously 5.258 Å). The model then returns the evGW
correction of the monomer exactly. The QP gap differs from the evGW gap only by the difference
between the PBE gap of the supplied MO file and that of ``MATERIAL_DB`` (5 meV).

Per-edge curves: HOMO/LUMO split
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each band edge gets its own two-anchor curve (``anchor_edge_curves``):

.. math::

   \delta_H(R) &= f_b\,\Delta_{\mathrm{bulk}} + \tfrac12 P(R) + A_H\,s(R) ,\\
   \delta_L(R) &= (1-f_b)\,\Delta_{\mathrm{bulk}} + \tfrac12 P(R) + A_L\,s(R) .

Here :math:`\delta_H` is the downward shift of the HOMO and :math:`\delta_L` the upward shift of the
LUMO.

* **Polarization.** Symmetric between electron and hole, since both see the same sphere.
* **Bulk opening.** Split by :math:`f_b` (the anchor-derived bulk HOMO fraction, :math:`f_b = d_{h0} / (d_{h0} + d_{l0}) \approx 41.2\%` for CdSe, matching first-principles bulk GW literature).
* **Residuals.** :math:`A_H` and :math:`A_L` are fitted so that each edge reproduces the evGW
  frontier shifts of the anchor: :math:`d_{h0} = 1.398` eV and :math:`d_{l0} = 1.997` eV for CdSe, i.e. 41.2% / 58.8%.
  For CdSe, :math:`A_H = -0.40` eV and :math:`A_L = -0.02` eV (summing to the gap residual :math:`A = -0.42` eV).

The asymmetry therefore sits in the non-classical residual. It is exact at :math:`R_0` and fades
toward the bulk split for large dots. The Resta and DIM models use the same split for absolute IP/EA by
default (``qp_edge_split: anchor``; see :doc:`/quasiparticles/models`).

The BSE kernel of the two-anchor model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The kernel ``resta-sphere`` (bulk Resta W plus the same sphere reaction field) is described in
:doc:`/excitons/screened_kernel`.

Legacy curve
~~~~~~~~~~~~

``--kernel resta`` restores the bulk-only kernel, and ``--qp-polarization legacy`` restores the old
curve:

.. math::

   \Delta_{\mathrm{GW}}^{\mathrm{legacy}}(R) = \Delta_{\mathrm{bulk}}
   + \frac{11.52\,\mathrm{eV\,\AA}\,(1/\epsilon_{\mathrm{out}} - 1/\epsilon_\infty)}{R+\ell}
   + A\left(\frac{R_0}{R}\right)^{p},\qquad \ell = 1\ \text{Å}.

Its prefactor is 0.8 of the Born term. That is 1.4 times smaller than the classical large-R limit
for CdSe in vacuum, and the 0.8 has no documented origin.

*CLI & YAML Invocation*:

.. code-block:: bash

   qdex --config config.yaml --qp_gap gw --eps-out 2.24          # sphere polarization, kernel resta-sphere
   qdex --config config.yaml --qp_gap gw --kernel resta           # bulk-only BSE kernel
   qdex --config config.yaml --qp_gap gw --qp-polarization legacy --kernel resta   # old behaviour
