Environment-consistent QP and BSE (``qp_gap: env``)
===================================================

Part of :doc:`/quasiparticles/index`.

.. figure:: /_static/figures/dielectric_sphere.svg
   :width: 100%
   :alt: dielectric sphere

   A charged excitation polarizes the nanocrystal–solvent interface (left). In a neutral exciton the
   electron and hole self-energies and the induced electron–hole attraction nearly cancel (right).

Why a separate model
--------------------

A nanocrystal with permittivity :math:`\epsilon_{\mathrm{in}}` sits in a medium with a smaller
:math:`\epsilon_{\mathrm{out}}`. Two excitations must be treated differently:

* **Charged excitation (QP gap).** An added electron or hole induces a surface polarization charge of
  its own sign. The charge is repelled by this polarization, so the ionization potential rises and the
  electron affinity falls. The fundamental gap opens by roughly
  :math:`(e^2/R)(1/\epsilon_{\mathrm{out}}-1/\epsilon_{\mathrm{in}})`.
* **Neutral excitation (optical gap).** The electron and the hole each carry that self-energy. Each
  also feels the surface charge induced by the other, which strengthens their attraction. For
  overlapping electron and hole densities the two contributions nearly cancel, so the absorption onset
  depends only weakly on the solvent (Brus 1984; Delerue, Lannoo and Allan 2003).

The scissor-type QP models (``gw``, ``sgw-*``) put the surface term into the QP gap. Their BSE
kernels (Resta, DIM, ``xs``) contain only interior screening, so the cancellation cannot happen and
the optical gap inherits the whole solvent shift (see :doc:`/validation/environment_cancellation`).
``qp_gap: env`` fixes this by using **one** screened interaction for both quantities.

Theory
------

For two points inside a dielectric sphere of radius :math:`R`, the electrostatic Green's function is
the interior term plus a reaction (image) term:

.. math::

   W(\mathbf r,\mathbf r') = \frac{e^2}{\epsilon_{\mathrm{in}}|\mathbf r-\mathbf r'|}
   + \underbrace{\frac{e^2}{\epsilon_{\mathrm{in}}}\sum_{l=0}^{\infty}
     c_l\,\frac{r^l r'^l}{R^{2l+1}}\,P_l(\cos\theta)}_{W^{\mathrm{refl}}(\mathbf r,\mathbf r')},
   \qquad
   c_l=\frac{(l+1)(\epsilon_{\mathrm{in}}-\epsilon_{\mathrm{out}})}
             {l\,\epsilon_{\mathrm{in}}+(l+1)\,\epsilon_{\mathrm{out}}}.

Its limits are:

* The :math:`l=0` term is the Born charging energy
  :math:`e^2(1/\epsilon_{\mathrm{out}}-1/\epsilon_{\mathrm{in}})/R`.
* :math:`W^{\mathrm{refl}}` vanishes identically when :math:`\epsilon_{\mathrm{out}}=\epsilon_{\mathrm{in}}`.
* For a conducting exterior the series gives Kelvin's image charge.

QDEX keeps its atomistic interior kernel, the Resta-MNOK :math:`W^{\mathrm{int}}_{AB}`, for the first
term. It evaluates :math:`W^{\mathrm{refl}}_{AB}` between atoms.

**Quasiparticle energies.** Each orbital :math:`p` with Mulliken atomic populations :math:`q_A(p)`
receives the polarization self-energy

.. math::

   \Sigma^{\mathrm{pol}}_p=\tfrac12\sum_{AB}q_A(p)\,W^{\mathrm{refl}}_{AB}\,q_B(p),

.. math::

   \varepsilon^{QP}_i=\varepsilon^{DFT}_i-\Sigma^{\mathrm{pol}}_i\ \ (\text{occupied}),\qquad
   \varepsilon^{QP}_a=\varepsilon^{DFT}_a+\Delta_{\mathrm{bulk}}+\Sigma^{\mathrm{pol}}_a\ \ (\text{virtual}),

where :math:`\Delta_{\mathrm{bulk}}=E_g^{GW,\mathrm{bulk}}-E_g^{PBE,\mathrm{bulk}}`. This is the
Delerue–Lannoo–Allan decomposition: bulk self-energy correction plus surface polarization. The shifts
are state resolved, so the QP levels are not a rigid scissor. The bulk opening is placed on the
virtual manifold (occupied-fixed gauge), as for the scissor models.

**BSE.** The direct kernel uses :math:`W=W^{\mathrm{int}}+W^{\mathrm{refl}}`. The exchange term keeps
the bare interaction. The environmental part of an exciton energy is then

.. math::

   \delta E_{\mathrm{opt}}=\tfrac12\,q_e^TW^{\mathrm{refl}}q_e+\tfrac12\,q_h^TW^{\mathrm{refl}}q_h
   -q_h^TW^{\mathrm{refl}}q_e=\tfrac12\,(q_e-q_h)^TW^{\mathrm{refl}}(q_e-q_h),

which is small for a 1S–1S exciton and large for charge-separated or surface-trapped states. Nothing is
tuned to obtain the cancellation; it follows from using the same :math:`W^{\mathrm{refl}}` twice.

Implementation
--------------

* Module: ``qdex.environment``. Functions ``sphere_cavity``, ``reaction_field_matrix``,
  ``environment_qp_energies``.
* CLI: ``--qp_gap env`` together with ``--eps-out``. Options: ``--env-cavity-buffer`` (default 1.0 Å) and
  ``--env-anchor-residual``.
* YAML: ``physics.qp_gap: env``, ``physics.eps_out``, ``physics.env_cavity_buffer``,
  ``physics.env_anchor_residual``.
* Kernels: ``resta`` and ``xs-resta``. For ``xs-resta`` the atom-pair reaction field is expanded to AO
  blocks. sBSE and DIM carry their own ``eps_out`` treatment and are rejected.
* The SOC spinor Hamiltonian is built from the QP energies, so the spinor BSE keeps the state-resolved
  shifts.

Parameters and conventions:

* :math:`\epsilon_{\mathrm{in}}` is the material :math:`\epsilon_\infty` from ``MATERIAL_DB``.
* :math:`\epsilon_{\mathrm{out}}` must be the **optical** permittivity :math:`n^2` of the solvent for
  vertical excitations. Examples: hexane 1.89, toluene 2.24, chloroform 2.09, polystyrene 2.5, vacuum 1.
  A static permittivity would describe a solvent that has fully relaxed around the charge.
* The cavity is centred on the atomic centroid. Its radius is the largest atomic distance plus the
  buffer, so that every atom (ligands included) lies inside the sphere, which the multipole series
  requires. The number of multipoles is chosen automatically for a 10\ :sup:`-6` eV tolerance.
* ``--env-anchor-residual`` adds the residual :math:`A(R_0/R)^p` of the anchor model
  (:doc:`anchor`). This is the part of the monomer GW anchor not explained by the dielectric term; it
  represents the size dependence of the short-range self-energy. Off by default.

CdSe 2 nm results
-----------------

Cd\ :sub:`68`\ Se\ :sub:`55`\ Cl\ :sub:`26`, spin-free, Resta interior kernel, 25 × 25 active space. The
cavity radius is 10.97 Å and :math:`\Sigma^{\mathrm{pol}}_{\mathrm{HOMO}}/\Sigma^{\mathrm{pol}}_{\mathrm{LUMO}}`
= 0.568/0.557 eV in vacuum.

.. list-table::
   :header-rows: 1

   * - :math:`\epsilon_{\mathrm{out}}`
     - QP gap (eV)
     - S\ :sub:`1` (eV)
     - binding (eV)
   * - 1.0
     - 3.852
     - 2.501
     - 1.351
   * - 1.9
     - 3.223
     - 2.503
     - 0.720
   * - 2.24
     - 3.115
     - 2.503
     - 0.612
   * - 6.2 (matched)
     - 2.727
     - 2.500
     - 0.227

The optical gap changes by 3 meV while the QP gap changes by 1.1 eV. With the scissor ``gw`` model the
optical gap moves by 0.66 eV between vacuum and :math:`\epsilon_{\mathrm{out}}=2.4`. With
``--env-anchor-residual`` the optical gap is 2.69 eV. The experimental first-exciton window for this
size is about 2.7–2.95 eV, using the Yu et al. (2003) sizing curve.

Limitations
-----------

* **Spherical cavity.** Adequate for near-spherical dots; this CdSe cluster is 15 % anisotropic.
  Platelets and rods need an ellipsoidal or boundary-element cavity.
* **Interior screening.** Uses bulk :math:`\epsilon_\infty`. The size dependence of the interior
  self-energy and screening is represented only by the optional anchor residual.
* **Exchange.** The exchange term is not environment screened, so the reaction-field (Onsager) shift
  of bright states, proportional to :math:`\mu^2/R^3`, is neglected.
* **Continuum boundary.** The reaction field of point-like atomic charges diverges as an atom
  approaches the cavity surface. The 1 Å buffer keeps every atom at least that far inside.
