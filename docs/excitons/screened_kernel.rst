The screened interaction in the BSE
===================================

Part of :doc:`/excitons/index`.

.. important::

   The direct term K\ :sup:`d` uses the W of the QP model, including its finite-size part ΔW. QDEX
   sets the kernel automatically: ``qp`` for the Resta and DIM models, ``resta-sphere`` for ``gw``.
   A kernel with a different W is rejected unless ``--allow-inconsistent-kernel`` is given.

1. Why the same W
-----------------

In many-body perturbation theory the BSE kernel is the functional derivative of the self-energy,
:math:`\Xi = \delta\Sigma/\delta G`. In the GW approximation with static W this gives

.. math::

   \Xi(1,2;3,4) = v(1,3)\,\delta(1,2)\,\delta(3,4) - W(1,2)\,\delta(1,3)\,\delta(2,4),

that is, the exchange term with the bare v and the direct term with **the same W as in Σ**. The
QP correction contains the self-image of each carrier (:doc:`/quasiparticles/gw`, section 6). The
direct term contains the image of the partner carrier. Both are parts of one polarization, so the
QP gap and the exciton binding must come from one W. What happens when they do not is shown in
:doc:`cancellation`.

2. Kernel of each QP model
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - ``qp_gap``
     - ``kernel``
     - W in K\ :sup:`d`
   * - ``brus``
     - ``resta`` (recommended)
     - :math:`W^{\mathrm{bulk}}_{AB} = S_{\epsilon_\infty}(r_{AB})\,\gamma_{AB}`: no ΔW, like the QP gap
   * - ``gw``
     - ``resta-sphere``
     - :math:`W^{\mathrm{bulk}}_{AB} + G(\mathbf r_A,\mathbf r_B)`: the sphere of the QP correction
   * - ``sgw-resta``, ``evgw-resta``, ``qsgw-resta``
     - ``qp``
     - :math:`W^{\mathrm{bulk}} + \bar Z\,(W^{\mathrm{QD}} + W^{\mathrm{add}} - W^{\mathrm{bulk}})`, Resta ε_in(R)
   * - ``sgw-dim``, ``evgw-dim``, ``qsgw-dim``
     - ``qp``
     - the same with the DIM :math:`W^{\mathrm{QD}}`

The terms are those of :doc:`/quasiparticles/models`. For ``brus`` the default kernel is still ``bse``
(uniform 1/ε∞), which also screens on-site interactions and underbinds. ``resta`` is the bulk W with
the correct short range.

**Z in the kernel.** In the Resta and DIM models, ΔW in the kernel is scaled by the mean quasiparticle
weight of the frontier orbitals,

.. math::

   W^{\mathrm{BSE}} = W^{\mathrm{bulk}} + \bar Z\,\Delta W^{\mathrm{kernel}},\qquad \bar Z = \tfrac12(Z_H + Z_L).

The same plasmon pole that reduces the QP shift by Z also makes the electron–hole interaction
dynamical. To first order the two reductions cancel in a neutral excitation (Bechstedt et al., PRL
78, 1528 (1997)). With Z̄ in the kernel the static BSE keeps that cancellation. For Z = 1 the kernel is
:math:`W^{\mathrm{QD}} + W^{\mathrm{add}}`.

``evgw`` and ``qsgw`` pass their converged ΔW; ``qsgw`` also passes its relaxed orbitals and energies.

3. K\ :sup:`x` and K\ :sup:`d` in the two representations
---------------------------------------------------------

With the representation of :doc:`/quasiparticles/representation`:

.. math::

   K^x_{ia,jb} = (ia|v|jb),\qquad K^d_{ia,jb} = (ij|W^{\mathrm{BSE}}|ab).

* ``mnok``:

  .. math::

     K^x_{ia,jb} = \sum_{AB} q^{ia}_A\,\gamma_{AB}\,q^{jb}_B,\qquad
     K^d_{ia,jb} = \sum_{AB} q^{ij}_A\,W^{\mathrm{BSE}}_{AB}\,q^{ab}_B .

* ``xs``:

  .. math::

     K^x_{ia,jb} = \sum_{\mu\nu} c_{\mu i}c_{\mu a}\,(\mu\mu|\nu\nu)\,c_{\nu j}c_{\nu b},\qquad
     K^d_{ia,jb} = \sum_{\mu\nu} c_{\mu i}c_{\mu j}\,W^{\mathrm{BSE}}_{\mu\nu}\,c_{\nu a}c_{\nu b}.

Exchange is always bare: it describes the annihilation and re-creation of the pair, which the
other electrons do not screen in the BSE kernel. It is short range for band-edge excitons and sets
the singlet–triplet (bright–dark) splitting.

4. How the models behave
------------------------

**Binding energy E_b = E_g^QP − S₁** at 2 nm (CdSe, 25 × 25, spin-free, calibrated):

.. list-table::
   :header-rows: 1

   * -
     - E_b vacuum (eV)
     - E_b toluene (eV)
   * - ``gw`` + ``resta-sphere``
     - 1.58
     - 0.69
   * - ``sgw-resta``
     - 1.55
     - 0.73
   * - ``sgw-dim``
     - 1.58
     - 0.76
   * - ``evgw-resta``
     - 1.67
     - 0.78
   * - ``qsgw-dim``
     - 1.58
     - 0.76

* **The binding is mostly the mutual image.** Without the surface term, bulk W binds by only
  0.2–0.3 eV at 2 nm. With ε_out = 1 the image roughly doubles the attraction, and in toluene it is
  smaller.
* **All shared-W models bind alike.** They share the sphere term, which dominates. The interior part
  of ΔW adds a little (Resta, DIM), and evGW slightly more because its converged ε_in is smaller.
* **The QP gap and E_b move together** with the solvent (about 0.9 eV at 2 nm), so S₁ moves by only
  0.04–0.12 eV (:doc:`cancellation`).
* **Resta vs DIM.** For spherical dots the difference is below 0.1 eV. DIM matters when the surface
  screening is inhomogeneous: shape, ligands, shells.
* **mnok vs xs.** xs binds more at short range: −0.6 eV in S₁ at 1.2 nm, −0.05 eV at 2 nm.

5. Kernels for models without W
-------------------------------

``pbe``, a numeric gap and ``gw`` with ``qp_polarization: legacy`` define no W. The kernel can then be
chosen freely. None of these kernels contains the ΔW of a QP model:

* ``resta`` / ``xs-resta``: bulk Resta W, :math:`S_{\epsilon_\infty}(r)\,\gamma`.
* ``dim`` / ``xs-dim``: DIM screening :math:`S_{\epsilon_{AB}}(r)\,\gamma`, without the environment term.
* ``rpa`` / ``xs-rpa``: :math:`W = (1 + \Gamma\Pi^0)^{-1}\Gamma` from the monopole (ZDO) independent-particle
  polarizability of the active space.
* ``sbse``: the simplified BSE screening of Cho, Bintrim and Berkelbach. It underscreens CdSe
  (ε_eff ≈ 1).
* ``bse``: uniform :math:`\alpha\,\gamma/\epsilon_\infty`.

``--allow-inconsistent-kernel`` combines a W-defining QP model with one of these; use it only to
reproduce old results.
