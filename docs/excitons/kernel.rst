The electron–hole interaction: K\ :sup:`x` and K\ :sup:`d`
==========================================================

Part of :doc:`/excitons/index`.

.. figure:: /_static/figures/bse_kernels.svg
   :width: 100%
   :alt: bse kernels

   The two interaction terms. Exchange uses the bare interaction and the transition densities; the
   direct term uses the screened W and the hole and electron densities.

.. important::

   The direct term K\ :sup:`d` uses the W of the QP model, including its finite-size part ΔW. QDEX
   sets the kernel automatically (``qp`` for the Resta and DIM models). A kernel with a different W is
   rejected unless ``allow_inconsistent_kernel`` is set.

1. The two terms
----------------

.. math::

   K^x_{ia,jb} = (ia|v|jb) = \iint \psi_i^*(\mathbf r)\psi_a(\mathbf r)\,v(\mathbf r,\mathbf r')\,
   \psi_j(\mathbf r')\psi_b^*(\mathbf r'),

.. math::

   K^d_{ia,jb} = (ij|W|ab) = \iint \psi_i^*(\mathbf r)\psi_j(\mathbf r)\,W(\mathbf r,\mathbf r')\,
   \psi_a(\mathbf r')\psi_b^*(\mathbf r').

* **Direct term K**\ :sup:`d` (attractive). The Coulomb attraction between the electron density
  :math:`\psi_a\psi_b^*` and the hole density :math:`\psi_i\psi_j^*`, screened by the rest of the
  system. It binds the exciton.
* **Exchange term K**\ :sup:`x` (repulsive). The annihilation of the pair ia and its re-creation as
  jb, through the transition densities :math:`\psi_i\psi_a`. It is unscreened: in the BSE kernel it
  comes from the bare Hartree term. It is short range for band-edge excitons and sets the
  singlet–triplet (dark–bright) splitting.

Spin enters through how the two terms combine:

.. math::

   \text{singlet: } 2K^x - K^d,\qquad \text{triplet: } -K^d,\qquad \text{spinor (SOC): } K^x - K^d .

The spin-adapted singlet couples to the pair annihilation with weight 2; the triplets have no net
transition density. With SOC the transitions are between Kramers spinors and the spin sum is already
in the orbitals.

``include_exchange: false`` or ``include_direct_eh: false`` switches one term off, to separate their
contributions.

2. Which W enters K\ :sup:`d`
-----------------------------

The BSE kernel is the functional derivative of the self-energy, :math:`\Xi = \delta\Sigma/\delta G`.
In the GW approximation with static W,

.. math::

   \Xi(1,2;3,4) = v(1,3)\,\delta(1,2)\,\delta(3,4) - W(1,2)\,\delta(1,3)\,\delta(2,4),

that is, exchange with the bare v and the direct term with **the same W as in Σ**. The QP energies
contain the self-image of each carrier; K\ :sup:`d` contains the image of the partner. Both are parts
of one polarization, so the QP energies and K\ :sup:`d` must come from one W.

.. list-table::
   :header-rows: 1
   :widths: 24 20 56

   * - ``quasiparticles.model``
     - ``kernel``
     - W in K\ :sup:`d`
   * - ``sgw-resta``, ``evgw-resta``, ``qsgw-resta``
     - ``qp``
     - :math:`W^{\mathrm{bulk}} + \bar Z\,(W^{\mathrm{QD}} + W^{\mathrm{add}} - W^{\mathrm{bulk}})`, Resta ε_in(R)
   * - ``sgw-dim``, ``evgw-dim``, ``qsgw-dim``
     - ``qp``
     - the same with the DIM :math:`W^{\mathrm{QD}}`
   * - ``bulk``, ``brus``
     - ``resta`` or ``dim``
     - :math:`W^{\mathrm{bulk}}`: no ΔW, like the QP energies (:doc:`sbse`)

The terms are those of :doc:`/quasiparticles/models`. ``evgw`` passes its converged ΔW; ``qsgw`` also
its relaxed orbitals and energies.

**Z in the kernel.** ΔW in the kernel is scaled by the mean quasiparticle weight of the frontier
orbitals, :math:`\bar Z = \tfrac12(Z_H + Z_L)`. The plasmon pole that reduces the QP shift by Z also
makes the electron–hole interaction dynamical; to first order the two reductions cancel in a neutral
excitation (Bechstedt et al., PRL 78, 1528 (1997)). With Z̄ in the kernel the static BSE keeps that
cancellation.

**Kernels for models without ΔW** (``none``, ``bulk``, ``brus``, a numeric gap):

* ``resta`` / ``xs-resta``: bulk Resta W, :math:`S_{\epsilon_\infty}(r)\,\gamma`.
* ``dim`` / ``xs-dim``: DIM screening :math:`S_{\epsilon_{AB}}(r)\,\gamma`, without the environment term.
* ``rpa`` / ``xs-rpa``: :math:`W = (1 + \Gamma\Pi^0)^{-1}\Gamma` from the monopole independent-particle
  polarizability of the active space.
* ``sbse``: the simplified-BSE screening of Cho, Bintrim and Berkelbach. It underscreens CdSe
  (ε_eff ≈ 1).
* ``bse``: uniform :math:`\alpha\,\gamma/\epsilon_\infty` (``kernel_scaling`` = α). It also screens
  on-site terms and underbinds.
* ``stda``: Grimme's γ\ :sup:`J` / γ\ :sup:`K`, set by ``mode: stda`` (:doc:`stda`).

3. Representation
-----------------

With the representation chosen in :doc:`/integrals/index`:

* ``mnok``:

  .. math::

     K^x_{ia,jb} = \sum_{AB} q^{ia}_A\,\gamma_{AB}\,q^{jb}_B,\qquad
     K^d_{ia,jb} = \sum_{AB} q^{ij}_A\,W_{AB}\,q^{ab}_B .

* ``xs``:

  .. math::

     K^x_{ia,jb} = \sum_{\mu\nu} c_{\mu i}c_{\mu a}\,(\mu\mu|\nu\nu)\,c_{\nu j}c_{\nu b},\qquad
     K^d_{ia,jb} = \sum_{\mu\nu} c_{\mu i}c_{\mu j}\,W_{\mu\nu}\,c_{\nu a}c_{\nu b}.

With the default on-site value γ_AA = IP − EA, MNOK reproduces xs within 0.15 eV at 1.2 nm and
0.01 eV at 2 nm (:doc:`/integrals/representation`). With the on-site hardness η instead, MNOK
underbinds by 0.4–0.5 eV at 1.2 nm.

4. Size of the binding
----------------------

Binding energy E_b = E_g^QP − S₁ at 2 nm (CdSe, 25 × 25, spin-free, coupled BSE):

.. list-table::
   :header-rows: 1

   * -
     - E_b vacuum (eV)
     - E_b toluene (eV)
   * - ``sgw-resta``
     - 1.52
     - 0.74
   * - ``sgw-dim``
     - 1.55
     - 0.77
   * - ``evgw-resta``
     - 1.65
     - 0.79
   * - ``qsgw-dim``
     - 1.55
     - 0.77

* **Most of the binding is the mutual image.** Bulk W alone binds by 0.2–0.3 eV at 2 nm.
* **All shared-W models bind alike.** The sphere term dominates; the interior part of ΔW adds a
  little.
* **The QP gap and E_b move together** with the solvent (about 0.85–1.00 eV), so S₁ moves by only
  0.07–0.12 eV (:doc:`cancellation`).
