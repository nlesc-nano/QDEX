Quasiparticle weight Z
======================

Part of :doc:`/quasiparticles/index`.

.. important::

   ``qp_z: derived`` is the default for every Resta and DIM model (``sgw-*``, ``evgw-*``, ``qsgw-*``, ``sgw``). Z is
   evaluated from one plasmon pole whose frequency follows from the same dielectric constant that
   builds the model's W. It is a linearized one-pole estimate, not the derivative of a computed
   frequency-dependent self-energy. A fixed value (``qp_z: 1.0``) is still accepted.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.hardness``
* Callable: ``qdex.hardness.compute_dynamic_z``
* CLI: ``--qp-z derived|<number>``
* YAML: ``physics.qp_z``

.. code-block:: python

   compute_dynamic_z(delta_sigma_stat_ev, gap_ev=None, eps_eff=None, material_name=None, omega_p_ev=None)


One-pole quasiparticle weight
-----------------------------

The finite-size correction of a Resta or DIM model is the static ΔCOHSEX self-energy
:math:`\sigma_p^{\mathrm{stat}} = \Delta\Sigma_p` (:doc:`gw`, section 5). Write the
frequency dependence of :math:`\Delta W` as a single pole of frequency :math:`\tilde\omega`:
:math:`\Delta\Sigma_p(\omega)` then has static value :math:`\sigma_p^{\mathrm{stat}}` and slope
:math:`\partial\Delta\Sigma_p/\partial\omega = -|\sigma_p^{\mathrm{stat}}|/\tilde\omega`. The linearized QP
equation gives

.. math::

   Z_p = \left(1 + \frac{|\sigma_p^{\mathrm{stat}}|}{\tilde\omega}\right)^{-1},
   \qquad \Delta\varepsilon_p = Z_p\,\sigma_p^{\mathrm{stat}}.

The pole is the generalized plasmon pole of Hybertsen and Louie (PRB 34, 5390 (1986)). For the single
Lorentz oscillator :math:`\epsilon(\omega) = 1 + \omega_p^2/(E_P^2-\omega^2)` with
:math:`\epsilon(0)=\epsilon_{\mathrm{eff}}`, it is the zero of :math:`\epsilon(\omega)`:

.. math::

   \tilde\omega = \sqrt{E_P^2+\omega_p^2} = \frac{\omega_p}{\sqrt{1-1/\epsilon_{\mathrm{eff}}}} .

It is also the pole of the Born reaction field :math:`1/\epsilon_{\mathrm{out}} - 1/\epsilon(\omega)` of a
charged sphere, so the interior and the surface (solvent) parts of :math:`\Delta W` share it.

Inputs, all taken from the model itself:

* :math:`\omega_p`: the valence (s, p) plasmon of the bulk crystal, computed from the lattice
  constant and valence electron count (``valence_plasmon_ev``; 14.1 eV for CdSe). Materials without
  an entry in ``PLASMON_DATA`` use 15 eV.
* :math:`\epsilon_{\mathrm{eff}}`:
  - the Penn :math:`\epsilon_{\mathrm{in}}(R)` for the Resta family;
  - the median inter-site screening for DIM;
  - the exciton screening of the sBSE response for ``sgw``.

Typical values for CdSe are Z ≈ 0.94–0.95 in vacuum and ≈ 0.98 in toluene. Z is clipped to [0.5, 1].


Z in the BSE kernel
-------------------

The same Z̄ = ½(Z_H + Z_L) scales ΔW in the BSE kernel, so the QP correction and the direct term
see the same dynamical reduction (:doc:`/excitons/screened_kernel`).

Limits
------

* **One pole.** The bulk part of the dynamical renormalization is already contained in the
  tabulated bulk GW gap and is not recomputed.
* **Z(ε_out).** Z depends on the solvent through :math:`\sigma^{\mathrm{stat}}`. The effect on S₁ is
  small: for sgw-resta at 1.2 nm, S₁ moves by 75 meV between vacuum and toluene with the derived Z
  and by 84 meV at Z = 1. That residual comes from the incomplete cancellation of the solvent term
  in a very small cluster, not from Z.
* **Earlier versions** used :math:`\tilde\omega=\sqrt{15^2/(\epsilon_{\mathrm{eff}}-1)+E_g^2}`. That
  expression combines the Penn gap with the band gap, not the plasmon pole, and used a fixed
  15 eV.
