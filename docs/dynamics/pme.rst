Pauli master equation (PME)
===========================

Part of :doc:`/dynamics/index`.

.. rubric:: QDEX implementation

* Callable: ``qdex.namd.master_equation.propagate_pme_tensor``
* CLI: ``--namd-run`` (``--namd-method pme``)
* YAML: ``namd.dynamics.method: master_equation``, ``pme_tau`` (``pairs``, default, ``pairs_nac`` or ``uniform``),
  ``tau_dec_fs``, ``temperature_k``

Idea
----

When the coherences decay much faster than the populations change, they can be eliminated: the density
matrix stays diagonal and only the populations :math:`P_I(t)` are propagated, with rates from second-order
perturbation theory in the coupling. This is the master-equation limit of the dynamics; it gives
ensemble-averaged populations directly, without trajectories.

Rates
-----

A coherence that decays as :math:`e^{-t/\tau_{IJ}}` turns the energy-conserving delta function of Fermi's
golden rule into a Lorentzian of width :math:`\hbar/\tau_{IJ}`:

.. math::

   k_{I \to J}(t) = \frac{2\, |d_{IJ}(t)|^2\, \tau_{IJ}}{1 + \big((E_J - E_I)\, \tau_{IJ}/\hbar\big)^2}\; B_{IJ}(T),
   \qquad
   B_{IJ} = \begin{cases} 1 & E_J \le E_I \\ e^{-(E_J - E_I)/k_BT} & E_J > E_I \end{cases},

the same rate that DISH produces on average (:doc:`dish`). With ``pme_tau: pairs`` (default when
``decoherence_times.npz`` exists) :math:`\tau_{IJ}` are the pair dephasing times of :doc:`decoherence`, as
in DISH; ``pme_tau: uniform`` uses one time, ``tau_dec_fs`` (the lowest-exciton cumulant time with
``decoherence: cumulant``) or the energy-based estimate with ``tau_dec_fs: edc``. :math:`d` are the
couplings of the step (:doc:`states_couplings`); the rates are rebuilt at every nuclear step.

Fluctuating couplings (``pme_tau: pairs_nac``)
----------------------------------------------

The rate above assumes that the coupling keeps its phase for the whole dephasing time. For a coupling
that fluctuates, the golden rule reads

.. math::

   k_{I \to J} = 2 \langle |d_{IJ}|^2 \rangle\, \operatorname{Re} \int_0^\infty C_d(s)\, D_{IJ}(s)\, e^{i (E_J - E_I) s/\hbar}\, ds ,

with :math:`C_d` the normalised autocorrelation of the coupling; with exponential decays it keeps the
Lorentzian form with :math:`\tau_{IJ} \to \tau^{\mathrm{eff}}_{IJ} = (1/\tau_{IJ} + 1/\tau_c)^{-1}`.
``--namd-nac`` stores in ``nac_correlation.npz`` the correlation of each channel (:math:`C_d(s)`,
:math:`\tau_c`) and, for every pair, the fraction of the constant-coupling rate that survives the
fluctuations,

.. math::

   r_{IJ} = \frac{\sum_s w_s\, C_{IJ}(s)\, D_{IJ}(s)}{\sum_s w_s\, D_{IJ}(s)}, \qquad
   C_{IJ}(s) = \frac{\operatorname{Re}\langle d_{IJ}^*(t)\, d_{IJ}(t+s) \rangle_t}{\langle |d_{IJ}|^2 \rangle_t},
   \qquad D_{IJ}(s) = e^{-s^2/2\tau_{IJ}^2},

(trapezoid weights :math:`w_s`, lags up to 60 fs, :math:`C_{IJ}` cut at its first zero), accumulated
over the trajectory with a buffer of the last 30 couplings. ``pme_tau: pairs_nac`` uses
:math:`\tau^{\mathrm{eff}}_{IJ} = r_{IJ}\, \tau_{IJ}` (``r = 1`` for a coupling that keeps its phase), and
the channel :math:`\tau_c` only when the pair values are missing.

In the CsPbX\ :sub:`3` dots the couplings decorrelate within one or two frames (:math:`C_d(2\,\mathrm{fs})`
= 0.17–0.32 for the holes, 0.45–0.81 for the electrons; :math:`\tau_c` = 1.6–3.9 fs against
:math:`\tau_{IJ}` = 6–100 fs): a pair keeps 5–15 % of its constant-coupling rate (:math:`|d|^2`-weighted
mean of :math:`r_{IJ}`). Time to lose :math:`1/e` of the exciton excess energy (spin-free, pump
2 E\ :sub:`g`, 2 ps run from frame 0; values beyond 2 ps from the exponential fit):

.. list-table::
   :header-rows: 1

   * -
     - PME ``pairs``
     - PME ``pairs_nac``
     - DISH
     - FSSH
   * - CsPbCl\ :sub:`3`
     - 0.51 ps
     - 1.27 ps
     - 0.99 ps
     - 1.38 ps
   * - CsPbBr\ :sub:`3`
     - 0.94 ps
     - 2.6 ps
     - 1.64 ps
     - 2.6 ps
   * - CsPbI\ :sub:`3`
     - 1.20 ps
     - 3.5 ps
     - 2.1 ps
     - 3.7 ps

The constant-coupling PME is 2× faster than DISH; with the pair-resolved coupling correlation the PME
coincides with decoherence-corrected FSSH, and both are 1.3–1.7× slower than DISH. The channel-averaged
:math:`\tau_c` (one value per channel) gave 1.6, 3.3 and 5.4 ps.

Exciton populations from electron and hole rates
------------------------------------------------

The exciton populations form a matrix :math:`P_{ia}` (hole :math:`i`, electron :math:`a`). As the
couplings move one carrier at a time (:doc:`states_couplings`), an exciton relaxes by an electron step
:math:`a \to b` at fixed hole or a hole step :math:`i \to j` at fixed electron:

.. math::

   \frac{dP_{ia}}{dt} = \sum_{b \neq a} \big[K^e_{b \to a} P_{ib} - K^e_{a \to b} P_{ia}\big]
   + \sum_{j \neq i} \big[K^h_{j \to i} P_{ja} - K^h_{i \to j} P_{ia}\big] .

:math:`K^e` (:math:`n_{\mathrm{virt}}^2`) and :math:`K^h` (:math:`n_{\mathrm{occ}}^2`) are the rates above with
the virtual and occupied couplings and orbital energy differences: :math:`\varepsilon_b - \varepsilon_a` for
the electron, :math:`\varepsilon_i - \varepsilon_j` for the hole. The pair rate matrix
(:math:`N_{\mathrm{pairs}}^2 \approx 10^{12}` for :math:`10^6` pairs) is never built: the equation is
:math:`\dot{\mathbf{P}} = \mathbf{P}\mathbf{K}^e - \mathbf{P}\operatorname{diag}\mathbf{L}^e +
\mathbf{K}^{h\mathsf{T}}\mathbf{P} - \operatorname{diag}\mathbf{L}^h\,\mathbf{P}` with the loss vectors
:math:`\mathbf{L}`, matrix products of size :math:`n_{\mathrm{occ}}^2 n_{\mathrm{virt}} + n_{\mathrm{occ}} n_{\mathrm{virt}}^2`.

Detailed balance uses the orbital gaps, while FSSH and DISH use the pair energies with :math:`K^d` and
:math:`K^x`: the PME equilibrium differs from theirs by the spread of the electron–hole binding over the
window, which matters for the final thermal distribution, not for the cooling rates.

Integration
-----------

Each nuclear step is divided into 20 sub-steps. Each sub-step applies a row-stochastic matrix
:math:`T = I + \delta t\, K` of the electron channel and one of the hole channel,
:math:`\mathbf{P} \leftarrow \mathbf{P}\, T_e`, then :math:`\mathbf{P} \leftarrow T_h^\mathsf{T}\, \mathbf{P}`; a row
whose leaving probability would exceed 1 is renormalised onto its outgoing transitions. Populations stay
non-negative and the total is conserved. Recombination (:math:`e^{-k_{\mathrm{loss}}\Delta t}`, with
recombination enabled) is applied once per nuclear step.

Output: the population matrix in time, the mean exciton energy and the excess energies of the electron and
of the hole (:doc:`analysis`). One step for :math:`10^6` pairs takes a fraction of a second.
