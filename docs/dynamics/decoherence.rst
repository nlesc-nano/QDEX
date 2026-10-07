Decoherence and dephasing times
===============================

Part of :doc:`/dynamics/index`.

.. rubric:: QDEX implementation

* Callable: ``qdex.namd.precompute.compute_trajectory_decoherence_times`` (pair times),
  ``qdex.namd.analysis.compute_band_gap_dynamics_and_spectral_density`` (lowest exciton)
* CLI: ``--namd-decoherence`` (after ``--namd-precompute``)
* YAML: ``namd.decoherence.method`` (``cumulant``, default, or ``gaussian``), ``max_lag_fs`` (1000),
  ``max_tau_fs`` (500)

Origin
------

In full quantum dynamics a superposition :math:`c_1|\psi_1\rangle + c_2|\psi_2\rangle` carries a nuclear
wavepacket on each surface,

.. math::

   |\Psi(t)\rangle = c_1 |\psi_1\rangle |\chi_1(t)\rangle + c_2 |\psi_2\rangle |\chi_2(t)\rangle ,
   \qquad \rho_{12}(t) \propto c_1 c_2^*\, \langle \chi_2(t) | \chi_1(t) \rangle .

The two surfaces exert different forces, the wavepackets separate, and the electronic coherence
:math:`\rho_{12}` decays with their overlap. A dot has thousands of vibrational modes, each modulating the
energy gap a little; the phases they add interfere destructively and the coherence is lost within
femtoseconds to tens of femtoseconds. In the classical path approximation all states share one classical
trajectory, so this loss does not happen by itself: it has to be added, as damping (FSSH), as stochastic
collapses (DISH) or as the line width of the golden-rule rates (PME), with a dephasing time for each pair
of states.

The dephasing function
----------------------

For a pair of states :math:`i, j` the gap :math:`\Delta\varepsilon_{ij}(t) = \varepsilon_i(t) - \varepsilon_j(t)`
fluctuates along the trajectory, :math:`\delta\Delta\varepsilon_{ij}(t) = \Delta\varepsilon_{ij}(t) - \langle
\Delta\varepsilon_{ij} \rangle`. In the second-order cumulant expansion the coherence between the two
states decays as

.. math::

   D_{ij}(t) = e^{-g_{ij}(t)}, \qquad
   g_{ij}(t) = \frac{1}{\hbar^2} \int_0^t dt_1 \int_0^{t_1} dt_2\; A_{ij}(t_2), \qquad
   A_{ij}(s) = \big\langle \delta\Delta\varepsilon_{ij}(t)\, \delta\Delta\varepsilon_{ij}(t+s) \big\rangle_t .

Two limits:

* **slow fluctuations** (the gap stays correlated for longer than :math:`\hbar/\sigma_{ij}`,
  :math:`\sigma_{ij}^2 = A_{ij}(0)`): :math:`A_{ij}(s) \approx \sigma_{ij}^2`, :math:`g_{ij} = \sigma_{ij}^2 t^2/2\hbar^2`,
  a Gaussian decay :math:`D_{ij}(t) = e^{-t^2/2\tau_{ij}^2}` with :math:`\tau_{ij} = \hbar/\sigma_{ij}`;
* **fast fluctuations** (correlation time :math:`\tau_c \ll \hbar/\sigma_{ij}`, motional narrowing): an
  exponential decay with rate :math:`\sigma_{ij}^2 \tau_c/\hbar^2`, slower than :math:`\hbar/\sigma_{ij}`
  suggests.

Pair dephasing times in QDEX
----------------------------

``--namd-decoherence`` evaluates the full cumulant for **every pair** of occupied states (hole channel)
and every pair of virtual states (electron channel) of the window, from the orbital energies of all
frames. The autocovariances of all pairs come from lagged covariance matrices of the energies,

.. math::

   C_{ij}(s) = \big\langle \delta\varepsilon_i(t)\, \delta\varepsilon_j(t+s) \big\rangle_t , \qquad
   A_{ij}(s) = C_{ii}(s) + C_{jj}(s) - C_{ij}(s) - C_{ji}(s),

one matrix product per lag; :math:`g_{ij}` is accumulated lag by lag (trapezoid rule) up to
``max_lag_fs``, and the time :math:`t^{1/e}_{ij}` where :math:`D_{ij}` falls to :math:`1/e` is located by
interpolation. The stored time keeps the Gaussian-width convention of DISH (Jaeger, Fischer & Prezhdo,
J. Chem. Phys. 137, 22A545 (2012)),

.. math::

   \tau_{ij} = t^{1/e}_{ij} / \sqrt{2} ,

which equals :math:`\hbar/\sigma_{ij}` in the slow-fluctuation limit and is longer for pairs whose gap
decorrelates fast. Pairs that do not reach :math:`1/e` within ``max_lag_fs``, and the diagonal, get
``max_tau_fs`` (500 fs); times are clipped to [1, 500] fs. ``method: gaussian`` stores
:math:`\hbar/\sigma_{ij}` directly; both are written to ``decoherence_times.npz`` (``tau_occ``,
``tau_virt``, ``tau_*_gaussian``, ``sigma_*``). The cost is about 1 s for 1000 states and 1000 frames.

**CsPbBr**\ :sub:`3` **(978 frames, 2 fs).** Hole pairs: median 12.4 fs (:math:`\hbar/\sigma` 12.4 fs),
electron pairs: 24.9 fs (24.2 fs); the cumulant and the short-time estimate agree within 0–10 % and every
pair reaches :math:`1/e` within 1 ps: the gap fluctuations of these soft lattices are large and slow
enough for the Gaussian limit. The pair times need the full trajectory: from 20 frames (38 fs) they come
out 5–10 times longer, because the variance misses the slow modes.

Where the times are used
------------------------

* **FSSH** (:doc:`fssh`): inactive amplitudes are damped by :math:`e^{-\Delta t/\tau_{IJ}}` per nuclear step.
* **DISH** (:doc:`dish`): a decoherence event of state :math:`J` occurs with probability
  :math:`1 - e^{-\Delta t/\tau_{IJ}}` per nuclear step.
* **PME** (:doc:`pme`): line width of the golden-rule rate of each pair (``pme_tau: pairs``).

Here :math:`I` is the active state of the channel. Without ``decoherence_times.npz`` a warning is
printed and all pairs get one time (``tau_dec_fs``, with ``decoherence: cumulant`` the lowest-exciton
cumulant time, otherwise 20 fs); the PME then uses that time or, with
``tau_dec_fs: edc``, the energy-based estimate of Granucci and Persico,
:math:`\tau_{IJ} = \hbar/|E_I - E_J|\,(1 + C/E_{\mathrm{kin}})`, :math:`C` = 0.1 Hartree,
:math:`E_{\mathrm{kin}} = \tfrac{3}{2} N_{\mathrm{atoms}} k_B T`.

The lowest exciton
------------------

``decoherence: cumulant`` also evaluates the cumulant of the lowest exciton energy along the trajectory
(normalised autocorrelation, :math:`D(t) = e^{-g(t)}`, time of :math:`D = 1/e`). It is a diagnostic of the
optical dephasing of the band-edge transition, reported with the spectral density of the gap
fluctuations (:doc:`analysis`); with fewer frames than the decay needs it is not reached and the log says
so.
