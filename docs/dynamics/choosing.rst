Choosing a method: FSSH, DISH or PME
====================================

Part of :doc:`/dynamics/index`.

All three methods use the same precomputed excitons, couplings and dephasing times; they differ in how
coherence and its loss are treated.

.. list-table::
   :widths: 20 27 27 26
   :header-rows: 1

   * -
     - CPA-FSSH (:doc:`fssh`)
     - DISH (:doc:`dish`)
     - PME (:doc:`pme`)
   * - Variable
     - amplitudes + active state per trajectory
     - amplitudes + active state per trajectory
     - populations :math:`P_{ia}`
   * - Transitions
     - fewest-switches flux :math:`\propto \operatorname{Re}(c_K^* c_J d_{KJ})`
     - collapse at decoherence events, probability :math:`|c_J|^2`
     - golden-rule rates :math:`\propto |d|^2` with Lorentzian width :math:`\hbar/\tau`
   * - Decoherence
     - damping :math:`c_J e^{-\Delta t/\tau_{KJ}}` per step
     - stochastic events, probability :math:`1 - e^{-\Delta t/\tau_{KJ}}`
     - assumed fast (coherences eliminated)
   * - Detailed balance
     - Boltzmann factor on upward hops (pair energies)
     - Boltzmann factor on upward collapses (pair energies)
     - Boltzmann factor on upward rates (orbital gaps)
   * - Limit
     - follows the coherent build-up of each hop; depends on how decoherence is added
     - golden-rule rates of the PME on average, with the same :math:`\tau`
     - golden rule; Markovian
   * - Gives
     - trajectories: dwell times, branching
     - trajectories: dwell times, branching
     - mean populations only
   * - Cost per nuclear step
     - :math:`N_{\mathrm{sub}}` products :math:`n^2 N_{\mathrm{traj}}` per channel
     - same as FSSH
     - one rate build and 20 products, :math:`< 1` s for :math:`10^6` pairs

When to use which
-----------------

* **Dense bands, fast dephasing** (hot-carrier cooling in dots of hundreds of atoms, dephasing times of
  10–30 fs): DISH and FSSH integrate the actual time dependence of the couplings and are the reference;
  they agree within a factor 1.3–1.7 on the CsPbX\ :sub:`3` dots. The PME agrees with them only if the
  couplings stay correlated over the dephasing time; when they decorrelate faster (dense valence bands,
  :doc:`pme`), the constant-coupling PME is 2× too fast, and the PME with the pair-resolved coupling
  correlation (``pme_tau: pairs_nac``) coincides with FSSH. FSSH without decoherence is 5–10× too fast
  (over-coherence).
* **Sparse levels** (the :math:`1P_e \to 1S_e` gap of a dot, a phonon bottleneck): transitions across
  gaps much larger than the phonon energies are rare and coherent effects matter; FSSH and DISH follow
  them explicitly, the PME only through the Lorentzian tails of its rates.
* **Surface traps**: trapping and de-trapping are rare, stochastic events with long dwell times; FSSH and
  DISH give their statistics, the PME a mean rate. Trapping needs trajectories much longer than the
  cooling (tens of picoseconds) and the trap states inside the active window.

For a new system, run DISH and FSSH, and the PME with ``pairs`` and ``pairs_nac``; check the coupling
correlation time in ``nac_correlation.npz`` against the dephasing times. If the PME variants agree with
DISH, the PME can be used for scans (pump energy, temperature, window).
