Introduction: why nonadiabatic dynamics
=======================================

Part of :doc:`/dynamics/index`.

.. figure:: /_static/figures/namd_workflow.svg
   :width: 100%
   :alt: namd workflow

   Non-adiabatic dynamics pipeline and its observables.

The problem
-----------

A photon above the gap creates an electron and a hole with excess energy. In a quantum dot that energy
is lost to the lattice within femtoseconds to picoseconds (**hot-carrier cooling**), carriers may be
captured by surface states (**trapping**), and the band-edge exciton finally recombines, radiatively or
not. These processes decide how bright a dot is, how fast it responds, and whether hot carriers can be
extracted before they cool.

None of them is visible in a static excited-state calculation (:doc:`/excitons/index`): they are
transitions *between* electronic states, driven by the motion of the nuclei. An electronic state that is
an eigenstate at one geometry is a superposition of eigenstates at the next; the rate at which the
nuclear motion mixes the states is the **non-adiabatic coupling**, and the energy released in each
transition goes into the vibrations. Following this requires solving the electronic time-dependent
Schrödinger equation along the nuclear motion, beyond the Born–Oppenheimer approximation, which assumes
the electrons always stay in one state. This is **non-adiabatic molecular dynamics (NAMD)**.

What QDEX computes
------------------

* carrier cooling curves: excess energy of the electron, of the hole and of the exciton above the band
  edges, cooling times;
* populations of the exciton states in time, arrival at the band edge, trapping in surface states;
* the vibrations that drive the transitions (spectral densities), dephasing times;
* with recombination enabled, radiative and non-radiative decay and the photoluminescence quantum yield
  (:doc:`/recombination/index`); transient-absorption maps (:doc:`/spectroscopy/index`).

The equations
-------------

**Electrons and nuclei.** The nuclei follow classical trajectories :math:`\mathbf{R}(t)`; the electrons
obey the time-dependent Schrödinger equation in the field of the moving nuclei,

.. math::

   i\hbar\, \frac{\partial}{\partial t} |\Psi(t)\rangle = \hat H_{\mathrm{el}}\big(\mathbf{R}(t)\big)\, |\Psi(t)\rangle .

**Adiabatic basis.** At every geometry the excited states :math:`|\psi_I(\mathbf{R})\rangle` with energies
:math:`E_I(\mathbf{R})` diagonalise :math:`\hat H_{\mathrm{el}}`. Expanding
:math:`|\Psi(t)\rangle = \sum_I c_I(t)\, |\psi_I(\mathbf{R}(t))\rangle` and projecting gives

.. math::

   i\hbar\, \dot c_I = E_I\, c_I - i\hbar \sum_J d_{IJ}\, c_J ,
   \qquad
   d_{IJ}(t) = \Big\langle \psi_I \Big| \frac{\partial \psi_J}{\partial t} \Big\rangle
   = \dot{\mathbf{R}} \cdot \langle \psi_I | \nabla_{\mathbf{R}}\, \psi_J \rangle .

:math:`\langle \psi_I | \nabla_{\mathbf{R}} \psi_J \rangle` is the non-adiabatic coupling *vector*; its
projection on the nuclear velocity, :math:`d_{IJ}`, is the time-derivative (scalar) coupling, which is all
the electrons need and which QDEX obtains directly from the overlaps of the states at consecutive frames
(:doc:`states_couplings`). :math:`d` is anti-Hermitian, so the amplitudes stay normalised.
In matrix form, with :math:`\mathbf{H}_{\mathrm{eff}} = \operatorname{diag}\mathbf{E} - i\hbar\,\mathbf{d}`
(Hermitian),

.. math::

   i\hbar\, \dot{\mathbf{c}} = \mathbf{H}_{\mathrm{eff}}(t)\, \mathbf{c} .

**What is still missing: the feedback on the nuclei and decoherence.** Solved as it stands (Ehrenfest
dynamics), the electrons would move the nuclei with the average force of a superposition, and the
superposition would never collapse: the carriers would not relax to a thermal distribution. In a real dot
the nuclear wavepackets that accompany different electronic states separate within femtoseconds and the
electronic superposition loses its phase (**decoherence**, :doc:`decoherence`). The three methods of QDEX
add this in different ways:

* **surface hopping** (:doc:`fssh`, :doc:`dish`): each trajectory occupies one *active* state; it jumps
  between states stochastically, with probabilities from the amplitudes :math:`c_I` and the couplings, and
  decoherence collapses the amplitudes;
* **master equation** (:doc:`pme`): coherences are assumed to decay fast; only populations
  :math:`P_I(t)` are propagated, with golden-rule rates from :math:`|d_{IJ}|^2`, the energy gaps and the
  dephasing times.

**Classical path approximation (CPA).** Exciting one electron–hole pair in a dot of hundreds of atoms
changes the forces on the nuclei very little, so the nuclear trajectory is taken from a single
ground-state MD run and reused for all electronic trajectories. The energy given to the vibrations in a
downward transition is not removed from the trajectory, and upward transitions, which in full dynamics need
kinetic energy from the nuclei, are weighted by the Boltzmann factor
:math:`B_{IJ} = \min\!\big(1, e^{-(E_J - E_I)/k_B T}\big)` (detailed balance at the temperature of the MD).
The heavy part, one DFT calculation per frame, is then done once, and any number of NAMD runs (methods,
pump energies, temperatures) are post-processing.

**Excitons and carriers.** The states :math:`|\psi_I\rangle` are excitons :math:`|ia\rangle`: a hole in
orbital :math:`i` and an electron in orbital :math:`a`, with the pair energy of the diagonal screened BSE
(:doc:`states_couplings`). The couplings between them are one-electron couplings of the electron or of
the hole, so the exciton dynamics separates into an **electron channel** and a **hole channel**, which is
what makes manifolds of :math:`10^6` exciton states tractable.

The workflow
------------

.. code-block:: text

   AIMD trajectory: frame_000001 ... frame_N (geometry + MOs of every frame)
         |
         v
   1. qdex --namd-precompute    excitons of every frame (diagonal sBSE), overlaps of consecutive frames,
                                Hungarian tracking of crossings, phases, Kramers pairs  -> step_*.npz
   2. qdex --namd-decoherence   dephasing time of every pair of states (cumulant)       -> decoherence_times.npz
   3. qdex --namd-nac           couplings d = log(U)/dt of every step                   -> nac_*.npz
         |
         v
   4. qdex --namd-run           FSSH, DISH or PME from a pump energy; one or many origins along the
                                trajectory  -> cooling curves, populations, lifetimes, dashboards
   5. analysis                  cooling times, spectral densities, gap statistics (:doc:`analysis`)

Steps 1–3 depend only on the trajectory and the excited-state model; step 4 is repeated for every method,
pump energy or temperature. :doc:`running` lists the commands and keys, :doc:`choosing` compares the
methods, and :doc:`gxtb` describes trajectories with g-xTB orbitals.

``--namd-run`` is a cooling calculation; Auger recombination, energy-conserving two-body hops and a
biexciton initial state are off unless ``--auger``, ``--namd-ecsh-auger`` or ``--namd-biexciton`` is set.
