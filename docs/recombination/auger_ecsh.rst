Auger ecsh
==========

Part of :doc:`/recombination/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.namd.surface_hopping``
* Callable: ``qdex.namd.surface_hopping.run_namd_dynamics``
* CLI: ``--auger, --auger-channel``
* YAML: ``auger.run, auger.channel``

.. code-block:: python

   run_namd_dynamics(config)


8. Energy-Conserving Surface Hopping (ECSH) for Auger in NAMD
-------------------------------------------------------------


Physical Foundations: Electron-Phonon vs. Coulomb Transitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In traditional non-adiabatic molecular dynamics (NAMD), electronic transitions are governed solely by **Non-Adiabatic Couplings (NACs)**:

.. math::

   d_{ij} = \langle \Phi_i | \frac{\partial}{\partial t} | \Phi_j \rangle = \sum_I \mathbf{v}_I \cdot \mathbf{d}_{ij}^I

Because the nuclear derivative :math:`\nabla_I` is a **single-particle operator**, non-adiabatic transitions occur only between many-body Slater determinants that differ by **exactly one orbital** (:math:`\Delta N_{\mathrm{orb}} = 1`). These transitions describe **electron-phonon scattering**, in which electronic energy is dissipated into the classical vibrational degrees of freedom (phonons). To satisfy detailed balance, upward hops are scaled by the Boltzmann factor:

.. math::

   P_{i \to j}^{\mathrm{NAC}} = \max(0, g_{ij}^{\mathrm{NAC}}) \times \begin{cases} 1, & E_j \le E_i \\ \exp\left(-\frac{E_j - E_i}{k_B T}\right), & E_j > E_i \end{cases}

In contrast, **Auger processes** (such as biexciton annihilation :math:`XX \to X + \text{carrier}^*` or Auger-assisted carrier cooling) are mediated by the **two-particle Coulomb operator**:

.. math::

   \hat{V} = \frac{1}{2} \sum_{k \neq l} \frac{e^2}{\epsilon |\mathbf{r}_k - \mathbf{r}_l|}

Coulomb matrix elements :math:`V_{ij} = \langle \Phi_i | \hat{V} | \Phi_j \rangle` couple many-body states that differ by **two orbitals** (:math:`\Delta N_{\mathrm{orb}} = 2`).


The Fundamental Flaw in Previous Approaches (The "Energy Leak" Error)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In earlier NAMD approaches (e.g. Zhou, Lu, & Prezhdo, *Nano Lett.* 2021, 21, 756), Coulomb matrix elements :math:`V_{ij}` were added to the off-diagonal Hamiltonian alongside NACs, but all hops were subjected to the same Boltzmann scaling / velocity readjustment.

**The physical failure**: When an Auger hop occurs, the recombination energy is transferred **entirely within the quantum electronic subsystem** to the spectator carrier. **Zero energy is transferred to lattice vibrations during the Auger hop itself.** Scaling Coulomb hops with the Boltzmann factor caused the huge recombination energy (:math:`\sim 2.5 - 3.5\text{ eV}`) to unphysically leak into classical phonons, leading to artificially accelerated decay rates (up to :math:`2\times` too fast) and distorted kinetics.


The Gumber-Prezhdo ECSH Methodology
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To resolve this issue, Gumber & Prezhdo (*J. Chem. Theory Comput.* 2024, 20, 13, 5408–5418) developed **Energy-Conserving Surface Hopping (ECSH)**, which partitions transitions by orbital permutation count:

1. **Single-orbital hops** (:math:`\Delta N_{\mathrm{orb}} = 1`):
   - Mediated by NAC :math:`d_{ij}` (electron-phonon coupling).
   - Electronic energy is transferred to nuclear vibrations.
   - Upward hops are scaled by the Boltzmann factor :math:`\exp(-\Delta E / k_B T)`.

2. **Two-orbital hops** (:math:`\Delta N_{\mathrm{orb}} = 2`):
   - Mediated by Coulomb coupling :math:`V_{ij}` (Auger processes).
   - Energy remains strictly within the electronic subsystem.
   - **NO Boltzmann factor is applied.**
   - The hop is allowed **only if** the initial and final states are energetically resonant within the thermal window:

.. math::

   P_{i \to j}^{\mathrm{Coul}} = \max(0, g_{ij}^{\mathrm{Coul}}) \times \Theta\left(\Delta E_{\mathrm{window}} - |E_j - E_i|\right)

where :math:`\Delta E_{\mathrm{window}} \approx k_B T` (or :math:`2 k_B T`). Once the hop occurs to the hot-carrier state, the newly created hot carrier subsequently relaxes down its band via standard single-particle NACs, correctly dissipating its energy to lattice heat as phonons.


Relation with Static Trions (:math:`eeh` and :math:`hhe`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

How does dynamic ECSH NAMD relate to the static trion picture used in Section 1?

* **Identical Microscopic Matrix Elements**:
  When ECSH or a biexciton clock is requested, the hop rate is the frame-0 golden-rule :math:`k_{XX}` from the static Auger module (or an explicit ``tau_auger_ps``). Resonant finals inside the window share that rate. The energy stays in the electronic subsystem: no Boltzmann factor and no velocity rescaling (Gumber and Prezhdo, JCTC 2024). A cooling run does not load this path. Per-frame :math:`V_{\mathrm{dir}} - V_{\mathrm{exch}}` is not stored on the compacted trajectory. The static channels are:
  
  .. math::
  
     \langle \Phi_{XX} | \hat{V} | \Phi_{X^*} \rangle = \begin{cases} V^{eeh}(e'), & \text{spectator electron promoted} \\ V^{hhe}(h'), & \text{spectator hole promoted} \end{cases}

* **What moves along the trajectory**:
  The static golden rule uses one geometry and a Gaussian of width :math:`\sigma`. In the dynamics the diagonal-BSE energies move with the frame, and a two-body hop is accepted only inside :math:`\max(k_B T, \sigma)`. The Coulomb weight shared by those resonant finals is the frame-0 :math:`k_{XX}`, unless ``tau_auger_ps`` or ``k_auger_fs`` is set. Per-frame :math:`V_{\mathrm{dir}}` is not recomputed. If no final state lies in the window, the trajectory stays where it is.
* **Complete Kinetic Cascade**:
  The static calculation outputs only the single instantaneous rate :math:`\Gamma_{XX}`. ECSH NAMD simulates the **full sequence of events**: the biexciton lives on the :math:`XX` surface, undergoes an Auger hop to a hot single exciton, and then emits phonons through the single-particle NAC manifold as it cools to the band edge. In Transient Absorption, this reproduces the bi-exponential bleach recovery observed in experiments.


Resolving the Timescale Mismatch (:math:`1 - 10\text{ ps}` MD vs. :math:`100\text{ ps} - 10\text{ ns}` Auger)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Ab initio* DFT MD trajectories for nanocrystals typically span only :math:`1 - 10\text{ ps}` due to the high computational cost of DFT forces. However, Auger lifetimes :math:`\tau_{XX}` often range from **tens of picoseconds to several nanoseconds**. If :math:`\tau_{XX} = 500\text{ ps}`, the probability of observing an Auger hop within a :math:`2\text{ ps}` trajectory is only :math:`\sim 0.4\%`.

``QDEX`` provides three complementary solutions to bridge this timescale gap:

1. **Trajectory Looping under the Classical Path Approximation (CPA)**:
   Under the CPA, the ground-state MD trajectory represents an ergodic thermal bath at equilibrium. The electronic Hamiltonian and surface hopping can be propagated over times :math:`t \gg T_{\mathrm{MD}}` (e.g. 50–100 ps) by looping the precomputed MD frames:

   .. code-block:: bash

      qdex --config test_auger.yaml --namd-run --namd-trajectory-loops 20

   The nuclear sequence is repeated. At the seam the geometry jumps from the last frame back to the first with no connecting overlap, and electronic phases are not randomized. The loop extends the bath by repetition. It is not a longer molecular-dynamics trajectory. ``--namd-trajectory-loops`` does not turn Auger on. A cooling run leaves ``ecsh_auger`` false and does not set ``initial_state: biexciton``.

2. **Integrated Survival Probability & Initial Linear Decay Slope**:
   Even without trajectory looping, for :math:`t \ll \tau_{XX}`, the decay of the biexciton population is linear:

   .. math::

      P_{XX}(t) = \exp\left( -\int_0^t \Gamma_{XX}(t') dt' \right) \approx 1 - \Gamma_{XX}^{(0)} t

   The lifetime is accurately determined from the initial rate of population transfer:

   .. math::

      \tau_{XX} = \frac{1}{\langle \Gamma_{XX} \rangle} = \left( -\left. \frac{d P_{XX}(t)}{dt} \right|_{t=0} \right)^{-1}

   The ``qdex.auger.extract_auger_kinetics_from_trajectory(...)`` function performs this analysis automatically.

3. **Ultrafast Intraband Auger-Assisted Carrier Cooling**:
   Unlike interband biexciton recombination, **intraband Auger cooling** (where a hot electron cools from :math:`1P_e \to 1S_e` by kicking a valence hole deep into the valence band, bypassing the phonon bottleneck) occurs on an ultrafast timescale of :math:`50 - 500\text{ fs}`. This ultrafast process naturally fits well within standard :math:`1 - 2\text{ ps}` AIMD trajectories.

