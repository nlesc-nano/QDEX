Pme
===

Part of :doc:`/dynamics/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.namd.master_equation``
* Callable: ``qdex.namd.master_equation.propagate_pme_tensor``
* CLI: ``--namd``
* YAML: ``namd.engine, namd.dt_fs``

.. code-block:: python

   propagate_pme_tensor(P_mat, E_mat, d_occ, d_virt, dt_fs, temp_k=300.0, tau_dec_fs=None, n_substeps=20, eps_occ=None, eps_virt=None, n_atoms=775, k_loss=None, return_flux=False)


4. Theoretical Foundations of the Dynamical Engines
---------------------------------------------------


1. Derivation of the Pauli Master Equation (PME)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The exact quantum dynamics of a coupled electron-nuclear system is governed by the Liouville-von Neumann equation for the total density operator :math:`\hat{\rho}(t)`:

.. math::

   i\hbar \frac{\partial \hat{\rho}(t)}{\partial t} = \left[ \hat{H}(t), \hat{\rho}(t) \right]

Applying the **Nakajima-Zwanzig projection operator technique**, the density matrix is partitioned into diagonal populations :math:`\mathcal{P} \hat{\rho} = \sum_I \rho_{II} |I\rangle \langle I|` and off-diagonal electronic coherences :math:`\mathcal{Q} \hat{\rho} = \sum_{I \neq J} \rho_{IJ} |I\rangle \langle J|`.

In a condensed-phase environment containing many nuclear degrees of freedom, thermal fluctuations of the nuclear bath induce rapid random phase fluctuations that cause the off-diagonal coherences :math:`\rho_{IJ}(t)` to decay exponentially with characteristic decoherence time :math:`\tau_{\mathrm{dec}}`:

.. math::

   \rho_{IJ}(t) \approx \rho_{IJ}(0) \, \exp\left( -\frac{t}{\tau_{\mathrm{dec}}} \right) \exp\left( -i \frac{\Delta E_{IJ}}{\hbar} t \right)

Under the **Markovian approximation** (where the bath correlation time is much shorter than the population relaxation time), integrating out the rapidly decaying coherences yields the closed **Pauli Master Equation** for populations :math:`P_I(t) \equiv \rho_{II}(t)`:

.. math::

   \frac{d P_I(t)}{dt} = \sum_{J \neq I} \left[ k_{J \to I}(t) P_J(t) - k_{I \to J}(t) P_I(t) \right]

Fermi's Golden Rule with Lorentzian Broadening
""""""""""""""""""""""""""""""""""""""""""""""

The instantaneous rate constant :math:`k_{I \to J}` is obtained from second-order time-dependent perturbation theory (Fermi's Golden Rule). The standard energy-conserving Dirac delta function :math:`\delta(E_I - E_J)` is broadened by the finite electronic decoherence time :math:`\tau_{\mathrm{dec}}`, yielding a Lorentzian line shape:

.. math::

   k_{I \to J}(t) = 2 |d_{IJ}(t)|^2 \, \left[ \frac{\tau_{\mathrm{dec}}}{1 + \left( \frac{(E_J(t) - E_I(t)) \tau_{\mathrm{dec}}}{\hbar} \right)^2} \right] \times B_{IJ}(T)

where :math:`d_{IJ}(t) = \langle \psi_I | \frac{\partial}{\partial t} | \psi_J \rangle` is the non-adiabatic coupling, and :math:`B_{IJ}(T)` enforces thermodynamic **detailed balance** at lattice temperature :math:`T`:

.. math::

   B_{IJ}(T) = \begin{cases}
     1 & \text{for downward transitions } (E_J \le E_I) \\
     \exp\left( -\frac{E_J - E_I}{k_B T} \right) & \text{for upward thermal activation } (E_J > E_I)
   \end{cases}

Vectorized Tensor Decomposition for Diagonal BSE
""""""""""""""""""""""""""""""""""""""""""""""""

In a two-particle excitonic manifold with :math:`N_{\mathrm{occ}}` occupied orbitals and :math:`N_{\mathrm{virt}}` virtual orbitals, the total number of electron-hole pairs is :math:`N_{\mathrm{pairs}} = N_{\mathrm{occ}} \times N_{\mathrm{virt}}`. Constructing and multiplying an :math:`(N_{\mathrm{pairs}} \times N_{\mathrm{pairs}})` rate matrix scales as :math:`O(N_{\mathrm{pairs}}^2) = O(N_{\mathrm{occ}}^2 N_{\mathrm{virt}}^2)`. For :math:`N_{\mathrm{occ}} = N_{\mathrm{virt}} = 500`, this corresponds to an intractable :math:`250,000 \times 250,000` dense matrix (:math:`500\text{ GB}` of RAM).

Under the **Diagonal BSE** representation, the exciton state :math:`|ia\rangle` factorizes into an independent occupied hole state :math:`i` and an independent virtual electron state :math:`a`. An exciton relaxes either via an electron transition (:math:`a \to b`) with rate :math:`K_e(a \to b)` or a hole transition (:math:`i \to j`) with rate :math:`K_h(i \to j)`:

.. math::

   \frac{d P_{ia}(t)}{dt} = \sum_{b \neq a} \left[ K_e(b \to a) P_{ib} - K_e(a \to b) P_{ia} \right] + \sum_{j \neq i} \left[ K_h(j \to i) P_{ja} - K_h(i \to j) P_{ia} \right]

In matrix notation, this decomposes into an exact **BLAS Level-3 tensor product**:

.. math::

   \frac{\partial \mathbf{P}}{\partial t} = \left( \mathbf{P} \, \mathbf{K}_e - \mathbf{P} \operatorname{diag}(\mathbf{L}_e) \right) + \left( \mathbf{K}_h^\mathsf{T} \, \mathbf{P} - \operatorname{diag}(\mathbf{L}_h) \, \mathbf{P} \right)

where :math:`\mathbf{L}_e = \sum_b K_e(a \to b)` and :math:`\mathbf{L}_h = \sum_j K_h(i \to j)` are the total state loss vectors.

This breakthrough reduces the computational scaling from :math:`O(N_{\mathrm{pairs}}^2)` to :math:`O(N_{\mathrm{occ}}^2 + N_{\mathrm{virt}}^2)`. A million exciton configurations are propagated in **less than 0.2 seconds per nuclear time step**.
