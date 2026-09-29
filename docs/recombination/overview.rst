Overview
========

Part of :doc:`/recombination/index`.

.. rubric:: QDEX implementation

Implementation entry point:

* Module: ``qdex.auger``
* Callable: ``qdex.auger.calculate_auger_rates``
* CLI: ``--auger, --auger-channel``
* YAML: ``auger.run, auger.channel``

.. code-block:: python

   calculate_auger_rates(C: np.ndarray, eps: np.ndarray, S: np.ndarray, atom_ao_ranges: List[Tuple[int, int]], coords: np.ndarray, atom_symbols: List[str], homo_idx: int, W_resta: Optional[np.ndarray]=None, material_name: Optional[str]='DEFAULT', eps_out: float=2.0, sigma_ev: float=0.05, broadening_mode: str='gaussian', channel: str='all', n_initial_elec: int=1, n_initial_hole: int=1, e_search_sigma_factor: float=4.0, spinor: bool=False, U_spinor_alpha: Optional[np.ndarray]=None, U_spinor_beta: Optional[np.ndarray]=None, lambda_reorg_ev: Optional[float]=None, temperature_k: float=300.0, eps_eff: Optional[float]=None, verbose: bool=True)


8. Radiative & Non-Radiative Recombination Mechanisms
-----------------------------------------------------

Once photoexcited hot carriers have cooled to the band edge (forming the lowest :math:`1S` exciton state), they recombine to the electronic ground state :math:`|S_0\rangle` on longer timescales (nanoseconds). In ``QDEX``, carrier recombination is coupled directly to the population dynamics:

.. math::

   \frac{d P_I(t)}{dt} = \sum_{J \neq I} \left[ k_{J \to I} P_J - k_{I \to J} P_I \right] - \left( k_{I \to 0}^{\mathrm{rad}} + k_{I \to 0}^{\mathrm{nr}} \right) P_I

where :math:`k_{I \to 0}^{\mathrm{rad}}` and :math:`k_{I \to 0}^{\mathrm{nr}}` represent the radiative and non-radiative recombination rates from excited state :math:`|I\rangle` to the ground state.


1. Microscopic Origin: From DFT and BSE to Recombination
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To avoid empirical approximations, all transition energies and optical matrix elements originate directly from the ground-state DFT and excited-state frameworks:

1. **Ground-State DFT (CP2K)**:
   Self-consistent field calculations yield the Kohn-Sham molecular orbitals :math:`\psi_i(\mathbf{r}) = \sum_\mu C_{\mu i} \chi_\mu(\mathbf{r})` and orbital eigenvalues :math:`\varepsilon_i`. The atomic-orbital transition dipole matrices are computed analytically via Libint:

   .. math::

      \boldsymbol{\mu}_{\mu\nu}^{\mathrm{AO}} = \langle \chi_\mu | e \mathbf{r} | \chi_\nu \rangle

2. **Spin-Orbit Coupling (SOC)**:
   When SOC is active, spatial orbitals are projected into two-component spinors:

   .. math::

      \Psi_n^{\mathrm{SP}}(\mathbf{r}) = U_{n\alpha} \psi_\alpha(\mathbf{r}) + U_{n\beta} \psi_\beta(\mathbf{r})

   The transition dipole matrix is rotated into the spinor basis:

   .. math::

      \boldsymbol{\mu}_{ia}^{\mathrm{SP}} = \mathbf{U}_{\mathrm{occ},\alpha}^\dagger \mathbf{M} \mathbf{U}_{\mathrm{virt},\alpha} + \mathbf{U}_{\mathrm{occ},\beta}^\dagger \mathbf{M} \mathbf{U}_{\mathrm{virt},\beta}

   where :math:`\mathbf{M} = \mathbf{C}_{\mathrm{act}}^T \boldsymbol{\mu}^{\mathrm{AO}} \mathbf{C}_{\mathrm{act}}`.

3. **Quasiparticle Corrections & Excitation Energies**:
   Applying the GW scissor shift :math:`\Delta_{\mathrm{scissor}}` yields quasiparticle eigenvalues :math:`\tilde{\varepsilon}_i = \varepsilon_i - \Delta_{\mathrm{scissor}} f_i^{\mathrm{HOMO}}` and :math:`\tilde{\varepsilon}_a = \varepsilon_a + \Delta_{\mathrm{scissor}} f_a^{\mathrm{LUMO}}`.

   * **In Diagonal BSE (``diagonal_bse``)**:
     The transition energy of an electron-hole pair :math:`|i \to a\rangle` incorporates the direct screened Coulomb attraction:

     .. math::

        E_{ia} = \tilde{\varepsilon}_a - \tilde{\varepsilon}_i - K_d(ia, ia)

     where :math:`K_d(ia, ia) = \sum_{AB} q_{ii}^A W_{AB} q_{aa}^B` is the Ohno-Klopman or Resta-MNOK screened kernel.
   * **In Full BSE / sTDA (``bse``, ``stda``)**:
     Exciton states are coherent superpositions of single-particle transitions:

     .. math::

        |\Psi_I\rangle = \sum_{ia} c_{ia}^I |i \to a\rangle

     with total transition dipole :math:`\boldsymbol{\mu}_I = \sum_{ia} c_{ia}^I \boldsymbol{\mu}_{ia}^{\mathrm{SP}}`.

4. **Oscillator Strengths**:
   The dimensionless oscillator strength for transition :math:`|I\rangle \to |S_0\rangle` is evaluated as:

   .. math::

      f = \begin{cases}
      \dfrac{4}{3}\, E_{\mathrm{Ha}} \, |\boldsymbol{\mu}|^2 & \text{closed-shell spatial orbital} \\
      \dfrac{2}{3}\, E_{\mathrm{Ha}} \, |\boldsymbol{\mu}|^2 & \text{spinor}
      \end{cases}

   with the excitation energy in Hartree and the dipole in :math:`ea_0`. This is the atomic-unit reduction of :math:`f = (2/3)(m_e/\hbar^2) E |\mu|^2` per spin-orbital, doubled for a singlet that uses a spatial orbital. Older precomputes that stored :math:`(2/3) E_{\mathrm{eV}} |\mu|^2` are rescaled when the dynamics read them.
