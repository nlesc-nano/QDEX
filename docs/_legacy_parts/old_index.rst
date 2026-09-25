.. QDEX documentation master file

Welcome to QDEX's Documentation!
================================

**QDEX** (*pronounced "qiu-di-ex"*, :math:`/\text{kjuː-diː-ɛks}/`, formerly `miniBSE`) is a high-performance, lightweight post-DFT exciton solver and non-adiabatic molecular dynamics (NAMD) engine designed for calculating, analyzing, and propagating the excited states of molecules and semiconductor nanoclusters.

By combining the flexibility of a Python interface with a fast C++ backend powered by `Libint2` and `Eigen3`, `QDEX` provides researchers with an end-to-end framework from ground-state Kohn-Sham orbitals to carrier cooling dynamics and photoluminescence quantum yields.

Pedagogical Structure (From Easiest to Most Complex)
----------------------------------------------------

The documentation is organized systematically into progressive parts:

1. **Part 1: Ground-State Electronic Structure**: Kohn-Sham molecular orbitals from CP2K, Lowdin/Mulliken population analysis, Projected Density of States (PDOS), Inverse Participation Ratio (IPR), Crystal Orbital Overlap Population (COOP), supercell unfolding (Fuzzy Bands), and 3D Gaussian ``.cube`` orbital visualization.
2. **Part 2: Relativistic Spin-Orbit Coupling (SOC)**: Dirac reduction, :math:`Z^4` scaling in heavy elements (Pb, Bi, I), GTH separable pseudopotential projectors, two-component spinor Hamiltonian, and sparse CSR BLAS DGEMM acceleration.
3. **Part 3: Quasiparticle Corrections & The Scaled GW Models**: Quasiparticle concepts, Hedin's :math:`G_0W_0`, DFT band gap underestimation, the Unified Two-Body Engine (``2e-integrals`` vs. ``kernel``), Scaled GW hierarchy (``sgw-anchor``, ``sgw-dim``, ``sgw-resta``, ``evGW``, ``qsGW``), Approach B (anchor-free microscopic wavefunction asymmetry :math:`f_H^{\mathrm{micro}}, f_L^{\mathrm{micro}}`), and spatial asymptotics converging to the bulk crystal limit (:math:`R_{\mathrm{QD}} \to \infty`).
4. **Part 4: Optical Excitations & Four Excited-State Frameworks**: Neutral two-particle excitation physics, Tamm-Dancoff Approximation (TDA), shared two-body engine (:math:`\Xi = v - W`), four excitation frameworks (``independent_dft``, ``independent_qp``, ``diagonal_bse``, ``bse``), relativistic 2-component spinor BSE, spatial asymptotics with Wannier-Mott bulk exciton convergence, and Davidson iterative eigensolver.
5. **Part 5: Excited-State Wavefunction Analysis**: Rigorous Plasser-Dreuw real-space descriptors (:math:`d_{eh}, R_{eh}, d_{\mathrm{CT}}, \sigma_h, \sigma_e`), Natural Transition Orbitals (NTOs), and volumetric 3D visualization.
6. **Part 6: Non-Adiabatic Dynamics & Photoluminescence (NAMD)**: Cross-frame overlaps and non-adiabatic couplings, gauge phase alignment, Hungarian crossing tracking, the tensorized Pauli Master Equation (PME), CPA-FSSH surface hopping, *ab initio* cumulant decoherence, Einstein radiative decay, and Photoluminescence Quantum Yield (PLQY).
7. **Part 7: Step-by-Step Tutorials**: Practical, end-to-end computational workflows mirroring the exact progression from Parts 1 through 6 on a model perovskite quantum dot.
8. **Part 8: Multi-Carrier Auger Recombination**: Non-radiative Auger decay in quantum dots, breakdown of translational symmetry and :math:`1/V^2` volume scaling, negative trion (:math:`eeh`) and positive trion (:math:`hhe`) channels, biexciton lifetimes (:math:`\tau_{XX}`), atom-centered Resta-screened Coulomb contractions, and static vs. trajectory-averaged NAMD dynamics.

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   getting_started/installation
   getting_started/quickstart
   getting_started/configuration

.. toctree::
   :maxdepth: 2
   :caption: Part 1: Ground-State Electronic Structure

   part1_ground_state/index

.. toctree::
   :maxdepth: 2
   :caption: Part 2: Spin-Orbit Coupling (SOC)

   part2_soc/index

.. toctree::
   :maxdepth: 2
   :caption: Part 3: Quasiparticle (GW) Corrections

   part3_gw_scissor/index

.. toctree::
   :maxdepth: 2
   :caption: Part 4: Optical Excitations & BSE

   part4_excited_states/index

.. toctree::
   :maxdepth: 2
   :caption: Part 5: Exciton Wavefunction Analysis

   part5_exciton_analysis/index

.. toctree::
   :maxdepth: 2
   :caption: Part 6: Non-Adiabatic Dynamics (NAMD)

   part6_namd/index

.. toctree::
   :maxdepth: 2
   :caption: Part 7: Step-by-Step Tutorials

   part7_examples/index

.. toctree::
   :maxdepth: 2
   :caption: Part 8: Multi-Carrier Auger Recombination

   part8_auger/index

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
