# Documentation migration inventory

13 RST files; 298 disjoint units; 5688 lines; every byte accounted for.

| Source | Lines | Current heading | Proposed destination |
|---|---:|---|---|
| `docs/api/index.rst` | 1–5 | API Reference | `docs/reference/api.rst` |
| `docs/api/index.rst` | 6–13 | qdex.solver | `docs/reference/api.rst` |
| `docs/api/index.rst` | 14–21 | qdex.cli | `docs/reference/api.rst` |
| `docs/api/index.rst` | 22–29 | qdex.pdos_coop | `docs/reference/api.rst` |
| `docs/api/index.rst` | 30–37 | qdex.fuzzy_bands | `docs/reference/api.rst` |
| `docs/api/index.rst` | 38–45 | qdex.soc_utils | `docs/reference/api.rst` |
| `docs/api/index.rst` | 46–53 | qdex.hardness | `docs/reference/api.rst` |
| `docs/api/index.rst` | 54–61 | qdex.exciton_analysis | `docs/reference/api.rst` |
| `docs/api/index.rst` | 62–69 | qdex.nto | `docs/reference/api.rst` |
| `docs/api/index.rst` | 70–77 | qdex.exciton_cube | `docs/reference/api.rst` |
| `docs/api/index.rst` | 78–85 | qdex.namd.precompute | `docs/reference/api.rst` |
| `docs/api/index.rst` | 86–93 | qdex.namd.master_equation | `docs/reference/api.rst` |
| `docs/api/index.rst` | 94–101 | qdex.namd.surface_hopping | `docs/reference/api.rst` |
| `docs/api/index.rst` | 102–109 | qdex.namd.analysis | `docs/reference/api.rst` |
| `docs/api/index.rst` | 110–117 | qdex.namd.transient_absorption | `docs/reference/api.rst` |
| `docs/api/index.rst` | 118–125 | qdex.io_utils | `docs/reference/api.rst` |
| `docs/api/index.rst` | 126–133 | qdex.device_utils | `docs/reference/api.rst` |
| `docs/api/index.rst` | 134–141 | qdex.profiler | `docs/reference/api.rst` |
| `docs/api/index.rst` | 142–148 | qdex.auger | `docs/reference/api.rst` |
| `docs/getting_started/configuration.rst` | 1–5 | Configuration Reference (YAML) | `docs/start/configuration.rst` |
| `docs/getting_started/configuration.rst` | 6–90 | Complete Example Configuration | `docs/start/configuration.rst` |
| `docs/getting_started/configuration.rst` | 91–93 | Detailed Keyword Reference | `docs/start/configuration.rst` |
| `docs/getting_started/configuration.rst` | 94–105 | system | `docs/start/configuration.rst` |
| `docs/getting_started/configuration.rst` | 106–138 | physics | `docs/start/configuration.rst` |
| `docs/getting_started/configuration.rst` | 139–151 | fuzzy | `docs/start/configuration.rst` |
| `docs/getting_started/configuration.rst` | 152–158 | cube | `docs/start/configuration.rst` |
| `docs/getting_started/configuration.rst` | 159–166 | analysis | `docs/start/configuration.rst` |
| `docs/getting_started/configuration.rst` | 167–174 | auger | `docs/start/configuration.rst` |
| `docs/getting_started/configuration.rst` | 175–185 | namd | `docs/start/configuration.rst` |
| `docs/getting_started/installation.rst` | 1–5 | Installation Guide | `docs/start/installation.rst` |
| `docs/getting_started/installation.rst` | 6–26 | Method 1: Conda / Micromamba (Recommended) | `docs/start/installation.rst` |
| `docs/getting_started/installation.rst` | 27–45 | Method 2: Manual Pip Installation | `docs/start/installation.rst` |
| `docs/getting_started/installation.rst` | 46–61 | Verifying the Installation | `docs/start/installation.rst` |
| `docs/getting_started/quickstart.rst` | 1–5 | Quickstart Guide | `docs/start/quickstart.rst` |
| `docs/getting_started/quickstart.rst` | 6–32 | 1. Static BSE Calculation via CLI | `docs/start/quickstart.rst` |
| `docs/getting_started/quickstart.rst` | 33–64 | 2. YAML Configuration Workflow | `docs/start/quickstart.rst` |
| `docs/getting_started/quickstart.rst` | 65–69 | 3. NAMD Carrier Cooling Workflow | `docs/start/quickstart.rst` |
| `docs/getting_started/quickstart.rst` | 70–82 | Step 1: Precompute overlaps and exciton states | `docs/start/quickstart.rst` |
| `docs/getting_started/quickstart.rst` | 83–91 | Step 2: Run Dynamics & Carrier Cooling | `docs/start/quickstart.rst` |
| `docs/index.rst` | 1–2 | Preamble | `docs/index.rst` |
| `docs/index.rst` | 3–9 | Welcome to QDEX's Documentation! | `docs/index.rst` |
| `docs/index.rst` | 10–85 | Pedagogical Structure (From Easiest to Most Complex) | `docs/index.rst` |
| `docs/index.rst` | 86–91 | Indices and tables | `docs/index.rst` |
| `docs/part1_ground_state/index.rst` | 1–9 | Part 1: Ground-State Electronic Structure Analysis | `docs/electronic_structure/orbitals.rst` |
| `docs/part1_ground_state/index.rst` | 10–31 | 1. Molecular Orbitals and the Atomic Orbital Basis | `docs/electronic_structure/orbitals.rst` |
| `docs/part1_ground_state/index.rst` | 32–50 | MO Orthonormality and Diagnostics | `docs/electronic_structure/orbitals.rst` |
| `docs/part1_ground_state/index.rst` | 51–55 | 2. Projected Density of States (PDOS) | `docs/electronic_structure/populations_pdos.rst` |
| `docs/part1_ground_state/index.rst` | 56–82 | Mulliken Population Weights | `docs/electronic_structure/populations_pdos.rst` |
| `docs/part1_ground_state/index.rst` | 83–93 | Energy Convolution & Smearing | `docs/electronic_structure/populations_pdos.rst` |
| `docs/part1_ground_state/index.rst` | 94–106 | Surface vs. Core Spatial Partitioning | `docs/electronic_structure/populations_pdos.rst` |
| `docs/part1_ground_state/index.rst` | 107–111 | 3. Inverse Participation Ratio (IPR) | `docs/electronic_structure/localization.rst` |
| `docs/part1_ground_state/index.rst` | 112–120 | Mathematical Formulation | `docs/electronic_structure/localization.rst` |
| `docs/part1_ground_state/index.rst` | 121–139 | Physical Interpretation | `docs/electronic_structure/localization.rst` |
| `docs/part1_ground_state/index.rst` | 140–144 | 4. Crystal Orbital Overlap Population (COOP) | `docs/electronic_structure/bonding_coop.rst` |
| `docs/part1_ground_state/index.rst` | 145–159 | Mathematical Definition | `docs/electronic_structure/bonding_coop.rst` |
| `docs/part1_ground_state/index.rst` | 160–170 | Bonding vs. Antibonding Character | `docs/electronic_structure/bonding_coop.rst` |
| `docs/part1_ground_state/index.rst` | 171–177 | 5. Fuzzy Band Structure (Supercell Unfolding) | `docs/electronic_structure/fuzzy_bands.rst` |
| `docs/part1_ground_state/index.rst` | 178–194 | Plane-Wave Fourier Projection | `docs/electronic_structure/fuzzy_bands.rst` |
| `docs/part1_ground_state/index.rst` | 195–208 | Brillouin Zone Folding | `docs/electronic_structure/fuzzy_bands.rst` |
| `docs/part1_ground_state/index.rst` | 209–218 | Automated High-Symmetry Paths & PCA Alignment | `docs/electronic_structure/fuzzy_bands.rst` |
| `docs/part1_ground_state/index.rst` | 219–223 | 6. Visualizing Orbitals via 3D Gaussian .cube Files | `docs/electronic_structure/orbital_cubes.rst` |
| `docs/part1_ground_state/index.rst` | 224–234 | Grid Discretization | `docs/electronic_structure/orbital_cubes.rst` |
| `docs/part1_ground_state/index.rst` | 235–251 | C++ Libint Acceleration | `docs/electronic_structure/orbital_cubes.rst` |
| `docs/part1_ground_state/index.rst` | 252–256 | 7. CLI Flags & YAML Configuration Reference | `docs/electronic_structure/configuration.rst` |
| `docs/part1_ground_state/index.rst` | 257–309 | Command-Line Arguments | `docs/electronic_structure/configuration.rst` |
| `docs/part1_ground_state/index.rst` | 310–336 | YAML Configuration Example | `docs/electronic_structure/configuration.rst` |
| `docs/part2_soc/index.rst` | 1–9 | Part 2: Relativistic Spin-Orbit Coupling (SOC) | `docs/relativity/foundations.rst` |
| `docs/part2_soc/index.rst` | 10–26 | 1. Physical Origin of Spin-Orbit Coupling | `docs/relativity/foundations.rst` |
| `docs/part2_soc/index.rst` | 27–39 | The :math:`Z^4` Scaling Law | `docs/relativity/foundations.rst` |
| `docs/part2_soc/index.rst` | 40–44 | 2. Impact of SOC on Semiconductor Nanocrystals | `docs/relativity/material_effects.rst` |
| `docs/part2_soc/index.rst` | 45–62 | Band Inversion & Giant Gap Contraction | `docs/relativity/material_effects.rst` |
| `docs/part2_soc/index.rst` | 63–75 | Rashba-Dresselhaus Splitting | `docs/relativity/material_effects.rst` |
| `docs/part2_soc/index.rst` | 76–100 | 3. Survey of Relativistic Methods & Rationale for GTH Pseudopotentials | `docs/relativity/pseudopotentials.rst` |
| `docs/part2_soc/index.rst` | 101–110 | Why QDEX Chooses GTH SOC Pseudopotentials | `docs/relativity/pseudopotentials.rst` |
| `docs/part2_soc/index.rst` | 111–131 | 4. Mathematical Formulation: Separable GTH Pseudopotential | `docs/relativity/spinor_hamiltonian.rst` |
| `docs/part2_soc/index.rst` | 132–150 | Angular Momentum Matrices in Spherical Harmonics | `docs/relativity/spinor_hamiltonian.rst` |
| `docs/part2_soc/index.rst` | 151–186 | Two-Component Spinor Hamiltonian Structure | `docs/relativity/spinor_hamiltonian.rst` |
| `docs/part2_soc/index.rst` | 187–204 | 5. High-Performance Sparse Assembly & DGEMM Optimization | `docs/relativity/assembly.rst` |
| `docs/part2_soc/index.rst` | 205–226 | 6. Spinor Representation & Unrestricted Kohn-Sham (UKS) | `docs/relativity/unrestricted.rst` |
| `docs/part2_soc/index.rst` | 227–242 | UKS Spin-Preserving Framework | `docs/relativity/unrestricted.rst` |
| `docs/part2_soc/index.rst` | 243–245 | 7. CLI Flags & YAML Configuration Reference | `docs/relativity/configuration.rst` |
| `docs/part2_soc/index.rst` | 246–271 | Command-Line Arguments | `docs/relativity/configuration.rst` |
| `docs/part2_soc/index.rst` | 272–286 | YAML Configuration Example | `docs/relativity/configuration.rst` |
| `docs/part3_gw_scissor/index.rst` | 1–9 | Part 3: Quasiparticle Corrections & The Scaled GW Models | `docs/quasiparticles/foundations.rst` |
| `docs/part3_gw_scissor/index.rst` | 10–23 | 1. Quasiparticle Theory and Hedin's GW Approximation | `docs/quasiparticles/foundations.rst` |
| `docs/part3_gw_scissor/index.rst` | 24–47 | Hedin's Equations & The :math:`G_0W_0` Approximation | `docs/quasiparticles/foundations.rst` |
| `docs/part3_gw_scissor/index.rst` | 48–80 | 2. Unified Two-Body Interaction Engine (``2e-integrals`` vs. ``kernel``) | `docs/interactions/architecture.rst` |
| `docs/part3_gw_scissor/index.rst` | 81–109 | Shared Two-Body Operator Architecture | `docs/interactions/architecture.rst` |
| `docs/part3_gw_scissor/index.rst` | 110–134 | Role in Quasiparticle Theory vs. Excited States | `docs/interactions/architecture.rst` |
| `docs/part3_gw_scissor/index.rst` | 135–153 | 3. The DFT Band Gap Problem | `docs/quasiparticles/dft_reference.rst` |
| `docs/part3_gw_scissor/index.rst` | 154–164 | 4. The Nanocrystal Scaling Bottleneck | `docs/quasiparticles/cost.rst` |
| `docs/part3_gw_scissor/index.rst` | 165–188 | 5. Hierarchy of Quasiparticle Models in QDEX | `docs/quasiparticles/model_selection.rst` |
| `docs/part3_gw_scissor/index.rst` | 189–207 | 6. Avenue 1: Two-Anchor Scaled GW (``sgw-anchor``) | `docs/quasiparticles/anchor.rst` |
| `docs/part3_gw_scissor/index.rst` | 208–241 | The Confinement Interpolation Formula | `docs/quasiparticles/anchor.rst` |
| `docs/part3_gw_scissor/index.rst` | 242–260 | Implementation in QDEX (``sgw-anchor``) | `docs/quasiparticles/anchor.rst` |
| `docs/part3_gw_scissor/index.rst` | 261–265 | 7. Avenue 2: Microscopic Dielectric Shift Model (:math:`\Delta W`) | `docs/quasiparticles/delta_w.rst` |
| `docs/part3_gw_scissor/index.rst` | 266–280 | The Physical Rationale: Cancellation of :math:`v_{xc}` | `docs/quasiparticles/delta_w.rst` |
| `docs/part3_gw_scissor/index.rst` | 281–317 | Screened COHSEX Operator | `docs/quasiparticles/delta_w.rst` |
| `docs/part3_gw_scissor/index.rst` | 318–355 | Screening Formulations for :math:`W^{\mathrm{QD}}` | `docs/interactions/qp_screening.rst` |
| `docs/part3_gw_scissor/index.rst` | 356–387 | Spatial Asymptotics of the Dielectric Kernel: Why Screening Fits Nanocrystals | `docs/interactions/asymptotics.rst` |
| `docs/part3_gw_scissor/index.rst` | 388–429 | Approach B: Microscopic Wavefunction Asymmetry (:math:`f_H^{\mathrm{micro}}, f_L^{\mathrm{micro}}`) | `docs/quasiparticles/edge_partition.rst` |
| `docs/part3_gw_scissor/index.rst` | 430–468 | Implementation in QDEX (``sgw-dim``, ``sgw-resta``, ``sgw``) | `docs/quasiparticles/delta_w_implementation.rst` |
| `docs/part3_gw_scissor/index.rst` | 469–471 | 8. Avenue 3: Dynamic Renormalization & Self-Consistency | `docs/quasiparticles/dynamic_z.rst` |
| `docs/part3_gw_scissor/index.rst` | 472–488 | Dynamic Renormalization Factor :math:`Z_p` (``--dynamic_z``) | `docs/quasiparticles/dynamic_z.rst` |
| `docs/part3_gw_scissor/index.rst` | 489–498 | Implementation in QDEX (``--dynamic_z``) | `docs/quasiparticles/dynamic_z.rst` |
| `docs/part3_gw_scissor/index.rst` | 499–509 | Eigenvalue Self-Consistent GW (``evgw-dim``, ``evgw-resta``) | `docs/quasiparticles/gap_iteration.rst` |
| `docs/part3_gw_scissor/index.rst` | 510–519 | Implementation in QDEX (``evgw-dim``, ``evgw-resta``) | `docs/quasiparticles/gap_iteration.rst` |
| `docs/part3_gw_scissor/index.rst` | 520–550 | Quasiparticle Self-Consistent GW (``qsgw-dim``, ``qsgw-resta``) | `docs/quasiparticles/orbital_iteration.rst` |
| `docs/part3_gw_scissor/index.rst` | 551–580 | Implementation in QDEX (``qsgw-dim``, ``qsgw-resta``, ``--update_orbitals``) | `docs/quasiparticles/orbital_iteration.rst` |
| `docs/part3_gw_scissor/index.rst` | 581–585 | 9. Environmental Dielectric Polarization | `docs/interactions/environment.rst` |
| `docs/part3_gw_scissor/index.rst` | 586–600 | Classical Image Charge Solvation | `docs/interactions/environment.rst` |
| `docs/part3_gw_scissor/index.rst` | 601–618 | 10. Rigid Scissor Operator & Frontier Splitting Strategies (Approach A vs. Approach B) | `docs/quasiparticles/alignment.rst` |
| `docs/part3_gw_scissor/index.rst` | 619–637 | Approach A: Database Monomer Anchor Frontier Splitting (Used by ``sgw-anchor``) | `docs/quasiparticles/alignment.rst` |
| `docs/part3_gw_scissor/index.rst` | 638–668 | Approach B: Microscopic Wavefunction Asymmetry (Used by ``sgw-dim``, ``sgw-resta``, ``evgw``, ``qsgw``) | `docs/quasiparticles/alignment.rst` |
| `docs/part3_gw_scissor/index.rst` | 669–699 | Implementation in QDEX (Frontier Alignment & CLI) | `docs/quasiparticles/alignment.rst` |
| `docs/part3_gw_scissor/index.rst` | 700–704 | 11. Material Database Reference | `docs/reference/materials.rst` |
| `docs/part3_gw_scissor/index.rst` | 705–782 | Perovskites | `docs/reference/materials.rst` |
| `docs/part3_gw_scissor/index.rst` | 783–890 | II–VI Semiconductors | `docs/reference/materials.rst` |
| `docs/part3_gw_scissor/index.rst` | 891–1000 | III–V Semiconductors | `docs/reference/materials.rst` |
| `docs/part3_gw_scissor/index.rst` | 1001–1062 | 12. Comprehensive Quasiparticle Options Matrix | `docs/quasiparticles/model_selection.rst` |
| `docs/part3_gw_scissor/index.rst` | 1063–1065 | 13. Recommended Workflow Presets | `docs/workflows/qp_presets.rst` |
| `docs/part3_gw_scissor/index.rst` | 1066–1080 | Preset 1: Ultra-Fast NAMD Trajectory Dynamics (``sgw-anchor``) | `docs/workflows/qp_presets.rst` |
| `docs/part3_gw_scissor/index.rst` | 1081–1096 | Preset 2: First-Principles Nanocrystal Screening (``sgw-dim``) | `docs/workflows/qp_presets.rst` |
| `docs/part3_gw_scissor/index.rst` | 1097–1111 | Preset 3: Eigenvalue Self-Consistency (``evgw-dim``) | `docs/workflows/qp_presets.rst` |
| `docs/part3_gw_scissor/index.rst` | 1112–1130 | Preset 4: Full Quasiparticle Self-Consistency (``qsgw-dim``) | `docs/workflows/qp_presets.rst` |
| `docs/part3_gw_scissor/index.rst` | 1131–1146 | 14. Combining Electronic Structure Analysis with QP Shifts | `docs/electronic_structure/qp_analysis.rst` |
| `docs/part3_gw_scissor/index.rst` | 1147–1149 | 15. CLI Flags & YAML Configuration Reference | `docs/quasiparticles/configuration.rst` |
| `docs/part3_gw_scissor/index.rst` | 1150–1187 | Command-Line Arguments | `docs/quasiparticles/configuration.rst` |
| `docs/part3_gw_scissor/index.rst` | 1188–1205 | YAML Configuration Reference | `docs/quasiparticles/configuration.rst` |
| `docs/part4_excited_states/index.rst` | 1–9 | Part 4: Optical Excitations & Four Excited-State Frameworks | `docs/excitons/foundations.rst` |
| `docs/part4_excited_states/index.rst` | 10–22 | 1. The Two-Particle Excitation Problem | `docs/excitons/foundations.rst` |
| `docs/part4_excited_states/index.rst` | 23–51 | Tamm-Dancoff Approximation (TDA) | `docs/excitons/foundations.rst` |
| `docs/part4_excited_states/index.rst` | 52–66 | Singlet vs. Triplet Matrix Elements | `docs/excitons/foundations.rst` |
| `docs/part4_excited_states/index.rst` | 67–79 | Relativistic 2-Component Spinor BSE | `docs/excitons/foundations.rst` |
| `docs/part4_excited_states/index.rst` | 80–142 | 2. Key Distinction: ``2e-integrals`` vs. ``kernel`` | `docs/interactions/architecture.rst` |
| `docs/part4_excited_states/index.rst` | 143–145 | 3. Two-Electron Integral Representations (``2e-integrals``) | `docs/interactions/representations.rst` |
| `docs/part4_excited_states/index.rst` | 146–179 | Semi-Empirical Atom-Centered Representation (``2e-integrals: mnok``) | `docs/interactions/representations.rst` |
| `docs/part4_excited_states/index.rst` | 180–206 | Exact Analytical Gaussian Representation (``2e-integrals: xs``) | `docs/interactions/representations.rst` |
| `docs/part4_excited_states/index.rst` | 207–224 | Implementation in QDEX (``2e-integrals``) | `docs/interactions/representations.rst` |
| `docs/part4_excited_states/index.rst` | 225–229 | 4. Dielectric Screening Kernels (``kernel``) | `docs/interactions/screening.rst` |
| `docs/part4_excited_states/index.rst` | 230–252 | 1. Resta Screened Dielectric Kernel (``kernel: resta``) | `docs/interactions/screening.rst` |
| `docs/part4_excited_states/index.rst` | 253–263 | 2. Atomistic Discrete Dipole Interaction Kernel (``kernel: dim``) | `docs/interactions/screening.rst` |
| `docs/part4_excited_states/index.rst` | 264–276 | 3. Parameter-Free Microscopic ZDO-RPA Kernel (``kernel: rpa`` / ``xs-rpa``) | `docs/interactions/screening.rst` |
| `docs/part4_excited_states/index.rst` | 277–285 | 4. Simplified BSE Kernel (``kernel: sbse``) | `docs/interactions/screening.rst` |
| `docs/part4_excited_states/index.rst` | 286–290 | 5. Uniform Dielectric Kernel (``kernel: bse``) | `docs/interactions/screening.rst` |
| `docs/part4_excited_states/index.rst` | 291–313 | Spatial Asymptotics & Wannier-Mott Bulk Limit | `docs/validation/bulk_exciton_limit.rst` |
| `docs/part4_excited_states/index.rst` | 314–338 | 5. The Four Excitation Frameworks (``excitation_mode``) | `docs/excitons/frameworks.rst` |
| `docs/part4_excited_states/index.rst` | 339–352 | Why Diagonal BSE Works in Nanocrystals | `docs/excitons/frameworks.rst` |
| `docs/part4_excited_states/index.rst` | 353–368 | Implementation in QDEX (``--excitation-mode``) | `docs/excitons/frameworks.rst` |
| `docs/part4_excited_states/index.rst` | 369–395 | 6. Transition Dipoles, Oscillator Strengths & Superradiance | `docs/spectroscopy/dipoles_oscillators.rst` |
| `docs/part4_excited_states/index.rst` | 396–410 | 7. Full BSE & The Davidson Iterative Solver | `docs/excitons/solvers.rst` |
| `docs/part4_excited_states/index.rst` | 411–460 | 8. Comprehensive Exciton Calculation Matrix | `docs/excitons/model_selection.rst` |
| `docs/part4_excited_states/index.rst` | 461–463 | 9. Recommended Workflow Presets | `docs/workflows/exciton_presets.rst` |
| `docs/part4_excited_states/index.rst` | 464–485 | Preset 1: Standard Colloidal QD Absorption Spectrum | `docs/workflows/exciton_presets.rst` |
| `docs/part4_excited_states/index.rst` | 486–500 | Preset 2: Ultrafast Non-Adiabatic Molecular Dynamics (NAMD) | `docs/workflows/exciton_presets.rst` |
| `docs/part4_excited_states/index.rst` | 501–520 | Preset 3: Benchmark First-Principles Calculation | `docs/workflows/exciton_presets.rst` |
| `docs/part4_excited_states/index.rst` | 521–543 | Preset 4: Spin-Orbit Coupling & Dark Excitons | `docs/workflows/exciton_presets.rst` |
| `docs/part4_excited_states/index.rst` | 544–546 | 10. CLI Flags & YAML Configuration Reference | `docs/excitons/configuration.rst` |
| `docs/part4_excited_states/index.rst` | 547–592 | Command-Line Arguments | `docs/excitons/configuration.rst` |
| `docs/part5_exciton_analysis/index.rst` | 1–9 | Part 5: Excited-State Wavefunction Analysis (Plasser-Dreuw) | `docs/exciton_analysis/transition_density.rst` |
| `docs/part5_exciton_analysis/index.rst` | 10–24 | 1. The Two-Particle Transition Density Matrix | `docs/exciton_analysis/transition_density.rst` |
| `docs/part5_exciton_analysis/index.rst` | 25–43 | Reduced Hole and Electron Densities | `docs/exciton_analysis/transition_density.rst` |
| `docs/part5_exciton_analysis/index.rst` | 44–60 | Vectorized Low-Memory Mulliken Population | `docs/exciton_analysis/transition_density.rst` |
| `docs/part5_exciton_analysis/index.rst` | 61–65 | 2. Rigorous Plasser-Dreuw Spatial Descriptors | `docs/exciton_analysis/descriptors.rst` |
| `docs/part5_exciton_analysis/index.rst` | 66–75 | 1. Spatial Centroids (:math:`\langle \mathbf{r}_h \rangle, \langle \mathbf{r}_e \rangle`) | `docs/exciton_analysis/descriptors.rst` |
| `docs/part5_exciton_analysis/index.rst` | 76–87 | 2. Charge-Transfer Distance (:math:`d_{\mathrm{CT}}`) | `docs/exciton_analysis/descriptors.rst` |
| `docs/part5_exciton_analysis/index.rst` | 88–102 | 3. Root-Mean-Square Particle Sizes (:math:`\sigma_h, \sigma_e`) | `docs/exciton_analysis/descriptors.rst` |
| `docs/part5_exciton_analysis/index.rst` | 103–113 | 4. Electron-Hole Spatial Covariance (:math:`\mathrm{Cov}_{eh}`) | `docs/exciton_analysis/descriptors.rst` |
| `docs/part5_exciton_analysis/index.rst` | 114–126 | 5. Pearson Correlation Coefficient (:math:`R_{eh}`) | `docs/exciton_analysis/descriptors.rst` |
| `docs/part5_exciton_analysis/index.rst` | 127–137 | 6. True Root-Mean-Square Exciton Size (:math:`d_{eh}`) | `docs/exciton_analysis/descriptors.rst` |
| `docs/part5_exciton_analysis/index.rst` | 138–146 | 7. Charge-Transfer Ratio (:math:`\mathrm{CT}`) | `docs/exciton_analysis/descriptors.rst` |
| `docs/part5_exciton_analysis/index.rst` | 147–160 | 8. Exciton Participation Ratio (:math:`\mathrm{PR}`) | `docs/exciton_analysis/descriptors.rst` |
| `docs/part5_exciton_analysis/index.rst` | 161–184 | 3. Natural Transition Orbitals (NTOs) | `docs/exciton_analysis/ntos.rst` |
| `docs/part5_exciton_analysis/index.rst` | 185–195 | NTO Compactness Metrics | `docs/exciton_analysis/ntos.rst` |
| `docs/part5_exciton_analysis/index.rst` | 196–208 | 4. 3D Volumetric Visualization (.cube Files) | `docs/exciton_analysis/cubes.rst` |
| `docs/part5_exciton_analysis/index.rst` | 209–222 | 5. Interactive Plotly 6-Panel Dashboard | `docs/exciton_analysis/dashboards.rst` |
| `docs/part5_exciton_analysis/index.rst` | 223–225 | 6. CLI Flags & YAML Configuration Reference | `docs/exciton_analysis/configuration.rst` |
| `docs/part5_exciton_analysis/index.rst` | 226–260 | Command-Line Arguments | `docs/exciton_analysis/configuration.rst` |
| `docs/part5_exciton_analysis/index.rst` | 261–276 | YAML Configuration Example | `docs/exciton_analysis/configuration.rst` |
| `docs/part6_namd/index.rst` | 1–9 | Part 6: Carrier Cooling Dynamics & Photoluminescence (NAMD) | `docs/dynamics/pipeline.rst` |
| `docs/part6_namd/index.rst` | 10–44 | 1. Overview & The NAMD Pipeline Architecture | `docs/dynamics/pipeline.rst` |
| `docs/part6_namd/index.rst` | 45–49 | 2. Multi-Timescale Integration: Separating Nuclear and Electronic Time Steps | `docs/dynamics/timesteps_cpa.rst` |
| `docs/part6_namd/index.rst` | 50–70 | The Timescale Mismatch | `docs/dynamics/timesteps_cpa.rst` |
| `docs/part6_namd/index.rst` | 71–81 | The Classical Path Approximation (CPA) | `docs/dynamics/timesteps_cpa.rst` |
| `docs/part6_namd/index.rst` | 82–132 | Electronic Sub-Stepping in CPA-FSSH | `docs/dynamics/timesteps_cpa.rst` |
| `docs/part6_namd/index.rst` | 133–153 | Electronic Sub-Stepping in the Pauli Master Equation (PME) | `docs/dynamics/timesteps_cpa.rst` |
| `docs/part6_namd/index.rst` | 154–187 | 3. Method Selection: PME vs. CPA-FSSH-EDC vs. DISH | `docs/dynamics/model_selection.rst` |
| `docs/part6_namd/index.rst` | 188–200 | The Phonon Bottleneck Case (:math:`1P_e \to 1S_e`) | `docs/dynamics/model_selection.rst` |
| `docs/part6_namd/index.rst` | 201–211 | Surface Trap States: Hopping & De-Hopping Kinetics | `docs/dynamics/model_selection.rst` |
| `docs/part6_namd/index.rst` | 212–221 | Summary Decision Rule | `docs/dynamics/model_selection.rst` |
| `docs/part6_namd/index.rst` | 222–224 | 4. Theoretical Foundations of the Dynamical Engines | `docs/dynamics/pme.rst` |
| `docs/part6_namd/index.rst` | 225–286 | 1. Derivation of the Pauli Master Equation (PME) | `docs/dynamics/pme.rst` |
| `docs/part6_namd/index.rst` | 287–413 | 2. Classical Path Approximation Surface Hopping with Energy-Based Decoherence (CPA-FSSH-EDC) | `docs/dynamics/fssh_edc.rst` |
| `docs/part6_namd/index.rst` | 414–522 | 3. Decoherence-Induced Surface Hopping (DISH) | `docs/dynamics/dish.rst` |
| `docs/part6_namd/index.rst` | 523–575 | 4. Comparative Synthesis: PME vs. CPA-FSSH-EDC vs. DISH | `docs/dynamics/comparison.rst` |
| `docs/part6_namd/index.rst` | 576–578 | 5. Trajectory Precomputation & Wavefunction Tracking | `docs/dynamics/nacs_tracking.rst` |
| `docs/part6_namd/index.rst` | 579–589 | Numerical Non-Adiabatic Couplings (NAC) | `docs/dynamics/nacs_tracking.rst` |
| `docs/part6_namd/index.rst` | 590–606 | Eliminating Gauge Phase Discontinuities | `docs/dynamics/nacs_tracking.rst` |
| `docs/part6_namd/index.rst` | 607–621 | Hungarian Matching for Trivial Crossings | `docs/dynamics/nacs_tracking.rst` |
| `docs/part6_namd/index.rst` | 622–624 | 6. Electronic Decoherence: Origin, Computation, and Rationale | `docs/dynamics/decoherence.rst` |
| `docs/part6_namd/index.rst` | 625–643 | Physical Origin of Electronic Decoherence | `docs/dynamics/decoherence.rst` |
| `docs/part6_namd/index.rst` | 644–648 | Why Dephasing is Ultrafast in Nanocrystals | `docs/dynamics/decoherence.rst` |
| `docs/part6_namd/index.rst` | 649–685 | Automated Ab Initio Cumulant Decoherence | `docs/dynamics/decoherence.rst` |
| `docs/part6_namd/index.rst` | 686–746 | State-Pair Pure-Dephasing Matrices (tau_ij) | `docs/dynamics/decoherence.rst` |
| `docs/part6_namd/index.rst` | 747–786 | Surface Hopping Schemes: FSSH-EDC vs. DISH | `docs/dynamics/decoherence.rst` |
| `docs/part6_namd/index.rst` | 787–789 | 7. Phonon Spectral Density J(ω): Mapping Electron-Phonon Coupling | `docs/spectroscopy/spectral_density.rst` |
| `docs/part6_namd/index.rst` | 790–800 | Mathematical Definition | `docs/spectroscopy/spectral_density.rst` |
| `docs/part6_namd/index.rst` | 801–809 | Physical Meaning | `docs/spectroscopy/spectral_density.rst` |
| `docs/part6_namd/index.rst` | 810–820 | Identifying Active Phonon Modes During Cooling | `docs/spectroscopy/spectral_density.rst` |
| `docs/part6_namd/index.rst` | 821–831 | 8. Radiative & Non-Radiative Recombination Mechanisms | `docs/recombination/overview.rst` |
| `docs/part6_namd/index.rst` | 832–890 | 1. Microscopic Origin: From DFT and BSE to Recombination | `docs/recombination/overview.rst` |
| `docs/part6_namd/index.rst` | 891–918 | 2. Einstein Radiative Rate: Single-Frame vs. NAMD Trajectory Averaging | `docs/recombination/radiative.rst` |
| `docs/part6_namd/index.rst` | 919–945 | 3. Non-Radiative Decay Across Large Gaps: Englman-Jortner Energy Gap Law | `docs/recombination/energy_gap_law.rst` |
| `docs/part6_namd/index.rst` | 946–982 | 4. Non-Empirical Extraction of Optical Phonon Energy from NAMD Spectral Density | `docs/recombination/trajectory_parameters.rst` |
| `docs/part6_namd/index.rst` | 983–1012 | 5. Derivation of Huang-Rhys Factor S and Reorganization Energy λ from Trajectory Data | `docs/recombination/trajectory_parameters.rst` |
| `docs/part6_namd/index.rst` | 1013–1044 | 6. Intermediate & Narrow Gap Decay: Franck-Condon Weighted Density of States (FCWD) | `docs/recombination/fcwd.rst` |
| `docs/part6_namd/index.rst` | 1045–1058 | 7. Defect Trap-Assisted Recombination (Shockley-Read-Hall) | `docs/recombination/traps.rst` |
| `docs/part6_namd/index.rst` | 1059–1072 | 8. Photoluminescence Quantum Yield (PLQY) | `docs/recombination/plqy.rst` |
| `docs/part6_namd/index.rst` | 1073–1077 | 9. In-Depth Analysis of NAMD Simulations | `docs/dynamics/analysis.rst` |
| `docs/part6_namd/index.rst` | 1078–1134 | 1. Carrier Cooling Curves, Lifetimes, and Band Edge Arrival Times | `docs/dynamics/analysis.rst` |
| `docs/part6_namd/index.rst` | 1135–1139 | 2. State-Resolved Population Kinetics | `docs/dynamics/analysis.rst` |
| `docs/part6_namd/index.rst` | 1140–1144 | 3. NAC vs. Energy Gap Distribution | `docs/dynamics/analysis.rst` |
| `docs/part6_namd/index.rst` | 1145–1157 | 4. 6-Panel Publication Figures & Dashboards | `docs/dynamics/analysis.rst` |
| `docs/part6_namd/index.rst` | 1158–1168 | 10. Ultrafast Pump-Probe Transient Absorption (TA) Spectroscopy | `docs/spectroscopy/transient_absorption.rst` |
| `docs/part6_namd/index.rst` | 1169–1182 | Physical Mechanisms in Nanocrystal Transient Absorption | `docs/spectroscopy/transient_absorption.rst` |
| `docs/part6_namd/index.rst` | 1183–1207 | Microscopic Formulation & State-Filling Factors | `docs/spectroscopy/transient_absorption.rst` |
| `docs/part6_namd/index.rst` | 1208–1228 | 1S Bleach Kinetic Profiling & Carrier Cooling Rates (:math:`k_C`) | `docs/spectroscopy/transient_absorption.rst` |
| `docs/part6_namd/index.rst` | 1229–1238 | Schematic Diagrams of Transient Absorption Processes | `docs/spectroscopy/transient_absorption.rst` |
| `docs/part6_namd/index.rst` | 1239–1279 | THREE OPTICAL MECHANISMS IN TRANSIENT ABSORPTION SPECTROSCOPY (ΔA) | `docs/spectroscopy/transient_absorption.rst` |
| `docs/part6_namd/index.rst` | 1280–1310 | HOT-CARRIER RELAXATION CASCADE & 1S BLEACH KINETIC RISE PROFILE | `docs/spectroscopy/transient_absorption.rst` |
| `docs/part6_namd/index.rst` | 1311–1350 | How QDEX Data Are Used to Compute Every Formula Term | `docs/spectroscopy/transient_absorption.rst` |
| `docs/part6_namd/index.rst` | 1351–1360 | Publication-Quality Visualizations | `docs/spectroscopy/transient_absorption.rst` |
| `docs/part6_namd/index.rst` | 1361–1363 | 11. CLI Flags & YAML Configuration Reference | `docs/dynamics/configuration.rst` |
| `docs/part6_namd/index.rst` | 1364–1416 | Command-Line Arguments | `docs/dynamics/configuration.rst` |
| `docs/part6_namd/index.rst` | 1417–1448 | YAML Configuration Example | `docs/dynamics/configuration.rst` |
| `docs/part7_examples/index.rst` | 1–14 | Part 7: Step-by-Step Tutorials & Workflows | `docs/workflows/index.rst` |
| `docs/part7_examples/index.rst` | 15–19 | Tutorial 1: Ground-State Electronic Structure Analysis | `docs/workflows/electronic_structure.rst` |
| `docs/part7_examples/index.rst` | 20–27 | Input Files Required | `docs/workflows/electronic_structure.rst` |
| `docs/part7_examples/index.rst` | 28–55 | Configuration File (``tutorial1_ground_state.yaml``) | `docs/workflows/electronic_structure.rst` |
| `docs/part7_examples/index.rst` | 56–64 | Execution | `docs/workflows/electronic_structure.rst` |
| `docs/part7_examples/index.rst` | 65–76 | Generated Output Files | `docs/workflows/electronic_structure.rst` |
| `docs/part7_examples/index.rst` | 77–81 | Tutorial 2: Incorporating Relativistic Spin-Orbit Coupling | `docs/workflows/soc.rst` |
| `docs/part7_examples/index.rst` | 82–107 | Configuration File (``tutorial2_soc.yaml``) | `docs/workflows/soc.rst` |
| `docs/part7_examples/index.rst` | 108–114 | Execution | `docs/workflows/soc.rst` |
| `docs/part7_examples/index.rst` | 115–123 | What to Observe | `docs/workflows/soc.rst` |
| `docs/part7_examples/index.rst` | 124–128 | Tutorial 3: Scaled GW Band Gap Correction & Absolute Edges | `docs/workflows/qp_edges.rst` |
| `docs/part7_examples/index.rst` | 129–153 | Configuration File (``tutorial3_gw.yaml``) | `docs/workflows/qp_edges.rst` |
| `docs/part7_examples/index.rst` | 154–160 | Execution | `docs/workflows/qp_edges.rst` |
| `docs/part7_examples/index.rst` | 161–178 | Key Output | `docs/workflows/qp_edges.rst` |
| `docs/part7_examples/index.rst` | 179–183 | Tutorial 4: Comparing the Four Excited-State Frameworks | `docs/workflows/four_frameworks.rst` |
| `docs/part7_examples/index.rst` | 184–210 | Compare via CLI | `docs/workflows/four_frameworks.rst` |
| `docs/part7_examples/index.rst` | 211–219 | Comparison Summary | `docs/workflows/four_frameworks.rst` |
| `docs/part7_examples/index.rst` | 220–224 | Tutorial 5: Exciton Wavefunction Descriptors & NTO Analysis | `docs/workflows/exciton_analysis.rst` |
| `docs/part7_examples/index.rst` | 225–258 | Configuration File (``tutorial5_analysis.yaml``) | `docs/workflows/exciton_analysis.rst` |
| `docs/part7_examples/index.rst` | 259–265 | Execution | `docs/workflows/exciton_analysis.rst` |
| `docs/part7_examples/index.rst` | 266–274 | Analysis Results | `docs/workflows/exciton_analysis.rst` |
| `docs/part7_examples/index.rst` | 275–279 | Tutorial 6: Carrier Cooling Dynamics & Photoluminescence (NAMD) | `docs/workflows/carrier_cooling.rst` |
| `docs/part7_examples/index.rst` | 280–290 | Directory Layout | `docs/workflows/carrier_cooling.rst` |
| `docs/part7_examples/index.rst` | 291–323 | Configuration File (``tutorial6_namd.yaml``) | `docs/workflows/carrier_cooling.rst` |
| `docs/part7_examples/index.rst` | 324–332 | Step 1: Trajectory Precomputation | `docs/workflows/carrier_cooling.rst` |
| `docs/part7_examples/index.rst` | 333–339 | Step 2: Carrier Cooling Simulation | `docs/workflows/carrier_cooling.rst` |
| `docs/part7_examples/index.rst` | 340–346 | Results & Visualizations | `docs/workflows/carrier_cooling.rst` |
| `docs/part8_auger/index.rst` | 1–11 | Part 8: Multi-Carrier Auger Recombination in Quantum Dots | `docs/recombination/auger_foundations.rst` |
| `docs/part8_auger/index.rst` | 12–33 | 1. Physical Foundations of Auger Scattering in Quantum Dots | `docs/recombination/auger_foundations.rst` |
| `docs/part8_auger/index.rst` | 34–72 | Breakdown of Momentum Conservation & :math:`1/V` Volume Scaling | `docs/recombination/auger_foundations.rst` |
| `docs/part8_auger/index.rst` | 73–81 | 2. Many-Body Formulation & Matrix Elements | `docs/recombination/auger_matrix_elements.rst` |
| `docs/part8_auger/index.rst` | 82–110 | The Two-Body Screened Coulomb Operator | `docs/recombination/auger_matrix_elements.rst` |
| `docs/part8_auger/index.rst` | 111–123 | Spin Summation & Multiplicity | `docs/recombination/auger_matrix_elements.rst` |
| `docs/part8_auger/index.rst` | 124–138 | Universal Multiexciton Statistical Scaling | `docs/recombination/auger_matrix_elements.rst` |
| `docs/part8_auger/index.rst` | 139–148 | Schematic Diagrams of the Auger Processes | `docs/recombination/auger_matrix_elements.rst` |
| `docs/part8_auger/index.rst` | 149–191 | NEGATIVE TRION / BIEXCITON: eeh (Electron Ejected) | `docs/recombination/auger_matrix_elements.rst` |
| `docs/part8_auger/index.rst` | 192–234 | POSITIVE TRION / BIEXCITON: hhe (Hole Ejected) | `docs/recombination/auger_matrix_elements.rst` |
| `docs/part8_auger/index.rst` | 235–262 | NEUTRAL BIEXCITON (XX) ANNIHILATION & ENERGY-CONSERVING SURFACE HOPPING (ECSH) CASCADE | `docs/recombination/auger_matrix_elements.rst` |
| `docs/part8_auger/index.rst` | 263–273 | 3. The Atom-Centered Resta Screening Contraction | `docs/recombination/auger_screening.rst` |
| `docs/part8_auger/index.rst` | 274–296 | Monopole Transition Charge Projection | `docs/recombination/auger_screening.rst` |
| `docs/part8_auger/index.rst` | 297–314 | The Microscopic Resta Screening Kernel | `docs/recombination/auger_screening.rst` |
| `docs/part8_auger/index.rst` | 315–335 | Local vs. Macroscopic Dielectric Screening in Auger Scattering | `docs/recombination/auger_screening.rst` |
| `docs/part8_auger/index.rst` | 336–383 | 4. How QDEX Data Are Used to Compute Every Formula Term | `docs/recombination/auger_implementation.rst` |
| `docs/part8_auger/index.rst` | 384–388 | 5. Energy Conservation Line Shapes: Gaussian vs. FCWD | `docs/recombination/auger_lineshapes.rst` |
| `docs/part8_auger/index.rst` | 389–399 | 1. Gaussian Line Shape Model | `docs/recombination/auger_lineshapes.rst` |
| `docs/part8_auger/index.rst` | 400–412 | 2. Multi-Phonon Marcus / Jortner Line Shape (FCWD) | `docs/recombination/auger_lineshapes.rst` |
| `docs/part8_auger/index.rst` | 413–431 | 6. Relativistic Spin-Orbit Coupling (SOC) in Auger Scattering | `docs/recombination/auger_soc.rst` |
| `docs/part8_auger/index.rst` | 432–434 | 7. Static Calculations vs. Trajectory Dynamics in NAMD | `docs/recombination/auger_workflows.rst` |
| `docs/part8_auger/index.rst` | 435–448 | Static Calculation | `docs/recombination/auger_workflows.rst` |
| `docs/part8_auger/index.rst` | 449–452 | QDEX AUGER RECOMBINATION REPORT | `docs/recombination/auger_workflows.rst` |
| `docs/part8_auger/index.rst` | 453–454 | Evaluated Active Pathways       : 177 (eeh)  /  48 (hhe) | `docs/recombination/auger_workflows.rst` |
| `docs/part8_auger/index.rst` | 455–458 | Channel                 Rate (s^-1)        Rate (ns^-1)        Lifetime (ns)          Lifetime (ps) | `docs/recombination/auger_workflows.rst` |
| `docs/part8_auger/index.rst` | 459–461 | Biexciton (XX)             1.8990e+11        1.8990e+02            0.0053 ns            5.27 ps | `docs/recombination/auger_workflows.rst` |
| `docs/part8_auger/index.rst` | 462–474 | Trajectory-Averaged Dynamic Auger Rates in NAMD | `docs/recombination/auger_workflows.rst` |
| `docs/part8_auger/index.rst` | 475–477 | 8. Energy-Conserving Surface Hopping (ECSH) for Auger in NAMD | `docs/recombination/auger_ecsh.rst` |
| `docs/part8_auger/index.rst` | 478–500 | Physical Foundations: Electron-Phonon vs. Coulomb Transitions | `docs/recombination/auger_ecsh.rst` |
| `docs/part8_auger/index.rst` | 501–507 | The Fundamental Flaw in Previous Approaches (The "Energy Leak" Error) | `docs/recombination/auger_ecsh.rst` |
| `docs/part8_auger/index.rst` | 508–529 | The Gumber-Prezhdo ECSH Methodology | `docs/recombination/auger_ecsh.rst` |
| `docs/part8_auger/index.rst` | 530–546 | Relation with Static Trions (:math:`eeh` and :math:`hhe`) | `docs/recombination/auger_ecsh.rst` |
| `docs/part8_auger/index.rst` | 547–582 | Resolving the Timescale Mismatch (:math:`1 - 10\text{ ps}` MD vs. :math:`100\text{ ps} - 10\text{ ns}` Auger) | `docs/recombination/auger_ecsh.rst` |
| `docs/part8_auger/index.rst` | 583–585 | 9. Configuration Reference (YAML & CLI) | `docs/recombination/auger_configuration.rst` |
| `docs/part8_auger/index.rst` | 586–607 | YAML Configuration Options | `docs/recombination/auger_configuration.rst` |
| `docs/part8_auger/index.rst` | 608–623 | Command-Line Arguments | `docs/recombination/auger_configuration.rst` |
