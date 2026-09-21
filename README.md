# QDEX: Quantum Dot Excitations & Dynamics

**QDEX** (pronounced *qiu-di-ex*, `/kjuː-diː-ɛks/`, formerly `miniBSE`) is a high-performance, lightweight post-DFT exciton solver and non-adiabatic molecular dynamics (NAMD) engine designed for calculating, analyzing, and propagating the excited states of semiconductor nanoclusters and molecules.

By combining the ease of a Python interface with a lightning-fast C++ backend powered by `Libint2` and `Eigen3`, `QDEX` provides researchers with an end-to-end framework—from ground-state Kohn-Sham orbitals to carrier cooling dynamics, trap-assisted recombination, and photoluminescence quantum yields (PLQY)—without the overhead of massive quantum chemistry suites.

---

## What QDEX Does

`QDEX` bridges the gap between static DFT calculations and real-time excited-state dynamics across seven integrated modules:

1. **Ground-State Electronic Structure**: 
   - Projects CP2K Kohn-Sham molecular orbitals (MOs) onto atomic sites via Lowdin and Mulliken population analysis.
   - Computes Projected Density of States (PDOS), Inverse Participation Ratio (IPR), and Crystal Orbital Overlap Population (COOP).
   - Unfolds discrete nanocluster states onto bulk crystal $k$-paths using the plane-wave **Fuzzy Bands** algorithm.
   - Generates 3D volumetric Gaussian `.cube` files with multi-threaded C++ evaluation.
2. **Relativistic Spin-Orbit Coupling (SOC)**:
   - Evaluates relativistic SOC using norm-conserving separable Goedecker-Teter-Hutter (GTH) pseudopotential projectors.
   - Fast sub-second spinor diagonalization via sparse CSR angular momentum matrices and purely real BLAS Level-3 DGEMM contractions.
3. **Scaled GW Quasiparticle Model**:
   - Resolves DFT band-gap underestimation via a two-anchor physical scaling model (vacuum monomer anchor $R_0$ and bulk ARPES limit).
   - Predicts absolute Ionization Potential (IP) and Electron Affinity (EA) with dielectric solvation screening.
4. **Optical Excitations & BSE**:
   - Solves the Bethe-Salpeter Equation (BSE) under the Tamm-Dancoff Approximation (TDA) using atom-centered MNOK transition charges and the Resta screened dielectric kernel.
   - Supports four progressive excitation approximations: `independent_dft`, `independent_qp`, `diagonal_bse`, and full `bse` (sTDA) via an optimized Davidson iterative subspace solver.
5. **Rigorous Wavefunction Analysis (Plasser-Dreuw)**:
   - Evaluates real-space exciton descriptors: exciton size ($d_{eh}$), charge-transfer distance ($d_{\mathrm{CT}}$), electron/hole cloud spread ($\sigma_e, \sigma_h$), and spatial Pearson correlation ($R_{eh}$).
   - Extracts Natural Transition Orbitals (NTOs) with participation ratios and Shannon entropy metrics.
   - Exports interactive 6-panel Plotly dashboards (`exciton_analysis.html`).
6. **Non-Adiabatic Molecular Dynamics (NAMD) & Carrier Cooling**:
   - Operates within the Classical Path Approximation (CPA) along *ab initio* molecular dynamics (AIMD) trajectories.
   - Computes analytic cross-frame overlaps $S_{IJ}(t, t+\Delta t)$ via Libint2, eliminating gauge phase flips and tracking trivial crossings via Hungarian matching.
   - Propagates carrier cooling deterministically via the tensorized Pauli Master Equation (PME) or stochastically via Fewest Switches Surface Hopping (CPA-FSSH).
   - Solves for *ab initio* decoherence times $\tau_{\mathrm{dec}}$ via second-order cumulant expansion of energy gap fluctuations.
   - Evaluates spontaneous emission rates ($k_{\mathrm{rad}}$), multi-phonon non-radiative rates ($k_{\mathrm{nonrad}}$ via Jortner, Marcus, and SRH defect models), carrier cooling lifetimes ($\tau_{\mathrm{cool}}$), band-edge arrival times, and photoluminescence quantum yields (PLQY).
7. **Multi-Carrier Auger Recombination**:
   - Solves Fermi's Golden Rule for negative trion ($eeh$) and positive trion ($hhe$) channels and neutral biexciton lifetimes ($\tau_{XX}$).
   - Contracts transition densities into atom-centered charges and uses the microscopic Resta dielectric kernel, completely avoiding $O(N_{\mathrm{ao}}^5)$ four-center integrals.
   - Supports both static single-geometry evaluation and trajectory-averaged NAMD dynamics.

---

## Installation

Because `QDEX` relies on C++ extensions compiled against `Libint2` and `Eigen3`, **Conda / Micromamba** is the recommended installation method.

### Method 1: Conda / Micromamba (Recommended)

1. Clone the repository:
   ```bash
   git clone https://github.com/nlesc-nano/miniBSE.git
   cd miniBSE
   ```
2. Create and activate the environment:
   ```bash
   micromamba env create -f environment.yml
   micromamba activate minibse_env
   ```
   *(Note: The environment file automatically installs `QDEX` in editable mode via pip).*

### Method 2: Standard Pip

If you already have `CMake` (>= 3.16), a C++17 compiler, `Eigen3`, and `Libint2` installed natively on your OS:

```bash
git clone https://github.com/nlesc-nano/miniBSE.git
cd miniBSE
pip install -r requirements.txt
pip install -e .
```

> **Backward Compatibility**: Any existing scripts importing `miniBSE` (`import miniBSE`) continue to work seamlessly via an automatic redirection hook that maps to `qdex`. Both `qdex` and `minibse` CLI commands are available globally.

---

## How to Use It

Once installed, the engine is accessible globally via the `qdex` CLI (or `minibse`).

### 1. Typical CLI Calculation

Below is a standard example for calculating the excited states of a semiconductor cluster using the scaled GW quasiparticle gap and plotting the results:

```bash
qdex \
  --mo_file MOs.mbse \
  --xyz structure.xyz \
  --basis_txt BASIS_MOLOPT \
  --basis_name DZVP-MOLOPT-SR-GTH \
  --e_thresh 2.0 \
  --qp_gap gw \
  --material CSPBBR3 \
  --eps-out 2.4 \
  --excitation-mode diagonal_bse \
  --plot \
  --nthreads 8
```

### 2. YAML Configuration Workflow

For complex workflows and NAMD trajectories, using a structured YAML configuration file is recommended:

```yaml
# config.yaml
system:
  mo_file: "CsPbBr3_MOs.mbse"
  xyz: "CsPbBr3_QD.xyz"
  basis_txt: "BASIS_MOLOPT"
  basis_name: "DZVP-MOLOPT-PBE-GTH"
  material: "CSPBBR3"
  nthreads: 8

physics:
  excitation_mode: "diagonal_bse"
  qp_gap: "gw"
  kernel: "resta"
  eps_out: 2.4
  soc: true
  gth_file: "GTH_SOC_POTENTIALS.txt"
```

Run static calculation:
```bash
qdex --config config.yaml
```

Run NAMD carrier cooling precomputation and propagation:
```bash
qdex --config config.yaml --namd-precompute
qdex --config config.yaml --namd-dynamics
```

---

## Outputs

A run of `QDEX` produces rich publication-ready data and interactive dashboards:

1. **Standard Output**: Formatted console tables listing states, energies, transition dipoles, oscillator strengths, cooling lifetimes, band-edge arrival times, and PLQY.
2. **`spectrum.png` & `spectrum_nm.png`**: UV-Vis absorption spectra with Gaussian broadening.
3. **`exciton_analysis.html`**: Interactive Plotly dashboard visualizing spatial correlations, exciton radii ($d_{eh}$), and charge-transfer metrics across the spectrum.
4. **`exciton_results.csv`**: Tabular data with spatial metrics ($d_{\mathrm{CT}}, \sigma_h, \sigma_e, R_{eh}$) for trajectory post-processing.
5. **`namd_cooling_dashboard.html`**: Coordinated 6-panel NAMD dashboard displaying cooling curves, state population cascades, non-adiabatic coupling heatmaps, gap-law distributions, phonon spectral densities $J(\omega)$, and cumulant decoherence decays.
6. **Volumetric Cube Files**: Gaussian `.cube` files of frontier orbitals or exciton hole, electron, and difference densities ready for VMD, PyMOL, or ChimeraX.

---

## Documentation

Full documentation with comprehensive mathematical formulations, physical explanations, and step-by-step tutorials is available in the `docs/` folder and can be built using Sphinx:

```bash
sphinx-build -b html docs docs/_build/html
```

---

## Citation & License

If you use **QDEX** in your research, please cite:
* Ivan Infante et al., *QDEX: Quantum Dot Excitations & Dynamics* (2026).
* Licensed under the Apache License 2.0.
