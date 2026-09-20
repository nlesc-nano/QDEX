# miniBSE

**miniBSE** is a high-performance, lightweight post-DFT exciton solver designed for calculating and analyzing the excited states of molecules and semiconductor nanoclusters. 

By combining the ease of a Python interface with a lightning-fast C++ backend powered by `Libint2` and `Eigen3`, `miniBSE` offers researchers a scalable tool to investigate light-matter interactions, exciton delocalization, and charge-transfer (CT) characteristics without the overhead of massive quantum chemistry suites.

## What it Does

At its core, `miniBSE` solves the **Bethe-Salpeter Equation (BSE)** under the Tamm-Dancoff Approximation (TDA). It constructs an active-space electron-hole Hamiltonian using:
1. **DFT Ground State Data**: Takes Molecular Orbitals (MOs) and orbital energies from a prior ground-state DFT calculation.
2. **Analytic Integrals**: Uses `Libint2` to instantly compute Gaussian basis set overlaps and dipole transition matrices in real-space.
3. **Separated BSE/TDA interactions**: Uses unscreened MNOK transition charges for the bare exchange/local-field term and an electronic Resta-MNOK interaction for the screened direct electron-hole attraction.
4. **Several excitation levels**: Solves the static BSE in the Tamm-Dancoff approximation with Davidson or dense diagonalization, or evaluates uncoupled one-electron transitions without diagonalization.

Beyond calculating energies, `miniBSE` performs **extensive wavefunction analysis** based on the Dreuw/Plasser framework, outputting physical descriptors such as exciton radii ($d_{eh}$), spatial correlation (Pearson $R$), and volumetric transition densities.

---

## Installation

Because `miniBSE` relies on C++ extensions, **Conda is the highly recommended installation method**. Our `environment.yml` handles the installation of C++ compilers, `CMake`, `Eigen3`, and `Libint2`, saving you the hassle of system-level configurations.

### Method 1: Conda (Recommended)

1. Clone the repository:
   ```bash
   git clone [https://github.com/nlesc-nano/miniBSE.git](https://github.com/nlesc-nano/miniBSE.git)
   cd miniBSE
   ```
2. Create and activate the environment:
   ```bash
   micromamba env create -f environment.yml
   micromamba activate minibse_env
   ```
   *(Note: The `environment.yml` automatically installs the `miniBSE` package in editable mode via pip at the end of the process).*

### Method 2: Standard Pip

If you already have `CMake` (>= 3.16), a C++17 compiler, `Eigen3`, and `Libint2` installed natively on your OS:

```bash
git clone [https://github.com/nlesc-nano/miniBSE.git](https://github.com/nlesc-nano/miniBSE.git)
cd miniBSE
pip install -r requirements.txt
pip install -e .
```

---

## How to Use It

Once installed, the solver is accessible globally via the `minibse` command-line interface. 

### Typical Calculation

Below is a standard example for calculating the excited states of an Indium Arsenide (InAs) semiconductor cluster, computing the lowest states within a 2 eV threshold, and plotting the results:

```bash
minibse \
  --mo_file MOs_cleaned.txt \
  --xyz last_opt.xyz \
  --basis_txt BASIS_MOLOPT \
  --basis_name DZVP-MOLOPT-SR-GTH \
  --e_thresh 2 \
  --qp_gap 3.0278 \
  --sigma 0.03 \
  --plot \
  --full-diag \
  --nthreads 8 \
  --material INAS \
  --exchange \
  --alpha 0.2
```

### CLI Argument Breakdown

**Inputs & Structure:**
* `--mo_file`: Path to your molecular orbitals (supports `.txt` or `.npz` arrays).
* `--xyz`: The Cartesian coordinates of your system.
* `--basis_txt` & `--basis_name`: The basis set file and the specific basis name (e.g., CP2K MOLOPT format) used to generate the C++ integrals.

**Physics & Truncation:**
* `--qp_gap`: The target quasi-particle gap (in eV). `miniBSE` uses this to apply a "scissor shift" to the raw DFT HOMO-LUMO gap.
* `--e_thresh`: Energy threshold (in eV). Truncates the active space by discarding electron-hole transitions that exceed this gap.
* `--material`: Uses a built-in material database to estimate dielectric screening.
* `--eps-out`: External dielectric used by the finite-size QP polarization correction. It does not replace the material electronic dielectric inside the microscopic Resta kernel.
* `--include-direct-eh` / `--no-direct-eh`: Enables or disables the screened attractive electron-hole direct term. The old `exchange: true` YAML key remains a deprecated compatibility alias.
* `--alpha`: Scales only the legacy non-Resta kernel. It is ignored, with a warning, when `kernel: resta` is selected.

**Excitation approximations:**

* `--excitation-mode bse`: Diagonalize the coupled BSE/TDA Hamiltonian, $D_{QP}+K_x^{bare}-K_d^{screened}$ (the default).
* `--excitation-mode independent_dft`: Do not diagonalize; use the underlying DFT occupied-to-virtual energy differences.
* `--excitation-mode independent_qp`: Do not diagonalize; use rigid-scissor or otherwise selected QP occupied-to-virtual energy differences.
* `--excitation-mode diagonal_bse`: Do not diagonalize; correct each QP transition by its diagonal bare $K_x$ and screened $K_d$ matrix elements. Off-diagonal configuration mixing is omitted.

**Solver & Output Controls:**
* `--full-diag`: Forces full dense diagonalization of the resonant BSE/TDA matrix. This is not a non-TDA BSE calculation. For an independent-transition mode it requests all retained transitions without diagonalizing a matrix.
* `--nthreads`: Number of CPU threads dedicated to C++ integral generation and PyTorch matrix contractions.
* `--sigma`: Broadening width (in eV) for the generated UV-Vis spectrum.
* `--plot`: Generates PNG spectra and an interactive HTML diagnostic dashboard.
* `--cube`: (Optional) Generates 3D volumetric `.cube` files of the brightest exciton's electron/hole densities.

---

## Outputs

A successful run of `miniBSE` will yield several outputs in your working directory:

1. **Standard Output**: A console table listing the excited states, energies, oscillator strengths, and primary orbital transitions (e.g., `HOMO -> LUMO+1`).
2. **`spectrum.png` & `spectrum_nm.png`**: UV-Vis absorption spectra utilizing your requested broadening (`--sigma`). 
3. **`exciton_analysis.html`**: An interactive Plotly dashboard. This visualizes exciton spatial correlations, sizes, and charge-transfer ratios across the energy spectrum.
4. **`exciton_results.csv`**: (If `--write-csv` is used) Tabular data containing detailed spatial metrics ($d_{CT}$, $\sigma_h$, $\sigma_e$) for post-processing or tracking across MD trajectories.
5. **Cube Files**: (If `--cube` is used) Volumetric densities ready to be visualized in software like VMD, PyMOL, or ChimeraX.
