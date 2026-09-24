Quickstart Guide
================

`QDEX` can be run either directly via the command-line interface (CLI) or through structured YAML configuration files.

1. Static BSE Calculation via CLI
---------------------------------

Below is a standard calculation for the excited states of an Indium Arsenide (InAs) semiconductor cluster:

.. code-block:: bash

   qdex \
     --mo_file MOs.mbse \
     --xyz structure.xyz \
     --basis_txt BASIS_MOLOPT \
     --basis_name DZVP-MOLOPT-SR-GTH \
     --material INAS \
     --e_thresh 2.0 \
     --qp_gap gw \
     --eps-out 2.4 \
     --plot \
     --nthreads 8

Key CLI Arguments:
* `--mo_file`: Path to binary `.mbse` or text MO file.
* `--xyz`: Cartesian coordinates of the system in standard XYZ format.
* `--basis_txt` & `--basis_name`: CP2K MOLOPT Gaussian basis set library and specific basis set name.
* `--material`: Built-in material identifier (e.g., `INAS`, `CSPBBR3`, `CDSE`) for dielectric screening.
* `--qp_gap`: Target quasiparticle gap or model (`gw`, `brus`, or numeric value in eV).
* `--plot`: Generates UV-Vis absorption spectra and interactive Plotly dashboards.

2. YAML Configuration Workflow
------------------------------

For reproducible and complex calculations, using a YAML file is recommended:

.. code-block:: yaml

   # config.yaml
   system:
     xyz_file: "structure.xyz"
     mo_file: "MOs.mbse"
     basis_txt: "BASIS_MOLOPT"
     basis_name: "DZVP-MOLOPT-PBE-GTH"
     material: "CSPBBR3"
     nthreads: 12

   physics:
     excitation_mode: "diagonal_bse"
     qp_gap: "gw"
     kernel: "resta"
     eps_out: 2.4
     nhomos: 1275
     nlumos: 522
     soc: true
     gth_file: "GTH_SOC_POTENTIALS.txt"

Run with:

.. code-block:: bash

   qdex --config config.yaml

3. NAMD Carrier Cooling Workflow
--------------------------------

For Non-Adiabatic Molecular Dynamics (NAMD) across an MD trajectory:

Step 1: Precompute overlaps and exciton states
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   qdex --config config.yaml --namd-precompute

This step computes:
* Quasiparticle energies and diagonal BSE exciton states.
* Cross-frame non-adiabatic overlaps :math:`S(t, t+\Delta t)`.
* Spinor phase alignment and Hungarian crossing tracking.
* Compact caching into `.npz` step files.

Step 2: Run Dynamics & Carrier Cooling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   qdex --config config.yaml --namd-run

Propagates the carrier cooling cascade deterministically via the Pauli Master Equation or stochastically via CPA-FSSH or DISH, exporting cooling curves, state populations, and photoluminescence summaries.

