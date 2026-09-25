> Archived output of `benchmarks/compare_models.py --profile full --soc --exp-gap 2.70 2.95`
> on the 2.0 nm cluster (inputs now in `tests/CdSe/2.0nm`), 25 September 2026, commit `6e50f97`.
> SOC was ON: the "S1" columns are spin-free, the "1st bright, SOC" column is from the SOC BSE.
> The two FAIL rows predate fixes: `no-exchange` crashed with SOC (fixed in `ed3f120`), and
> `as:50x50` was killed for memory with dense SOC diagonalization (now run with Davidson).

# QDEX model comparison: /home/user/QDEX/tests/CdSe

profile `full`, eps_solvent = 2.24, eps_inf = 6.2, bulk GW-PBE opening = 1.270 eV, spin-free

## A. QP models, each with its compatible BSE kernel (vacuum)

| case | QP gap (eV) | S1 (eV) | f(S1) | QP - S1 (eV) | 1st bright, SOC (eV) | wall (s) | status |
|---|---|---|---|---|---|---|---|
| `pbe` | 1.457 | 1.230 | 0.133 | 0.227 | 1.168 | 12.8 | ok |
| `brus` | 3.856 | 3.629 | 0.392 | 0.227 | 3.539 | 11.9 | ok |
| `gw` | 3.861 | 3.634 | 0.393 | 0.227 | 3.544 | 11.8 | ok |
| `sgw-resta` | 3.762 | 2.240 | 0.330 | 1.522 | 2.150 | 14.5 | ok |
| `sgw-resta-pure` | 3.498 | 2.321 | 0.317 | 1.177 | 2.245 | 15.2 | ok |
| `sgw-dim` | 3.593 | 2.312 | 0.335 | 1.281 | 2.232 | 15.4 | ok |
| `evgw-resta` | 3.851 | 2.440 | 0.361 | 1.411 | 2.354 | 14.0 | ok |
| `evgw-dim` | 3.732 | 2.465 | 0.351 | 1.267 | 2.386 | 18.3 | ok |
| `qsgw-resta` | 4.234 | 2.738 | 0.477 | 1.496 | 2.634 | 52.1 | ok |
| `qsgw-dim` | 3.955 | 2.689 | 0.427 | 1.266 | 2.592 | 21.3 | ok |
| `sgw(sbse)` | 4.598 | 3.219 | 0.412 | 1.379 | 3.125 | 15.1 | ok |

## B. Solvent: same model at eps_out = solvent vs. vacuum (group A)

| case | QP gap (eV) | S1 (eV) | f(S1) | QP - S1 (eV) | 1st bright, SOC (eV) | wall (s) | status |
|---|---|---|---|---|---|---|---|
| `gw@2.24` | 3.237 | 3.010 | 0.325 | 0.227 | 2.920 | 11.4 | ok |
| `sgw-resta@2.24` | 3.253 | 2.383 | 0.342 | 0.870 | 2.302 | 13.9 | ok |
| `sgw-dim@2.24` | 3.084 | 2.441 | 0.311 | 0.643 | 2.370 | 13.9 | ok |
| `sgw(sbse)@2.24` | 4.297 | 3.285 | 0.415 | 1.012 | 3.192 | 16.0 | ok |
| `sgw-resta[Z=1]@2.24` | 3.385 | 2.515 | 0.361 | 0.870 | 2.434 | 13.2 | ok |
| `sgw-dim[Z=1]@2.24` | 3.173 | 2.531 | 0.323 | 0.642 | 2.433 | 14.2 | ok |

## C. Quasiparticle weight Z (shared W, vacuum)

| case | QP gap (eV) | S1 (eV) | f(S1) | QP - S1 (eV) | 1st bright, SOC (eV) | wall (s) | status |
|---|---|---|---|---|---|---|---|
| `sgw-resta` | 3.762 | 2.240 | 0.330 | 1.522 | 2.150 | 14.5 | ok |
| `sgw-dim` | 3.593 | 2.312 | 0.335 | 1.281 | 2.232 | 15.4 | ok |
| `sgw-resta[Z=1]` | 4.021 | 2.499 | 0.369 | 1.522 | 2.409 | 13.7 | ok |
| `sgw-resta[Z=derived]` | 3.949 | 2.427 | 0.358 | 1.522 | 2.336 | 13.5 | ok |
| `sgw-dim[Z=1]` | 3.809 | 2.529 | 0.366 | 1.280 | 2.435 | 14.5 | ok |
| `sgw-dim[Z=derived]` | 3.744 | 2.463 | 0.357 | 1.281 | 2.383 | 13.7 | ok |

## D. Independent BSE kernels for the gap-only gw model

| case | QP gap (eV) | S1 (eV) | f(S1) | QP - S1 (eV) | 1st bright, SOC (eV) | wall (s) | status |
|---|---|---|---|---|---|---|---|
| `gw` | 3.861 | 3.634 | 0.393 | 0.227 | 3.544 | 11.8 | ok |
| `gw+dim` | 3.861 | 3.540 | 0.412 | 0.321 | 3.447 | 12.1 | ok |
| `gw+xs-resta` | 3.861 | 3.594 | 0.349 | 0.267 | 3.502 | 154.2 | ok |
| `gw+sbse` | 3.861 | 2.482 | 0.317 | 1.379 | 2.388 | 16.2 | ok |
| `gw+mnok-bare` | 3.861 | 2.184 | 0.240 | 1.677 | 2.087 | 10.9 | ok |

## E. Excitation framework, spin and charge partition (gw + Resta)

| case | QP gap (eV) | S1 (eV) | f(S1) | QP - S1 (eV) | 1st bright, SOC (eV) | wall (s) | status |
|---|---|---|---|---|---|---|---|
| `gw` | 3.861 | 3.634 | 0.393 | 0.227 | 3.544 | 11.8 | ok |
| `mode:independent_dft` | 3.861 | 1.457 | 0.284 | 2.404 | 1.380 | 7.5 | ok |
| `mode:independent_qp` | 3.861 | 3.861 | 0.752 | 0.000 | 3.784 | 7.5 | ok |
| `mode:diagonal_bse` | 3.861 | 3.653 | 0.711 | 0.208 | 3.537 | 7.5 | ok |
| `spin:triplet` | 3.861 | 3.578 | 0.000 | 0.283 | — | 11.8 | ok |
| `charges:lowdin` | 3.861 | 3.629 | 0.414 | 0.232 | 3.538 | 51.8 | ok |
| `no-exchange` | 3.861 | 3.578 | 1.090 | 0.283 | — | 8.6 | FAILED |

## F. Active-space convergence (sgw-resta, shared W, Z = 1)

| case | QP gap (eV) | S1 (eV) | f(S1) | QP - S1 (eV) | 1st bright, SOC (eV) | wall (s) | status |
|---|---|---|---|---|---|---|---|
| `as:50x50` | 4.021 | 2.475 | 0.311 | 1.546 | — | 93.0 | FAILED |
| `as:100x100` | 4.021 | 2.443 | 0.210 | 1.578 | 2.354 | 232.8 | ok |

## L. Legacy mismatched combination (QP W differs from BSE W)

| case | QP gap (eV) | S1 (eV) | f(S1) | QP - S1 (eV) | 1st bright, SOC (eV) | wall (s) | status |
|---|---|---|---|---|---|---|---|
| `legacy sgw-resta+resta` | 3.762 | 3.535 | 0.382 | 0.227 | 3.445 | 14.2 | ok |
| `legacy sgw-resta+resta@2.24` | 3.253 | 3.026 | 0.327 | 0.227 | 2.936 | 13.4 | ok |

## Checks

| status | check | detail |
|---|---|---|
| **FAIL** | no-exchange ran | AttributeError: 'ExcitonHamiltonian' object has no attribute 'W_elec_spinor'. Did you mean: 'n_occ_spinor'? |
| **FAIL** | as:50x50 ran | return code -9; see run.out |
| **WARN** | gw: S1 insensitive to solvent | eps_out 1 -> 2.24: dQP = -0.624 eV, dS1 = -0.624 eV, ratio 1.00  (no shared W: the BSE does not see the solvent term of the QP model) |
| **WARN** | sgw-resta: S1 insensitive to solvent | eps_out 1 -> 2.24: dQP = -0.509 eV, dS1 = +0.143 eV, ratio 0.28  (shared W but Z < 1: the QP shift is scaled by Z, the BSE attraction is not) |
| **WARN** | sgw-dim: S1 insensitive to solvent | eps_out 1 -> 2.24: dQP = -0.509 eV, dS1 = +0.129 eV, ratio 0.25  (shared W but Z < 1: the QP shift is scaled by Z, the BSE attraction is not) |
| **WARN** | sgw(sbse): S1 insensitive to solvent | eps_out 1 -> 2.24: dQP = -0.301 eV, dS1 = +0.067 eV, ratio 0.22  (shared W but Z < 1: the QP shift is scaled by Z, the BSE attraction is not) |
| **WARN** | legacy sgw-resta+resta: S1 insensitive to solvent | eps_out 1 -> 2.24: dQP = -0.509 eV, dS1 = -0.509 eV, ratio 1.00  (no shared W: the BSE does not see the solvent term of the QP model) |
| **WARN** | gw: binding robust to the independent kernel | resta 0.23 eV; dim 0.32 eV; xs-resta 0.27 eV; sbse 1.38 eV; mnok-bare 1.68 eV  (spread 1.45 eV) |
| **WARN** | S1 converged in the active space | 25x25: 2.499 -> 50x50: 2.475 -> 100x100: 2.443 eV (last change -31.9 meV) |
| **WARN** | sgw-resta@2.24: first bright state within experiment [2.7, 2.95] eV | 2.302 eV (-0.40 eV outside) |
| **WARN** | sgw-dim@2.24: first bright state within experiment [2.7, 2.95] eV | 2.370 eV (-0.33 eV outside) |
| **WARN** | sgw(sbse)@2.24: first bright state within experiment [2.7, 2.95] eV | 3.192 eV (+0.24 eV outside) |
| **WARN** | sgw-resta[Z=1]@2.24: first bright state within experiment [2.7, 2.95] eV | 2.434 eV (-0.27 eV outside) |
| **WARN** | sgw-dim[Z=1]@2.24: first bright state within experiment [2.7, 2.95] eV | 2.433 eV (-0.27 eV outside) |
| **WARN** | sgw-resta[Z=1]: first bright state within experiment [2.7, 2.95] eV | 2.409 eV (-0.29 eV outside) |
| **WARN** | sgw-dim[Z=1]: first bright state within experiment [2.7, 2.95] eV | 2.435 eV (-0.26 eV outside) |
| **PASS** | independent_dft lowest = DFT gap | difference -0.0 meV |
| **PASS** | independent_qp lowest = QP gap | difference -0.1 meV |
| **PASS** | full BSE S1 <= lowest diagonal element | diagonal - full = +19.1 meV (variational principle) |
| **PASS** | T1 <= S1 (exchange is repulsive) | S1 - T1 = +56.7 meV |
| **PASS** | singlet without exchange = triplet | difference +0.0 meV |
| **PASS** | sgw-resta[Z=1]: S1 insensitive to solvent | eps_out 1 -> 2.24: dQP = -0.636 eV, dS1 = +0.016 eV, ratio 0.02 |
| **PASS** | sgw-dim[Z=1]: S1 insensitive to solvent | eps_out 1 -> 2.24: dQP = -0.636 eV, dS1 = +0.002 eV, ratio 0.00 |
| **PASS** | sgw-resta: S1 vs Z | [Z=0.8] 2.240 eV, [Z=1] 2.499 eV, [Z=derived] 2.427 eV  (informative: Z changes S1 through the QP shift only) |
| **PASS** | sgw-dim: S1 vs Z | [Z=0.8] 2.312 eV, [Z=1] 2.529 eV, [Z=derived] 2.463 eV  (informative: Z changes S1 through the QP shift only) |
| **PASS** | gw@2.24: first bright state within experiment [2.7, 2.95] eV | 2.920 eV |
| **PASS** | legacy sgw-resta+resta@2.24: first bright state within experiment [2.7, 2.95] eV | 2.936 eV |

How to read this report: `docs/validation/model_comparison.rst`.
