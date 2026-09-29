> Archived output of `benchmarks/compare_models.py --profile full --only A` on
> `tests/CdSe/1.2nm` (Cd16Se13Cl6, spin-free, vacuum), 25 September 2026. The evGW comparison
> and the Z = 1 runs are in `docs/validation/anchor_evgw_benchmark.rst` and in audit report §8.

# QDEX model comparison: /home/user/QDEX/tests/CdSe/1.2nm

profile `full`, eps_solvent = 2.24, eps_inf = 6.2, bulk GW-PBE opening = 1.270 eV, spin-free

## A. QP models, each with its compatible BSE kernel (vacuum)

| case | QP gap (eV) | S1 (eV) | f(S1) | QP - S1 (eV) | 1st bright, SOC (eV) | wall (s) | status |
|---|---|---|---|---|---|---|---|
| `pbe` | 2.639 | 2.186 | 0.000 | 0.453 | — | 2.7 | ok |
| `brus` | 6.220 | 5.766 | 0.000 | 0.454 | — | 2.3 | ok |
| `gw` | 6.009 | 5.555 | 0.000 | 0.454 | — | 2.4 | ok |
| `sgw-resta` | 6.215 | 2.609 | 0.000 | 3.606 | — | 2.4 | ok |
| `sgw-resta-pure` | 5.201 | 3.043 | 0.000 | 2.158 | — | 2.4 | ok |
| `sgw-dim` | 5.316 | 3.010 | 0.000 | 2.306 | — | 2.4 | ok |
| `evgw-resta` | 6.338 | 3.026 | 0.000 | 3.312 | — | 2.7 | ok |
| `evgw-dim` | 5.514 | 3.218 | 0.000 | 2.296 | — | 2.9 | ok |
| `qsgw-resta` | 6.726 | 3.395 | 0.385 | 3.331 | — | 3.4 | ok |
| `qsgw-dim` | 5.810 | 3.530 | 0.000 | 2.280 | — | 3.2 | ok |
| `sgw(sbse)` | 6.182 | 3.977 | 0.264 | 2.205 | — | 2.6 | ok |

## C. Quasiparticle weight Z (shared W, vacuum)

| case | QP gap (eV) | S1 (eV) | f(S1) | QP - S1 (eV) | 1st bright, SOC (eV) | wall (s) | status |
|---|---|---|---|---|---|---|---|
| `sgw-resta` | 6.215 | 2.609 | 0.000 | 3.606 | — | 2.4 | ok |
| `sgw-dim` | 5.316 | 3.010 | 0.000 | 2.306 | — | 2.4 | ok |

## D. Independent BSE kernels for the gap-only gw model

| case | QP gap (eV) | S1 (eV) | f(S1) | QP - S1 (eV) | 1st bright, SOC (eV) | wall (s) | status |
|---|---|---|---|---|---|---|---|
| `gw` | 6.009 | 5.555 | 0.000 | 0.454 | — | 2.4 | ok |

## E. Excitation framework, spin and charge partition (gw + Resta)

| case | QP gap (eV) | S1 (eV) | f(S1) | QP - S1 (eV) | 1st bright, SOC (eV) | wall (s) | status |
|---|---|---|---|---|---|---|---|
| `gw` | 6.009 | 5.555 | 0.000 | 0.454 | — | 2.4 | ok |

## Checks

| status | check | detail |
|---|---|---|

How to read this report: `docs/validation/model_comparison.rst`.
