# Post-fix verification and documentation migration

25 September 2026. This report checks the current working tree against the
earlier `SCIENTIFIC_AUDIT.md`. It is a code-path verification, not a validation
of predictive accuracy for a physical CdSe nanocrystal. The user's eight CdSe
audit observations came from an external example folder supplied after the
initial verification. Its files were inspected read-only; the large AO
calculation was not rerun for this report.

## Verified repairs

| Earlier concern | Current evidence | Result |
|---|---|---|
| Complex frontier populations | `hardness.py` uses absolute-square populations; synthetic global phase gives identical 3.0403134975 eV shift. | Repaired for the tested QP path. |
| Polar solvent fraction singularity | Signed exterior reaction shift is separated from allocation. In the synthetic Resta probe at exterior permittivity 10, the scissor is 0.99311965 eV and both fractions are 0.5. | No fraction divergence. |
| Zero-contrast bulk offset | `qsgw-resta` with zero screening contrast returns the 1.27 eV CdSe database bulk opening. | Repaired in the probe. |
| Relaxed QP eigenvalues lost at CLI boundary | `cli.py` retains `eps_qp_active`, passes it as `eps_shifted`, sets solver scissor to zero, and supplies DFT eigenvalues separately for filtering. | Wiring repaired. A full physical QP/BSE benchmark remains separate. |
| Direct and exchange switches coupled | CLI, solver, and Hamiltonian now carry independent flags. This review found matrix-free `kernel_actions` still applied exchange when disabled; that CPU/GPU and spinor/nonspinor gate was repaired and regression-tested. | Repaired in covered paths. |
| Diagonal versus full Löwdin charges | This review aligned the full atom-resolved Löwdin charge builder with diagonal mode. The synthetic full-matrix diagonal and diagonal-mode energies now agree within 9e-16 eV for direct enabled and disabled. | Repaired in the tested restricted-space case. |
| Incomplete MO input to orbital relaxation | A rectangular 3×2 MO matrix now raises an informative `ValueError` instead of a broadcast failure. | Input failure clarified. |

The repository's `unittest` discovery ran **100 tests, all passing** after
these changes, including the new matrix-free switch regression. The synthetic
checks are in `reproduce_findings.py` and current results in
`probe_results.json`. The now-supplied CdSe folder contains a 149-atom
Cd68Se55Cl26 structure, 2,753-AO calculation records and the audit script.
The saved result file reports 0.2774 eV diagonal and 0.4756 eV coupled BSE
binding for its chosen setup. It is not an independent GW/BSE benchmark.
The script has no eight pass/fail assertions: its solvent-pathology result is
a stale fixed string, and its `confirmed_discarded_in_cli` value is hard-coded
`true` even though the current CLI passes relaxed eigenvalues. See
`docs/validation/cdse_benchmark.rst` for provenance and comparison details.

## Documentation delivered

The former eight sequential chapters are now **86 migrated topic pages plus
three new validation/style pages in 12 groups**,
with implementation entry points, signatures, CLI flags and YAML pointers
adjacent to theory. The original chapter bytes are preserved in
`docs/_legacy_parts`; their SHA-256 hashes match `documentation_inventory.json`.
The 255 migration units partition all **5,112 original chapter lines** without
gaps or overlaps. Historical chapter URLs remain as stubs linking to the new
topics and original source download. The Sphinx HTML build passes with
warnings treated as errors. A rendered-HTML scan found and corrected dozens of
literal `:math:` fragments caused by math roles nested inside strong or
emphasis markup; no literal `:math:` remains in the reader-facing HTML pages.

The topical pages correct several historical overclaims while the byte-exact
archive retains their provenance: `xs` is restricted to analytical
`(mu mu|nu nu)` density-pair integrals; `dynamic_z` is an empirical damping
formula, not an f-sum-rule calculation; `qsgw-*` is a static COHSEX-like
orbital-relaxation model; and an unchanged BSE direct kernel is not proof of
universal solvent cancellation. See `docs/validation/environment_cancellation.rst`
and `docs/validation/cdse_benchmark.rst` for the physical conditions and
benchmark provenance.

## Scientific limits still open

The fixes close the listed implementation defects but do not turn the scaled
models into conventional GW/QSGW or a full four-index BSE. A quantitative
environment calculation still needs a consistent charged/neutral reaction
operator, and the CdSe target ranges need a geometry and reference data set.
Diagonal BSE still omits off-diagonal transition mixing; its binding error is
system- and active-space-dependent. The original audit's separate model-level
concerns remain applicable unless explicitly addressed above.
