# Bulk reference band structures and bulk PBE gaps

Spin-free and spin-orbit bulk bands of the QDSpaceWebApp materials for the fuzzy-band dashboards
(`qdex/data/bulk_bands`) and the bulk PBE gaps of `MATERIAL_DB` (`qdex/hardness.py`), at the level of
theory of the dot calculations (CP2K 2026.2, PBE, DZVP-MOLOPT-PBE-GTH, GTH-PBE of POTENTIAL_UZH).
See `docs/electronic_structure/bulk_bands.rst`.

1. `python make_cp2k_inputs.py CIF_ROOT OUT_DIR [BASIS_DIR]` -- one input per CIF at the CIF lattice
   (`inputs/`).
2. `python make_lattice_scan.py OUT_DIR` -- each input scaled to the experimental lattice
   (MATERIAL_DB index 2) and to 0.97 ... 1.03 of the CIF lattice (energy-volume scan).
   `python analyse_lattice_scan.py RUNS_DIR BASIS GTH` -- Birch-Murnaghan PBE lattice constants and
   the gaps at both lattices (`lattice_scan.json`).
3. `python make_pbe_lattice_inputs.py OUT_DIR` -- band-structure inputs at the PBE lattice
   (`inputs_pbe_lattice/`, the ones behind qdex/data/bulk_bands).
4. Run CP2K in every folder, then
   `./process_all.sh RUNS_DIR CIF_ROOT BASIS_MOLOPT_UZH GTH_SOC_POTENTIALS` rebuilds the spin-free
   bands from the real-space matrices (checked against CP2K's band structure), adds spin-orbit
   coupling with the QDEX GTH-SOC operator and writes the band files and metadata.

5. Bulk QP reference (`MATERIAL_DB` index 8, `MATERIAL_DB_BULK_PBE`): `python make_exp_lattice_inputs.py OUT_DIR`
   (inputs at the experimental lattice with the real-space matrices), run CP2K and `qdex.bulk_soc` (as in 4,
   with `-d` to a scratch folder), `python summarize_soc_gaps.py JSON_DIR exp_lattice_soc_gaps.json` (signed
   spin-free and PBE+SOC gaps; `pbe_lattice_soc_gaps.json` is the same for `qdex/data/bulk_bands`), then
   `python build_qsgw_reference.py pbe_lattice_soc_gaps.json`: Delta_Sigma = literature QSGW+SOC
   (`QSGW_SOC_LITERATURE`) - PBE+SOC at a_exp, and the table entries.

6. Room-temperature structures (`inputs_rt_structures/`, gaps in `rt_structure_soc_gaps.json`): CsPbBr3 Pnma
   (Stoumpos et al. 2013, COD 4510745), gamma-CsPbI3 Pnma (Straus et al. 2019, COD 4127359), CsPbCl3 Pnma at
   the measured volume with PBE-relaxed tilts (`CsPbCl3_ortho_relax`: fixed-cell GEO_OPT; `CsPbCl3_ortho`:
   bands of the relaxed structure), generated with `make_cp2k_inputs.py`; `qdex.bulk_soc` gives the PBE and PBE+SOC gaps of the measured phase
   (`gap_exp_structure`, `gap_soc_exp_structure` in `MATERIAL_DB_BULK_PBE`), the partner of the experimental
   gap in `bulk_vertex_factor: material`.

7. Wurtzite CdSe and CdS (`MATERIAL_DB` CDSE_WZ, CDS_WZ): `inputs/CdS_wz` (COD 9011663),
   `inputs_exp_lattice_wz/` (room-temperature lattices), the CdS scan with `make_lattice_scan.scaled_input`
   (a_PBE 4.215 A) and `inputs_pbe_lattice/CdS_wz` (overlay bands). Run on NHR (CP2K 2026.1, no TREXIO:
   `qdex.bulk_soc` reads the structure from the input).

8. Unfolded overlay of the tilted perovskites (`qdex/data/bulk_bands/CsPbX3_cubic_unfolded.npz`): the Pnma
   inputs of step 6 scaled to the PBE pseudo-cubic volume, run, then
   `python -m qdex.bulk_unfold RUN_DIR --basis BASIS --gth GTH --path-bs qdex/data/bulk_bands/CsPbX3_cubic.bs.gz
   -o qdex/data/bulk_bands/CsPbX3_cubic_unfolded.npz` (validated on untilted cubic CsPbBr3 in the same supercell).

`gw_test/`: CP2K periodic G0W0 test for bulk CdSe (sets up, runs out of memory on 191 GB nodes).
