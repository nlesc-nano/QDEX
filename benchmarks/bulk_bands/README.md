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

`gw_test/`: CP2K periodic G0W0 test for bulk CdSe (sets up, runs out of memory on 191 GB nodes).
