# Bulk reference band structures

Spin-free and spin-orbit bulk bands of the QDSpaceWebApp materials, used by the fuzzy-band
dashboards (`qdex/data/bulk_bands`, see `docs/electronic_structure/bulk_bands.rst`).

1. `python make_cp2k_inputs.py CIF_ROOT OUT_DIR [BASIS_DIR]` writes one CP2K input per CIF
   (`inputs/` holds the inputs used, CP2K 2026.2; each runs in 15 s to 2 min on 48 cores).
2. Run CP2K in every folder (`cp2k.popt -i cp2k_job.in -o cp2k_job.out`).
3. `./process_all.sh RUNS_DIR CIF_ROOT BASIS_MOLOPT_UZH GTH_SOC_POTENTIALS` rebuilds the
   spin-free bands from the real-space matrices (checked against CP2K's band structure), adds
   spin-orbit coupling with the QDEX GTH-SOC operator and writes the band files and metadata.
