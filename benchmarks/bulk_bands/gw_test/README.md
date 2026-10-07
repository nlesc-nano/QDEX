# CP2K periodic G0W0 for bulk zinc-blende CdSe (test, October 2026)

Can CP2K 2026.2 give a bulk GW correction at the level of theory of the dots (PBE, DZVP-MOLOPT-PBE-GTH)?
`&PROPERTIES &BANDSTRUCTURE &GW` (with `&SOC`) on a k-point SCF:

- the input here sets up and runs: SCF on a 4x4x4 Monkhorst-Pack mesh (`PARALLEL_GROUP_SIZE -1`;
  at least 4 k-points per periodic direction are required), RI basis from BASIS_RI_MOLOPT given
  explicitly in `&KIND` (`BASIS_SET RI_AUX ...`; `AUTO_BASIS RI_AUX` is not accepted by the GW code),
  `SPECIAL_POINT` lines need a label (`SPECIAL_POINT GAMMA 0 0 0`);
- CP2K reports 577 cells / 72,919 cell pairs for the three-center integrals (the diffuse MOLOPT
  functions reach 10.9 A) and 45.5 GB for them;
- every run was killed for memory after 20-50 min in the chi/W stage: 48 x 2 GB, 12 x 7.5 GB,
  4 x 45 GB (one 191 GB node), 8 ranks on two 191 GB nodes; with the default W mesh (16^3, 32^3
  extrapolation) and with `KPOINTS_W 6 6 6`.

The CP2K documentation lists 3D periodic GW as work in progress (2D is supported). Next steps, if
wanted: the 1-2 TB nodes (topmat, hitz-exclusive), a more compact basis (SZV/DZVP without the most
diffuse functions) or a later CP2K. The CP2K DFT+SOC band structure (`&BANDSTRUCTURE &SOC` without
`&GW`) works and agrees with qdex.bulk_soc to 0.5 meV.
