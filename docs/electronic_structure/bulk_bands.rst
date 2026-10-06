Bulk reference band structures
==============================

Part of :doc:`/electronic_structure/index`.

The fuzzy-band dashboards draw the bulk band structure of the material over the quantum-dot
map. ``qdex/data/bulk_bands`` holds, for every bulk CIF of the QDSpaceWebApp (II-VI, III-V,
IV-VI and the cubic lead-halide perovskites):

* ``<name>.bs.gz`` -- spin-free PBE bands along the pymatgen high-symmetry path of the CIF
  (the path of the fuzzy bands);
* ``<name>_soc.bs.gz`` -- the same bands with spin-orbit coupling;
* ``<name>.json`` -- band edges and gaps, band characters at :math:`\Gamma`, the s-p band order
  and spin-orbit splittings, and the semicore manifold that places the bulk bands on a dot's
  energy axis.

``<name>`` is the CIF stem (``CdSe_zb``, ``CdSe_wz``, ``PbSe_rs``, ``CsPbBr3_cubic``, ...). A run
uses the files of its ``analysis.cif``; otherwise those of ``system.material`` (zinc blende,
rock salt or cubic first).

Level of theory
---------------

The bands are computed at the level of theory of the dot calculations: CP2K PBE, GPW (400 Ry,
EPS_DEFAULT 1e-12), DZVP-MOLOPT-PBE-GTH basis sets and GTH-PBE pseudopotentials (POTENTIAL_UZH)
with the valence charges QDEX uses (Cd/Zn/Hg q12, Ga/In q13, Pb q4, Cs q9, Al q3, P/As/Sb q5,
S/Se/Te q6, Cl/Br/I q7), at the lattice of the CIF, with a Monkhorst-Pack grid of about
0.22 :math:`2\pi`/Å spacing (8x8x8 for zinc-blende CdSe; InSb with 300 K Fermi smearing).
CP2K writes the real-space Kohn-Sham and overlap matrices, and ``python -m qdex.bulk_soc``
rebuilds the spin-free bands (they agree with CP2K's own band structure to 0.3 meV or
better) and adds the QDEX GTH-SOC operator in the full AO basis, as for the dots
(:doc:`fuzzy_bands`). The SOC constants are those of GTH_SOC_POTENTIALS (PBE). For Pb the
scalar part of the q4 pseudopotential in POTENTIAL_UZH differs from the one the SOC constants
were fitted with (the same approximation as in the dot calculations); with the GTH_SOC_POTENTIALS
scalar part instead, the PBE+SOC gap of CsPbBr3 is 0.620 instead of 0.602 eV.

Validation against CP2K's own spin-orbit band structure (``&PROPERTIES &BANDSTRUCTURE &SOC``,
CP2K 2026.2): on the same SCF (GTH_SOC_POTENTIALS), the QDEX operator reproduces every spinor
band at 81 k-points of the path within 0.5 meV (CP2K prints 1 meV) for zinc-blende CdSe (76
bands, PBE+SOC gap 0.425 eV) and cubic CsPbBr3 (132 bands, 0.620 eV).

The spin-free gaps are the bulk PBE gaps of ``MATERIAL_DB`` (index 7). For the zinc-blende
semimetals of PBE (HgS, HgSe, HgTe, InAs, InSb) the gap is zero and the table and the database
give the s-p band order :math:`E(\Gamma_1) - E(\Gamma_{15})` (spin-free) and
:math:`E(\Gamma_6) - E(\Gamma_8)` (with SOC), negative when inverted; their spin-orbit
splittings are well defined and shown in the overlay.

Results
-------

============= =========== ============ ================ ========== ========================= ========== ======
Material      Space group E_g PBE (eV) E_g PBE+SOC (eV) Edges (SF) Γ s-p order SF / SOC (eV) Δ_so (eV)  Anchor
============= =========== ============ ================ ========== ========================= ========== ======
AlAs_zb       F-43m       1.541        1.443            Γ-X        +2.101 / +2.003           0.297      As s  
AlP_zb        F-43m       1.678        1.658            Γ-X        +3.437 / +3.417           0.059      P s   
AlSb_zb       F-43m       1.245        1.033            Γ-path     +1.526 / +1.314           0.654      Sb s  
CdS_zb        F-43m       1.073        1.056            Γ-Γ        +1.073 / +1.056           0.046      Cd d  
CdSe_wz       P6_3mc      0.603        0.482            Γ-Γ        --                        --         Cd d  
CdSe_zb       F-43m       0.551        0.429            Γ-Γ        +0.551 / +0.430           0.361      Cd d  
CdTe_zb       F-43m       0.621        0.334            Γ-Γ        +0.621 / +0.334           0.875      Cd d  
CsPbBr3_cubic Pm-3m       1.533        0.602            R-R        --                        1.327 (CB) Cs p  
CsPbCl3_cubic Pm-3m       1.955        1.053            R-R        --                        1.317 (CB) Cs p  
CsPbI3_cubic  Pm-3m       0.977        0.052            R-R        --                        1.300 (CB) Cs p  
GaAs_zb       F-43m       0.110        0.003            Γ-Γ        +0.110 / +0.003           0.325      Ga d  
GaP_zb        F-43m       1.634        1.606            Γ-path     +1.780 / +1.752           0.085      Ga d  
GaSb_zb       F-43m       0.024        0.001            Γ-Γ        +0.024 / -0.203           0.696      Ga d  
HgS_zb        F-43m       0.000        0.040            Γ-Γ        -0.483 / -0.521           0.052      Hg d  
HgSe_zb       F-43m       0.000        -0.004           Γ-Γ        -0.946 / -1.042           0.230      Hg d  
HgTe_zb       F-43m       0.000        -0.002           Γ-Γ        -0.803 / -1.074           0.776      Hg d  
InAs_zb       F-43m       0.000        0.001            Γ-Γ        -0.365 / -0.477           0.340      In d  
InP_zb        F-43m       0.583        0.553            Γ-Γ        +0.584 / +0.553           0.097      In d  
InSb_zb       F-43m       0.000        0.000            Γ-Γ        -0.552 / -0.786           0.710      In d  
PbS_rs        Fm-3m       0.317        0.028            L-L        --                        --         S s   
PbSe_rs       Fm-3m       0.238        0.080            L-L        --                        --         Se s  
PbTe_rs       Fm-3m       0.580        0.027            L-L        --                        --         Te s  
ZnS_zb        F-43m       2.105        2.084            Γ-Γ        +2.105 / +2.084           0.060      Zn d  
ZnSe_zb       F-43m       1.263        1.136            Γ-Γ        +1.263 / +1.136           0.382      Zn d  
ZnTe_zb       F-43m       1.199        0.901            Γ-Γ        +1.199 / +0.901           0.911      Zn d  
============= =========== ============ ================ ========== ========================= ========== ======

:math:`\Delta_{so}` is the :math:`\Gamma_8`-:math:`\Gamma_7` splitting of the valence-band
triplet at :math:`\Gamma` (zinc blende) or the :math:`j=3/2`-:math:`j=1/2` splitting of the
Pb 6p conduction-band triplet at R (perovskites, "CB"). The zinc-blende splittings follow the
measured ones (AlP 0.06, GaAs 0.33, CdSe 0.36, GaSb 0.70, CdTe 0.88 eV, against about 0.07,
0.34, 0.42, 0.75, 0.95 eV).

Energy anchor
-------------

The bulk bands are placed on a dot's energy axis with a semicore manifold, measured on the
interior bulk-like atoms of the dot (full crystal coordination, no ligand neighbour, inner half by
radius) as the Mulliken-weighted energy of its l-character. ``<name>.json`` lists the manifolds
in order of preference (``semicore_candidates``): the cation d band where it is in the valence,
Cs 5p in the perovskites, the anion s band. The table shows the first one. A run uses the
first manifold its orbitals reach (MO files often hold only the orbitals near the gap; the InP
example file stops 9 eV below the HOMO, above both In 4d and P 3s) and aligns the bulk mid-gap
at zero when none is reached. The dashboard overlay uses the same manifold.

The anchor assumes the dot and the bulk bands share the level of theory. When the anchored bulk
VBM lies more than 0.2 eV below the dot's HOMO (a confined hole lies below it), the run warns: a
CdS dot computed with another basis set and a larger gap than PBE gives a Cd 4d to HOMO
separation 0.9 eV larger than bulk PBE.

Recomputing
-----------

The CP2K input of a material prints the band path, the real-space matrices and a TREXIO file:

.. code-block:: text

   &PRINT
     &BAND_STRUCTURE
       FILE_NAME CdSe_zb.bs
       &KPOINT_SET ... &END KPOINT_SET          # pymatgen HighSymmKpath of the CIF
     &END BAND_STRUCTURE
     &KS_CSR_WRITE
       REAL_SPACE .TRUE.
       UPPER_TRIANGULAR .FALSE.
     &END KS_CSR_WRITE
     &S_CSR_WRITE
       REAL_SPACE .TRUE.
       UPPER_TRIANGULAR .FALSE.
     &END S_CSR_WRITE
     &TREXIO
     &END TREXIO
   &END PRINT

.. code-block:: bash

   python -m qdex.bulk_soc RUN_DIR --basis BASIS_MOLOPT_UZH --gth GTH_SOC_POTENTIALS \
       --name CdSe_zb --cif CdSe_zb.cif -d qdex/data/bulk_bands --gzip
