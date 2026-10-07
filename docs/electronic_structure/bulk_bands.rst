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
S/Se/Te q6, Cl/Br/I q7), with a Monkhorst-Pack grid of about 0.22 :math:`2\pi`/Å spacing
(8x8x8 for zinc-blende CdSe; 300 K Fermi smearing for InSb, GaAs and GaSb, gapless or inverted
in PBE at this lattice), at the PBE equilibrium lattice constant of the same level of theory,
since the dots are relaxed with PBE as well. The structure is the CIF of the QDSpaceWebApp scaled
isotropically to the minimum of an energy-volume scan (five lattices from -3 to +3 % of the CIF,
third-order Birch-Murnaghan fit; wurtzite keeps c/a and u of the CIF). The PBE lattice constants
are 0.4-2.3 % larger than experiment, as usual for PBE (CdSe 6.217 against 6.05 Å, GaAs 5.774
against 5.65 Å). A cell optimization with the analytic stress tensor was tried for CdSe: with
k-points it breaks the cubic symmetry of zinc blende, so the energy-volume scan is used.
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

For the zinc-blende semimetals of PBE (HgS, HgSe, HgTe, InAs, InSb; GaSb at its PBE lattice)
the gap is zero; the table gives the s-p band order :math:`E(\Gamma_1) - E(\Gamma_{15})`
(spin-free) and :math:`E(\Gamma_6) - E(\Gamma_8)` (with SOC), negative when inverted. Their
spin-orbit splittings are well defined and shown in the overlay.

Bulk PBE gaps for the quasiparticle shift
------------------------------------------

The bulk QP shift of the dots is :math:`\Delta_\text{bulk} = E_g^\text{QSGW} - E_g^\text{PBE}`
(``MATERIAL_DB`` index 8 minus index 7). The QSGW gaps are literature values at the experimental
lattice constant; QSGW is self-consistent, so they do not depend on the DFT starting point, basis
or pseudopotential. The PBE gap has to be computed at the level of theory of the dots and at the
same lattice as the QSGW gap: gaps change with volume by several eV per unit strain, and a PBE gap
at the PBE lattice would add a strain effect to :math:`\Delta_\text{bulk}` (GaAs: 0.50 eV at the
experimental lattice, 0.02 eV at its 2 % larger PBE lattice). ``MATERIAL_DB`` index 7 is therefore
the spin-free PBE gap at the experimental lattice (``MATERIAL_DB`` index 2), computed with the same
CP2K settings; :math:`\Delta_\text{bulk}` is then a self-energy correction at fixed geometry,
which transfers to the dots and their PBE-relaxed geometries. ``MATERIAL_DB_BULK_PBE`` (next to
``MATERIAL_DB`` in ``qdex/hardness.py``) keeps the gaps at both lattices; negative values are the
inverted s-p order.

======== ========= ================ ========= ================ ======== ======
Material a_exp (Å) E_g PBE at a_exp a_PBE (Å) E_g PBE at a_PBE E_g QSGW Δ_bulk
======== ========= ================ ========= ================ ======== ======
CSPBCL3  5.600     +1.771           5.759     +2.113           4.20     2.43  
CSPBBR3  5.830     +1.278           6.024     +1.671           3.40     2.12  
CSPBI3   6.200     +0.832           6.417     +1.207           2.80     1.97  
ZNS      5.410     +2.051           5.424     +2.018           4.17     2.12  
ZNSE     5.670     +1.252           5.752     +1.085           3.26     2.01  
ZNTE     6.100     +1.225           6.206     +0.984           2.96     1.73  
CDS      5.820     +1.146           5.947     +1.004           2.78     1.63  
CDSE     6.050     +0.644           6.217     +0.474           2.19     1.55  
CDTE     6.480     +0.740           6.630     +0.530           2.12     1.38  
HGS      5.850     -0.421           6.013     -0.534           0.31     0.73  
HGSE     6.080     -0.872           6.285     -1.004           0.05     0.92  
HGTE     6.460     -0.665           6.679     -0.908           0.23     0.90  
ALP      5.460     +1.667           5.511     +1.724           2.74     1.07  
ALAS     5.660     +1.522           5.738     +1.592           2.46     0.94  
ALSB     6.140     +1.208           6.243     +1.212           1.80     0.59  
GAP      5.450     +1.632           5.433     +1.614           2.49     0.86  
GAAS     5.650     +0.496           5.774     +0.022           1.89     1.39  
GASB     6.100     +0.156           6.270     -0.425           1.20     1.04  
INP      5.870     +0.680           5.978     +0.377           1.65     0.97  
INAS     6.060     -0.247           6.204     -0.603           0.80     1.05  
INSB     6.480     -0.153           6.661     -0.622           0.77     0.92  
PBS      5.940     +0.253           6.040     +0.421           0.73     0.48  
PBSE     6.120     +0.139           6.241     +0.325           0.65     0.51  
======== ========= ================ ========= ================ ======== ======

The gaps at the experimental lattice are close to the bulk PBE values QDEX used before (CdSe
0.62, GaAs 0.49 eV), which were also experimental-lattice PBE gaps.

Results
-------

============= =========== ========= ============ ================ ========== ========================= ========== ======
Material      Space group a_PBE (Å) E_g PBE (eV) E_g PBE+SOC (eV) Edges (SF) Γ s-p order SF / SOC (eV) Δ_so (eV)  Anchor
============= =========== ========= ============ ================ ========== ========================= ========== ======
AlAs_zb       F-43m       5.738     1.592        1.495            Γ-X        +1.849 / +1.752           0.293      As s  
AlP_zb        F-43m       5.511     1.724        1.704            Γ-X        +3.232 / +3.212           0.059      P s   
AlSb_zb       F-43m       6.243     1.212        1.004            Γ-L        +1.320 / +1.111           0.647      Sb s  
CdS_zb        F-43m       5.947     1.005        0.986            Γ-Γ        +1.005 / +0.987           0.049      Cd d  
CdSe_wz       P6_3mc      4.399     0.522        0.401            Γ-Γ        --                        --         Cd d  
CdSe_zb       F-43m       6.217     0.474        0.352            Γ-Γ        +0.474 / +0.353           0.362      Cd d  
CdTe_zb       F-43m       6.630     0.530        0.245            Γ-Γ        +0.530 / +0.245           0.871      Cd d  
CsPbBr3_cubic Pm-3m       6.024     1.671        0.738            R-R        --                        1.324 (CB) Cs p  
CsPbCl3_cubic Pm-3m       5.759     2.113        1.211            R-R        --                        1.319 (CB) Cs p  
CsPbI3_cubic  Pm-3m       6.417     1.207        0.155            R-R        --                        1.290 (CB) Cs p  
GaAs_zb       F-43m       5.774     0.022        0.000            Γ-Γ        +0.022 / -0.085           0.324      Ga d  
GaP_zb        F-43m       5.433     1.614        1.586            Γ-path     +1.866 / +1.838           0.085      Ga d  
GaSb_zb       F-43m       6.270     0.000        0.002            Γ-Γ        -0.425 / -0.646           0.679      Ga d  
HgS_zb        F-43m       6.013     0.000        0.036            Γ-Γ        -0.534 / -0.570           0.047      Hg d  
HgSe_zb       F-43m       6.285     0.000        -0.003           Γ-Γ        -1.003 / -1.102           0.239      Hg d  
HgTe_zb       F-43m       6.679     0.000        -0.002           Γ-Γ        -0.908 / -1.179           0.776      Hg d  
InAs_zb       F-43m       6.204     0.000        0.001            Γ-Γ        -0.603 / -0.713           0.334      In d  
InP_zb        F-43m       5.978     0.376        0.346            Γ-Γ        +0.377 / +0.346           0.096      In d  
InSb_zb       F-43m       6.661     0.000        0.000            Γ-Γ        -0.622 / -0.854           0.706      In d  
PbS_rs        Fm-3m       6.040     0.421        0.134            L-L        --                        --         S s   
PbSe_rs       Fm-3m       6.241     0.325        0.007            L-L        --                        --         Se s  
PbTe_rs       Fm-3m       6.584     0.616        0.012            L-L        --                        --         Te s  
ZnS_zb        F-43m       5.424     2.018        1.997            Γ-Γ        +2.018 / +1.997           0.061      Zn d  
ZnSe_zb       F-43m       5.752     1.085        0.958            Γ-Γ        +1.085 / +0.958           0.380      Zn d  
ZnTe_zb       F-43m       6.206     0.984        0.689            Γ-Γ        +0.984 / +0.689           0.902      Zn d  
============= =========== ========= ============ ================ ========== ========================= ========== ======

:math:`\Delta_{so}` is the :math:`\Gamma_8`-:math:`\Gamma_7` splitting of the valence-band
triplet at :math:`\Gamma` (zinc blende) or the :math:`j=3/2`-:math:`j=1/2` splitting of the
Pb 6p conduction-band triplet at R (perovskites, "CB"). The zinc-blende splittings follow the
measured ones (AlP 0.06, GaAs 0.32, CdSe 0.36, GaSb 0.68, CdTe 0.87 eV, against about 0.07,
0.34, 0.42, 0.75, 0.95 eV). At its PBE lattice GaSb is inverted as well.

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

``benchmarks/bulk_bands`` holds the input generators, the CP2K inputs (CIF lattice, lattice scan,
PBE lattice) and the analysis scripts (``README.md`` there). The CP2K input of a material prints the band path, the real-space matrices and a TREXIO file:

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
