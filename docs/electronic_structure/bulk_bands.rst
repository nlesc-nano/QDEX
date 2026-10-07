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

Bulk quasiparticle reference
----------------------------

The bulk QP shift of a dot computed with PBE has two parts,

.. math::

   \Delta_\text{bulk} = \Delta_\Sigma + \Delta_\text{geom},\qquad
   \Delta_\Sigma = E_g^\text{QSGW+SOC}(a_\text{exp}) - E_g^\text{PBE+SOC}(a_\text{exp}),\qquad
   \Delta_\text{geom} = E_g^\text{PBE}(a_\text{exp}) - E_g^\text{PBE}(a_\text{dot}).

**Self-energy part.** :math:`\Delta_\Sigma` compares QSGW and PBE at the same lattice and with the
same spin-orbit treatment. The QSGW gaps are literature values (100 % QSGW, RPA screening) that
include spin-orbit coupling, added to the converged QSGW Hamiltonian; QSGW is self-consistent, so they
do not depend on a DFT starting point, basis or pseudopotential. The PBE+SOC gaps are ours, at the level
of theory of the dots (the CP2K settings above and the QDEX GTH-SOC operator), at the experimental
lattice constant ``MATERIAL_DB`` index 2 (inputs: ``benchmarks/bulk_bands/make_exp_lattice_inputs.py``).
Gaps change with volume by several eV per unit strain, so both sides must be at the same lattice (GaAs
PBE: 0.50 eV at :math:`a_\text{exp}`, 0.02 eV at its 2 % larger PBE lattice). ``MATERIAL_DB`` index 7
is the spin-free PBE gap at :math:`a_\text{exp}` and index 8 the spin-free equivalent
:math:`E_g^\text{PBE}(a_\text{exp}) + \Delta_\Sigma`, so index 8 − index 7 = :math:`\Delta_\Sigma` corrects
spin-free and SOC dots alike, and the spin-orbit lowering of the gap is the PBE/GTH-SOC one in the bulk
reference and in the dots.

Gaps are signed band-edge transitions (negative = inverted): :math:`E(\Gamma_6) - E(\Gamma_8)` for the
direct zinc blendes and the Hg compounds (:math:`E(\Gamma_1) - E(\Gamma_{15})` spin-free), the
fundamental gap (X) for the indirect AlP, AlAs, AlSb and GaP, L in PbS and PbSe, R in the
perovskites. PBE+SOC inverts PbS, PbSe and CsPbI\ :sub:`3` at :math:`a_\text{exp}` (as LDA/PBE+SOC is
known to do for the lead chalcogenides); there the gap is signed by the character of the edge states
(``edge_band_order`` in the ``qdex.bulk_soc`` metadata), since counting bands only gives the size of
the anticrossing. HgS has a negative :math:`\Delta_\text{so}` (:math:`\Gamma_7` above :math:`\Gamma_8`,
-0.115 eV; LDA -0.11 eV in Svane et al.).

Sources (``QSGW_SOC_LITERATURE`` in ``qdex/hardness.py``): Deguchi16 = D. Deguchi, K. Sato, H. Kino,
T. Kotani, *Jpn. J. Appl. Phys.* 55, 051201 (2016), Table II QSGW+SO (GaSb: :math:`\Gamma_6` of
Table III; their minimum gap is at L); Svane11 = A. Svane et al., *Phys. Rev. B* 84, 205205 (2011),
Table I (QSGW :math:`E_0`); Huang16 = L.-y. Huang and W. R. L. Lambrecht, *Phys. Rev. B* 93, 195211
(2016), Table II, at the cubic lattices 6.289 / 5.874 / 5.605 Å, carried to :math:`a_\text{exp}` with
their QSGW :math:`dE_g/d\ln V` (6.4 / 7.5 / 7.7 eV); Svane10 = A. Svane et al., *Phys. Rev. B* 81,
245120 (2010), Table I (QSGW, SOC added after self-consistency) at their low-temperature lattices
5.909 / 6.098 Å, carried to :math:`a_\text{exp}` with their QSGW :math:`dE_g/d\ln V` (Table V: 5.3 /
4.9 eV); Deguchi16 gives 0.49 eV for PbS at 5.936 Å, about 0.1 eV higher. The other sources used
lattices within 0.2 % of :math:`a_\text{exp}` (carried with the PBE deformation potential, ≤ 0.01 eV).
``old Δ`` is :math:`\Delta_\text{bulk}` before this table: the III–V values
were the spin-free QSGW of Deguchi16 (hence nearly unchanged), while the II–VI, Hg and perovskite
values matched no QSGW source (the perovskite ones were 0.5–0.6 eV above Huang16). Energies in eV,
lattices in Å; table built by ``benchmarks/bulk_bands/build_qsgw_reference.py``.

======== ===== ====== ======= ========================= ================= ===== ======= ===== ===== ============ ============
Material a_exp PBE    PBE+SOC QSGW+SOC lit. (a, source) QSGW+SOC at a_exp Δ_Σ   index 8 old Δ a_PBE PBE at a_PBE Δ_geom (s=1)
======== ===== ====== ======= ========================= ================= ===== ======= ===== ===== ============ ============
CSPBCL3  5.600 +1.771 +0.870  +2.678 (5.605, Huang16)   +2.657            1.787 +3.558  2.43  5.759 +2.113       -0.342
CSPBBR3  5.830 +1.278 +0.349  +1.868 (5.874, Huang16)   +1.699            1.350 +2.628  2.12  6.024 +1.671       -0.393
CSPBI3   6.200 +0.832 -0.187  +1.331 (6.289, Huang16)   +1.057            1.244 +2.076  1.97  6.417 +1.207       -0.375
ZNS      5.410 +2.051 +2.030  +4.100 (5.413, Deguchi16) +4.107            2.077 +4.128  2.12  5.424 +2.018       +0.033
ZNSE     5.670 +1.252 +1.125  +3.100 (5.667, Deguchi16) +3.094            1.969 +3.221  2.01  5.752 +1.085       +0.167
ZNTE     6.100 +1.225 +0.927  +2.640 (6.101, Deguchi16) +2.642            1.715 +2.940  1.73  6.206 +0.984       +0.241
CDS      5.820 +1.146 +1.130  +2.840 (5.826, Deguchi16) +2.847            1.717 +2.863  1.63  5.947 +1.004       +0.142
CDSE     6.050 +0.644 +0.522  +2.160 (6.054, Deguchi16) +2.164            1.642 +2.286  1.55  6.217 +0.474       +0.170
CDTE     6.480 +0.740 +0.452  +1.970 (6.482, Deguchi16) +1.973            1.521 +2.261  1.38  6.630 +0.530       +0.210
HGS      5.850 -0.421 -0.403  +0.610 (5.840, Svane11)   +0.603            1.006 +0.585  0.73  6.013 -0.534       +0.113
HGSE     6.080 -0.872 -0.964  -0.110 (6.080, Svane11)   -0.110            0.854 -0.018  0.92  6.285 -1.004       +0.132
HGTE     6.460 -0.665 -0.937  +0.090 (6.470, Svane11)   +0.101            1.038 +0.373  0.90  6.679 -0.908       +0.243
ALP      5.460 +1.667 +1.647  +2.720 (5.467, Deguchi16) +2.712            1.065 +2.732  1.07  5.511 +1.724       -0.057
ALAS     5.660 +1.522 +1.423  +2.360 (5.661, Deguchi16) +2.359            0.936 +2.458  0.94  5.738 +1.592       -0.070
ALSB     6.140 +1.208 +0.994  +1.590 (6.136, Deguchi16) +1.590            0.596 +1.804  0.59  6.243 +1.212       -0.004
GAP      5.450 +1.632 +1.604  +2.460 (5.451, Deguchi16) +2.459            0.855 +2.487  0.86  5.433 +1.614       +0.018
GAAS     5.650 +0.496 +0.386  +1.770 (5.653, Deguchi16) +1.782            1.395 +1.891  1.39  5.774 +0.022       +0.474
GASB     6.100 +0.156 -0.071  +1.090 (6.096, Deguchi16) +1.076            1.148 +1.304  1.04  6.270 -0.425       +0.581
INP      5.870 +0.680 +0.649  +1.620 (5.870, Deguchi16) +1.620            0.971 +1.651  0.97  5.978 +0.377       +0.303
INAS     6.060 -0.247 -0.360  +0.680 (6.058, Deguchi16) +0.675            1.035 +0.788  1.05  6.204 -0.603       +0.356
INSB     6.480 -0.153 -0.394  +0.540 (6.479, Deguchi16) +0.537            0.931 +0.778  0.92  6.661 -0.622       +0.469
PBS      5.940 +0.253 -0.037  +0.310 (5.909, Svane10)   +0.393            0.430 +0.683  0.48  6.040 +0.421       -0.168
PBSE     6.120 +0.139 -0.180  +0.210 (6.098, Svane10)   +0.263            0.443 +0.582  0.51  6.241 +0.325       -0.186
======== ===== ====== ======= ========================= ================= ===== ======= ===== ===== ============ ============

**Geometry part.** A dot relaxed with PBE has (almost) the PBE lattice, and its PBE gap carries the
bulk PBE gap change between :math:`a_\text{exp}` and its own lattice. :math:`\Delta_\text{geom}`
removes it, so that the QP gap is that of the dot at the experimental lattice, the one measured. No
QSGW at :math:`a_\text{PBE}` is needed:
:math:`\Delta_\Sigma(a_\text{dot}) + [E^\text{QSGW}(a_\text{exp}) - E^\text{QSGW}(a_\text{dot})] =
\Delta_\Sigma(a_\text{exp}) + \Delta_\text{geom}` exactly. This matters because the self-energy
correction itself depends on the lattice, strongly in the lead halide perovskites (QSGW
:math:`dE_g/d\ln V` about twice the PBE one). The lattice of the dot is measured on its interior
cation–anion bonds (``quasiparticles.bulk_geometry: strain``, default): with the strain fraction
:math:`s = (d_\text{dot}/g - a_\text{exp})/(a_\text{PBE} - a_\text{exp})` (:math:`g = \sqrt{3}/4` in zinc
blende, 1/2 for the rock-salt and the Pb–X bond of the perovskites),
:math:`\Delta_\text{geom} = s\,[E_g^\text{PBE}(a_\text{exp}) - E_g^\text{PBE}(a_\text{PBE})]`. ``full`` takes
:math:`s = 1` and ``none`` :math:`s = 0` (the dot taken at :math:`a_\text{exp}`). The vertex correction
(``quasiparticles.bulk_vertex``) scales :math:`\Delta_\Sigma` only. Every other use of the bulk PBE gap
as the bulk limit of the dot's PBE gap (confinement energy, Penn screening, the anchor models) uses
the gap at the dot's lattice, :math:`E_g^\text{PBE}(a_\text{exp}) - \Delta_\text{geom}`; for the
perovskites the confinement energy is measured from the tilted structure (below).

For the Cd\ :sub:`156`\ Se\ :sub:`111`\ Cl\ :sub:`90` dot (2.62 nm) the interior Cd–Se bond is 2.690 Å
(bulk 2.620 Å at :math:`a_\text{exp}`, 2.692 Å at :math:`a_\text{PBE}`), :math:`s` = 0.975 and
:math:`\Delta_\text{geom}` = +0.166 eV; with ``bulk_vertex: scaled``
:math:`\Delta_\text{bulk}` = 1.379 + 0.166 = 1.545 eV (1.286 eV before), the spin-free QP gap
2.748 eV (2.489) and the first bright SOC exciton 2.499 eV (2.240), against 2.36–2.74 eV from the
sizing curves at this diameter (:doc:`/validation/cdse_experiment`).

Experimental anchor: material vertex factor
-------------------------------------------

QSGW (100 %, RPA screening) misses the electron–hole vertex, the zero-point and the thermal
renormalisation of the gap, and the sum differs between materials. ``bulk_vertex_factor: material``
scales :math:`\Delta_\Sigma` by the factor that puts the bulk limit on the measured gap,
:math:`a_m = (E_g^\text{exp} - E_g^\text{PBE+SOC})/\Delta_\Sigma`, with the room-temperature gap (the
temperature of the sizing experiments) and the PBE+SOC gap of the measured structure. Experimental
gaps: ``MATERIAL_DB`` index 3; zinc-blende CdSe 1.675 eV and wurtzite CdSe 1.751 eV from Aubert et al.,
*Nano Lett.* 22, 1778 (2022) (the bulk gaps of their sizing curves); the Hg compounds the experimental
:math:`E_0` of Svane11; CdS 2.42 eV for both polytypes (the wurtzite value; room-temperature zinc-blende
CdS is not well established, 2.37–2.50 eV in films). "QSGW − exp" is the error of QSGW for the measured
structure; "residual at a = 0.8" is :math:`E_g^\text{exp} - (E_g^\text{PBE+SOC} + 0.8\,\Delta_\Sigma)`, what
a universal factor of 0.8 leaves.

======== ========== ======= ======== ========== ==== ===================
Material E_exp (RT) PBE+SOC QSGW+SOC QSGW − exp a_m  residual at a = 0.8
======== ========== ======= ======== ========== ==== ===================
CSPBCL3  +3.000     +1.549  +3.336   +0.34      0.81 +0.02
CSPBBR3  +2.300     +0.797  +2.147   -0.15      1.11 +0.42
CSPBI3   +1.730     +0.514  +1.758   +0.03      0.98 +0.22
ZNS      +3.600     +2.030  +4.107   +0.51      0.76 -0.09
ZNSE     +2.700     +1.125  +3.094   +0.39      0.80 -0.00
ZNTE     +2.260     +0.927  +2.642   +0.38      0.78 -0.04
CDS      +2.420     +1.130  +2.847   +0.43      0.75 -0.08
CDSE     +1.675     +0.522  +2.164   +0.49      0.70 -0.16
CDTE     +1.440     +0.452  +1.973   +0.53      0.65 -0.23
HGS      -0.110     -0.403  +0.603   +0.71      0.29 -0.51
HGSE     -0.200     -0.964  -0.110   +0.09      0.89 +0.08
HGTE     -0.300     -0.937  +0.101   +0.40      0.61 -0.19
ALP      +2.450     +1.647  +2.712   +0.26      0.75 -0.05
ALAS     +2.160     +1.423  +2.359   +0.20      0.79 -0.01
ALSB     +1.620     +0.994  +1.590   -0.03      1.05 +0.15
GAP      +2.260     +1.604  +2.459   +0.20      0.77 -0.03
GAAS     +1.420     +0.386  +1.782   +0.36      0.74 -0.08
GASB     +0.730     -0.071  +1.076   +0.35      0.70 -0.12
INP      +1.340     +0.649  +1.620   +0.28      0.71 -0.09
INAS     +0.350     -0.360  +0.675   +0.33      0.69 -0.12
INSB     +0.170     -0.394  +0.537   +0.37      0.61 -0.18
PBS      +0.410     -0.037  +0.393   -0.02      1.04 +0.10
PBSE     +0.270     -0.180  +0.263   -0.01      1.02 +0.10
CDSE_WZ  +1.751     +0.554  +2.196   +0.45      0.73 -0.12
CDS_WZ   +2.420     +1.189  +2.920   +0.50      0.71 -0.15
======== ========== ======= ======== ========== ==== ===================

**Polytypes.** Zinc-blende and wurtzite CdSe and CdS converge to different bulk limits:
``CDSE``/``CDS`` (aliases ``CDSE_ZB``/``CDS_ZB``) are zinc blende, ``CDSE_WZ``/``CDS_WZ`` wurtzite. The
wurtzite PBE and PBE+SOC gaps are computed at the room-temperature lattice (CdSe a = 4.299, c = 7.01 Å;
CdS a = 4.137, c = 6.714 Å, COD 9011663) and at the PBE lattice of an isotropic energy-volume scan
(c/a and u of the CIF; CdS a = 4.215 Å). Wurtzite CdS has a QSGW+SOC gap (Deguchi16, 2.88 eV at
a = 4.160 Å, carried to the experimental lattice with the PBE deformation potential); wurtzite CdSe
has none, and takes :math:`\Delta_\Sigma` of zinc blende: in CdS the two polytypes give 1.717 and 1.731 eV,
so the self-energy correction is the same within 0.02 eV. ``MATERIAL_DB`` index 2 of the wurtzite
entries is the zinc-blende-equivalent lattice constant (same volume per formula unit), used for the
cluster sizes; the lattice of the dot is measured on its Cd–X bonds with the wurtzite relation
:math:`d = \sqrt{3/8}\,a`.

**Perovskites: structure of the reference.** The QSGW (Huang16) and PBE reference gaps are for the cubic
cell, while the lead halide perovskites are measured in their tilted room-temperature phases and the
dots relax with tilted octahedra (Pb–X–Pb about 150–170°). PBE and PBE+SOC of the measured structures
(``benchmarks/bulk_bands/inputs_rt_structures``, ``rt_structure_soc_gaps.json``):

* CsPbBr\ :sub:`3`: Pnma of Stoumpos et al., *Cryst. Growth Des.* 13, 2722 (2013), the volume of the
  cubic cell at :math:`a_\text{exp}`: 1.657 / 0.797 eV against 1.279 / 0.349 eV for the cubic cell; the
  tilts open the gap by 0.38 eV (0.45 eV with SOC).
* γ-CsPbI\ :sub:`3`: Pnma of Straus et al., *J. Am. Chem. Soc.* 141, 11435 (2019), 295 K (pseudo-cubic
  a = 6.19 Å): 1.402 / 0.514 eV against 0.832 / −0.187 eV (cubic PBE+SOC is inverted at R), +0.57 /
  +0.70 eV.
* CsPbCl\ :sub:`3`: no room-temperature structure is in COD; Pnma with the measured volume (pseudo-cubic
  5.60 Å) and the tilts relaxed in PBE (Pb–Cl–Pb 153–158°): 2.389 / 1.549 eV against 1.771 / 0.870 eV,
  +0.62 / +0.68 eV. This is an upper bound: at room temperature CsPbCl\ :sub:`3` is within about 20 K of
  its cubic transition and its tilts are smaller than the 0 K ones, so :math:`a_m` lies between 0.81
  (relaxed tilts) and 1.19 (cubic).

With the tilts, cubic QSGW+SOC plus the tilt opening is within 0.15 eV of the measured gaps of CsPbBr\ :sub:`3`
(2.15 against 2.30 eV) and CsPbI\ :sub:`3` (1.76 against 1.73 eV), instead of 0.6–0.7 eV below them against
the cubic gaps. The self-energy correction :math:`\Delta_\Sigma` is taken from the cubic cell; the tilts are
in the PBE gap of the dot itself.

**Confinement of perovskite dots.** Because the dots are tilted, their confinement energy
:math:`\Delta E_\text{conf} = E_g^\text{PBE}(\text{dot}) - E_g^\text{ref}` (Penn :math:`\epsilon_\text{eff}`
of the ΔW models and of the split vertex) is measured from the bulk PBE gap of the measured, tilted
structure at the dot's lattice,
:math:`E_g^\text{ref} = E_g^\text{PBE,cubic}(a_\text{dot}) + T_\text{exp} + s\,(T_\text{PBE} - T_\text{exp})`,
with :math:`T` the tilt opening at :math:`a_\text{exp}` (``gap_exp_structure`` − ``gap_exp_lattice``) and at
the PBE volume (``gap_pbe_structure``, the unfolded overlay data, − ``gap_pbe_lattice``): 0.62/0.48 eV
(Cl), 0.38/0.29 eV (Br), 0.57/0.45 eV (I). Against the cubic gap the tilt opening was counted as
confinement: the 4.1 nm CsPbBr\ :sub:`3` cube had :math:`\Delta E_\text{conf}` = 0.69 instead of 0.40 eV
and a vertex fraction f = 0.83 instead of 0.89. The self-energy shift keeps the cubic reference.

**Lattice of perovskite dots.** Tilts lengthen the Pb–X bond at fixed volume (Pnma CsPbBr\ :sub:`3`: Pb–Br
2.962 Å against 2.916 Å for the cubic cell of the same volume) while the Pb–Pb distance stays the
pseudo-cubic lattice constant (5.831 Å), so the strain fraction of perovskite dots is measured on Pb–Pb
distances. Frames of a 300 K trajectory can be expanded beyond the 0 K PBE lattice; that part is
thermal expansion, which the room-temperature experiment has too, so :math:`s` is capped at 1.

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
HgS_zb        F-43m       6.013     0.000        0.036            Γ-Γ        -0.534 / -0.524           -0.091     Hg d  
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
0.34, 0.42, 0.75, 0.95 eV). At its PBE lattice GaSb is inverted as well. HgS has the doublet
:math:`\Gamma_7` above the quartet (negative :math:`\Delta_{so}`), as in LDA (Svane et al. 2011).

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
