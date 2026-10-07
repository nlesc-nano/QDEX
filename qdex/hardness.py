import functools
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import ConvexHull, distance_matrix
from qdex.constants import (
    HA_TO_EV, ANG_PER_BOHR, BOHR_PER_ANG,
    IMAGE_CHARGE_CONST_EV_ANG, BRUS_KINETIC_EV_ANG2, VDW_SURFACE_ANG
)
import logging

logger = logging.getLogger(__name__)

# Atomic Hardness values (eta) in eV. 
# Used for the Ohno-Klopman damping in the Coulomb kernel.
HARDNESS_DICT = {
    'h': 6.4299, 'he': 12.5449, 'li': 2.3746, 'be': 3.4968, 'b': 4.619, 'c': 5.7410,
    'n': 6.8624, 'o': 7.9854, 'f': 9.1065, 'ne': 10.2303, 'na': 2.4441, 'mg': 3.0146,
    'al': 3.5849, 'si': 4.1551, 'p': 4.7258, 's': 5.2960, 'cl': 5.8662, 'ar': 6.4366,
    'k': 2.3273, 'ca': 2.7587, 'sc': 2.8582, 'ti': 2.9578, 'v': 3.0573, 'cr': 3.1567,
    'mn': 3.2564, 'fe': 3.3559, 'co': 3.4556, 'ni': 3.555, 'cu': 3.6544, 'zn': 3.7542,
    'ga': 4.1855, 'ge': 4.6166, 'as': 5.0662, 'se': 5.4795, 'br': 5.9111, 'kr': 6.3418,
    'rb': 2.1204, 'sr': 2.5374, 'y': 2.6335, 'zr': 2.7297, 'nb': 2.8260, 'mo': 2.9221,
    'tc': 3.0184, 'ru': 3.1146, 'rh': 3.2107, 'pd': 3.3069, 'ag': 3.4032, 'cd': 3.4994,
    'in': 3.9164, 'sn': 4.3332, 'sb': 4.7501, 'te': 5.167, 'i': 5.5839, 'xe': 6.0009,
    'cs': 0.6829, 'ba': 0.9201, 'la': 1.1571, 'ce': 1.3943, 'pr': 1.6315, 'nd': 1.8686,
    'pm': 2.1056, 'sm': 2.3427, 'eu': 2.5798, 'gd': 2.8170, 'tb': 3.0540, 'dy': 3.2912,
    'ho': 3.5283, 'er': 3.7655, 'tm': 4.0026, 'yb': 4.2395, 'lu': 4.4766, 'hf': 4.7065,
    'ta': 4.9508, 'w': 5.1879, 're': 5.4256, 'os': 5.6619, 'ir': 5.900, 'pt': 6.1367,
    'au': 6.3741, 'hg': 6.6103, 'tl': 1.7043, 'pb': 1.9435, 'bi': 2.1785, 'po': 2.4158,
    'at': 2.6528, 'rn': 2.8899, 'fr': 0.9882, 'ra': 1.2819, 'ac': 1.3497, 'th': 1.4175,
    'pa': 1.9368, 'u': 2.2305, 'np': 2.5241, 'pu': 3.0436, 'am': 3.4169, 'cm': 3.4050,
    'bk': 3.9244, 'cf': 4.2181, 'es': 4.5116, 'fm': 4.8051, 'md': 5.0100, 'no': 5.3926,
    'lr': 5.4607
}

# =====================================================================
# MATERIAL DATABASE (BULK QP REFERENCE: LITERATURE 1.0 QSGW + SOC, PBE AT THE SAME LATTICE)
# =====================================================================
# Updated Format: (eps_inf, Bohr_diam_A, lattice_A, gap_exp_eV, eps_static, m_eff, E_LO_eV, 
#                  gap_pbe_bulk_eV, gap_gw_bulk_eV, 
#                  R_mono_A, pbe_homo_mono, pbe_lumo_mono, gw_homo_mono, gw_lumo_mono)
#
# PHYSICAL RATIONALE: 1.0 QSGW vs. 0.8 QSGW IN COLLOIDAL QUANTUM DOTS (2-6 nm)
# ----------------------------------------------------------------------------
# In extended 3D bulk semiconductors, standard QSGW (in the random phase
# approximation, RPA) neglects attractive electron-hole ladder vertex corrections
# in the polarizability / screening W. In macroscopic solids, these vertex corrections
# increase macroscopic screening (raising eps_inf by ~10-20%), which compresses the
# quasiparticle self-energy Sigma = iGW and reduces the fundamental gap. Consequently,
# bulk electronic-structure practitioners often scale the QSGW self-energy by 0.8
# (e.g., "0.8 QSGW" / QSGW80) to mimic ladder diagrams and match low-temperature
# experimental bulk band gaps.
#
# However, in finite colloidal quantum dots (typically 2 - 6 nm diameter), quantum
# confinement fundamentally suppresses these vertex corrections:
# 1. Discrete Level Spectra: The continuous electron-hole continuum is replaced by
#    discrete particle-in-a-box states with large level spacings (Delta epsilon >> k_B T),
#    strongly quenching and detuning virtual electron-hole polarization pairs.
# 2. Screening Cloud Truncation: In bulk, the excitonic ladder correction requires a
#    long-range correlation volume spanning multiple unit cells (often R >= a_B). In
#    nanocrystals with R_QD <= a_B, the physical dielectric boundary cuts off this
#    spatial correlation cloud.
# 3. Dielectric Confinement: The large dielectric mismatch (eps_core ~ 6-10 vs.
#    eps_ligand ~ 2) strongly amplifies the bare classical image-charge self-energy,
#    which operates without bulk ladder reduction.
# Therefore, scaling by 0.8 inappropriately underestimates the quasiparticle opening
# in 2-6 nm QDs (e.g. predicting ~2.2 eV for a 3 nm CdSe dot instead of the experimental
# ~2.6 eV). Pure 1.0 QSGW provides the correct high-energy, vertex-quenched asymptote
# matching molecular/vacuum evGW anchors and accurately reproducing experimental sizing
# curves (such as Hens & Rodina, Nano Lett. 2022).
#
# BULK QP SHIFT: INDEX 7 (PBE), INDEX 8 (QSGW) AND SPIN-ORBIT COUPLING
# -------------------------------------------------------------------
# The bulk QP shift applied to the PBE orbitals of a dot is
#     Delta_bulk = Delta_Sigma + Delta_geom,
#     Delta_Sigma = E_g^QSGW+SOC(a_exp) - E_g^PBE+SOC(a_exp)      = index 8 - index 7,
#     Delta_geom  = E_g^PBE(a_exp) - E_g^PBE(a_dot)                (bulk_geometry_shift, below).
# Delta_Sigma is the bulk self-energy correction: QSGW and PBE at the same lattice (a_exp, index 2)
# and with the same spin-orbit treatment. The QSGW gaps are literature values WITH SOC
# (QSGW_SOC_LITERATURE, after MATERIAL_DB_BULK_PBE; SOC added to the converged QSGW Hamiltonian),
# carried to a_exp where the source used another lattice; the PBE+SOC gaps are ours at the level of
# theory of the dots (CP2K PBE, DZVP-MOLOPT-PBE-GTH, GTH-PBE of POTENTIAL_UZH, the QDEX GTH-SOC
# operator). Gaps change with volume by several eV per unit strain (GaAs PBE: 0.50 eV at a_exp,
# 0.02 eV at its 2 % larger PBE lattice), so both must be at the same lattice.
#
# Index 7 is the spin-free PBE gap at a_exp; index 8 is the spin-free QSGW-equivalent gap
#     index 8 = index 7 + Delta_Sigma = E_g^QSGW+SOC(a_exp) + [E_g^PBE(a_exp) - E_g^PBE+SOC(a_exp)],
# so the same Delta_Sigma corrects the spin-free and the SOC dots, and the SOC lowering of the gap is
# the PBE/GTH-SOC one in both the bulk reference and the dots (the QDEX spinors). Gaps are band-edge
# transitions, signed (negative = inverted): E(Gamma1/Gamma6) - E(Gamma15/Gamma8) for direct zinc blende
# (inverted in PBE for HgS, HgSe, HgTe, InAs, InSb; with SOC also GaSb), the fundamental gap for the
# indirect zinc blendes (X), L in rock salt and R in the cubic perovskites (PBS, PBSE and CSPBI3 are
# inverted in PBE+SOC at a_exp). Earlier versions converted literature gaps to spin-free ones with
# Delta_so/3 (zinc blende) or 2/3 Delta_so(CB) (perovskites); the III-V values were in fact the
# spin-free QSGW of Deguchi et al., and the II-VI, Hg and perovskite values matched no QSGW source
# (the perovskite ones were 0.5-0.6 eV above Huang and Lambrecht's QSGW).
#
# A dot relaxed with PBE is (almost) at the PBE lattice, a_PBE > a_exp, so its PBE gap carries the
# bulk PBE gap change between a_exp and its own lattice; Delta_geom (quasiparticles.bulk_geometry)
# removes it, and the QP gap is that of the dot at the experimental lattice. This needs no QSGW at
# a_PBE: Delta_Sigma(a_dot) + [E^QSGW(a_exp) - E^QSGW(a_dot)] = Delta_Sigma(a_exp) + Delta_geom exactly
# (which matters: the self-energy correction itself depends on the lattice, strongly so in the lead
# halide perovskites, where QSGW's dEg/dlnV is about twice PBE's). The vertex correction
# (quasiparticles.bulk_vertex) scales Delta_Sigma only.
#
# MATERIAL_DB_BULK_PBE (below) keeps the spin-free and PBE+SOC gaps at a_exp and at a_PBE, the QSGW+SOC
# gap at a_exp and Delta_Sigma; the bulk bands drawn over the fuzzy bands (qdex/data/bulk_bands) are
# at the PBE lattice, like the dots. Computed with qdex.bulk_soc; inputs, runs and the construction of
# the table in benchmarks/bulk_bands; see docs/electronic_structure/bulk_bands.rst.
#
# Monomer values are computed with CP2K (DZVP-RI/GTH/PBE) then evGW.

MATERIAL_DB = {
    "CS3BI2BR9": (3.9, 24.0, 8.01, 2.55, 11.7, 0.57, 0.021, 3.33, 4.35, 7.850, -5.5000, -1.8500, -7.2000, 0.1500),
    # Halide Perovskites (cubic; index 8 from Huang & Lambrecht QSGW+SO, see QSGW_SOC_LITERATURE)
    "CSPBCL3": (4.0, 25.0, 5.6, 3.0, 20.0, 0.15, 0.022, 1.771, 3.558, 7.604, -4.7512, -1.5076, -6.5706, 0.3314),
    "CSPBBR3": (4.8, 35.0, 5.83, 2.3, 25.0, 0.12, 0.018, 1.278, 2.628, 7.9149, -4.7779, -1.6961, -6.3775, 0.4555),
    "CSPBI3": (5.1, 60.0, 6.2, 1.73, 30.0, 0.1, 0.015, 0.832, 2.076, 8.3726, -4.4874, -1.7988, -5.6323, 0.3893),
    "MAPBI3": (6.5, 45.0, 6.27, 1.55, 30.0, 0.1, 0.015, 1.55, 2.73),
    "FAPBI3": (6.2, 50.0, 6.36, 1.48, 28.0, 0.1, 0.015, 1.45, 2.60),

    # II–VI Semiconductors (zinc blende; index 8 from Deguchi et al. / Svane et al. QSGW+SO)
    "ZNS": (5.1, 25.0, 5.41, 3.6, 8.9, 0.28, 0.043, 2.051, 4.128, 4.8204, -6.6591, -3.208, -8.2883, -1.0303),
    "ZNSE": (5.9, 38.0, 5.67, 2.7, 8.6, 0.17, 0.031, 1.252, 3.221, 4.9394, -6.4914, -3.1993, -7.9032, -1.063),
    "ZNTE": (6.7, 52.0, 6.1, 2.26, 9.8, 0.12, 0.026, 1.225, 2.940, 5.1907, -5.789, -2.8062, -7.0003, -0.901),
    "CDS": (5.4, 30.0, 5.82, 2.42, 8.9, 0.16, 0.037, 1.146, 2.863, 5.1533, -6.598, -3.88, -8.2128, -1.8481),
    "CDSE": (6.2, 56.0, 6.05, 1.675, 9.5, 0.13, 0.026, 0.644, 2.286, 5.31331, -6.4196, -3.7852, -7.8177, -1.7884),
    "CDSE_WZ": (6.2, 56.0, 6.077, 1.751, 9.5, 0.13, 0.026, 0.675, 2.317, 5.31331, -6.4196, -3.7852, -7.8177, -1.7884),   # wurtzite; index 2: zinc-blende-equivalent a (volume per formula unit)
    "CDS_WZ": (5.4, 30.0, 5.839, 2.42, 8.9, 0.16, 0.037, 1.205, 2.936, 5.1533, -6.598, -3.88, -8.2128, -1.8481),
    "CDTE": (7.1, 73.0, 6.48, 1.44, 10.2, 0.1, 0.021, 0.740, 2.261, 5.4433, -5.7978, -2.9767, -6.9951, -0.9143),
    "HGS": (11.3, 50.0, 5.85, 0.5, 13.0, 0.2, 0.03, -0.421, 0.585, 5.0759, -6.8468, -4.6185, -8.0873, -2.6135),
    "HGSE": (14.0, 460.0, 6.08, 0.0, 18.0, 0.04, 0.017, -0.872, -0.018, 5.2017, -6.5897, -4.4284, -7.6523, -2.4936),
    "HGTE": (15.0, 400.0, 6.46, 0.0, 20.0, 0.03, 0.015, -0.665, 0.373, 5.4091, -5.9729, -3.7913, -6.9007, -1.9993),

    # III–V Semiconductors (zinc blende; index 8 from Deguchi et al. QSGW+SO)
    "ALP": (7.5, 15.0, 5.46, 2.45, 10.0, 0.2, 0.05, 1.667, 2.732, 5.1769, -6.4313, -4.4894, -8.2744, -2.4365),
    "ALAS": (8.2, 30.0, 5.66, 2.16, 10.1, 0.15, 0.049, 1.522, 2.458, 5.2267, -6.2214, -4.1183, -7.8386, -2.1083),
    "ALSB": (10.2, 60.0, 6.14, 1.62, 12.0, 0.14, 0.036, 1.208, 1.804, 5.5593, -5.9163, -3.9868, -7.301, -2.1319),
    "GAP": (9.1, 15.0, 5.45, 2.26, 11.1, 0.15, 0.049, 1.632, 2.487, 5.3208, -6.562, -4.5789, -8.2171, -2.8192),
    "GAAS": (10.9, 100.0, 5.65, 1.42, 13.1, 0.067, 0.036, 0.496, 1.891, 5.4503, -6.3334, -4.4405, -7.7475, -2.7198),
    "GASB": (14.4, 200.0, 6.1, 0.73, 15.7, 0.04, 0.028, 0.156, 1.304, 5.7221, -6.0195, -4.5359, -7.2596, -2.9446),
    "INP": (9.6, 150.0, 5.87, 1.34, 12.4, 0.08, 0.042, 0.680, 1.651, 5.6309, -6.5944, -4.8483, -8.0539, -3.0678),
    "INAS": (11.8, 340.0, 6.06, 0.35, 15.0, 0.023, 0.029, -0.247, 0.788, 5.697, -6.4406, -4.8382, -7.7358, -3.1218),
    "INSB": (15.7, 650.0, 6.48, 0.17, 17.9, 0.014, 0.023, -0.153, 0.778, 5.9465, -6.184, -4.6223, -7.3335, -3.0342),

    # IV–VI Semiconductors
    "PBS": (17.2, 200.0, 5.94, 0.41, 23.0, 0.09, 0.027, 0.253, 0.683, 7.2378, -6.638, -4.3253, -8.128, -2.3724),
    "PBSE": (22.9, 460.0, 6.12, 0.27, 30.0, 0.07, 0.017, 0.139, 0.582, 7.3935, -6.5104, -4.2415, -7.8234, -2.3177),

    "DEFAULT": (1.0, 1.0, 5.0, 0.0, 1.0, 1.0, 0.02, 0.00, 0.00)

}

# Bulk gaps (eV) at the level of theory of the dots, at the experimental lattice constant a_exp
# (= MATERIAL_DB index 2) and at the PBE equilibrium lattice a_pbe (energy-volume minimum with the same
# basis, pseudopotentials and cutoff; the bulk bands of the fuzzy-band overlay are at a_pbe):
# spin-free PBE (gap_exp_lattice = MATERIAL_DB index 7, gap_pbe_lattice), PBE + GTH-SOC
# (gap_soc_exp_lattice, gap_soc_pbe_lattice), the literature QSGW+SOC gap carried to a_exp
# (qsgw_soc_exp_lattice) and Delta_Sigma = qsgw_soc_exp_lattice - gap_soc_exp_lattice
# (= MATERIAL_DB index 8 - index 7). gap_exp_structure / gap_soc_exp_structure: PBE and PBE+SOC gaps of
# the measured room-temperature structure when it is not the reference one (CsPbBr3: orthorhombic Pnma
# of Stoumpos et al., Cryst. Growth Des. 13, 2722 (2013), the same volume as the cubic cell at a_exp;
# the octahedral tilts open the gap by 0.38 eV, 0.45 eV with SOC; gamma-CsPbI3: Pnma of Straus et al.,
# J. Am. Chem. Soc. 141, 11435 (2019), 295 K, +0.57 / +0.70 eV; CsPbCl3: Pnma at the measured volume
# with the tilts relaxed in PBE (no room-temperature structure in COD), +0.62 / +0.68 eV, an upper
# bound: CsPbCl3 is within ~20 K of its cubic transition at room temperature and the 0 K tilts are
# larger). Wurtzite: a_exp/c_exp (hexagonal) of the room-temperature structure, isotropic PBE scan.
# Signed band-edge transitions (negative = inverted; see the header
# of MATERIAL_DB). Lattice constants in A (conventional cell; wurtzite: hexagonal a, isotropic scan
# with c/a and u of the CIF). B0: bulk modulus of the energy-volume fit (GPa). PBTE and CDSE_WZ: overlay
# bands only.
MATERIAL_DB_BULK_PBE = {
    "CSPBCL3": dict(a_exp=5.600, gap_exp_lattice=1.771, gap_soc_exp_lattice=0.870, a_pbe=5.759, gap_pbe_lattice=2.113, gap_soc_pbe_lattice=1.211, qsgw_soc_exp_lattice=2.657, delta_sigma=1.787, B0=21.6,
                    gap_exp_structure=2.389, gap_soc_exp_structure=1.549),   # Pnma, PBE tilts: upper bound
    "CSPBBR3": dict(a_exp=5.830, gap_exp_lattice=1.278, gap_soc_exp_lattice=0.349, a_pbe=6.024, gap_pbe_lattice=1.671, gap_soc_pbe_lattice=0.738, qsgw_soc_exp_lattice=1.699, delta_sigma=1.350, B0=18.3,
                    gap_exp_structure=1.657, gap_soc_exp_structure=0.797),   # Pnma (RT) structure, see below
    "CSPBI3": dict(a_exp=6.200, gap_exp_lattice=0.832, gap_soc_exp_lattice=-0.187, a_pbe=6.417, gap_pbe_lattice=1.207, gap_soc_pbe_lattice=0.155, qsgw_soc_exp_lattice=1.057, delta_sigma=1.244, B0=15.2,
                   gap_exp_structure=1.402, gap_soc_exp_structure=0.514),   # gamma-CsPbI3 Pnma (295 K)
    "ZNS": dict(a_exp=5.410, gap_exp_lattice=2.051, gap_soc_exp_lattice=2.030, a_pbe=5.424, gap_pbe_lattice=2.018, gap_soc_pbe_lattice=1.997, qsgw_soc_exp_lattice=4.107, delta_sigma=2.077, B0=75.2),
    "ZNSE": dict(a_exp=5.670, gap_exp_lattice=1.252, gap_soc_exp_lattice=1.125, a_pbe=5.752, gap_pbe_lattice=1.085, gap_soc_pbe_lattice=0.958, qsgw_soc_exp_lattice=3.094, delta_sigma=1.969, B0=56.7),
    "ZNTE": dict(a_exp=6.100, gap_exp_lattice=1.225, gap_soc_exp_lattice=0.927, a_pbe=6.206, gap_pbe_lattice=0.984, gap_soc_pbe_lattice=0.689, qsgw_soc_exp_lattice=2.642, delta_sigma=1.715, B0=42.4),
    "CDS": dict(a_exp=5.820, gap_exp_lattice=1.146, gap_soc_exp_lattice=1.130, a_pbe=5.947, gap_pbe_lattice=1.004, gap_soc_pbe_lattice=0.987, qsgw_soc_exp_lattice=2.847, delta_sigma=1.717, B0=53.7),
    "CDSE": dict(a_exp=6.050, gap_exp_lattice=0.644, gap_soc_exp_lattice=0.522, a_pbe=6.217, gap_pbe_lattice=0.474, gap_soc_pbe_lattice=0.353, qsgw_soc_exp_lattice=2.164, delta_sigma=1.642, B0=45.3),
    "CDTE": dict(a_exp=6.480, gap_exp_lattice=0.740, gap_soc_exp_lattice=0.452, a_pbe=6.630, gap_pbe_lattice=0.530, gap_soc_pbe_lattice=0.245, qsgw_soc_exp_lattice=1.973, delta_sigma=1.521, B0=37.0),
    "HGS": dict(a_exp=5.850, gap_exp_lattice=-0.421, gap_soc_exp_lattice=-0.403, a_pbe=6.013, gap_pbe_lattice=-0.534, gap_soc_pbe_lattice=-0.524, qsgw_soc_exp_lattice=0.603, delta_sigma=1.006, B0=50.0),
    "HGSE": dict(a_exp=6.080, gap_exp_lattice=-0.872, gap_soc_exp_lattice=-0.964, a_pbe=6.285, gap_pbe_lattice=-1.004, gap_soc_pbe_lattice=-1.102, qsgw_soc_exp_lattice=-0.110, delta_sigma=0.854, B0=42.8),
    "HGTE": dict(a_exp=6.460, gap_exp_lattice=-0.665, gap_soc_exp_lattice=-0.937, a_pbe=6.679, gap_pbe_lattice=-0.908, gap_soc_pbe_lattice=-1.179, qsgw_soc_exp_lattice=0.101, delta_sigma=1.038, B0=35.2),
    "ALP": dict(a_exp=5.460, gap_exp_lattice=1.667, gap_soc_exp_lattice=1.647, a_pbe=5.511, gap_pbe_lattice=1.724, gap_soc_pbe_lattice=1.704, qsgw_soc_exp_lattice=2.712, delta_sigma=1.065, B0=81.3),
    "ALAS": dict(a_exp=5.660, gap_exp_lattice=1.522, gap_soc_exp_lattice=1.423, a_pbe=5.738, gap_pbe_lattice=1.592, gap_soc_pbe_lattice=1.495, qsgw_soc_exp_lattice=2.359, delta_sigma=0.936, B0=66.6),
    "ALSB": dict(a_exp=6.140, gap_exp_lattice=1.208, gap_soc_exp_lattice=0.994, a_pbe=6.243, gap_pbe_lattice=1.212, gap_soc_pbe_lattice=1.004, qsgw_soc_exp_lattice=1.590, delta_sigma=0.596, B0=48.6),
    "GAP": dict(a_exp=5.450, gap_exp_lattice=1.632, gap_soc_exp_lattice=1.604, a_pbe=5.433, gap_pbe_lattice=1.614, gap_soc_pbe_lattice=1.586, qsgw_soc_exp_lattice=2.459, delta_sigma=0.855, B0=85.9),
    "GAAS": dict(a_exp=5.650, gap_exp_lattice=0.496, gap_soc_exp_lattice=0.386, a_pbe=5.774, gap_pbe_lattice=0.022, gap_soc_pbe_lattice=-0.086, qsgw_soc_exp_lattice=1.782, delta_sigma=1.395, B0=55.0),
    "GASB": dict(a_exp=6.100, gap_exp_lattice=0.156, gap_soc_exp_lattice=-0.071, a_pbe=6.270, gap_pbe_lattice=-0.425, gap_soc_pbe_lattice=-0.646, qsgw_soc_exp_lattice=1.076, delta_sigma=1.148, B0=39.7),
    "INP": dict(a_exp=5.870, gap_exp_lattice=0.680, gap_soc_exp_lattice=0.649, a_pbe=5.978, gap_pbe_lattice=0.377, gap_soc_pbe_lattice=0.346, qsgw_soc_exp_lattice=1.620, delta_sigma=0.971, B0=59.1),
    "INAS": dict(a_exp=6.060, gap_exp_lattice=-0.247, gap_soc_exp_lattice=-0.360, a_pbe=6.204, gap_pbe_lattice=-0.603, gap_soc_pbe_lattice=-0.713, qsgw_soc_exp_lattice=0.675, delta_sigma=1.035, B0=48.7),
    "INSB": dict(a_exp=6.480, gap_exp_lattice=-0.153, gap_soc_exp_lattice=-0.394, a_pbe=6.661, gap_pbe_lattice=-0.622, gap_soc_pbe_lattice=-0.855, qsgw_soc_exp_lattice=0.537, delta_sigma=0.931, B0=36.6),
    "PBS": dict(a_exp=5.940, gap_exp_lattice=0.253, gap_soc_exp_lattice=-0.037, a_pbe=6.040, gap_pbe_lattice=0.421, gap_soc_pbe_lattice=0.134, qsgw_soc_exp_lattice=0.393, delta_sigma=0.430, B0=52.7),
    "PBSE": dict(a_exp=6.120, gap_exp_lattice=0.139, gap_soc_exp_lattice=-0.180, a_pbe=6.241, gap_pbe_lattice=0.325, gap_soc_pbe_lattice=0.006, qsgw_soc_exp_lattice=0.263, delta_sigma=0.443, B0=46.6),
    "PBTE": dict(a_exp=None, gap_exp_lattice=None, a_pbe=6.584, gap_pbe_lattice=0.616, B0=39.8),
    "CDSE_WZ": dict(a_exp=4.299, c_exp=7.013, gap_exp_lattice=0.675, gap_soc_exp_lattice=0.554, a_pbe=4.399,
                    c_pbe=7.177, gap_pbe_lattice=0.522, gap_soc_pbe_lattice=0.401, qsgw_soc_exp_lattice=2.196,
                    delta_sigma=1.642, B0=45.1),   # Delta_Sigma of zinc blende (no QSGW+SOC for wurtzite CdSe)
    "CDS_WZ": dict(a_exp=4.137, c_exp=6.714, gap_exp_lattice=1.205, gap_soc_exp_lattice=1.189, a_pbe=4.215,
                   c_pbe=6.841, gap_pbe_lattice=1.072, gap_soc_pbe_lattice=1.054, qsgw_soc_exp_lattice=2.920,
                   delta_sigma=1.731, B0=53.5),
}

# Literature QSGW gaps WITH spin-orbit coupling (100 % QSGW, RPA W; SOC added to the converged QSGW
# Hamiltonian), the reference of the bulk self-energy correction (see "BULK QP SHIFT" above MATERIAL_DB):
#   material: (gap_eV, a_lit_A, transition, source, dEg/dlnV_eV or None)
# The gap is the band-edge transition, signed (negative = inverted): E(Gamma6) - E(Gamma8) for the
# direct zinc blendes and the Hg compounds (E0), the fundamental gap of the indirect zinc blendes (X), L
# in rock salt, R in the cubic perovskites. dEg/dlnV, when the source gives it, carries the gap from
# a_lit to a_exp (MATERIAL_DB index 2); otherwise the PBE deformation potential does.
#   Deguchi16: D. Deguchi, K. Sato, H. Kino, T. Kotani, Jpn. J. Appl. Phys. 55, 051201 (2016),
#              Table II (QSGW+SO), Table III (Gamma6c of GaSb); lattice constants of its Table I.
#   Svane11:   A. Svane et al., Phys. Rev. B 84, 205205 (2011), Table I (QSGW E0).
#   Huang16:   L.-y. Huang and W. R. L. Lambrecht, Phys. Rev. B 93, 195211 (2016), Tables I-II.
#   Svane10:   A. Svane et al., Phys. Rev. B 81, 245120 (2010), Table I (QSGW) at the low-temperature
#              lattices of its Ref. 43, Table V (QSGW deformation potentials). Deguchi16 gives 0.49 eV for
#              PbS at 5.936 A, about 0.1 eV above Svane10 carried to the same lattice.
QSGW_SOC_LITERATURE = {
    "ALP": (2.72, 5.467, "Gamma-X (indirect)", "Deguchi16", None),
    "ALAS": (2.36, 5.661, "Gamma-X (indirect)", "Deguchi16", None),
    "ALSB": (1.59, 6.136, "Gamma-X (indirect)", "Deguchi16", None),
    "GAP": (2.46, 5.451, "Gamma-X (indirect)", "Deguchi16", None),
    "GAAS": (1.77, 5.653, "Gamma6-Gamma8", "Deguchi16", None),
    "GASB": (1.09, 6.096, "Gamma6-Gamma8 (Table III; the QSGW minimum is at L)", "Deguchi16", None),
    "INP": (1.62, 5.870, "Gamma6-Gamma8", "Deguchi16", None),
    "INAS": (0.68, 6.058, "Gamma6-Gamma8", "Deguchi16", None),
    "INSB": (0.54, 6.479, "Gamma6-Gamma8", "Deguchi16", None),
    "ZNS": (4.10, 5.413, "Gamma6-Gamma8", "Deguchi16", None),
    "ZNSE": (3.10, 5.667, "Gamma6-Gamma8", "Deguchi16", None),
    "ZNTE": (2.64, 6.101, "Gamma6-Gamma8", "Deguchi16", None),
    "CDS": (2.84, 5.826, "Gamma6-Gamma8", "Deguchi16", None),
    "CDS_WZ": (2.88, 4.160, "Gamma7c-Gamma9v (wurtzite; c = 6.756 A)", "Deguchi16", None),
    "CDSE": (2.16, 6.054, "Gamma6-Gamma8", "Deguchi16", None),
    "CDTE": (1.97, 6.482, "Gamma6-Gamma8", "Deguchi16", None),
    "HGS": (0.61, 5.84, "E0 = Gamma6-Gamma8", "Svane11", None),
    "HGSE": (-0.11, 6.08, "E0 = Gamma6-Gamma8", "Svane11", None),
    "HGTE": (0.09, 6.47, "E0 = Gamma6-Gamma8", "Svane11", None),
    "PBS": (0.31, 5.909, "L6+-L6-", "Svane10", 5.3),
    "PBSE": (0.21, 6.098, "L6+-L6-", "Svane10", 4.9),
    "CSPBCL3": (2.678, 5.605, "R", "Huang16", 7.7),
    "CSPBBR3": (1.868, 5.874, "R", "Huang16", 7.5),
    "CSPBI3": (1.331, 6.289, "R", "Huang16", 6.4),
}

# Core inorganic elements for each material (ignores organic ligands like MA/FA)
MATERIAL_ELEMENTS = {
    "CSPBCL3": ["Cs", "Pb", "Cl"], "CSPBBR3": ["Cs", "Pb", "Br"], "CSPBI3":  ["Cs", "Pb", "I"],
    "MAPBI3":  ["Pb", "I"], "FAPBI3":  ["Pb", "I"], 
    "ZNS":     ["Zn", "S"], "ZNSE":    ["Zn", "Se"], "ZNTE":    ["Zn", "Te"],
    "CDS":     ["Cd", "S"], "CDSE":    ["Cd", "Se"], "CDTE":    ["Cd", "Te"],
    "HGS":     ["Hg", "S"], "HGSE":    ["Hg", "Se"], "HGTE":    ["Hg", "Te"],
    "ALP":     ["Al", "P"], "ALAS":    ["Al", "As"], "ALSB":    ["Al", "Sb"],
    "GAP":     ["Ga", "P"], "GAAS":    ["Ga", "As"], "GASB":    ["Ga", "Sb"],
    "INP":     ["In", "P"], "INAS":    ["In", "As"], "INSB":    ["In", "Sb"],
    "PBS":     ["Pb", "S"], "PBSE":    ["Pb", "Se"]
}

# Optical Refractive Indices (at optical / band-edge frequencies)
# Extracted from experimental literature and bulk eps_inf (n_r ≈ sqrt(eps_inf))
REFRACTIVE_INDEX_DICT = {
    "CS3BI2BR9": 1.97,  # sqrt(3.9)
    "CSPBCL3": 2.00,    # sqrt(4.0)
    "CSPBBR3": 2.19,    # sqrt(4.8)
    "CSPBI3": 2.26,     # sqrt(5.1)
    "MAPBI3": 2.55,     # sqrt(6.5)
    "FAPBI3": 2.49,     # sqrt(6.2)
    "ZNS": 2.26,        # sqrt(5.1)
    "ZNSE": 2.43,       # sqrt(5.9)
    "ZNTE": 2.59,       # sqrt(6.7)
    "CDS": 2.32,        # sqrt(5.4)
    "CDSE": 2.49,       # sqrt(6.2)
    "CDTE": 2.66,       # sqrt(7.1)
    "HGS": 3.36,        # sqrt(11.3)
    "HGSE": 3.74,       # sqrt(14.0)
    "HGTE": 3.87,       # sqrt(15.0)
    "ALP": 2.74,        # sqrt(7.5)
    "ALAS": 2.86,       # sqrt(8.2)
    "ALSB": 3.19,       # sqrt(10.2)
    "GAP": 3.02,        # sqrt(9.1)
    "GAAS": 3.30,       # sqrt(10.9)
    "GASB": 3.79,       # sqrt(14.4)
    "INP": 3.10,        # sqrt(9.6)
    "INAS": 3.44,       # sqrt(11.8)
    "INSB": 3.96,       # sqrt(15.7)
    "PBS": 4.15,        # sqrt(17.2)
    "PBSE": 4.79,       # sqrt(22.9)
    "DEFAULT": 2.0
}

def get_refractive_index(material_name):
    """Returns the optical refractive index n_r for a given semiconductor material."""
    mat_key = str(material_name).upper().strip() if material_name else "DEFAULT"
    if mat_key in REFRACTIVE_INDEX_DICT:
        return REFRACTIVE_INDEX_DICT[mat_key]
    if mat_key in MATERIAL_DB:
        return float(np.sqrt(MATERIAL_DB[mat_key][0]))
    return REFRACTIVE_INDEX_DICT["DEFAULT"]


def compute_radiative_rates(E_ev, f_osc, refractive_index=2.0):
    """
    Computes Einstein A spontaneous emission rates (s^-1 and fs^-1) from exciton energies (eV)
    and oscillator strengths f_osc in a dielectric medium of refractive index n_r::

      k_rad = (2 * n_r * e^4 * E^2 * f) / (4 * pi * eps_0 * m_e * c^3 * hbar^2)
            = C_rad * n_r * E^2 * f

    where C_rad = 4.3391988e7 s^-1 eV^-2 (4.3391988e-8 fs^-1 eV^-2).
    """
    C_RAD_S = 4.3391988e7  # s^-1 * eV^-2
    E_arr = np.asarray(E_ev, dtype=np.float64)
    f_arr = np.asarray(f_osc, dtype=np.float64)
    k_rad_s = C_RAD_S * refractive_index * (np.maximum(E_arr, 0.0) ** 2) * np.maximum(f_arr, 0.0)
    k_rad_fs = k_rad_s * 1e-15  # fs^-1
    return k_rad_s, k_rad_fs


def compute_energy_gap_law_rate(E_gap_ev, E_LO_ev=0.018, S_hr=1.0, A_nr=1e13):
    """
    Computes multi-phonon non-radiative recombination rate via the Englman-Jortner Energy Gap Law::

      k_nr = A_nr * exp(-gamma * (E_gap / E_LO))

    where gamma = ln(E_gap / (S * E_LO)) - 1.
    """
    if E_gap_ev <= 0.0 or E_LO_ev <= 0.0:
        return 0.0, 0.0
    num_phonons = E_gap_ev / E_LO_ev
    gamma = max(0.1, np.log(max(num_phonons / max(S_hr, 1e-3), 1.01)) - 1.0)
    exponent = -gamma * num_phonons
    exponent = max(-100.0, min(0.0, exponent))
    k_nr_s = A_nr * np.exp(exponent)
    k_nr_fs = k_nr_s * 1e-15
    return k_nr_s, k_nr_fs


def compute_fcwd_rate(E_gap_ev, V_el_ev, lambda_ev, sigma_ev):
    r"""
    Computes non-radiative recombination rate via Fermi's Golden Rule with
    Franck-Condon Weighted Density of States (FCWD)::

      k_nr = (2 * pi / hbar) * |V_el|^2 * FCWD(E_gap)
      FCWD(E_gap) = (1 / sqrt(2 * pi * sigma^2)) * exp(-(E_gap - lambda)^2 / (2 * sigma^2))

    where:
      - V_el is the effective electronic coupling (eV), e.g. :math:`\hbar \langle |d_{10}| \rangle` from NAMD
      - lambda is the nuclear reorganization energy (eV)
      - sigma is the thermal Gaussian broadening (eV), sigma = sqrt(2 * lambda * k_B * T) = sigma_E
      - E_gap is the transition energy (eV)

    Returns:
      (k_nr_s, k_nr_fs) in s^-1 and fs^-1.
    """
    HBAR_EV_FS = 0.6582119569  # eV * fs
    sigma = max(float(sigma_ev), 1e-6)
    V_el = float(V_el_ev)
    E_g = float(E_gap_ev)
    lam = float(lambda_ev)

    fcwd = (1.0 / (np.sqrt(2.0 * np.pi) * sigma)) * np.exp(-((E_g - lam) ** 2) / (2.0 * (sigma ** 2)))
    k_nr_fs = (2.0 * np.pi / HBAR_EV_FS) * (V_el ** 2) * fcwd
    k_nr_s = k_nr_fs * 1e15
    return k_nr_s, k_nr_fs


def extract_recombination_parameters_from_namd(
    var_E_gap_ev2,
    dominant_freq_cm1=None,
    temp_k=300.0,
    mean_nac_fs=None,
    material_name=None
):
    r"""
    Extracts non-empirical physical parameters for radiative and non-radiative
    recombination from NAMD trajectory data:

    1. Dominant optical phonon energy :math:`\hbar \omega_{\mathrm{LO}}` from spectral density peak:
       E_LO = h * c * nu_peak  (eV)
    2. Nuclear reorganization energy :math:`\lambda` from classical fluctuation-dissipation:
       :math:`\lambda = \sigma_E^2 / (2 k_B T)`  (eV)
    3. Dimensionless Huang-Rhys factor S:
       :math:`S = \lambda / (\hbar \omega_{\mathrm{LO}}) = \sigma_E^2 / (2 k_B T \hbar \omega_{\mathrm{LO}})`
    4. Gaussian broadening parameter :math:`\sigma`:
       :math:`\sigma = \sqrt{\sigma_E^2}`  (eV)
    5. Effective electronic coupling :math:`V_{\mathrm{el}}` from mean non-adiabatic coupling:
       :math:`V_{\mathrm{el}} = \hbar \langle |d_{10}| \rangle`  (eV)

    Parameters:
      var_E_gap_ev2: float, variance of the band edge gap fluctuations :math:`\sigma_E^2` in eV^2.
      dominant_freq_cm1: float, dominant vibrational mode wavenumber in cm^-1 from J(omega).
      temp_k: float, temperature in Kelvin.
      mean_nac_fs: float, mean non-adiabatic coupling magnitude in fs^-1.
      material_name: str, material name fallback if dominant_freq_cm1 is not available.

    Returns:
      dict with keys: 'E_LO_ev', 'dominant_freq_cm1', 'lambda_ev', 'S_hr', 'sigma_ev', 'V_el_ev'
    """
    KB_EV = 8.617333262e-5  # eV / K
    HBAR_EV_FS = 0.6582119569  # eV * fs
    HC_EV_CM = 1.239841984e-4  # eV * cm

    # 1. Optical phonon energy \hbar\omega_LO
    if dominant_freq_cm1 is not None and dominant_freq_cm1 > 0:
        E_LO_ev = float(dominant_freq_cm1 * HC_EV_CM)
        freq_cm1 = float(dominant_freq_cm1)
    elif material_name and str(material_name).upper() in MATERIAL_DB:
        mat_entry = MATERIAL_DB[str(material_name).upper()]
        E_LO_ev = float(mat_entry[6]) if len(mat_entry) > 6 else 0.018
        freq_cm1 = E_LO_ev / HC_EV_CM
    else:
        E_LO_ev = 0.018  # default ~145 cm^-1 (typical perovskite Pb-X LO mode)
        freq_cm1 = E_LO_ev / HC_EV_CM

    # 2. Nuclear Reorganization Energy \lambda = \sigma_E^2 / (2 * k_B * T)
    var_g = max(float(var_E_gap_ev2), 1e-12)
    T = max(float(temp_k), 1.0)
    lambda_ev = var_g / (2.0 * KB_EV * T)

    # 3. Huang-Rhys factor S = \lambda / (\hbar\omega_LO)
    S_hr = lambda_ev / max(E_LO_ev, 1e-6)

    # 4. Thermal Gaussian broadening \sigma = \sqrt{\sigma_E^2}
    sigma_ev = np.sqrt(var_g)

    # 5. Effective electronic coupling V_el = \hbar * \langle |d_10| \rangle
    V_el_ev = float(HBAR_EV_FS * mean_nac_fs) if mean_nac_fs is not None else None

    return {
        "E_LO_ev": E_LO_ev,
        "dominant_freq_cm1": freq_cm1,
        "lambda_ev": lambda_ev,
        "S_hr": S_hr,
        "sigma_ev": sigma_ev,
        "V_el_ev": V_el_ev
    }


def estimate_gw_qp_gap(
    coords, atom_symbols, material_name, eps_out, return_details=False,
    regularization_length_ang=1.0, residual_power=2.0, strict=False,
    polarization_model="sphere", dft_gap=None, residual_scaling="econf",
    anchor_residual="on",
):
    """
    Estimate the PBE-to-QP gap correction from a finite vacuum anchor and a bulk limit.

    ``polarization_model="sphere"`` (default): the finite-size term is the
    classical surface polarization of the electron and the hole in a dielectric
    sphere with the bulk eps_inf inside and eps_out outside, averaged over the
    1S envelope, P(R) = F(eps_inf, eps_out) e^2/R (``sphere_polarization_factor``).
    This is the leading self-energy correction found in tight-binding GW for
    nanocrystals (Delerue, Lannoo, Allan, PRL 84, 2457 (2000); PRL 90, 076803
    (2003)).  A residual A (R0/R)^p fixed at the vacuum anchor makes the curve
    reproduce the evGW anchor exactly at R0.  Each band edge has its own curve
    (``anchor_edge_curves``), so the HOMO/LUMO split is also exact at R0.

    ``polarization_model="legacy"``: the earlier kappa/(R + ell) form with
    kappa = 11.52 eV A (1 - 1/eps_inf), kept to reproduce old results.
    """
    if material_name is None:
        logger.warning("  [Warning] Material not specified. Cannot compute GW scaling.")
        return (None, None) if return_details else None
        
    m_name = material_name.upper()
    if m_name not in MATERIAL_DB:
        logger.warning(f"  [Warning] Material {m_name} not found in MATERIAL_DB. GW estimation failed.")
        return (None, None) if return_details else None
        
    entry = MATERIAL_DB[m_name]
    if len(entry) < 9:
        logger.warning(f"  [Warning] MATERIAL_DB entry for {m_name} is outdated. GW estimation failed.")
        return (None, None) if return_details else None
        
    # Extract Bulk Data
    eps_inf = entry[0]
    gap_pbe_bulk = bulk_pbe_gap_dot(m_name)     # bulk PBE gap at the dot's lattice (geometry correction)
    gap_gw_bulk = entry[8]
    
    if gap_gw_bulk == 0.0:
        logger.warning(f"  [Warning] Missing GW bulk gap for {m_name}. GW estimation failed.")
        return (None, None) if return_details else None
        
    # Extract Monomer Data (if available in the tuple)
    has_monomer_data = len(entry) >= 14
    R_mono, gap_pbe_mono, gap_gw_mono = 0.0, 0.0, 0.0
    if has_monomer_data:
        R_mono = entry[9]
        pbe_h, pbe_l = entry[10], entry[11]
        gw_h, gw_l = entry[12], entry[13]
        
        # Compute the gaps on the fly!
        gap_pbe_mono = pbe_l - pbe_h
        gap_gw_mono = gw_l - gw_h

    # Get physical size of the cluster
    metrics = get_cluster_size_metrics(coords, atom_symbols, material_name)
    R_QD_ang = metrics['R_eff_hull']
    
    # --- 1. Bulk GW Correction (Infinite Limit) ---
    delta_bulk_qp = gap_gw_bulk - gap_pbe_bulk
    
    logger.info(f"\n  [Scaled GW Model] Quasiparticle Correction for {m_name}:")
    logger.info(f"    Cluster Radius (R_QD) : {R_QD_ang:.3f} Å")
    anisotropy = float(metrics.get("anisotropy", 1.0))
    if anisotropy > 2.0:
        logger.warning(
            f"    [Size Warning] Principal-axis anisotropy is {anisotropy:.2f}; "
            "the scalar equivalent-volume radius may be a coarse approximation."
        )
    logger.info(f"    Bulk PBE Gap          : {gap_pbe_bulk:.3f} eV")
    logger.info(f"    Bulk GW Gap           : {gap_gw_bulk:.3f} eV")
    logger.info(f"    -> Bulk Shift         : {delta_bulk_qp:+.3f} eV")
   
    ell = float(regularization_length_ang)
    p = float(residual_power)
    if ell < 0.0:
        raise ValueError("regularization_length_ang must be non-negative")
    if p <= 1.0:
        raise ValueError("residual_power must be greater than 1")
    pol_model = str(polarization_model).lower()
    if pol_model not in ("sphere", "legacy"):
        raise ValueError(f"polarization_model must be 'sphere' or 'legacy', got '{polarization_model}'")

    # --- 2. Anchor-constrained finite-size correction ---
    kappa_vac = IMAGE_CHARGE_CONST_EV_ANG * (1.0 - (1.0 / eps_inf))
    kappa_solvent = IMAGE_CHARGE_CONST_EV_ANG * ((1.0 / eps_out) - (1.0 / eps_inf))
    anchor_residual = 0.0
    anchor_gap_shift = None
    radius_used = R_QD_ang
    edges = None

    if pol_model == "sphere":
        F_vac = sphere_polarization_factor(eps_inf, 1.0)
        F_out = sphere_polarization_factor(eps_inf, max(1.0, float(eps_out)))

        def pol_term(F, R):
            return F * COULOMB_EV_ANG / R
    else:
        F_vac = F_out = None

        def pol_term(kappa, R):
            return kappa / (R + ell)
    c_vac = F_vac if pol_model == "sphere" else kappa_vac
    c_out = F_out if pol_model == "sphere" else kappa_solvent

    if has_monomer_data and gap_gw_mono > 0.0:
        anchor_gap_shift = gap_gw_mono - gap_pbe_mono
        anchor_residual = anchor_gap_shift - delta_bulk_qp - pol_term(c_vac, R_mono)
        if R_QD_ang < R_mono - 1.0e-8:
            message = (
                f"Target radius {R_QD_ang:.3f} Å is below the finite-anchor radius "
                f"{R_mono:.3f} Å; extrapolation is disabled."
            )
            if strict:
                raise ValueError(message)
            logger.warning(f"    [Warning] {message} Using R=R0 for the QP model.")
            radius_used = R_mono

        use_residual = str(anchor_residual).lower() not in ("off", "false", "0", "no")
        if not use_residual:
            anchor_residual = 0.0
            residual = 0.0
            res_scale = 0.0
            res_mode = "off"
        else:
            res_scale, res_mode = anchor_residual_scale(m_name, radius_used, dft_gap, residual_scaling, p)
            residual = anchor_residual * res_scale
        sigma_pol = pol_term(c_out, radius_used) + residual
        sigma_pol_vac = pol_term(c_vac, radius_used) + residual
        if pol_model == "sphere":
            edges = anchor_edge_curves(m_name, radius_used, eps_out, residual_power=p, decay=res_scale,
                                       residual_enabled=use_residual)
        logger.info(f"    Finite Vacuum Anchor : R0={R_mono:.3f} Å, gap shift={anchor_gap_shift:+.3f} eV")
        if pol_model == "sphere":
            logger.info(f"    Sphere Polarization  : F = {F_out:.4f} (vacuum {F_vac:.4f}), "
                  f"P(R) = {pol_term(c_out, radius_used):+.3f} eV")
            logger.info(f"    Anchor Residual A    : {anchor_residual:+.3f} eV x {res_scale:.3f} "
                  f"({'E_conf(R)/E_conf(R0)' if res_mode == 'econf' else (f'(R0/R)^p, p={p:.2f}' if use_residual else 'off')})")
        else:
            logger.info(f"    Anchor Residual A    : {anchor_residual:+.3f} eV (ell={ell:.3f} Å, p={p:.3f})")
    else:
        # No finite anchor is available; retain the correct bulk and dielectric limits.
        sigma_pol = pol_term(c_out, radius_used)
        sigma_pol_vac = pol_term(c_vac, radius_used)
        logger.warning("    [Warning] No finite-QD GW anchor; using only the classical polarization term.")

    logger.info(f"    Solvent Dielectric    : eps_out = {eps_out:.2f}, eps_inf = {eps_inf:.2f}")
    logger.info(f"    -> Polarization Shift : {sigma_pol:+.3f} eV")
 
    # --- 4. Total Quasiparticle Scissor ---
    total_scissor = delta_bulk_qp + sigma_pol
    logger.info(f"    ==> Total GW Scissor  : {total_scissor:+.3f} eV\n")

    details = {
        "qp_model": "anchor_scaled_pbe_to_qp_model",
        "material": m_name,
        "cluster_radius_ang": float(R_QD_ang),
        "radius_definition_version": metrics.get("radius_definition_version"),
        "selected_atom_indices": metrics.get("selected_atom_indices"),
        "principal_extents_ang": metrics.get("principal_extents_ang"),
        "anisotropy_ratio": metrics.get("anisotropy_ratio"),
        "surface_offset_ang": metrics.get("surface_offset_ang"),
        "eps_out": float(eps_out),
        "eps_inf": float(eps_inf),
        "bulk_pbe_gap_ev": float(gap_pbe_bulk),
        "bulk_gw_gap_ev": float(gap_gw_bulk),
        "bulk_gw_shift_ev": float(delta_bulk_qp),
        "has_monomer_anchor": bool(has_monomer_data and gap_gw_mono > 0),
        "monomer_radius_ang": float(R_mono) if has_monomer_data else None,
        "monomer_pbe_gap_ev": float(gap_pbe_mono) if has_monomer_data else None,
        "monomer_gw_gap_ev": float(gap_gw_mono) if has_monomer_data else None,
        "anchor_pbe_homo_ev": float(pbe_h) if has_monomer_data else None,
        "anchor_pbe_lumo_ev": float(pbe_l) if has_monomer_data else None,
        "anchor_qp_homo_ev": float(gw_h) if has_monomer_data else None,
        "anchor_qp_lumo_ev": float(gw_l) if has_monomer_data else None,
        "pbe_to_qp_homo_shift_ev": float(gw_h - pbe_h) if has_monomer_data else None,
        "pbe_to_qp_lumo_shift_ev": float(gw_l - pbe_l) if has_monomer_data else None,
        "anchor_gap_shift_ev": float(anchor_gap_shift) if anchor_gap_shift is not None else None,
        "polarization_model": pol_model,
        "polarization_factor_vacuum": float(F_vac) if F_vac is not None else None,
        "polarization_factor_solvent": float(F_out) if F_out is not None else None,
        "regularization_length_ang": ell if pol_model == "legacy" else None,
        "residual_power": p,
        "residual_scaling": res_mode if (has_monomer_data and gap_gw_mono > 0.0) else None,
        "residual_scale": float(res_scale) if (has_monomer_data and gap_gw_mono > 0.0) else None,
        "anchor_residual_ev": float(anchor_residual),
        "radius_used_ang": float(radius_used),
        "kappa_vacuum_ev_ang": float(kappa_vac),
        "kappa_solvent_ev_ang": float(kappa_solvent),
        "finite_size_shift_vacuum_ev": float(sigma_pol_vac),
        "finite_size_shift_solvent_ev": float(sigma_pol),
        "total_scissor_vacuum_ev": float(delta_bulk_qp + sigma_pol_vac),
        "total_scissor_solvent_ev": float(total_scissor),
    }

    if edges is not None:
        # Per-edge two-anchor curves: exact HOMO/LUMO split at R0.
        details.update({
            "f_homo": float(edges["f_homo"]),
            "f_lumo": float(edges["f_lumo"]),
            "edge_split_source": "anchor_edge_curves",
            "edge_residual_homo_ev": float(edges["residual_homo_ev"]),
            "edge_residual_lumo_ev": float(edges["residual_lumo_ev"]),
            "bulk_homo_fraction": float(edges["bulk_homo_fraction"]),
        })

    if return_details:
        return total_scissor, details
    return total_scissor

def compute_delta_xc(material):

    entry = MATERIAL_DB.get(material.upper(), None)

    if entry is None:
        return 0.0

    # Experimental bulk gap
    gap_exp = entry[3]

    # If PBE gap is not present return zero correction
    if len(entry) < 8:
        logger.warning(f"[Δxc] Warning: no PBE gap for {material}. Using Δxc = 0.")
        return 0.0

    gap_pbe = bulk_pbe_gap_dot(material)

    delta_xc = gap_exp - gap_pbe

    logger.info(f"[Δxc] Bulk correction for {material}: {delta_xc:.3f} eV")

    return max(delta_xc, 0.0)

# =====================================================================
# Single-oscillator (Penn) dielectric model shared by the QP layer
# =====================================================================
# Valence (s, p) electrons per formula unit and the formula-unit volume as a
# fraction of the cubic lattice constant cubed.  Semicore d electrons are not
# counted (Phillips-Van Vechten convention).  Zincblende/wurtzite and rocksalt
# have 4 formula units per cubic cell; cubic ABX3 perovskites have one.
PLASMON_DATA = {
    "ZNS": (8, 0.25), "ZNSE": (8, 0.25), "ZNTE": (8, 0.25),
    "CDS": (8, 0.25), "CDSE": (8, 0.25), "CDTE": (8, 0.25),
    "HGS": (8, 0.25), "HGSE": (8, 0.25), "HGTE": (8, 0.25),
    "ALP": (8, 0.25), "ALAS": (8, 0.25), "ALSB": (8, 0.25),
    "GAP": (8, 0.25), "GAAS": (8, 0.25), "GASB": (8, 0.25),
    "INP": (8, 0.25), "INAS": (8, 0.25), "INSB": (8, 0.25),
    "PBS": (10, 0.25), "PBSE": (10, 0.25),
    "CSPBCL3": (26, 1.0), "CSPBBR3": (26, 1.0), "CSPBI3": (26, 1.0),
}

# Polytypes: wurtzite CdSe/CdS (MATERIAL_DB "CDSE_WZ", "CDS_WZ") share the composition-level data of the
# zinc-blende entries; "CDSE_ZB"/"CDS_ZB" are aliases of the zinc-blende "CDSE"/"CDS" (bulk limits differ
# between polytypes: experimental gap, PBE gaps and lattice, MATERIAL_DB_BULK_PBE).
for _wz, _zb in (("CDSE_WZ", "CDSE"), ("CDS_WZ", "CDS")):
    for _tab in (MATERIAL_ELEMENTS, REFRACTIVE_INDEX_DICT, PLASMON_DATA):
        _tab.setdefault(_wz, _tab[_zb])
for _alias, _zb in (("CDSE_ZB", "CDSE"), ("CDS_ZB", "CDS")):
    for _tab in (MATERIAL_DB, MATERIAL_DB_BULK_PBE, QSGW_SOC_LITERATURE, MATERIAL_ELEMENTS, REFRACTIVE_INDEX_DICT,
                 PLASMON_DATA):
        if _zb in _tab:
            _tab.setdefault(_alias, _tab[_zb])
DEFAULT_PLASMON_EV = 15.0
COULOMB_EV_ANG = 14.3996


def valence_plasmon_ev(material_name):
    """Free-electron valence plasmon energy hbar*omega_p = sqrt(4 pi n) (Hartree units).

    n is the valence (s, p) electron density of the bulk crystal from
    ``PLASMON_DATA`` and the lattice constant in ``MATERIAL_DB``.  Materials
    without an entry fall back to ``DEFAULT_PLASMON_EV``.
    """
    m_name = str(material_name).upper() if material_name else "DEFAULT"
    if m_name not in PLASMON_DATA or m_name not in MATERIAL_DB:
        return DEFAULT_PLASMON_EV
    n_val, vol_factor = PLASMON_DATA[m_name]
    a_ang = float(MATERIAL_DB[m_name][2])
    n_bohr = n_val / (vol_factor * a_ang ** 3) * ANG_PER_BOHR ** 3
    return float(np.sqrt(4.0 * np.pi * n_bohr) * HA_TO_EV)


def penn_gap_ev(material_name, eps_inf=None):
    """Penn gap E_P = hbar*omega_p / sqrt(eps_inf - 1) of the bulk material."""
    m_name = str(material_name).upper() if material_name else "DEFAULT"
    if eps_inf is None:
        eps_inf = float(MATERIAL_DB.get(m_name, MATERIAL_DB["DEFAULT"])[0])
    return float(valence_plasmon_ev(m_name) / np.sqrt(max(1.0e-3, float(eps_inf) - 1.0)))


def penn_eps_eff(eps_inf, delta_gap_ev, penn_gap):
    """Size-dependent interior dielectric constant of the Penn model.

    In the Penn model eps - 1 = (hbar omega_p / E_P)^2.  Confinement opens the
    average gap, E_P -> E_P + Delta E, so

        eps(R) = 1 + (eps_inf - 1) * [E_P / (E_P + Delta E)]^2

    (Penn, Phys. Rev. 128, 2093 (1962); Tsu, Babic, Ioriatti, J. Appl. Phys.
    82, 1327 (1997)).  Delta E is the confinement opening of the gap that
    enters the polarizability.
    """
    dE = max(0.0, float(delta_gap_ev))
    E_P = max(1.0e-3, float(penn_gap))
    return float(1.0 + (float(eps_inf) - 1.0) * (E_P / (E_P + dE)) ** 2)


@functools.lru_cache(maxsize=256)
def sphere_polarization_factor(eps_in, eps_out, n_terms=4000, n_grid=4000):
    """Dimensionless surface-polarization factor F of a dielectric sphere.

    A charge at radius r inside a sphere of radius R (eps_in inside, eps_out
    outside) has the self-image energy (Boettcher; Brus, J. Chem. Phys. 80,
    4403 (1984))

        Sigma(r) = e^2/(2R) sum_n (eps_in - eps_out)(n + 1)
                   / [eps_in (n eps_in + (n + 1) eps_out)] (r/R)^(2n).

    Averaged over the 1S envelope j0(pi r/R)^2 for the electron and for the
    hole, the gap opening is Sigma_e + Sigma_h = F e^2 / R.  The n = 0 term is
    the Born term (1/eps_out - 1/eps_in); the higher multipoles make F larger.
    F = 0 for eps_in = eps_out.
    """
    ei = float(eps_in)
    eo = float(eps_out)
    if abs(ei - eo) < 1.0e-12:
        return 0.0
    x = (np.arange(n_grid) + 0.5) / n_grid
    rho = np.sin(np.pi * x) ** 2
    rho /= rho.sum()
    n = np.arange(n_terms, dtype=float)
    c = (ei - eo) * (n + 1.0) / (ei * (n * ei + (n + 1.0) * eo))
    # sum_n c_n x^(2n) evaluated as a power series in y = x^2
    y = x ** 2
    series = np.zeros_like(x)
    for cn in c[::-1]:
        series = series * y + cn
    return float(np.sum(rho * series))


def anchor_residual_scale(material_name, radius_ang, dft_gap=None, mode="econf", residual_power=2.0):
    """Size scaling of the non-classical anchor residual (1 at the anchor, 0 in bulk).

    mode 'econf' (default): E_conf(R) / E_conf(R0), with E_conf the PBE confinement
    energy E_g^PBE(cluster) - E_g^PBE(bulk).  If the residual is band stretching
    (an energy-dependent bulk GW correction), it is proportional to how far the
    confined levels lie from the band edges, which each cluster's own PBE gap
    measures; no radius or power law enters.  Clipped to [0, 1].
    mode 'power': (R0/R)^p with R clamped to R >= R0 (the earlier form); also the
    fallback when no DFT gap is available.
    """
    m_name = str(material_name).upper()
    entry = MATERIAL_DB.get(m_name)
    if entry is None or len(entry) < 14:
        return 0.0, "none"
    R0 = float(entry[9])
    if str(mode).lower() == "econf" and dft_gap is not None:
        eg_bulk = bulk_pbe_gap_dot(m_name)
        e0 = (float(entry[11]) - float(entry[10])) - eg_bulk
        if e0 > 1.0e-6:
            return float(np.clip((float(dft_gap) - eg_bulk) / e0, 0.0, 1.0)), "econf"
    R = max(float(radius_ang), R0)
    return float((R0 / R) ** float(residual_power)), "power"


def anchor_bulk_homo_fraction(material_name, default=0.5):
    """Return the HOMO share of the bulk GW opening from the anchor monomer.

    In the literature (e.g. Hinuma et al., Phys. Rev. B 90, 155405 (2014);
    Schleife et al., Phys. Rev. B 73, 245212 (2006)), the bulk GW gap opening
    is asymmetric between valence and conduction bands: the VBM typically takes
    ~42-44% and the CBM ~56-58% for II-VI semiconductors like CdSe, owing to the
    localized Se 4p / Cd 4d character of the valence states versus the diffuse
    Cd 5s conduction states.

    Rather than imposing an arbitrary 50:50 midpoint split, we take the fraction
    directly from the anchor monomer's PBE->GW shift ratio::

        f_b = d_h0 / (d_h0 + d_l0)

    where d_h0 = -(GW_HOMO - PBE_HOMO) and d_l0 = GW_LUMO - PBE_LUMO.
    For CdSe: 1.3981 / 3.3949 = 0.4118 (41.2% HOMO / 58.8% LUMO), which matches
    first-principles bulk GW band-alignment benchmarks (Hinuma 2014, Schleife 2006:
    42.5-43.3% HOMO) to within 1-2%.
    Falls back to `default` (0.5) if no anchor data is available.
    """
    m_name = str(material_name).upper() if material_name else "DEFAULT"
    entry = MATERIAL_DB.get(m_name)
    if entry is None or len(entry) < 14:
        return float(default)
    d_h0 = -(float(entry[12]) - float(entry[10]))
    d_l0 = float(entry[13]) - float(entry[11])
    if d_h0 + d_l0 <= 0.0:
        return float(default)
    return float(d_h0 / (d_h0 + d_l0))


def anchor_edge_curves(material_name, radius_ang, eps_out, residual_power=2.0,
                       bulk_homo_fraction=None, decay=None, residual_enabled=True):
    """Per-edge two-anchor curves of the PBE-to-QP shift (eV, both positive).

    Each band edge is interpolated between the bulk GW limit and the finite
    evGW anchor (the monomer in ``MATERIAL_DB``)::

        dE_h(R) = f_b dGW_bulk + P(R)/2 + A_h (R0/R)^p      (HOMO moves down)
        dE_l(R) = (1 - f_b) dGW_bulk + P(R)/2 + A_l (R0/R)^p (LUMO moves up)

    P(R) = F(eps_inf, eps_out) e^2 / R is the classical surface polarization of
    the electron plus the hole (``sphere_polarization_factor``), which is
    symmetric in electron and hole.  A_h and A_l are fixed so that each edge
    reproduces the vacuum anchor exactly at R0; they carry the non-classical
    HOMO/LUMO asymmetry and decay as (R0/R)^p.  f_b is the HOMO share of the
    bulk GW opening (defaults to the anchor monomer's fraction, ~41.2% for CdSe).
    Returns None without anchor data.
    """
    m_name = str(material_name).upper()
    entry = MATERIAL_DB.get(m_name)
    if entry is None or len(entry) < 14:
        return None
    eps_inf = float(entry[0])
    d_bulk = float(entry[8]) - bulk_pbe_gap_dot(m_name)
    R0 = float(entry[9])
    d_h0 = -(float(entry[12]) - float(entry[10]))
    d_l0 = float(entry[13]) - float(entry[11])
    if d_h0 + d_l0 <= 0.0:
        return None
    R = max(float(radius_ang), R0)
    F_vac = sphere_polarization_factor(eps_inf, 1.0)
    F_out = sphere_polarization_factor(eps_inf, max(1.0, float(eps_out)))
    P0 = F_vac * COULOMB_EV_ANG / R0
    P = F_out * COULOMB_EV_ANG / R
    fb = anchor_bulk_homo_fraction(m_name) if bulk_homo_fraction is None else float(bulk_homo_fraction)
    if not residual_enabled:
        A_h = 0.0
        A_l = 0.0
        decay = 0.0
    else:
        A_h = d_h0 - fb * d_bulk - 0.5 * P0
        A_l = d_l0 - (1.0 - fb) * d_bulk - 0.5 * P0
        if decay is None:
            decay = (R0 / R) ** float(residual_power)
    d_h = fb * d_bulk + 0.5 * P + A_h * decay
    d_l = (1.0 - fb) * d_bulk + 0.5 * P + A_l * decay
    return {
        "radius_used_ang": R, "anchor_radius_ang": R0, "bulk_shift_ev": d_bulk,
        "polarization_factor_vacuum": F_vac, "polarization_factor_solvent": F_out,
        "polarization_ev": P, "polarization_anchor_ev": P0,
        "residual_homo_ev": A_h, "residual_lumo_ev": A_l,
        "shift_homo_ev": d_h, "shift_lumo_ev": d_l,
        "f_homo": d_h / (d_h + d_l), "f_lumo": d_l / (d_h + d_l),
        "bulk_homo_fraction": fb,
    }


def get_cluster_size_metrics(coords_ang, atom_symbols=None, material_name=None):
    """Return rotation/translation-invariant size metrics for selected inorganic atoms."""
    coords = np.asarray(coords_ang, dtype=float)
    selected_indices = np.arange(len(coords), dtype=int)

    if atom_symbols is not None and material_name is not None:
        m_name = material_name.upper()
        if 'MATERIAL_ELEMENTS' in globals() and m_name in MATERIAL_ELEMENTS:
            core_elements = [el.lower() for el in MATERIAL_ELEMENTS[m_name]]
            selected_indices = np.array(
                [i for i, sym in enumerate(atom_symbols) if sym.lower() in core_elements],
                dtype=int,
            )
            core_coords = coords[selected_indices]
            if len(core_coords) > 0:
                coords = np.array(core_coords)

    if len(coords) < 2:
        return {
            'R_eff_hull': 1.0, 'diameter_hull': 2.0,
            'selected_atom_indices': selected_indices.tolist(),
            'radius_definition_version': 'convex_hull_equivalent_volume_v1',
            'principal_extents_ang': [0.0, 0.0, 0.0], 'anisotropy_ratio': 1.0,
        }

    centered = coords - np.mean(coords, axis=0)
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    principal_coords = centered @ vh.T
    extents = np.ptp(principal_coords, axis=0)
    positive_extents = extents[extents > 1.0e-12]
    anisotropy = float(np.max(positive_extents) / np.min(positive_extents)) if len(positive_extents) else 1.0

    try:
        hull = ConvexHull(coords)
        volume = float(hull.volume)
        radius_geom = (3.0 * volume / (4.0 * np.pi)) ** (1.0 / 3.0)
    except Exception:
        # Degenerate (linear/coplanar) selections have no 3-D hull volume.
        radius_geom = float(np.max(np.linalg.norm(centered, axis=1)))

    # VDW_SURFACE_ANG is a diameter allowance in the legacy model; apply half once.
    surface_offset = 0.5 * VDW_SURFACE_ANG
    radius_eff = radius_geom + surface_offset
    return {
        'R_eff_hull': float(radius_eff),
        'diameter_hull': float(2.0 * radius_eff),
        'hull_volume_ang3': float(volume) if 'volume' in locals() else 0.0,
        'surface_offset_ang': float(surface_offset),
        'selected_atom_indices': selected_indices.tolist(),
        'radius_definition_version': 'convex_hull_equivalent_volume_v1',
        'principal_extents_ang': [float(x) for x in extents],
        'anisotropy_ratio': anisotropy,
    }

# Radius of the dielectric sphere and of the confinement models: 'saxs' (default; SAXS-equivalent
# sphere of the inorganic electron density, qdex.cluster_size) or 'hull' (core convex hull + 1.25 A).
RADIUS_DEFINITION = "saxs"
_RADIUS_CACHE = {}


def qd_radius(coords, atom_symbols, material_name=None, definition=None):
    """Radius (A) of the dot used by the sphere reaction field and the confinement models."""
    how = str(definition or RADIUS_DEFINITION).lower()
    coords = np.asarray(coords, dtype=float)
    if how == "hull":
        return float(get_cluster_size_metrics(coords, atom_symbols, material_name)["R_eff_hull"])
    key = (coords.tobytes(), tuple(atom_symbols), str(material_name))
    if key not in _RADIUS_CACHE:
        from qdex.cluster_size import cluster_size
        size = cluster_size(coords, atom_symbols, material_name)
        _RADIUS_CACHE.clear()
        _RADIUS_CACHE[key] = 5.0 * float(size.get("d_saxs_nm", size["d_hull_nm"]))
    return _RADIUS_CACHE[key]


# MNOK interaction gamma_AB = (r^beta + a_AB^beta)^(-1/beta), a_AB = (1/(s eta_A) + 1/(s eta_B)) / 2.
# Set from integrals.mnok_exponent, mnok_exponent_exchange and mnok_onsite.
MNOK_EXPONENT = 2.0        # beta for all MNOK interactions: 2 = Ohno-Klopman (default), 1 = Mataga-Nishimoto
MNOK_EXPONENT_K = None     # beta of the exchange gamma (K^x); None = MNOK_EXPONENT
MNOK_ONSITE_SCALE = 2.0    # s: 2 -> gamma_AA = IP - EA (default, benchmarked against xs); 1 -> eta = (IP-EA)/2


def set_mnok_options(exponent=2.0, exponent_exchange=None, onsite="ip_ea", verbose=True):
    """Set the MNOK exponent(s) and on-site convention used by every MNOK builder."""
    global MNOK_EXPONENT, MNOK_EXPONENT_K, MNOK_ONSITE_SCALE
    MNOK_EXPONENT = float(exponent if exponent not in (None, "") else 2.0)
    MNOK_EXPONENT_K = None if exponent_exchange in (None, "", "none") else float(exponent_exchange)
    onsite = str(onsite or "ip_ea").lower()
    if onsite not in ("eta", "ip_ea"):
        raise ValueError(f"mnok_onsite must be 'eta' or 'ip_ea', not '{onsite}'")
    MNOK_ONSITE_SCALE = 2.0 if onsite == "ip_ea" else 1.0
    if verbose and (MNOK_EXPONENT, MNOK_EXPONENT_K, MNOK_ONSITE_SCALE) != (2.0, None, 2.0):
        k = MNOK_EXPONENT_K if MNOK_EXPONENT_K is not None else MNOK_EXPONENT
        logger.info(f"  [MNOK] exponent {MNOK_EXPONENT:g}, exchange exponent {k:g}, "
              f"on-site {'IP - EA' if MNOK_ONSITE_SCALE == 2.0 else 'eta'}")


_KERNEL_DESCRIPTION = {
    "qp": "W of the QP model (the same W as in the QP correction)",
    "resta": "bulk Resta W_AB = S_eps_inf(R_AB) gamma_AB",
    "resta-sphere": "bulk Resta W + dielectric-sphere reaction field",
    "dim": "DIM W_AB = S_eps_AB(R_AB) gamma_AB",
    "rpa": "RPA W = (1 + gamma Pi0)^-1 gamma of the active space",
    "sbse": "sBSE screening (Cho, Bintrim, Berkelbach)",
    "bse": "uniform alpha gamma / eps_inf",
}


def format_integrals_block(representation, charges, kernel, symbols, stda_info=None, stda_ax_source=None,
                           include_direct=True, include_exchange=True, hubbard_beta=0.0,
                           eta_dict=HARDNESS_DICT):
    """Printable description of the two-electron integrals actually used in the run."""
    rep = str(representation).lower()
    kern = str(kernel).lower() if kernel is not None else "none"
    elements = sorted({str(s).capitalize() for s in symbols})
    lines = ["", "--- Two-electron integrals ---"]
    if kern == "stda":
        ax = stda_info["ax"]
        src = f" ({stda_ax_source})" if stda_ax_source else ""
        a_rule = (" = 1.42 + 0.48 a_x" if abs(stda_info["alpha_K"] - (STDA_ALPHA1 + STDA_ALPHA2 * ax)) < 1e-12
                  else " (set)")
        b_rule = (" = 0.20 + 1.83 a_x" if abs(stda_info["beta_J"] - (STDA_BETA1 + STDA_BETA2 * ax)) < 1e-12
                  else " (set)")
        lines += [
            "  Representation : sTDA (Grimme, J. Chem. Phys. 138, 244104 (2013)), Loewdin transition charges",
            f"  Exchange  K^x  : gamma^K_AB = (R^a + eta_AB^-a)^(-1/a),        a{a_rule}"
            f" = {stda_info['alpha_K']:.3f}",
            f"  Direct    K^d  : gamma^J_AB = (R^b + (a_x eta_AB)^-b)^(-1/b),  b{b_rule}"
            f" = {stda_info['beta_J']:.3f}",
            f"  a_x            : {ax:.3f}{src}",
            "  eta_AB         : (eta_A + eta_B)/2 with eta = IP - EA (R, eta in atomic units)",
        ]
        onsite = {e: 2.0 * eta_dict.get(e.lower(), 5.0) for e in elements}
        lines.append("  eta per element: " + ", ".join(f"{e} {v:.2f} eV" for e, v in onsite.items()))
        return "\n".join(lines)
    k_exp = MNOK_EXPONENT_K if MNOK_EXPONENT_K is not None else MNOK_EXPONENT
    names = {1.0: "Mataga-Nishimoto", 2.0: "Ohno-Klopman"}
    if rep.startswith("xs"):
        lines.append("  Representation : xs, exact AO density-pair integrals (mu mu|nu nu) (Libint2, ZDO)")
    else:
        s = MNOK_ONSITE_SCALE
        lines += [
            f"  Representation : MNOK, atom pairs, {str(charges).capitalize()} transition charges",
            "  gamma_AB       = (R^beta + a_AB^beta)^(-1/beta),  a_AB = (1/gamma_AA + 1/gamma_BB)/2",
            f"  beta (direct)  : {MNOK_EXPONENT:g}" + (f" ({names[MNOK_EXPONENT]})" if MNOK_EXPONENT in names else ""),
            f"  beta (exchange): {k_exp:g}" + (f" ({names[k_exp]})" if k_exp in names else ""),
            "  On-site        : " + ("gamma_AA = IP - EA = 2 eta_A" if s == 2.0 else "gamma_AA = eta_A = (IP - EA)/2")
            + "  (integrals.mnok_onsite: " + ("ip_ea" if s == 2.0 else "eta") + ")",
            "  gamma_AA       : " + ", ".join(f"{e} {s * eta_dict.get(e.lower(), 5.0):.2f} eV" for e in elements),
        ]
        if hubbard_beta:
            lines.append(f"  Hubbard stiffening of the diagonal: beta = {hubbard_beta:g}")
    lines.append("  Exchange  K^x  : " + ("bare interaction (unscreened)" if include_exchange else "off (triplets or include_exchange: false)"))
    if include_direct:
        lines.append(f"  Direct    K^d  : {_KERNEL_DESCRIPTION.get(kern, kern)}  (excitations.kernel: {kern})")
    else:
        lines.append("  Direct    K^d  : off")
    return "\n".join(lines)


# Geometry correction of the bulk QP shift (quasiparticles.bulk_geometry).
# The bulk shift splits into a self-energy and a geometry part,
#     Delta_bulk = [E_g^QSGW(a_exp) - E_g^PBE(a_exp)] + [E_g^PBE(a_exp) - E_g^PBE(a_dot)]
#                =  Delta_Sigma                        +  Delta_geom.
# A dot relaxed with PBE has (almost) the PBE lattice: its PBE gap carries the gap change of bulk PBE
# between a_exp and its own lattice a_dot, which Delta_geom removes. Delta_geom is a geometry effect,
# not a self-energy: the bulk vertex correction scales Delta_Sigma only. With
# E_g^PBE(a_dot) = E_g^PBE(a_exp) - s * [E_g^PBE(a_exp) - E_g^PBE(a_PBE)], the strain fraction s of the dot
# is measured on its interior: s = (d_dot/g - a_exp) / (a_PBE - a_exp), with d the cation-anion bond and
# g the bond per lattice constant (zinc blende sqrt(3)/4, rock salt 1/2); in the perovskites d is the
# B-B (Pb-Pb) distance, g = 1, because octahedral tilts lengthen the Pb-X bond at fixed volume (by 1.6 %
# in orthorhombic CsPbBr3) while the Pb-Pb distance stays the pseudo-cubic lattice constant. s is capped
# at 1: frames of a 300 K trajectory are expanded beyond the 0 K PBE lattice by thermal expansion, which
# the room-temperature experiment has as well.
#   strain : s measured on the dot (default; s = 1 when it cannot be measured)
#   full   : s = 1, the dot at the bulk PBE lattice
#   none   : s = 0, no geometry correction (the dot taken at a_exp)
BULK_GEOMETRY = "none"
DOT_STRAIN = {}            # material -> strain fraction s of the dot of this run (set_dot_strain)
# bond per lattice constant: zinc blende sqrt(3)/4 a; wurtzite sqrt(3/8) a_hex (ideal c/a and u); rock salt
# a/2; perovskite: the B-B distance, a
_BOND_PER_LATTICE = {"zb": np.sqrt(3.0) / 4.0, "wz": np.sqrt(3.0 / 8.0), "rs": 0.5, "perovskite": 1.0}


def _structure_of(m_name):
    if m_name.startswith("CSPB") or m_name in ("MAPBI3", "FAPBI3"):
        return "perovskite"
    if m_name.endswith("_WZ"):
        return "wz"
    if m_name in ("PBS", "PBSE", "PBTE"):
        return "rs"
    return "zb"


def set_bulk_geometry(mode="strain", verbose=True):
    """Set the geometry correction of the bulk QP shift ('strain', 'full' or 'none')."""
    global BULK_GEOMETRY
    mode = str(mode or "none").lower()
    if mode not in ("strain", "full", "none"):
        raise ValueError(f"bulk_geometry must be 'strain', 'full' or 'none', not '{mode}'")
    BULK_GEOMETRY = mode
    DOT_STRAIN.clear()
    if verbose and mode != "none":
        logger.info(f"  [QP] Bulk geometry correction: {mode}")


def dot_lattice_strain(material_name, atom_symbols, coords_ang, interior=0.6):
    """Strain fraction s of a dot between a_exp (s = 0) and the bulk PBE lattice (s = 1).

    Median of the nearest cation-anion bonds whose atoms lie within `interior` x the radius of the
    inorganic core (widened when there are fewer than 8 such bonds). Returns (s, info); s is None
    when the material has no bulk PBE lattice data or the dot has no cation-anion bonds."""
    from scipy.spatial import cKDTree
    m_name = str(material_name).upper()
    data = MATERIAL_DB_BULK_PBE.get(m_name)
    elements = MATERIAL_ELEMENTS.get(m_name)
    if not data or data.get("a_exp") is None or not elements:
        return None, {}
    structure = _structure_of(m_name)
    # perovskites: B-B (Pb-Pb) distances, insensitive to the octahedral tilts
    cation, anion = (elements[1], elements[1]) if structure == "perovskite" else (elements[0], elements[1])
    syms = np.array([str(s).capitalize() for s in atom_symbols])
    X = np.asarray(coords_ang, float)
    ic, ia = np.flatnonzero(syms == cation), np.flatnonzero(syms == anion)
    if ic.size == 0 or ia.size == 0:
        return None, {}
    g = _BOND_PER_LATTICE[structure]
    d_exp, d_pbe = g * data["a_exp"], g * data["a_pbe"]
    core = np.concatenate([ic, ia])
    centre = X[core].mean(axis=0)
    r = np.linalg.norm(X - centre, axis=1)
    R = float(r[core].max())
    tree = cKDTree(X[ia])
    same = cation == anion
    k = min(8 + same, ia.size)
    dist, nb = tree.query(X[ic], k=k)
    dist, nb = np.atleast_2d(dist), np.atleast_2d(nb)
    keep = (dist < 1.25 * d_pbe) & (dist > 0.5 * d_exp)            # (B-B: not the atom itself)
    bond_d = dist[keep]
    bond_r = np.maximum(np.repeat(r[ic][:, None], k, axis=1)[keep], r[ia][nb[keep]])
    if bond_d.size == 0:
        return None, {}
    frac_used = interior
    while frac_used < 1.0 and np.count_nonzero(bond_r <= frac_used * R) < 8:
        frac_used = min(1.0, frac_used + 0.1)
    sel = bond_r <= frac_used * R + 1e-9
    d_dot = float(np.median(bond_d[sel]))
    s = (d_dot - d_exp) / (d_pbe - d_exp)
    info = dict(bond=f"{cation}-{anion}", bond_dot_ang=d_dot, bond_exp_ang=float(d_exp), bond_pbe_ang=float(d_pbe),
                n_bonds=int(np.count_nonzero(sel)), interior_fraction=float(frac_used), strain_fraction=float(s))
    return float(s), info


def structure_is_bb(m_name):
    return _structure_of(str(m_name).upper()) == "perovskite"


def set_dot_strain(material_name, atom_symbols, coords_ang):
    """Measure the strain fraction of the dot of this run (bulk_geometry 'strain'); returns (s, info)."""
    m_name = str(material_name).upper() if material_name else "DEFAULT"
    if BULK_GEOMETRY != "strain" or m_name not in MATERIAL_DB_BULK_PBE:
        return None, {}
    s, info = dot_lattice_strain(m_name, atom_symbols, coords_ang)
    if s is None:
        logger.warning(f"  [QP] Bulk geometry: no {m_name} cation-anion bonds found; the dot is taken at the bulk "
                       f"PBE lattice (s = 1).")
        return None, {}
    if s > 1.0:
        logger.info(f"  [QP] Bulk geometry: strain fraction {s:.2f} > 1 (a thermal frame?); the expansion beyond "
                    f"the PBE lattice is taken as thermal, as in the room-temperature experiment: s = 1.")
    elif s < -0.25:
        logger.warning(f"  [QP] Bulk geometry: strain fraction {s:.2f} < 0: the dot is smaller than at a_exp. "
                       f"Was it relaxed with PBE? s = 0.")
    DOT_STRAIN[m_name] = float(np.clip(s, 0.0, 1.0))
    logger.info(f"  [QP] Dot lattice: interior {info['bond']} {'distance' if structure_is_bb(m_name) else 'bond'} {info['bond_dot_ang']:.4f} A over {info['n_bonds']} "
                f"bonds (bulk {info['bond_exp_ang']:.4f} A at a_exp, {info['bond_pbe_ang']:.4f} A at a_PBE): "
                f"strain fraction s = {s:.3f}")
    return DOT_STRAIN[m_name], info


def bulk_geometry_shift(material_name):
    """Delta_geom = s * [E_g^PBE(a_exp) - E_g^PBE(a_PBE)] (spin-free) for the dot of this run. Returns (shift, s)."""
    m_name = str(material_name).upper() if material_name else "DEFAULT"
    data = MATERIAL_DB_BULK_PBE.get(m_name)
    if BULK_GEOMETRY == "none" or not data or data.get("gap_exp_lattice") is None:
        return 0.0, 0.0
    s = 1.0 if BULK_GEOMETRY == "full" else DOT_STRAIN.get(m_name, 1.0)
    return float(s * (data["gap_exp_lattice"] - data["gap_pbe_lattice"])), float(s)


def bulk_pbe_gap_dot(material_name):
    """Bulk PBE gap at the lattice of the dot of this run: the bulk limit of its PBE gap
    (MATERIAL_DB index 7, at a_exp, minus the geometry correction)."""
    m_name = str(material_name).upper() if material_name else "DEFAULT"
    entry = MATERIAL_DB.get(m_name, MATERIAL_DB["DEFAULT"])
    return float(entry[7]) - bulk_geometry_shift(m_name)[0]


# Vertex correction of the bulk QSGW gap (quasiparticles.bulk_vertex).
# QSGW overestimates bulk gaps because its W lacks the electron-hole (ladder) vertex, which makes the
# polarizability too small; the bulk fix is Delta_bulk -> factor * Delta_bulk (factor 0.8). In a dot the
# polarizability is reduced by confinement, and with it the missing vertex part.
#   none   : pure QSGW correction (default)
#   full   : the bulk vertex correction at every size, Delta = factor * Delta_QSGW
#   scaled : the correction times the fraction of bulk screening the dot keeps,
#            f = (eps_eff - 1)/(eps_inf - 1), eps_eff from the Penn model at the DFT gap
#   factor: a number (0.8: the usual bulk QSGW80 scaling) or 'material': the factor that puts the bulk
#   limit on the experimental gap of the material,
#       a_m = (E_exp - E_g^PBE+SOC(a_exp)) / Delta_Sigma,
#   so that PBE+SOC + a_m Delta_Sigma = E_exp in the bulk. It absorbs every difference between 0 K QSGW
#   and the measured gap (missing vertex, zero-point and thermal renormalisation), so E_exp must be at the
#   temperature, and E_g^PBE+SOC for the structure, of the experiments the dots are compared with
#   (room temperature; EXPERIMENTAL_BULK_GAP). With 'scaled' the dot goes from QSGW (molecular limit)
#   to the experimental gap (bulk) with the Penn fraction f. Same bulk limit as the g-xTB route.
BULK_VERTEX = "none"
BULK_VERTEX_FACTOR = 0.8

# Experimental bulk gaps for bulk_vertex_factor: material, the same transition as the QSGW and PBE gaps
# of MATERIAL_DB_BULK_PBE (signed). Room temperature, MATERIAL_DB index 3, except the Hg compounds:
# E0 = E(Gamma6) - E(Gamma8) from the representative experimental values in Svane et al., Phys. Rev. B 84,
# 205205 (2011), Table I (HgS uncertain: experiments disagree on its sign). Notes on the phase:
# CdS 2.42 and CdSe 1.74 eV are the wurtzite values (zinc blende about 0.05-0.07 eV lower); the lead
# halide perovskites are measured in their room-temperature (tilted) phases, while the PBE and QSGW gaps
# are for the cubic structure: the PBE+SOC partner of E_exp is then the gap of the measured structure
# (gap_soc_exp_structure in MATERIAL_DB_BULK_PBE; CsPbBr3 only so far), with Delta_Sigma from the cubic one.
EXPERIMENTAL_E0 = {"HGS": -0.11, "HGSE": -0.20, "HGTE": -0.30}


def experimental_bulk_gap(material_name):
    m_name = str(material_name).upper()
    if m_name in EXPERIMENTAL_E0:
        return float(EXPERIMENTAL_E0[m_name])
    return float(MATERIAL_DB[m_name][3])


def material_vertex_factor(material_name):
    """a_m = (E_exp - E_PBE+SOC(a_exp)) / Delta_Sigma; None without the bulk data. Returns (a_m, info)."""
    m_name = str(material_name).upper()
    b = MATERIAL_DB_BULK_PBE.get(m_name, {})
    if b.get("delta_sigma") is None or b.get("gap_soc_exp_lattice") is None:
        return None, {}
    e_exp = experimental_bulk_gap(m_name)
    partner = float(b.get("gap_soc_exp_structure", b["gap_soc_exp_lattice"]))   # PBE+SOC of the measured structure
    a = (e_exp - partner) / b["delta_sigma"]
    return float(a), dict(bulk_gap_exp_ev=e_exp, bulk_gap_pbe_soc_ev=partner,
                          bulk_delta_sigma_ev=float(b["delta_sigma"]))


def set_bulk_vertex(mode="none", factor=0.8, verbose=True):
    """Set the vertex correction of the bulk QP shift ('none', 'full' or 'scaled') and the bulk factor
    (a number in (0, 1], or 'material' for the per-material factor a_m)."""
    global BULK_VERTEX, BULK_VERTEX_FACTOR
    mode = str(mode or "none").lower()
    if mode not in ("none", "full", "scaled"):
        raise ValueError(f"bulk_vertex must be 'none', 'full' or 'scaled', not '{mode}'")
    if isinstance(factor, str) and factor.strip().lower() == "material":
        factor = "material"
    else:
        factor = float(0.8 if factor in (None, "") else factor)
        if not 0.0 < factor <= 1.0:
            raise ValueError(f"bulk_vertex_factor must be in (0, 1] or 'material', not {factor}")
    BULK_VERTEX, BULK_VERTEX_FACTOR = mode, factor
    if verbose and mode != "none":
        logger.info(f"  [QP] Bulk vertex correction: {mode} (bulk factor {factor if factor == 'material' else f'{factor:g}'})")


# Residual of the bulk reference (quasiparticles.bulk_residual).
#   none         : the bulk limit is PBE+SOC + a Delta_Sigma (0 K, frozen lattice; QSGW with the vertex factor)
#   experimental : plus a constant, size-independent residual that puts the bulk limit on the measured gap,
#                  delta_res = E_exp - (E_g^PBE+SOC + a_bulk Delta_Sigma), a_bulk the vertex factor in the bulk
#                  (1 without vertex correction). The split model: the vertex part fades with the dot's
#                  screening (bulk_vertex: scaled), the residual - thermal and zero-point renormalisation and
#                  what else separates 0 K QSGW80 from the room-temperature gap - does not.
BULK_RESIDUAL = "none"


def set_bulk_residual(mode="none", verbose=True):
    """Set the residual of the bulk reference: 'none' or 'experimental' (split model)."""
    global BULK_RESIDUAL
    mode = str(mode or "none").lower()
    if mode not in ("none", "experimental"):
        raise ValueError(f"bulk_residual must be 'none' or 'experimental', not '{mode}'")
    BULK_RESIDUAL = mode
    if verbose and mode != "none":
        logger.info(f"  [QP] Bulk residual: {mode} (bulk limit on the room-temperature experimental gap)")


def bulk_residual_shift(material_name, a_bulk):
    """delta_res = E_exp - (E_PBE+SOC of the measured structure + a_bulk Delta_Sigma); (0, {}) without data."""
    m_name = str(material_name).upper()
    b = MATERIAL_DB_BULK_PBE.get(m_name, {})
    if BULK_RESIDUAL == "none" or b.get("delta_sigma") is None or b.get("gap_soc_exp_lattice") is None:
        return 0.0, {}
    e_exp = experimental_bulk_gap(m_name)
    partner = float(b.get("gap_soc_exp_structure", b["gap_soc_exp_lattice"]))
    res = e_exp - (partner + a_bulk * float(b["delta_sigma"]))
    return float(res), dict(bulk_gap_exp_ev=e_exp, bulk_gap_pbe_soc_ev=partner)


def _vertex_factor(m_name):
    """The bulk factor of this run for a material: (a, info)."""
    if BULK_VERTEX_FACTOR != "material":
        return float(BULK_VERTEX_FACTOR), {}
    a, info = material_vertex_factor(m_name)
    if a is None:
        logger.warning(f"  [QP] bulk_vertex_factor: material has no experimental/PBE+SOC data for {m_name}; using 0.8.")
        return 0.8, {}
    if a > 1.0:
        logger.warning(f"  [QP] Material vertex factor of {m_name} is {a:.3f} > 1: the bulk QSGW gap is below the "
                       f"experimental one, so the factor is not a vertex correction (check the structure and "
                       f"temperature of the reference).")
    return a, info


def bulk_qp_shift(material_name, dft_gap=None):
    """Bulk QP correction added to the PBE orbital energies, with the vertex correction of BULK_VERTEX.

    Returns (shift_ev, info). Without bulk GW data the shift is 0.
    """
    m_name = str(material_name).upper() if material_name else "DEFAULT"
    entry = MATERIAL_DB.get(m_name)
    a_vertex, a_info = _vertex_factor(m_name) if BULK_VERTEX != "none" else (
        (0.8 if BULK_VERTEX_FACTOR == "material" else float(BULK_VERTEX_FACTOR)), {})
    info = {"bulk_vertex": BULK_VERTEX, "bulk_vertex_factor": float(a_vertex),
            "bulk_vertex_factor_mode": "material" if BULK_VERTEX_FACTOR == "material" else "fixed", **a_info}
    if entry is None or len(entry) < 9:
        info.update(bulk_qsgw_shift_ev=0.0, bulk_vertex_fraction=0.0, bulk_vertex_correction_ev=0.0)
        return 0.0, info
    eps_inf, pbe_gap, gw_gap = float(entry[0]), float(entry[7]), float(entry[8])
    d_qsgw = gw_gap - pbe_gap                        # Delta_Sigma: QSGW - PBE, both at a_exp
    d_geom, s = bulk_geometry_shift(m_name)          # PBE gap at a_exp - PBE gap at the dot's lattice
    pbe_gap_dot = pbe_gap - d_geom
    delta_v = (1.0 - a_vertex) * d_qsgw
    eps_eff = None
    if BULK_VERTEX == "none":
        frac = 0.0
    elif BULK_VERTEX == "full" or dft_gap is None:
        frac = 1.0
    else:
        d_e_conf = max(0.0, float(dft_gap) - pbe_gap_dot)
        eps_eff = penn_eps_eff(eps_inf, d_e_conf, penn_gap_ev(m_name, eps_inf))
        frac = float(np.clip((eps_eff - 1.0) / max(eps_inf - 1.0, 1e-12), 0.0, 1.0))
    d_sigma = d_qsgw - delta_v * frac if frac else d_qsgw
    d_res, r_info = bulk_residual_shift(m_name, a_vertex if BULK_VERTEX != "none" else 1.0)
    if d_res:
        info.update(r_info)
    shift = d_sigma + d_res + d_geom
    info.update(bulk_qsgw_shift_ev=float(d_qsgw), bulk_vertex_fraction=float(frac),
                bulk_vertex_correction_ev=float(d_sigma - d_qsgw), bulk_selfenergy_shift_ev=float(d_sigma),
                bulk_geometry=BULK_GEOMETRY, bulk_geometry_shift_ev=float(d_geom), bulk_strain_fraction=float(s),
                bulk_pbe_gap_dot_lattice_ev=float(pbe_gap_dot), bulk_residual=BULK_RESIDUAL,
                bulk_residual_shift_ev=float(d_res))
    if eps_eff is not None:
        info["bulk_vertex_eps_eff"] = float(eps_eff)
    if BULK_VERTEX != "none":
        extra = f", Penn eps_eff {eps_eff:.2f} at the DFT gap" if eps_eff is not None else ""
        logger.info(f"  [QP] Bulk vertex correction ({BULK_VERTEX}, factor {a_vertex:.3f}"
                    f"{' (material)' if BULK_VERTEX_FACTOR == 'material' else ''}): fraction {frac:.3f}{extra}; "
                    f"Delta_Sigma {d_qsgw:.3f} -> {d_sigma:.3f} eV")
    if BULK_GEOMETRY != "none" and d_geom != 0.0:
        logger.info(f"  [QP] Bulk geometry correction ({BULK_GEOMETRY}, s = {s:.3f}): Delta_geom {d_geom:+.3f} eV "
                    f"(bulk PBE gap {pbe_gap:.3f} eV at a_exp, {pbe_gap_dot:.3f} eV at the dot's lattice)")
    if BULK_RESIDUAL != "none":
        logger.info(f"  [QP] Bulk residual (experimental): {d_res:+.3f} eV (E_exp {r_info.get('bulk_gap_exp_ev', float('nan')):.3f} eV)")
    logger.info(f"  [QP] Bulk shift = Delta_Sigma {d_sigma:.3f} + residual {d_res:+.3f} + Delta_geom {d_geom:+.3f} "
                f"= {shift:.3f} eV")
    return float(shift), info


# g-xTB route (qdex.xtb): bulk QP shift for g-xTB orbitals on g-xTB-relaxed geometries.
#
#   Delta_bulk^gxtb = E_g^ref,SF - E_g^gxtb,bulk
#   E_g^ref,SF      = gap_exp + [E_g^PBE(a_exp) - E_g^PBE+SOC(a_exp)]   (spin-free experimental gap)
#   E_g^gxtb,bulk  ~= E_g^PBE(a_PBE) + <E_g^gxtb(dot_i; g-xTB geometry) - E_g^PBE(dot_i; PBE geometry)>_i
#
# The spin-free reference adds back the spin-orbit lowering of the gap of PBE with the QDEX GTH-SOC
# operator (MATERIAL_DB_BULK_PBE, 0.121 eV for CdSe), the one the SOC spinors of the dots apply, so
# that a SOC run lands on gap_exp in the bulk (as the PBE route does with Delta_Sigma). The PBE dots of
# delta were relaxed with PBE, so their bulk partner is the PBE gap at the PBE lattice (0.474 eV for
# CdSe, not 0.644 at a_exp): E_gxtb,bulk is then the g-xTB bulk gap at the g-xTB lattice, and the shift
# carries the dot from its own lattice to the experimental one, like Delta_geom of the PBE route.
# Equivalently Delta_bulk^gxtb = (E_g^ref,SF - E_g^PBE(a_PBE)) - delta: on average over the dots the
# g-xTB route reproduces the PBE route with a bulk limit at the experimental gap (bulk_vertex_factor
# = (E_exp - E_PBE+SOC)/Delta_Sigma). PROVISIONAL (2026-10-02): periodic g-xTB (xtb-bleed, macOS arm64)
# did not converge reliably, so delta comes from the CdSe dot series 1.2-3.4 nm (CP2K PBE gaps; g-xTB
# gaps with the Cd-Se bonds rescaled to the g-xTB value 2.60 A). Replace it by a periodic g-xTB bulk
# gap when one is available (key 'gap_gxtb_bulk').
GXTB_BULK = {
    # delta (eV), sample std over the dots (eV), number of dots; Delta_so (eV, experimental): used only
    # when MATERIAL_DB_BULK_PBE has no PBE+SOC gap for the material
    "CDSE": dict(delta=3.731, delta_std=0.520, n_dots=6, delta_so=0.42),
}


def gxtb_bulk_shift(material_name):
    """Bulk QP shift for g-xTB orbitals (negative: g-xTB overestimates the gap). Returns (shift_ev, info)."""
    m_name = str(material_name).upper() if material_name else "DEFAULT"
    data, entry = GXTB_BULK.get(m_name), MATERIAL_DB.get(m_name)
    if data is None or entry is None:
        raise ValueError(f"quasiparticles.reference: gxtb has no g-xTB bulk data for material '{material_name}'. "
                         f"Available: {', '.join(sorted(GXTB_BULK))}")
    bulk = MATERIAL_DB_BULK_PBE.get(m_name, {})
    gap_exp = float(entry[3])
    if bulk.get("gap_soc_exp_lattice") is not None:
        soc_lowering = float(bulk["gap_exp_lattice"] - bulk["gap_soc_exp_lattice"])
        soc_source = "PBE/GTH-SOC at a_exp"
    else:
        soc_lowering, soc_source = data["delta_so"] / 3.0, "experimental Delta_so/3"
    gap_ref_sf = gap_exp + soc_lowering
    gap_pbe_bulk = float(bulk.get("gap_pbe_lattice", entry[7]))   # partner of the PBE-relaxed dots of delta
    if "gap_gxtb_bulk" in data:
        gap_gxtb_bulk, source = float(data["gap_gxtb_bulk"]), "periodic g-xTB"
    else:
        gap_gxtb_bulk, source = gap_pbe_bulk + data["delta"], "PBE bulk at a_PBE + dot-averaged g-xTB/PBE offset"
    shift = gap_ref_sf - gap_gxtb_bulk
    info = dict(qp_reference="gxtb", gap_ref_spin_free_ev=gap_ref_sf, gap_exp_ev=gap_exp,
                soc_lowering_ev=soc_lowering, soc_lowering_source=soc_source, gap_pbe_bulk_ev=gap_pbe_bulk,
                gap_gxtb_bulk_ev=gap_gxtb_bulk, gxtb_bulk_source=source, gxtb_pbe_offset_ev=data.get("delta"),
                gxtb_pbe_offset_std_ev=data.get("delta_std"), bulk_shift_ev=shift)
    logger.info(f"  [QP] g-xTB bulk shift: E_ref(SF) {gap_exp:.3f} + {soc_lowering:.3f} ({soc_source}) - "
                f"E_gxtb(bulk) {gap_gxtb_bulk:.3f} = {shift:+.3f} eV ({source})")
    return float(shift), info


def _mnok_denom(r_mat_au, damp_mat_au, exponent=None):
    """MNOK denominator (r^beta + damp^beta)^(1/beta). Default beta = MNOK_EXPONENT (2, Ohno-Klopman)."""
    beta = float(MNOK_EXPONENT if exponent is None else exponent)
    if beta == 2.0:
        return np.sqrt(r_mat_au**2 + damp_mat_au**2)
    elif beta == 1.0:
        return r_mat_au + damp_mat_au
    return np.power(r_mat_au**beta + damp_mat_au**beta, 1.0 / beta)


# =====================================================================
# sTDA (Grimme, J. Chem. Phys. 138, 244104 (2013)) interaction matrices
# =====================================================================
# Fock-exchange fraction a_x of common functionals (the functional of the MO file).
STDA_FUNCTIONAL_AX = {
    "pbe": 0.0, "blyp": 0.0, "bp86": 0.0, "lda": 0.0, "tpss": 0.0, "r2scan": 0.0,
    "tpssh": 0.10, "b3lyp": 0.20, "b3pw91": 0.20, "pbe0": 0.25, "pw6b95": 0.28,
    "m06": 0.27, "bhlyp": 0.50, "bhandhlyp": 0.50, "m06-2x": 0.54, "hf": 1.0,
}
# Global parameters of sTDA (std2 source, stda.f): beta = b1 + b2 a_x, alpha = a1 + a2 a_x.
STDA_BETA1, STDA_BETA2 = 0.20, 1.83
STDA_ALPHA1, STDA_ALPHA2 = 1.42, 0.48
# Range-separated hybrids: fitted (a_x, alpha, beta) of Risthaus, Hansen, Grimme, PCCP 16, 14408 (2014),
# as in the stda source (main.f: -CAMB3LYP, -wB97XD2, -wB97XD3, -wB97MV). The XsTD range-separation
# terms of that program are not part of plain sTDA. g-xTB is fitted to omegaB97M-V and uses its set.
STDA_RSH_PARAMS = {
    "cam-b3lyp": (0.38, 0.90, 1.86),
    "wb97x-d2": (0.51, 4.51, 8.0),
    "wb97x-d": (0.51, 4.51, 8.0),
    "wb97x-d3": (0.51, 4.51, 8.0),
    "wb97m-v": (0.51, 4.51, 8.0),
    "gxtb": (0.51, 4.51, 8.0),
    "g-xtb": (0.51, 4.51, 8.0),
}
# sTDA-xTB set of the stda program (main.f, -xtb): a_x 0.50, alpha 2.0, beta 4.0 (Grimme, Bannwarth,
# J. Chem. Phys. 145, 054103 (2016)). It was fitted for the sTDA-xTB Hamiltonian (xtb4stda orbitals,
# with a +3.1 eV shift of the virtual levels and a K_ia diagonal shift, neither applied here). No sTDA
# set has been published for GFN2-xTB orbitals; 'gfn2' takes this one as the closest xTB values.
STDA_XTB_PARAMS = {
    "stda-xtb": (0.50, 2.0, 4.0),
    "gfn2": (0.50, 2.0, 4.0),
    "gfn2-xtb": (0.50, 2.0, 4.0),
}


def stda_parameters(functional=None, ax=None, alpha=None, beta=None, material_name=None):
    """
    sTDA parameters (a_x, alpha_K, beta_J, source).

    A fitted set (STDA_RSH_PARAMS, STDA_XTB_PARAMS) gives all three; otherwise a_x comes from
    ``ax`` or the functional (``stda_ax``) and alpha, beta from the global-hybrid formulas.
    Explicit ``ax``, ``alpha`` or ``beta`` override the preset.
    """
    key = str(functional).lower() if functional else None
    if key in STDA_RSH_PARAMS or key in STDA_XTB_PARAMS:
        if key in STDA_RSH_PARAMS:
            ax0, alpha0, beta0 = STDA_RSH_PARAMS[key]
            source = f"{functional} (range-separated set)"
        else:
            ax0, alpha0, beta0 = STDA_XTB_PARAMS[key]
            source = f"{functional} (sTDA-xTB set)"
            if key != "stda-xtb":
                logger.warning(f"  [sTDA] '{functional}': no sTDA parameters exist for GFN2-xTB orbitals; using the "
                               "sTDA-xTB set (a_x 0.50, alpha 2.0, beta 4.0), fitted for xtb4stda orbitals.")
        if ax is not None and str(ax).strip() != "":
            ax0, _ = stda_ax(None, ax, material_name)
            source += ", explicit a_x"
    else:
        ax0, source = stda_ax(functional, ax, material_name)
        alpha0 = STDA_ALPHA1 + STDA_ALPHA2 * ax0
        beta0 = STDA_BETA1 + STDA_BETA2 * ax0
    if alpha is not None:
        alpha0, source = float(alpha), source + ", explicit alpha"
    if beta is not None:
        beta0, source = float(beta), source + ", explicit beta"
    return float(ax0), float(alpha0), float(beta0), source


def stda_ax(functional=None, ax=None, material_name=None):
    """Resolve a_x from an explicit value, 'dielectric' (1/eps_inf of the material) or a functional name."""
    if ax is not None and str(ax).strip() != "":
        if str(ax).lower() == "dielectric":
            entry = MATERIAL_DB.get(str(material_name).upper()) if material_name else None
            if entry is None:
                raise ValueError("ax: dielectric needs a material with eps_inf in MATERIAL_DB")
            return 1.0 / float(entry[0]), f"1/eps_inf ({material_name})"
        return float(ax), "explicit"
    if functional is None:
        raise ValueError("sTDA needs excitations.functional (functional of the MO file) or excitations.ax")
    key = str(functional).lower()
    if key not in STDA_FUNCTIONAL_AX:
        raise ValueError(f"Unknown functional '{functional}' for sTDA; give excitations.ax explicitly. "
                         f"Known: {', '.join(sorted(STDA_FUNCTIONAL_AX))}")
    return STDA_FUNCTIONAL_AX[key], str(functional)


def build_stda_gammas(atom_symbols, coords, ax, eta_dict=HARDNESS_DICT, alpha=None, beta=None):
    """Coulomb (J) and exchange (K) interaction matrices of sTDA, in eV.

    gamma^J_AB = (R^beta  + (a_x eta_AB)^-beta)^(-1/beta),   beta  = 0.20 + 1.83 a_x
    gamma^K_AB = (R^alpha + eta_AB^-alpha)^(-1/alpha),        alpha = 1.42 + 0.48 a_x
    with R in bohr and eta_AB = (eta_A + eta_B)/2 in hartree.  Grimme's hardness is twice the
    Ghosh-Islam value stored in HARDNESS_DICT ((ii|ii) = IP - EA).  gamma^J enters the direct term
    without a further a_x prefactor (std2, rtdamat); for a_x = 0 it vanishes.
    ``alpha`` and ``beta`` replace the global-hybrid formulas (range-separated sets, stda_parameters).
    """
    ax = float(ax)
    beta = STDA_BETA1 + STDA_BETA2 * ax if beta is None else float(beta)
    alpha = STDA_ALPHA1 + STDA_ALPHA2 * ax if alpha is None else float(alpha)
    r = squareform(pdist(np.asarray(coords, dtype=float) / ANG_PER_BOHR))
    eta = 2.0 * np.array([eta_dict.get(s.lower(), 5.0) for s in atom_symbols]) / HA_TO_EV
    eta_ab = 0.5 * (eta[:, None] + eta[None, :])
    gam_k = (r ** alpha + eta_ab ** (-alpha)) ** (-1.0 / alpha)
    if ax > 0.0:
        gam_j = (r ** beta + (ax * eta_ab) ** (-beta)) ** (-1.0 / beta)
    else:
        gam_j = np.zeros_like(gam_k)
    return gam_j * HA_TO_EV, gam_k * HA_TO_EV, {"ax": ax, "alpha_K": alpha, "beta_J": beta}


# =====================================================================
# Model 1: Classic MNOK Kernel (sTDA style)
# =====================================================================
def build_gamma(atom_symbols, coords, alpha, beta=0.0, eta_dict=HARDNESS_DICT, exponent=None):
    """
    Standard Ohno-Klopman kernel. 
    `alpha` scales the entire matrix (macroscopic screening, e.g., 1/eps_inf).
    `beta` [0.0 to 1.0] applies exact-exchange Hubbard stiffening to the diagonal.
    """
    BOHR_TO_ANG = ANG_PER_BOHR
    HA_TO_EV_val = HA_TO_EV
    coords_au = coords / BOHR_TO_ANG
    r_mat_au = squareform(pdist(coords_au))
    
    etas_ev = np.array([eta_dict.get(s.lower(), 5.0) for s in atom_symbols])
    etas_au = etas_ev / HA_TO_EV_val
    a_au = 1.0 / (MNOK_ONSITE_SCALE * etas_au)
    damp_mat_au = 0.5 * (a_au[:, np.newaxis] + a_au[np.newaxis, :])
    
    gamma_au = 1.0 / _mnok_denom(r_mat_au, damp_mat_au, exponent)
    gamma_ev = gamma_au * HA_TO_EV
    
    # Apply macroscopic screening
    gamma_mat = alpha * gamma_ev
    
    # Apply short-range beta stiffening to the diagonal
    if beta > 0.0:
        raise ValueError(
            "beta > 0 is disabled: the bare on-site U table has not been defined and validated"
        )
            
    return gamma_mat

# =====================================================================
# Model 2: Universal Tunable Resta-MNOK Kernel (Dual-Screening)
# =====================================================================
def build_resta_mnok(atom_symbols, coords, alpha, material_name, eps_out=2.0, eta_dict=HARDNESS_DICT):

    BOHR_TO_ANG = ANG_PER_BOHR
    HA_TO_EV_val = HA_TO_EV

    # --------------------------------------------------
    # 1. Geometry and Size Metrics
    # --------------------------------------------------
    metrics = get_cluster_size_metrics(coords, atom_symbols, material_name)
    R_QD_ang = metrics['R_eff_hull']

    m_name = material_name.upper() if material_name else "DEFAULT"
    entry = MATERIAL_DB.get(m_name, MATERIAL_DB["DEFAULT"])

    eps_inf_bulk = entry[0]   # electronic dielectric

    # --------------------------------------------------
    # 2. Core atoms for geometry
    # --------------------------------------------------
    core_coords = coords
    if m_name in MATERIAL_ELEMENTS:
        core_elements = [el.lower() for el in MATERIAL_ELEMENTS[m_name]]
        core_coords = np.array(
            [coords[i] for i, sym in enumerate(atom_symbols) if sym.lower() in core_elements]
        )

    if len(core_coords) > 1:
        r_core_ang = squareform(pdist(core_coords))
        np.fill_diagonal(r_core_ang, np.inf)
        d_NN_ang = np.median(np.min(r_core_ang, axis=1))
    else:
        d_NN_ang = 2.5

    # --------------------------------------------------
    # 3. Screening length (Thomas-Fermi / Resta Model)
    # --------------------------------------------------
    d_NN_au = d_NN_ang / BOHR_TO_ANG

    # Tie screening length to bulk dielectric response (Microscopic)
    k_s_au = np.sqrt(max(0.0, eps_inf_bulk - 1.0)) / d_NN_au
    
    # Calculate the physical screening length in Angstroms
    lambda_s_ang = (1.0 / k_s_au) * BOHR_TO_ANG if k_s_au > 0 else 0.0

    # --------------------------------------------------
    confinement_ratio = R_QD_ang / lambda_s_ang if lambda_s_ang > 0 else np.inf

    logger.info("\n    [Kernel: Electronic Resta-MNOK]")
    logger.info(f"    Material            = {m_name}")
    logger.info(f"    R_QD (hull_eff)     = {R_QD_ang:.3f} Å")
    logger.info(f"    R / lambda_s        = {confinement_ratio:.3f}")
    logger.info(f"    epsilon_in           = {eps_inf_bulk:.3f} (electronic/high-frequency)")
    logger.info(f"    epsilon_out          = {eps_out:.3f} (QP model only; not used in RESTA)")
    logger.info(f"    Screening length    = {lambda_s_ang/BOHR_TO_ANG:.3f} a.u. ({lambda_s_ang:.3f} Å)")
    logger.info("")

    # --------------------------------------------------
    # 5. Build MNOK matrices
    # --------------------------------------------------
    coords_au = coords / BOHR_TO_ANG
    r_mat_au = squareform(pdist(coords_au))

    etas_au = np.array([eta_dict[s.lower()] for s in atom_symbols]) / HA_TO_EV
    a_au = 1.0 / (MNOK_ONSITE_SCALE * etas_au)

    damp_mat_au = 0.5 * (a_au[:, np.newaxis] + a_au[np.newaxis, :])
    mnok_denom_au = _mnok_denom(r_mat_au, damp_mat_au)

    # --------------------------------------------------
    # 6. Electronic screened direct kernel.  The environment is intentionally
    # absent here; eps_out is handled only by the QP polarization model.
    # --------------------------------------------------
    c_inf = 1.0 / eps_inf_bulk
    w_resta_au = (c_inf + (1.0 - c_inf) * np.exp(-k_s_au * r_mat_au)) / mnok_denom_au
    w_resta_ev = w_resta_au * HA_TO_EV
    # Preserve the historical two-return-value API.  Both are electronic W;
    # the first is used only by the explicitly experimental COHSEX path.
    return w_resta_ev, w_resta_ev


def build_sphere_reaction_field(coords, atom_symbols, material_name, eps_out, eps_in=None,
                                n_terms=300, x_max=0.9, radius_definition=None):
    """Atom-pair reaction field of a dielectric sphere (eV).

    G(r, r') = (e^2/R) sum_n (eps_in - eps_out)(n + 1) / [eps_in (n eps_in + (n + 1) eps_out)]
               (r r'/R^2)^n P_n(cos theta)

    is the potential at r of the surface charge induced by a unit charge at r'
    (Boettcher; Brus 1984).  Its 1S-averaged diagonal, G(r, r)/2 per carrier,
    is the polarization term P(R) of the anchor-scaled QP model; adding G to
    the static BSE kernel gives the electron-hole image interaction of the same
    sphere, so the surface polarization cancels in the neutral excitation as
    it does in the QP gap.  The sphere is centred on the core atoms with
    R = qd_radius (SAXS by default) and eps_in = eps_inf by default.  Radial fractions r/R are
    capped at ``x_max`` because the image series diverges at the boundary for
    atoms (ligands) at or beyond R; the cap is a regularization choice.
    """
    m_name = str(material_name).upper() if material_name else "DEFAULT"
    entry = MATERIAL_DB.get(m_name, MATERIAL_DB["DEFAULT"])
    ei = float(entry[0]) if eps_in is None else float(eps_in)
    eo = max(1.0, float(eps_out))
    coords = np.asarray(coords, dtype=float)
    n_at = len(coords)
    if abs(ei - eo) < 1.0e-12:
        return np.zeros((n_at, n_at))
    metrics = get_cluster_size_metrics(coords, atom_symbols, m_name)
    R = qd_radius(coords, atom_symbols, m_name, radius_definition)
    sel = metrics.get("selected_atom_indices") or list(range(n_at))
    center = np.mean(coords[np.asarray(sel, dtype=int)], axis=0)
    rel = coords - center
    r = np.linalg.norm(rel, axis=1)
    x = np.minimum(r / R, float(x_max))
    unit = rel / np.maximum(r, 1.0e-12)[:, None]
    cos_t = np.clip(unit @ unit.T, -1.0, 1.0)
    cos_t[r < 1.0e-12, :] = 1.0
    cos_t[:, r < 1.0e-12] = 1.0
    xx = np.outer(x, x)
    G = np.zeros((n_at, n_at))
    P_prev = np.ones_like(cos_t)
    P_curr = cos_t.copy()
    pow_n = np.ones_like(xx)
    xmax2 = float(x.max()) ** 2
    tmp = np.empty_like(xx)
    for n in range(int(n_terms)):
        if n == 0:
            Pn = P_prev
        elif n == 1:
            Pn = P_curr
        else:
            # P_n = ((2n-1) x P_{n-1} - (n-1) P_{n-2}) / n, in place
            np.multiply(cos_t, P_curr, out=tmp)
            tmp *= (2 * n - 1) / n
            P_prev *= (n - 1) / n
            tmp -= P_prev
            P_prev, P_curr, tmp = P_curr, tmp, P_prev
            Pn = P_curr
        cn = (ei - eo) * (n + 1.0) / (ei * (n * ei + (n + 1.0) * eo))
        G += cn * pow_n * Pn
        # |P_n| <= 1 and (r r'/R^2)^n <= xmax^(2n): stop once the tail is negligible
        if n > 2 and abs(cn) * xmax2 ** (n + 1) / max(1.0e-300, 1.0 - xmax2) < 1.0e-10 * abs(G[0, 0] if G[0, 0] else 1.0):
            break
        pow_n *= xx
    return G * COULOMB_EV_ANG / R


# =========================================================================
# ATOMISTIC POLARIZABLE DIPOLE INTERACTION MODEL (DIM / THOLE MODEL)
# =========================================================================

POLARIZABILITY_TABLE_AU = {
    # Atomic polarizabilities in Bohr^3 (from standard CRC/Miller/Thole tables)
    "H": 4.5, "HE": 1.4,
    "LI": 164.0, "BE": 38.0, "B": 21.0, "C": 11.8, "N": 7.4, "O": 5.4, "F": 3.8, "NE": 2.7,
    "NA": 163.0, "MG": 71.0, "AL": 58.0, "SI": 37.3, "P": 25.0, "S": 19.6, "CL": 15.0, "AR": 11.1,
    "K": 290.0, "CA": 160.0, "SC": 97.0, "TI": 80.0, "V": 64.0, "CR": 52.0, "MN": 46.0, "FE": 40.0,
    "CO": 35.0, "NI": 32.0, "CU": 42.0, "ZN": 38.0, "GA": 49.0, "GE": 41.0, "AS": 29.0, "SE": 31.0,
    "BR": 21.0, "KR": 16.8,
    "RB": 319.0, "SR": 197.0, "Y": 120.0, "ZR": 95.0, "NB": 80.0, "MO": 68.0, "TC": 58.0, "RU": 50.0,
    "RH": 43.0, "PD": 38.0, "AG": 55.0, "CD": 48.0, "IN": 65.0, "SN": 53.0, "SB": 43.0, "TE": 38.0,
    "I": 35.0, "XE": 27.3,
    "CS": 400.0, "BA": 275.0, "LA": 150.0, "HF": 90.0, "TA": 75.0, "W": 64.0, "RE": 54.0, "OS": 48.0,
    "IR": 42.0, "PT": 40.0, "AU": 40.0, "HG": 34.0, "TL": 50.0, "PB": 47.0, "BI": 45.0,
    "DEFAULT": 25.0
}


def build_dim_screening_factors(coords, atom_symbols, material_name=None, eps_out=2.4, alpha=1.0):
    """
    Computes real-space atomistic dielectric screening factors S_AB via the
    Atomistic Polarizable Dipole Interaction Model (DIM / Thole Model).

    References:
      - J. Applequist, J. R. Carl, K.-K. Fung, J. Am. Chem. Soc. 94, 2956 (1972).
      - B. T. Thole, Chem. Phys. 59, 341 (1981).
      - M. Lannoo, C. Delerue, G. Allan, Phys. Rev. Lett. 74, 3415 (1995).
      - C. Delerue, M. Lannoo, G. Allan, Phys. Rev. B 53, 15837 (1996).
    """
    n_atoms = len(atom_symbols)
    m_name = material_name.upper() if material_name else "DEFAULT"
    entry = MATERIAL_DB.get(m_name, MATERIAL_DB["DEFAULT"])
    eps_inf_bulk = float(entry[0])

    coords_bohr = coords / ANG_PER_BOHR

    # 1. Atomic polarizabilities in Bohr^3
    pol_raw = np.array([POLARIZABILITY_TABLE_AU.get(s.upper(), POLARIZABILITY_TABLE_AU["DEFAULT"]) for s in atom_symbols], dtype=np.float64)
    pol = pol_raw * float(alpha)

    # 2. Assemble 3N x 3N Thole dipole interaction matrix M = alpha^-1 + T
    M = np.zeros((3 * n_atoms, 3 * n_atoms), dtype=np.float64)
    for i in range(n_atoms):
        M[3 * i : 3 * i + 3, 3 * i : 3 * i + 3] = np.eye(3) / max(1e-4, pol[i])

    a_thole = 2.1304
    # Thole-damped dipole tensor for all pairs at once (vectorized; same formula as before)
    rij = coords_bohr[:, None, :] - coords_bohr[None, :, :]
    R = np.linalg.norm(rij, axis=-1)
    off = ~np.eye(n_atoms, dtype=bool) & (R >= 1e-6)
    Rs = np.where(off, R, 1.0)
    u = Rs / np.maximum(1e-4, np.outer(pol, pol) ** (1.0 / 6.0))
    au = a_thole * u
    exp_au = np.exp(-au)
    fT = 1.0 - (1.0 + au + 0.5 * au ** 2) * exp_au
    fD = 1.0 - (1.0 + au + 0.5 * au ** 2 + (1.0 / 6.0) * au ** 3) * exp_au
    T = (fT / Rs ** 3)[:, :, None, None] * np.eye(3)[None, None, :, :] \
        - (3.0 * fD / Rs ** 5)[:, :, None, None] * rij[:, :, :, None] * rij[:, :, None, :]
    T[~off] = 0.0
    M += T.transpose(0, 2, 1, 3).reshape(3 * n_atoms, 3 * n_atoms)
    del T, rij

    try:
        invM = np.linalg.inv(M)
    except np.linalg.LinAlgError:
        invM = np.linalg.pinv(M, rcond=1e-8)

    # 3. Determine atom-specific effective polarizability under uniform electric fields (x, y, z)
    E_ext = np.zeros((3 * n_atoms, 3), dtype=np.float64)
    for dim in range(3):
        E_ext[dim::3, dim] = 1.0
    p_resp = invM @ E_ext
    p_eff = sum(p_resp[dim::3, dim] for dim in range(3)) / 3.0

    p_max = np.max(p_eff) if np.max(p_eff) > 0 else 1.0
    eta_atom = np.clip(p_eff / p_max, 0.05, 1.0)

    # 4. Nearest-neighbor distance for screening length
    core_coords = coords
    if m_name in MATERIAL_ELEMENTS:
        core_elements = [el.lower() for el in MATERIAL_ELEMENTS[m_name]]
        core_coords = np.array(
            [coords[i] for i, sym in enumerate(atom_symbols) if sym.lower() in core_elements]
        )

    if len(core_coords) > 1:
        r_core_ang = squareform(pdist(core_coords))
        np.fill_diagonal(r_core_ang, np.inf)
        d_NN_ang = np.median(np.min(r_core_ang, axis=1))
    else:
        d_NN_ang = 2.5
    d_NN_au = d_NN_ang / ANG_PER_BOHR

    # 5. Distance matrix in Angstrom
    R_mat_ang = squareform(pdist(coords))

    # 6. Pairwise screening factor S_AB:
    # Pure electronic internal screening (eps_out is intentionally absent;
    # solvent polarization is handled exclusively by the QP model).
    # Short-range: S -> 1.0 (unscreened atomic core).
    # Long-range core: S -> 1/eps_inf_bulk.
    eta_pair = np.sqrt(np.outer(eta_atom, eta_atom))
    eps_pair = 1.0 + (eps_inf_bulk - 1.0) * eta_pair
    k_s_ang = np.sqrt(np.maximum(0.0, eps_pair - 1.0)) / max(1e-4, d_NN_au) / ANG_PER_BOHR
    c_inf = 1.0 / np.maximum(1.0, eps_pair)
    S_atom = c_inf + (1.0 - c_inf) * np.exp(-k_s_ang * R_mat_ang)
    np.fill_diagonal(S_atom, 1.0)

    return S_atom, eta_atom, eps_inf_bulk, d_NN_ang


def build_dim_mnok(atom_symbols, coords, material_name=None, alpha=1.0, eta_dict=HARDNESS_DICT, **kwargs):
    """
    Constructs the MNOK two-electron interaction matrix screened by the
    Atomistic Polarizable Dipole Interaction Model (DIM / Thole Model).
    """
    if eta_dict is None:
        eta_dict = HARDNESS_DICT

    S_atom, eta_atom, eps_inf_bulk, d_NN_ang = build_dim_screening_factors(
        coords=coords, atom_symbols=atom_symbols, material_name=material_name, alpha=alpha
    )

    BOHR_TO_ANG = ANG_PER_BOHR
    coords_au = coords / BOHR_TO_ANG
    r_mat_au = squareform(pdist(coords_au))

    etas_au = np.array([eta_dict[s.lower()] for s in atom_symbols]) / HA_TO_EV
    a_au = 1.0 / (MNOK_ONSITE_SCALE * etas_au)

    damp_mat_au = 0.5 * (a_au[:, np.newaxis] + a_au[np.newaxis, :])
    mnok_denom_au = _mnok_denom(r_mat_au, damp_mat_au)
    gamma_mnok_bare_ev = (1.0 / mnok_denom_au) * HA_TO_EV

    w_dim_ev = S_atom * gamma_mnok_bare_ev

    n_atoms = len(atom_symbols)
    m_name = material_name.upper() if material_name else "DEFAULT"
    inter_mask = ~np.eye(n_atoms, dtype=bool)
    eps_eff_median = float(1.0 / np.median(S_atom[inter_mask])) if np.any(inter_mask) else 1.0

    logger.info(f"\n    [Kernel: Atomistic Polarizable Dipole Model (DIM-MNOK)]")
    logger.info(f"    Material            = {m_name}")
    logger.info(f"    epsilon_in (bulk)   = {eps_inf_bulk:.3f}")
    logger.info(f"    Nearest neighbor    = {d_NN_ang:.3f} Å")
    logger.info(f"    Median interatomic ε= {eps_eff_median:.3f}")
    logger.info(f"    Dipole Matrix (3N)  = {3*n_atoms} x {3*n_atoms}")

    eps_info = {
        "eps_eff_exciton": eps_eff_median,
        "eps_bulk": eps_inf_bulk,
        "eps_interatomic": eps_eff_median,
        "kernel_mode": "dim"
    }
    return w_dim_ev, w_dim_ev, gamma_mnok_bare_ev, eps_info


def build_sbse_kernel(atom_symbols, coords, atom_ao_ranges=None, shells=None,
                      C_occ_low=None, C_virt_low=None, eps_occ=None, eps_virt=None,
                      C_occ_b_low=None, C_virt_b_low=None, eps_occ_b=None, eps_virt_b=None,
                      mode="atom", eps_out=1.0, material_name=None, alpha=1.0,
                      nthreads=1, eta_dict=HARDNESS_DICT, return_eps_info=False):
    """
    Constructs the simplified Bethe-Salpeter Equation (sBSE) screened interaction kernel W
    following Cho, Bintrim, and Berkelbach [J. Chem. Theory Comput. 18, 3438 (2022), doi:10.1021/acs.jctc.2c00087]::

        W = (I + J_solv * Pi^0)^{-1} * J_solv

    Features:
      - mode in ['atom', 'sbse-atom']: Atom-resolved sBSE (matrix dimension N_atom x N_atom)
      - mode in ['ao', 'sbse-ao']: AO-resolved sBSE (matrix dimension N_ao x N_ao)
      - Solvent screening: if eps_out > 1.0, J_solv includes the asymptotic dielectric screening
        J_solv = J_bare - (1 - 1/eps_out) / sqrt(R_{AB}^2 + (2 R_{QD})^2)
      - Parameter-free: metric S' cancels identically; screening emerges purely from transition
        polarizability Pi^0 and Coulomb repulsion J.
    """
    from qdex.constants import HA_TO_EV, ANG_PER_BOHR
    from scipy.spatial.distance import pdist, squareform

    if C_occ_low is None or C_virt_low is None or eps_occ is None or eps_virt is None:
        raise ValueError(
            "sBSE kernel requires active Löwdin MOs (C_occ_low, C_virt_low) and eigenvalues (eps_occ, eps_virt)."
        )

    n_atoms = len(atom_symbols)
    m_name = material_name.upper() if material_name else "DEFAULT"
    entry = MATERIAL_DB.get(m_name, MATERIAL_DB["DEFAULT"])
    eps_bulk = float(entry[0])

    # 1. Cluster size metrics for solvent asymptotic boundary
    try:
        size_metrics = get_cluster_size_metrics(coords, atom_symbols, m_name)
        r_qd_ang = float(size_metrics.get("R_eff_hull", 0.0))
    except Exception:
        r_qd_ang = 0.0

    if r_qd_ang <= 0.0:
        r_qd_ang = float(0.5 * np.max(pdist(coords))) if len(coords) > 1 else 5.0
    r_qd_au = r_qd_ang / ANG_PER_BOHR

    coords_bohr = coords / ANG_PER_BOHR
    r_mat_au = squareform(pdist(coords_bohr))
    etas_au = np.array([eta_dict.get(s.lower(), 7.0) for s in atom_symbols]) / HA_TO_EV
    a_au = 1.0 / (MNOK_ONSITE_SCALE * etas_au)
    damp_mat_au = 0.5 * (a_au[:, np.newaxis] + a_au[np.newaxis, :])
    J_bare_atom_au = 1.0 / _mnok_denom(r_mat_au, damp_mat_au)

    # Solvent screening correction on J
    eps_out_val = float(eps_out) if eps_out is not None else 1.0
    if eps_out_val > 1.0:
        corr_atom_au = (1.0 - 1.0 / eps_out_val) / np.sqrt(r_mat_au**2 + (2.0 * r_qd_au)**2)
        J_solv_atom_au = J_bare_atom_au - corr_atom_au
    else:
        corr_atom_au = np.zeros_like(J_bare_atom_au)
        J_solv_atom_au = J_bare_atom_au.copy()

    # Determine mode: atom-resolved vs AO-resolved
    mode_str = str(mode).lower()
    is_ao_mode = mode_str in ["ao", "sbse-ao", "sbse_ao"]

    n_occ_a = C_occ_low.shape[1]
    n_virt_a = C_virt_low.shape[1]
    d_alpha_au = (eps_virt[:, None] - eps_occ[None, :]) / HA_TO_EV
    inv_sqrt_deps_a = np.sqrt(4.0 / np.maximum(d_alpha_au.T, 1e-6)).reshape(-1)

    if not is_ao_mode:
        # =====================================================================
        # Atom-Resolved sBSE
        # =====================================================================
        if atom_ao_ranges is None:
            raise ValueError("atom_ao_ranges required for atom-resolved sBSE.")

        n_trans_a = n_occ_a * n_virt_a
        V_atom_a = np.zeros((n_atoms, n_trans_a), dtype=np.float64)
        for A, (a0, a1) in enumerate(atom_ao_ranges):
            Q_A = (C_occ_low[a0:a1, :].T @ C_virt_low[a0:a1, :]).reshape(-1)
            V_atom_a[A, :] = Q_A * inv_sqrt_deps_a

        if C_occ_b_low is not None and C_virt_b_low is not None and eps_occ_b is not None and eps_virt_b is not None:
            n_occ_b = C_occ_b_low.shape[1]
            n_virt_b = C_virt_b_low.shape[1]
            d_beta_au = (eps_virt_b[:, None] - eps_occ_b[None, :]) / HA_TO_EV
            inv_sqrt_deps_b = np.sqrt(2.0 / np.maximum(d_beta_au.T, 1e-6)).reshape(-1)
            # Re-scale alpha channel by sqrt(2) instead of 2 for open-shell
            V_atom_a = V_atom_a * (np.sqrt(2.0) / 2.0)
            V_atom_b = np.zeros((n_atoms, n_occ_b * n_virt_b), dtype=np.float64)
            for A, (a0, a1) in enumerate(atom_ao_ranges):
                Q_B = (C_occ_b_low[a0:a1, :].T @ C_virt_b_low[a0:a1, :]).reshape(-1)
                V_atom_b[A, :] = Q_B * inv_sqrt_deps_b
            V_atom = np.hstack([V_atom_a, V_atom_b])
            n_trans_tot = n_trans_a + (n_occ_b * n_virt_b)
        else:
            V_atom = V_atom_a
            n_trans_tot = n_trans_a

        Pi_atom = float(alpha) * (V_atom @ V_atom.T)
        eps_mat = np.eye(n_atoms) + J_solv_atom_au @ Pi_atom
        try:
            W_au = np.linalg.solve(eps_mat, J_solv_atom_au)
        except np.linalg.LinAlgError:
            W_au = np.linalg.pinv(eps_mat) @ J_solv_atom_au
        W_au = 0.5 * (W_au + W_au.T)

        W_ev = W_au * HA_TO_EV
        J_solv_ev = J_solv_atom_au * HA_TO_EV
        J_bare_ev = J_bare_atom_au * HA_TO_EV

        # Diagnostics
        q_h = np.zeros(n_atoms)
        q_l = np.zeros(n_atoms)
        for A, (a0, a1) in enumerate(atom_ao_ranges):
            q_h[A] = np.sum(np.abs(C_occ_low[a0:a1, -1])**2)
            q_l[A] = np.sum(np.abs(C_virt_low[a0:a1, 0])**2)

        v_eh_solv = float(q_h @ J_solv_ev @ q_l)
        v_eh_bare = float(q_h @ J_bare_ev @ q_l)
        w_eh_screened = float(q_h @ W_ev @ q_l)
        eps_exciton = (v_eh_bare / w_eh_screened) if abs(w_eh_screened) > 1e-12 else 1.0

        mask = ~np.eye(n_atoms, dtype=bool)
        eps_inter = float(1.0 / np.median(W_ev[mask] / J_bare_ev[mask])) if np.any(mask) else 1.0
        mat_dim_str = f"{n_atoms} x {n_atoms} (Atom-resolved)"

    else:
        # =====================================================================
        # AO-Resolved sBSE
        # =====================================================================
        n_ao = C_occ_low.shape[0]
        n_trans_a = n_occ_a * n_virt_a
        V_ao_a = np.empty((n_ao, n_trans_a), dtype=np.float64)
        for i in range(n_occ_a):
            c_i = C_occ_low[:, i, None]
            w_i = inv_sqrt_deps_a[i * n_virt_a : (i + 1) * n_virt_a]
            V_ao_a[:, i * n_virt_a : (i + 1) * n_virt_a] = (c_i * C_virt_low) * w_i[None, :]

        if C_occ_b_low is not None and C_virt_b_low is not None and eps_occ_b is not None and eps_virt_b is not None:
            n_occ_b = C_occ_b_low.shape[1]
            n_virt_b = C_virt_b_low.shape[1]
            d_beta_au = (eps_virt_b[:, None] - eps_occ_b[None, :]) / HA_TO_EV
            inv_sqrt_deps_b = np.sqrt(2.0 / np.maximum(d_beta_au.T, 1e-6)).reshape(-1)
            V_ao_a = V_ao_a * (np.sqrt(2.0) / 2.0)
            V_ao_b = np.empty((n_ao, n_occ_b * n_virt_b), dtype=np.float64)
            for i in range(n_occ_b):
                c_ib = C_occ_b_low[:, i, None]
                w_ib = inv_sqrt_deps_b[i * n_virt_b : (i + 1) * n_virt_b]
                V_ao_b[:, i * n_virt_b : (i + 1) * n_virt_b] = (c_ib * C_virt_b_low) * w_ib[None, :]
            V_ao = np.hstack([V_ao_a, V_ao_b])
            n_trans_tot = n_trans_a + (n_occ_b * n_virt_b)
        else:
            V_ao = V_ao_a
            n_trans_tot = n_trans_a

        Pi_ao = float(alpha) * (V_ao @ V_ao.T)

        # Build J_ao: exact one-center if shells provided, MNOK off-diagonal
        from qdex.integrals import compute_two_electron_ao
        if shells is not None:
            J_bare_ao_au = compute_two_electron_ao(shells, nthreads=nthreads)
            # Replace off-diagonal blocks with MNOK
            if atom_ao_ranges is not None:
                for A, (a0, a1) in enumerate(atom_ao_ranges):
                    for B, (b0, b1) in enumerate(atom_ao_ranges):
                        if A != B:
                            J_bare_ao_au[a0:a1, b0:b1] = J_bare_atom_au[A, B]
        else:
            J_bare_ao_au = np.zeros((n_ao, n_ao), dtype=np.float64)
            if atom_ao_ranges is not None:
                for A, (a0, a1) in enumerate(atom_ao_ranges):
                    for B, (b0, b1) in enumerate(atom_ao_ranges):
                        J_bare_ao_au[a0:a1, b0:b1] = J_bare_atom_au[A, B]

        # Solvent correction on AO matrix
        J_solv_ao_au = J_bare_ao_au.copy()
        if eps_out_val > 1.0 and atom_ao_ranges is not None:
            for A, (a0, a1) in enumerate(atom_ao_ranges):
                for B, (b0, b1) in enumerate(atom_ao_ranges):
                    J_solv_ao_au[a0:a1, b0:b1] -= corr_atom_au[A, B]

        eps_mat = np.eye(n_ao) + J_solv_ao_au @ Pi_ao
        try:
            W_au = np.linalg.solve(eps_mat, J_solv_ao_au)
        except np.linalg.LinAlgError:
            W_au = np.linalg.pinv(eps_mat) @ J_solv_ao_au
        W_au = 0.5 * (W_au + W_au.T)

        W_ev = W_au * HA_TO_EV
        J_solv_ev = J_solv_ao_au * HA_TO_EV
        J_bare_ev = J_bare_ao_au * HA_TO_EV

        # Diagnostics in AO basis
        q_h = np.abs(C_occ_low[:, -1])**2
        q_l = np.abs(C_virt_low[:, 0])**2
        v_eh_solv = float(q_h @ J_solv_ev @ q_l)
        v_eh_bare = float(q_h @ J_bare_ev @ q_l)
        w_eh_screened = float(q_h @ W_ev @ q_l)
        eps_exciton = (v_eh_bare / w_eh_screened) if abs(w_eh_screened) > 1e-12 else 1.0

        mask = ~np.eye(n_ao, dtype=bool)
        eps_inter = float(1.0 / np.median(W_ev[mask] / J_bare_ev[mask])) if np.any(mask) else 1.0
        mat_dim_str = f"{n_ao} x {n_ao} (AO-resolved)"

    logger.info(f"\n    ==========================================================================")
    logger.info(f"    [Kernel: sBSE (Simplified Bethe-Salpeter Equation, mode='{mode_str}')]")
    logger.info(f"    ==========================================================================")
    logger.info(f"    Reference Theory        : Cho, Bintrim, Berkelbach [JCTC 18, 3438 (2022)]")
    logger.info(f"    Material                = {m_name}")
    logger.info(f"    epsilon_in (bulk)       = {eps_bulk:.3f}")
    logger.info(f"    epsilon_out (solvent)   = {eps_out_val:.3f}")
    logger.info(f"    Cluster Radius (R_QD)   = {r_qd_ang:.3f} Å")
    logger.info(f"    Active Transitions      = {n_trans_tot}")
    logger.info(f"    Matrix Dimension        = {mat_dim_str}")
    logger.info(f"    RPA Scaling (alpha)     = {float(alpha):.3f}")
    logger.info(f"    --------------------------------------------------------------------------")
    logger.info(f"    Computed Microscopic Dielectric Constants (ε_eff):")
    logger.info(f"      ε_eff (1S Exciton e-h) : {eps_exciton:8.3f}   [Lowest exciton screening]")
    logger.info(f"      ε_eff (Inter-atomic)   : {eps_inter:8.3f}   [Median inter-site screening]")
    logger.info(f"    ==========================================================================\n")

    eps_info = {
        "eps_eff_exciton": eps_exciton,
        "eps_interatomic": eps_inter,
        "eps_bulk": eps_bulk,
        "eps_out": eps_out_val,
        "cluster_radius_ang": r_qd_ang,
        "kernel_mode": "sbse",
        "sbse_mode": mode_str
    }
    if return_eps_info:
        return W_ev, W_ev, J_bare_ev, eps_info
    return W_ev, W_ev, J_bare_ev


def estimate_sgw_qp_gap(coords, atom_symbols, material_name=None, eps_out=1.0,
                        C_occ_low=None, C_virt_low=None, eps_occ=None, eps_virt=None,
                        atom_ao_ranges=None, shells=None, mode="atom", alpha=1.0, Z=0.8,
                        nthreads=1, return_details=False, dynamic_z=False):
    """
    Computes the Quasiparticle (GW) gap shift for a quantum dot via the microscopic
    screening difference (Delta W / Delta COHSEX) approach::

        Delta E_p = Delta E_p^{bulk} + (Z / 2) * <psi_p | W^{QD} - W^{bulk} | psi_p>

    Features:
      - Uses W^{QD} computed directly from the microscopic sGW/sBSE kernel
      - Solves the exchange-correlation challenge: local v^{xc} cancels out against bulk
      - Requires zero empirical parameters or arbitrary alpha_K fudge factors
      - Incorporates solvent dielectric screening via W^{QD}(eps_out)
    """
    from qdex.constants import HA_TO_EV, ANG_PER_BOHR
    m_name = material_name.upper() if material_name else "DEFAULT"
    entry = MATERIAL_DB.get(m_name, None)
    if entry is None or len(entry) < 4:
        raise ValueError(f"Material '{m_name}' not found in database or lacks bulk GW data.")

    eps_bulk = float(entry[0])
    pbe_bulk_gap = bulk_pbe_gap_dot(m_name)    # at the dot's lattice (bulk_geometry)
    gw_bulk_gap = float(entry[8])
    bulk_shift = gw_bulk_gap - pbe_bulk_gap

    # 1. Compute microscopic W^{QD} using sBSE kernel builder
    res = build_sbse_kernel(
        atom_symbols=atom_symbols,
        coords=coords,
        atom_ao_ranges=atom_ao_ranges,
        shells=shells,
        C_occ_low=C_occ_low,
        C_virt_low=C_virt_low,
        eps_occ=eps_occ,
        eps_virt=eps_virt,
        mode=mode,
        eps_out=eps_out,
        material_name=m_name,
        alpha=alpha,
        nthreads=nthreads,
        return_eps_info=True
    )
    W_ev, _, J_bare_ev, eps_info = res

    # 2. Compute Delta W on diagonal sites: W_{AA}^{QD} - W_{AA}^{bulk}
    # In bulk, W_{AA}^{bulk} = J_{AA} / eps_bulk
    n_dim = W_ev.shape[0]
    W_diag_qd = np.diag(W_ev)
    J_diag_bare = np.diag(J_bare_ev)
    W_diag_bulk = J_diag_bare / max(1.0, eps_bulk)
    delta_W_diag = np.maximum(0.0, W_diag_qd - W_diag_bulk)

    # 3. Project Delta W onto frontier HOMO and LUMO states
    n_atoms = len(atom_symbols)
    is_ao_mode = (n_dim != n_atoms)

    if not is_ao_mode and atom_ao_ranges is not None:
        q_h = np.zeros(n_atoms)
        q_l = np.zeros(n_atoms)
        for A, (a0, a1) in enumerate(atom_ao_ranges):
            q_h[A] = np.sum(np.abs(C_occ_low[a0:a1, -1])**2)
            q_l[A] = np.sum(np.abs(C_virt_low[a0:a1, 0])**2)
        sig_h_stat = 0.5 * float(np.sum(q_h * delta_W_diag))
        sig_l_stat = 0.5 * float(np.sum(q_l * delta_W_diag))
    else:
        q_h_ao = np.abs(C_occ_low[:, -1])**2
        q_l_ao = np.abs(C_virt_low[:, 0])**2
        sig_h_stat = 0.5 * float(np.sum(q_h_ao * delta_W_diag))
        sig_l_stat = 0.5 * float(np.sum(q_l_ao * delta_W_diag))

    if dynamic_z:
        eps_z = float(eps_info.get("eps_eff_exciton", eps_bulk)) if isinstance(eps_info, dict) else eps_bulk
        Z_h = compute_dynamic_z(sig_h_stat, None, eps_z, m_name)
        Z_l = compute_dynamic_z(sig_l_stat, None, eps_z, m_name)
    else:
        Z_h = Z_l = float(Z)
    Z_eh = 0.5 * (Z_h + Z_l)
    delta_sigma_h = Z_h * sig_h_stat
    delta_sigma_l = Z_l * sig_l_stat

    confinement_shift = delta_sigma_h + delta_sigma_l
    total_sgw_scissor = bulk_shift + confinement_shift
    if confinement_shift > 1e-8:
        f_h_micro = float(np.clip(delta_sigma_h / confinement_shift, 0.0, 1.0))
        f_l_micro = 1.0 - f_h_micro
    else:
        f_h_micro = 0.5
        f_l_micro = 0.5

    logger.info(f"  [sGW Model] Microscopic Quasiparticle Correction for {m_name}:")
    logger.info(f"    Bulk PBE -> GW Gap   : {pbe_bulk_gap:.3f} -> {gw_bulk_gap:.3f} eV (Shift: +{bulk_shift:.3f} eV)")
    logger.info(f"    Microscopic Delta W  : HOMO shift = +{delta_sigma_h:.3f} eV, LUMO shift = +{delta_sigma_l:.3f} eV")
    logger.info(f"    Microscopic Split    : HOMO takes {f_h_micro*100:.1f}%, LUMO takes {f_l_micro*100:.1f}%")
    logger.info(f"    Confinement Opening  : +{confinement_shift:.3f} eV (Z_h = {Z_h:.3f}, Z_l = {Z_l:.3f}"
          f"{', derived' if dynamic_z else ''})")
    logger.info(f"    Solvent Dielectric   : eps_out = {float(eps_out):.2f}")
    logger.info(f"    ==> Total sGW Scissor: +{total_sgw_scissor:.3f} eV\n")

    provenance = {
        "material": m_name,
        "bulk_pbe_gap_ev": pbe_bulk_gap,
        "bulk_gw_gap_ev": gw_bulk_gap,
        "bulk_shift_ev": bulk_shift,
        "bulk_gw_shift_ev": bulk_shift,
        "confinement_shift_ev": confinement_shift,
        "delta_sigma_homo_ev": delta_sigma_h,
        "delta_sigma_lumo_ev": delta_sigma_l,
        "f_homo": f_h_micro,
        "f_lumo": f_l_micro,
        "f_homo_micro": f_h_micro,
        "f_lumo_micro": f_l_micro,
        "z_factor": float(Z_eh),
        "z_homo": float(Z_h),
        "z_lumo": float(Z_l),
        "dynamic_z": bool(dynamic_z),
        "eps_out": float(eps_out),
        "total_scissor_ev": total_sgw_scissor,
        "total_scissor_solvent_ev": total_sgw_scissor,
        "total_scissor_vacuum_ev": total_sgw_scissor,
        "qp_model": "sgw_dw"
    }
    # The screened interaction this QP model is built on.  The CLI passes it
    # to the BSE (kernel "qp") so that GW and BSE share the same W.
    # Delta W enters the kernel with the same Z as the QP shift (scale_w_difference);
    # the bulk reference is J/eps_inf, as in the site-diagonal QP term.
    provenance["w_bse_ev"] = scale_w_difference(W_ev, J_bare_ev / max(1.0, eps_bulk), Z_eh)
    provenance["w_bse_z"] = float(Z_eh)
    provenance["w_bse_label"] = "sBSE RPA W (solvent inside J)"

    if return_details:
        return total_sgw_scissor, provenance
    return total_sgw_scissor


def compute_dynamic_z(delta_sigma_stat_ev, gap_ev=None, eps_eff=None, material_name=None, omega_p_ev=None):
    """
    Quasiparticle weight of the finite-size correction from one plasmon pole::

        Z_p = [1 + abs(Delta Sigma_p) / omega_tilde]^-1,
        omega_tilde = omega_p / sqrt(1 - 1/eps_eff)

    The finite-size self-energy is modelled as coupling to one bosonic mode of
    frequency omega_tilde with static strength abs(Delta Sigma_p), so
    dSigma/domega = -abs(Delta Sigma_p)/omega_tilde and Z_p follows from the
    linearized QP equation.  omega_tilde is the generalized plasmon-pole
    frequency of Hybertsen and Louie (PRB 34, 5390 (1986)) for the same
    eps_eff that builds the model's W: with a single Lorentz oscillator,
    eps(omega) = 1 + omega_p^2 / (E_P^2 - omega^2) and eps(0) = eps_eff, it is
    the zero of eps(omega), sqrt(E_P^2 + omega_p^2).  It is also the pole of the
    Born reaction field (1/eps_out - 1/eps(omega)).  omega_p is the valence
    plasmon of the bulk material (``valence_plasmon_ev``) unless given.
    ``gap_ev`` is not used; it is kept for call compatibility.

    Limits: one pole only; the bulk part of the dynamical renormalization is
    already inside the tabulated bulk GW gap and is not recomputed here.
    """
    wp = float(omega_p_ev) if omega_p_ev is not None else valence_plasmon_ev(material_name)
    eps_val = max(1.01, float(eps_eff) if eps_eff is not None else 1.01)
    omega_tilde = wp / np.sqrt(1.0 - 1.0 / eps_val)
    sig_stat = abs(float(delta_sigma_stat_ev))
    Z = 1.0 / (1.0 + sig_stat / omega_tilde)
    return float(np.clip(Z, 0.50, 1.0))


def solvent_reaction_term(coords, atom_symbols, material_name, eps_out, eps_bulk, R_mat_ang, R_QD_ang,
                          mode="sphere"):
    """Environment part of Delta W for the Delta-W models (eV, atom pairs).

    mode 'sphere' (default): reaction field of a dielectric sphere with eps_bulk inside and eps_out
    outside, the full multipole Green function (``build_sphere_reaction_field``), the same term the
    two-anchor gw model uses.  mode 'born': the earlier softened Born form
    (1/eps_out - 1/eps_bulk) e^2 / sqrt(r_AB^2 + R^2).
    """
    eps_out_val = max(1.0, float(eps_out))
    if str(mode).lower() == "born":
        return ((1.0 / eps_out_val - 1.0 / eps_bulk) * COULOMB_EV_ANG) / np.sqrt(R_mat_ang ** 2 + R_QD_ang ** 2)
    return build_sphere_reaction_field(np.asarray(coords, dtype=float), atom_symbols, material_name,
                                       eps_out_val, eps_in=eps_bulk)


def expand_atom_to_ao(M_atom, atom_ao_ranges, n_ao):
    """Expand an atom-pair matrix to AO blocks: M_ao[mu, nu] = M[A(mu), B(nu)]."""
    owner = np.empty(n_ao, dtype=int)
    for A, (a0, a1) in enumerate(atom_ao_ranges):
        owner[a0:a1] = A
    M_atom = np.asarray(M_atom, dtype=float)
    return M_atom[np.ix_(owner, owner)]


def ao_screening_ratio(W_atom, gamma_atom, atom_ao_ranges, n_ao):
    """AO expansion of the atom-pair screening ratio W_AB / gamma_AB (dimensionless)."""
    ratio = np.asarray(W_atom, dtype=float) / np.asarray(gamma_atom, dtype=float)
    return expand_atom_to_ao(ratio, atom_ao_ranges, n_ao)


def ao_delta_w(dW_screen_atom, gamma_atom, dW_add_atom, gamma_ao, atom_ao_ranges):
    """Delta W in the AO density-pair (xs) representation.

    The screening part of Delta W is a change of the atom-pair screening ratio,
    Delta S_AB = Delta W_AB / gamma_AB (gamma = MNOK), and is applied to the
    exact AO integrals (mu mu | nu nu).  Additive classical terms (solvent
    reaction field, sphere image) are smooth and are expanded to AO blocks.
    """
    n_ao = gamma_ao.shape[0]
    dW = ao_screening_ratio(dW_screen_atom, gamma_atom, atom_ao_ranges, n_ao) * gamma_ao
    if dW_add_atom is not None:
        dW = dW + expand_atom_to_ao(dW_add_atom, atom_ao_ranges, n_ao)
    return dW


def scale_w_difference(W_qd_ev, W_bulk_ev, Z_eh):
    """Kernel counterpart of the Z-renormalized QP correction.

    Returns W_bulk + Z_eh (W_QD - W_bulk).  The same pole that renormalizes
    the QP shift by Z also screens the electron-hole interaction dynamically;
    to first order in 1/omega_tilde the two effects cancel in the neutral
    excitation (Bechstedt et al., PRL 78, 1528 (1997)).  Scaling Delta W by the
    same Z in the static kernel keeps that cancellation, so the solvent
    (surface-polarization) term still drops out of S1 when Z < 1.
    Z_eh = 1 returns W_QD unchanged.
    """
    return np.asarray(W_bulk_ev, dtype=float) + float(Z_eh) * (
        np.asarray(W_qd_ev, dtype=float) - np.asarray(W_bulk_ev, dtype=float)
    )


# ---------------------------------------------------------------------
# Delta-W quasiparticle estimators (sGW / evGW / qsGW) on the Resta or DIM
# screened interaction.  The two dielectric models differ only in how
# W^QD is built at a given gap (_ScreeningModel); everything else is shared.
# ---------------------------------------------------------------------
class _ScreeningModel:
    """W^QD of the Resta (Penn-scaled) or DIM (atomistic polarizable dipole) model on a cluster."""

    def __init__(self, kind, coords, atom_symbols, material_name=None, eps_out=2.4, solvent_term="sphere",
                 alpha=1.0, penn_scaling=True):
        from qdex.constants import HA_TO_EV, ANG_PER_BOHR
        from scipy.spatial.distance import pdist, squareform

        self.kind = kind
        self.label = {"dim": "DIM", "resta": "Resta"}[kind]
        self.coords = coords
        self.atom_symbols = atom_symbols
        self.alpha = alpha
        self.penn_scaling = penn_scaling
        self.solvent_term = solvent_term
        self.eps_out = eps_out
        self.m_name = m_name = material_name.upper() if material_name else "DEFAULT"
        entry = MATERIAL_DB.get(m_name, None)
        if entry is None or len(entry) < 9:
            raise ValueError(f"Material '{m_name}' not found in database or lacks bulk GW data.")

        self.eps_bulk = eps_bulk = float(entry[0])
        self.pbe_bulk_gap = bulk_pbe_gap_dot(m_name)   # bulk limit of the dot's PBE gap (its lattice)
        self.gw_bulk_gap = float(entry[8])
        self.bulk_shift = self.gw_bulk_gap - self.pbe_bulk_gap
        self.qp_bulk_ref = self.gw_bulk_gap   # bulk QP gap used as the screening reference
        self.bulk_vertex_info = {}
        self.n_atoms = len(atom_symbols)
        self.R_QD_ang = qd_radius(coords, atom_symbols, m_name)

        # 1. Bare Ohno-Klopman interaction matrix (eV)
        r_mat_au = squareform(pdist(coords / ANG_PER_BOHR))
        etas_au = np.array([HARDNESS_DICT.get(s.lower(), 5.0) for s in atom_symbols]) / HA_TO_EV
        a_au = 1.0 / (MNOK_ONSITE_SCALE * etas_au)
        damp_mat_au = 0.5 * (a_au[:, None] + a_au[None, :])
        self.gamma_bare_ev = (1.0 / _mnok_denom(r_mat_au, damp_mat_au)) * HA_TO_EV

        # 2. Bulk reference W^bulk (Resta, eps_inf) on the core nearest-neighbour distance
        core_coords = coords
        if m_name in MATERIAL_ELEMENTS:
            core_elements = [el.lower() for el in MATERIAL_ELEMENTS[m_name]]
            core_coords = np.array([coords[i] for i, sym in enumerate(atom_symbols) if sym.lower() in core_elements])
        if len(core_coords) > 1:
            r_core_ang = squareform(pdist(core_coords))
            np.fill_diagonal(r_core_ang, np.inf)
            d_NN_ang = float(np.median(np.min(r_core_ang, axis=1)))
        else:
            d_NN_ang = 2.5
        self.d_NN_au = d_NN_ang / ANG_PER_BOHR
        self.R_mat_ang = squareform(pdist(coords))
        self.W_bulk_ev = self._resta_w(eps_bulk)

        # 3. Solvent reaction field Delta W^solv
        self.eps_out_val = max(1.0, float(eps_out))
        self.delta_W_solv = solvent_reaction_term(coords, atom_symbols, m_name, self.eps_out_val, eps_bulk,
                                                  self.R_mat_ang, self.R_QD_ang, solvent_term)

    def set_dft_gap(self, dft_gap):
        """Bulk QP shift (with the bulk vertex correction) for a cluster with this DFT gap."""
        if BULK_VERTEX == "none":          # pure QSGW: record the setting, keep the shift as it is
            _, self.bulk_vertex_info = bulk_qp_shift(self.m_name, dft_gap)
            return
        self.bulk_shift, self.bulk_vertex_info = bulk_qp_shift(self.m_name, dft_gap)
        self.qp_bulk_ref = self.pbe_bulk_gap + self.bulk_shift

    def _resta_w(self, eps):
        from qdex.constants import ANG_PER_BOHR
        k_s_ang = (np.sqrt(max(0.0, eps - 1.0)) / self.d_NN_au) / ANG_PER_BOHR
        c_inf = 1.0 / max(1.0, eps)
        return (c_inf + (1.0 - c_inf) * np.exp(-k_s_ang * self.R_mat_ang)) * self.gamma_bare_ev

    def _penn(self, d_e_conf):
        return penn_eps_eff(self.eps_bulk, d_e_conf, penn_gap_ev(self.m_name, self.eps_bulk))

    def _dim(self, alpha):
        S_dim, _, _, _ = build_dim_screening_factors(coords=self.coords, atom_symbols=self.atom_symbols,
                                                     material_name=self.m_name, eps_out=self.eps_out,
                                                     alpha=alpha)
        inter_mask = ~np.eye(self.n_atoms, dtype=bool)
        eps_z = float(1.0 / np.median(S_dim[inter_mask])) if np.any(inter_mask) else self.eps_bulk
        return S_dim * self.gamma_bare_ev, eps_z

    def _dim_scaled_alpha(self, gap):
        return float(np.clip(self.alpha * (self.qp_bulk_ref / max(0.5, gap)), 0.20, 1.0))

    def one_shot(self, gap_dft):
        """W^QD and the eps entering Z for a one-shot correction at the DFT gap."""
        if self.kind == "dim":
            return self._dim(self.alpha)
        if self.penn_scaling and gap_dft > self.pbe_bulk_gap:
            eps = self._penn(float(gap_dft - self.pbe_bulk_gap))
            return self._resta_w(eps), eps
        return self.W_bulk_ev, self.eps_bulk

    def at_gap(self, gap, strict_penn=False):
        """W^QD, eps for Z and an iteration note at the current QP gap (evGW, qsGW).

        strict_penn: Penn scaling only above the bulk GW gap (qsGW) instead of max(0, gap - E_g^GW) (evGW).
        """
        if self.kind == "dim":
            a = self._dim_scaled_alpha(gap)
            w, eps_z = self._dim(a)
            return w, eps_z, f"alpha_scale = {a:.3f}"
        if self.penn_scaling and (gap > self.qp_bulk_ref or not strict_penn):
            eps = self._penn(max(0.0, gap - self.qp_bulk_ref))
            return self._resta_w(eps), eps, f"eps_eff = {eps:.3f}"
        return self.W_bulk_ev, self.eps_bulk, f"eps_eff = {self.eps_bulk:.3f}"

    def vacuum_solvent(self):
        return solvent_reaction_term(self.coords, self.atom_symbols, self.m_name, 1.0, self.eps_bulk,
                                     self.R_mat_ang, self.R_QD_ang, self.solvent_term)

    def w_parts(self, w_qd, eps_z):
        return {"w_qd": np.array(w_qd, dtype=float), "w_bulk": np.array(self.W_bulk_ev, dtype=float),
                "w_add": np.array(self.delta_W_solv, dtype=float), "gamma": np.array(self.gamma_bare_ev, dtype=float),
                "eps_z": float(eps_z), "bulk_shift": float(self.bulk_shift),
                "bulk_homo_fraction": anchor_bulk_homo_fraction(self.m_name)}


def _estimate_sgw_delta_w(model, dft_gap=None, C_occ_low=None, C_virt_low=None, eps_occ=None, eps_virt=None,
                          atom_ao_ranges=None, Z=0.8, dynamic_z=False, self_consistent=False, max_iter=25,
                          tol=1e-4, damping=0.5, return_details=False, frontier_pops=None):
    """One-shot sGW or eigenvalue self-consistent evGW scissor on the Delta W of ``model``."""
    m = model
    m_name, eps_bulk, bulk_shift = m.m_name, m.eps_bulk, m.bulk_shift
    pbe_bulk_gap, gw_bulk_gap, R_QD_ang = m.pbe_bulk_gap, m.gw_bulk_gap, m.R_QD_ang
    W_bulk_ev, delta_W_solv = m.W_bulk_ev, m.delta_W_solv
    coords, n_atoms = m.coords, m.n_atoms

    # Frontier densities
    if frontier_pops is not None:
        # atomic populations of HOMO and LUMO supplied by the caller (Mulliken or Löwdin)
        q_h, q_l = (np.asarray(x, dtype=float) for x in frontier_pops)
    elif C_occ_low is not None and C_virt_low is not None and atom_ao_ranges is not None:
        q_h = np.zeros(n_atoms)
        q_l = np.zeros(n_atoms)
        for A, (a0, a1) in enumerate(atom_ao_ranges):
            q_h[A] = np.sum(np.abs(C_occ_low[a0:a1, -1])**2)
            q_l[A] = np.sum(np.abs(C_virt_low[a0:a1, 0])**2)
    else:
        center = np.mean(coords, axis=0)
        dist = np.linalg.norm(coords - center, axis=1)
        psi_env = np.maximum(0.0, np.cos(np.pi * dist / (2.0 * max(1.0, R_QD_ang))))
        q_h = (psi_env**2) / np.sum(psi_env**2)
        q_l = q_h.copy()

    if dft_gap is not None:
        gap_dft = float(dft_gap)
    elif eps_occ is not None and eps_virt is not None:
        gap_dft = float(eps_virt[0] - eps_occ[-1])
    else:
        gap_dft = pbe_bulk_gap
    m.set_dft_gap(gap_dft)
    bulk_shift = m.bulk_shift

    def sigmas(w_qd):
        delta_W_conf = np.maximum(0.0, w_qd - W_bulk_ev)
        return (0.5 * float(q_h @ delta_W_conf @ q_h), 0.5 * float(q_l @ delta_W_conf @ q_l),
                0.5 * float(q_h @ delta_W_solv @ q_h), 0.5 * float(q_l @ delta_W_solv @ q_l))

    def z_factors(sig_h_stat, sig_l_stat, gap, eps_z):
        if dynamic_z:
            return (compute_dynamic_z(sig_h_stat, gap, eps_z, m_name),
                    compute_dynamic_z(sig_l_stat, gap, eps_z, m_name))
        return float(Z), float(Z)

    model_name = f"{'evGW' if self_consistent else 'sGW'}-{m.label}"
    if not self_consistent:
        W_qd_ev, eps_z = m.one_shot(gap_dft)
        sig_h_conf, sig_l_conf, sig_h_solv, sig_l_solv = sigmas(W_qd_ev)
        sig_h_stat = sig_h_conf + sig_h_solv
        sig_l_stat = sig_l_conf + sig_l_solv
        Z_h, Z_l = z_factors(sig_h_stat, sig_l_stat, gap_dft + bulk_shift, eps_z)
        delta_sigma_h = Z_h * sig_h_stat
        delta_sigma_l = Z_l * sig_l_stat
        confinement_shift = delta_sigma_h + delta_sigma_l
        confinement_shift_internal = Z_h * sig_h_conf + Z_l * sig_l_conf
        confinement_shift_solv = Z_h * sig_h_solv + Z_l * sig_l_solv
        total_sgw_scissor = bulk_shift + confinement_shift
        n_iters = 1
        converged = True
    else:
        gap_curr = gap_dft + bulk_shift
        logger.info(f"\n  [{model_name}] Starting Eigenvalue Self-Consistent Loop (Initial Gap = {gap_curr:.4f} eV):")
        n_iters = 0
        converged = False
        scissor_next = bulk_shift
        delta_sigma_h = delta_sigma_l = 0.0
        Z_h = Z_l = float(Z)
        eps_z = eps_bulk
        for it in range(max_iter):
            n_iters += 1
            W_qd_ev, eps_z, note = m.at_gap(gap_curr)
            sig_h_conf, sig_l_conf, sig_h_solv, sig_l_solv = sigmas(W_qd_ev)
            sig_h_stat = sig_h_conf + sig_h_solv
            sig_l_stat = sig_l_conf + sig_l_solv
            Z_h, Z_l = z_factors(sig_h_stat, sig_l_stat, gap_curr, eps_z)
            delta_sigma_h = Z_h * sig_h_stat
            delta_sigma_l = Z_l * sig_l_stat
            confinement_shift = delta_sigma_h + delta_sigma_l
            scissor_next = bulk_shift + confinement_shift
            gap_next = gap_dft + scissor_next

            diff = abs(gap_next - gap_curr)
            logger.debug(f"    Iter {it+1:2d}: Gap = {gap_curr:.4f} eV, Scissor = +{scissor_next:.4f} eV, {note}, Z_h = {Z_h:.3f}, Z_l = {Z_l:.3f}, Diff = {diff:.5f} eV")
            if diff < tol:
                converged = True
                logger.info(f"    -> {model_name} converged in {it+1} iterations! Final QP Gap = {gap_next:.4f} eV")
                break
            gap_curr = (1.0 - damping) * gap_curr + damping * gap_next
        total_sgw_scissor = scissor_next
        confinement_shift_internal = Z_h * sig_h_conf + Z_l * sig_l_conf
        confinement_shift_solv = Z_h * sig_h_solv + Z_l * sig_l_solv

    int_shift = confinement_shift_internal
    if int_shift > 1e-8:
        f_h_micro = float(np.clip(Z_h * sig_h_conf / int_shift, 0.0, 1.0))
        f_l_micro = 1.0 - f_h_micro
    elif delta_sigma_h + delta_sigma_l > 1e-8:
        f_h_micro = float(np.clip(delta_sigma_h / (delta_sigma_h + delta_sigma_l), 0.0, 1.0))
        f_l_micro = 1.0 - f_h_micro
    else:
        f_h_micro = 0.5
        f_l_micro = 0.5

    penn_active = m.kind == "resta" and m.penn_scaling and eps_z < eps_bulk
    z_str = f"Z_h={Z_h:.3f}, Z_l={Z_l:.3f} (dynamic)" if dynamic_z else f"Z={float(Z):.2f} (fixed)"
    if m.kind == "dim":
        logger.info(f"\n  [{model_name} Model] Quasiparticle Correction for {m_name}:")
    else:
        logger.info(f"\n  [{model_name} Model ({'Penn-scaled' if penn_active else 'Pure Boundary'})] "
              f"Quasiparticle Correction for {m_name}:")
    logger.info(f"    Cluster Radius (R_QD)    : {R_QD_ang:.3f} Å")
    logger.info(f"    Bulk PBE -> GW Gap      : {pbe_bulk_gap:.3f} -> {gw_bulk_gap:.3f} eV (Shift: +{bulk_shift:.3f} eV)")
    if m.kind == "dim":
        logger.info(f"    Internal DIM Contrast    : +{confinement_shift_internal:.3f} eV (surface coordination under-screening)")
        logger.info(f"    Solvent Reaction Field   : +{confinement_shift_solv:.3f} eV (eps_out = {m.eps_out_val:.2f}, eps_bulk = {eps_bulk:.2f})")
    else:
        if penn_active:
            logger.info(f"    Penn Dielectric eps_eff  : {eps_z:.3f} (bulk eps_inf = {eps_bulk:.3f})")
            logger.info(f"    Internal Resta Contrast  : +{confinement_shift_internal:.3f} eV")
        logger.info(f"    Solvent Reaction Field   : +{confinement_shift_solv:.3f} eV (eps_out = {m.eps_out_val:.2f})")
    logger.info(f"    HOMO Quasiparticle Shift : +{delta_sigma_h:.3f} eV ({z_str})")
    logger.info(f"    LUMO Quasiparticle Shift : +{delta_sigma_l:.3f} eV ({z_str})")
    logger.info(f"    Microscopic Split        : HOMO takes {f_h_micro*100:.1f}%, LUMO takes {f_l_micro*100:.1f}%")
    logger.info(f"    Total Confinement Opening: +{confinement_shift:.3f} eV")
    if self_consistent:
        logger.info(f"    evGW Iterations          : {n_iters} (converged: {converged})")
    logger.info(f"    ==> Total {model_name} Scissor: +{total_sgw_scissor:.3f} eV\n")

    w_solv_vac = m.vacuum_solvent()
    provenance = {
        "material": m_name,
        "qp_model": f"{'evgw' if self_consistent else 'sgw'}_{m.kind}",
        "cluster_radius_ang": R_QD_ang,
        "bulk_pbe_gap_ev": pbe_bulk_gap,
        "bulk_gw_gap_ev": gw_bulk_gap,
        "bulk_shift_ev": bulk_shift,
        "bulk_gw_shift_ev": bulk_shift,
        "confinement_shift_ev": confinement_shift,
        "confinement_shift_internal_ev": confinement_shift_internal,
        "confinement_shift_solvent_ev": confinement_shift_solv,
        "delta_sigma_homo_ev": delta_sigma_h,
        "delta_sigma_lumo_ev": delta_sigma_l,
        "f_homo": f_h_micro,
        "f_lumo": f_l_micro,
        "f_homo_micro": f_h_micro,
        "f_lumo_micro": f_l_micro,
        "z_factor": float(0.5 * (Z_h + Z_l)),
        "z_homo": float(Z_h),
        "z_lumo": float(Z_l),
        "dynamic_z": bool(dynamic_z),
        "self_consistent": bool(self_consistent),
        "evgw_converged": bool(converged),
        "evgw_iterations": int(n_iters),
        "eps_out": float(m.eps_out_val),
        "eps_bulk": float(eps_bulk),
    }
    if m.kind == "resta":
        provenance["eps_eff_qd"] = float(eps_z)
        provenance["penn_scaling"] = bool(penn_active)
    provenance.update(m.bulk_vertex_info)
    provenance.update({
        "total_scissor_ev": total_sgw_scissor,
        "total_scissor_solvent_ev": total_sgw_scissor,
        "total_scissor_vacuum_ev": bulk_shift + confinement_shift_internal
        + (0.5 * float(Z_h * (q_h @ w_solv_vac @ q_h) + Z_l * (q_l @ w_solv_vac @ q_l))),
    })
    # The screened interaction this QP model is built on.  The CLI passes it
    # to the BSE (kernel "qp") so that GW and BSE share the same W.
    # Delta W enters the kernel with the same Z as the QP shift (scale_w_difference).
    provenance["w_bse_ev"] = scale_w_difference(W_qd_ev + delta_W_solv, W_bulk_ev, 0.5 * (Z_h + Z_l))
    provenance["w_bse_z"] = float(0.5 * (Z_h + Z_l))
    # Components for the xs representation and for orbital-resolved QP levels.
    provenance["w_parts"] = m.w_parts(W_qd_ev, eps_z)
    provenance["w_bse_label"] = (f"DIM W_QD + {m.solvent_term} solvent term" if m.kind == "dim"
                                 else f"Resta W_QD(eps_eff) + {m.solvent_term} solvent term")

    if return_details:
        return total_sgw_scissor, provenance
    return total_sgw_scissor


def estimate_sgw_dim_qp_gap(coords, atom_symbols, material_name=None, eps_out=2.4,
                            C_occ_low=None, C_virt_low=None, eps_occ=None, eps_virt=None,
                            atom_ao_ranges=None, alpha=1.0, Z=0.8, dynamic_z=False,
                            self_consistent=False, max_iter=25, tol=1e-4, damping=0.5,
                            return_details=False, frontier_pops=None, solvent_term="sphere"):
    """
    Computes the Quasiparticle (GW) gap shift for a quantum dot via the microscopic
    screening difference (Delta W) approach using the Atomistic Polarizable Dipole Interaction
    Model (DIM / Thole Model)::

        Delta Sigma_p = (Z_p / 2) * <psi_p | W^{QD, DIM} - W^{bulk} + W^{solv} | psi_p>
        Scissor_{sGW-DIM} = (E_g^{bulk, GW} - E_g^{bulk, PBE}) + Delta Sigma_{HOMO} + Delta Sigma_{LUMO}

    Options:
      - dynamic_z=True: Dynamically computes Z_p from the microscopic plasmon-pole f-sum rule.
      - self_consistent=True (evGW): Solves the eigenvalue self-consistent Dyson equation
        E_g^{(k+1)} = E_g^DFT + Scissor(E_g^{(k)}) until convergence; the DIM polarizabilities
        are scaled by alpha * E_g^{bulk,GW} / E_g (clipped to [0.2, 1]).
    """
    model = _ScreeningModel("dim", coords, atom_symbols, material_name, eps_out, solvent_term, alpha=alpha)
    return _estimate_sgw_delta_w(model, None, C_occ_low, C_virt_low, eps_occ, eps_virt, atom_ao_ranges, Z,
                                 dynamic_z, self_consistent, max_iter, tol, damping, return_details, frontier_pops)


def estimate_sgw_resta_qp_gap(coords, atom_symbols, material_name=None, eps_out=2.4,
                              dft_gap=None, C_occ_low=None, C_virt_low=None,
                              eps_occ=None, eps_virt=None, atom_ao_ranges=None,
                              alpha=1.0, Z=0.8, penn_scaling=True, dynamic_z=False,
                              self_consistent=False, max_iter=25, tol=1e-4, damping=0.5,
                              return_details=False, frontier_pops=None, solvent_term="sphere"):
    """
    Computes the Quasiparticle (GW) gap shift for a quantum dot via the microscopic
    screening difference (Delta W) approach using the Resta dielectric screening model::

        Delta Sigma_p = (Z_p / 2) * <psi_p | W^{QD, Resta} - W^{bulk} + W^{solv} | psi_p>
        Scissor_{sGW-Resta} = (E_g^{bulk, GW} - E_g^{bulk, PBE}) + Delta Sigma_{HOMO} + Delta Sigma_{LUMO}

    Options:
      - penn_scaling=True (default): Scales the QD dielectric constant via the Penn model
        eps_eff = 1 + (eps_inf - 1) / (1 + (Delta E_conf / E_g^bulk)^2).
      - dynamic_z=True: Dynamically computes Z_p from the microscopic plasmon-pole f-sum rule.
      - self_consistent=True (evGW): Solves the eigenvalue self-consistent Dyson equation.
    """
    model = _ScreeningModel("resta", coords, atom_symbols, material_name, eps_out, solvent_term, alpha=alpha,
                            penn_scaling=penn_scaling)
    return _estimate_sgw_delta_w(model, dft_gap, C_occ_low, C_virt_low, eps_occ, eps_virt, atom_ao_ranges, Z,
                                 dynamic_z, self_consistent, max_iter, tol, damping, return_details, frontier_pops)


def estimate_evgw_dim_qp_gap(coords, atom_symbols, material_name=None, eps_out=2.4,
                             C_occ_low=None, C_virt_low=None, eps_occ=None, eps_virt=None,
                             atom_ao_ranges=None, alpha=1.0, max_iter=25, tol=1e-4,
                             damping=0.5, return_details=False):
    """Convenience wrapper for Eigenvalue Self-Consistent evGW using the Atomistic DIM kernel and dynamic Z."""
    return estimate_sgw_dim_qp_gap(
        coords=coords, atom_symbols=atom_symbols, material_name=material_name, eps_out=eps_out,
        C_occ_low=C_occ_low, C_virt_low=C_virt_low, eps_occ=eps_occ, eps_virt=eps_virt,
        atom_ao_ranges=atom_ao_ranges, alpha=alpha, dynamic_z=True, self_consistent=True,
        max_iter=max_iter, tol=tol, damping=damping, return_details=return_details
    )


def estimate_evgw_resta_qp_gap(coords, atom_symbols, material_name=None, eps_out=2.4,
                               dft_gap=None, C_occ_low=None, C_virt_low=None,
                               eps_occ=None, eps_virt=None, atom_ao_ranges=None,
                               alpha=1.0, max_iter=25, tol=1e-4,
                               damping=0.5, return_details=False):
    """Convenience wrapper for Eigenvalue Self-Consistent evGW using the Resta-Penn kernel and dynamic Z."""
    return estimate_sgw_resta_qp_gap(
        coords=coords, atom_symbols=atom_symbols, material_name=material_name, eps_out=eps_out,
        dft_gap=dft_gap, C_occ_low=C_occ_low, C_virt_low=C_virt_low,
        eps_occ=eps_occ, eps_virt=eps_virt, atom_ao_ranges=atom_ao_ranges,
        alpha=alpha, penn_scaling=True, dynamic_z=True, self_consistent=True,
        max_iter=max_iter, tol=tol, damping=damping, return_details=return_details
    )


def _estimate_qsgw_delta_w(model, C, eps, S, atom_ao_ranges, homo_index, dynamic_z=True, max_iter=25, tol=1e-4,
                           damping=0.5, return_details=False, Z=0.8, gamma_ao=None):
    """Full-AO qsGW (eigenvalues and orbitals) on the Delta W of ``model``."""
    m = model
    m_name, eps_bulk, bulk_shift = m.m_name, m.eps_bulk, m.bulk_shift
    W_bulk_ev, delta_W_solv, gamma_bare_ev = m.W_bulk_ev, m.delta_W_solv, m.gamma_bare_ev
    tag = f"qsGW-{m.label}"
    n_ao = C.shape[0]

    # Löwdin orthogonalization operators
    S_dense = S.toarray() if hasattr(S, "toarray") else np.asarray(S, dtype=np.float64)
    from qdex.lowdin import lowdin_factor
    eigvals_S, U_S = lowdin_factor(S_dense)
    eigvals_S = np.maximum(eigvals_S, 1e-12)
    S_half = (U_S * np.sqrt(eigvals_S)[None, :]) @ U_S.T
    S_inv_half = (U_S * (1.0 / np.sqrt(eigvals_S))[None, :]) @ U_S.T

    C_dense = C.toarray() if hasattr(C, "toarray") else np.asarray(C, dtype=np.float64)
    if C_dense.shape[1] != n_ao:
        raise ValueError(
            f"Full-AO qsGW requires complete square MO coefficients (got shape {C_dense.shape} for {n_ao} AOs). "
            f"Please provide full MOs or use non-orbital QP methods (e.g. sgw-{m.kind}, evgw-{m.kind})."
        )
    C_low_init = S_half @ C_dense
    eps_dft = np.asarray(eps, dtype=np.float64)
    H_dft_low = (C_low_init * eps_dft[None, :]) @ C_low_init.T

    n_occ = homo_index + 1
    lumo_index = homo_index + 1
    dft_gap = float(eps_dft[lumo_index] - eps_dft[homo_index])
    m.set_dft_gap(dft_gap)
    bulk_shift = m.bulk_shift
    gap_curr = dft_gap + bulk_shift

    C_curr = C_low_init.copy()
    H_eff_prev = None
    n_iters = 0
    converged = False

    logger.info(f"\n  [{tag}] Starting Full AO Quasiparticle Self-Consistent Loop ({n_ao}x{n_ao} AOs, Initial Gap = {gap_curr:.4f} eV):")
    for it in range(max_iter):
        n_iters += 1
        P_low = 2.0 * (C_curr[:, :n_occ] @ C_curr[:, :n_occ].conj().T)
        Q_low = np.eye(n_ao) - 0.5 * P_low

        W_qd_ev, eps_z, note = m.at_gap(gap_curr, strict_penn=True)
        delta_W_atom = np.maximum(0.0, W_qd_ev - W_bulk_ev) + delta_W_solv

        if gamma_ao is None:
            delta_W_ao = expand_atom_to_ao(delta_W_atom, atom_ao_ranges, n_ao)
            q_h = np.array([np.sum(np.abs(C_curr[a0:a1, homo_index])**2) for a0, a1 in atom_ao_ranges])
            q_l = np.array([np.sum(np.abs(C_curr[a0:a1, lumo_index])**2) for a0, a1 in atom_ao_ranges])
            sig_h_stat = 0.5 * float(q_h @ delta_W_atom @ q_h)
            sig_l_stat = 0.5 * float(q_l @ delta_W_atom @ q_l)
        else:
            # xs: the same screening ratio applied to exact AO density-pair integrals
            delta_W_ao = ao_delta_w(np.maximum(0.0, W_qd_ev - W_bulk_ev), gamma_bare_ev, delta_W_solv,
                                    gamma_ao, atom_ao_ranges)
            q_h = np.abs(C_curr[:, homo_index])**2
            q_l = np.abs(C_curr[:, lumo_index])**2
            sig_h_stat = 0.5 * float(q_h @ delta_W_ao @ q_h)
            sig_l_stat = 0.5 * float(q_l @ delta_W_ao @ q_l)

        tot_sig = sig_h_stat + sig_l_stat
        if tot_sig > 1e-8:
            f_h_iter = float(np.clip(sig_h_stat / tot_sig, 0.0, 1.0))
            f_l_iter = 1.0 - f_h_iter
        else:
            f_h_iter = 0.5
            f_l_iter = 0.5
        H_bulk_low = - (f_h_iter * bulk_shift) * (0.5 * P_low) + (f_l_iter * bulk_shift) * Q_low

        if dynamic_z:
            Z_h = compute_dynamic_z(sig_h_stat, gap_curr, eps_z, m_name)
            Z_l = compute_dynamic_z(sig_l_stat, gap_curr, eps_z, m_name)
            Z_avg = 0.5 * (Z_h + Z_l)
        else:
            Z_h = float(Z)
            Z_l = float(Z)
            Z_avg = float(Z)

        Sigma_sex = -0.5 * P_low * delta_W_ao
        Sigma_coh = 0.5 * np.diag(np.diag(delta_W_ao))

        H_eff_target = H_dft_low + H_bulk_low + Z_avg * (Sigma_sex + Sigma_coh)
        if H_eff_prev is None:
            H_eff = H_eff_target
        else:
            H_eff = (1.0 - damping) * H_eff_prev + damping * H_eff_target
        H_eff_prev = H_eff.copy()

        eigvals_qp, C_new = np.linalg.eigh(H_eff)

        # Enforce phase consistency
        signs = np.sign(np.sum(C_curr * C_new, axis=0))
        signs[signs == 0] = 1.0
        C_new = C_new * signs[None, :]

        gap_next = float(eigvals_qp[lumo_index] - eigvals_qp[homo_index])
        diff = abs(gap_next - gap_curr)

        fid_h = float((C_low_init[:, homo_index] @ C_new[:, homo_index])**2)
        fid_l = float((C_low_init[:, lumo_index] @ C_new[:, lumo_index])**2)

        note_str = f", {note}" if m.kind == "resta" else ""
        logger.debug(f"    Iter {it+1:2d}: Gap = {gap_next:.4f} eV, Scissor = +{gap_next - dft_gap:.4f} eV{note_str}, Diff = {diff:.5f} eV, Zh={Z_h:.3f}, Zl={Z_l:.3f}, Fid_H={fid_h:.5f}, Fid_L={fid_l:.5f}")

        if diff < tol:
            converged = True
            logger.info(f"    -> {tag} converged in {it+1} iterations! Final QP Gap = {gap_next:.4f} eV")
            break
        C_curr = C_new
        gap_curr = gap_next
    else:
        logger.info(f"    -> {tag} reached max iterations ({max_iter}). Final QP Gap = {gap_curr:.4f} eV")

    C_qp = S_inv_half @ C_new
    eps_qp = eigvals_qp.copy()
    final_scissor = float(eps_qp[lumo_index] - eps_qp[homo_index] - dft_gap)
    delta_sigma_h = float(eps_qp[homo_index] - eps_dft[homo_index])
    delta_sigma_l = float(eps_qp[lumo_index] - eps_dft[lumo_index] - bulk_shift)

    h_shift_tot = abs(float(eps_qp[homo_index] - eps_dft[homo_index]))
    l_shift_tot = abs(float(eps_qp[lumo_index] - eps_dft[lumo_index]))
    tot_shift = max(1e-12, h_shift_tot + l_shift_tot)
    f_h_final = float(h_shift_tot / tot_shift)
    f_l_final = float(l_shift_tot / tot_shift)

    provenance = {
        "material": m_name,
        "qp_model": f"qsgw_{m.kind}",
        "cluster_radius_ang": m.R_QD_ang,
        "bulk_pbe_gap_ev": m.pbe_bulk_gap,
        "bulk_gw_gap_ev": m.gw_bulk_gap,
        "bulk_shift_ev": bulk_shift,
        "confinement_shift_ev": final_scissor - bulk_shift,
        "delta_sigma_homo_ev": delta_sigma_h,
        "delta_sigma_lumo_ev": delta_sigma_l,
        "f_homo": f_h_final,
        "f_lumo": f_l_final,
        "f_homo_micro": f_h_final,
        "f_lumo_micro": f_l_final,
        "z_factor": float(Z_avg),
        "z_homo": float(Z_h),
        "z_lumo": float(Z_l),
        "dynamic_z": bool(dynamic_z),
        "self_consistent": True,
        "orbital_update": True,
        "qsgw_converged": bool(converged),
        "qsgw_iterations": int(n_iters),
        "homo_fidelity": float(fid_h),
        "lumo_fidelity": float(fid_l),
        "eps_out": float(m.eps_out_val),
        "eps_bulk": float(eps_bulk),
    }
    if m.kind == "resta":
        provenance["eps_eff_qd"] = float(eps_z)
    provenance.update(m.bulk_vertex_info)
    provenance["total_scissor_ev"] = final_scissor
    # The screened interaction this QP model is built on.  The CLI passes it
    # to the BSE (kernel "qp") so that GW and BSE share the same W.
    # Delta W enters the kernel with the same Z as the QP shift (scale_w_difference).
    provenance["w_bse_ev"] = scale_w_difference(W_qd_ev + delta_W_solv, W_bulk_ev, Z_avg)
    provenance["w_bse_z"] = float(Z_avg)
    # Components for the xs representation and for orbital-resolved QP levels.
    provenance["w_parts"] = m.w_parts(W_qd_ev, eps_z)
    provenance["w_bse_label"] = ("DIM W_QD (final iteration) + solvent term" if m.kind == "dim"
                                 else "Resta W_QD(eps_eff, final iteration) + solvent term")

    if return_details:
        return final_scissor, provenance, C_qp, eps_qp
    return final_scissor


def estimate_qsgw_dim_qp_gap(coords, atom_symbols, C, eps, S, atom_ao_ranges, homo_index,
                             material_name=None, eps_out=2.4, alpha=1.0, dynamic_z=True,
                             max_iter=25, tol=1e-4, damping=0.5, return_details=False, Z=0.8,
                             gamma_ao=None, solvent_term="sphere"):
    """
    Computes Quasiparticle Self-Consistent GW (qsGW) by updating BOTH eigenvalues
    and molecular orbitals across the full AO basis using the Atomistic Polarizable
    Dipole Interaction Model (DIM / Thole Model)::

        H_eff = H_DFT + Delta H_bulk + Z (Sigma^SEX + Sigma^COH)
        H_eff C_new = S C_new E_new

    Returns:
        total_scissor, provenance, C_qp, eps_qp
    """
    model = _ScreeningModel("dim", coords, atom_symbols, material_name, eps_out, solvent_term, alpha=alpha)
    return _estimate_qsgw_delta_w(model, C, eps, S, atom_ao_ranges, homo_index, dynamic_z, max_iter, tol, damping,
                                  return_details, Z, gamma_ao)


def estimate_qsgw_resta_qp_gap(coords, atom_symbols, C, eps, S, atom_ao_ranges, homo_index,
                               material_name=None, eps_out=2.4, alpha=1.0, penn_scaling=True,
                               dynamic_z=True, max_iter=25, tol=1e-4, damping=0.5,
                               return_details=False, Z=0.8, gamma_ao=None, solvent_term="sphere"):
    """
    Computes Quasiparticle Self-Consistent GW (qsGW) by updating BOTH eigenvalues
    and molecular orbitals across the full AO basis using the Resta dielectric model::

        H_eff = H_DFT + Delta H_bulk + Z (Sigma^SEX + Sigma^COH)
        H_eff C_new = S C_new E_new

    Returns:
        total_scissor, provenance, C_qp, eps_qp
    """
    model = _ScreeningModel("resta", coords, atom_symbols, material_name, eps_out, solvent_term, alpha=alpha,
                            penn_scaling=penn_scaling)
    return _estimate_qsgw_delta_w(model, C, eps, S, atom_ao_ranges, homo_index, dynamic_z, max_iter, tol, damping,
                                  return_details, Z, gamma_ao)


def build_xs_kernel(shells, atom_symbols, coords, atom_ao_ranges, material_name=None, kernel_mode="bse", alpha=1.0, nthreads=1,
                    C_occ_low=None, C_virt_low=None, eps_occ=None, eps_virt=None,
                    C_occ_b_low=None, C_virt_b_low=None, eps_occ_b=None, eps_virt_b=None,
                    eps_out=2.4, return_eps_info=False):
    """
    Builds the exact AO two-electron integral kernel (Xs-QDEX) for BSE/TDA.
    Gamma_{mu, nu} = (mu mu | nu nu) evaluated analytically via Libint2.

    Screening Modes:
      1. kernel_mode in ['dim', 'dipole', 'xs-dim', 'xs-dipole']:
         Atomistic Polarizable Dipole Interaction Model (DIM / Thole Model),
         solving M * p = E_ext with M = alpha^-1 + T.
         S_AB computed from atom-specific polarizabilities and local dielectric response.
         W_{mu, nu} = S_{AB} * Gamma_{mu, nu}.
      2. kernel_mode in ['rpa', 'xs-rpa', 'xs_rpa', 'zdo-rpa']:
         Parameter-free microscopic Random Phase Approximation (RPA) screening
         under the Zero Differential Overlap (ZDO) approximation.
      3. kernel_mode in ['sbse', 'sbse-ao', 'sbse_ao']:
         Simplified Bethe-Salpeter Equation (sBSE) AO-resolved kernel (Cho et al. 2022).
      4. kernel_mode in ['resta', 'xs-resta']:
         Screens Gamma using the microscopic Resta screening profile.
      5. kernel_mode == 'bse' / 'uniform':
         W_{mu, nu} = alpha * Gamma_{mu, nu}

    Returns:
        (W_ev, W_ev, Gamma_bare_ev) or (W_ev, W_ev, Gamma_bare_ev, eps_info) if return_eps_info=True
    """
    km = str(kernel_mode).lower()
    if km in ["sbse", "sbse-ao", "sbse_ao"]:
        return build_sbse_kernel(
            atom_symbols=atom_symbols,
            coords=coords,
            atom_ao_ranges=atom_ao_ranges,
            shells=shells,
            C_occ_low=C_occ_low,
            C_virt_low=C_virt_low,
            eps_occ=eps_occ,
            eps_virt=eps_virt,
            C_occ_b_low=C_occ_b_low,
            C_virt_b_low=C_virt_b_low,
            eps_occ_b=eps_occ_b,
            eps_virt_b=eps_virt_b,
            mode="ao",
            eps_out=eps_out,
            material_name=material_name,
            alpha=alpha,
            nthreads=nthreads,
            return_eps_info=return_eps_info
        )

    from qdex.integrals import compute_two_electron_ao
    from qdex.constants import HA_TO_EV, ANG_PER_BOHR
    from scipy.spatial.distance import pdist, squareform

    # 1. Compute exact bare AO repulsion matrix (Hartree -> eV)
    gamma_bare_au = compute_two_electron_ao(shells, nthreads=nthreads)
    gamma_bare_ev = gamma_bare_au * HA_TO_EV

    n_ao = gamma_bare_ev.shape[0]
    m_name = material_name.upper() if material_name else "DEFAULT"
    km = str(kernel_mode).lower()

    if km in ["rpa", "xs-rpa", "xs_rpa", "zdo-rpa"]:
        if C_occ_low is None or C_virt_low is None or eps_occ is None or eps_virt is None:
            raise ValueError(
                "ZDO-RPA screening requires active Löwdin MOs (C_occ_low, C_virt_low) and eigenvalues (eps_occ, eps_virt)."
            )

        # 1. Non-interacting polarizability and microscopic dielectric screening
        coords_bohr = (coords / ANG_PER_BOHR) if coords is not None else np.zeros((len(atom_ao_ranges), 3))
        n_occ_a = C_occ_low.shape[1]
        n_virt_a = C_virt_low.shape[1]
        d_alpha_au = (eps_virt[:, None] - eps_occ[None, :]) / HA_TO_EV
        d_alpha_au = np.maximum(d_alpha_au, 1e-6)

        # Quantum mechanical transition dipoles: d_ia = sum_A R_A * sum_{mu in A} C_occ[mu, i] * C_virt[mu, a]
        d_ia_x_a = np.zeros((n_occ_a, n_virt_a), dtype=np.float64)
        d_ia_y_a = np.zeros((n_occ_a, n_virt_a), dtype=np.float64)
        d_ia_z_a = np.zeros((n_occ_a, n_virt_a), dtype=np.float64)
        for A, (a0, a1) in enumerate(atom_ao_ranges):
            q_ia_A = np.einsum('mi,ma->ia', C_occ_low[a0:a1, :], C_virt_low[a0:a1, :])
            d_ia_x_a += coords_bohr[A, 0] * q_ia_A
            d_ia_y_a += coords_bohr[A, 1] * q_ia_A
            d_ia_z_a += coords_bohr[A, 2] * q_ia_A

        ax_a = np.sum(d_ia_x_a**2 / d_alpha_au.T)
        ay_a = np.sum(d_ia_y_a**2 / d_alpha_au.T)
        az_a = np.sum(d_ia_z_a**2 / d_alpha_au.T)
        alpha_iso_a = (ax_a + ay_a + az_a) / 3.0

        if C_occ_b_low is not None and C_virt_b_low is not None and eps_occ_b is not None and eps_virt_b is not None:
            n_occ_b = C_occ_b_low.shape[1]
            n_virt_b = C_virt_b_low.shape[1]
            d_beta_au = (eps_virt_b[:, None] - eps_occ_b[None, :]) / HA_TO_EV
            d_beta_au = np.maximum(d_beta_au, 1e-6)
            d_ia_x_b = np.zeros((n_occ_b, n_virt_b), dtype=np.float64)
            d_ia_y_b = np.zeros((n_occ_b, n_virt_b), dtype=np.float64)
            d_ia_z_b = np.zeros((n_occ_b, n_virt_b), dtype=np.float64)
            for A, (a0, a1) in enumerate(atom_ao_ranges):
                q_ia_B = np.einsum('mi,ma->ia', C_occ_b_low[a0:a1, :], C_virt_b_low[a0:a1, :])
                d_ia_x_b += coords_bohr[A, 0] * q_ia_B
                d_ia_y_b += coords_bohr[A, 1] * q_ia_B
                d_ia_z_b += coords_bohr[A, 2] * q_ia_B
            ax_b = np.sum(d_ia_x_b**2 / d_beta_au.T)
            ay_b = np.sum(d_ia_y_b**2 / d_beta_au.T)
            az_b = np.sum(d_ia_z_b**2 / d_beta_au.T)
            alpha_iso_b = (ax_b + ay_b + az_b) / 3.0
            # Open-shell UKS: sum of alpha and beta channels
            alpha_cluster_au = alpha_iso_a + alpha_iso_b
            n_trans_tot = (n_occ_a * n_virt_a) + (n_occ_b * n_virt_b)
        else:
            # Closed-shell singlet: factor of 2.0 for spin
            alpha_cluster_au = 2.0 * alpha_iso_a
            n_trans_tot = n_occ_a * n_virt_a

        # 2. Extract equivalent cluster radius and evaluate microscopic dielectric constant
        entry = MATERIAL_DB.get(m_name, MATERIAL_DB["DEFAULT"])
        eps_bulk = float(entry[0])
        r_qd_ang = None
        if coords is not None and len(coords) > 1 and atom_symbols is not None:
            try:
                size_metrics = get_cluster_size_metrics(coords, atom_symbols, m_name)
                r_qd_ang = float(size_metrics.get("R_eff_hull", 0.0))
            except Exception:
                r_qd_ang = None

        if r_qd_ang is None or r_qd_ang <= 0.0:
            r_qd_ang = 5.0
        r_qd_bohr = r_qd_ang / ANG_PER_BOHR
        v_sphere_bohr = r_qd_bohr ** 3

        eta = float(alpha) * (alpha_cluster_au / v_sphere_bohr)
        if eps_bulk > 1.0:
            bulk_eta = (eps_bulk - 1.0) / (eps_bulk + 2.0)
            eta_eff = min(0.98 * bulk_eta, max(0.01, eta))
        else:
            eta_eff = min(0.85, max(0.01, eta))

        eps_eff_micro = (1.0 + 2.0 * eta_eff) / max(0.01, 1.0 - eta_eff)

        # 3. Construct distance-dependent microscopic screening matrix W:
        # Core atoms for nearest-neighbor distance d_NN
        core_coords = coords
        if m_name in MATERIAL_ELEMENTS:
            core_elements = [el.lower() for el in MATERIAL_ELEMENTS[m_name]]
            core_coords = np.array(
                [coords[i] for i, sym in enumerate(atom_symbols) if sym.lower() in core_elements]
            )

        if len(core_coords) > 1:
            r_core_ang = squareform(pdist(core_coords))
            np.fill_diagonal(r_core_ang, np.inf)
            d_NN_ang = np.median(np.min(r_core_ang, axis=1))
        else:
            d_NN_ang = 2.5

        d_NN_au = d_NN_ang / ANG_PER_BOHR
        k_s_au = np.sqrt(max(0.0, eps_eff_micro - 1.0)) / max(1e-4, d_NN_au)
        lambda_s_ang = (1.0 / k_s_au) * ANG_PER_BOHR if k_s_au > 0 else 0.0

        r_mat_au = squareform(pdist(coords_bohr))
        c_inf = 1.0 / eps_eff_micro
        # S_AB = c_inf + (1 - c_inf)*exp(-k_s * r_AB)
        s_atom = c_inf + (1.0 - c_inf) * np.exp(-k_s_au * r_mat_au)

        s_ao = np.zeros((n_ao, n_ao), dtype=np.float64)
        for A, (a0, a1) in enumerate(atom_ao_ranges):
            for B, (b0, b1) in enumerate(atom_ao_ranges):
                s_ao[a0:a1, b0:b1] = s_atom[A, B]

        w_rpa_ev = s_ao * gamma_bare_ev

        # 4. Compute Direct & Self-Energy Screening Diagnostics:
        # Hole density (HOMO) and Electron density (LUMO)
        q_h = np.abs(C_occ_low[:, -1])**2
        q_l = np.abs(C_virt_low[:, 0])**2
        v_eh_bare = float(q_h @ gamma_bare_ev @ q_l)
        w_eh_screened = float(q_h @ w_rpa_ev @ q_l)
        eps_exciton = (v_eh_bare / w_eh_screened) if abs(w_eh_screened) > 1e-12 else eps_eff_micro

        # Hole self-energy (HOMO)
        v_hh = float(q_h @ gamma_bare_ev @ q_h)
        w_hh = float(q_h @ w_rpa_ev @ q_h)
        eps_hole = (v_hh / w_hh) if abs(w_hh) > 1e-12 else 1.0

        # Electron self-energy (LUMO)
        v_ll = float(q_l @ gamma_bare_ev @ q_l)
        w_ll = float(q_l @ w_rpa_ev @ q_l)
        eps_elec = (v_ll / w_ll) if abs(w_ll) > 1e-12 else 1.0

        # On-site / core screening
        diag_bare = np.diag(gamma_bare_ev)
        diag_w = np.diag(w_rpa_ev)
        eps_onsite = float(np.mean(diag_bare / np.maximum(diag_w, 1e-12)))

        # Inter-atomic long-range screening
        inter_mask = np.ones((n_ao, n_ao), dtype=bool)
        for a0, a1 in atom_ao_ranges:
            inter_mask[a0:a1, a0:a1] = False
        if np.any(inter_mask):
            eps_inter = float(np.median(gamma_bare_ev[inter_mask] / np.maximum(w_rpa_ev[inter_mask], 1e-12)))
        else:
            eps_inter = eps_onsite

        alpha_cluster_ang3 = alpha_cluster_au * (ANG_PER_BOHR ** 3)

        logger.info(f"\n    ==========================================================================")
        logger.info(f"    [Kernel: Xs-QDEX (Exact 2-Electron AO Integrals + Microscopic ZDO-RPA)]")
        logger.info(f"    ==========================================================================")
        logger.info(f"    Screening Mode          : Parameter-free Microscopic ZDO-RPA (W = ε⁻¹ v)")
        logger.info(f"    Active Polarizability   : {alpha_cluster_au:8.2f} a.u. ({alpha_cluster_ang3:.2f} Å³) across {n_trans_tot} transitions")
        logger.info(f"    Cluster Radius (R_QD)   : {r_qd_ang:8.3f} Å (R_bohr = {r_qd_bohr:.2f})")
        logger.info(f"    AO Matrix Dimension     : {n_ao} x {n_ao}")
        logger.info(f"    RPA Scaling (alpha)     : {alpha:.3f}")
        logger.info(f"    --------------------------------------------------------------------------")
        logger.info(f"    Computed Microscopic Dielectric Constants (ε_eff):")
        logger.info(f"      ε_eff (1S Exciton e-h) : {eps_exciton:8.3f}   [Direct e-h attraction screening]")
        logger.info(f"      ε_eff (HOMO Hole self) : {eps_hole:8.3f}   [Hole self-energy screening]")
        logger.info(f"      ε_eff (LUMO Elec self) : {eps_elec:8.3f}   [Electron self-energy screening]")
        logger.info(f"      ε_eff (Inter-atomic)   : {eps_inter:8.3f}   [Asymptotic inter-atomic screening]")
        logger.info(f"      ε_eff (On-site core)   : {eps_onsite:8.3f}   [Atomic core screening (~1.0 = unscreened)]")
        logger.info(f"    --------------------------------------------------------------------------")
        if eps_bulk > 1.0:
            pct_bulk = (eps_eff_micro / eps_bulk) * 100.0
            logger.info(f"    Nanoscale Confinement & Size Comparison:")
            logger.info(f"      Material Bulk ε_∞      : {eps_bulk:8.3f}   [Experimental bulk limit]")
            logger.info(f"      Dielectric Retention   : {pct_bulk:7.1f}% of bulk (confinement suppression: {100.0 - pct_bulk:.1f}%)")
        logger.info(f"    ==========================================================================\n")

        eps_info = {
            "eps_eff_exciton": eps_exciton,
            "eps_eff_hole": eps_hole,
            "eps_eff_elec": eps_elec,
            "eps_interatomic": eps_inter,
            "eps_onsite": eps_onsite,
            "eps_bulk": eps_bulk,
            "cluster_radius_ang": r_qd_ang,
            "dielectric_retention_pct": (eps_eff_micro / eps_bulk * 100.0) if eps_bulk > 1.0 else 100.0,
            "alpha_cluster_bohr3": float(alpha_cluster_au),
            "kernel_mode": "rpa"
        }

        if return_eps_info:
            return w_rpa_ev, w_rpa_ev, gamma_bare_ev, eps_info
        return w_rpa_ev, w_rpa_ev, gamma_bare_ev

    elif km in ["dim", "dipole", "xs-dim", "xs_dim", "xs-dipole", "xs_dipole", "polarizable_dipole"]:
        S_atom, eta_atom, eps_inf_bulk, d_NN_ang = build_dim_screening_factors(
            coords=coords, atom_symbols=atom_symbols, material_name=m_name, alpha=alpha
        )

        s_ao = np.zeros((n_ao, n_ao), dtype=np.float64)
        for A, (a0, a1) in enumerate(atom_ao_ranges):
            for B, (b0, b1) in enumerate(atom_ao_ranges):
                s_ao[a0:a1, b0:b1] = S_atom[A, B]

        w_dim_ev = s_ao * gamma_bare_ev

        n_atoms = len(atom_symbols)
        inter_mask = ~np.eye(n_atoms, dtype=bool)
        eps_eff_median = float(1.0 / np.median(S_atom[inter_mask])) if np.any(inter_mask) else 1.0

        eps_exciton = eps_eff_median
        if C_occ_low is not None and C_virt_low is not None:
            q_h = np.abs(C_occ_low[:, -1])**2
            q_l = np.abs(C_virt_low[:, 0])**2
            v_eh_bare = float(q_h @ gamma_bare_ev @ q_l)
            w_eh_screened = float(q_h @ w_dim_ev @ q_l)
            if abs(w_eh_screened) > 1e-12:
                eps_exciton = v_eh_bare / w_eh_screened

        logger.info(f"\n    ==========================================================================")
        logger.info(f"    [Kernel: Xs-QDEX (Exact 2-Electron AO Integrals + Polarizable Dipole DIM)]")
        logger.info(f"    ==========================================================================")
        logger.info(f"    Material                = {m_name}")
        logger.info(f"    epsilon_in (bulk)       = {eps_inf_bulk:.3f}")
        logger.info(f"    Nearest neighbor d_NN   = {d_NN_ang:.3f} Å")
        logger.info(f"    Computed ε_eff (1S)     = {eps_exciton:.3f}")
        logger.info(f"    Median interatomic ε    = {eps_eff_median:.3f}")
        logger.info(f"    Dipole Matrix (3N)      = {3*n_atoms} x {3*n_atoms}")
        logger.info(f"    AO matrix dimension     = {n_ao} x {n_ao}")
        logger.info(f"    ==========================================================================\n")

        eps_info = {
            "eps_eff_exciton": eps_exciton,
            "eps_bulk": eps_inf_bulk,
            "eps_interatomic": eps_eff_median,
            "kernel_mode": "dim"
        }
        if return_eps_info:
            return w_dim_ev, w_dim_ev, gamma_bare_ev, eps_info
        return w_dim_ev, w_dim_ev, gamma_bare_ev

    elif km in ["resta", "xs_resta", "xs-resta"]:
        BOHR_TO_ANG = ANG_PER_BOHR
        entry = MATERIAL_DB.get(m_name, MATERIAL_DB["DEFAULT"])
        eps_inf_bulk = entry[0]

        core_coords = coords
        if m_name in MATERIAL_ELEMENTS:
            core_elements = [el.lower() for el in MATERIAL_ELEMENTS[m_name]]
            core_coords = np.array(
                [coords[i] for i, sym in enumerate(atom_symbols) if sym.lower() in core_elements]
            )

        if len(core_coords) > 1:
            r_core_ang = squareform(pdist(core_coords))
            np.fill_diagonal(r_core_ang, np.inf)
            d_NN_ang = np.median(np.min(r_core_ang, axis=1))
        else:
            d_NN_ang = 2.5

        d_NN_au = d_NN_ang / BOHR_TO_ANG
        k_s_au = np.sqrt(max(0.0, eps_inf_bulk - 1.0)) / d_NN_au
        lambda_s_ang = (1.0 / k_s_au) * BOHR_TO_ANG if k_s_au > 0 else 0.0

        coords_au = coords / BOHR_TO_ANG
        r_mat_au = squareform(pdist(coords_au))

        c_inf = 1.0 / eps_inf_bulk
        # Atom-atom Resta screening factor: S_AB = c_inf + (1 - c_inf)*exp(-k_s * r_AB)
        # On-site (A == B, r_AB == 0) gives exactly 1.0 (unscreened atomic core)
        s_atom = c_inf + (1.0 - c_inf) * np.exp(-k_s_au * r_mat_au)

        # Expand atom-atom screening factor to AO-AO matrix
        s_ao = np.zeros((n_ao, n_ao), dtype=np.float64)
        for A, (a0, a1) in enumerate(atom_ao_ranges):
            for B, (b0, b1) in enumerate(atom_ao_ranges):
                s_ao[a0:a1, b0:b1] = s_atom[A, B]

        w_resta_ev = s_ao * gamma_bare_ev
        logger.info(f"\n    [Kernel: Xs-QDEX (Exact 2-Electron AO Integrals + Resta Screening)]")
        logger.info(f"    Material            = {m_name}")
        logger.info(f"    epsilon_in          = {eps_inf_bulk:.3f}")
        logger.info(f"    Screening length    = {lambda_s_ang/BOHR_TO_ANG:.3f} a.u. ({lambda_s_ang:.3f} Å)")
        logger.info(f"    AO matrix dimension = {n_ao} x {n_ao}")

        eps_info = {
            "eps_eff_exciton": eps_inf_bulk,
            "eps_bulk": eps_inf_bulk,
            "kernel_mode": "resta"
        }
        if return_eps_info:
            return w_resta_ev, w_resta_ev, gamma_bare_ev, eps_info
        return w_resta_ev, w_resta_ev, gamma_bare_ev

    else:
        # Uniform screening
        w_bse_ev = alpha * gamma_bare_ev
        logger.info(f"\n    [Kernel: Xs-QDEX (Exact 2-Electron AO Integrals + Uniform Screening)]")
        logger.info(f"    alpha               = {alpha:.3f}")
        logger.info(f"    AO matrix dimension = {n_ao} x {n_ao}")

        eps_eff_uni = 1.0 / alpha if alpha > 0 else 1.0
        eps_info = {
            "eps_eff_exciton": eps_eff_uni,
            "kernel_mode": "uniform"
        }
        if return_eps_info:
            return w_bse_ev, w_bse_ev, gamma_bare_ev, eps_info
        return w_bse_ev, w_bse_ev, gamma_bare_ev


def estimate_brus_qp_gap(material_name, coords, atom_symbols):
    """Analytically estimates the Quantum Confined QP Gap using the Brus equation with Non-Parabolic corrections."""
    m_name = material_name.upper() if material_name else "DEFAULT"
    entry = MATERIAL_DB.get(m_name, MATERIAL_DB["DEFAULT"])
    
    E_bulk = entry[3]
    m_eff = entry[5]
    
    if E_bulk == 0.0 or m_eff == 0.0:
        logger.warning(f"  [Warning] Missing bulk gap or effective mass for {m_name}. Brus estimation failed.")
        return None
        
    R_QD_ang = qd_radius(coords, atom_symbols, material_name)
    
    # 1. Standard Parabolic Kinetic Confinement Energy
    E_conf_parabolic = (BRUS_KINETIC_EV_ANG2 * np.pi**2) / (m_eff * (R_QD_ang ** 2))
    
    # 2. Non-Parabolicity Correction (Hyperbolic Band Model)
    # Highly necessary for narrow gap materials like InAs, PbS, PbSe (E_bulk < 1.0 eV)
    if E_bulk < 2.0:
        predicted_gap = np.sqrt(E_bulk**2 + 2 * E_bulk * E_conf_parabolic)
        is_non_parabolic = True
    else:
        predicted_gap = E_bulk + E_conf_parabolic
        is_non_parabolic = False
        
    logger.info(f"\n  [Brus Model] Estimating Confinement for {m_name}:")
    logger.info(f"    Radius (R_QD)    : {R_QD_ang:.2f} Å")
    logger.info(f"    Bulk Gap         : {E_bulk:.3f} eV")
    logger.info(f"    Effective Mass   : {m_eff:.3f} m_e")
    
    if is_non_parabolic:
        logger.info(f"    Model Used       : Hyperbolic (Non-Parabolic)")
        logger.info(f"    Raw Parabolic dE : +{E_conf_parabolic:.3f} eV (Unphysical, applying correction...)")
    else:
        logger.info(f"    Model Used       : Standard Parabolic")
        logger.info(f"    Confinement (dE) : +{E_conf_parabolic:.3f} eV")
        
    logger.info(f"    Predicted QP Gap : {predicted_gap:.3f} eV")
    
    return predicted_gap
