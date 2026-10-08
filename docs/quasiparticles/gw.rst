GW quasiparticles in a quantum dot
==================================

Part of :doc:`/quasiparticles/index`.

This page starts from the full GW approximation, lists what each term is and how it is computed
in a standard GW code, and then shows how the quantum-dot problem reduces to the finite-size
correction :math:`\Delta W`. Excitons are not discussed here; they are the subject of
:doc:`/excitons/index`, which uses the same W.

1. Quasiparticle energies
-------------------------

A quasiparticle (QP) energy is the energy needed to add an electron to the system (empty levels)
or to remove one (occupied levels). It is the solution of the Dyson equation for the one-particle
Green's function. Expanding around Kohn–Sham (KS) orbitals :math:`\psi_n` with energies
:math:`\varepsilon_n`, and keeping only the diagonal:

.. math::

   \varepsilon_n^{\mathrm{QP}} = \varepsilon_n
   + \langle n|\,\Sigma(\varepsilon_n^{\mathrm{QP}}) - v_{xc}\,|n\rangle
   \;\approx\; \varepsilon_n + Z_n\,\langle n|\,\Sigma(\varepsilon_n) - v_{xc}\,|n\rangle .

The terms are:

* :math:`\varepsilon_n`, :math:`\psi_n`: the KS energies and orbitals of the starting point (PBE in
  QDEX).
* :math:`v_{xc}`: the KS exchange–correlation potential. It is removed because Σ replaces it.
* :math:`\Sigma(\mathbf r,\mathbf r';\omega)`: the self-energy. It is non-local and energy dependent,
  and contains all exchange and correlation beyond the Hartree term.
* :math:`Z_n = [1 - \partial_\omega\langle n|\Sigma(\omega)|n\rangle|_{\varepsilon_n}]^{-1}`: the
  quasiparticle weight. It linearizes the energy dependence of Σ and is between 0 and 1 (0.7–0.8 for
  bulk semiconductors).

2. The GW self-energy
---------------------

In the GW approximation (Hedin 1965) the self-energy is the product of the Green's function G and
the screened Coulomb interaction W:

.. math::

   \Sigma(\mathbf r,\mathbf r';\omega) = \frac{i}{2\pi}\int d\omega'\,
   G(\mathbf r,\mathbf r';\omega+\omega')\,W(\mathbf r,\mathbf r';\omega')\,e^{i\eta\omega'} .

**Green's function.** Built from the orbitals and energies,

.. math::

   G(\mathbf r,\mathbf r';\omega) = \sum_m \frac{\psi_m(\mathbf r)\psi_m^*(\mathbf r')}
   {\omega - \varepsilon_m + i\eta\,\mathrm{sgn}(\varepsilon_m - \mu)} .

**Screened interaction.** The bare Coulomb interaction :math:`v = 1/|\mathbf r-\mathbf r'|`
screened by the inverse dielectric function:

.. math::

   W(\mathbf r,\mathbf r';\omega) = \int d\mathbf r''\,\epsilon^{-1}(\mathbf r,\mathbf r'';\omega)\,v(\mathbf r'',\mathbf r'),
   \qquad
   \epsilon = 1 - v\,\chi_0 .

**Independent-particle response** (RPA):

.. math::

   \chi_0(\mathbf r,\mathbf r';\omega) = 2\sum_{i\in\mathrm{occ}}\sum_{a\in\mathrm{empty}}
   \psi_i(\mathbf r)\psi_a^*(\mathbf r)\psi_a(\mathbf r')\psi_i^*(\mathbf r')
   \left[\frac{1}{\omega - (\varepsilon_a-\varepsilon_i) + i\eta}
        - \frac{1}{\omega + (\varepsilon_a-\varepsilon_i) - i\eta}\right].

Everything that distinguishes a QP energy from a KS energy is in W: how strongly the rest of the
system screens an added charge. It is convenient to write :math:`W = v + W^{\mathrm{p}}`, where
:math:`W^{\mathrm{p}}` is the potential of the polarization induced by the charge.

3. How Σ is computed in a GW code
---------------------------------

A standard molecular or periodic GW calculation (for example evGW in CP2K) does the following:

1. **Response.** χ₀ is built from all occupied and many empty states, usually in an auxiliary
   (resolution-of-identity) basis of size :math:`N_{\mathrm{aux}}`.
2. **Screening.** ε is formed and inverted at each frequency.
3. **Frequency integral.** The ω′ integral is done with a plasmon-pole model, by contour
   deformation, or by analytic continuation from the imaginary axis.
4. **Self-consistency.**

   * G₀W₀: one shot on the KS orbitals;
   * evGW: the energies in G and W are iterated, the orbitals are kept;
   * qsGW: orbitals and energies are iterated to a static, Hermitian Σ.

The cost is :math:`\mathcal O(N^4)` in common implementations, with a large prefactor from the empty
states and the frequency grid. One calculation on a 3 nm CdSe dot (about 10³ atoms, 10⁴ basis
functions) is possible on a large machine. Size series, ligand shells, solvents and molecular
dynamics are not. The result also depends on the starting point (PBE, PBE0) and on the level of
self-consistency.

4. The static COHSEX limit
--------------------------

If the characteristic frequencies of W (plasmons, 10–20 eV) are large compared with the energies
of interest, W can be taken at ω = 0. The frequency integral then gives two static terms
(Hedin 1965):

.. math::

   \Sigma^{\mathrm{SEX}}(\mathbf r,\mathbf r') = -\sum_{m\in\mathrm{occ}}
   \psi_m(\mathbf r)\psi_m^*(\mathbf r')\,W(\mathbf r,\mathbf r';0),
   \qquad
   \Sigma^{\mathrm{COH}}(\mathbf r,\mathbf r') = \tfrac12\,\delta(\mathbf r-\mathbf r')\,
   W^{\mathrm{p}}(\mathbf r,\mathbf r;0).

Their diagonal elements are

.. math::

   \mathrm{SEX}_n = -\sum_{m\in\mathrm{occ}} (nm|W|mn),
   \qquad
   \mathrm{COH}_n = \tfrac12\int d\mathbf r\,|\psi_n(\mathbf r)|^2\,W^{\mathrm{p}}(\mathbf r,\mathbf r).

* **Screened exchange (SEX)** is the Fock exchange with the occupied states, with v replaced by W.
  It is non-local, it acts only through the occupied manifold, and it depends on how orbital n
  overlaps that manifold.
* **Coulomb hole (COH)** is the interaction of a charge with the polarization it induces at its
  own position. It is local and the same for occupied and empty states.

Static COHSEX overestimates bulk gaps by a few tenths of an eV because it neglects the dynamics of
W. It is used below only for the finite-size part of W, where the static limit is better justified
(section 5).

5. Quantum dots converge to the bulk: the ΔW split
--------------------------------------------------

A quantum dot is a piece of a crystal. Far from the surface, and for a large enough dot, its W
becomes the W of the bulk crystal. Following Delerue, Lannoo and Allan (PRL 84, 2457 (2000);
PRL 90, 076803 (2003)), write

.. math::

   W_{\mathrm{QD}}(\mathbf r,\mathbf r') = W_{\mathrm{bulk}}(\mathbf r,\mathbf r') + \Delta W(\mathbf r,\mathbf r'),
   \qquad \Delta W \to 0\quad (R\to\infty).

The self-energy splits the same way:

.. math::

   \langle n|\Sigma_{\mathrm{QD}} - v_{xc}|n\rangle =
   \underbrace{\langle n|\Sigma_{\mathrm{bulk}} - v_{xc}|n\rangle}_{\text{bulk GW correction}}
   + \underbrace{\Delta\Sigma_n[\Delta W]}_{\text{finite-size correction}} .

**Bulk part.** The GW correction of the crystal is transferable. QDEX takes it from bulk QSGW as
the gap opening

.. math::

   \Delta_{\mathrm{bulk}} = \Delta_\Sigma + \Delta_{\mathrm{geom}},\qquad
   \Delta_\Sigma = E_g^{\mathrm{QSGW+SOC}}(a_{\mathrm{exp}}) - E_g^{\mathrm{PBE+SOC}}(a_{\mathrm{exp}})
   \qquad (1.64\ \mathrm{eV\ for\ CdSe}),

with the literature QSGW gap (with spin-orbit coupling) and our PBE+SOC gap at the same, experimental
lattice (``MATERIAL_DB`` index 8 − index 7). :math:`\Delta_{\mathrm{geom}} = E_g^{\mathrm{PBE}}(a_{\mathrm{exp}}) -
E_g^{\mathrm{PBE}}(a_{\mathrm{dot}})` (``quasiparticles.bulk_geometry``) removes the gap change of a
PBE-relaxed dot, whose lattice is close to the PBE one (+0.166 eV for a 2.62 nm CdSe dot; GaAs +0.47,
CsPbBr\ :sub:`3` −0.39 eV at the full PBE lattice). Sources, lattices and the measurement of the dot
lattice: :doc:`/electronic_structure/bulk_bands`. :math:`\Delta_{\mathrm{bulk}}` is

split between the valence and conduction band: :math:`-f_b\,\Delta_{\mathrm{bulk}}` for occupied, :math:`+(1 - f_b)\,\Delta_{\mathrm{bulk}}` for empty states, where :math:`f_b` is the valence share of the bulk opening. The split does not change the gap; it places the levels on the absolute scale (IP, EA). It is not known
to better than about ±0.5 in nanocrystals: bulk GW with vertex corrections and experiment put 75–98 % of
the opening on the valence band of II-VI and III-V crystals, while evGW of 1–2 nm clusters puts almost
none of it there for the Zn, Cd, Ga and In compounds. QDEX therefore splits it evenly by default,
:math:`f_b = 1/2` (``quasiparticles.bulk_edge_split: symmetric``, as Biffi et al., arXiv:2210.01324);
``cluster`` and ``bulk`` select the two calibrations. Evidence and tables: :doc:`/workflows/qp_edges`. All dynamical effects of the bulk are inside this number.

**Vertex correction of the bulk part** (``quasiparticles.bulk_vertex``). QSGW overestimates bulk gaps
by 10–20 %: its W lacks the electron–hole (ladder) vertex, so the polarizability and ε are too small.
The usual bulk remedy scales the opening, :math:`\Delta_{\mathrm{bulk}} \to a\,\Delta_{\mathrm{bulk}}`
with :math:`a = 0.8` (``bulk_vertex_factor``); for CdSe this gives a bulk gap of 1.84 eV with SOC
(0.522 + 0.8 × 1.642), 0.10 eV above the experimental 1.74 eV. The vertex correction scales
:math:`\Delta_\Sigma` only, not :math:`\Delta_{\mathrm{geom}}`.

The error of QSGW is not the same for every material: against the room-temperature gaps it ranges from
+0.5 eV (ZnS, CdTe) through about zero (PbS, PbSe, AlSb) to −0.15 eV (CsPbBr\ :sub:`3` in its measured structure), and it contains more than the missing vertex (zero-point
and thermal renormalisation). ``bulk_vertex_factor: material`` therefore uses, for each material, the
factor that puts the bulk limit on its experimental gap,

.. math::

   a_m = \frac{E_g^{\mathrm{exp}} - E_g^{\mathrm{PBE+SOC}}}{\Delta_\Sigma}
   \qquad (\text{CdSe zb } 0.70,\ \text{wz } 0.73,\ \text{GaAs } 0.74,\ \text{InP } 0.71,\
   \text{PbS } 1.04,\ \text{CsPbBr}_3\ 1.11,\ \text{CsPbI}_3\ 0.98),

with :math:`E_g^{\mathrm{exp}}` at room temperature, the temperature of the sizing experiments, and the
PBE+SOC gap of the measured structure (for the lead halide perovskites the orthorhombic Pnma phases,
whose tilts open the gap by 0.45–0.70 eV with respect to the cubic cell of the same volume). With ``scaled`` the dot then goes
from QSGW in the molecular limit to the experimental gap in the bulk, the bulk limit of the g-xTB route.
A factor above 1 (Pb chalcogenides, perovskites) means QSGW is below the measured gap and is not a vertex
correction; the run says so. Values per material: :doc:`/electronic_structure/bulk_bands`. The missing vertex part is proportional to the
polarizability, which is reduced in a dot. With ``scaled`` the correction follows the fraction of bulk
screening the dot keeps,

.. math::

   \Delta_\Sigma(R) = \Delta_{\mathrm{QSGW}} - (1 - a)\,\Delta_{\mathrm{QSGW}}\,f(R),\qquad
   f(R) = \frac{\epsilon_{\mathrm{eff}}(R) - 1}{\epsilon_\infty - 1},

with :math:`\epsilon_{\mathrm{eff}}` the Penn value at the DFT gap of the cluster (the same Penn model
as the Resta W; the confinement energy is taken from the bulk PBE gap at the dot's lattice, of the
measured tilted structure for the perovskites, :doc:`/electronic_structure/bulk_bands`). f → 1 in the bulk (experimental gap), f → 0 in the molecular limit (QSGW). f uses
the DFT gap, so Δ_bulk is the same for every QP model of a given cluster; the screening of the ΔW
models is referenced to the same corrected bulk gap, so the vertex correction changes the bulk shift
and nothing else.

.. list-table::
   :header-rows: 1

   * - ``bulk_vertex``
     - Δ_Σ
     - CdSe: 1.2 nm / 2.0 nm / bulk (eV)
   * - ``none`` (default)
     - pure QSGW
     - 1.64 / 1.64 / 1.64
   * - ``full``
     - :math:`a\,\Delta_{\mathrm{QSGW}}` at every size
     - 1.31 / 1.31 / 1.31
   * - ``scaled``
     - :math:`\Delta_{\mathrm{QSGW}}\,[1 - (1-a) f(R)]`
     - ≈ 1.45 / 1.39 / 1.31 (2.62 nm: 1.38)

**Split model** (``quasiparticles.bulk_residual: experimental``). The residual of a fixed factor,
:math:`\delta_{\mathrm{res}} = E_g^{\mathrm{exp}} - (E_g^{\mathrm{PBE+SOC}} + a\,\Delta_\Sigma)`, is not a
vertex effect: for most II–VI and III–V semiconductors it is −0.08 to −0.23 eV (CdSe −0.16 eV, about
the thermal narrowing to room temperature plus the zero-point renormalisation), for the Pb compounds
+0.1 to +0.4 eV, where QSGW is already below the measured gap. The split model keeps the vertex part
size-dependent and adds the residual at every size,

.. math::

   \Delta_{\mathrm{bulk}}(R) = \Delta_\Sigma\,[1 - (1 - a) f(R)] + \delta_{\mathrm{res}} + \Delta_{\mathrm{geom}},

so the bulk limit is the measured gap, like ``bulk_vertex_factor: material``, but the temperature and
zero-point parts do not fade in small dots (with the material factor everything is scaled by
:math:`f(R)` as if it were vertex). For the 2.62 nm CdSe dot the two give shifts 0.03 eV apart; for
small dots the difference grows to about 0.07 eV. Recommended settings: ``bulk_vertex: scaled``,
``bulk_vertex_factor: 0.8``, ``bulk_residual: experimental``. The residuals per material are in
:doc:`/electronic_structure/bulk_bands`; a large positive one (CsPbBr\ :sub:`3` +0.42 eV) says that the 0.8
vertex scaling does not suit that material.

It applies to ``bulk`` (the sBSE), the ΔW models and the NAMD precompute, not to the legacy ``gw``
and ``sgw`` models.

**Finite-size part.** ΔW is computed for each dot and inserted in the static COHSEX form:

.. math::

   \Delta\Sigma_n = \underbrace{\tfrac12\int d\mathbf r\,|\psi_n(\mathbf r)|^2\,\Delta W(\mathbf r,\mathbf r)}_{\Delta\mathrm{COH}_n}
   \;\underbrace{-\sum_{m\in\mathrm{occ}}\iint d\mathbf r\,d\mathbf r'\,
   \psi_n^*(\mathbf r)\psi_m(\mathbf r)\,\Delta W(\mathbf r,\mathbf r')\,\psi_m^*(\mathbf r')\psi_n(\mathbf r')}_{\Delta\mathrm{SEX}_n} .

Written with the density matrix :math:`P(\mathbf r,\mathbf r') = 2\sum_{m\in\mathrm{occ}}\psi_m(\mathbf r)\psi_m^*(\mathbf r')`,

.. math::

   \Delta\mathrm{SEX}_n = -\tfrac12\iint d\mathbf r\,d\mathbf r'\,
   \psi_n^*(\mathbf r)\,P(\mathbf r,\mathbf r')\,\Delta W(\mathbf r,\mathbf r')\,\psi_n(\mathbf r') .

The QP energies of all orbitals are then

.. math::

   \varepsilon_n^{\mathrm{QP}} = \begin{cases}
   \varepsilon_n - f_b\,\Delta_{\mathrm{bulk}} + Z_n\,\Delta\Sigma_n, & n \in \mathrm{occ}, \\
   \varepsilon_n + (1 - f_b)\,\Delta_{\mathrm{bulk}} + Z_n\,\Delta\Sigma_n, & n \in \mathrm{virt},
   \end{cases}

where :math:`f_b` is the valence share of the bulk opening.

**Why the static limit is good for ΔW.** ΔW is dominated by the long-range polarization of the
dot surface and of the environment. Its characteristic frequency is the valence plasmon of the
dot (about 15 eV for CdSe), far above the band-edge level spacings. The remaining energy dependence
is kept to first order by :math:`Z_n` (:doc:`dynamic_z`).

**What ΔW contains.**

1. **Surface polarization.** The dot has ε∞ ≈ 6 inside and a smaller ε_out outside (1 in vacuum,
   2.24 in toluene). An added charge polarizes the interface, and the induced charge acts back on
   it.
2. **Reduced interior screening.** Confinement opens the gaps of the states that screen, so a small
   dot is less polarizable per volume than the bulk.
3. **Local fields at the surface.** Under-coordinated surface atoms and ligands screen differently
   from bulk atoms.
4. **The environment,** through ε_out.

6. The classical limit
----------------------

Take ΔW smooth on the scale of an atom, and let the occupied states act as a complete set for it,
:math:`\tfrac12 P(\mathbf r,\mathbf r') \approx \delta(\mathbf r-\mathbf r')` within the valence
manifold. Then

.. math::

   \Delta\mathrm{SEX}_n \to -\int|\psi_n|^2\,\Delta W(\mathbf r,\mathbf r)\ \ (n\ \mathrm{occupied}),
   \qquad \Delta\mathrm{SEX}_n \to 0\ \ (n\ \mathrm{empty}),

so that

.. math::

   \Delta\Sigma_{\mathrm{empty}} = +\tfrac12\langle\Delta W(\mathbf r,\mathbf r)\rangle_e,
   \qquad
   \Delta\Sigma_{\mathrm{occ}} = -\tfrac12\langle\Delta W(\mathbf r,\mathbf r)\rangle_h,
   \qquad
   \Delta E_g^{\mathrm{QP}} = \Sigma^{\mathrm{pol}}_e + \Sigma^{\mathrm{pol}}_h .

The gap opens by the **self-image** energies of the added electron and hole: the classical result
of Brus (J. Chem. Phys. 80, 4403 (1984)). A constant ΔW = c opens the gap by exactly c.

The completeness assumption fails in a small dot. The actual ΔSEX then depends on how each orbital
overlaps the occupied states. This **non-classical screened exchange** is where a microscopic ΔW
and the actual orbitals matter.

7. What a model has to supply
-----------------------------

With the bulk correction tabulated, a QP model for a quantum dot needs:

1. **ΔW(r, r′)**, or at least its classical self-image. This is where the models differ
   (:doc:`models`).
2. **Z**, from the same dielectric model (:doc:`dynamic_z`).
3. **A representation** of the integrals (nm|ΔW|mn) on the atomic basis: MNOK or ZDO xs
   (:doc:`/integrals/index`).

The same W, including ΔW, is the screened interaction of the BSE (:doc:`/excitons/kernel`).
