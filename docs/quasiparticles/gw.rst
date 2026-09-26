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

A standard molecular or periodic GW calculation (for example the CP2K evGW that gives the QDEX
anchor) does the following:

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
dynamics are not. The result also depends on the starting point: the QDEX anchor is evGW on PBE0
orbitals (PBE → PBE0 → evGW), not on PBE.

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

**Bulk part.** The GW correction of the crystal is transferable. QDEX takes it from a bulk GW
calculation as the gap opening

.. math::

   \Delta_{\mathrm{bulk}} = E_g^{\mathrm{GW}}(\mathrm{bulk}) - E_g^{\mathrm{PBE}}(\mathrm{bulk})
   \qquad (1.27\ \mathrm{eV\ for\ CdSe}),

split between the valence and conduction band: −½Δ_bulk for occupied, +½Δ_bulk for empty states.
All dynamical effects of the bulk are inside this number.

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

   \varepsilon_n^{\mathrm{QP}} = \varepsilon_n \mp \tfrac12\Delta_{\mathrm{bulk}} + Z_n\,\Delta\Sigma_n ,

with − for occupied and + for empty states.

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
3. **A calibration** of what the model misses, on one evGW calculation of a small cluster
   (:doc:`anchor`).
4. **A representation** of the integrals (nm|ΔW|mn) on the atomic basis: MNOK or ZDO xs
   (:doc:`representation`).

The same W, including ΔW, is the screened interaction of the BSE (:doc:`/excitons/screened_kernel`).
