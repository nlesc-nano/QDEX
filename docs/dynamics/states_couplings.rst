Excitons and couplings along the trajectory
===========================================

Part of :doc:`/dynamics/index`.

.. rubric:: QDEX implementation

* Modules: ``qdex.namd.precompute`` (excitons, overlaps, tracking), ``qdex.namd.nac`` (couplings)
* CLI: ``--namd-precompute``, then ``--namd-nac``
* YAML: ``excitations.*``, ``quasiparticles.*``, ``namd.tracking.*``, ``namd.nac.*``,
  ``namd.dynamics.nac_scheme``

For every MD frame the precompute builds the exciton states and, for every pair of consecutive frames,
the overlaps from which the couplings follow. This page derives the couplings, first in the one-electron
basis, then in the basis of excitons, and describes how the states are followed from frame to frame.

1. Exciton states of each frame
-------------------------------

The states are electron–hole pairs :math:`|ia\rangle = \hat a_a^\dagger \hat a_i |\Phi_0\rangle` of the
active window (``excitations.nhomos`` occupied and ``nlumos`` virtual orbitals; with SOC, spinors), with
the energies of the **diagonal** screened BSE (one energy per pair, no mixing between pairs), so that the
same transitions can be followed along the trajectory:

.. math::

   E_{ia}(t) = \varepsilon_a(t) - \varepsilon_i(t) + \Delta_{\mathrm{bulk}}
   + k_x\, K^x_{ia,ia}(t) - K^d_{ia,ia}(t),

with :math:`k_x` = 2 (spin-free singlets) or 1 (spinors), Mulliken charges and the MNOK representation
(:doc:`/excitons/sbse`). :math:`\Delta_{\mathrm{bulk}}` is the bulk QP correction of the material
(``quasiparticles.model: bulk``, :doc:`/quasiparticles/gw`), computed once on frame 0 and kept for all
frames: the energy fluctuations along the trajectory then come only from the Kohn–Sham levels and the
interaction terms, and the surface polarization is absent from both the orbital energies and the kernel.
:math:`K^d` (direct, screened electron–hole attraction) and :math:`K^x` (exchange) are rebuilt from every
frame's geometry and orbitals. :math:`K^x` is a few meV for the pairs of a 4 nm CsPbX\ :sub:`3` dot and
costs 2.7× the time of a frame; ``include_exchange: false`` leaves it out. With SOC, :math:`K^x`,
:math:`K^d` and :math:`|\mu|^2` inside a Kramers pair are averaged over the pair (section 4.4).

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Key
     - Default
     - Meaning
   * - ``quasiparticles.model``
     - ``bulk``
     - ``bulk`` (one constant correction for all frames), ``none``, ``brus``, ``gw`` or a gap in eV.
       Models that need a QP step per frame (``sgw-*``, ``evgw-*``, ``qsgw-*``) are rejected.
   * - ``quasiparticles.reference``
     - ``pbe``
     - orbitals the ``bulk`` shift corrects: ``pbe`` or ``gxtb`` (:doc:`gxtb`)
   * - ``excitations.mode``
     - ``diagonal_sbse``
     - ``diagonal_sbse`` / ``diagonal_bse`` (same solver), ``diagonal_stda``, ``independent_qp``,
       ``independent_dft``
   * - ``excitations.kernel``
     - ``resta``
     - ``resta`` (bulk Resta W), ``dim`` or ``bse``; not used by ``diagonal_stda``
   * - ``excitations.include_exchange``, ``include_direct_eh``
     - true
     - :math:`K^x` and :math:`K^d` in every frame
   * - ``excitations.nhomos``, ``nlumos``
     - —
     - active window; it must extend about 1 eV beyond the largest excess energy of the pump (section 5)

With ``diagonal_stda`` the energy is Grimme's diagonal sTDA element
:math:`\varepsilon_a - \varepsilon_i + 2(ia|ia)_K - (ii|aa)_J` (:doc:`/excitons/stda`) on the orbital
energies as they are (``quasiparticles.model: none``), the setting for orbitals of a hybrid functional or
of g-xTB (:doc:`gxtb`). ``gw`` with the bulk kernel is accepted for old runs with a warning: its QP gap
contains the surface polarization but the kernel lacks the matching electron–hole image, which places the
excitons too high.

2. Couplings in the one-electron basis
--------------------------------------

The orbitals of the active window depend on time through the geometry. Their coupling is

.. math::

   d_{pq}(t) = \langle \phi_p(t) | \partial_t\, \phi_q(t) \rangle , \qquad d = -d^\dagger ,

separately for the occupied (:math:`d^{\mathrm{occ}}_{ij}`) and the virtual (:math:`d^{\mathrm{virt}}_{ab}`)
orbitals; couplings between occupied and virtual orbitals connect :math:`|ia\rangle` to the ground state
and to double excitations, which are outside the exciton basis (no recombination through them).
:math:`d` is not differentiated analytically: it follows from the overlap of the orbitals at consecutive
frames (section 4).

3. Couplings in the exciton basis
---------------------------------

**Single pairs.** Differentiating :math:`|ia\rangle = \hat a_a^\dagger \hat a_i |\Phi_0\rangle` with the
orbitals moving gives, to first order in the orbital couplings,

.. math::

   D_{ia,jb} \equiv \langle ia | \partial_t\, jb \rangle
   = \delta_{ij}\, d^{\mathrm{virt}}_{ab} \;-\; \delta_{ab}\, d^{\mathrm{occ}}_{ji}
   = \delta_{ij}\, d^{\mathrm{virt}}_{ab} \;+\; \delta_{ab}\, \big(d^{\mathrm{occ}}_{ij}\big)^* ,

up to a diagonal phase term that the phase alignment removes. Two pairs are coupled only when they share
the hole (the electron moves, :math:`a \to b`) or the electron (the hole moves, :math:`i \to j`): one
carrier changes state at a time. The hole sees the complex conjugate of the orbital coupling, because a
hole in :math:`i` is the absence of an electron there; for real orbitals
:math:`(d^{\mathrm{occ}}_{ij})^* = d^{\mathrm{occ}}_{ij}`, for spinors it is not, and the dynamics uses
:math:`d^{\mathrm{occ}*}` in the hole channel.

**Correlated excitons (full basis).** An exciton of the full BSE (or TDA, CIS) is a superposition
:math:`|A\rangle = \sum_{ia} X^A_{ia} |ia\rangle`. Its coupling to :math:`|B\rangle` has two parts,

.. math::

   \langle A | \partial_t B \rangle
   = \sum_{i,ab} X^{A*}_{ia} X^B_{ib}\, d^{\mathrm{virt}}_{ab}
   + \sum_{ij,a} X^{A*}_{ia} X^B_{ja}\, \big(d^{\mathrm{occ}}_{ij}\big)^*
   + \sum_{ia} X^{A*}_{ia}\, \dot X^B_{ia} ,

the orbital couplings dressed by the exciton coefficients, and the change of the coefficients themselves
along the trajectory. In practice it is obtained, like the orbital couplings, from the overlap of the
excitons of consecutive frames, :math:`\langle A(t) | B(t+\Delta t) \rangle`, which needs the overlaps of
the configurations (determinant overlaps, or their one-electron approximation
:math:`\langle ia(t) | jb(t+\Delta t)\rangle \approx S^{\mathrm{occ}*}_{ij} S^{\mathrm{virt}}_{ab}` with the
other occupied orbitals following themselves) and a tracking of the excitons, whose ordering and mixing
change from frame to frame in a dense manifold.

**What QDEX uses.** In the diagonal sBSE every exciton is a single pair, :math:`X = \mathbb{1}`, and the
coupling is the single-pair :math:`D` above. Its structure, couplings only through one carrier, means that
the exciton amplitudes :math:`c_{ia}` never need the :math:`N_{\mathrm{pairs}}^2` coupling matrix:

.. math::

   \dot c_{ia} = -\frac{i}{\hbar} E_{ia}\, c_{ia}
   - \sum_b d^{\mathrm{virt}}_{ab}\, c_{ib}
   - \sum_j \big(d^{\mathrm{occ}}_{ij}\big)^* c_{ja} ,

an **electron channel** (:math:`n_{\mathrm{virt}} \times n_{\mathrm{virt}}`, at the hole of the pair) and a
**hole channel** (:math:`n_{\mathrm{occ}} \times n_{\mathrm{occ}}`, at the electron). Surface hopping
propagates the two channels of each trajectory with the pair energies at the active partner carrier, the
master equation propagates the population matrix :math:`P_{ia}` with electron and hole rates
(:doc:`fssh`, :doc:`pme`). For the CsPbX\ :sub:`3` dots the window of 1300 × 800 orbitals (2600 × 1600
spinors) holds :math:`10^6` (:math:`4 \times 10^6`) pairs. The full-basis coupling, with the mixing of
the excitons, is not used in the dynamics: the diagonal approximation keeps one label per transition,
which can be followed along the trajectory, while mixed excitons would have to be re-identified at every
frame in a manifold of :math:`10^6` states and the coupling matrix would no longer factorise.

4. Overlaps and tracking of the states
--------------------------------------

For every pair of consecutive frames :math:`t_k, t_{k+1} = t_k + \Delta t` (the MD time step,
``namd.trajectory.dt_nuc_fs``) the precompute performs, in this order:

1. overlaps of the occupied and of the virtual active states;
2. Hungarian relabelling of trivial crossings;
3. phase alignment;
4. SOC only: parallel transport inside each Kramers pair;
5. energies, :math:`K^d` and :math:`K^x` permuted with the states; pair energies rebuilt;

and stores ``step_<k>_to_<k+1>.npz``. Frame :math:`t_{k+1}`, aligned, is the reference of the next step:
labels and phases are carried from frame 0 through the whole trajectory.

4.1 Cross-frame overlaps
~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   S_{ij}(t, t+\Delta t) = \langle \phi_i(t) | \phi_j(t+\Delta t) \rangle
   = \sum_{\mu\nu} C_{\mu i}(t)\, S^{\mathrm{AO}}_{\mu\nu}(t, t+\Delta t)\, C_{\nu j}(t+\Delta t),
   \qquad
   S^{\mathrm{AO}}_{\mu\nu} = \int \chi_\mu(\mathbf{r}; \mathbf{R}(t))\, \chi_\nu(\mathbf{r}; \mathbf{R}(t+\Delta t))\, d\mathbf{r},

with the basis functions of both geometries (Libint2). For spinors
:math:`\psi_j = \sum_n (U^\alpha_{nj}\, \phi_n \alpha + U^\beta_{nj}\, \phi_n \beta)` over the active MOs,

.. math::

   S^{\mathrm{spinor}} = U^{\alpha\dagger}(t)\, M\, U^\alpha(t+\Delta t) + U^{\beta\dagger}(t)\, M\, U^\beta(t+\Delta t),
   \qquad M = C^{\mathrm{act}}(t)^\mathsf{T} S^{\mathrm{AO}} C^{\mathrm{act}}(t+\Delta t).

4.2 Hungarian matching of trivial crossings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A trivial crossing is a label swap between two states that do not interact: within one step the diagonal
overlap collapses and the character moves to an off-diagonal element; following the energy label would
put a spike :math:`|S_{ij}|/\Delta t` into the coupling. A state is **locked** to its label when
:math:`|S_{ii}| \ge 0.5` (avoided crossings stay in the adiabatic basis that the dynamics propagates); the
unlocked states are matched by the Hungarian algorithm (``scipy.optimize.linear_sum_assignment``),

.. math::

   \pi = \arg\min_{\pi} \sum_{i\ \mathrm{unlocked}} \left(1 - |S_{i,\pi(i)}|^2\right),

over the unlocked block only (a locked state keeps its column, so no other state can take it). The
permutation is applied to the columns of :math:`S`, to the coefficients, energies, :math:`K^d`,
:math:`K^x` and dipoles of frame :math:`t+\Delta t`. The number of states whose label changed energy rank
in a step is stored (``n_swap_occ``, ``n_swap_virt``) and printed. The assignment is not a global
diabatisation.

4.3 Phase alignment
~~~~~~~~~~~~~~~~~~~

Eigenvectors come with an arbitrary phase; after the matching each state of :math:`t+\Delta t` is
multiplied by :math:`\operatorname{sign} S_{jj}` (real orbitals) or :math:`e^{-i \arg S_{jj}}` (spinors), so
that :math:`S_{jj} \ge 0`.

4.4 Kramers pairs of the SOC spinors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A Kramers pair is exactly degenerate and comes out of the diagonaliser in an arbitrary SU(2) frame, a
different one at every frame. A U(1) phase per spinor cannot undo that rotation, which appears as a
coupling :math:`|S_{12}|/\Delta t` between the partners. Each degenerate group :math:`G` (energies within
``tracking.degeneracy_tol_ev``, default 10\ :sup:`−5` eV) is parallel-transported: from the SVD
:math:`S_{GG} = W \Sigma V^\dagger` its columns are rotated by :math:`R = V W^\dagger`, so that
:math:`S_{GG} R = W \Sigma W^\dagger` is Hermitian positive semidefinite. Pair-diagonal quantities inside a
group (:math:`K^x`, :math:`K^d`, :math:`|\mu|^2`) depend on the frame of the pair; only their trace does
not, so they are replaced by the group average.

On 20 frames of the CsPbBr\ :sub:`3` and CsPbCl\ :sub:`3` cubes (400 × 200 spinors) the U(1)-only
alignment left a coupling between Kramers partners of 0.21 and 0.18 fs\ :sup:`−1` (median; up to
0.43 fs\ :sup:`−1`, :math:`\hbar d` ≈ 0.14 eV) in every step; with the parallel transport it is zero, the
other couplings are unchanged and the median :math:`|S_{jj}|` rises from 0.76 to 0.92 (Br) and from 0.66
to 0.87 (Cl). SOC precomputes made before this correction (2026-10-07) carry the spurious couplings.

4.5 The couplings from the overlaps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The coupling is taken constant over the step (the integrators hold it fixed between frames,
:doc:`propagation`). Two schemes (``dynamics.nac_scheme``):

* ``hst`` (Hammes-Schiffer & Tully, J. Chem. Phys. 101, 4657 (1994)), first order in the rotation:

  .. math::

     d \approx \frac{S - S^\dagger}{2\Delta t};

  a pair rotating by an angle :math:`\theta` within the step gets :math:`\sin\theta/\Delta t` instead of
  :math:`\theta/\Delta t`;

* ``logm`` (default): :math:`S` is made unitary (Löwdin; this also absorbs the norm that leaks out of the
  active window at its edges) and the coupling is the generator of that unitary,

  .. math::

     U = S\,(S^\dagger S)^{-1/2} = W V^\dagger \quad (S = W \Sigma V^\dagger), \qquad
     d = \frac{\log U}{\Delta t}.

  Propagating the amplitudes with :math:`e^{-d\,\Delta t}`, as the integrators do, reproduces the overlap
  exactly, :math:`c(t+\Delta t) = U^\dagger c(t)`. :math:`\log U` is evaluated from the eigenvectors of the
  Hermitian :math:`(U - U^\dagger)/2i` (which :math:`U` shares), with a residual check and
  ``scipy.linalg.logm`` as fallback for rotations beyond 90°; it is real for real orbitals and
  anti-Hermitian for spinors.

In the valence band of the CsPbX\ :sub:`3` dots (1000 occupied orbitals over 2 eV, ~5 meV apart) states
rotate by up to ~50° in a 2 fs step: the logm coupling of the hole states is 16 % (CsPbBr\ :sub:`3`) and
37 % (CsPbCl\ :sub:`3`) larger than HST (median per state; 49 % deep in the band of CsPbCl\ :sub:`3`), so
HST rates (:math:`\propto |d|^2`) are 35–120 % too low. At the band edges and for the electrons the two
agree within 2–6 %.

The logm couplings cost an SVD and an eigendecomposition per step (≈ 1–4 s for 1300 real, 10–35 s for
2600 complex states), so ``--namd-nac`` computes them once after the precompute, in parallel over the
steps (``namd.nac.workers`` × ``threads_per_worker``; restartable), and stores ``nac_<k>_to_<k+1>.npz``
next to the step files. Without them the dynamics computes the chosen scheme on the fly.

5. The active window
--------------------

States at the edges of the window lose norm to states outside it (:math:`1 - \sum_j |S_{ij}|^2` up to
0.1–0.8 for the outermost few states, below 10\ :sup:`−3` near the gap); the warning in the log refers to
them. The window must reach about 1 eV beyond the largest excess energy the pump gives to the electron or
the hole. For the 4 nm CsPbX\ :sub:`3` cubes, 1000 occupied orbitals cover only 2.0–2.3 eV of the valence
band and 350 virtual orbitals 2.1–2.4 eV of the conduction band, so a pump at 2 E\ :sub:`g` puts carriers
at the window edges; 1300 × 800 covers 3.1–3.35 eV (the whole halide-p band) and 3.35–3.65 eV.
