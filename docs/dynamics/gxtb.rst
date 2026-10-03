Dynamics on g-xTB trajectories
==============================

Part of :doc:`/dynamics/index`.

.. rubric:: QDEX implementation

Implementation entry points:

* Modules: ``qdex.xtb.calculator``, ``qdex.xtb.md``, ``qdex.xtb.molden``, ``qdex.namd.precompute``
* Callables: ``qdex.xtb.md.main``, ``qdex.xtb.molden.convert_molden``, ``qdex.namd.precompute.precompute_namd_data``
* CLI: ``python -m qdex.xtb.md``, ``--namd-precompute``
* YAML: ``namd.trajectory.basis_file``, ``excitations.mode: diagonal_stda``, ``excitations.functional: gxtb``,
  ``quasiparticles.reference: gxtb``, ``system.basis_name: per-atom``

The standard NAMD workflow of QDEX (:doc:`pipeline`) starts from an *ab initio* MD trajectory with DFT
orbitals at every frame. Both halves are expensive: every 2 fs step needs a converged DFT calculation
for the forces and another one, with a good basis, for the orbitals. For a nanocrystal of a few
hundred atoms this limits trajectories to a few picoseconds.

This page describes a route in which both the nuclear dynamics and the orbitals come from the
tight-binding method **g-xTB** (Froitzheim, Müller, Hansen, Grimme, ChemRxiv 2025). g-xTB is fitted to
ωB97M-V/def2-TZVPPD, covers the elements H–Lr, and has analytic gradients. For the 149-atom
Cd₆₈Se₅₅Cl₂₆ dot (2 nm) one energy and gradient takes about 2 s on a laptop, so a 2 ps trajectory at
2 fs with orbitals at every frame takes about 1.5 h.

The price is that g-xTB orbital energies are not PBE orbital energies. The quasiparticle corrections
and kernels of QDEX were built for PBE, so the excited-state part of the route has to be reconsidered.
The page explains the choices, shows what was validated, and lists what did not work.


The route at a glance
---------------------

.. code-block:: text

   start.xyz (relaxed structure)
        │
        ▼
   [1] ASE molecular dynamics, forces from g-xTB         qdex.xtb.md
        ├─ Maxwell–Boltzmann velocities, Langevin equilibration
        └─ NVE velocity Verlet production; each step also writes the molden file
        │
        ▼
   [2] Orbitals of each frame                             qdex.xtb.molden
        ├─ per-atom basis (the g-xTB basis depends on the atomic charge)
        └─ frame_XXXXXX/{frame.xyz, BASIS_GXTB, MOs.mbse}
        │
        ▼
   [3] NAMD precompute                                    qdex.namd.precompute
        ├─ cross-frame overlaps with each frame's own basis
        └─ excited states per frame: diagonal_stda (g-xTB as ωB97M-V)
                                     or diagonal_sbse + g-xTB bulk shift
        │
        ▼
   [4] Dynamics and analysis, unchanged (PME, DISH, CPA-FSSH; see the other pages of this section)


1. Setting up g-xTB
-------------------

At present g-xTB is available as the ``--gxtb`` method of the development ("bleed") build of xtb,
distributed by the g-xTB repository (https://github.com/grimme-lab/g-xtb). It is not yet part of a
released ``tblite``, so it has no Python API. QDEX therefore drives the command-line program.

The calculator reads the location of the program and of its runtime libraries from the environment.
A small file sourced before the runs is convenient:

.. code-block:: bash

   # env.sh
   export GXTB=/path/to/xtb-bleed-macos-arm64/bin/xtb      # binary with --gxtb
   export XTB_LIBDIR=$HOME/miniforge3/envs/xtb/lib         # runtime libraries (libgfortran, libgcc_s)
   export XTBPATH=$HOME/miniforge3/envs/xtb/share/xtb
   export OMP_STACKSIZE=1G

Two practical notes for macOS:

* The bleed binary needs ``libgcc_s`` from a conda ``xtb`` environment. QDEX puts ``XTB_LIBDIR`` on
  the dynamic-library path of the child process. System binaries such as ``/usr/bin/time`` or
  ``nohup`` strip ``DYLD_*`` variables, so do not wrap the binary in them.
* ``OMP_STACKSIZE=2G`` with 8 threads made xtb crash at start-up in some runs. QDEX uses 1 G.

The MD driver uses ASE (``pip install ase``); it is an optional dependency of QDEX.


2. Molecular dynamics with ASE
------------------------------

**Why ASE.** xtb has its own MD driver (``--md``), but with ``--gxtb`` it stops at the first step
("MD is unstable"). ASE provides tested integrators and thermostats, and QDEX adds a thin calculator,
``qdex.xtb.calculator.GXTB``, that runs ``xtb geom.xyz --gxtb --grad`` and reads energy and forces
from the ``.engrad`` file.

**The protocol** follows the usual practice for NAMD in the classical path approximation:

1. Velocities from a Maxwell–Boltzmann distribution at the target temperature, with centre-of-mass
   translation and rotation removed.
2. Equilibration with a Langevin thermostat (friction 0.01 fs⁻¹ by default, a time constant of
   100 fs).
3. Production in the microcanonical ensemble with velocity Verlet. The trajectory is not perturbed
   by a thermostat, so the energy fluctuations that drive decoherence and couplings are physical.

During production, every force calculation also writes the molden file of the new geometry. It is
converted at once into a QDEX frame (section 3) and then deleted. Forces and orbitals of a frame
therefore come from the same g-xTB calculation, and no second pass over the trajectory is needed.

.. code-block:: bash

   source env.sh
   python -m qdex.xtb.md start.xyz --outdir run --temp 300 --dt 2 \
          --n-equil 500 --n-prod 1000 --nthreads 8

Options: ``--friction`` (fs⁻¹), ``--seed``, ``--charge``, ``--acc`` and ``--etemp`` (passed to xtb),
``--binary`` and ``--libdir`` (instead of the environment), ``--keep-molden``. The start file may be
an ASE ``.traj`` with velocities, which continues a run without a new equilibration.

**Output** in ``run/``:

* ``frames/frame_000001 ...``: the QDEX frames (section 3), plus ``velocities.npy``;
* ``equil_log.csv``, ``prod_log.csv``: time, potential, kinetic and total energy, temperature, wall time
  per step;
* ``equil.traj``, ``prod.traj``, ``equil_final.traj``, ``prod_final.traj``: ASE trajectories;
* ``calc/``: the working directory of xtb. Its ``xtbrestart`` file is reused as the SCF guess of the
  next step, which reduces the SCF to about 5 cycles.

**Test run** (Cd₆₈Se₅₅Cl₂₆, 1 ps equilibration and 50 production frames at 300 K, 8 threads):

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Quantity
     - Value
   * - wall time per production step (forces, molden, conversion)
     - 2.7 s
   * - mean temperature, last 250 fs of equilibration / production
     - 304 K / 302 K
   * - total-energy drift in NVE
     - 6 meV over 98 fs (std 2 meV)
   * - g-xTB gap along the trajectory
     - 4.65 ± 0.09 eV (4.86 eV at the minimum)

For longer runs the drift can be reduced with a tighter SCF (``--acc 0.1``).

**Custom workflows.** The calculator can be used directly with any ASE dynamics or optimizer:

.. code-block:: python

   from ase.io import read
   from qdex.xtb.calculator import GXTB

   atoms = read("start.xyz")
   atoms.calc = GXTB(directory="calc", nthreads=8)
   forces = atoms.get_forces()               # eV/Å
   atoms.calc.molden_dest = "frame/molden.input"
   atoms.calc.reset(); atoms.get_forces()    # this calculation also keeps the molden file


3. Orbitals of each frame
-------------------------

**A basis that changes with the charges.** The g-xTB basis (q-vSZP) is a minimal valence basis whose
contraction coefficients depend on the effective charge of each atom. In the 2 nm dot, for example,
the leading Cd s coefficient is 0.1677, 0.1705, 0.1725 and 0.1692 on four different Cd atoms of the
same frame. The basis therefore differs from atom to atom and from frame to frame. Every QDEX frame
carries its own basis, and the cross-frame overlap S(t, t + Δt) is computed between two different
basis sets.

**Conversion.** ``qdex.xtb.molden.convert_molden`` turns the molden file into the QDEX conventions
(spherical functions in libint order, unit-normalized contractions):

1. The contraction coefficients in xtb molden files refer to unnormalized primitives. They are
   converted to normalized primitives, the convention libint expects.
2. xtb writes Cartesian d functions. The two xtb writers normalize them differently: the g-xTB
   (tblite) writer gives all components the norm of the axis component, the classic xtb writer (GFN2)
   normalizes every component. The converter tries both conventions and keeps the one for which the
   orbitals are orthonormal.
3. The orbitals are projected onto spherical functions with libint overlaps,

   .. math::

      C_{\mathrm{sph}} = S_{ss}^{-1}\, S_{sc}\, C_{\mathrm{cart}},

   which is exact for orbitals without s/p contamination of the Cartesian shells.
4. A spin-unrestricted molden file of a closed shell (identical α and β sets) is merged into one
   restricted set.

The acceptance test is orthonormality, :math:`\max|C^{\mathrm{T}} S C - 1|`, in the libint metric (default tolerance 10⁻⁶).

**Files of a frame.**

* ``frame.xyz``: the geometry;
* ``BASIS_GXTB``: CP2K-style basis file with one block per atom, named ``GXTB-A<index>``:

  .. code-block:: text

     Cd GXTB-A1
       2
       1 0 0 5 1
           1.6110073208000000e+00   1.6456352095383311e-01
           ...

* ``MOs.mbse``: all orbitals, energies (hartree) and occupations, in the QDEX binary format;
* ``conversion.json``: detected conventions, orthonormality error, gap.

A frame of the 2 nm dot (1001 basis functions) takes 7.8 MB.

**Validation** on the 2 nm dot: orthonormality ≤ 1.2·10⁻⁹ on all 50 frames; 648 electrons; Mulliken
charges from the converted orbitals identical to those printed by xtb (Cd +0.673, Se −0.589,
Cl −0.515).

The converter also reads molden files of GFN2-xTB and of xtb4stda (sTDA-xTB). It can be used on its
own:

.. code-block:: python

   from qdex.xtb.molden import convert_molden
   info = convert_molden("molden.input", "frame_000001", nthreads=8)
   print(info["gap_ev"], info["orthonormality_error"])


4. Precomputing the couplings
-----------------------------

The NAMD precompute (:doc:`configuration`) needs one extra key, the name of the per-frame basis file:

.. code-block:: yaml

   system:
     material: CDSE
     nthreads: 8

   namd:
     trajectory:
       dir: run/frames
       frame_pattern: frame_*
       xyz_file: frame.xyz
       mo_file: MOs.mbse
       basis_file: BASIS_GXTB      # per-frame g-xTB basis; replaces system.basis_txt
       dt_nuc_fs: 2.0
     storage:
       precompute_dir: namd_precomputed

With ``basis_file`` set, the shells of every frame are built from that frame's file, and the
cross-frame overlap uses the bases of both frames. Phase alignment and crossing tracking are those of
:doc:`nacs_tracking`.

**How well consecutive frames overlap.** On the 50-frame test the change of basis between frames
loses at most 1.7·10⁻³ of the orbital norm over the full orbital space. The band-edge orbitals
(HOMO−4 … LUMO+4) lose less than 2·10⁻⁴ inside the active window. The precompute does report larger
"active-space norm loss" values (up to 0.86 with a 40 + 40 window), but these come only from the
orbitals at the window edges (HOMO−40, LUMO+39), which mix with orbitals outside the window. They do
not affect the band-edge dynamics; a larger window moves them further away.


5. Excited states of each frame
-------------------------------

This is where the g-xTB route differs most from the PBE route, so it is worth looking at what the
g-xTB orbital energies are before choosing a model.

**What the g-xTB gap is.** g-xTB reproduces ωB97M-V, a range-separated hybrid with 100 %
long-range exact exchange. Its HOMO–LUMO gap behaves like the generalized Kohn–Sham gap of such a
functional:

* For a small cluster it is close to the fundamental gap. For Cd₁₆Se₁₃Cl₆ (1.2 nm) g-xTB gives
  5.96 eV and evGW\@PBE0 6.03 eV (:doc:`/validation/evgw_cluster`).
* For a larger dot it is too large, because the long-range exchange is not screened by the dot's
  polarizability. For the 2 nm dot g-xTB gives 4.47 eV at the PBE geometry, against a QP gap of about
  2.8–3.0 eV from the QDEX models.
* Absolute levels are about 1.8 eV too deep (HOMO −9.6 eV against −7.8 eV for evGW at 1.2 nm); only
  differences are meaningful.

So the g-xTB gap is neither a PBE gap, to which the bulk GW correction of QDEX could be added, nor a
quasiparticle gap. QDEX offers two consistent ways to handle it.


5.1 sTDA with g-xTB as the hybrid functional (``diagonal_stda``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first way takes g-xTB at its word: if its orbitals behave like those of ωB97M-V, the excitations
can be computed as Grimme's sTDA computes them for ωB97M-V orbitals (:doc:`/excitons/stda`).

.. code-block:: yaml

   quasiparticles:
     model: none               # g-xTB orbital energies as they are

   excitations:
     mode: diagonal_stda       # stda in a single-point run
     functional: gxtb          # Grimme's omegaB97M-V set: a_x 0.51, alpha 4.51, beta 8.0
     nhomos: 40
     nlumos: 40

Each transition energy is

.. math::

   E_{ia}(t) = \varepsilon_a(t) - \varepsilon_i(t) + 2\,(ia|ia)_K - (ii|aa)_J ,

with Löwdin transition charges and Grimme's interactions γ\ :sup:`K` (exchange) and γ\ :sup:`J`
(direct), rebuilt from every frame's geometry. ``functional: gxtb`` selects the parameters Grimme
fitted for ωB97M-V (stda program, ``-wB97MV``). They can be changed with ``ax``, ``stda_alpha`` and
``stda_beta``.

**Why this is consistent.** The electron–hole attraction of sTDA tends to the bare 1/R at long range:
like the g-xTB gap, it contains no screening by the dot. An unscreened gap and an unscreened
attraction belong together, as in TDHF or range-separated TDDFT, and their errors largely cancel. No
quasiparticle model and no material parameter (ε∞, bulk gaps) enter.

**Results** (2 nm dot):

* along the 50-frame trajectory at 300 K, the lowest transition is at 2.99 ± 0.09 eV;
* the experimental window for this size is 2.77–3.31 eV (:doc:`/validation/cdse_experiment`).

Part of this agreement comes from the diagonal approximation. The full sTDA (``mode: stda``) on the
same orbitals gives S₁ = 2.54 eV at the PBE geometry, against 2.86 eV for ``diagonal_stda``: the
coupling between transitions lowers S₁ by 0.32 eV at 2 nm (0.46 eV at 1.2 nm). The energy spacings
between transitions, which drive the cooling, are less affected than the absolute energy of S₁.

**Open point.** Whether the cancellation holds for larger dots has not been tested. The screening of a
real dot grows with size while the sTDA kernel keeps its bare tail; the overestimate of the g-xTB gap
also grows. Dots of 3–4 nm (700–1300 atoms) would answer the question.


5.2 sBSE with a g-xTB bulk shift (``quasiparticles.reference: gxtb``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The second way keeps the sBSE of the PBE route (:doc:`/excitons/sbse`): a constant shift of the
orbital energies plus the bulk Resta kernel. Only the shift changes, because it must now correct
g-xTB instead of PBE.

.. code-block:: yaml

   quasiparticles:
     model: bulk
     reference: gxtb           # g-xTB bulk shift instead of bulk QSGW - bulk PBE

   excitations:
     mode: diagonal_sbse
     kernel: resta

**The shift.** As for PBE, the shift is the difference between a reference bulk gap and the bulk gap
of the method. The reference is the spin-free experimental gap (g-xTB has no SOC):

.. math::

   \Delta_{\mathrm{bulk}}^{\mathrm{gxtb}} = E_g^{\mathrm{ref,SF}} - E_g^{\mathrm{gxtb,bulk}},
   \qquad E_g^{\mathrm{ref,SF}} = E_g^{\mathrm{exp}} + \tfrac13\Delta_{\mathrm{so}} = 1.74 + 0.14 = 1.88\ \mathrm{eV}.

Since g-xTB overestimates the gap, the shift is negative.

**The g-xTB bulk gap.** Periodic g-xTB did not give a reliable bulk gap (section 8), so the bulk gap is
estimated from the dot series. The PBE bulk gap is known (0.62 eV, ``MATERIAL_DB``), and on every dot
where both methods were run one measures how much higher the g-xTB gap is:

.. math::

   \delta_i = E_g^{\mathrm{gxtb}}(\text{dot } i;\ \text{g-xTB geometry}) - E_g^{\mathrm{PBE}}(\text{dot } i;\ \text{PBE geometry}),
   \qquad E_g^{\mathrm{gxtb,bulk}} \approx E_g^{\mathrm{PBE,bulk}} + \langle\delta\rangle .

Each method is taken on its own geometry, because the NAMD frames have g-xTB bonds (Cd–Se 2.60 Å)
and the PBE route uses PBE structures (2.67–2.69 Å). The g-xTB gaps were computed on the PBE
structures with all coordinates scaled about the centre to a mean Cd–Se bond of 2.60 Å; on the 2 nm
dot this reproduces the gap of the fully relaxed g-xTB structure (4.869 against 4.858 eV).

.. list-table::
   :header-rows: 1

   * - Dot
     - E_g PBE
     - E_g g-xTB, PBE geometry
     - E_g g-xTB, 2.60 Å bonds
     - δ
   * - 1.2 nm
     - 2.639
     - 5.956
     - 6.173
     - 3.53
   * - 1.8 nm
     - 1.318
     - 5.221
     - 5.841
     - 4.52
   * - 2.0 nm
     - 1.457
     - 4.472
     - 4.869
     - 3.41
   * - 2.4 nm
     - 0.680
     - 3.255
     - 3.721
     - 3.04
   * - 2.8 nm
     - 1.454
     - 4.873
     - 5.501
     - 4.05
   * - 3.4 nm
     - 1.151
     - 4.387
     - 4.980
     - 3.83

⟨δ⟩ = 3.73 eV, with a spread of 0.52 eV between dots and no trend with size. Hence

.. math::

   E_g^{\mathrm{gxtb,bulk}} \approx 0.62 + 3.73 = 4.35\ \mathrm{eV},\qquad
   \Delta_{\mathrm{bulk}}^{\mathrm{gxtb}} = 1.88 - 4.35 = -2.47\ \mathrm{eV}.

Written as :math:`\Delta^{\mathrm{gxtb}}_{\mathrm{bulk}} = (E_g^{\mathrm{ref,SF}} - E_g^{\mathrm{PBE,bulk}}) - \langle\delta\rangle`,
the shift is chosen so that, on average over the dots, the g-xTB route gives the same QP gap as the
PBE route with the same experimental reference. The data and the formula are in
``qdex.hardness.GXTB_BULK`` and ``gxtb_bulk_shift``. A periodic g-xTB gap, once available, is entered
as ``gap_gxtb_bulk`` and replaces the estimate.

**Results** (2 nm dot, 50 frames at 300 K): QP gap 2.18 ± 0.09 eV, lowest transition
1.90 ± 0.05 eV. This is lower than the PBE route (sBSE S₁ = 2.78 eV) for three reasons: this dot's δ
(3.41 eV) is 0.32 eV below the average, the gap at 300 K is 0.21 eV lower than at the minimum, and the
PBE route uses the QSGW instead of the experimental reference. The first reason shows the main
limitation: any single dot can deviate by about ±0.5 eV from the estimate.


5.3 Choosing between them
~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * -
     - ``diagonal_stda``, ``functional: gxtb``
     - ``diagonal_sbse``, ``reference: gxtb``
   * - orbital energies
     - g-xTB as they are
     - g-xTB plus one constant shift
   * - electron–hole interaction
     - Grimme's γ\ :sup:`J`, γ\ :sup:`K` (unscreened tail)
     - bulk Resta W and bare MNOK exchange
   * - material parameters
     - none
     - experimental gap, Δso, PBE bulk gap, ⟨δ⟩
   * - lowest transition, 2 nm, 300 K
     - 2.99 eV
     - 1.90 eV
   * - main uncertainty
     - size dependence of the cancellation
     - ±0.5 eV dot-to-dot spread of δ

For the 2 nm CdSe test, ``diagonal_stda`` is the simpler and better-performing choice. The sBSE
variant is useful as a cross-check and for comparison with PBE-based runs. Raw g-xTB transitions
(``independent_dft``) remain available for diagnostics.


6. Single-point QDEX runs on g-xTB orbitals
-------------------------------------------

Any converted frame can also be used in an ordinary QDEX run, for example to compare the full and the
diagonal sTDA or to look at a relaxed structure. The basis is read with ``basis_name: per-atom``:

.. code-block:: yaml

   system:
     mo_file: "MOs.mbse"
     xyz: "frame.xyz"
     basis_txt: "BASIS_GXTB"
     basis_name: "per-atom"
     material: "CDSE"

   quasiparticles:
     model: none

   excitations:
     mode: stda
     functional: gxtb
     nhomos: 25
     nlumos: 25
     full_diag: true

Results at the PBE geometries (vacuum, 25 × 25; S₁ / first bright state in eV, spin-free unless noted):

.. list-table::
   :header-rows: 1

   * - Setting
     - 1.2 nm
     - 2.0 nm
   * - g-xTB gap
     - 5.956
     - 4.472
   * - ``stda``, ``functional: gxtb`` (a_x 0.51, α 4.51, β 8.0)
     - 2.795 / 3.077
     - 2.541 / 2.541
   * - same, with SOC (first bright)
     - 2.879
     - 2.510
   * - ``stda``, ``ax: 1.0`` (global-hybrid formulas)
     - 2.297 / 2.693
     - 2.512 / 2.512
   * - ``diagonal_stda``, ``functional: gxtb``
     - 3.258
     - 2.858
   * - for comparison: PBE route, sBSE (``bulk`` + Resta)
     - 3.323 / 3.617
     - 2.779
   * - experiment (2.0 nm)
     -
     - 2.77–3.31

At 2 nm the two parameter sets agree, because the exciton extends over the whole dot and only the
long-range (bare) part of the attraction matters. At 1.2 nm the exciton is compact and the short-range
damping of the ωB97M-V set weakens the binding by 0.5 eV. At 1.2 nm there is no measured reference;
the binding energy of sTDA on g-xTB (3.2–3.7 eV) is in the range of the ΔW models of QDEX
(3.1–3.8 eV), and the lower S₁ comes from the smaller (evGW-like) gap.


7. What transfers from PBE and what does not
--------------------------------------------

**Interaction terms.** On the same structure (2 nm, 10 × 10 transitions, Resta W, Mulliken charges):

.. list-table::
   :header-rows: 1

   * - Term (HOMO → LUMO)
     - g-xTB
     - PBE (DZVP)
   * - K\ :sup:`d`, direct (mean over 100 transitions)
     - 0.26 eV (0.28 eV)
     - 0.30 eV (0.30 eV)
   * - 2K\ :sup:`x`, exchange
     - 0.006 eV
     - 0.084 eV
   * - oscillator strength
     - 0.046
     - 0.284

The direct term transfers well: it depends on where each orbital's density sits, which a minimal basis
describes well. Exchange and oscillator strengths depend on the overlap density of hole and electron
and are several times too small in the compact minimal basis. For carrier cooling this matters little;
oscillator strengths serve to separate bright from dark states. Absorption spectra would need a larger
basis on selected frames.

**Level structure.** What drives carrier cooling is the spacing of levels inside the bands, which a
constant shift cannot repair:

.. list-table::
   :header-rows: 1

   * - 2 nm dot
     - PBE
     - g-xTB
     - GFN2-xTB
     - sTDA-xTB
   * - LUMO → LUMO+1 spacing
     - 0.28 eV
     - 0.41 eV
     - 0.23 eV
     - 0.01 eV
   * - conduction states within 1 eV of the LUMO
     - 7
     - 6
     - 17
     - 50

g-xTB keeps the discrete conduction band of PBE, slightly stretched, as a hybrid functional would.

**Spin–orbit coupling.** The GTH spin–orbit operator of QDEX (:doc:`/relativity/index`) is built from
atomic projectors integrated over whatever basis is given, so it also runs on g-xTB orbitals. It gives
about two thirds of the PBE strength: the gap of the 2 nm dot shrinks by 0.052 eV with g-xTB orbitals
against 0.078 eV with PBE. CdSe has weak spin–orbit coupling at the band edges, so this is used as is.

**Geometry.** g-xTB contracts the dot: mean Cd–Se 2.60 Å, against 2.67 Å for PBE and 2.63 Å for bulk
zinc-blende CdSe. At fixed method this raises the gap by about 0.4 eV relative to the PBE structure.


8. Alternatives that were tested
--------------------------------

**GFN2-xTB.** Gap 2.21 eV on the 2 nm dot (2.68 eV at the g-xTB geometry), closer to PBE than g-xTB.
But its conduction band is too dense (17 states within 1 eV, against 7 for PBE), which would speed up
electron cooling artificially. It is also slower here (87 SCF cycles; 9.8 s per energy and gradient
against 2.1 s).

**IPEA-xTB (``--vipea``).** Vertical IP − EA from ΔSCF, 2.63 eV for the 2 nm dot in vacuum. It uses
the GFN1-type IPEA1 Hamiltonian, not the Hamiltonian of the orbitals, and the charging correction it
adds hardly changes with geometry at fixed size. A constant shift is cheaper and equivalent.

**sTDA-xTB (xtb4stda).** Grimme's tight-binding method made for sTDA (J. Chem. Phys. 145, 054103
(2016); Mol. Phys. 117, 1104 (2019)). It is its own Hamiltonian, not GFN2: charges from a minimal-basis
step, then one diagonalization in an extended basis with diffuse sp shells, and a +3.1 eV shift of all
virtual orbitals. In the CdSe dots the diffuse Cl and Se shells turn the bottom of the conduction band
into a quasi-continuum (50 states within 1 eV, 22 % Cl character in the LUMO), which would remove the
electron-cooling bottleneck. It has no periodic mode. Its molden files can be converted, but it is not
recommended for dots.

**Periodic g-xTB for the bulk gap.** With the macOS arm64 bleed build:

* a 64-atom zinc-blende cell (Γ point) collapses in the SCF;
* a 216-atom cell reaches a density error of 4·10⁻⁶ with Broyden damping 0.10 (``$scc broydamp``), but
  the energy keeps fluctuating at 10⁻⁵ Eh, above the 10⁻⁶ Eh criterion, so xtb prints no orbital
  energies;
* any non-default ``--acc`` crashes the periodic run at start-up, and some runs crash during the SCF;
* damping 0.05 with 1000 K smearing diverges.

The Linux build on a cluster may behave better. Production periodic xTB calculations normally use
k-point sampling in DFTB+ or CP2K (GFN1, GFN2), which do not offer g-xTB.


9. Limitations and open points
------------------------------

* The g-xTB bulk shift (section 5.2) is provisional: ±0.5 eV from dot to dot.
* The size dependence of the sTDA route (section 5.1) is untested beyond 2 nm.
* The diagonal approximation raises S₁ by 0.3–0.5 eV relative to the full sTDA in these dots.
* Exchange and oscillator strengths are underestimated in the minimal basis.
* g-xTB has no Python API yet; QDEX drives the command-line program and parses molden files. When
  g-xTB is released in ``tblite``, the orbitals could be taken directly from its Python interface.
* The xtb ``--md`` driver does not work with ``--gxtb``; use ``qdex.xtb.md``.


10. Files of the test case
--------------------------

The 2 nm test lives in ``CdSe/2.0nm/gxtb_namd`` of the CdSe project data:

* ``env.sh``: environment;
* ``start.xyz``: g-xTB-relaxed structure;
* ``frames/``: 50 frames; ``equil_log.csv``, ``prod_log.csv``: MD logs;
* ``config_namd_gxtb_stda.yaml``: precompute with ``diagonal_stda`` (section 5.1);
* ``config_namd_gxtb.yaml``: precompute with ``diagonal_sbse`` and the g-xTB bulk shift (section 5.2);
* ``bulk/``: dot series and periodic tests.

The single-point sTDA tests (section 6) are in ``CdSe/gxtb_stda_test/{1.2nm,2.0nm}``.
