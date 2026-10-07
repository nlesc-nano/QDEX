Nacs tracking
=============

Part of :doc:`/dynamics/index`.

.. rubric:: QDEX implementation

* Modules: ``qdex.namd.precompute`` (overlaps, tracking), ``qdex.namd.nac`` (couplings)
* Callables: ``align_phases_and_crossings``, ``align_spinor_phases_and_crossings``,
  ``align_degenerate_blocks``, ``qdex.namd.nac.nac_logm``, ``qdex.namd.nac.compute_nac_files``
* CLI: ``--namd-precompute``, then ``--namd-nac`` and ``--namd-decoherence``
* YAML: ``namd.tracking`` (``phase_correction``, ``hungarian_tracking``, ``degeneracy_tol_ev``),
  ``namd.nac`` (``scheme``, ``workers``, ``threads_per_worker``), ``namd.dynamics.nac_scheme``

5. Trajectory Precomputation & Wavefunction Tracking
-----------------------------------------------------

For every pair of consecutive MD frames :math:`t_k, t_{k+1} = t_k + \Delta t` (the MD time step,
``namd.trajectory.dt_nuc_fs``) the precompute builds the overlaps of the active orbitals (or spinors),
fixes their labels and phases, and stores them in ``step_<k>_to_<k+1>.npz``. The couplings are made from
these overlaps, in the precompute stage ``--namd-nac`` (stored) or in the dynamics.

Per step the operations are, in this order:

1. overlap :math:`S = C(t)^\dagger\, S^{\mathrm{AO}}(t, t+\Delta t)\, C(t+\Delta t)` of the occupied and of
   the virtual active states;
2. Hungarian relabelling of trivial crossings (states whose overlap with their own label is below 0.5);
3. phase alignment (sign of each real orbital, U(1) phase of each spinor);
4. SOC only: parallel transport inside each Kramers pair;
5. energies, :math:`K^d` and :math:`K^x` permuted with the states; pair energies rebuilt;
6. couplings :math:`d = \log(U)/\Delta t` from the Löwdin-orthonormalised overlap (``--namd-nac``).

Cross-frame overlaps
~~~~~~~~~~~~~~~~~~~~

.. math::

   S_{ij}(t, t+\Delta t) = \langle \phi_i(t) | \phi_j(t+\Delta t) \rangle
   = \sum_{\mu\nu} C_{\mu i}(t)\, S^{\mathrm{AO}}_{\mu\nu}(t, t+\Delta t)\, C_{\nu j}(t+\Delta t),
   \qquad
   S^{\mathrm{AO}}_{\mu\nu} = \int \chi_\mu(\mathbf{r}; \mathbf{R}(t))\, \chi_\nu(\mathbf{r}; \mathbf{R}(t+\Delta t))\, d\mathbf{r},

with the basis functions of both geometries (Libint2). For spinors :math:`\psi_j = \sum_n (U^\alpha_{nj}\,
\phi_n \alpha + U^\beta_{nj}\, \phi_n \beta)` over the active MOs,

.. math::

   S^{\mathrm{spinor}} = U^{\alpha\dagger}(t)\, M\, U^\alpha(t+\Delta t) + U^{\beta\dagger}(t)\, M\, U^\beta(t+\Delta t),
   \qquad M = C^{\mathrm{act}}(t)^\mathsf{T} S^{\mathrm{AO}} C^{\mathrm{act}}(t+\Delta t).

Electrons and holes are tracked separately (``S_occ``, ``S_virt``): in the diagonal BSE an exciton
:math:`|ia\rangle` moves by one-body steps of its hole :math:`i` or of its electron :math:`a`.

Hungarian matching for trivial crossings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A trivial crossing is a label swap between two states that do not interact: within one step the diagonal
overlap collapses and the character moves to an off-diagonal element; following the energy label would
put a spike :math:`|S_{ij}|/\Delta t` into the coupling. A state is **locked** to its label when
:math:`|S_{ii}| \ge 0.5` (avoided crossings stay in the adiabatic basis that surface hopping propagates);
the unlocked states are matched by the Hungarian algorithm (``scipy.optimize.linear_sum_assignment``) with

.. math::

   \pi = \arg\min_{\pi} \sum_{i\ \mathrm{unlocked}} \left(1 - |S_{i,\pi(i)}|^2\right),

over the unlocked block only (a locked state keeps its column, so no other state can take it). The
permutation is applied to the columns of :math:`S`, to the coefficients and to the energies,
:math:`K^d`, :math:`K^x` and dipoles of frame :math:`t+\Delta t`, which then becomes the reference of the
next step: labels and phases are carried from frame 0 through the whole trajectory. The number of states
whose label changed rank in a step is stored (``n_swap_occ``, ``n_swap_virt``) and printed.

Phase alignment
~~~~~~~~~~~~~~~

Eigenvectors come with an arbitrary phase; after the matching each state of :math:`t+\Delta t` is
multiplied by :math:`e^{-i\theta_j}` with :math:`\theta_j = \arg S_{jj}` (spinors) or by
:math:`\operatorname{sign}\, S_{jj}` (real orbitals), so that :math:`S_{jj} \ge 0`.

Kramers pairs of the SOC spinors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
0.43 fs\ :sup:`−1`, ħ\ *d* ≈ 0.14 eV) in every step; with the parallel transport it is zero, the other
couplings are unchanged and the median :math:`|S_{jj}|` rises from 0.76 to 0.92 (Br) and from 0.66 to
0.87 (Cl). SOC precomputes made before this correction (2026-10-07) carry the spurious couplings.

Non-adiabatic couplings
~~~~~~~~~~~~~~~~~~~~~~~

The coupling :math:`d_{ij} = \langle \phi_i | \partial_t \phi_j \rangle` is taken constant over the step
(the integrators hold it fixed between frames). Two schemes (``dynamics.nac_scheme``):

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

  Propagating the amplitudes with :math:`e^{-d\,\Delta t}`, as the integrators do, then reproduces the
  overlap exactly, :math:`c(t+\Delta t) = U^\dagger c(t)`. :math:`\log U` is evaluated from the
  eigenvectors of the Hermitian :math:`(U - U^\dagger)/2i` (which :math:`U` shares), with a residual check
  and ``scipy.linalg.logm`` as fallback for rotations beyond 90°; it is real for real orbitals and
  anti-Hermitian for spinors.

In the valence band of the CsPbX\ :sub:`3` dots (1000 occupied orbitals over 2 eV, ~5 meV apart) states
rotate by up to ~50° in a 2 fs step: the logm coupling of the hole states is 16 % (CsPbBr\ :sub:`3`) and
37 % (CsPbCl\ :sub:`3`) larger than HST (median per state; 49 % deep in the band of CsPbCl\ :sub:`3`), i.e.
HST rates (∝ :math:`|d|^2`) are 35–120 % too low. At the band edges and for the electrons the two agree within 2–6 %.

The logm couplings cost an SVD and an eigendecomposition per step (≈ 1–4 s for 1300 real, 10–35 s for
2600 complex states), so ``qdex --config ... --namd-nac`` computes them once after the precompute, in
parallel over the steps (``namd.nac.workers`` processes × ``threads_per_worker``; restartable), and stores
``nac_<k>_to_<k+1>.npz`` next to the step files. Without them the dynamics computes the chosen scheme on
the fly.

**Hole channel.** The hole states are :math:`a_i|\Phi_0\rangle`; their coupling is
:math:`\langle \Phi_0 a_i^\dagger | \partial_t\, a_j \Phi_0 \rangle = d_{ij}^*`. The dynamics uses
:math:`d^*` for the hole channel of spinors (for real orbitals :math:`d^* = d`).

**Window edges.** States at the edges of the active window lose norm to states outside it
(:math:`1 - \sum_j |S_{ij}|^2` up to 0.1–0.8 for the outermost few states of a 1300 × 800 window, < 10\ :sup:`−3`
near the gap); the warning in the log refers to them. Choose the window so that the carriers never come
near its edges (about 1 eV beyond the largest excess energy of the pump).
