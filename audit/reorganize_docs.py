"""Lossless topic migration from the audited RST snapshot.

Creates topical pages from disjoint manifest units, retaining byte-identical copies of
all eight source parts in docs/_legacy_parts. Historical paths become redirect stubs.
"""
from __future__ import annotations
from collections import defaultdict
from pathlib import Path
import hashlib
import json
import re
import shutil

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / 'docs'
manifest = json.loads((ROOT/'audit/documentation_inventory.json').read_text())
api = json.loads((ROOT/'audit/api_signatures.json').read_text())
by_key = {(x['file'].removesuffix('.py').replace('/','.'),x['name']):x['signature'] for x in api}
parts = [f for f in manifest['files'] if '/part' in f['path']]
legacy = DOCS/'_legacy_parts'
legacy.mkdir(exist_ok=True)
for f in parts:
    original=ROOT/f['path']
    archived=legacy/(original.parent.name+'.rst')
    if archived.exists():
        assert hashlib.sha256(archived.read_bytes()).hexdigest()==f['sha256'], archived
    else:
        assert hashlib.sha256(original.read_bytes()).hexdigest()==f['sha256'], original
        shutil.copyfile(original,archived)
old_index=legacy/'old_index.rst'
if not old_index.exists():
    shutil.copyfile(DOCS/'index.rst',old_index)

# Exact functions, CLI flags and YAML keys for each topic family; page overrides
# narrow the pointer to the relevant implementation without changing the source.
group_contract={
 'electronic_structure':('qdex.pdos_coop','compute_pdos_and_coop','--charge_type, --run_fuzzy','physics.charge_type, fuzzy.run'),
 'relativity':('qdex.soc_utils','compute_spinor_subspace','--soc_flag, --gth_file','physics.soc, physics.soc_window_ev'),
 'interactions':('qdex.hardness','build_gamma','--2e-integrals, --kernel, --eps-out','physics.2e-integrals, physics.kernel, physics.eps_out'),
 'quasiparticles':('qdex.hardness','estimate_gw_qp_gap','--qp_gap, --dynamic_z','physics.qp_gap, physics.dynamic_z'),
 'excitons':('qdex.solver','solve','--excitation-mode, --include-direct-eh, --include-exchange','physics.excitation_mode, physics.include_direct_eh, physics.include_exchange'),
 'exciton_analysis':('qdex.nto','analyze_nto_state','--nto, --nto-states','analysis.nto, analysis.nto_states'),
 'dynamics':('qdex.namd.precompute','precompute_namd_data','--namd','namd.engine, namd.dt_fs'),
 'spectroscopy':('qdex.namd.transient_absorption','compute_transient_absorption','--sigma','namd.transient_absorption'),
 'recombination':('qdex.auger','calculate_auger_rates','--auger, --auger-channel','auger.run, auger.channel'),
 'workflows':('qdex.cli','main','--config','system.*, physics.*, namd.*'),
 'reference':('qdex.hardness','get_cluster_size_metrics','--material','system.material, physics.material'),
 'validation':('qdex.solver','solve','--nroots, --tol','bse.nroots, bse.tol'),
}
overrides={
 'interactions/representations':('qdex.integrals','compute_two_electron_ao'),
 'interactions/screening':('qdex.hardness','build_xs_kernel'),
 'interactions/qp_screening':('qdex.hardness','build_dim_screening_factors'),
 'quasiparticles/anchor':('qdex.hardness','estimate_gw_qp_gap'),
 'quasiparticles/delta_w':('qdex.hardness','estimate_sgw_resta_qp_gap'),
 'quasiparticles/delta_w_implementation':('qdex.hardness','estimate_sgw_dim_qp_gap'),
 'quasiparticles/dynamic_z':('qdex.hardness','compute_dynamic_z'),
 'quasiparticles/gap_iteration':('qdex.hardness','estimate_evgw_dim_qp_gap'),
 'quasiparticles/orbital_iteration':('qdex.hardness','estimate_qsgw_dim_qp_gap'),
 'excitons/frameworks':('qdex.exciton_hamiltonian','independent_transition_energies'),
 'excitons/solvers':('qdex.davidson','davidson'),
 'exciton_analysis/descriptors':('qdex.exciton_analysis','ExcitonAnalyzer'),
 'exciton_analysis/ntos':('qdex.nto','compute_nto_pairs'),
 'dynamics/decoherence':('qdex.namd.integrator','apply_edc_decoherence'),
 'dynamics/nacs_tracking':('qdex.namd.precompute','align_phases_and_crossings'),
 'dynamics/pme':('qdex.namd.master_equation','propagate_pme_tensor'),
 'dynamics/fssh_edc':('qdex.namd.surface_hopping','run_namd_dynamics'),
 'dynamics/dish':('qdex.namd.integrator','step_dish_batch'),
 'spectroscopy/spectral_density':('qdex.namd.analysis','compute_spectral_density'),
 'recombination/radiative':('qdex.hardness','compute_radiative_rates'),
 'recombination/energy_gap_law':('qdex.hardness','compute_energy_gap_law_rate'),
 'recombination/fcwd':('qdex.hardness','compute_fcwd_rate'),
 'recombination/auger_matrix_elements':('qdex.auger','compute_auger_matrix_element'),
 'recombination/auger_ecsh':('qdex.namd.surface_hopping','run_namd_dynamics'),
 'reference/materials':('qdex.hardness','estimate_gw_qp_gap'),
}
titles={
 'electronic_structure':'Electronic structure and orbital analysis',
 'relativity':'Relativity and spin',
 'interactions':'Two-body interactions',
 'quasiparticles':'Quasiparticle models',
 'excitons':'Excitons and BSE',
 'exciton_analysis':'Exciton wavefunctions',
 'dynamics':'Nonadiabatic dynamics',
 'spectroscopy':'Spectroscopy',
 'recombination':'Recombination and Auger processes',
 'workflows':'Worked workflows',
 'reference':'Reference data',
 'validation':'Validation and limiting cases',
}
method_notes={
 'interactions/representations': 'The ``xs`` backend evaluates analytical Gaussian integrals of the restricted form ``(mu mu | nu nu)``. It does not store the full four-index electron-repulsion tensor; the MO interactions are density-pair approximations. The MNOK heteronuclear damping currently uses the mean of inverse hardnesses, not ``2/(eta_A+eta_B)``.',
 'interactions/screening': 'Resta and DIM screening builders omit external-medium screening from their BSE direct kernel. This is a model choice; their QP paths include an external reaction term. The two pieces have not been shown to cancel for arbitrary electron and hole densities.',
 'interactions/architecture': 'The ``--2e-integrals`` axis chooses the representation of the bare interaction; ``--kernel`` chooses a screening builder. Some kernel names force AO or atom resolution in ``qdex.solver.ExcitonSolver``. Inspect that dispatcher for supported combinations.',
 'interactions/environment': 'The solvent reaction term is a softened dielectric-boundary model. Its sign reverses when the exterior electronic permittivity exceeds the bulk value. It is not the complete Green function of a dielectric sphere.',
 'quasiparticles/anchor': '``sgw-anchor`` interpolates a bulk gap correction and one finite vacuum anchor. Its exponent and regularization length are model parameters, not uniquely fixed by two endpoints.',
 'quasiparticles/delta_w': 'These routines use static atom-centered screened-interaction contrasts and charging energies. They do not evaluate a dynamical GW self-energy; ``v_xc`` cancellation is an approximation and is least reliable near reconstructed or trapped surfaces.',
 'quasiparticles/delta_w_implementation': 'DIM maps a three-field polarizability response to pairwise screening with additional normalization. The mapping is heuristic and must be benchmarked against a microscopic dielectric response.',
 'quasiparticles/dynamic_z': '``compute_dynamic_z`` is an empirical scalar damping formula. It does not evaluate the derivative of a computed frequency-dependent self-energy or enforce an f-sum rule. The default pole energy is 15 eV; material-specific plasmon energies are not read from a database.',
 'quasiparticles/gap_iteration': 'The ``evgw-*`` implementation iterates an effective gap and screening model; it does not update a complete Green function or the state-resolved GW self-energy.',
 'quasiparticles/orbital_iteration': 'The ``qsgw-*`` implementation relaxes orbitals with a static atom-block COHSEX-like operator. It is distinct from conventional QSGW. The CLI now passes its returned QP eigenvalues to the BSE solver with zero additional scissor, while retaining DFT eigenvalues separately for the DFT framework and filtering.',
 'quasiparticles/edge_partition': 'Frontier fractions are model allocations of a gap correction. Signed solvent contributions are evaluated separately; the fractions are bounded and default to 50:50 when the relevant contrast vanishes. Absolute IP/EA still require a vacuum reference or a declared monomer-anchor reconstruction.',
 'quasiparticles/alignment': 'The current CLI may reconstruct absolute edges from monomer entries in ``MATERIAL_DB`` even for an anchor-free gap model. Distinguish a computed gap from calibrated vacuum-level alignment.',
 'quasiparticles/model_selection': 'The labels ``sgw``, ``evgw`` and ``qsgw`` in QDEX name reduced models. They should not be interpreted as a claim of numerical equivalence to conventional GW/QSGW calculations. Validate each against a common reference set.',
 'excitons/frameworks': '``diagonal_bse`` omits all off-diagonal transition mixing. It can miss a substantial part of binding; it cannot generally be assumed to reproduce the bulk Wannier exciton. Exchange and direct attraction have independent switches.',
 'excitons/solvers': '``bse`` solves the resonant/Tamm-Dancoff matrix. Davidson returns selected roots; ``diagonal_bse`` uses the Hamiltonian diagonal without an iterative eigensolve.',
 'validation/bulk_exciton_limit': 'The Wannier–Mott limit requires converged coupled BSE transitions, effective masses and a macroscopic screened tail. It is not guaranteed by a pair potential alone, and diagonal BSE does not generally recover finite bulk binding.',
 'reference/materials': 'Material entries are model inputs, with scalar-relativistic PBE and GW gap conventions where given. Current ``MATERIAL_DB["CDSE"]`` uses 0.64 eV PBE, 1.91 eV GW, and epsilon_inf=6.2; alternative CdSe numbers must be separate, cited data sets.',
}

# Correct stale names in the topical copies. The archived bytes remain untouched.
fix_names={
 'estimate_anchor_scaled_qp_gap':'estimate_gw_qp_gap',
 'build_damped_mnok_matrix':'build_gamma',
 'qdex.solver.solve_diagonal_bse':'qdex.solver.ExcitonSolver.solve',
 'qdex.solver.davidson':'qdex.davidson.davidson',
 '--ci_threshold':'--e_thresh',
}
# The archived text is exact. These corrections apply only to reader-facing
# topical copies where the historical prose makes a stronger claim than the code.
fix_claims={
 'from the **Plasmon-Pole Model (PPM) constrained by the :math:`f`-sum rule**': 'from an **empirical plasmon-energy damping ansatz** (no frequency-dependent self-energy or :math:`f`-sum-rule constraint is evaluated)',
 'Retrieves material-specific valence plasmon energy :math:`\\Omega_p \\approx 15 - 20\\text{ eV}` from ``MATERIAL_PLASMON_ENERGIES``.': 'Uses the optional ``omega_p_ev`` argument (15 eV by default); no material plasmon-energy table is queried.',
 'Compute state-dependent Z_p via PPM f-sum rule': 'Apply empirical state-dependent Z_p damping',
 'State-dependent Z_p via PPM f-sum rule': 'Empirical state-dependent Z_p damping',
 'State-dependent Z_p via Plasmon-Pole f-sum rule': 'Empirical state-dependent Z_p damping',
 'Dynamically compute Z_p via PPM f-sum rule': 'Apply empirical state-dependent Z_p damping',
 'Local :math:`v_{xc}` cancels identically; parameter-free.': 'Local :math:`v_{xc}` differences are neglected as a model approximation; bulk reference parameters remain required.',
 'dynamic plasmon-pole renormalization :math:`Z_p`': 'empirical damping factor :math:`Z_p`',
 'full AO-basis orbital relaxation (:math:`qsGW`)': 'static AO-basis orbital relaxation (``qsgw-*``)',
 'full AO-basis quasiparticle self-consistent GW (:math:`qsGW`)': 'the static AO-basis orbital-relaxation model (``qsgw-*``)',
 'Perform full AO-basis Quasiparticle Self-Consistent GW (:math:`qsGW`) orbital update.': 'Run the static AO-basis COHSEX-like orbital-relaxation model.',
 'Z_p = \\left( 1 - \\left. \\frac{\\partial \\operatorname{Re}\\Sigma_p}{\\partial \\omega} \\right|_{\\varepsilon_p} \\right)^{-1} = \\left( 1 + \\frac{\\Delta \\Sigma_p^{\\mathrm{stat}}}{\\tilde{\\Omega}_p} \\right)^{-1}': 'Z_p^{\\mathrm{model}} = \\left( 1 + \\frac{\\max(0,\\Delta \\Sigma_p^{\\mathrm{stat}})}{\\tilde{\\Omega}_p} \\right)^{-1}',
 'Z \\in [0.50, 0.98]': 'Z \\in [0.50, 0.99]',
 '**exact four-center two-electron Gaussian repulsion integrals**': '**analytical AO density-pair Coulomb integrals**',
 'exact AO-to-MO contraction': 'density-pair AO-to-MO contraction',
 'Exact Analytical Gaussian Representation': 'Analytical AO Density-Pair Representation',
 'Exact four-center integrals': 'Analytical AO density-pair integrals',
 'full AO-basis **Quasiparticle Self-Consistent GW (:math:`qsGW`)**': 'full AO-basis **static COHSEX-like orbital relaxation (``qsgw-*``)**',
 'Full AO-basis Quasiparticle Self-Consistent GW': 'Full AO-basis static orbital relaxation model',
 'Eigenvalue self-consistent GW': 'Effective-gap self-consistent screening model',
 'Eigenvalue Self-Consistent GW': 'Effective-gap Self-Consistent Screening Model',
 'exact same ratio partitions the bulk reference Hamiltonian': 'same model allocation is used to partition the bulk reference Hamiltonian',
 '**Rigorous qsGW Background**': '**Model bulk projector allocation**',
 'Scales strictly as :math:`O(N_{\\mathrm{atoms}}^2)`. Memory requirements are minimal (:math:`< 10\\text{ MB}`).': 'The atom-pair kernel stores :math:`O(N_{\\mathrm{atoms}}^2)` elements; actual memory depends on system size and solver intermediates.',
 'accounting for **over 90% of the total exciton binding energy** :math:`E_b`': 'whose share of the coupled exciton binding energy is system-dependent; the supplied CdSe example gives 58.3% for diagonal BSE',
 'Diagonal BSE provides an accurate, energy-conserving potential energy surface at a fraction of the cost.': 'Diagonal BSE provides a less expensive approximate surface; its error must be checked against coupled BSE for the chosen active space.',
 'gives the true exciton binding energy :math:`E_b \\approx 280\\text{ meV}`': 'gives a model-dependent estimate of binding relative to the chosen QP reference, not a universal 280 meV value',
 'exact analytical 4-center Gaussian integrals': 'analytical Gaussian AO density-pair integrals',
 'exact analytical Gaussian integrals': 'analytical Gaussian AO density-pair integrals',
 'Parameter-free exact Gaussian integrals and microscopic RPA screening.': 'Analytical AO density-pair Coulomb integrals and model RPA screening; benchmark accuracy requires validation.',
 'two rigorously computed physical limits': 'a finite-cluster calibration and a bulk reference',
 'these exact same microscopic fractions': 'the same model-dependent microscopic fractions',
 'In low-permittivity solvents, dielectric confinement strongly enhances the direct electron-hole attraction :math:`K^d`.': 'In a consistent charged/neutral treatment this boundary can change optical energies; QDEX currently does not include ``eps_out`` in its Resta/DIM BSE direct kernel.',
 'the BSE optical transition energy :math:`\\Omega_1` smoothly converges to the bulk band edge minus the Wannier-Mott binding energy': 'a converged bulk BSE optical transition should approach the bulk band edge minus the Wannier-Mott binding energy',
}

role_pattern=re.compile(r':[a-z][a-z0-9_-]*:`[^`]+`')
def render_roles_outside_strong(line):
    """RST roles are not parsed when nested inside ``**strong**`` markup."""
    if line.count('**')%2:
        return line  # E.g. a literal exponent ``x**2`` is not strong markup.
    line=re.sub(r'\*\*([^*`]+?) \(:math:`([^`]+)`\)\*\*',
                r'**\1** (:math:`\2`)',line)
    line=re.sub(r'\*\*([^*`]+?) :math:`([^`]+)`\*\*',
                r'**\1** :math:`\2`',line)
    pieces=line.split('**')
    out=pieces[0]
    for i in range(1,len(pieces)):
        out+=(pieces[i] if i%2==0 or role_pattern.search(pieces[i]) else '**'+pieces[i]+'**')
    stripped=out.strip()
    if (':math:' in stripped and stripped.startswith('*')
            and not stripped.startswith('**') and stripped.endswith('*')
            and not stripped.endswith('**')):
        prefix=out[:len(out)-len(out.lstrip())]
        suffix='\n' if out.endswith('\n') else ''
        return prefix+stripped[1:-1]+suffix
    return out
units=defaultdict(list)
for u in manifest['units']:
    if not any(u['source']==f['path'] for f in parts): continue
    destination=u['proposed_destination']
    if destination=='docs/index.rst': continue
    if destination=='docs/workflows/index.rst': destination='docs/workflows/overview.rst'
    units[destination].append(u)

for destination, records in sorted(units.items()):
    relative=Path(destination).relative_to('docs')
    key=relative.with_suffix('').as_posix()
    group=relative.parts[0]
    page=DOCS/relative
    page.parent.mkdir(parents=True,exist_ok=True)
    title=key.split('/')[-1].replace('_',' ').capitalize()
    chunks=[title,'='*len(title),'',f'Part of :doc:`/{group}/index`.', '']
    if key in method_notes:
        chunks+=['.. important::','', '   '+method_notes[key], '']
    module,name=overrides.get(key,group_contract[group][:2])
    default_module,default_name,flags,yaml=group_contract[group]
    signature=by_key.get((module,name))
    chunks+=['.. rubric:: Theory and QDEX implementation','','The detailed theory and worked equations follow below. The corresponding entry point is:', '', f'* Module: ``{module}``',f'* Callable: ``{module}.{name}``',f'* CLI: ``{flags}``',f'* YAML: ``{yaml}``','']
    if signature:
        chunks+=['.. code-block:: python','', '   '+signature,'']
    chunks+=['.. rubric:: Detailed derivations and reference material','']
    for u in records:
        source=legacy/(Path(u['source']).parent.name+'.rst')
        lines=source.read_text().splitlines(keepends=True)
        body=''.join(lines[u['start_line']-1:u['end_line']])
        assert hashlib.sha256(body.encode()).hexdigest()==u['sha256']
        if u['start_line']==1:
            bl=body.splitlines(keepends=True)
            if len(bl)>=2 and re.fullmatch(r'[=~^\-]{4,}\s*',bl[1]): body=''.join(bl[2:])
        for old,new in fix_names.items(): body=body.replace(old,new)
        for old,new in fix_claims.items(): body=body.replace(old,new)
        corrected=[render_roles_outside_strong(line) for line in body.splitlines(keepends=True)]
        for i in range(1,len(corrected)):
            underline=corrected[i].rstrip('\n')
            heading=corrected[i-1].rstrip('\n')
            if (heading and not heading[0].isspace()
                    and re.fullmatch(r'[=~^\-]{4,}',underline)
                    and len(underline)<len(heading)):
                corrected[i]=underline[0]*len(heading)+'\n'
        body=''.join(corrected)
        if not lines[u['start_line']-1][:1].isspace():
            chunks+=[f'.. rubric:: From ``{u["source"]}:{u["start_line"]}-{u["end_line"]}``','']
        chunks+=body.splitlines()
        chunks+=['']
    page.write_text('\n'.join(chunks).rstrip()+'\n')

# Create per-family indexes and topical cross-links.
by_group=defaultdict(list)
for destination in units:
    relative=Path(destination).relative_to('docs')
    by_group[relative.parts[0]].append(relative.stem)
by_group['validation'] += ['cdse_benchmark', 'environment_cancellation']
by_group['reference'] += ['documentation_math']
for group,pages in sorted(by_group.items()):
    title=titles[group]
    lines=[title,'='*len(title),'',f'Topics in {title.lower()} pair each derivation with the QDEX callable, CLI flags and YAML keys used by the implementation.','', '.. toctree::','   :maxdepth: 1','']
    lines += ['   '+p for p in sorted(pages)]
    (DOCS/group/'index.rst').write_text('\n'.join(lines)+'\n')

# Preserve old URLs as stubs; exact bytes remain downloadable from _legacy_parts.
for f in parts:
    original=ROOT/f['path']
    group=original.parent.name
    title=(legacy/(group+'.rst')).read_text().splitlines()[0]
    links=[u['proposed_destination'] for u in manifest['units'] if u['source']==f['path']]
    groups=sorted({Path(p).parts[1] for p in links if p!='docs/index.rst'})
    lines=[':orphan:','',title,'='*len(title),'',
        'This historical chapter has been reorganized by scientific topic. Its full original text is retained as a download, and all sections appear in the topic pages.',
        '', f'Original RST: :download:`{group} source <../_legacy_parts/{group}.rst>`.', '', 'Start with:']
    lines += [f'* :doc:`/{g}/index`' for g in groups]
    original.write_text('\n'.join(lines)+'\n')

main=['QDEX documentation','==================','','QDEX connects post-DFT electronic structure, quasiparticle models, BSE/TDA excitations and carrier dynamics for molecules and semiconductor nanoclusters.','',
'The topic pages place the physical equations beside the QDEX function, flags and YAML configuration that implement them. Model-status notes distinguish established results from QDEX approximations.','',
'Quick entry points: :doc:`getting_started/installation`, :doc:`getting_started/quickstart`, :doc:`getting_started/configuration`, and :doc:`api/index`.','']
for caption,refs in [
 ('Start here',['getting_started/installation','getting_started/quickstart']),
 ('Electronic structure and relativity',['electronic_structure/index','relativity/index']),
 ('Interactions, quasiparticles and excitons',['interactions/index','quasiparticles/index','excitons/index','validation/index']),
 ('Analysis, dynamics and spectroscopy',['exciton_analysis/index','dynamics/index','spectroscopy/index','recombination/index']),
 ('Workflows and reference',['workflows/index','reference/index','getting_started/configuration','api/index'])]:
    main+=['.. toctree::','   :maxdepth: 2',f'   :caption: {caption}','']+['   '+x for x in refs]+['']
main+=['Indices and tables','==================','','* :ref:`genindex`','* :ref:`modindex`','* :ref:`search`','']
(DOCS/'index.rst').write_text('\n'.join(main))
print(f'Created {len(units)} topical pages in {len(by_group)} sections; archived {len(parts)} historical chapters.')
