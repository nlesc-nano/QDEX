"""Inventory all RST content in disjoint, hash-addressed migration units."""
from pathlib import Path
import ast, hashlib, json, re, subprocess
root=Path(__file__).resolve().parents[1]
out=root/'audit'
routes={
 'part1_ground_state':[(1,'electronic_structure/orbitals'),(51,'electronic_structure/populations_pdos'),(107,'electronic_structure/localization'),(140,'electronic_structure/bonding_coop'),(171,'electronic_structure/fuzzy_bands'),(219,'electronic_structure/orbital_cubes'),(252,'electronic_structure/configuration')],
 'part2_soc':[(1,'relativity/foundations'),(40,'relativity/material_effects'),(76,'relativity/pseudopotentials'),(111,'relativity/spinor_hamiltonian'),(187,'relativity/assembly'),(205,'relativity/unrestricted'),(243,'relativity/configuration')],
 'part3_gw_scissor':[(1,'quasiparticles/foundations'),(48,'interactions/architecture'),(135,'quasiparticles/dft_reference'),(154,'quasiparticles/cost'),(165,'quasiparticles/model_selection'),(189,'quasiparticles/anchor'),(261,'quasiparticles/delta_w'),(318,'interactions/qp_screening'),(356,'interactions/asymptotics'),(388,'quasiparticles/edge_partition'),(430,'quasiparticles/delta_w_implementation'),(469,'quasiparticles/dynamic_z'),(499,'quasiparticles/gap_iteration'),(520,'quasiparticles/orbital_iteration'),(581,'interactions/environment'),(601,'quasiparticles/alignment'),(700,'reference/materials'),(1001,'quasiparticles/model_selection'),(1063,'workflows/qp_presets'),(1131,'electronic_structure/qp_analysis'),(1147,'quasiparticles/configuration')],
 'part4_excited_states':[(1,'excitons/foundations'),(80,'interactions/architecture'),(143,'interactions/representations'),(225,'interactions/screening'),(291,'validation/bulk_exciton_limit'),(314,'excitons/frameworks'),(369,'spectroscopy/dipoles_oscillators'),(396,'excitons/solvers'),(411,'excitons/model_selection'),(461,'workflows/exciton_presets'),(544,'excitons/configuration')],
 'part5_exciton_analysis':[(1,'exciton_analysis/transition_density'),(61,'exciton_analysis/descriptors'),(161,'exciton_analysis/ntos'),(196,'exciton_analysis/cubes'),(209,'exciton_analysis/dashboards'),(223,'exciton_analysis/configuration')],
 'part6_namd':[(1,'dynamics/pipeline'),(45,'dynamics/timesteps_cpa'),(154,'dynamics/model_selection'),(222,'dynamics/pme'),(287,'dynamics/fssh_edc'),(414,'dynamics/dish'),(523,'dynamics/comparison'),(576,'dynamics/nacs_tracking'),(622,'dynamics/decoherence'),(787,'spectroscopy/spectral_density'),(821,'recombination/overview'),(891,'recombination/radiative'),(919,'recombination/energy_gap_law'),(946,'recombination/trajectory_parameters'),(1013,'recombination/fcwd'),(1045,'recombination/traps'),(1059,'recombination/plqy'),(1073,'dynamics/analysis'),(1158,'spectroscopy/transient_absorption'),(1361,'dynamics/configuration')],
 'part7_examples':[(1,'workflows/index'),(15,'workflows/electronic_structure'),(77,'workflows/soc'),(124,'workflows/qp_edges'),(179,'workflows/four_frameworks'),(220,'workflows/exciton_analysis'),(275,'workflows/carrier_cooling')],
 'part8_auger':[(1,'recombination/auger_foundations'),(73,'recombination/auger_matrix_elements'),(263,'recombination/auger_screening'),(336,'recombination/auger_implementation'),(384,'recombination/auger_lineshapes'),(413,'recombination/auger_soc'),(432,'recombination/auger_workflows'),(475,'recombination/auger_ecsh'),(583,'recombination/auger_configuration')],
}
files=[]; units=[]
for path in sorted((root/'docs').rglob('*.rst')):
 rel=path.relative_to(root).as_posix(); data=path.read_bytes(); lines=data.splitlines(keepends=True)
 heads=[]
 for i,line in enumerate(lines):
  if i and re.fullmatch(rb'[=~^\-]{4,}',line.strip()) and lines[i-1].strip():
   heads.append((i,lines[i-1].decode().strip()))
 starts=sorted(set([1]+[x[0] for x in heads])); titles=dict(heads)
 for j,start in enumerate(starts):
  end=starts[j+1]-1 if j+1<len(starts) else len(lines)
  chunk=b''.join(lines[start-1:end]); parent=path.parent.name
  if parent in routes:
   dest=[dst for ln,dst in routes[parent] if ln<=start][-1]+'.rst'
  elif parent=='getting_started': dest='start/'+path.name
  elif parent=='api': dest='reference/api.rst'
  else: dest='index.rst'
  units.append(dict(source=rel,start_line=start,end_line=end,title=titles.get(start,'Preamble'),proposed_destination='docs/'+dest,sha256=hashlib.sha256(chunk).hexdigest(),bytes=len(chunk)))
 assert sum(x['bytes'] for x in units if x['source']==rel)==len(data)
 files.append(dict(path=rel,sha256=hashlib.sha256(data).hexdigest(),lines=len(lines),bytes=len(data)))
manifest=dict(note='Proposed routing only. RST units partition every source byte exactly once. Labels/directives must be checked before actual moves; preserve duplicate derivations with provenance, do not deduplicate during migration.',files=files,units=units)
(out/'documentation_inventory.json').write_text(json.dumps(manifest,indent=2)+'\n')
rows=['# Documentation migration inventory','',f'{len(files)} RST files; {len(units)} disjoint units; {sum(x["lines"] for x in files)} lines; every byte accounted for.','', '| Source | Lines | Current heading | Proposed destination |','|---|---:|---|---|']
for u in units:
 rows.append(f'| `{u["source"]}` | {u["start_line"]}–{u["end_line"]} | {u["title"].replace("|", " / ")} | `{u["proposed_destination"]}` |')
(out/'documentation_inventory.md').write_text('\n'.join(rows)+'\n')
# Exact signatures for topic pages, without importing extension modules.
sigs=[]
for path in sorted((root/'qdex').rglob('*.py')):
 source=path.read_text(); tree=ast.parse(source)
 for node in ast.walk(tree):
  if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef)):
   sigs.append(dict(file=str(path.relative_to(root)),line=node.lineno,name=node.name,signature=f'{node.name}({ast.unparse(node.args)})'))
(out/'api_signatures.json').write_text(json.dumps(sigs,indent=2)+'\n')
# Snapshot fingerprint includes dirty files, not just HEAD.
paths=list((root/'qdex').rglob('*.py'))+list((root/'docs').rglob('*.rst'))+list((root/'libint').glob('*.*'))
(out/'source_fingerprints.json').write_text(json.dumps({str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file()},indent=2)+'\n')
print(f'{len(files)} files, {len(units)} units, {sum(x["lines"] for x in files)} lines; {len(sigs)} signatures')
