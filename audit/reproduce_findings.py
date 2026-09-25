"""Small synthetic probes, not a physical CdSe nanocrystal benchmark.
Run with the repository root as cwd and its NumPy/SciPy Python environment.
"""
import contextlib
import io
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from qdex import hardness as h

def quiet(fn, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args, **kwargs)

coords = np.array([[0.,0.,0.],[2.6,0.,0.]])
common = dict(coords=coords, atom_symbols=['Cd','Se'], material_name='CDSE',
              C_occ_low=np.array([[1.],[0.]]), C_virt_low=np.array([[0.],[1.]]),
              eps_occ=np.array([-4.]), eps_virt=np.array([-1.4]),
              atom_ao_ranges=[(0,1),(1,2)], return_details=True)
out = {}
for e in [1., 2.4, 6.2, 10.]:
    shift, p = quiet(h.estimate_sgw_resta_qp_gap, **common, eps_out=e, penn_scaling=False)
    out[f'pure_resta_eps_{e}'] = {k:p[k] for k in ['total_scissor_ev','confinement_shift_ev','f_homo','f_lumo']}
qargs = {k:common[k] for k in ['coords','atom_symbols','material_name','atom_ao_ranges','return_details']}
shift,p,c,eps = quiet(h.estimate_qsgw_resta_qp_gap, **qargs, C=np.eye(2), S=np.eye(2),
                    eps=np.array([-4.,-1.4]), homo_index=0, eps_out=6.2, penn_scaling=False)
out['qsgw_zero_delta_w'] = dict(scissor=shift, expected_bulk_shift=1.27, eps=eps.tolist(), converged=p['qsgw_converged'])
out['dynamic_z'] = {str(s):h.compute_dynamic_z(s, 3., 6.2, 'CDSE') for s in [0., .2, .5, 1., 2.]}
out['dynamic_z_low_eps'] = h.compute_dynamic_z(2.,3.,1.5)
out['dynamic_z_documented_low_eps'] = 1/(1+2/np.sqrt(15**2/max(1,1.5-1)+3**2))
# A two-site delocalized frontier state changes only its global complex phase.
co=np.array([[1.],[1.]])/np.sqrt(2); cv=np.array([[1.],[-1.]])/np.sqrt(2)
a={**common,'C_occ_low':co,'C_virt_low':cv}
base=quiet(h.estimate_sgw_resta_qp_gap,**a,eps_out=2.4)[0]
try:
    phased=quiet(h.estimate_sgw_resta_qp_gap,**{**a,'C_occ_low':co*np.exp(1j*np.pi/4)},eps_out=2.4)[0]
    out['global_phase_probe']=dict(real_shift=base, phased_shift=phased)
except Exception as exc:
    out['global_phase_probe']=dict(real_shift=base,error=repr(exc))
# Zero-filled complement from a rectangular MO coefficient matrix.
try:
    shift,p,c,eps=quiet(h.estimate_qsgw_resta_qp_gap,**qargs,C=np.eye(3)[:,:2], S=np.eye(3),
       eps=np.array([-4.,-1.4]),homo_index=0,eps_out=6.2,penn_scaling=False)
    out['rectangular_mos']=dict(input_shape=[3,2],output_shape=list(c.shape),eps=eps.tolist())
except Exception as exc: out['rectangular_mos']=dict(error=repr(exc))
# DIM response normalization persists in a separated heteronuclear pair.
s,eta,eb,d=quiet(h.build_dim_screening_factors, np.array([[0.,0.,0.],[1e4,0.,0.]]),['Cd','Se'],'CDSE')
out['dim_separated_pair']=dict(eta=eta.tolist(),screening=float(s[0,1]),claimed_bulk_screening=1/eb)
# Conditional arithmetic at exactly R_model = 10 A, independent of any geometry input.
e=h.MATERIAL_DB['CDSE']; R=10.; ell=1.; p=2.
b=e[8]-e[7]; anchor=(e[13]-e[12])-(e[11]-e[10]); kv=11.52*(1-1/e[0]); ko=11.52*(1/2.4-1/e[0]); A=anchor-b-kv/(e[9]+ell)
out['cdse_anchor_at_model_radius_10A']=dict(bulk_pbe=e[7],bulk_gw=e[8],bulk_shift=b,anchor_shift=anchor,A=A,kappa_out=ko,scissor=b+ko/(R+ell)+A*(e[9]/R)**p)
out['cdse_consistency']=dict(volume_A3=4*np.pi*10**3/3,mu_for_15meV=.015*6.2**2/13.605693,
   radius_nm_for_15meV=0.052917721*6.2/(.015*6.2**2/13.605693),binding_eV_for_radius_5p6nm=14.3996/(2*6.2*56))
text=json.dumps(out,indent=2)
Path(__file__).with_name('probe_results.json').write_text(text+'\n')
print(text)
from qdex.exciton_hamiltonian import ExcitonHamiltonian
S=np.array([[1.,.25,.05],[.25,1.,.2],[.05,.2,1.]])
sv,su=np.linalg.eigh(S)
rng=np.random.default_rng(17)
u,_=np.linalg.qr(rng.normal(size=(3,3)))
C=(su/np.sqrt(sv))@su.T@u
baseargs=dict(C=C,eps=np.array([-4.,-1.4,0.]),overlap=S,atom_ao_ranges=[(0,1),(1,3)],
 homo_index=0,n_occ=1,n_virt=2,scissor_ev=1.27,gamma_qp=np.array([[2.,1.],[1.,3.]]),
 gamma_bse=np.array([[2.,1.],[1.,3.]]),gamma_bare=np.array([[4.,2.],[2.,5.]]),
 charge_type='lowdin',device='numpy',kernel_type='mnok')
for direct in [False,True]:
 full=quiet(ExcitonHamiltonian,**baseargs,include_exchange=direct,excitation_mode='bse')
 diag=quiet(ExcitonHamiltonian,**baseargs,include_exchange=direct,excitation_mode='diagonal_bse')
 A=np.column_stack([full.matvec(x) for x in np.eye(full.dim)])
 diag_e=diag.independent_transition_energies('diagonal_bse')[0]
 out[f'diagonal_full_consistency_direct_{direct}']=dict(full_diagonal=np.diag(A).tolist(),diagonal_mode=diag_e.tolist(),max_error=float(np.max(abs(np.diag(A)-diag_e))))
Path(__file__).with_name('probe_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k.startswith('diagonal_full')},indent=2))
