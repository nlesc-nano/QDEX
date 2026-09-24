import argparse
import atexit
import json
import numpy as np
import sys
import os
import time 
import gc
import platform
import yaml 

import libint_cpp

if __package__ is None or __package__ == "":
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from qdex.io_utils import (
    read_xyz, parse_basis, build_shell_dicts,
    count_ao_from_shells, build_atom_ao_ranges, read_mos_auto, read_mos_uks
)
from qdex.solver import ExcitonSolver
from qdex.constants import HA_TO_EV, BOHR_PER_ANG
from qdex.exciton_analysis import ExcitonAnalyzer, plot_analysis_summary
from qdex.integrals import compute_dipole_ao
from qdex.oscillator import compute_oscillator_strengths
from qdex.hardness import MATERIAL_DB, estimate_brus_qp_gap, estimate_gw_qp_gap
from qdex.orbital_analysis import (
    compute_spin_character, compute_uks_soc_spin_free_channels,
    compute_uks_spin_free_channels, format_uks_soc_spin_free_character,
    format_uks_spin_free_character, infer_reference_spin,
    print_orbital_summary, spin_multiplicity_name
)
from qdex.fuzzy_bands import run_fuzzy_bands_and_pdos, build_qp_energies, build_qp_energies_vacuum
from qdex.nto import run_nto_analysis
from qdex.profiler import ResourceTracker


class TeeStream:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)
            stream.flush()

    def flush(self):
        for stream in self.streams:
            stream.flush()


from qdex.device_utils import is_gpu, resolve_device

def transform_ao_operator(mu_ao, c_left, c_right, device="numpy"):
    """Transform an AO operator without silently reducing reference precision.

    Apple MPS does not implement float64/complex128 matrix multiplication.  The
    transition-dipole transform is small compared with the BSE solve, so MPS
    runs deliberately use the NumPy reference path here instead of changing
    the physical result to float32.
    """
    dtype = np.complex128 if any(np.iscomplexobj(x) for x in (mu_ao, c_left, c_right)) else np.float64
    mu = np.asarray(mu_ao, dtype=dtype)
    left = np.asarray(c_left, dtype=dtype)
    right = np.asarray(c_right, dtype=dtype)
    half = mu @ right

    if is_gpu(device):
        import torch
        dev = torch.device(device) if isinstance(device, str) else device
        return (
            torch.as_tensor(left, device=dev).conj().T
            @ torch.as_tensor(half, device=dev)
        ).cpu().numpy().astype(dtype, copy=False)

    return (left.conj().T @ half).astype(dtype, copy=False)


def setup_run_logging(log_file):
    if log_file in (None, "", "none", "None", False):
        return None

    log_handle = open(log_file, "w", encoding="utf-8")
    sys.stdout = TeeStream(sys.__stdout__, log_handle)
    sys.stderr = TeeStream(sys.__stderr__, log_handle)

    def close_log():
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        log_handle.close()

    atexit.register(close_log)
    print(f"  [Log] Writing full run log to {log_file}")
    return log_handle


def print_qp_provenance(details, dft_gap=None, target_qp_gap=None, output_file=None):
    if not details:
        return

    print("\n  [QP Provenance] Anchor-scaled PBE-to-QP model")
    print(f"    Material                 : {details['material']}")
    if "cluster_radius_ang" in details:
        print(f"    Cluster radius           : {details['cluster_radius_ang']:.3f} Å")
    print(f"    Bulk PBE -> GW gap       : {details['bulk_pbe_gap_ev']:.3f} -> {details['bulk_gw_gap_ev']:.3f} eV")
    print(f"    Bulk GW shift            : {details['bulk_gw_shift_ev']:+.3f} eV")
    if details.get("periodic_bulk_limit"):
        print("    Periodic mode            : using tabulated bulk GW-PBE scissor only")
    if details.get("has_monomer_anchor"):
        print(f"    Monomer anchor radius    : {details['monomer_radius_ang']:.3f} Å")
        print(f"    Monomer PBE -> GW gap    : {details['monomer_pbe_gap_ev']:.3f} -> {details['monomer_gw_gap_ev']:.3f} eV")
        print(f"    Anchor residual A        : {details['anchor_residual_ev']:+.3f} eV")
        print(f"    ell, p                   : {details['regularization_length_ang']:.3f} Å, {details['residual_power']:.3f}")
    if details.get("principal_extents_ang"):
        extents = ", ".join(f"{x:.3f}" for x in details["principal_extents_ang"])
        print(f"    Principal extents        : [{extents}] Å")
        print(f"    Anisotropy ratio         : {details['anisotropy_ratio']:.3f}")
    print(f"    Vacuum finite-size shift : {details['finite_size_shift_vacuum_ev']:+.3f} eV")
    print(f"    Solvent finite-size shift: {details['finite_size_shift_solvent_ev']:+.3f} eV")
    print(f"    Vacuum scissor           : {details['total_scissor_vacuum_ev']:+.3f} eV")
    print(f"    Solvent scissor used     : {details['total_scissor_solvent_ev']:+.3f} eV")
    if dft_gap is not None and target_qp_gap is not None:
        print(f"    Final gap                : {dft_gap:.3f} -> {target_qp_gap:.3f} eV")
    if output_file:
        print(f"    JSON                     : {output_file}")


def write_qp_provenance(details, dft_gap, target_qp_gap, scissor, args, filename="qp_provenance.json"):
    if not details:
        return None

    payload = dict(details)
    payload.update({
        "dft_gap_ev": float(dft_gap),
        "target_qp_gap_ev": float(target_qp_gap),
        "scissor_used_ev": float(scissor),
        "estimate_qp": bool(getattr(args, "estimate_qp", False)),
        "experimental_cohsex_tb_enabled": bool(getattr(args, "estimate_qp", False)),
        "use_cohsex_gap": bool(getattr(args, "use_cohsex_gap", False)),
        "notes": [
            "The scaled-GW hardness dictionary model is the recommended QP correction path.",
            "The COHSEX/TB Mulliken correction is experimental and not recommended for production QP provenance.",
        ],
    })

    with open(filename, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
        f.write("\n")
    return filename


def estimate_periodic_bulk_gw_scissor(material_name):
    if material_name is None:
        return None, None
    m_name = material_name.upper()
    entry = MATERIAL_DB.get(m_name)
    if entry is None or len(entry) < 9:
        return None, None

    gap_pbe_bulk = float(entry[7])
    gap_gw_bulk = float(entry[8])
    scissor = gap_gw_bulk - gap_pbe_bulk
    details = {
        "qp_model": "bulk_gw_hardness_dictionary",
        "periodic_bulk_limit": True,
        "material": m_name,
        "bulk_pbe_gap_ev": gap_pbe_bulk,
        "bulk_gw_gap_ev": gap_gw_bulk,
        "bulk_gw_shift_ev": scissor,
        "has_monomer_anchor": False,
        "finite_size_shift_vacuum_ev": 0.0,
        "finite_size_shift_solvent_ev": 0.0,
        "total_scissor_vacuum_ev": scissor,
        "total_scissor_solvent_ev": scissor,
    }
    return scissor, details

def run_solver_and_analysis(solver, coords_ang, syms, shells, mu_ia_x, mu_ia_y, mu_ia_z, dft_gap, scissor, confinement_energy, args, suffix="", soc_gap=None, soc_U=None, soc_E=None):
    """Encapsulates the solving, printing, analysis, and exporting."""
    spin = getattr(solver.ham, 'spin', 'singlet')
    if solver.soc_flag:
        label = "SOC"
    elif spin == 'uks_spin_preserving':
        label = "UKS SPIN-PRESERVING"
    elif spin == 'triplet':
        label = "TRIPLET (SPIN-FLIP)"
    else:
        label = "SPIN-FREE"
    print(f"\n===================================================")
    print(f" [ {label} ] EXCITON CALCULATION ")
    print(f"===================================================")

    start_solve = time.time()
    energies_ev, vectors = solver.solve(
        nroots=args.nroots, full_diag=args.full_diag, tol=args.tol,
        excitation_mode=args.excitation_mode,
    )
    
    if args.soc != 0.0:
        energies_ev = energies_ev - args.soc
        print(f"  [SOC Shift] Applied empirical energy shift: -{args.soc:.3f} eV")
    print(f"  {label} Solver converged in {time.time() - start_solve:2.2f} s")

    mu_ia = solver.ham.get_transition_dipoles(mu_ia_x, mu_ia_y, mu_ia_z)
    f_strengths = compute_oscillator_strengths(
        energies_ev, vectors, mu_ia,
        is_spinor=solver.soc_flag,
        spin_resolved=(spin == "uks_spin_preserving"),
    )

    def project_spinor_vec_to_spatial(vec):
        X_IA = np.zeros(solver.ham.dim_spinor_full, dtype=complex)
        if hasattr(solver.ham, 'valid_spinor_idx'):
            X_IA[solver.ham.valid_spinor_idx] = vec
        else:
            X_IA = vec
        X_IA = X_IA.reshape(solver.ham.n_occ_spinor, solver.ham.n_virt_spinor)

        if getattr(solver.ham, 'spin', 'singlet') == 'uks_spin_preserving':
            n_alpha = solver.ham.n_occ_act + solver.ham.n_virt_act
            occ_cols = slice(0, solver.ham.n_occ_spinor)
            virt_cols = slice(solver.ham.n_occ_spinor, solver.ham.n_occ_spinor + solver.ham.n_virt_spinor)
            U_occ_a = soc_U[0:solver.ham.n_occ_act, occ_cols]
            U_virt_a = soc_U[solver.ham.n_occ_act:n_alpha, virt_cols]
            U_occ_b = soc_U[n_alpha:n_alpha + solver.ham.n_occ_act_b, occ_cols]
            U_virt_b = soc_U[n_alpha + solver.ham.n_occ_act_b:n_alpha + solver.ham.n_occ_act_b + solver.ham.n_virt_act_b, virt_cols]
            X_a = U_occ_a.conj() @ X_IA @ U_virt_a.T
            X_b = U_occ_b.conj() @ X_IA @ U_virt_b.T
            return np.concatenate([
                np.abs(X_a[solver.ham.vi_a, solver.ham.va_a]),
                np.abs(X_b[solver.ham.vi_b, solver.ham.va_b]),
            ])

        n_mo = soc_U.shape[0] // 2
        n_occ_sp = solver.ham.n_occ_spinor
        U_occ_a = soc_U[:n_mo, :n_occ_sp]
        U_virt_a = soc_U[:n_mo, n_occ_sp:]
        U_occ_b = soc_U[n_mo:, :n_occ_sp]
        U_virt_b = soc_U[n_mo:, n_occ_sp:]
        X_ia_a = U_occ_a.conj() @ X_IA @ U_virt_a.T
        X_ia_b = U_occ_b.conj() @ X_IA @ U_virt_b.T
        # Do not add alpha and beta amplitudes before forming densities: valid
        # spin channels can cancel at amplitude level.  This norm preserves the
        # spin trace for the approximate spatial analysis path.
        return np.sqrt(
            np.abs(X_ia_a[solver.ham.valid_i, solver.ham.valid_a]) ** 2
            + np.abs(X_ia_b[solver.ham.valid_i, solver.ham.valid_a]) ** 2
        )

    print("\n" + "-"*60)
    print(f" SYSTEM ENERGY SUMMARY ({label})")
    print("-" * 60)
    print(f"  Raw DFT Gap           : {dft_gap:8.4f} eV")
    if solver.soc_flag and soc_gap is not None:
        print(f"  SOC Gap               : {soc_gap:8.4f} eV")
    print(f"  QP Correction (Shift) : {scissor:8.4f} eV")
    print(f"  Confinement Energy    : {confinement_energy:8.4f} eV")
    print(f"  Excitation Mode        : {args.excitation_mode}")
    if args.kernel != "resta":
        print(f"  Legacy Kernel Scaling  : {args.alpha:8.4f}")
    print("-" * 60)

    print("\n" + "="*172)
    print(f"{'State':>5} {'Energy':>10} {'Main Trans':>12} {'Weight':>8} {'f_osc':>10} | {'PR':>5} | {'D(eV)':>8} {'Kx(eV)':>8} {'-Kd(eV)':>8} | {'Spin-Free Character':>62}")
    print("-" * 172)

    n_print = min(100, len(energies_ev))

    for n in range(n_print):
        vec = vectors[:, n]
        hole_idx, elec_idx, weight = solver.main_transition(vec)
        
        if solver.soc_flag:
            n_occ_sp = solver.ham.n_occ_spinor
            abs_h, abs_e = hole_idx + 1, elec_idx + 1 + n_occ_sp
            h_lbl = f"spH" if abs_h == n_occ_sp else f"spH-{n_occ_sp - abs_h}"
            e_lbl = f"spL" if abs_e == n_occ_sp + 1 else f"spL+{abs_e - (n_occ_sp + 1)}"
            trans_str = f"{h_lbl}->{e_lbl}"
        elif getattr(solver.ham, 'spin', 'singlet') == 'uks_spin_preserving':
            # Determine which channel the dominant transition belongs to
            idx_dom = np.argmax(np.abs(vectors[:, n]))
            if idx_dom < solver.ham.dim_a:
                # Alpha channel
                abs_h = (solver.ham.homo_index_alpha - solver.ham.n_occ_act + 1) + hole_idx
                abs_e = (solver.ham.homo_index_alpha + 1) + elec_idx
                trans_str = f"{abs_h:3d}->{abs_e:3d}(α)"
            else:
                # Beta channel
                abs_h = (solver.ham.homo_index_beta - solver.ham.n_occ_act_b + 1) + hole_idx
                abs_e = (solver.ham.homo_index_beta + 1) + elec_idx
                trans_str = f"{abs_h:3d}->{abs_e:3d}(β)"
        else:
            abs_h = (solver.homo_index - solver.ham.n_occ_act + 1) + hole_idx
            abs_e = (solver.homo_index + 1) + elec_idx
            trans_str = f"{abs_h:3d}->{abs_e:3d}"

        vec_conj = vec.conj() if np.iscomplexobj(vec) else vec
        dE_val, Kx_val, minus_Kd_val = solver.expectation_components(vec)

        pr = 1.0 / np.sum(np.abs(vec)**4)

        if solver.soc_flag and getattr(solver.ham, 'spin', 'singlet') == 'uks_spin_preserving':
            character = compute_uks_soc_spin_free_channels(
                vec, solver.ham, soc_U,
                getattr(args, "n_alpha_ref", 0.0),
                getattr(args, "n_beta_ref", 0.0)
            )
            spin_str = format_uks_soc_spin_free_character(character)
        elif solver.soc_flag:
            soc_mask = None
            if args.e_thresh is not None and soc_E is not None:
                n_occ_sp = solver.ham.n_occ_spinor
                n_virt_sp = solver.ham.n_virt_spinor
                soc_occ_E = soc_E[:n_occ_sp]
                soc_virt_E = soc_E[n_occ_sp:n_occ_sp + n_virt_sp]
                soc_mask = (soc_virt_E.reshape(1, -1) - soc_occ_E.reshape(-1, 1)) <= args.e_thresh
                
            s_pct, t_pct = compute_spin_character(vec, soc_U, solver.ham.n_occ_spinor, solver.ham.n_virt_spinor, valid_mask=soc_mask)
            spin_str = f"{s_pct:5.1f}% S / {t_pct:5.1f}% T"
        elif getattr(solver.ham, 'spin', 'singlet') == 'uks_spin_preserving':
            character = compute_uks_spin_free_channels(
                vec, solver.ham,
                getattr(args, "n_alpha_ref", 0.0),
                getattr(args, "n_beta_ref", 0.0)
            )
            spin_str = format_uks_spin_free_character(character)
        else:
            spin_str = "100.0% S /   0.0% T"

        print(f"{n+1:5d} {energies_ev[n]:10.4f}  {trans_str:>12}  {weight**2:8.3f}  {f_strengths[n]:10.5f} | {pr:5.1f} | {dE_val:8.4f} {Kx_val:8.4f} {minus_Kd_val:8.4f} | {spin_str:>62}")

    if len(energies_ev) > 100: print(f" ... {len(energies_ev) - 100} additional states computed (output truncated) ...")
    print("="*172)

    print(f"\n--- Performing Dreuw/Plasser Analysis ({label}) ---")
    analyzer = ExcitonAnalyzer(solver, np.array(coords_ang), syms)
    analysis_results = []
    analysis_t0 = time.perf_counter()

    print(f"{'State':>5} {'Energy':>8} {'f_osc':>8} | {'PR':>5} {'d_eh~(A)':>8} {'d_CT~(A)':>8} {'sig_h~':>7} {'sig_e~':>7} | {'Type':>8}")
    print("-" * 95)

    for n in range(n_print):
        vec = vectors[:, n]
        
        # Project spinor back to spatial for Dreuw-Plasser natively
        if solver.soc_flag:
            vec_spatial = project_spinor_vec_to_spatial(vec)
        else:
            vec_spatial = vec

        res = analyzer.analyze_state(vec_spatial, energies_ev[n], f_strengths[n])
        analysis_results.append(res)
        
        ct_ratio = res['CT_Character']
        if ct_ratio > 0.6: ex_type = "CT"
        elif res['d_eh'] < 3.0 and ct_ratio < 0.2: ex_type = "Frenkel"
        else: ex_type = "Wannier"

        print(f"{n+1:5d} {res['energy']:8.3f} {res['f_osc']:8.4f} | {res['PR']:5.1f} {res['d_eh']:7.2f} {res['d_CT']:7.2f} {res['sigma_h']:6.1f} {res['sigma_e']:6.1f} | {ex_type:>8}")

    print(
        f"  [Analyzer] Completed {len(analysis_results)} coherent-state analyses "
        f"in {time.perf_counter() - analysis_t0:.2f} s"
    )

    if getattr(args, 'nto', False):
        run_nto_analysis(
            solver=solver,
            vectors=vectors,
            energies_ev=energies_ev,
            f_strengths=f_strengths,
            coords=np.array(coords_ang),
            symbols=syms,
            mu_ia=mu_ia,
            args=args,
            suffix=suffix,
            soc_U=soc_U if solver.soc_flag else None,
        )

    # =========================================================================
    # EXCITON CUBE GENERATION (MOs are now generated globally before this)
    # =========================================================================
    if getattr(args, 'cube', False):
        from qdex.exciton_cube import generate_cubes
        
        bse_states_arg = getattr(args, 'bse_states', None)
        if bse_states_arg:
            top_indices = [i - 1 for i in bse_states_arg if (i - 1) < len(f_strengths)]
            n_bse = len(top_indices)
        else:
            n_bse = min(getattr(args, 'nbse', 3), len(f_strengths))
            top_indices = np.argsort(f_strengths)[-n_bse:][::-1]
        
        bse_states = {}
        for idx in top_indices:
            raw_vec = vectors[:, idx]
            if solver.soc_flag:
                cube_vec = project_spinor_vec_to_spatial(raw_vec)
            else:
                cube_vec = raw_vec
            
            bse_states[f"exciton_state_{idx + 1}{suffix}"] = cube_vec
            
        print(f"\n--- Generating Cubes ({n_bse} Excitons) ---")
        use_cpp_writer = not getattr(args, 'disable_cpp_cube', False)
        
        generate_cubes(
            solver=solver, 
            bse_states_dict=bse_states, 
            mo_list=[],     # MOs are generated earlier
            spinor_list=[], # Spinors are generated earlier
            soc_U=soc_U if solver.soc_flag else None,
            shells=shells, symbols=syms, coords=coords_ang, 
            spacing_ang=args.cube_spacing, nthreads=args.nthreads,
            use_cpp=use_cpp_writer
        )

    if args.write_csv:
        import csv
        csv_file = f"exciton_results{suffix}.csv"
        n_to_write = min(args.csv_roots, len(analysis_results))
        with open(csv_file, mode='w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Time", "State", "Energy_eV", "f_osc", "mu_x", "mu_y", "mu_z", "PR", "d_eh_A", "d_CT_A", "sigma_h_A", "sigma_e_A", "Type"])
            for n in range(n_to_write):
                res = analysis_results[n]
                mu_state = np.sum(vectors[:, n][:, None] * mu_ia, axis=0)
                ct_ratio = res['CT_Character']
                ex_type = "CT" if ct_ratio > 0.6 else "Frenkel" if res['d_eh'] < 3.0 and ct_ratio < 0.2 else "Wannier"
                writer.writerow([
                    args.time, n + 1, f"{energies_ev[n]:.6f}", f"{f_strengths[n]:.6e}", f"{np.real(mu_state[0]):.6f}", f"{np.real(mu_state[1]):.6f}", f"{np.real(mu_state[2]):.6f}",
                    f"{res['PR']:.3f}", f"{res['d_eh']:.4f}", f"{res['d_CT']:.4f}", f"{res['sigma_h']:.4f}", f"{res['sigma_e']:.4f}", ex_type
                ])

    if args.save_xia:
        npz_filename = f"xia{suffix}_{args.time:08.2f}fs.npz"
        n_to_write = min(args.csv_roots, len(energies_ev))
        np.savez_compressed(
            npz_filename, time=args.time, energies=energies_ev[:n_to_write], X_ia=vectors[:, :n_to_write], 
            valid_i=solver.ham.valid_i, valid_a=solver.ham.valid_a, homo_index=solver.homo_index, n_occ=solver.n_occ
        )
        print(f"  Saved X_ia coefficients to {npz_filename}")

    if args.broadening != "none":
        from qdex.spectrum import generate_spectrum, plot_spectrum
        e_min, e_max = max(0.0, np.min(energies_ev) - 2.5), np.max(energies_ev) + 2.5
        x_grid, y_grid = generate_spectrum(energies_ev, f_strengths, e_min=e_min, e_max=e_max, sigma=args.sigma, profile=args.broadening)
        spec_file = f"spectrum{suffix}.dat"
        np.savetxt(spec_file, np.column_stack((x_grid, y_grid)), header=f"Energy(eV) Intensity(arb.u.) | {args.broadening}, sigma={args.sigma}")
        if args.plot or args.show:
            plot_file = f"spectrum{suffix}.png" if args.plot else None
            plot_spectrum(x_grid, y_grid, energies_ev, f_strengths, filename=plot_file, show=args.show)

    if args.plot or args.show:
        ref_gap = soc_gap + scissor if (solver.soc_flag and soc_gap is not None) else args.qp_gap_num
        metrics = {
            "dft_gap": dft_gap, "qp_correction": scissor, "confinement_energy": confinement_energy, 
            "binding_energy": ref_gap - energies_ev[0] if len(energies_ev) > 0 else 0.0, 
            "first_exc_energy": energies_ev[0] if len(energies_ev) > 0 else 0.0, "is_soc": solver.soc_flag
        }
        if solver.soc_flag and soc_gap is not None: metrics["soc_gap"] = soc_gap
        plot_file = f"exciton_analysis{suffix}.html" if args.plot else None
        plot_analysis_summary(analysis_results, physics_metrics=metrics, filename=plot_file, show=args.show, broadening=args.broadening, sigma=args.sigma)

    if getattr(args, "auger", False):
        from qdex.auger import calculate_auger_rates
        eps_eval = solver.eps.copy()
        if scissor is not None and scissor != 0.0:
            eps_eval[solver.homo_index + 1:] += scissor
        if solver.soc_flag and soc_E is not None and soc_U is not None:
            k = soc_U.shape[0] // 2
            act_start = solver.homo_index - solver.n_occ + 1
            act_end = solver.homo_index + 1 + solver.n_virt
            C_act = solver.C[:, act_start:act_end]
            U_spinor_alpha = C_act @ soc_U[:k, :]
            U_spinor_beta = C_act @ soc_U[k:, :]
            eps_eval = soc_E
            homo_idx_eval = (solver.n_occ * 2) - 1
            is_spinor = True
        else:
            U_spinor_alpha = None
            U_spinor_beta = None
            homo_idx_eval = solver.homo_index
            is_spinor = False

        calculate_auger_rates(
            C=solver.C,
            eps=eps_eval,
            S=solver.overlap,
            atom_ao_ranges=solver.atom_ao_ranges,
            coords=getattr(solver, "coords", coords_ang),
            atom_symbols=getattr(solver, "atom_symbols", syms),
            homo_idx=homo_idx_eval,
            W_resta=getattr(solver.ham, "W_resta", None),
            material_name=getattr(args, "material", "DEFAULT"),
            eps_out=getattr(args, "eps_out", 2.0),
            sigma_ev=getattr(args, "auger_sigma", 0.05),
            broadening_mode=getattr(args, "auger_lineshape", "gaussian"),
            channel=getattr(args, "auger_channel", "all"),
            n_initial_elec=getattr(args, "auger_states", 1),
            n_initial_hole=getattr(args, "auger_states", 1),
            spinor=is_spinor,
            U_spinor_alpha=U_spinor_alpha,
            U_spinor_beta=U_spinor_beta,
            eps_eff=getattr(args, "auger_eps_eff", None),
            verbose=True,
        )


def validate_args(args, parser):
    if args.e_thresh is not None and (args.nhomos is not None or args.nlumos is not None):
        parser.error("Conflict: Please specify either energy threshold selection (--e_thresh) OR direct active space count (--nhomos / --nlumos), but not both.")
        
    if args.soc_flag and args.gth_file is None:
        parser.error("Validation error: --soc_flag requires a GTH pseudopotential file via --gth_file.")
        
    if args.run_fuzzy and args.cif is None:
        parser.error("Validation error: --run_fuzzy requires a CIF crystal structure file via --cif.")

    if args.beta > 0.0:
        parser.error("Validation error: beta > 0 is unavailable until a validated on-site U table is provided.")

    if args.qp_regularization_length < 0.0:
        parser.error("Validation error: --qp-regularization-length must be non-negative.")

    if args.qp_residual_power <= 1.0:
        parser.error("Validation error: --qp-residual-power must be greater than 1.")

    if getattr(args, "periodic_enabled", False):
        lattice_vectors = getattr(args, "lattice_vectors", None)
        if lattice_vectors is None:
            parser.error("Validation error: periodic.enabled requires periodic.lattice_vectors in the YAML config.")
        lattice = np.asarray(lattice_vectors, dtype=float)
        if lattice.shape != (3, 3):
            parser.error("Validation error: periodic.lattice_vectors must be a 3x3 list of vectors in angstrom.")


def _apply_config(args, config_data):
    for section, parameters in config_data.items():
        if section == "periodic" and isinstance(parameters, dict):
            setattr(args, "periodic_enabled", bool(parameters.get("enabled", False)))
            if "lattice_vectors" in parameters:
                setattr(args, "lattice_vectors", parameters["lattice_vectors"])
            if "overlap_cutoff" in parameters:
                setattr(args, "overlap_cutoff", parameters["overlap_cutoff"])
            continue

        if section == "namd" and isinstance(parameters, dict):
            setattr(args, "namd_cfg", parameters)
            continue

        if section == "auger" and isinstance(parameters, dict):
            setattr(args, "auger", bool(parameters.get("run", parameters.get("enabled", True))))
            if "sigma" in parameters:
                setattr(args, "auger_sigma", float(parameters["sigma"]))
            if "channel" in parameters:
                setattr(args, "auger_channel", str(parameters["channel"]))
            if "n_initial_states" in parameters:
                setattr(args, "auger_states", int(parameters["n_initial_states"]))
            if "lineshape" in parameters:
                setattr(args, "auger_lineshape", str(parameters["lineshape"]))
            if "eps_eff" in parameters and parameters["eps_eff"] is not None:
                setattr(args, "auger_eps_eff", float(parameters["eps_eff"]))
            continue

        if isinstance(parameters, dict):
            for key, value in parameters.items():
                if not hasattr(args, key):
                    raise ValueError(f"Unknown YAML key '{section}.{key}'")
                setattr(args, key, value)
        else:
            if not hasattr(args, section):
                raise ValueError(f"Unknown YAML key '{section}'")
            setattr(args, section, parameters)


def _select_soc_window_indices(eps_shifted, homo_index, soc_window):
    eps_shifted = np.asarray(eps_shifted)
    if soc_window is None or soc_window <= 0:
        return np.arange(len(eps_shifted), dtype=int)

    mask_idx = np.where(np.abs(eps_shifted) <= float(soc_window))[0]
    if mask_idx.size == 0:
        start = max(0, homo_index)
        stop = min(len(eps_shifted), homo_index + 2)
    else:
        start = min(int(mask_idx[0]), max(0, homo_index))
        stop = max(int(mask_idx[-1]) + 1, min(len(eps_shifted), homo_index + 2))

    return np.arange(start, stop, dtype=int)
        

def main():
    parser = argparse.ArgumentParser(description="QDEX - Quantum Dot Excitations & Dynamics exciton solver")

    parser.add_argument("--config", type=str, help="Path to a YAML configuration file.")
    parser.add_argument("--mo_file")
    parser.add_argument("--xyz") 
    parser.add_argument("--basis_txt")
    parser.add_argument("--basis_name")
    parser.add_argument(
        "--cache-mos", dest="cache_mos", action="store_true",
        help="Cache parsed text MOs as a validated uncompressed binary NPZ sidecar for faster repeated runs.",
    )

    parser.add_argument("--n-occ", type=int, default=50)
    parser.add_argument("--n-virt", type=int, default=50)
    parser.add_argument("--e_thresh", type=float, default=None)
    parser.add_argument("--f_thresh", type=float, default=0.0)

    parser.add_argument("--qp_gap", type=str, default="brus")
    parser.add_argument("--soc", type=float, default=0.0)
    parser.add_argument("--soc_flag", action="store_true")
    parser.add_argument("--gth_file", type=str, default=None)

    parser.add_argument("--kernel", choices=["bse", "resta"], default="bse")
    parser.add_argument("--alpha", type=float, default=1.0, help="Scaling for the legacy non-RESTA kernel; ignored by RESTA.")
    parser.add_argument("--beta", type=float, default=0.0, help="Reserved on-site stiffening parameter; beta > 0 is currently rejected.")
    parser.add_argument("--exchange", action="store_true", default=None, help="Deprecated alias for --include-direct-eh.")
    direct_group = parser.add_mutually_exclusive_group()
    direct_group.add_argument("--include-direct-eh", dest="include_direct_eh", action="store_true", default=None,
                              help="Include the Resta-screened attractive electron-hole direct term (default).")
    direct_group.add_argument("--no-direct-eh", dest="include_direct_eh", action="store_false",
                              help="Disable the attractive electron-hole direct term.")
    parser.add_argument("--estimate_qp", action="store_true", help="Compute G0W0-lite Quasiparticle corrections via COHSEX")
    parser.add_argument("--use_cohsex_gap", action="store_true", help="Override the tabulated GW gap with the pure COHSEX computed gap")
    parser.add_argument("--vxc_ao", type=str, default=None, help="Path to cleaned CP2K AO-basis Vxc matrix text file")
    parser.add_argument("--material", type=str, default="DEFAULT")
    parser.add_argument("--eps-out", type=float, default=2.0)
    parser.add_argument("--qp-regularization-length", dest="qp_regularization_length", type=float, default=1.0,
                        help="Regularization length ell in angstrom for the anchor-scaled QP model.")
    parser.add_argument("--qp-residual-power", dest="qp_residual_power", type=float, default=2.0,
                        help="Power p > 1 controlling decay of the anchor residual.")
    parser.add_argument("--qp-strict", action="store_true",
                        help="Reject clusters smaller than the finite-size anchor instead of clamping to it.")

    parser.add_argument("--broadening", choices=["gaussian", "lorentzian", "none"], default="gaussian")
    parser.add_argument("--sigma", type=float, default=0.1)
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--show", action="store_true")
    
    # --- Cube Arguments ---
    parser.add_argument("--cube", action="store_true")
    parser.add_argument("--cube-spacing", type=float, default=0.5)
    parser.add_argument("--disable_cpp_cube", action="store_true")
    parser.add_argument("--cube-nhomos", type=int, default=2, dest="cube_nhomos", help="Number of HOMO states (spatial or spinor) to export")
    parser.add_argument("--cube-nlumos", type=int, default=2, dest="cube_nlumos", help="Number of LUMO states (spatial or spinor) to export")
    parser.add_argument("--nbse", type=int, default=3, help="Number of Top Exciton states to export if bse_states is not provided")
    parser.add_argument("--bse_states", type=int, nargs='+', help="Specific exciton states to plot (1-indexed, e.g., 2 9 19)")
    
    # --- BSE Active Space Arguments ---
    parser.add_argument("--nhomos", type=int, default=None, help="Number of occupied MOs to include in the BSE active space")
    parser.add_argument("--nlumos", type=int, default=None, help="Number of virtual MOs to include in the BSE active space")
    parser.add_argument("--charge_type", choices=["mulliken", "lowdin"], default="mulliken", help="Method to compute transition charges (Mulliken or Lowdin).")

    # --- Triplet State / UKS Arguments ---
    parser.add_argument("--triplet", action="store_true", help="Perform triplet state BSE calculations")
    parser.add_argument("--mo_file_beta", type=str, default=None, help="Path to molecular orbitals file for beta spin channel (UKS triplet)")
    # ----------------------
    
    parser.add_argument("--write-csv", action="store_true")
    parser.add_argument("--csv-roots", type=int, default=10)
    parser.add_argument("--time", type=float, default=0.0)
    parser.add_argument("--save-xia", action="store_true")
    parser.add_argument("--log_file", type=str, default="minibse.log", help="Write terminal output to this log file as well as stdout; use 'none' to disable.")

    parser.add_argument("--nto", action="store_true", help="Run Natural Transition Orbital analysis after solving excitons.")
    parser.add_argument("--nto-states", type=int, nargs='+', help="Specific exciton states for NTO analysis, 1-indexed.")
    parser.add_argument("--nto-top", type=int, default=3, help="Number of dominant NTO pairs to print per state.")
    parser.add_argument("--nto-csv", action="store_true", help="Write detailed NTO descriptors to nto_results*.csv.")
 
    parser.add_argument("--nroots", type=int, default=10)
    parser.add_argument("--full-diag", action="store_true")
    parser.add_argument(
        "--excitation-mode",
        choices=["bse", "independent_dft", "independent_qp", "diagonal_bse"],
        default="bse",
        help=("Excitation model: diagonalize BSE/TDA, or use uncoupled transitions with "
              "DFT gaps, QP gaps, or QP gaps plus diagonal Kx/Kd corrections."),
    )
    parser.add_argument("--tol", type=float, default=1e-5)
    parser.add_argument("--orthonormality-tol", type=float, default=1e-5,
                        help="Abort when max|C^dagger S C-I| exceeds this tolerance.")
    parser.add_argument("--nthreads", type=int, default=1)
    parser.add_argument("--device", choices=["auto", "cpu", "cuda", "mps", "numpy"], default="auto")

    # Auger arguments
    parser.add_argument("--auger", action="store_true", help="Compute non-radiative Auger recombination rates and biexciton lifetimes.")
    parser.add_argument("--auger-sigma", type=float, default=0.05, help="Energy conservation Gaussian broadening for Auger recombination in eV (default: 0.05).")
    parser.add_argument("--auger-channel", choices=["all", "eeh", "hhe"], default="all", help="Auger channel to compute: all, eeh, or hhe (default: all).")
    parser.add_argument("--auger-states", type=int, default=1, help="Number of band-edge frontier states to consider as initial carriers (default: 1).")
    parser.add_argument("--auger-lineshape", choices=["gaussian", "fcwd"], default="gaussian", help="Energy conservation line shape for Auger rates (default: gaussian).")
    parser.add_argument("--auger-eps-eff", type=float, default=None, help="Effective dielectric constant for dynamic screening at energy transfer hbar*omega=Eg (e.g. 1.8; default: Resta eps_inf).")

    # Fuzzy arguments
    parser.add_argument("--run_fuzzy", action="store_true")
    parser.add_argument("--cif", type=str)
    parser.add_argument("--soc_window", type=float, default=10.0)
    parser.add_argument("--pdos_atoms", type=str, nargs='+')
    parser.add_argument("--coop_pairs", type=str, nargs='+')
    parser.add_argument("--population_print_range", type=int, default=15)
    parser.add_argument("--fuzzy_sigma", type=float, default=0.03)
    parser.add_argument("--pdos_sigma", type=float, default=0.10)
    parser.add_argument("--ewin", type=float, nargs=2, default=[-5.0, 5.0])
    parser.add_argument("--fold_to_bz", action="store_true", help="Fold fuzzy-band weights into the first Brillouin zone by summing reciprocal replicas.")
    parser.add_argument("--g_shell", type=int, default=0, help="Reciprocal-vector shell for folded fuzzy-band weights; 0, 1, and 2 use 1, 27, and 125 replicas.")
    parser.add_argument("--dashboard_energy_mode", choices=["dft", "qp", "both"], default="dft", help="Generate fuzzy dashboards on DFT, QP-corrected, or both energy axes.")
    parser.add_argument("--qp_energy_reference", choices=["vacuum", "fermi"], default="vacuum", help="Energy reference for QP fuzzy dashboards.")

    # Periodic/Gamma-only arguments. YAML is the preferred interface.
    parser.add_argument("--periodic", action="store_true", dest="periodic_enabled", help="Use Gamma-only periodic AO overlap.")
    parser.add_argument("--lattice_vectors", type=float, nargs=9, default=None, help="3x3 lattice vectors in angstrom, row-major.")
    parser.add_argument("--overlap_cutoff", type=float, default=-1.0, help="Periodic overlap image cutoff in angstrom; <=0 uses minimum image only.")

    # NAMD arguments
    parser.add_argument("--namd-precompute", action="store_true", help="Run NAMD precomputation (cross-overlaps, tracking, caching).")
    parser.add_argument("--namd-run", action="store_true", help="Run NAMD carrier cooling simulation from precomputed data.")
    parser.add_argument("--namd-compact", type=str, nargs="?", const="default", default=None, help="Compact precomputed NAMD directory (compresses and removes redundant arrays).")
    parser.add_argument("--namd-decoherence", type=str, nargs="?", const="default", default=None, help="Compute state-pair pure-dephasing decoherence times from trajectory fluctuations.")
    parser.add_argument("--namd-soc", action="store_true", help="Enable Spin-Orbit Coupling (SOC) for NAMD precomputation.")
    parser.add_argument("--namd-ta", action="store_true", help="Compute ultrafast pump-probe transient absorption (TA) spectra from NAMD dynamics.")
    parser.add_argument("--namd-ta-sigma", type=float, default=0.03, help="Gaussian line broadening in eV for transient absorption probe spectra (default: 0.03).")
    parser.add_argument("--namd-ta-plot", action="store_true", help="Generate 2D false-color TA map and 1S bleach rise kinetics plot.")
    parser.add_argument("--namd-ecsh-auger", action="store_true", help="Enable Energy-Conserving Surface Hopping (ECSH) for Auger processes in NAMD.")
    parser.add_argument("--namd-ecsh-window", type=float, default=None, help="Resonance energy window in eV for ECSH Auger transitions (default: k_B * T).")
    parser.add_argument("--namd-trajectory-loops", type=int, default=1, help="Number of times to loop precomputed MD trajectory to reach long Auger timescales (default: 1).")
    parser.add_argument("--namd-biexciton", action="store_true", help="Initialize NAMD from a biexciton state (XX) to simulate Auger annihilation dynamics.")
    parser.add_argument("--namd-method", type=str, default=None, choices=["cpa_fssh", "fssh", "dish", "cpa_fssh_edc", "cpa_fssh_gdc", "fssh_gdc", "pme", "master_equation"], help="NAMD simulation method (cpa_fssh, cpa_fssh_gdc, dish, pme).")
    parser.add_argument("--namd-tr-sd", action="store_true", help="Compute time-resolved vibrational action spectrum / dynamical phonon spectrogram J(omega, t).")
    parser.add_argument("--namd-tr-sd-plot", action="store_true", help="Generate multi-panel 2D false-color heatmap of frequency vs time.")
    parser.add_argument("--namd-tr-sd-wmax", type=float, default=400.0, help="Maximum phonon frequency in cm^-1 for time-resolved spectral density (default: 400.0).")
    parser.add_argument("--namd-tr-sd-sigma", type=float, default=15.0, help="Gaussian time broadening in fs for discrete surface hops in J(omega, t) (default: 15.0).")
    parser.add_argument("--namd-2d-map", action="store_true", help="Compute 2D non-adiabatic vibronic action map S(omega_acc, Omega_prom).")
    parser.add_argument("--namd-2d-plot", action="store_true", help="Generate publication-grade 2D contour map of promoting vs accepting modes.")
    parser.add_argument("--namd-initial-conditions", type=str, default=None, choices=["single", "multiple", "multi"], help="Initial condition sampling mode: 'single' (t0 = 0) or 'multiple' (automated multi-origin ensemble).")
    parser.add_argument("--namd-multi-init", action="store_true", help="Enable automated multi-origin ensemble sampling across the MD trajectory.")
    parser.add_argument("--namd-origins", type=int, default=None, help="Number of independent initial condition origins for multi-origin NAMD (default: auto-calibrated).")
    parser.add_argument("--namd-window-fs", type=float, default=None, help="Simulation window duration in fs for each multi-origin sub-trajectory (default: auto-calibrated).")

    args = parser.parse_args()

    config_path = args.config
    config_data = {}
    if args.config:
        with open(args.config, 'r') as f:
            config_data = yaml.safe_load(f) or {}

        _apply_config(args, config_data)

    if getattr(args, "namd_soc", False):
        config_data.setdefault("physics", {})["soc"] = True
    if getattr(args, "gth_file", None):
        config_data.setdefault("system", {})["gth_file"] = args.gth_file

    if getattr(args, "namd_method", None):
        config_data.setdefault("namd", {}).setdefault("dynamics", {})["method"] = args.namd_method

    if getattr(args, "namd_multi_init", False):
        config_data.setdefault("namd", {}).setdefault("dynamics", {})["initial_conditions"] = "multiple"
    elif getattr(args, "namd_initial_conditions", None):
        config_data.setdefault("namd", {}).setdefault("dynamics", {})["initial_conditions"] = args.namd_initial_conditions

    if getattr(args, "namd_origins", None) is not None:
        config_data.setdefault("namd", {}).setdefault("dynamics", {})["n_origins"] = args.namd_origins
    if getattr(args, "namd_window-fs", None) is not None:
        config_data.setdefault("namd", {}).setdefault("dynamics", {})["window_fs"] = args.namd_window_fs
    elif getattr(args, "namd_window_fs", None) is not None:
        config_data.setdefault("namd", {}).setdefault("dynamics", {})["window_fs"] = args.namd_window_fs

    if getattr(args, "namd_tr_sd", False) or getattr(args, "namd_tr_sd_plot", False):
        tr_dict = config_data.setdefault("namd", {}).setdefault("time_resolved_spectral_density", {})
        tr_dict["run"] = True
        if getattr(args, "namd_tr_sd_plot", False):
            tr_dict["plot"] = True
        if getattr(args, "namd_tr_sd_wmax", None) is not None:
            tr_dict["w_max_cm"] = args.namd_tr_sd_wmax
        if getattr(args, "namd_tr_sd_sigma", None) is not None:
            tr_dict["sigma_t_fs"] = args.namd_tr_sd_sigma
        if getattr(args, "namd_2d_map", False):
            tr_dict["compute_2d_vibronic"] = True
        if getattr(args, "namd_2d_plot", False):
            tr_dict["plot_2d"] = True

    if getattr(args, "namd_ta", False):
        ta_dict = config_data.setdefault("namd", {}).setdefault("transient_absorption", {})
        ta_dict["run"] = True
        if getattr(args, "namd_ta_sigma", None) is not None:
            ta_dict["sigma"] = args.namd_ta_sigma
        if getattr(args, "namd_ta_plot", False):
            ta_dict["plot"] = True

    if getattr(args, "namd_ecsh_auger", False) or getattr(args, "namd_biexciton", False):
        dyn_dict = config_data.setdefault("namd", {}).setdefault("dynamics", {})
        if getattr(args, "namd_ecsh_auger", False):
            dyn_dict["ecsh_auger"] = True
        if getattr(args, "namd_ecsh_window", None) is not None:
            dyn_dict["ecsh_window_ev"] = args.namd_ecsh_window
        if getattr(args, "namd_trajectory_loops", 1) > 1:
            dyn_dict["trajectory_loops"] = args.namd_trajectory_loops
        if getattr(args, "namd_biexciton", False):
            dyn_dict["initial_state"] = "biexciton"
    elif getattr(args, "namd_trajectory_loops", 1) > 1:
        dyn_dict = config_data.setdefault("namd", {}).setdefault("dynamics", {})
        dyn_dict["trajectory_loops"] = args.namd_trajectory_loops

    if getattr(args, "namd_compact", None) is not None:
        from qdex.namd import compact_precomputed_data
        compact_dir = args.namd_compact
        if compact_dir == "default":
            compact_dir = config_data.get("namd", {}).get("storage", {}).get("precompute_dir", "namd_precomputed")
        compact_precomputed_data(compact_dir)
        return

    if getattr(args, "namd_decoherence", None) is not None:
        from qdex.namd.precompute import compute_trajectory_decoherence_times
        dec_dir = args.namd_decoherence
        if dec_dir == "default":
            dec_dir = config_data.get("namd", {}).get("storage", {}).get("precompute_dir", "namd_precomputed")
        compute_trajectory_decoherence_times(dec_dir)
        return

    if getattr(args, "namd_precompute", False):
        from qdex.namd import precompute_namd_data
        setup_run_logging(getattr(args, "log_file", "minibse.log"))
        precompute_namd_data(config_data)
        return

    if getattr(args, "namd_run", False):
        from qdex.namd import run_namd_dynamics
        setup_run_logging(getattr(args, "log_file", "minibse.log"))
        run_namd_dynamics(config_data)
        return

    if getattr(args, "lattice_vectors", None) is not None:
        args.lattice_vectors = np.asarray(args.lattice_vectors, dtype=float).reshape(3, 3).tolist()

    if getattr(args, "include_direct_eh", None) is None:
        args.include_direct_eh = True if args.exchange is None else bool(args.exchange)
    if args.exchange is not None:
        print("  [Deprecated] 'exchange' now maps to include_direct_eh; use include_direct_eh instead.")

    setup_run_logging(getattr(args, "log_file", "minibse.log"))
    validate_args(args, parser)
    if config_path:
        print(f"Loading configuration from {config_path}...")

    required_args = ['mo_file', 'xyz', 'basis_txt', 'basis_name', 'qp_gap']
    missing = [arg for arg in required_args if getattr(args, arg) is None]
    if missing: parser.error(f"Missing required arguments: {', '.join(missing)}")

    print("\n===================================================")
    print(" QDEX - Post-DFT Exciton Solver")
    print("===================================================")

    tracker = ResourceTracker()

    tracker.start_stage("Geometry & Basis Parsing")
    print("\n--- Parsing Geometry and Basis Set ---")
    t0_parse = time.time()
    syms, coords_ang = read_xyz(args.xyz)
    basis_dict = parse_basis(args.basis_txt, args.basis_name, required_elements=set(syms))
    shells = build_shell_dicts(syms, coords_ang, basis_dict)
    shells = [{**sh, 'pure': True} for sh in shells] # Use Sphericals 
    n_ao = count_ao_from_shells(shells)
    atom_ao_ranges = build_atom_ao_ranges(shells)
    print(f"  -> Parsed in {time.time() - t0_parse:.2f} s | Total AOs: {n_ao}")

    tracker.start_stage("AO Overlap Matrix (S)")
    print("\n--- Computing AO overlap ---")
    t0_s = time.time()
    if getattr(args, "periodic_enabled", False):
        lattice_ang = np.asarray(args.lattice_vectors, dtype=float)
        lattice_bohr = lattice_ang * BOHR_PER_ANG
        cutoff = float(getattr(args, "overlap_cutoff", -1.0))
        if cutoff <= 0.0:
            print("  [Warning] periodic.overlap_cutoff <= 0 uses minimum-image overlap only.")
            print("  [Warning] Gamma-periodic MO orthonormality usually requires summing neighboring cell images; try overlap_cutoff: 6.0.")
        S = libint_cpp.overlap_pbc(shells, lattice_bohr, cutoff, args.nthreads)
        print(f"  ->  PBC overlap computed in {time.time() - t0_s:.2f} s")
        print(f"  ->  Lattice vectors read from YAML in angstrom; Libint lattice passed in bohr. cutoff={cutoff:.3f} A")
    else:
        S = libint_cpp.overlap(shells, args.nthreads)
        print(f"  ->  Finite-system overlap computed in {time.time() - t0_s:.2f} s")

    C_beta, eps_beta, occ_beta = None, None, None
    homo_index_beta = None

    tracker.start_stage("Molecular Orbitals (MOs)")
    if args.mo_file_beta is not None:
        print(f"\n--- Reading Alpha Molecular Orbitals from {args.mo_file} ---")
        t0_mos = time.time()
        C, eps, occ = read_mos_auto(args.mo_file, n_ao, verbose=True, cache=args.cache_mos)
        print(f"  -> Alpha MOs parsed in {time.time() - t0_mos:.2f} s | C shape {C.shape}")
        
        print(f"\n--- Reading Beta Molecular Orbitals from {args.mo_file_beta} ---")
        t0_beta = time.time()
        C_beta, eps_beta, occ_beta = read_mos_auto(args.mo_file_beta, n_ao, verbose=True, cache=args.cache_mos)
        print(f"  -> Beta MOs parsed in {time.time() - t0_beta:.2f} s | C shape {C_beta.shape}")
    else:
        print(f"\n--- Reading Molecular Orbitals from {args.mo_file} ---")
        t0_mos = time.time()
        C, eps, occ = read_mos_auto(args.mo_file, n_ao, verbose=True, cache=args.cache_mos)
        print(f"  -> MOs parsed in {time.time() - t0_mos:.2f} s | C shape {C.shape}")

    t0_gap = time.time()
    homo_index = np.where(occ > 0.0)[0].max()
    eps = eps * HA_TO_EV
    
    if eps_beta is not None:
        homo_index_beta = np.where(occ_beta > 0.0)[0].max()
        eps_beta = eps_beta * HA_TO_EV
    else:
        homo_index_beta = homo_index
        eps_beta = eps

    args.n_alpha_ref = float(np.sum(occ))
    args.n_beta_ref = float(np.sum(occ_beta)) if occ_beta is not None else float(np.sum(occ))
    args.S_ref = infer_reference_spin(args.n_alpha_ref, args.n_beta_ref)
    if C_beta is not None:
        ref_spin_name = spin_multiplicity_name(args.S_ref)
        print(
            f"  [Spin] UKS reference: N_alpha={args.n_alpha_ref:.1f}, "
            f"N_beta={args.n_beta_ref:.1f}, S_ref~{args.S_ref:.1f}, "
            f"multiplicity~{int(round(2 * args.S_ref + 1))} ({ref_spin_name})"
        )

    # Define Shifted Energy Axis
    e_fermi_raw = (eps[homo_index] + eps[homo_index + 1]) / 2.0
    eps_shifted = eps - e_fermi_raw
    eps_beta_shifted = eps_beta - e_fermi_raw
    e_homo = eps_shifted[homo_index]
    e_lumo = eps_beta_shifted[homo_index_beta + 1]

    tracker.start_stage("Quasiparticle (GW) Model")
    dft_gap = eps_beta[homo_index_beta + 1] - eps[homo_index]
    target_qp_gap = dft_gap
    confinement_energy = 0.0
    qp_provenance = None
    
    if isinstance(args.qp_gap, str):
        if args.qp_gap.lower() == "brus":
            target_qp_gap = estimate_brus_qp_gap(np.array(coords_ang), syms, args.material)
            if target_qp_gap is not None:
                scissor = target_qp_gap - dft_gap
                confinement_energy = target_qp_gap - MATERIAL_DB.get(args.material.upper(), [0]*9)[3]
            else:
                raise ValueError(
                    "Brus QP model requested, but the required material data are missing. "
                    "Select qp_gap: pbe explicitly for an uncorrected calculation."
                )
                
        elif args.qp_gap.lower() == "gw":
            if getattr(args, "periodic_enabled", False):
                print("\n  [Bulk GW Model] Periodic mode enabled: using tabulated bulk GW-PBE scissor.")
                gw_scissor, qp_provenance = estimate_periodic_bulk_gw_scissor(args.material)
            else:
                # Compute the Scaled GW Scissor directly
                gw_scissor, qp_provenance = estimate_gw_qp_gap(
                    np.array(coords_ang), syms, args.material, args.eps_out,
                    regularization_length_ang=args.qp_regularization_length,
                    residual_power=args.qp_residual_power,
                    strict=args.qp_strict,
                    return_details=True,
                )
            
            if gw_scissor is not None:
                scissor = gw_scissor
                target_qp_gap = dft_gap + scissor
            else:
                raise ValueError(
                    "GW QP model requested, but the required anchor/bulk data are missing. "
                    "Select qp_gap: pbe explicitly for an uncorrected calculation."
                )
                
            if args.estimate_qp:
                print("  [QP Warning] qp_gap is set to 'gw', which already uses the recommended scaled-GW hardness model.")
                print("  [QP Warning] estimate_qp enables the experimental COHSEX/TB Mulliken correction and can double-count QP shifts.")
                print("  [QP Warning] Production runs should use estimate_qp: false unless you explicitly want this experimental path.")
        elif args.qp_gap.lower() == "pbe":
            scissor = 0.0
            target_qp_gap = dft_gap
            print("  [QP] Explicit uncorrected PBE mode selected.")
        else:
            raise ValueError(
                f"Unknown qp_gap mode '{args.qp_gap}'. Use 'gw', 'brus', 'pbe', or a numeric gap."
            )
    else:
        # Numeric explicit gap provided
        target_qp_gap = float(args.qp_gap)
        scissor = target_qp_gap - dft_gap
    
    print(f"\n  [DFT] Initial Gap  : {dft_gap:.4f} eV")
    print(f"  [QP]  Target Gap   : {target_qp_gap:.4f} eV")
    print(f"  [QP]  Scissor Shift: {scissor:.4f} eV")
    qp_provenance_file = write_qp_provenance(qp_provenance, dft_gap, target_qp_gap, scissor, args)
    print_qp_provenance(qp_provenance, dft_gap=dft_gap, target_qp_gap=target_qp_gap, output_file=qp_provenance_file)
   
    # =========================================================================
    # DYNAMIC IP & EA PREDICTION (Vacuum-Anchored Projection Method)
    # =========================================================================
    dft_homo_raw = eps[homo_index]
    dft_lumo_raw = eps[homo_index + 1]
    
    entry = MATERIAL_DB.get(args.material.upper(), None)

    if getattr(args, "periodic_enabled", False):
        f_homo, f_lumo = 0.0, 1.0
        qp_homo = dft_homo_raw
        qp_lumo = dft_lumo_raw + scissor

        print(f"\n  [Periodic Band Edges]")
        print(f"    Raw CP2K HOMO    : {dft_homo_raw:8.4f} eV (arbitrary periodic eigenvalue zero)")
        print(f"    Raw CP2K LUMO    : {dft_lumo_raw:8.4f} eV")
        print(f"    Bulk GW scissor  : virtual manifold shifted by {scissor:+.4f} eV")
        print("    Note             : absolute IP/EA levels are not assigned in periodic mode.")

    elif entry is not None and len(entry) >= 14:
        pbe_h_mono, pbe_l_mono, gw_h_mono, gw_l_mono = entry[10], entry[11], entry[12], entry[13]
        gap_pbe_mono = pbe_l_mono - pbe_h_mono
        
        # 1. Asymmetry Fractions from the GW Anchor
        delta_h = gw_h_mono - pbe_h_mono
        delta_l = gw_l_mono - pbe_l_mono
        anchor_gap_opening = delta_l - delta_h

        if anchor_gap_opening > 1.0e-12 and delta_h <= 0.0 and delta_l >= 0.0:
            f_homo = -delta_h / anchor_gap_opening
            f_lumo = delta_l / anchor_gap_opening
        else:
            f_homo = f_lumo = 0.5
            print("  [QP Warning] Anchor frontier shifts do not bracket the PBE gap; using a symmetric edge split.")
        
        # 2. Project Absolute PBE Levels (Bypassing CP2K floating vacuum)
        # We use the computed intermediate PBE gap (dft_gap) as the physical truth
        shrinkage_pbe = gap_pbe_mono - dft_gap 
        
        true_pbe_homo = pbe_h_mono + (shrinkage_pbe * f_homo)
        true_pbe_lumo = pbe_l_mono - (shrinkage_pbe * f_lumo)
        
        # 3. Apply the Dielectric Scissor to get QP levels
        qp_homo = true_pbe_homo - (scissor * f_homo)
        qp_lumo = true_pbe_lumo + (scissor * f_lumo)
        
        print(f"\n  [Absolute Band Edges (IP & EA)]")
        print(f"    Raw CP2K HOMO    : {dft_homo_raw:8.4f} eV (Floating Vacuum)")
        print(f"    Modeled PBE HOMO : {true_pbe_homo:8.4f} eV (anchor-reconstructed)")
        print(f"    -> Shift Split   : HOMO takes {f_homo*100:.1f}%, LUMO takes {f_lumo*100:.1f}%")

    else:
        # Fallback if no 14-item monomer data is available
        f_homo, f_lumo = 0.5, 0.5
        qp_homo = dft_homo_raw - (scissor * f_homo)
        qp_lumo = dft_lumo_raw + (scissor * f_lumo)
        
        print(f"\n  [Absolute Band Edges (IP & EA)]")
        print(f"    Raw CP2K HOMO    : {dft_homo_raw:8.4f} eV")
        print(f"    -> Shift Split   : HOMO takes 50.0%, LUMO takes 50.0% (Default)")

    if getattr(args, "periodic_enabled", False):
        print(f"    QP HOMO-like     : {qp_homo:8.4f} eV (relative eigenvalue)")
        print(f"    QP LUMO-like     : {qp_lumo:8.4f} eV (relative eigenvalue)")
    else:
        print(f"    QP HOMO (IP)     : {qp_homo:8.4f} eV   -> IP = {-qp_homo:8.4f} eV")
        print(f"    QP LUMO (EA)     : {qp_lumo:8.4f} eV   -> EA = {-qp_lumo:8.4f} eV")
    if qp_provenance is not None:
        qp_provenance.update({
            "target_pbe_gap_ev": float(dft_gap),
            "target_qp_gap_ev": float(target_qp_gap),
            "modeled_qp_homo_ev": float(qp_homo),
            "modeled_qp_lumo_ev": float(qp_lumo),
            "modeled_ip_ev": float(-qp_homo),
            "modeled_ea_ev": float(-qp_lumo),
            "absolute_edge_model": "anchor_reconstructed",
        })
        write_qp_provenance(qp_provenance, dft_gap, target_qp_gap, scissor, args)
    # =========================================================================

    print(f"  -> Energy axis shifted and target gap resolved in {time.time() - t0_gap:.4f} s")
 
    # -----------------------------------------------------------------
    # Unified S@C Computation
    # -----------------------------------------------------------------
    tracker.start_stage("MO Orthonormality & Populations")
    print("\n--- Computing Unified S@C Population Analysis ---")
    t0_pop = time.time()
    C_dense = C.toarray() if hasattr(C, 'toarray') else C
    SC_dense = S @ C_dense 
    C_dense_beta_pop, SC_dense_beta_pop, pops_beta = None, None, None
    
    # === DIAGNOSTIC: STRICT C^T S C ORTHONORMALITY CHECK ===
    overlap_label = "PBC" if getattr(args, "periodic_enabled", False) else "finite"
    print(f"\n  [Diag] Testing MO Orthonormality with {overlap_label} overlap (C^T S C = I) ...")
    norm_matrix = C_dense.conj().T @ SC_dense
    orth_delta = norm_matrix - np.eye(C_dense.shape[1])
    orth_err = np.linalg.norm(orth_delta)
    orth_max = np.max(np.abs(orth_delta))
    print(f"[CHECK] alpha ||C†SC - I||_F = {orth_err:.3e}")
    print(f"[CHECK] alpha max|C†SC - I|  = {orth_max:.3e}")
    if orth_max > args.orthonormality_tol:
        raise ValueError(
            f"MO orthonormality failure: max|C†SC-I|={orth_max:.3e} exceeds "
            f"{args.orthonormality_tol:.3e}. Check AO ordering, normalization, and spherical conventions."
        )
    soc_assume_orthonormal = orth_max < 1.0e-6

    if C_beta is not None:
        C_dense_beta_pop = C_beta.toarray() if hasattr(C_beta, 'toarray') else np.asarray(C_beta)
        SC_dense_beta_pop = S @ C_dense_beta_pop
        norm_matrix_beta = C_dense_beta_pop.conj().T @ SC_dense_beta_pop
        orth_delta_beta = norm_matrix_beta - np.eye(C_dense_beta_pop.shape[1])
        orth_err_beta = np.linalg.norm(orth_delta_beta)
        orth_max_beta = np.max(np.abs(orth_delta_beta))
        print(f"[CHECK] beta  ||C†SC - I||_F = {orth_err_beta:.3e}")
        print(f"[CHECK] beta  max|C†SC - I|  = {orth_max_beta:.3e}")
        if orth_max_beta > args.orthonormality_tol:
            raise ValueError(
                f"Beta MO orthonormality failure: max|C†SC-I|={orth_max_beta:.3e} exceeds "
                f"{args.orthonormality_tol:.3e}."
            )
        soc_assume_orthonormal = soc_assume_orthonormal and orth_max_beta < 1.0e-6
        cross_err = np.linalg.norm(C_dense.conj().T @ SC_dense_beta_pop)
        print(f"[CHECK] alpha/beta ||Caᵀ S Cb|| = {cross_err:.3e} (diagnostic)")
    # ===========================================================

    pops_sf = np.real(C_dense.conj() * SC_dense)
    if C_dense_beta_pop is not None:
        pops_beta = np.real(C_dense_beta_pop.conj() * SC_dense_beta_pop)
    print(f"  -> S@C projection and populations computed in {time.time() - t0_pop:.2f} s")

    # -----------------------------------------------------------------
    # BSE Active Space Setup
    # -----------------------------------------------------------------
    is_uks_sp = (C_beta is not None) and (not args.triplet)
    if args.e_thresh is not None:
        raw_gap = eps_beta_shifted[homo_index_beta + 1:].reshape(1, -1) - eps_shifted[:homo_index + 1].reshape(-1, 1)
        valid_pairs = raw_gap <= args.e_thresh
        bse_n_occ = homo_index - np.where(np.any(valid_pairs, axis=1))[0][0] + 1
        bse_n_virt = np.where(np.any(valid_pairs, axis=0))[0][-1] + 1
    else:
        bse_n_occ = args.nhomos if args.nhomos is not None else 25
        bse_n_virt = args.nlumos if args.nlumos is not None else 25

    # Safe bound active space sizes to avoid out of bounds
    bse_n_occ = min(bse_n_occ, homo_index + 1)
    bse_n_virt = min(bse_n_virt, len(eps) - (homo_index + 1), len(eps_beta) - (homo_index_beta + 1))
    bse_n_occ_beta = min(bse_n_occ, homo_index_beta + 1)
    bse_n_virt_beta = min(bse_n_virt, len(eps_beta) - (homo_index_beta + 1))

    bse_active_indices = np.arange(homo_index - bse_n_occ + 1, homo_index + 1 + bse_n_virt)
    bse_active_indices_beta = np.arange(homo_index_beta - bse_n_occ_beta + 1, homo_index_beta + 1 + bse_n_virt_beta)
    bse_soc_E, bse_soc_U, bse_spinor_homo_idx = None, None, None
    calculated_soc_gap = None
    soc_overlap_cache = None
    
    if args.soc_flag:
        tracker.start_stage("SOC Active Space (BSE)")
        print(f"\n--- Computing SOC Spinor Subspace for BSE (Small Window) ---")
        if is_uks_sp:
            from qdex.soc_utils import compute_spinor_subspace_uks
            C_dense_beta = C_beta.toarray() if hasattr(C_beta, 'toarray') else np.asarray(C_beta)
            bse_soc_E, bse_soc_U, soc_overlap_cache = compute_spinor_subspace_uks(
                atom_symbols=syms, coords_ang=coords_ang, shells=shells,
                C_alpha_AO=C_dense, eps_alpha_Ha=eps / HA_TO_EV, active_alpha_indices=bse_active_indices,
                C_beta_AO=C_dense_beta, eps_beta_Ha=eps_beta / HA_TO_EV, active_beta_indices=bse_active_indices_beta,
                S_AO=S, gth_file=args.gth_file, nthreads=args.nthreads,
                assume_orthonormal=soc_assume_orthonormal,
                SC_alpha_AO=SC_dense, SC_beta_AO=SC_dense_beta_pop,
                device=args.device,
            )
            bse_spinor_homo_idx = bse_n_occ + bse_n_occ_beta - 1
        else:
            from qdex.soc_utils import compute_spinor_subspace
            bse_soc_E, bse_soc_U, soc_overlap_cache = compute_spinor_subspace(
                atom_symbols=syms, coords_ang=coords_ang, shells=shells, 
                C_AO=C_dense, eps_Ha=eps / HA_TO_EV, S_AO=S, 
                active_indices=bse_active_indices, gth_file=args.gth_file,
                nthreads=args.nthreads, assume_orthonormal=soc_assume_orthonormal,
                SC_AO=SC_dense, device=args.device,
            )
            bse_spinor_homo_idx = (bse_n_occ * 2) - 1
        bse_soc_E = (bse_soc_E * HA_TO_EV) - e_fermi_raw
        bse_soc_E -= (bse_soc_E[bse_spinor_homo_idx] + bse_soc_E[bse_spinor_homo_idx + 1]) / 2.0

    if C_beta is not None:
        pop_range = max(0, int(getattr(args, "population_print_range", 15)))
        pop_tags = getattr(args, "population_bars", None)
        print("\n--- Spin-Free Alpha MO Population Analysis ---")
        print_orbital_summary(eps_shifted, occ, homo_index, pops_sf, syms, shells, is_soc=False, print_range=pop_range, population_bars=pop_tags)
        print("\n--- Spin-Free Beta MO Population Analysis ---")
        print_orbital_summary(eps_beta_shifted, occ_beta, homo_index_beta, pops_beta, syms, shells, is_soc=False, print_range=pop_range, population_bars=pop_tags)
    else:
        print("\n--- Spin-Free MO Population Analysis ---")
        pop_range = max(0, int(getattr(args, "population_print_range", 15)))
        pop_tags = getattr(args, "population_bars", None)
        print_orbital_summary(eps_shifted, occ, homo_index, pops_sf, syms, shells, is_soc=False, print_range=pop_range, population_bars=pop_tags)

    if args.soc_flag:
        if args.gth_file is None: sys.exit("ERROR: --gth_file is required when --soc_flag is enabled.")
        print("\n--- SOC Spinor Population Analysis (BSE Active Space) ---")
        t_pop = time.time()
        
        if is_uks_sp:
            C_dense_beta = C_beta.toarray() if hasattr(C_beta, 'toarray') else np.asarray(C_beta)
            SC_dense_beta = S @ C_dense_beta
            C_act = C_dense[:, bse_active_indices]
            C_act_beta = C_dense_beta[:, bse_active_indices_beta]
            SC_act = SC_dense[:, bse_active_indices]
            SC_act_beta = SC_dense_beta[:, bse_active_indices_beta]
            n_mo_act = len(bse_active_indices)
            n_mo_act_beta = len(bse_active_indices_beta)
            C_spinor_act_alpha = C_act @ bse_soc_U[:n_mo_act, :]
            C_spinor_act_beta  = C_act_beta @ bse_soc_U[n_mo_act:n_mo_act + n_mo_act_beta, :]
            SC_spinor_act_alpha = SC_act @ bse_soc_U[:n_mo_act, :]
            SC_spinor_act_beta  = SC_act_beta @ bse_soc_U[n_mo_act:n_mo_act + n_mo_act_beta, :]
        else:
            C_act = C_dense[:, bse_active_indices]
            SC_act = SC_dense[:, bse_active_indices]
            n_mo_act = len(bse_active_indices)
            C_spinor_act_alpha = C_act @ bse_soc_U[:n_mo_act, :]
            C_spinor_act_beta  = C_act @ bse_soc_U[n_mo_act:, :]
            SC_spinor_act_alpha = SC_act @ bse_soc_U[:n_mo_act, :]
            SC_spinor_act_beta  = SC_act @ bse_soc_U[n_mo_act:, :]
        
        pops_soc_act = np.real(C_spinor_act_alpha.conj() * SC_spinor_act_alpha) + \
                       np.real(C_spinor_act_beta.conj() * SC_spinor_act_beta)
                       
        print(f"  -> Projected populations in {time.time() - t_pop:.2f}s")
        
        soc_occ = np.zeros_like(bse_soc_E)
        soc_occ[:bse_spinor_homo_idx + 1] = 1.0
        soc_offset = (homo_index - bse_n_occ + 1) * 2 
        print_orbital_summary(
            bse_soc_E, soc_occ, bse_spinor_homo_idx, pops_soc_act, syms, shells,
            is_soc=True, offset=soc_offset, print_range=pop_range, population_bars=pop_tags
        )
        calculated_soc_gap = bse_soc_E[bse_spinor_homo_idx + 1] - bse_soc_E[bse_spinor_homo_idx]

        # --- RIGID SCISSOR APPLICATION ---
        print("\n--- Scissor Operator Application ---")
        print(f"  Rigid Scissor (Computed above)  : {scissor:+8.4f} eV")
        print(f"  Spin-Free DFT Gap               : {dft_gap:8.4f} eV")
        print(f"  SOC-Shrunken DFT Gap            : {calculated_soc_gap:8.4f} eV")
        
        final_sf_qp_gap = dft_gap + scissor
        final_soc_qp_gap = calculated_soc_gap + scissor
        
        print(f"  -> Final Spin-Free QP Gap       : {final_sf_qp_gap:8.4f} eV")
        print(f"  -> Final SOC QP Gap             : {final_soc_qp_gap:8.4f} eV (D_SOC = {dft_gap - calculated_soc_gap:.4f} eV)")

    else:
        # --- RIGID SCISSOR FOR SPIN-FREE ONLY ---
        print("\n--- Scissor Operator Application (Spin-Free) ---")
        print(f"  Rigid Scissor (Computed above)  : {scissor:+8.4f} eV")
        print(f"  Spin-Free DFT Gap               : {dft_gap:8.4f} eV")
        print(f"  -> Final Spin-Free QP Gap       : {dft_gap + scissor:8.4f} eV")

    # Update Confinement Energy (Always relative to bulk)
    db_gap = MATERIAL_DB.get(args.material.upper(), MATERIAL_DB["DEFAULT"])[3]
    confinement_energy = target_qp_gap - db_gap

    # -----------------------------------------------------------------
    # EXCITON CUBE GENERATION: EXECUTED BEFORE FUZZY PLOTTING 
    # -----------------------------------------------------------------
    if getattr(args, 'cube', False):
        tracker.start_stage("Exciton Cube Generation")
        from qdex.exciton_cube import generate_cubes
        print("\n--- Generating Cubes for MOs / Spinors ---")
        
        class DummySolver:
            def __init__(self):
                self.C = C_dense
                self.homo_index = homo_index
                self.n_occ = bse_n_occ
                class DummyHam: pass
                self.ham = DummyHam()
                self.soc_flag = args.soc_flag
                if args.soc_flag:
                    self.ham.n_occ_spinor = bse_spinor_homo_idx + 1
                    
        dummy = DummySolver()
        mo_list = []
        spinor_list = []
        n_h = getattr(args, 'cube_nhomos', 2)
        n_l = getattr(args, 'cube_nlumos', 2)
       
        # ALWAYS generate Spin-Free MOs
        mo_homo = dummy.homo_index
        mo_list = [mo_homo - i for i in range(n_h) if (mo_homo - i) >= 0]
        mo_list += [mo_homo + 1 + i for i in range(n_l) if (mo_homo + 1 + i) < dummy.C.shape[1]]
        
        # ADDITIONALLY generate Spinors if SOC is enabled
        if args.soc_flag:
            sp_homo = dummy.ham.n_occ_spinor - 1
            spinor_list = [sp_homo - i for i in range(n_h) if (sp_homo - i) >= 0]
            spinor_list += [sp_homo + 1 + i for i in range(n_l) if (sp_homo + 1 + i) < bse_soc_U.shape[1]]
 
        generate_cubes(
            solver=dummy, 
            bse_states_dict={}, 
            mo_list=mo_list,
            spinor_list=spinor_list,
            soc_U=bse_soc_U if args.soc_flag else None,
            shells=shells, symbols=syms, coords=coords_ang, 
            spacing_ang=args.cube_spacing, nthreads=args.nthreads,
            use_cpp=not getattr(args, 'disable_cpp_cube', False)
        )

    # -----------------------------------------------------------------
    # MODULE DELEGATION: FUZZY BANDS & PDOS (Large Window)
    # -----------------------------------------------------------------
    if getattr(args, 'run_fuzzy', False):
        tracker.start_stage("Fuzzy Bands & PDOS/COOP")
        fuzzy_active_indices = _select_soc_window_indices(eps_shifted, homo_index, args.soc_window)
        fuzzy_active_indices_beta = _select_soc_window_indices(eps_beta_shifted, homo_index_beta, args.soc_window)

        qp_occ_shift_abs = qp_homo - eps[homo_index]
        qp_virt_shift_abs = qp_lumo - eps_beta[homo_index_beta + 1]
        qp_energies_abs = build_qp_energies_vacuum(eps, homo_index, occ_shift=qp_occ_shift_abs, virt_shift=qp_virt_shift_abs)
        qp_energies_beta_abs = build_qp_energies_vacuum(eps_beta, homo_index_beta, occ_shift=qp_occ_shift_abs, virt_shift=qp_virt_shift_abs)
        qp_energies_rel = build_qp_energies(eps_shifted, homo_index, scissor_ev=scissor)
        qp_energies_beta_rel = build_qp_energies(eps_beta_shifted, homo_index_beta, scissor_ev=scissor)
        
        fuzzy_soc_E, fuzzy_soc_E_abs, fuzzy_soc_U, fuzzy_spinor_homo_idx = None, None, None, None
        if args.soc_flag:
            print(f"\n--- Computing SOC Spinor Subspace for Fuzzy Bands (|E-Ef| <= {args.soc_window:.3f} eV) ---")
            print(f"  -> Alpha fuzzy SOC active MOs: {len(fuzzy_active_indices)} / {len(eps)}")
            if is_uks_sp:
                print(f"  -> Beta fuzzy SOC active MOs : {len(fuzzy_active_indices_beta)} / {len(eps_beta)}")
            if is_uks_sp:
                from qdex.soc_utils import compute_spinor_subspace_uks
                C_dense_beta = C_beta.toarray() if hasattr(C_beta, 'toarray') else np.asarray(C_beta)
                fuzzy_soc_E, fuzzy_soc_U, _ = compute_spinor_subspace_uks(
                    atom_symbols=syms, coords_ang=coords_ang, shells=shells,
                    C_alpha_AO=C_dense, eps_alpha_Ha=eps / HA_TO_EV, active_alpha_indices=fuzzy_active_indices,
                    C_beta_AO=C_dense_beta, eps_beta_Ha=eps_beta / HA_TO_EV, active_beta_indices=fuzzy_active_indices_beta,
                    S_AO=S, gth_file=args.gth_file, nthreads=args.nthreads,
                    soc_cache=soc_overlap_cache, assume_orthonormal=soc_assume_orthonormal,
                    SC_alpha_AO=SC_dense, SC_beta_AO=SC_dense_beta_pop,
                    device=args.device,
                )
                f_n_occ = np.sum(fuzzy_active_indices <= homo_index)
                f_n_occ_beta = np.sum(fuzzy_active_indices_beta <= homo_index_beta)
                fuzzy_spinor_homo_idx = f_n_occ + f_n_occ_beta - 1
            else:
                from qdex.soc_utils import compute_spinor_subspace
                fuzzy_soc_E, fuzzy_soc_U, _ = compute_spinor_subspace(
                    atom_symbols=syms, coords_ang=coords_ang, shells=shells,
                    C_AO=C_dense, eps_Ha=eps / HA_TO_EV, S_AO=S,
                    active_indices=fuzzy_active_indices, gth_file=args.gth_file, nthreads=args.nthreads,
                    soc_cache=soc_overlap_cache, assume_orthonormal=soc_assume_orthonormal,
                    SC_AO=SC_dense, device=args.device,
                )
                f_n_occ = np.sum(fuzzy_active_indices <= homo_index)
                fuzzy_spinor_homo_idx = (f_n_occ * 2) - 1
            fuzzy_soc_E_abs = fuzzy_soc_E * HA_TO_EV
            fuzzy_soc_E = fuzzy_soc_E_abs - e_fermi_raw
            fuzzy_soc_E -= (fuzzy_soc_E[fuzzy_spinor_homo_idx] + fuzzy_soc_E[fuzzy_spinor_homo_idx + 1]) / 2.0

        run_fuzzy_bands_and_pdos(
            args, C_dense, S, eps_shifted, occ, homo_index, e_homo, e_lumo, e_fermi_raw, 
            syms, coords_ang, shells, pops_sf, 
            soc_active_indices=fuzzy_active_indices, soc_E_act=fuzzy_soc_E, soc_U_act=fuzzy_soc_U, spinor_homo_idx=fuzzy_spinor_homo_idx,
            qp_energies=qp_energies_rel,
            eps_abs=eps, qp_energies_abs=qp_energies_abs, soc_E_abs_act=fuzzy_soc_E_abs,
            C_beta_dense=(C_beta.toarray() if hasattr(C_beta, 'toarray') else np.asarray(C_beta)) if is_uks_sp else None,
            eps_beta_shifted=eps_beta_shifted if is_uks_sp else None,
            eps_beta_abs=eps_beta if is_uks_sp else None,
            homo_index_beta=homo_index_beta if is_uks_sp else None,
            qp_energies_beta=qp_energies_beta_rel if is_uks_sp else None,
            qp_energies_beta_abs=qp_energies_beta_abs if is_uks_sp else None,
            soc_active_indices_beta=fuzzy_active_indices_beta if is_uks_sp else None
        )


    # -----------------------------------------------------------------
    # EARLY EXIT LOGIC (If run_bse is False)
    # -----------------------------------------------------------------
    run_bse = getattr(args, 'run_bse', True)
    if not run_bse:
        tracker.end_stage()
        tracker.print_summary(device=args.device, nthreads=args.nthreads)
        print("\n--- BSE Calculation Skipped (run_bse: false) ---")
        print("\nAll requested tasks finished successfully.")
        return

    # -----------------------------------------------------------------
    # BSE CONTINUATION (Only if run_bse is True)
    # -----------------------------------------------------------------
    tracker.start_stage("Transition Dipoles")
    print("\n--- Computing Transition Dipoles ---")
    if getattr(args, "periodic_enabled", False):
        print("  [Warning] Periodic mode is enabled, but transition dipoles use the finite-cell AO position operator.")
        print("  [Warning] Oscillator strengths should be treated as Gamma-only finite-cell approximations.")
    t0_dip = time.time()
    mu_ao_x, mu_ao_y, mu_ao_z = compute_dipole_ao(shells, nthreads=args.nthreads)
    print(f"  ->  Dipoles computed in {time.time() - t0_dip:.2f} s")

    compute_device, dev_obj = resolve_device(args.device, verbose=True)
    if dev_obj is not None:
        try:
            import torch
            torch.set_num_threads(args.nthreads)
        except Exception:
            pass

    # ---------------------------------------------------------------
    # Determine spin mode and compute transition dipoles accordingly
    # ---------------------------------------------------------------
    is_uks_sp = (C_beta is not None) and (not args.triplet)
    spin_mode = 'uks_spin_preserving' if is_uks_sp else ('triplet' if args.triplet else 'singlet')

    if is_uks_sp:
        # Manifold B: compute alpha and beta dipoles separately and pass as tuples
        print("  [UKS Spin-Preserving] Computing alpha-alpha and beta-beta transition dipoles separately...")
        C_dense_alpha = C.toarray() if hasattr(C, 'toarray') else np.asarray(C)
        C_dense_beta  = C_beta.toarray() if hasattr(C_beta, 'toarray') else np.asarray(C_beta)

        C_occ_a  = C_dense_alpha[:, homo_index - bse_n_occ + 1 : homo_index + 1].astype(np.float64)
        C_virt_a = C_dense_alpha[:, homo_index + 1 : homo_index + 1 + bse_n_virt].astype(np.float64)
        C_occ_b  = C_dense_beta[:, homo_index_beta - bse_n_occ + 1 : homo_index_beta + 1].astype(np.float64)
        C_virt_b = C_dense_beta[:, homo_index_beta + 1 : homo_index_beta + 1 + bse_n_virt].astype(np.float64)

        def transform_dipole_chan(mu_ao, C_occ_ch, C_virt_ch):
            return transform_ao_operator(mu_ao, C_occ_ch, C_virt_ch, compute_device)

        mu_ia_x_a = transform_dipole_chan(mu_ao_x, C_occ_a, C_virt_a)
        mu_ia_y_a = transform_dipole_chan(mu_ao_y, C_occ_a, C_virt_a)
        mu_ia_z_a = transform_dipole_chan(mu_ao_z, C_occ_a, C_virt_a)
        mu_ia_x_b = transform_dipole_chan(mu_ao_x, C_occ_b, C_virt_b)
        mu_ia_y_b = transform_dipole_chan(mu_ao_y, C_occ_b, C_virt_b)
        mu_ia_z_b = transform_dipole_chan(mu_ao_z, C_occ_b, C_virt_b)

        # Pass as tuples: (alpha_block, beta_block)
        mu_ia_x = (mu_ia_x_a, mu_ia_x_b)
        mu_ia_y = (mu_ia_y_a, mu_ia_y_b)
        mu_ia_z = (mu_ia_z_a, mu_ia_z_b)

        # Sizes for the beta active space (may differ if bse_n_occ/bse_n_virt are clamped differently)
        bse_n_occ_beta = min(bse_n_occ, homo_index_beta + 1)
        bse_n_virt_beta = min(bse_n_virt, len(eps_beta) - (homo_index_beta + 1))
    else:
        # Singlet or spin-flip triplet: standard single-channel dipoles
        C_occ = C[:, homo_index - bse_n_occ + 1 : homo_index + 1].astype(np.float64)
        if C_beta is not None:
            C_virt = C_beta[:, homo_index_beta + 1 : homo_index_beta + 1 + bse_n_virt].astype(np.float64)
        else:
            C_virt = C[:, homo_index + 1 : homo_index + 1 + bse_n_virt].astype(np.float64)

        def transform_dipole(mu_ao):
            return transform_ao_operator(mu_ao, C_occ, C_virt, compute_device)

        mu_ia_x, mu_ia_y, mu_ia_z = transform_dipole(mu_ao_x), transform_dipole(mu_ao_y), transform_dipole(mu_ao_z)
        bse_n_occ_beta  = bse_n_occ
        bse_n_virt_beta = bse_n_virt

    del mu_ao_x, mu_ao_y, mu_ao_z; gc.collect()

    # Store this for the analysis printouts later
    args.qp_gap_num = target_qp_gap

    tracker.start_stage("BSE Exciton Solver (Spin-Free)")
    solver_sf = ExcitonSolver(
        C=C, eps=eps_shifted, occ=occ, overlap=S, atom_symbols=syms, atom_coords=np.array(coords_ang),
        atom_ao_ranges=atom_ao_ranges, homo_index=homo_index, n_occ=bse_n_occ, n_virt=bse_n_virt, 
        scissor_ev=scissor, kernel=args.kernel, alpha=args.alpha, beta=args.beta, include_exchange=args.include_direct_eh,
        include_direct_eh=args.include_direct_eh,
        estimate_qp=args.estimate_qp, material=args.material, e_thresh=args.e_thresh, f_thresh=args.f_thresh,
        mu_ia_x=mu_ia_x, mu_ia_y=mu_ia_y, mu_ia_z=mu_ia_z, eps_out=args.eps_out,
        soc_U=None, soc_E=None, device=compute_device,
        vxc_ao_path=args.vxc_ao,
        nthreads=args.nthreads,
        spin=spin_mode,
        C_beta=C_beta, eps_beta=eps_beta_shifted, homo_index_beta=homo_index_beta,
        charge_type=args.charge_type,
        n_occ_beta=bse_n_occ_beta, n_virt_beta=bse_n_virt_beta,
        excitation_mode=args.excitation_mode
    )

    if args.estimate_qp:
        # Pull the dynamically calculated QP gap correction from the solver
        calculated_gap_shift = solver_sf.ham.sigma_virt[0] - solver_sf.ham.sigma_occ[-1]
        
        if getattr(args, 'use_cohsex_gap', False):
            print(f"\n  [QP] OVERRIDE: Using Pure COHSEX Gap Correction ({calculated_gap_shift:.4f} eV) instead of Tabulated GW.")
            scissor = calculated_gap_shift
        else:
            print(f"\n  [QP] Note: Tabulated GW Scissor ({scissor:.4f} eV) was used as the anchor. COHSEX provided orbital dispersion.")
            
        # Update confinement energy based on the final total gap
        confinement_energy = (dft_gap + scissor) - db_gap

    suffix_sf = "_sf" if args.soc_flag else ""
    run_solver_and_analysis(solver_sf, np.array(coords_ang), syms, shells, mu_ia_x, mu_ia_y, mu_ia_z, 
                            dft_gap, scissor, confinement_energy, args, suffix=suffix_sf)

    # Extract the computed shifts to avoid redundant calculation in SOC
    precalc_sigma = None
    if args.estimate_qp and hasattr(solver_sf.ham, 'sigma_occ'):
        precalc_sigma = (solver_sf.ham.sigma_occ, solver_sf.ham.sigma_virt)

    if args.soc_flag:
        tracker.start_stage("BSE Exciton Solver (SOC)")
        solver_soc = ExcitonSolver(
            C=C, eps=eps_shifted, occ=occ, overlap=S, atom_symbols=syms, atom_coords=np.array(coords_ang),
            atom_ao_ranges=atom_ao_ranges, homo_index=homo_index, n_occ=bse_n_occ, n_virt=bse_n_virt, 
            scissor_ev=scissor, kernel=args.kernel, alpha=args.alpha, beta=args.beta, include_exchange=args.include_direct_eh,
            include_direct_eh=args.include_direct_eh,
            estimate_qp=args.estimate_qp, material=args.material, e_thresh=args.e_thresh, f_thresh=args.f_thresh, 
            mu_ia_x=mu_ia_x, mu_ia_y=mu_ia_y, mu_ia_z=mu_ia_z, eps_out=args.eps_out,
            soc_U=bse_soc_U, soc_E=bse_soc_E, device=compute_device, 
            vxc_ao_path=args.vxc_ao,
            nthreads=args.nthreads,
            spin=spin_mode,
            C_beta=C_beta, eps_beta=eps_beta_shifted, homo_index_beta=homo_index_beta,
            charge_type=args.charge_type,
            n_occ_beta=bse_n_occ_beta, n_virt_beta=bse_n_virt_beta,
            excitation_mode=args.excitation_mode
        )
 
        run_solver_and_analysis(solver_soc, np.array(coords_ang), syms, shells, mu_ia_x, mu_ia_y, mu_ia_z, 
                                dft_gap, scissor, confinement_energy, args, suffix="_soc", soc_gap=calculated_soc_gap, soc_U=bse_soc_U, soc_E=bse_soc_E)

    tracker.end_stage()
    tracker.print_summary(device=args.device, nthreads=args.nthreads)
    print("\nAll calculations finished successfully.")

if __name__ == "__main__":
    main()
