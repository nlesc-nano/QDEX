import argparse
import atexit
import json
import numpy as np
import sys
import os
import time 
import gc
import inspect
import platform
import yaml 
from qdex.config_schema import SECTIONS, section_dest

import libint_cpp

if __package__ is None or __package__ == "":
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from qdex.io_utils import (
    read_xyz, parse_basis, build_shell_dicts,
    count_ao_from_shells, build_atom_ao_ranges, read_mos_auto, read_mos_uks,
    read_geometry_h5, is_h5_file, geometry_source
)
from qdex.solver import ExcitonSolver
from qdex.constants import HA_TO_EV, BOHR_PER_ANG
from qdex.exciton_analysis import ExcitonAnalyzer, plot_analysis_summary
from qdex.integrals import compute_dipole_ao
from qdex.oscillator import compute_oscillator_strengths
from qdex.hardness import MATERIAL_DB, estimate_brus_qp_gap, estimate_gw_qp_gap, build_gamma, anchor_bulk_homo_fraction
from qdex.qp_levels import (xs_shared_w, atom_delta_w, orbital_qp_energies, orbital_populations,
                             cohsex_diagonal, cohsex_qp_energies, anchor_key, load_anchor_table,
                             save_anchor_entry)
from qdex.hardness import get_cluster_size_metrics
from qdex.orbital_analysis import (
    compute_spin_character, compute_uks_soc_spin_free_channels,
    compute_uks_spin_free_channels, format_uks_soc_spin_free_character,
    format_uks_spin_free_character, infer_reference_spin,
    print_orbital_summary, spin_multiplicity_name
)
from qdex.fuzzy_bands import run_fuzzy_bands_and_pdos, build_qp_energies, build_qp_energies_vacuum
from qdex.nto import run_nto_analysis
from qdex.profiler import ResourceTracker
import logging

logger = logging.getLogger("qdex.cli")


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


def setup_run_logging(log_file, verbosity="full"):
    """Console output at ``verbosity``; with a log file, everything is also written to it."""
    from qdex.logging_setup import attach_run_log, detach_run_log, set_verbosity
    set_verbosity(verbosity)
    if log_file in (None, "", "none", "None", False):
        return None

    log_handle = open(log_file, "w", encoding="utf-8")
    # Logged messages: console at the chosen verbosity, the log file in full.
    # Anything written directly to stdout/stderr (C++ extensions, libraries) is copied to both.
    attach_run_log(log_handle, sys.__stdout__, verbosity)
    sys.stdout = TeeStream(sys.__stdout__, log_handle)
    sys.stderr = TeeStream(sys.__stderr__, log_handle)

    def close_log():
        detach_run_log()
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        log_handle.close()

    atexit.register(close_log)
    logger.info(f"  [Log] Writing full run log to {log_file}")
    return log_handle


def print_qp_provenance(details, dft_gap=None, target_qp_gap=None, output_file=None):
    if not details:
        return

    if details.get("qp_model") in ["sgw_dw", "sgw_dim", "sgw_resta", "evgw_dim", "evgw_resta", "qsgw_dim", "qsgw_resta"]:
        if details.get("qp_model") in ["qsgw_dim", "evgw_dim", "sgw_dim"]:
            model_name = "qsGW-DIM" if details.get("orbital_update") else ("evGW-DIM" if details.get("self_consistent") else "sGW-DIM")
        elif details.get("qp_model") in ["qsgw_resta", "evgw_resta", "sgw_resta"]:
            model_name = "qsGW-Resta" if details.get("orbital_update") else ("evGW-Resta" if details.get("self_consistent") else "sGW-Resta")
        else:
            model_name = "Microscopic sGW / Delta-W"
        logger.info(f"\n  [QP Provenance] {model_name} Model")
        logger.info(f"    Material                 : {details['material']}")
        logger.info(f"    Bulk PBE -> GW gap       : {details['bulk_pbe_gap_ev']:.3f} -> {details['bulk_gw_gap_ev']:.3f} eV")
        logger.info(f"    Bulk GW shift            : {details['bulk_shift_ev']:+.3f} eV")
        if details.get("dynamic_z"):
            z_info = f"dynamic: Zh={details.get('z_homo', details['z_factor']):.3f}, Zl={details.get('z_lumo', details['z_factor']):.3f}"
        else:
            z_info = f"Z = {details['z_factor']:.2f}"
        logger.info(f"    Confinement shift        : {details['confinement_shift_ev']:+.3f} eV ({z_info})")
        if "confinement_shift_internal_ev" in details and "confinement_shift_solvent_ev" in details:
            logger.info(f"      - Internal Contrast    : {details['confinement_shift_internal_ev']:+.3f} eV")
            logger.info(f"      - Solvent Reaction Field: {details['confinement_shift_solvent_ev']:+.3f} eV")
        logger.info(f"    HOMO shift (Delta Sigma) : {details['delta_sigma_homo_ev']:+.3f} eV")
        logger.info(f"    LUMO shift (Delta Sigma) : {details['delta_sigma_lumo_ev']:+.3f} eV")
        if details.get("orbital_update"):
            logger.info(f"    qsGW Iterations          : {details.get('qsgw_iterations', 1)} (converged: {details.get('qsgw_converged', False)})")
            logger.info(f"    HOMO Orbital Fidelity    : {details.get('homo_fidelity', 1.0):.5f}  (|<psi_PBE|psi_QP>|^2)")
            logger.info(f"    LUMO Orbital Fidelity    : {details.get('lumo_fidelity', 1.0):.5f}  (|<psi_PBE|psi_QP>|^2)")
        elif details.get("self_consistent"):
            logger.info(f"    evGW Iterations          : {details.get('evgw_iterations', 1)} (converged: {details.get('evgw_converged', False)})")
        logger.info(f"    Solvent epsilon_out      : {details['eps_out']:.2f}")
        logger.info(f"    Total scissor            : {details['total_scissor_ev']:+.3f} eV")
        if dft_gap is not None and target_qp_gap is not None:
            logger.info(f"    Final gap                : {dft_gap:.3f} -> {target_qp_gap:.3f} eV")
        if output_file:
            logger.info(f"    JSON                     : {output_file}")
        return

    if details.get("qp_model") == "bulk_gw_scissor":
        logger.info("\n  [QP Provenance] Bulk GW Scissor Model (Bypassing Finite-Size QP)")
        logger.info(f"    Material                 : {details.get('material', 'N/A')}")
        logger.info(f"    Bulk PBE -> GW gap       : {details.get('bulk_pbe_gap_ev', 0.0):.3f} -> {details.get('bulk_gw_gap_ev', 0.0):.3f} eV")
        logger.info(f"    Bulk GW scissor          : {details.get('bulk_gw_shift_ev', 0.0):+.3f} eV")
        logger.info("    Finite-size Delta-W      : none (pure bulk limit)")
        if dft_gap is not None and target_qp_gap is not None:
            logger.info(f"    Final gap                : {dft_gap:.3f} -> {target_qp_gap:.3f} eV")
        if output_file:
            logger.info(f"    JSON                     : {output_file}")
        return

    logger.info("\n  [QP Provenance] Anchor-scaled PBE-to-QP model")
    logger.info(f"    Material                 : {details['material']}")
    if "cluster_radius_ang" in details:
        logger.info(f"    Cluster radius           : {details['cluster_radius_ang']:.3f} Å")
    logger.info(f"    Bulk PBE -> GW gap       : {details['bulk_pbe_gap_ev']:.3f} -> {details['bulk_gw_gap_ev']:.3f} eV")
    logger.info(f"    Bulk GW shift            : {details['bulk_gw_shift_ev']:+.3f} eV")
    if details.get("periodic_bulk_limit"):
        logger.info("    Periodic mode            : using tabulated bulk GW-PBE scissor only")
    if details.get("has_monomer_anchor"):
        logger.info(f"    Monomer anchor radius    : {details['monomer_radius_ang']:.3f} Å")
        logger.info(f"    Monomer PBE -> GW gap    : {details['monomer_pbe_gap_ev']:.3f} -> {details['monomer_gw_gap_ev']:.3f} eV")
        logger.info(f"    Anchor residual A        : {details['anchor_residual_ev']:+.3f} eV")
        if details.get("polarization_model", "legacy") == "sphere":
            logger.info(f"    Polarization             : sphere, F = {details['polarization_factor_solvent']:.4f} "
                  f"(vacuum {details['polarization_factor_vacuum']:.4f}), p = {details['residual_power']:.3f}")
        else:
            logger.info(f"    ell, p                   : {details['regularization_length_ang']:.3f} Å, {details['residual_power']:.3f}")
        if details.get("edge_split_source") == "anchor_edge_curves":
            logger.info(f"    Edge residuals A_h, A_l  : {details['edge_residual_homo_ev']:+.3f}, "
                  f"{details['edge_residual_lumo_ev']:+.3f} eV (HOMO share of bulk shift {details['bulk_homo_fraction']:.2f})")
    if details.get("principal_extents_ang"):
        extents = ", ".join(f"{x:.3f}" for x in details["principal_extents_ang"])
        logger.info(f"    Principal extents        : [{extents}] Å")
        logger.info(f"    Anisotropy ratio         : {details['anisotropy_ratio']:.3f}")
    logger.info(f"    Vacuum finite-size shift : {details['finite_size_shift_vacuum_ev']:+.3f} eV")
    logger.info(f"    Solvent finite-size shift: {details['finite_size_shift_solvent_ev']:+.3f} eV")
    logger.info(f"    Vacuum scissor           : {details['total_scissor_vacuum_ev']:+.3f} eV")
    logger.info(f"    Solvent scissor used     : {details['total_scissor_solvent_ev']:+.3f} eV")
    if dft_gap is not None and target_qp_gap is not None:
        logger.info(f"    Final gap                : {dft_gap:.3f} -> {target_qp_gap:.3f} eV")
    if output_file:
        logger.info(f"    JSON                     : {output_file}")


def write_qp_provenance(details, dft_gap, target_qp_gap, scissor, args, filename="qp_provenance.json"):
    if not details:
        return None

    payload = {}
    if getattr(args, "cluster_size_info", None):
        payload["cluster_size"] = args.cluster_size_info
    for k, v in details.items():
        if isinstance(v, np.ndarray):
            if v.size <= 10:
                payload[k] = v.tolist()
        else:
            payload[k] = v
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
    logger.info(f"\n===================================================")
    logger.info(f" [ {label} ] EXCITON CALCULATION ")
    logger.info(f"===================================================")

    start_solve = time.time()
    energies_ev, vectors = solver.solve(
        nroots=args.nroots, full_diag=args.full_diag, tol=args.tol,
        excitation_mode=args.excitation_mode,
    )
    
    if args.soc != 0.0:
        energies_ev = energies_ev - args.soc
        logger.info(f"  [SOC Shift] Applied empirical energy shift: -{args.soc:.3f} eV")
    logger.debug(f"  {label} Solver converged in {time.time() - start_solve:2.2f} s")

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

    logger.info("\n" + "-"*60)
    logger.info(f" SYSTEM ENERGY SUMMARY ({label})")
    logger.info("-" * 60)
    logger.info(f"  Raw DFT Gap           : {dft_gap:8.4f} eV")
    if solver.soc_flag and soc_gap is not None:
        logger.info(f"  SOC Gap               : {soc_gap:8.4f} eV")
    logger.info(f"  QP Correction (Shift) : {scissor:8.4f} eV")
    logger.info(f"  Confinement Energy    : {confinement_energy:8.4f} eV")
    logger.info(f"  Excitation Mode        : {args.excitation_mode}")
    if hasattr(solver, 'eps_info') and solver.eps_info:
        eps_info = solver.eps_info
        logger.info(f"  Microscopic ε_eff (1S) : {eps_info.get('eps_eff_exciton', 1.0):8.3f} (Lowest exciton screening)")
        if eps_info.get('eps_bulk'):
            logger.info(f"  Bulk Dielectric (ε_∞)  : {eps_info.get('eps_bulk'):8.3f} (Retention: {eps_info.get('dielectric_retention_pct', 100.0):.1f}%)")
        dielectric_file = f"dielectric_summary{suffix}.json"
        try:
            with open(dielectric_file, "w", encoding="utf-8") as f:
                json.dump(eps_info, f, indent=2, sort_keys=True)
            logger.info(f"  Dielectric Log File    : {dielectric_file}")
        except Exception:
            pass
    elif args.kernel != "resta":
        logger.info(f"  Legacy Kernel Scaling  : {args.alpha:8.4f}")
    logger.info("-" * 60)

    logger.info("\n" + "="*172)
    logger.info(f"{'State':>5} {'Energy':>10} {'Main Trans':>12} {'Weight':>8} {'f_osc':>10} | {'PR':>5} | {'D(eV)':>8} {'Kx(eV)':>8} {'-Kd(eV)':>8} | {'Spin-Free Character':>62}")
    logger.info("-" * 172)

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

        logger.info(f"{n+1:5d} {energies_ev[n]:10.4f}  {trans_str:>12}  {weight**2:8.3f}  {f_strengths[n]:10.5f} | {pr:5.1f} | {dE_val:8.4f} {Kx_val:8.4f} {minus_Kd_val:8.4f} | {spin_str:>62}")

    if len(energies_ev) > 100: logger.info(f" ... {len(energies_ev) - 100} additional states computed (output truncated) ...")
    logger.info("="*172)

    logger.info(f"\n--- Performing Dreuw/Plasser Analysis ({label}) ---")
    analyzer = ExcitonAnalyzer(solver, np.array(coords_ang), syms)
    analysis_results = []
    analysis_t0 = time.perf_counter()

    logger.info(f"{'State':>5} {'Energy':>8} {'f_osc':>8} | {'PR':>5} {'d_eh~(A)':>8} {'d_CT~(A)':>8} {'sig_h~':>7} {'sig_e~':>7} | {'Type':>8}")
    logger.info("-" * 95)

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

        logger.info(f"{n+1:5d} {res['energy']:8.3f} {res['f_osc']:8.4f} | {res['PR']:5.1f} {res['d_eh']:7.2f} {res['d_CT']:7.2f} {res['sigma_h']:6.1f} {res['sigma_e']:6.1f} | {ex_type:>8}")

    logger.debug(
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
            
        logger.info(f"\n--- Generating Cubes ({n_bse} Excitons) ---")
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
        logger.info(f"  Saved X_ia coefficients to {npz_filename}")

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
        eff_scissor = getattr(solver, "scissor_ev", None)
        if eff_scissor is None:
            eff_scissor = scissor
        if eff_scissor is not None and eff_scissor != 0.0:
            eps_eval[solver.homo_index + 1:] += eff_scissor
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


def _apply_config(args, config_data, explicit_cli_args=None):
    if explicit_cli_args is None:
        explicit_cli_args = set()
    for section, parameters in config_data.items():
        if section == "periodic" and isinstance(parameters, dict):
            if "periodic_enabled" not in explicit_cli_args and "periodic" not in explicit_cli_args:
                setattr(args, "periodic_enabled", bool(parameters.get("enabled", False)))
            if "lattice_vectors" in parameters and "lattice_vectors" not in explicit_cli_args:
                setattr(args, "lattice_vectors", parameters["lattice_vectors"])
            if "overlap_cutoff" in parameters and "overlap_cutoff" not in explicit_cli_args:
                setattr(args, "overlap_cutoff", parameters["overlap_cutoff"])
            continue

        if section == "namd" and isinstance(parameters, dict):
            setattr(args, "namd_cfg", parameters)
            continue

        if section == "auger" and isinstance(parameters, dict):
            if "auger" not in explicit_cli_args:
                setattr(args, "auger", bool(parameters.get("run", parameters.get("enabled", True))))
            if "sigma" in parameters and "auger_sigma" not in explicit_cli_args:
                setattr(args, "auger_sigma", float(parameters["sigma"]))
            if "channel" in parameters and "auger_channel" not in explicit_cli_args:
                setattr(args, "auger_channel", str(parameters["channel"]))
            if "n_initial_states" in parameters and "auger_states" not in explicit_cli_args:
                setattr(args, "auger_states", int(parameters["n_initial_states"]))
            if "lineshape" in parameters and "auger_lineshape" not in explicit_cli_args:
                setattr(args, "auger_lineshape", str(parameters["lineshape"]))
            if "eps_eff" in parameters and parameters["eps_eff"] is not None and "auger_eps_eff" not in explicit_cli_args:
                setattr(args, "auger_eps_eff", float(parameters["eps_eff"]))
            continue

        if isinstance(parameters, dict):
            if section == "physics":
                logger.info("  [Config] The 'physics' section is the old layout; it still works. New layout: "
                      "quasiparticles, environment, integrals, excitations (see docs, Configuration).")
            for key, value in parameters.items():
                dest = section_dest(section, key)
                if dest is None:
                    valid = ", ".join(sorted(SECTIONS[section]))
                    raise ValueError(f"Unknown YAML key '{section}.{key}'. Valid keys: {valid}")
                if (
                    key in explicit_cli_args
                    or dest in explicit_cli_args
                    or (dest == "kernel_type" and any(k in explicit_cli_args for k in ["kernel_type", "2e_integrals", "two_electron_integrals"]))
                ):
                    continue
                if hasattr(args, dest):
                    setattr(args, dest, value)
                else:
                    raise ValueError(f"Unknown YAML key '{section}.{key}'")
        else:
            norm_sec = section.replace("-", "_")
            if section in explicit_cli_args or norm_sec in explicit_cli_args:
                continue
            if hasattr(args, section):
                setattr(args, section, parameters)
            elif hasattr(args, norm_sec):
                setattr(args, norm_sec, parameters)
            else:
                raise ValueError(f"Unknown YAML key '{section}'")


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
        

W_BASED_QP_MODELS = {
    "sgw-dim", "sgw_dim", "sgw-dim-dz", "evgw", "evgw-dim", "evgw_dim", "qsgw", "qsgw-dim", "qsgw_dim",
    "scgw", "scgw-dim", "scgw_dim",
    "sgw-resta", "sgw_resta", "sgw-resta-penn", "sgw-resta-pure", "sgw-resta-bulk", "sgw-resta-dz",
    "evgw-resta", "evgw_resta", "qsgw-resta", "qsgw_resta", "scgw-resta", "scgw_resta",
    "sgw", "sgw-dw", "sgw_dw", "sgw-atom", "sgw-ao",
}


def resolve_qp_z(args, model_default_dynamic=True):
    """Return (dynamic_z, Z) for a Delta-W QP model from --qp-z / --dynamic_z.

    Without an explicit qp_z, Z is derived from the model's own screening
    (plasmon pole, see ``compute_dynamic_z``) for every Delta-W model.
    ``model_default_dynamic`` is kept for callers that need the old default.
    """
    qz = getattr(args, "qp_z", None)
    if qz is None:
        if getattr(args, "dynamic_z", False):
            return True, 0.8
        return bool(model_default_dynamic), 0.8
    if str(qz).strip().lower() == "derived":
        return True, 0.8
    try:
        z = float(qz)
    except ValueError:
        raise ValueError(f"qp_z must be 'derived' or a number, got '{qz}'")
    if not 0.0 < z <= 1.0:
        raise ValueError(f"qp_z = {z} is outside (0, 1]; a quasiparticle weight must lie in that range")
    return False, z


def resolve_bse_kernel(args, qp_w):
    """Choose the BSE direct kernel so that GW and BSE use the same W.

    Delta-W QP models (``sgw-*``, ``evgw-*``, ``qsgw-*``, ``sgw``) define a screened interaction
    W; the BSE must use that same W (kernel 'qp').  Gap-only models (pbe, brus,
    gw, numeric) define no W, so any interior kernel may be chosen.
    """
    qp_name = str(args.qp_gap).lower()
    kernel = args.kernel
    if str(getattr(args, "excitation_mode", "")).lower() in ("stda", "diagonal_stda"):
        if qp_w is not None:
            raise ValueError(f"excitation_mode '{args.excitation_mode}' uses the sTDA interactions; qp_gap "
                             f"'{args.qp_gap}' defines its own W. Use quasiparticles.model: none (faithful sTDA) "
                             "or bulk (dielectric variant).")
        if kernel not in (None, "stda"):
            raise ValueError(f"excitation_mode '{args.excitation_mode}' sets the kernel itself; remove kernel '{kernel}'.")
        if str(args.kernel_type).lower() != "mnok":
            raise ValueError("sTDA is a monopole (MNOK-type) method: use integrals.representation: mnok.")
        return "stda"
    if qp_w is not None:
        label = qp_w[1]
        equivalent = kernel in ("sbse", "sbse-atom") and qp_name in ("sgw", "sgw-dw", "sgw_dw", "sgw-atom")
        if kernel in (None, "qp") or equivalent:
            if str(args.kernel_type).lower() not in ("mnok", "xs", "xs-qdex"):
                raise ValueError(f"kernel 'qp' supports two_electron_integrals mnok or xs, not '{args.kernel_type}'.")
            logger.info(f"  [Consistency] BSE direct kernel = W of the QP model ({label}).")
            return "qp"
        if args.allow_inconsistent_kernel:
            logger.warning(f"  [Consistency WARNING] qp_gap '{args.qp_gap}' uses W = {label}, but the BSE uses kernel "
                  f"'{kernel}'. GW and BSE do not share W (allowed by --allow-inconsistent-kernel).")
            return kernel
        raise ValueError(
            f"qp_gap '{args.qp_gap}' is built on its own screened interaction ({label}). The BSE must use the same "
            f"W: set kernel: qp (the default for this model). kernel '{kernel}' would use a different W in GW and "
            "BSE. Use --allow-inconsistent-kernel only to reproduce legacy results."
        )
    if kernel == "qp":
        raise ValueError(
            f"kernel 'qp' needs a QP model that defines W (sgw-dim, sgw-resta, evgw-*, qsgw-*, sgw); "
            f"qp_gap '{args.qp_gap}' is a gap-only model. Choose resta, xs-resta, dim, sbse or bse."
        )
    if kernel is None:
        kernel = "bse"
    logger.info(f"  [Consistency] qp_gap '{args.qp_gap}' is gap-only (no W); the BSE uses the independent kernel '{kernel}'.")
    return kernel


def _build_parser():
    """Command-line interface of QDEX."""
    parser = argparse.ArgumentParser(description="QDEX - Quantum Dot Excitations & Dynamics exciton solver")

    parser.add_argument("--config", type=str, help="Path to a YAML configuration file.")
    parser.add_argument("--mo_file")
    parser.add_argument("--xyz") 
    parser.add_argument("--basis_txt")
    parser.add_argument("--basis_name", help="CP2K basis name, or 'per-atom' for a per-atom basis file (g-xTB frames, qdex.xtb.molden).")
    parser.add_argument(
        "--cache-mos", dest="cache_mos", action="store_true",
        help="Cache parsed text MOs as a validated uncompressed binary NPZ sidecar for faster repeated runs.",
    )

    parser.add_argument("--n-occ", type=int, default=50)
    parser.add_argument("--n-virt", type=int, default=50)
    parser.add_argument("--e_thresh", type=float, default=None)
    parser.add_argument("--f_thresh", type=float, default=0.0)

    parser.add_argument("--qp_gap", type=str, default="brus",
                        help="Quasiparticle gap model: 'sgw-anchor' / 'gw' (anchor-scaled), 'sgw-dim' (atomistic polarizable dipole Delta-W), 'evgw-dim' / 'evgw' (DIM gap/screening fixed-point iteration), 'qsgw-dim' / 'qsgw' (static DIM Delta-COHSEX orbital-relaxation model, full AO update), 'sgw-resta' (Resta Penn-scaled Delta-W), 'evgw-resta' (Resta gap/screening fixed-point iteration), 'qsgw-resta' (static Resta Delta-COHSEX orbital-relaxation model, full AO update), 'sgw-resta-pure' (Resta boundary Delta-W), 'sgw' (site-diagonal Delta-W on the sBSE-screened kernel), 'brus' (bulk experimental gap + effective-mass kinetic confinement), 'bulk' (PBE orbitals + bulk GW correction, PBE only), 'none' / 'pbe' / 'dft' (uncorrected DFT energies), or explicit gap in eV.")
    parser.add_argument("--bulk-vertex", dest="bulk_vertex", choices=["none", "full", "scaled"], default="none",
                        help="Vertex correction of the bulk QSGW shift: 'none' (pure QSGW, default), 'full' (bulk factor "
                             "at every size) or 'scaled' (times the Penn fraction of bulk screening the dot keeps).")
    parser.add_argument("--bulk-vertex-factor", dest="bulk_vertex_factor", type=float, default=0.8,
                        help="Bulk vertex factor: Delta_bulk -> factor * Delta_bulk in the bulk (default 0.8).")
    parser.add_argument("--qp-reference", dest="qp_reference", choices=["pbe", "gxtb"], default="pbe",
                        help="Orbitals the 'bulk' QP shift corrects: 'pbe' (default, bulk QSGW - bulk PBE) or 'gxtb' "
                             "(g-xTB frames: spin-free experimental gap - g-xTB bulk gap; NAMD precompute only).")
    parser.add_argument("--qp-z", dest="qp_z", type=str, default=None,
                        help="Quasiparticle renormalization Z for Delta-W models: 'derived' (default; one plasmon pole whose "
                             "frequency follows from the same eps as the model's W, see compute_dynamic_z) or a fixed "
                             "number, e.g. 1.0. Not applicable to gap-only models.")
    parser.add_argument("--allow-inconsistent-kernel", action="store_true", default=False,
                        help="Allow a BSE kernel whose W differs from the W of the QP model (legacy behaviour; "
                             "only for reproducing old results).")
    parser.add_argument("--dynamic_z", action="store_true", default=False,
                        help="Deprecated alias of --qp-z derived (now the default for Delta-W models).")
    parser.add_argument("--update_orbitals", "--qsgw", dest="update_orbitals", action="store_true", default=False,
                        help="Perform full AO-basis Quasiparticle Self-Consistent GW (qsGW) orbital update.")
    parser.add_argument("--soc", type=float, default=0.0)
    parser.add_argument("--soc_flag", action="store_true")
    parser.add_argument("--gth_file", type=str, default=None)

    parser.add_argument("--kernel", choices=["qp", "stda", "resta-sphere", "bse", "resta", "mnok", "xs", "xs-resta", "xs-qdex", "xs-rpa", "rpa", "dim", "xs-dim", "dipole", "xs-dipole", "sbse", "sbse-atom", "sbse-ao", "xs-sbse"], default=None,
                        help="Exciton interaction kernel. 'qp': the screened W built by the QP model (default and only allowed choice for "
                             "sgw-*, evgw-*, qsgw-* and sgw, so that GW and BSE share one W). For gap-only QP models (pbe, brus, gw, "
                             "numeric) choose: 'bse' (MNOK uniform, legacy default), 'resta' (MNOK Resta-screened), 'dim' / 'dipole' (MNOK atomistic polarizable dipole model), 'sbse' / 'sbse-atom' (Simplified BSE atom-resolved kernel, Cho et al. 2022), 'sbse-ao' / 'xs-sbse' (Simplified BSE AO-resolved kernel), 'xs' (Xs-QDEX uniform), 'xs-resta' (Xs-QDEX Resta-screened), 'xs-rpa' (Xs-QDEX microscopic RPA screening), 'xs-dim' (Xs-QDEX atomistic polarizable dipole model).")
    parser.add_argument("--two-electron-integrals", "--two_electron_integrals", "--2e-integrals", "--2e_integrals", "--kernel_type", "--kernel-type",
                        dest="kernel_type", choices=["mnok", "xs", "xs-qdex"], default="mnok",
                        help="Two-electron integral representation: 'mnok' (semi-empirical atom-centered damped Coulomb) or 'xs' / 'xs-qdex' (exact analytical Gaussian AO four-center integrals via Libint2).")
    parser.add_argument("--alpha", type=float, default=1.0, help="Scaling for the legacy non-RESTA kernel; ignored by RESTA.")
    parser.add_argument("--beta", type=float, default=0.0, help="Reserved on-site stiffening parameter; beta > 0 is currently rejected.")
    parser.add_argument("--exchange", action="store_true", default=None, help="Deprecated alias for --include-direct-eh.")
    direct_group = parser.add_mutually_exclusive_group()
    direct_group.add_argument("--include-direct-eh", dest="include_direct_eh", action="store_true", default=None,
                              help="Include the Resta-screened attractive electron-hole direct term (default).")
    direct_group.add_argument("--no-direct-eh", dest="include_direct_eh", action="store_false",
                              help="Disable the attractive electron-hole direct term.")
    exchange_group = parser.add_mutually_exclusive_group()
    exchange_group.add_argument("--include-exchange", dest="include_exchange", action="store_true", default=None,
                                help="Include the repulsive electron-hole exchange term Kx (default).")
    exchange_group.add_argument("--no-exchange", dest="include_exchange", action="store_false",
                                help="Disable the repulsive electron-hole exchange term Kx.")
    parser.add_argument("--estimate_qp", action="store_true", help="Compute G0W0-lite Quasiparticle corrections via COHSEX")
    parser.add_argument("--use_cohsex_gap", action="store_true", help="Override the tabulated GW gap with the pure COHSEX computed gap")
    parser.add_argument("--vxc_ao", type=str, default=None, help="Path to cleaned CP2K AO-basis Vxc matrix text file")
    parser.add_argument("--material", type=str, default="DEFAULT")
    parser.add_argument("--eps-out", type=float, default=2.0)
    parser.add_argument("--qp-regularization-length", dest="qp_regularization_length", type=float, default=1.0,
                        help="Regularization length ell in angstrom for the anchor-scaled QP model.")
    parser.add_argument("--qp-residual-power", dest="qp_residual_power", type=float, default=2.0,
                        help="Power p > 1 controlling decay of the anchor residual.")
    parser.add_argument("--qp-polarization", dest="qp_polarization", choices=["sphere", "legacy"], default="sphere",
                        help="Finite-size term of the anchor-scaled QP model: 'sphere' (classical surface polarization "
                             "of a dielectric sphere, 1S-averaged, Delerue-Lannoo-Allan; default) or 'legacy' "
                             "(11.52 eV A (1/eps_out - 1/eps_inf)/(R + ell)).")
    parser.add_argument("--qp-edge-split", dest="qp_edge_split", choices=["anchor", "model"], default="model",
                        help="HOMO/LUMO split of the QP correction for absolute IP/EA: 'model' (the model's own "
                             "levels, default) or 'anchor' (legacy per-edge curves fitted to the monomer evGW).")
    parser.add_argument("--qp-levels", dest="qp_levels", choices=["orbital", "rigid"], default="orbital",
                        help="QP energies of models that define W: 'orbital' (default; every molecular orbital gets its "
                             "own Z_p * sigma_p) or 'rigid' (one scissor on all virtual orbitals).")
    parser.add_argument("--qp-selfenergy", dest="qp_selfenergy", choices=["cohsex", "classical"], default="cohsex",
                        help="Orbital QP levels of the Delta-W models: 'cohsex' (default; one-shot static Delta-COHSEX "
                             "diagonal incl. non-classical screened exchange) or 'classical' (1/2 q^T dW q).")
    parser.add_argument("--qp-window", dest="qp_window", choices=["active", "all"], default="active",
                        help="Orbital space evaluated in orbital Delta-W / COHSEX self-energy: 'active' "
                             "(default; evaluates states covering the BSE active space and clamps edges, fast) "
                             "or 'all' (evaluates every orbital in the basis).")
    parser.add_argument("--qp-window-size", dest="qp_window_size", type=int, default=None,
                        help="Optional override for the number of active occupied/virtual states around the "
                             "Fermi level evaluated in orbital QP mode (default: matches BSE active space, min 25).")
    parser.add_argument("--qp-solvent-term", dest="qp_solvent_term", choices=["sphere", "born"], default="sphere",
                        help="Environment part of Delta W in the Delta-W models: 'sphere' (default; reaction field "
                             "of a dielectric sphere, as in gw) or 'born' (earlier softened Born form).")
    parser.add_argument("--qp-residual-scaling", dest="qp_residual_scaling", choices=["econf", "power"],
                        default="econf",
                        help="Size scaling of the non-classical anchor residual: 'econf' (default; "
                             "E_conf^PBE(R)/E_conf^PBE(R0)) or 'power' ((R0/R)^p).")
    parser.add_argument("--qp-anchor-residual", dest="qp_anchor_residual", choices=["on", "off"], default=None,
                        help="Add the residual calibrated on the monomer evGW: default off for the Delta-W models "
                             "(the approximation is used as it is); the two-anchor gw model always has it.")
    parser.add_argument("--qp-anchor-calibrate", dest="qp_anchor_calibrate", action="store_true",
                        help="Run on the anchor cluster (vacuum): store the per-edge residual of this Delta-W "
                             "model against the evGW anchor in the residual table.")
    parser.add_argument("--qp-anchor-table", dest="qp_anchor_table", default=None,
                        help="Path of the anchor residual table (default: qdex/data/dw_anchor_residuals.json).")
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
    parser.add_argument("--verbosity", choices=["full", "normal", "quiet"], default="full",
                        help="Console output: 'full' (default, everything), 'normal' (no iteration traces, timings "
                             "or diagnostics) or 'quiet' (warnings and errors). The log file always gets everything.")

    parser.add_argument("--nto", action="store_true", help="Run Natural Transition Orbital analysis after solving excitons.")
    parser.add_argument("--nto-states", type=int, nargs='+', help="Specific exciton states for NTO analysis, 1-indexed.")
    parser.add_argument("--nto-top", type=int, default=3, help="Number of dominant NTO pairs to print per state.")
    parser.add_argument("--nto-csv", action="store_true", help="Write detailed NTO descriptors to nto_results*.csv.")
 
    parser.add_argument("--nroots", type=int, default=10)
    parser.add_argument("--full-diag", action="store_true")
    parser.add_argument(
        "--excitation-mode",
        choices=["bse", "independent_dft", "independent_qp", "diagonal_bse", "sbse", "diagonal_sbse",
                 "stda", "diagonal_stda"],
        default="bse",
        help=("Excitation model: diagonalize BSE/TDA ('bse', 'sbse' or Grimme's 'stda'), or use uncoupled transitions "
              "with DFT gaps ('independent_dft'), QP gaps ('independent_qp'), or diagonal Kx/Kd corrections "
              "('diagonal_bse', 'diagonal_sbse', 'diagonal_stda')."),
    )
    parser.add_argument("--selection", dest="selection", choices=["none", "perturbative"], default="none",
                        help="Transition selection for the coupled solvers: 'perturbative' keeps the transitions "
                             "below selection_energy and adds the higher ones strongly coupled to them (Grimme).")
    parser.add_argument("--selection-energy", dest="selection_energy", type=float, default=7.0,
                        help="Primary-space threshold E_thr in eV on the diagonal A_ia,ia (default 7.0, as std2).")
    parser.add_argument("--selection-pt", dest="selection_pt", type=float, default=1e-4,
                        help="Second-order selection threshold t in hartree (default 1e-4, as std2).")
    parser.add_argument("--selection-shift", dest="selection_shift", choices=["on", "off"], default="on",
                        help="Lower the primary diagonal by the second-order contributions of the rejected "
                             "transitions (std2 behaviour, default on).")
    parser.add_argument("--mnok-exponent", dest="mnok_exponent", type=float, default=2.0,
                        help="Exponent beta of the MNOK interaction (r^beta + a^beta)^(-1/beta) for all interactions "
                             "(2 = Ohno-Klopman, default; 1 = Mataga-Nishimoto).")
    parser.add_argument("--mnok-exponent-exchange", dest="mnok_exponent_exchange", type=float, default=None,
                        help="Exponent of the exchange interaction K^x only (default: mnok_exponent).")
    parser.add_argument("--mnok-onsite", dest="mnok_onsite", choices=["eta", "ip_ea"], default="ip_ea",
                        help="On-site value of the MNOK interaction: 'ip_ea' (IP - EA = 2 eta, the monopole (ii|ii); "
                             "default, reproduces the exact xs exchange) or 'eta' (Ghosh-Islam hardness, (IP-EA)/2; "
                             "the earlier convention).")
    parser.add_argument("--qp-radius", dest="qp_radius", choices=["saxs", "hull"], default="saxs",
                        help="Radius of the dielectric sphere and of the Brus confinement: 'saxs' (default, SAXS-"
                             "equivalent sphere of the inorganic electron density) or 'hull' (core hull + 1.25 A). "
                             "The two-anchor gw model always uses the hull radius.")
    parser.add_argument("--inorganic-elements", dest="inorganic_elements", nargs="+", default=None,
                        help="Elements seen by SAXS for the reported size (default: all but H, C, N, O, P, B, Si, F).")
    parser.add_argument("--stda-functional", dest="stda_functional", type=str, default=None,
                        help="sTDA: functional of the MO file, sets a_x (pbe 0, b3lyp 0.20, pbe0 0.25, ...); a "
                             "range-separated one (cam-b3lyp, wb97x-d3, wb97m-v, gxtb) or stda-xtb/gfn2 sets a_x, alpha and beta.")
    parser.add_argument("--stda-ax", dest="stda_ax", type=str, default=None,
                        help="sTDA: explicit Fock-exchange fraction a_x, or 'dielectric' for 1/eps_inf of the material.")
    parser.add_argument("--stda-alpha", dest="stda_alpha", type=float, default=None,
                        help="sTDA: explicit alpha of gamma^K (default 1.42 + 0.48 a_x, or the range-separated set).")
    parser.add_argument("--stda-beta", dest="stda_beta", type=float, default=None,
                        help="sTDA: explicit beta of gamma^J (default 0.20 + 1.83 a_x, or the range-separated set).")
    parser.add_argument("--tol", type=float, default=1e-5)
    parser.add_argument("--skip-orthonormality-check", dest="skip_orthonormality_check", action="store_true",
                        help="Skip the C^T S C = I check (a full n_ao^3 product); use only for MO files already "
                             "known to be orthonormal (SOC still re-orthonormalizes its small active window). YAML: "
                             "system.skip_orthonormality_check: true")
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
    return parser


def _load_arguments(parser):
    """Parse the command line and the YAML config; run NAMD/utility modes (returns None then)."""
    args = parser.parse_args()

    config_path = args.config
    config_data = {}
    if args.config:
        with open(args.config, 'r') as f:
            config_data = yaml.safe_load(f) or {}

        explicit_cli_args = set()
        for arg_str in sys.argv[1:]:
            if arg_str.startswith("--"):
                flag_name = arg_str.lstrip("-").split("=")[0].replace("-", "_")
                explicit_cli_args.add(flag_name)

        _apply_config(args, config_data, explicit_cli_args=explicit_cli_args)

    if getattr(args, "namd_soc", False):
        config_data.setdefault("soc", {})["enabled"] = True
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
        setup_run_logging(getattr(args, "log_file", "minibse.log"), getattr(args, "verbosity", "full"))
        precompute_namd_data(config_data)
        return

    if getattr(args, "namd_run", False):
        from qdex.namd import run_namd_dynamics
        setup_run_logging(getattr(args, "log_file", "minibse.log"), getattr(args, "verbosity", "full"))
        run_namd_dynamics(config_data)
        return

    if getattr(args, "lattice_vectors", None) is not None:
        args.lattice_vectors = np.asarray(args.lattice_vectors, dtype=float).reshape(3, 3).tolist()

    if getattr(args, "include_direct_eh", None) is None:
        args.include_direct_eh = True if args.exchange is None else bool(args.exchange)
    if args.exchange is not None:
        logger.warning("  [Deprecated] 'exchange' now maps to include_direct_eh; use include_direct_eh instead.")
    if getattr(args, "include_exchange", None) is None:
        args.include_exchange = True
    return {"args": args, "config_path": config_path}


def _prepare_run(args, *, config_path, parser):
    """Logging, argument validation, MNOK/radius settings and the compute device."""
    setup_run_logging(getattr(args, "log_file", "minibse.log"), getattr(args, "verbosity", "full"))
    validate_args(args, parser)
    import qdex.hardness as _hardness
    _hardness.RADIUS_DEFINITION = str(getattr(args, "qp_radius", "saxs") or "saxs").lower()
    _hardness.set_bulk_vertex(getattr(args, "bulk_vertex", "none"), getattr(args, "bulk_vertex_factor", 0.8))
    _hardness.set_mnok_options(getattr(args, "mnok_exponent", 2.0), getattr(args, "mnok_exponent_exchange", None),
                               getattr(args, "mnok_onsite", "ip_ea"))
    compute_device, dev_obj = resolve_device(args.device, verbose=True)
    if config_path:
        logger.info(f"Loading configuration from {config_path}...")

    if not getattr(args, "xyz", None) and geometry_source(None, getattr(args, "mo_file", None)):
        args.xyz = args.mo_file   # geometry from the TREXIO 'nucleus' group of the HDF5 MO file
        logger.info(f"  [Geometry] system.xyz not set; reading the geometry from '{args.mo_file}'.")
    required_args = ['mo_file', 'xyz', 'basis_txt', 'basis_name', 'qp_gap']
    missing = [arg for arg in required_args if getattr(args, arg) is None]
    if missing: parser.error(f"Missing required arguments: {', '.join(missing)}")

    logger.info("\n===================================================")
    logger.info(" QDEX - Post-DFT Exciton Solver")
    logger.info("===================================================")

    tracker = ResourceTracker()
    return _export(locals(), (
        "compute_device", "dev_obj", "tracker"
    ))


def _read_geometry_and_basis(args, *, tracker):
    """Geometry, basis set and the size of the dot (always reported)."""
    tracker.start_stage("Geometry & Basis Parsing")
    logger.info("\n--- Parsing Geometry and Basis Set ---")
    t0_parse = time.time()
    mo_path = str(getattr(args, "mo_file", "") or "")
    is_h5_mo = os.path.exists(mo_path) and (mo_path.lower().endswith((".h5", ".hdf5")) or is_h5_file(mo_path))

    if getattr(args, "xyz", None) and not os.path.exists(args.xyz) and is_h5_mo:
        logger.info(f"  [Geometry] '{args.xyz}' not found on disk; extracting nuclear geometry directly from HDF5 '{mo_path}'")
        syms, coords_ang = read_geometry_h5(mo_path)
    else:
        syms, coords_ang = read_xyz(args.xyz)
        if is_h5_mo:
            try:
                syms_h5, coords_h5 = read_geometry_h5(mo_path)
                if len(syms_h5) == len(syms):
                    max_dev = float(np.max(np.abs(np.asarray(coords_ang) - coords_h5)))
                    if max_dev > 1e-4:
                        logger.warning(
                            f"  [Geometry:Warning] Coordinates in '{args.xyz}' differ from '{mo_path}' "
                            f"by up to {max_dev:.4f} Å. Using the self-consistent nuclear coordinates "
                            f"from '{mo_path}' to ensure exact basis set centering and MO orthonormality."
                        )
                        syms, coords_ang = syms_h5, coords_h5
            except Exception:
                pass
    if str(args.basis_name).lower() == "per-atom":
        # One basis block per atom (qdex.xtb.molden): the g-xTB basis depends on the atomic charge
        from qdex.xtb.molden import build_frame_shells
        shells = build_frame_shells(args.basis_txt, syms, coords_ang)
        logger.info(f"  [Basis] per-atom basis file '{args.basis_txt}'")
    else:
        basis_dict = parse_basis(args.basis_txt, args.basis_name, required_elements=set(syms))
        shells = build_shell_dicts(syms, coords_ang, basis_dict)
    shells = [{**sh, 'pure': True} for sh in shells] # Use Sphericals
    n_ao = count_ao_from_shells(shells)
    atom_ao_ranges = build_atom_ao_ranges(shells)
    logger.debug(f"  -> Parsed in {time.time() - t0_parse:.2f} s | Total AOs: {n_ao}")

    # Size of the dot, always reported (SAXS definition), whether or not the model uses it.
    cluster_size_info = None
    try:
        from qdex.cluster_size import cluster_size, format_cluster_size
        cluster_size_info = cluster_size(np.array(coords_ang), syms, args.material,
                                         getattr(args, "inorganic_elements", None))
        qp_name_sz = str(args.qp_gap).lower()
        uses_radius = not (qp_name_sz in ("none", "pbe", "dft", "bulk")
                           or qp_name_sz.replace(".", "", 1).isdigit())
        if not uses_radius:
            radius_note = None
        elif qp_name_sz in ("gw", "sgw-anchor") or str(getattr(args, "qp_radius", "saxs")).lower() == "hull":
            radius_note = (f"hull radius {cluster_size_info['hull_radius_ang']:.3f} A "
                           f"(sphere reaction field / confinement of '{args.qp_gap}')")
        else:
            radius_note = (f"SAXS radius {5.0 * cluster_size_info['d_saxs_nm']:.3f} A "
                           f"(sphere reaction field / confinement of '{args.qp_gap}')")
        logger.info(format_cluster_size(cluster_size_info, radius_note))
        with open("cluster_size.json", "w", encoding="utf-8") as fh:
            json.dump(cluster_size_info, fh, indent=2)
        args.cluster_size_info = cluster_size_info
    except Exception as exc:  # the size report must never stop a calculation
        logger.info(f"  [Size] Could not evaluate the cluster size: {exc}")
    return _export(locals(), (
        "atom_ao_ranges", "coords_ang", "n_ao", "shells", "syms"
    ))


def _read_molecular_orbitals(args, *, n_ao, shells, tracker):
    """AO overlap and molecular orbitals (restricted or unrestricted), energy axis."""
    tracker.start_stage("AO Overlap Matrix (S)")
    logger.info("\n--- Computing AO overlap ---")
    t0_s = time.time()
    if getattr(args, "periodic_enabled", False):
        lattice_ang = np.asarray(args.lattice_vectors, dtype=float)
        lattice_bohr = lattice_ang * BOHR_PER_ANG
        cutoff = float(getattr(args, "overlap_cutoff", -1.0))
        if cutoff <= 0.0:
            logger.warning("  [Warning] periodic.overlap_cutoff <= 0 uses minimum-image overlap only.")
            logger.warning("  [Warning] Gamma-periodic MO orthonormality usually requires summing neighboring cell images; try overlap_cutoff: 6.0.")
        S = libint_cpp.overlap_pbc(shells, lattice_bohr, cutoff, args.nthreads)
        logger.debug(f"  ->  PBC overlap computed in {time.time() - t0_s:.2f} s")
        logger.info(f"  ->  Lattice vectors read from YAML in angstrom; Libint lattice passed in bohr. cutoff={cutoff:.3f} A")
    else:
        S = libint_cpp.overlap(shells, args.nthreads)
        logger.debug(f"  ->  Finite-system overlap computed in {time.time() - t0_s:.2f} s")

    C_beta, eps_beta, occ_beta = None, None, None
    homo_index_beta = None

    tracker.start_stage("Molecular Orbitals (MOs)")
    if args.mo_file_beta is not None:
        logger.info(f"\n--- Reading Alpha Molecular Orbitals from {args.mo_file} ---")
        t0_mos = time.time()
        C, eps, occ = read_mos_auto(args.mo_file, n_ao, verbose=True, cache=args.cache_mos)
        logger.debug(f"  -> Alpha MOs parsed in {time.time() - t0_mos:.2f} s | C shape {C.shape}")
        
        logger.info(f"\n--- Reading Beta Molecular Orbitals from {args.mo_file_beta} ---")
        t0_beta = time.time()
        C_beta, eps_beta, occ_beta = read_mos_auto(args.mo_file_beta, n_ao, verbose=True, cache=args.cache_mos)
        logger.debug(f"  -> Beta MOs parsed in {time.time() - t0_beta:.2f} s | C shape {C_beta.shape}")
    else:
        logger.info(f"\n--- Reading Molecular Orbitals from {args.mo_file} ---")
        t0_mos = time.time()
        C, eps, occ = read_mos_auto(args.mo_file, n_ao, verbose=True, cache=args.cache_mos)
        logger.debug(f"  -> MOs parsed in {time.time() - t0_mos:.2f} s | C shape {C.shape}")

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
        logger.info(
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
    return _export(locals(), (
        "C", "C_beta", "S", "e_fermi_raw", "e_homo", "e_lumo", "eps", "eps_beta", "eps_beta_shifted",
        "eps_shifted", "homo_index", "homo_index_beta", "occ", "occ_beta", "t0_gap"
    ))


def _quasiparticle_correction(args, *,
        C, S, atom_ao_ranges, coords_ang, eps, eps_beta, homo_index, homo_index_beta, shells, syms, tracker):
    """QP model: gap correction and the screened interaction W it is built on."""
    tracker.start_stage("Quasiparticle (GW) Model")
    dft_gap = eps_beta[homo_index_beta + 1] - eps[homo_index]
    target_qp_gap = dft_gap
    confinement_energy = 0.0
    qp_provenance = None
    eps_qp_active = None
    _xs_cache = {}
    # QP populations follow the BSE charge partition: Löwdin for charge_type lowdin and for xs
    # (the xs Hamiltonian is Löwdin-based), Mulliken otherwise (no diagonalization of S needed).
    qp_pop_mode = ("lowdin" if (str(getattr(args, "charge_type", "mulliken")).lower() == "lowdin"
                                or str(args.kernel_type).lower() in ("xs", "xs-qdex")) else "mulliken")
    _pop_cache = {}
    _sc_cache = {}

    def _get_SC():
        """S @ C for the current C, computed once and shared with the population analysis stage."""
        key = id(C)
        if key not in _sc_cache:
            _sc_cache.clear()
            C_d = C.toarray() if hasattr(C, "toarray") else C
            _sc_cache[key] = S @ C_d
        return _sc_cache[key]

    def _all_populations(representation):
        """Populations of ALL molecular orbitals (one vectorized pass, cached per representation)."""
        if representation not in _pop_cache:
            t0p = time.time()
            _pop_cache[representation] = orbital_populations(
                C, S, atom_ao_ranges, qp_pop_mode, representation,
                SC=_get_SC() if qp_pop_mode == "mulliken" else None)
            logger.debug(f"  [QP] {qp_pop_mode.capitalize()} populations of all {C.shape[1]} orbitals ({representation}) "
                  f"in {time.time() - t0p:.1f} s")
        return _pop_cache[representation]

    def _frontier_populations(representation="atom"):
        """Populations of HOMO and LUMO only (fast: avoids calculating all columns)."""
        key = ("frontier", representation)
        if key not in _pop_cache:
            t0p = time.time()
            h_cols = np.array([homo_index, homo_index + 1])
            C_front = C[:, h_cols]
            SC_front = None
            if qp_pop_mode == "mulliken":
                if id(C) in _sc_cache:
                    SC_front = _sc_cache[id(C)][:, h_cols]
                else:
                    C_front_d = C_front.toarray() if hasattr(C_front, "toarray") else C_front
                    SC_front = S @ C_front_d
            _pop_cache[key] = orbital_populations(
                C_front, S, atom_ao_ranges, qp_pop_mode, representation,
                SC=SC_front
            )
            logger.debug(f"  [QP] Frontier {qp_pop_mode} populations (HOMO/LUMO) in {time.time() - t0p:.3f} s")
        return _pop_cache[key]

    def _xs_gamma_ao():
        """Exact AO density-pair integrals (mu mu|nu nu) in eV when two_electron_integrals: xs."""
        if str(args.kernel_type).lower() not in ("xs", "xs-qdex"):
            return None
        if "gamma" not in _xs_cache:
            from qdex.integrals import compute_two_electron_ao
            logger.info("  [xs] Computing exact AO density-pair integrals for the shared W ...")
            _xs_cache["gamma"] = compute_two_electron_ao(shells, nthreads=args.nthreads) * HA_TO_EV
        return _xs_cache["gamma"]
    
    entry = MATERIAL_DB.get(args.material.upper(), None)
    if isinstance(args.qp_gap, str):
        if args.qp_gap.lower() == "brus":
            target_qp_gap = estimate_brus_qp_gap(material_name=args.material, coords=np.array(coords_ang), atom_symbols=syms)
            if target_qp_gap is not None:
                scissor = target_qp_gap - dft_gap
                confinement_energy = target_qp_gap - MATERIAL_DB.get(args.material.upper(), [0]*9)[3]
            else:
                raise ValueError(
                    "Brus QP model requested, but the required material data are missing. "
                    "Select qp_gap: pbe explicitly for an uncorrected calculation."
                )
                
        elif args.qp_gap.lower() in ["gw", "sgw-anchor", "sgw_anchor"]:
            if getattr(args, "periodic_enabled", False):
                logger.info("\n  [Bulk GW Model] Periodic mode enabled: using tabulated bulk GW-PBE scissor.")
                gw_scissor, qp_provenance = estimate_periodic_bulk_gw_scissor(args.material)
            else:
                # Compute the Scaled GW Scissor directly
                gw_scissor, qp_provenance = estimate_gw_qp_gap(
                    np.array(coords_ang), syms, args.material, args.eps_out,
                    regularization_length_ang=args.qp_regularization_length,
                    residual_power=args.qp_residual_power,
                    strict=args.qp_strict,
                    polarization_model=getattr(args, "qp_polarization", "sphere"),
                    return_details=True,
                    dft_gap=dft_gap,
                    residual_scaling=getattr(args, "qp_residual_scaling", "econf"),
                    anchor_residual=getattr(args, "qp_anchor_residual", None) or "on",
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
                logger.warning("  [QP Warning] qp_gap is set to 'sgw-anchor', which already uses the recommended scaled-GW hardness model.")
                logger.warning("  [QP Warning] estimate_qp enables the experimental COHSEX/TB Mulliken correction and can double-count QP shifts.")
                logger.warning("  [QP Warning] Production runs should use estimate_qp: false unless you explicitly want this experimental path.")
        elif args.qp_gap.lower() in ["sgw-dim", "sgw_dim", "sgw-dim-dz", "evgw", "evgw-dim", "evgw_dim", "qsgw", "qsgw-dim", "qsgw_dim", "scgw", "scgw-dim", "scgw_dim"]:
            from qdex.hardness import estimate_sgw_dim_qp_gap, estimate_qsgw_dim_qp_gap
            from qdex.lowdin import lowdin_sqrt

            is_qsgw = ("qsgw" in args.qp_gap.lower()) or ("scgw" in args.qp_gap.lower()) or getattr(args, "update_orbitals", False)
            is_evgw = "evgw" in args.qp_gap.lower()
            use_dz, z_val = resolve_qp_z(args, True)
            model_tag = "qsGW-DIM" if is_qsgw else ("evGW-DIM" if is_evgw else ("sGW-DIM (Dynamic Z)" if use_dz else "sGW-DIM"))
            logger.info(f"\n--- Estimating Quasiparticle Gap using {model_tag} (Atomistic Delta-W) ---")
            t0_sgw = time.time()

            if is_qsgw:
                sgw_scissor, qp_provenance, C_qp, eps_qp = estimate_qsgw_dim_qp_gap(
                    coords=np.array(coords_ang),
                    atom_symbols=syms,
                    C=C,
                    eps=eps,
                    S=S,
                    atom_ao_ranges=atom_ao_ranges,
                    homo_index=homo_index,
                    material_name=args.material,
                    eps_out=args.eps_out,
                    alpha=args.alpha,
                    dynamic_z=use_dz,
                    Z=z_val,
                    return_details=True,
                    gamma_ao=_xs_gamma_ao(),
                    solvent_term=getattr(args, "qp_solvent_term", "sphere"),
                )
                C = C_qp
                eps_qp_active = eps_qp
                scissor = sgw_scissor
                target_qp_gap = dft_gap + scissor
            else:
                occ_idx_a = np.array([homo_index])
                virt_idx_a = np.array([homo_index + 1])
                q_edge = _frontier_populations("atom")
                C_occ_low = C_virt_low = None

                sgw_scissor, qp_provenance = estimate_sgw_dim_qp_gap(
                    coords=np.array(coords_ang),
                    atom_symbols=syms,
                    material_name=args.material,
                    eps_out=args.eps_out,
                    C_occ_low=C_occ_low,
                    C_virt_low=C_virt_low,
                    eps_occ=eps[occ_idx_a],
                    eps_virt=eps[virt_idx_a],
                    atom_ao_ranges=atom_ao_ranges,
                    alpha=args.alpha,
                    dynamic_z=use_dz,
                    Z=z_val,
                    self_consistent=is_evgw,
                    return_details=True,
                    frontier_pops=(q_edge[:, 0], q_edge[:, 1]),
                    solvent_term=getattr(args, "qp_solvent_term", "sphere"),
                )
                scissor = sgw_scissor
                target_qp_gap = dft_gap + scissor

            confinement_energy = qp_provenance.get("confinement_shift_ev", 0.0)
            logger.debug(f"  -> {model_tag} Quasiparticle shift computed in {time.time() - t0_sgw:.2f} s")

            if args.estimate_qp:
                logger.warning(f"  [QP Warning] qp_gap is set to '{args.qp_gap}', which already applies the microscopic Delta-W quasiparticle correction.")
                logger.warning("  [QP Warning] estimate_qp enables the experimental COHSEX/TB Mulliken correction and can double-count QP shifts.")
                logger.warning("  [QP Warning] Production runs should use estimate_qp: false unless you explicitly want this experimental path.")

        elif args.qp_gap.lower() in ["sgw-resta", "sgw_resta", "sgw-resta-penn", "sgw-resta-pure", "sgw-resta-bulk", "sgw-resta-dz", "evgw-resta", "evgw_resta", "qsgw-resta", "qsgw_resta", "scgw-resta", "scgw_resta"]:
            from qdex.hardness import estimate_sgw_resta_qp_gap, estimate_qsgw_resta_qp_gap
            from qdex.lowdin import lowdin_sqrt

            is_qsgw = ("qsgw" in args.qp_gap.lower()) or ("scgw" in args.qp_gap.lower()) or getattr(args, "update_orbitals", False)
            is_evgw = "evgw" in args.qp_gap.lower()
            use_dz, z_val = resolve_qp_z(args, True)
            use_penn = not any(k in args.qp_gap.lower() for k in ["pure", "bulk"])
            penn_label = "Penn-scaled" if use_penn else "Pure Boundary"
            model_tag = f"qsGW-Resta ({penn_label})" if is_qsgw else (f"evGW-Resta ({penn_label})" if is_evgw else (f"sGW-Resta ({penn_label}, Dynamic Z)" if use_dz else f"sGW-Resta ({penn_label})"))
            logger.info(f"\n--- Estimating Quasiparticle Gap using {model_tag} (Delta-W) ---")
            t0_sgw = time.time()

            if is_qsgw:
                sgw_scissor, qp_provenance, C_qp, eps_qp = estimate_qsgw_resta_qp_gap(
                    coords=np.array(coords_ang),
                    atom_symbols=syms,
                    C=C,
                    eps=eps,
                    S=S,
                    atom_ao_ranges=atom_ao_ranges,
                    homo_index=homo_index,
                    material_name=args.material,
                    eps_out=args.eps_out,
                    alpha=args.alpha,
                    penn_scaling=use_penn,
                    dynamic_z=use_dz,
                    Z=z_val,
                    return_details=True,
                    gamma_ao=_xs_gamma_ao(),
                    solvent_term=getattr(args, "qp_solvent_term", "sphere"),
                )
                C = C_qp
                eps_qp_active = eps_qp
                scissor = sgw_scissor
                target_qp_gap = dft_gap + scissor
            else:
                occ_idx_a = np.array([homo_index])
                virt_idx_a = np.array([homo_index + 1])
                q_edge = _frontier_populations("atom")
                C_occ_low = C_virt_low = None

                sgw_scissor, qp_provenance = estimate_sgw_resta_qp_gap(
                    coords=np.array(coords_ang),
                    atom_symbols=syms,
                    material_name=args.material,
                    eps_out=args.eps_out,
                    dft_gap=dft_gap,
                    C_occ_low=C_occ_low,
                    C_virt_low=C_virt_low,
                    eps_occ=eps[occ_idx_a],
                    eps_virt=eps[virt_idx_a],
                    atom_ao_ranges=atom_ao_ranges,
                    alpha=args.alpha,
                    penn_scaling=use_penn,
                    dynamic_z=use_dz,
                    Z=z_val,
                    self_consistent=is_evgw,
                    return_details=True,
                    frontier_pops=(q_edge[:, 0], q_edge[:, 1]),
                    solvent_term=getattr(args, "qp_solvent_term", "sphere"),
                )
                scissor = sgw_scissor
                target_qp_gap = dft_gap + scissor

            confinement_energy = qp_provenance.get("confinement_shift_ev", 0.0)
            logger.debug(f"  -> {model_tag} Quasiparticle shift computed in {time.time() - t0_sgw:.2f} s")

            if args.estimate_qp:
                logger.warning(f"  [QP Warning] qp_gap is set to '{args.qp_gap}', which already applies the microscopic Delta-W quasiparticle correction.")
                logger.warning("  [QP Warning] estimate_qp enables the experimental COHSEX/TB Mulliken correction and can double-count QP shifts.")
                logger.warning("  [QP Warning] Production runs should use estimate_qp: false unless you explicitly want this experimental path.")

        elif args.qp_gap.lower() in ["sgw", "sgw-dw", "sgw_dw", "sgw-atom", "sgw-ao"]:
            from qdex.hardness import estimate_sgw_qp_gap
            from qdex.lowdin import lowdin_sqrt

            logger.info(f"\n--- Estimating Quasiparticle Gap using Microscopic sGW (Delta-W formulation) ---")
            t0_sgw = time.time()
            from qdex.lowdin import lowdin_apply
            n_occ_tot = homo_index + 1
            n_virt_tot = min(1000, len(eps) - n_occ_tot)
            occ_idx_a = np.arange(0, n_occ_tot)
            virt_idx_a = np.arange(n_occ_tot, n_occ_tot + n_virt_tot)
            C_occ_act = C[:, occ_idx_a].toarray() if hasattr(C, 'toarray') else C[:, occ_idx_a]
            C_virt_act = C[:, virt_idx_a].toarray() if hasattr(C, 'toarray') else C[:, virt_idx_a]
            C_occ_low = lowdin_apply(S, C_occ_act)   # sBSE needs all occupied (RPA response)
            C_virt_low = lowdin_apply(S, C_virt_act)
            eps_occ_act = eps[occ_idx_a]
            eps_virt_act = eps[virt_idx_a]

            sgw_mode = "ao" if "ao" in args.qp_gap.lower() else "atom"
            sgw_dz, sgw_z = resolve_qp_z(args, True)
            sgw_scissor, qp_provenance = estimate_sgw_qp_gap(
                coords=np.array(coords_ang),
                atom_symbols=syms,
                material_name=args.material,
                eps_out=args.eps_out,
                C_occ_low=C_occ_low,
                C_virt_low=C_virt_low,
                eps_occ=eps_occ_act,
                eps_virt=eps_virt_act,
                atom_ao_ranges=atom_ao_ranges,
                shells=shells,
                mode=sgw_mode,
                alpha=args.alpha,
                Z=sgw_z,
                dynamic_z=sgw_dz,
                nthreads=args.nthreads,
                return_details=True
            )
            scissor = sgw_scissor
            target_qp_gap = dft_gap + scissor
            confinement_energy = qp_provenance.get("confinement_shift_ev", 0.0)
            logger.debug(f"  -> sGW Quasiparticle shift computed in {time.time() - t0_sgw:.2f} s")

            if args.estimate_qp:
                logger.warning("  [QP Warning] qp_gap is set to 'sgw', which already applies the microscopic Delta-W quasiparticle correction.")
                logger.warning("  [QP Warning] estimate_qp enables the experimental COHSEX/TB Mulliken correction and can double-count QP shifts.")
                logger.warning("  [QP Warning] Production runs should use estimate_qp: false unless you explicitly want this experimental path.")
        elif args.qp_gap.lower() in ["none", "pbe", "dft"]:
            scissor = 0.0
            target_qp_gap = dft_gap
            logger.info("  [QP] No QP correction: DFT orbital energies used as they are (scissor = 0.0000 eV).")
            if str(getattr(args, "excitation_mode", "")).lower() in ("sbse", "diagonal_sbse"):
                logger.info("  [QP Notice] The sBSE adds the bulk GW correction to PBE orbitals: use quasiparticles.model: bulk. "
                      "'none' no longer adds it.")
        elif args.qp_gap.lower() == "bulk":
            if str(getattr(args, "qp_reference", "pbe")).lower() == "gxtb":
                raise ValueError("quasiparticles.reference: gxtb is only supported by the NAMD precompute (g-xTB frames).")
            logger.info("  [QP] Bulk GW correction (bulk QSGW - bulk PBE). Valid only for PBE orbitals.")
            from qdex.hardness import bulk_qp_shift
            scissor, bulk_vertex_info = bulk_qp_shift(args.material, dft_gap)
            if entry is not None and len(entry) >= 9:
                pbe_bulk_gap = float(entry[7])
                gw_bulk_gap = float(entry[8])
            target_qp_gap = dft_gap + scissor
            f_homo = anchor_bulk_homo_fraction(args.material)
            f_lumo = 1.0 - f_homo
            qp_provenance = {
                "qp_model": "bulk_gw_scissor",
                "material": str(args.material).upper(),
                "bulk_pbe_gap_ev": float(pbe_bulk_gap) if entry is not None and len(entry) >= 9 else 0.0,
                "bulk_gw_gap_ev": float(gw_bulk_gap) if entry is not None and len(entry) >= 9 else 0.0,
                "bulk_gw_shift_ev": float(scissor),
                "bulk_shift_ev": float(scissor),
                **bulk_vertex_info,
                "f_homo": float(f_homo),
                "f_lumo": float(f_lumo),
                "finite_size_shift_vacuum_ev": 0.0,
                "finite_size_shift_solvent_ev": 0.0,
                "anchor_residual_homo_ev": 0.0,
                "anchor_residual_lumo_ev": 0.0,
                "edge_residual_homo_ev": 0.0,
                "edge_residual_lumo_ev": 0.0,
                "residual_scale": 0.0,
                "anchor_residual_scale": 0.0,
                "bulk_homo_fraction": float(f_homo),
            }
            logger.info(f"  [QP] Pure Bulk GW Scissor mode selected ('{args.qp_gap}'):")
            logger.info(f"       -> Scissor = {scissor:+.4f} eV (PBE bulk {pbe_bulk_gap:.3f} -> GW bulk {gw_bulk_gap:.3f} eV)")
            logger.info(f"       -> Zero finite-size Delta-W or boundary polarization applied.")
        else:
            raise ValueError(
                f"Unknown qp_gap mode '{args.qp_gap}'. Use 'none' (uncorrected), 'bulk' (PBE + bulk GW), 'gw', 'sgw-dim', 'evgw-dim', 'qsgw-dim', 'sgw-resta', 'evgw-resta', 'qsgw-resta', 'sgw', 'brus', 'pbe', or a numeric gap."
            )
    else:
        # Numeric explicit gap provided
        target_qp_gap = float(args.qp_gap)
        scissor = target_qp_gap - dft_gap
    
    qp_w = None
    if qp_provenance is not None and "w_bse_ev" in qp_provenance:
        qp_w = (qp_provenance.pop("w_bse_ev"), qp_provenance.pop("w_bse_label", "QP-model W"))
    elif getattr(args, "qp_z", None) is not None or getattr(args, "dynamic_z", False):
        raise ValueError(
            f"qp_z / dynamic_z apply only to Delta-W QP models (sgw-*, evgw-*, qsgw-*, sgw); "
            f"qp_gap '{args.qp_gap}' has no Z."
        )
    if (qp_w is None and qp_provenance is not None
            and qp_provenance.get("polarization_model") == "sphere"
            and args.kernel in (None, "resta-sphere")):
        # Anchor-scaled gw model: the BSE sees the same dielectric sphere as the
        # QP polarization term, W = Resta(eps_inf) + sphere reaction field.
        from qdex.hardness import build_resta_mnok, build_sphere_reaction_field
        w_resta, _ = build_resta_mnok(syms, np.array(coords_ang), args.alpha, args.material, eps_out=args.eps_out)
        # The two-anchor gw model is defined with the hull radius (its anchor R0 uses it).
        w_sphere = build_sphere_reaction_field(np.array(coords_ang), syms, args.material, args.eps_out,
                                               radius_definition="hull")
        qp_w = (w_resta + w_sphere, f"Resta(eps_inf) + sphere reaction field (eps_out = {args.eps_out:.2f})")
        qp_provenance["bse_kernel_model"] = "resta_plus_sphere_reaction_field"
        qp_provenance["w_parts"] = {"w_qd": w_resta, "w_bulk": w_resta, "w_add": w_sphere,
                                    "gamma": build_gamma(syms, np.array(coords_ang), 1.0, 0.0),
                                    "eps_z": None, "bulk_shift": qp_provenance.get("bulk_gw_shift_ev", 0.0),
                                    "bulk_homo_fraction": anchor_bulk_homo_fraction(args.material),
                                    "anchor_edges": True}
        args.kernel = None
    elif args.kernel == "resta-sphere":
        raise ValueError("kernel 'resta-sphere' requires qp_gap: gw with qp_polarization: sphere.")
    return _export(locals(), (
        "C", "_all_populations", "_get_SC", "_xs_gamma_ao", "confinement_energy", "dft_gap", "entry",
        "eps_qp", "eps_qp_active", "f_homo", "f_lumo", "qp_pop_mode", "qp_provenance", "qp_w", "scissor",
        "target_qp_gap"
    ))


def _qp_levels_and_kernel(args, *,
        C, S, _all_populations, _xs_gamma_ao, atom_ao_ranges, coords_ang, dft_gap, eps, eps_qp,
        eps_qp_active, homo_index, qp_pop_mode, qp_provenance, qp_w, scissor, syms, target_qp_gap):
    """Orbital-resolved QP levels, the BSE kernel and the sTDA interactions; QP provenance."""
    w_parts = qp_provenance.pop("w_parts", None) if qp_provenance is not None else None
    use_xs = str(args.kernel_type).lower() in ("xs", "xs-qdex")
    if qp_w is not None and use_xs and w_parts is None:
        raise ValueError(f"qp_gap '{args.qp_gap}' provides W only in the atom (mnok) representation; "
                         "use two_electron_integrals: mnok with this model.")
    shared_gamma_ao = None
    if w_parts is not None:
        z_eh = float(qp_provenance.get("z_factor", 1.0)) if not w_parts.get("anchor_edges") else 1.0
        if use_xs:
            gamma_ao = _xs_gamma_ao()
            W_kernel_ao, _, dW_levels = xs_shared_w(w_parts, z_eh, gamma_ao, atom_ao_ranges)
            qp_w = (W_kernel_ao, qp_w[1] + " [xs: same screening on exact AO integrals]")
            shared_gamma_ao = gamma_ao
            representation = "ao"
        else:
            dW_levels = atom_delta_w(w_parts)
            representation = "atom"
        qp_provenance["two_electron_representation"] = "xs" if use_xs else "mnok"

        levels_mode = str(getattr(args, "qp_levels", "orbital")).lower()
        if levels_mode == "orbital" and eps_qp_active is None:
            qp_win_mode = str(getattr(args, "qp_window", "active")).lower()
            if qp_win_mode == "all":
                occ_w = np.arange(0, homo_index + 1)
                virt_w = np.arange(homo_index + 1, len(eps))
                eval_indices = None
            else:
                n_win_override = getattr(args, "qp_window_size", None)
                if n_win_override is not None:
                    n_occ_req = int(n_win_override)
                    n_virt_req = int(n_win_override)
                elif args.nhomos is not None or args.nlumos is not None:
                    n_occ_req = args.nhomos if args.nhomos is not None else 25
                    n_virt_req = args.nlumos if args.nlumos is not None else 25
                elif args.e_thresh is not None:
                    dft_pairs = (eps[homo_index + 1:] - eps[homo_index]) <= float(args.e_thresh)
                    n_virt_req = int(np.where(dft_pairs)[0][-1] + 1) if np.any(dft_pairs) else 25
                    dft_pairs_occ = (eps[homo_index + 1] - eps[:homo_index + 1]) <= float(args.e_thresh)
                    n_occ_req = int(homo_index - np.where(dft_pairs_occ)[0][0] + 1) if np.any(dft_pairs_occ) else 25
                else:
                    n_occ_req = 25
                    n_virt_req = 25
                n_occ_act = min(homo_index + 1, max(1, n_occ_req))
                n_virt_act = min(len(eps) - (homo_index + 1), max(1, n_virt_req))
                occ_w = np.arange(homo_index - n_occ_act + 1, homo_index + 1)
                virt_w = np.arange(homo_index + 1, homo_index + 1 + n_virt_act)
                eval_indices = np.concatenate([occ_w, virt_w])

            is_anchor_model = bool(w_parts.get("anchor_edges"))
            z_mode = "derived" if qp_provenance.get("dynamic_z") else "fixed"
            z_fixed = float(qp_provenance.get("z_factor", 1.0))
            fb = float(w_parts.get("bulk_homo_fraction", 0.5))
            selfenergy = "classical" if is_anchor_model else str(getattr(args, "qp_selfenergy", "cohsex")).lower()
            if selfenergy == "cohsex":
                t0c = time.time()
                coh, sex = cohsex_diagonal(C, S, homo_index, atom_ao_ranges,
                                           dW_atom=None if use_xs else dW_levels,
                                           dW_ao=dW_levels if use_xs else None,
                                           eval_indices=eval_indices)
                eps_qp, lev_info = cohsex_qp_energies(
                    eps, coh, sex, homo_index, float(w_parts.get("bulk_shift", 0.0)), z_mode, z_fixed,
                    w_parts.get("eps_z"), args.material, homo_fraction_bulk=fb,
                    occ_idx=occ_w, virt_idx=virt_w)
                n_eval_str = f"{len(eval_indices)} active" if eval_indices is not None else f"all {len(eps)}"
                logger.info(f"\n  [QP Levels] One-shot Delta-COHSEX for {n_eval_str} orbitals in {time.time() - t0c:.1f} s: "
                      f"HOMO COH {coh[homo_index]:+.3f} SEX {sex[homo_index]:+.3f} eV; "
                      f"LUMO COH {coh[homo_index + 1]:+.3f} SEX {sex[homo_index + 1]:+.3f} eV")
            else:
                if eval_indices is not None:
                    C_act = C[:, eval_indices]
                    q_act = orbital_populations(C_act, S, atom_ao_ranges, qp_pop_mode, representation)
                    q_occ = q_act[:, :len(occ_w)]
                    q_virt = q_act[:, len(occ_w):]
                else:
                    q_w = _all_populations(representation)
                    q_occ = q_w[:, :len(occ_w)]
                    q_virt = q_w[:, len(occ_w):]
                edges = None
                if is_anchor_model:
                    f_h = float(qp_provenance.get("f_homo", 0.5))
                    edges = (f_h * scissor, (1.0 - f_h) * scissor)
                eps_qp, lev_info = orbital_qp_energies(
                    eps, q_occ, q_virt, occ_w, virt_w, dW_levels,
                    float(w_parts.get("bulk_shift", 0.0)), z_mode, z_fixed,
                    w_parts.get("eps_z"), args.material, homo_fraction_bulk=fb, edge_shifts=edges,
                    representation=representation)
                lev_info["qp_selfenergy"] = "classical"
            lev_info["qp_populations"] = qp_pop_mode

            # The BSE kernel carries the same Z as the QP levels: W_BSE = W_bulk + Zbar (W_QD + W_add - W_bulk),
            # Zbar = (Z_HOMO + Z_LUMO)/2 of the Delta-COHSEX levels (the estimator's classical Z is replaced).
            if selfenergy == "cohsex" and not is_anchor_model:
                z_new = 0.5 * (float(lev_info["z_homo_orbital"]) + float(lev_info["z_lumo_orbital"]))
                if abs(z_new - z_eh) > 1.0e-12:
                    from qdex.hardness import scale_w_difference
                    if use_xs:
                        W_new = xs_shared_w(w_parts, z_new, _xs_gamma_ao(), atom_ao_ranges)[0]
                    else:
                        add = w_parts.get("w_add")
                        W_new = scale_w_difference(w_parts["w_qd"] + (0.0 if add is None else add),
                                                   w_parts["w_bulk"], z_new)
                    qp_w = (W_new, qp_w[1])
                    logger.info(f"  [Consistency] BSE kernel Zbar = {z_new:.4f} (Delta-COHSEX levels; estimator gave {z_eh:.4f}).")
                lev_info["w_bse_z"] = float(z_new)

            # Non-classical anchor residual of the Delta-W models (calibrated on the evGW anchor)
            if not is_anchor_model:
                from qdex.hardness import anchor_residual_scale
                z_label = "derived" if z_mode == "derived" else f"{z_fixed:g}"
                a_key = anchor_key(args.material, str(args.qp_gap).lower(), selfenergy,
                                   "xs" if use_xs else "mnok", qp_pop_mode, z_label,
                                   getattr(args, "qp_solvent_term", "sphere"))
                table_path = getattr(args, "qp_anchor_table", None)
                ent = MATERIAL_DB.get(str(args.material).upper(), ())
                if getattr(args, "qp_anchor_calibrate", False):
                    if len(ent) < 14:
                        raise ValueError("qp_anchor_calibrate needs monomer evGW data in MATERIAL_DB.")
                    if abs(float(args.eps_out) - 1.0) > 1e-9:
                        logger.warning("  [Anchor] Warning: calibrating with eps_out != 1; the evGW anchor is in vacuum.")
                    d_h = float(eps_qp[homo_index] - eps[homo_index])
                    d_l = float(eps_qp[homo_index + 1] - eps[homo_index + 1])
                    res_h = (float(ent[12]) - float(ent[10])) - d_h
                    res_l = (float(ent[13]) - float(ent[11])) - d_l
                    path_w = save_anchor_entry(a_key, {"residual_homo_ev": res_h, "residual_lumo_ev": res_l,
                                                       "anchor_dft_gap_ev": float(dft_gap),
                                                       "model_homo_shift_ev": d_h, "model_lumo_shift_ev": d_l},
                                               table_path)
                    logger.info(f"  [Anchor] Calibrated '{a_key}': residual HOMO {res_h:+.3f} eV, LUMO {res_l:+.3f} eV "
                          f"-> {path_w}")
                    lev_info.update({"anchor_residual_homo_ev": res_h, "anchor_residual_lumo_ev": res_l,
                                     "anchor_residual_scale": 1.0, "anchor_key": a_key})
                    eps_qp[:homo_index + 1] += res_h
                    eps_qp[homo_index + 1:] += res_l
                elif str(getattr(args, "qp_anchor_residual", None) or "off").lower() == "on" and len(ent) >= 14:
                    table = load_anchor_table(table_path)
                    if a_key in table:
                        r_cl = qp_provenance.get("cluster_radius_ang") or get_cluster_size_metrics(
                            np.array(coords_ang), syms, args.material)["R_eff_hull"]
                        scale, smode = anchor_residual_scale(args.material, float(r_cl), dft_gap,
                                                             getattr(args, "qp_residual_scaling", "econf"),
                                                             args.qp_residual_power)
                        res_h = float(table[a_key]["residual_homo_ev"]) * scale
                        res_l = float(table[a_key]["residual_lumo_ev"]) * scale
                        eps_qp[:homo_index + 1] += res_h
                        eps_qp[homo_index + 1:] += res_l
                        logger.info(f"  [Anchor] Residual '{a_key}': HOMO {res_h:+.3f}, LUMO {res_l:+.3f} eV "
                              f"(anchor values x {scale:.3f}, scaling {smode})")
                        lev_info.update({"anchor_residual_homo_ev": res_h, "anchor_residual_lumo_ev": res_l,
                                         "anchor_residual_scale": scale, "anchor_key": a_key})
                    else:
                        logger.info(f"  [Anchor] No calibrated residual for '{a_key}' (run the anchor cluster with "
                              f"--qp-anchor-calibrate); none applied.")
                        lev_info["anchor_key"] = a_key
            new_scissor = float(eps_qp[homo_index + 1] - eps_qp[homo_index]) - dft_gap
            logger.info(f"\n  [QP Levels] Orbital-resolved ({representation}): {len(occ_w)} occ + {len(virt_w)} virt levels; "
                  f"HOMO {lev_info['qp_homo_shift_ev']:+.3f} eV, LUMO {lev_info['qp_lumo_shift_ev']:+.3f} eV; "
                  f"spread occ {lev_info['qp_shift_spread_occ_ev']:.3f} / virt {lev_info['qp_shift_spread_virt_ev']:.3f} eV; "
                  f"Z in [{lev_info['z_min_window']:.3f}, {lev_info['z_max_window']:.3f}]")
            if abs(new_scissor - scissor) > 1e-4:
                logger.info(f"  [QP Levels] Gap correction {scissor:+.4f} -> {new_scissor:+.4f} eV "
                      f"({'xs integrals' if use_xs else 'same formula, all orbitals'}).")
            lev_info["model_scissor_ev"] = float(scissor)
            qp_provenance.update(lev_info)
            scissor = new_scissor
            target_qp_gap = dft_gap + scissor
            eps_qp_active = eps_qp
        else:
            qp_provenance["qp_levels"] = "orbital (qsGW)" if eps_qp_active is not None else "rigid"
    args.kernel = resolve_bse_kernel(args, qp_w)
    stda_gamma_j = stda_gamma_k = None
    if args.kernel == "stda":
        from qdex.hardness import build_stda_gammas, stda_parameters
        ax_val, alpha_k, beta_j, ax_src = stda_parameters(
            getattr(args, "stda_functional", None), getattr(args, "stda_ax", None),
            getattr(args, "stda_alpha", None), getattr(args, "stda_beta", None), args.material)
        stda_gamma_j, stda_gamma_k, stda_info = build_stda_gammas(syms, np.array(coords_ang), ax_val,
                                                                  alpha=alpha_k, beta=beta_j)
        if str(args.charge_type).lower() != "lowdin":
            logger.info("  [sTDA] Transition charges set to Loewdin, as in sTDA.")
            args.charge_type = "lowdin"
        logger.info(f"  [sTDA] a_x = {ax_val:.3f} ({ax_src}); beta(J) = {stda_info['beta_J']:.3f}, "
              f"alpha(K) = {stda_info['alpha_K']:.3f}")
        qp_name = str(args.qp_gap).lower()
        functional = str(getattr(args, "stda_functional", "") or "").lower()
        if ax_val == 0.0:
            logger.info("  [sTDA Notice] a_x = 0: no electron-hole attraction (gamma^J = 0); only K^x acts on the DFT gap.")
        if qp_name == "bulk" and functional not in ("", "pbe"):
            logger.warning("  [sTDA Warning] quasiparticles.model: bulk corrects PBE orbitals only; the MO file is "
                  f"declared as '{functional}'.")
        if qp_name not in ("none", "pbe", "dft", "bulk"):
            logger.warning(f"  [sTDA Warning] qp_gap '{args.qp_gap}' changes the orbital energies; faithful sTDA uses 'none'.")
        if qp_name == "bulk":
            logger.info("  [sTDA] PBE orbitals + bulk GW correction: dielectric variant, not standard sTDA.")

    if str(getattr(args, "excitation_mode", "")).lower() not in ("independent_qp", "independent_dft"):
        from qdex.hardness import format_integrals_block
        logger.info(format_integrals_block(
            args.kernel_type, args.charge_type, args.kernel, syms,
            stda_info=stda_info if args.kernel == "stda" else None,
            stda_ax_source=ax_src if args.kernel == "stda" else None,
            include_direct=getattr(args, "include_direct_eh", True) is not False,
            include_exchange=bool(getattr(args, "include_exchange", True)) and not getattr(args, "triplet", False),
            hubbard_beta=float(getattr(args, "beta", 0.0) or 0.0)))

    logger.info(f"\n  [DFT] Initial Gap  : {dft_gap:.4f} eV")
    logger.info(f"  [QP]  Target Gap   : {target_qp_gap:.4f} eV")
    logger.info(f"  [QP]  Scissor Shift: {scissor:.4f} eV")
    qp_provenance_file = write_qp_provenance(qp_provenance, dft_gap, target_qp_gap, scissor, args)
    print_qp_provenance(qp_provenance, dft_gap=dft_gap, target_qp_gap=target_qp_gap, output_file=qp_provenance_file)
    return _export(locals(), (
        "C_act", "anchor_residual_scale", "edges", "eps_qp_active", "f_h", "fb", "qp_w", "scissor",
        "shared_gamma_ao", "stda_gamma_j", "stda_gamma_k", "target_qp_gap"
    ))


def _ip_ea_and_energy_axis(args, *,
        C_beta, anchor_residual_scale, coords_ang, dft_gap, e_fermi_raw, edges, entry, eps, eps_qp_active,
        eps_shifted, f_homo, f_lumo, homo_index, qp_provenance, scissor, syms, t0_gap, target_qp_gap):
    """IP/EA prediction and the shifted energy axis."""
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

        logger.info(f"\n  [Periodic Band Edges]")
        logger.info(f"    Raw CP2K HOMO    : {dft_homo_raw:8.4f} eV (arbitrary periodic eigenvalue zero)")
        logger.info(f"    Raw CP2K LUMO    : {dft_lumo_raw:8.4f} eV")
        logger.info(f"    Bulk GW scissor  : virtual manifold shifted by {scissor:+.4f} eV")
        logger.info("    Note             : absolute IP/EA levels are not assigned in periodic mode.")

    elif qp_provenance is not None and "f_homo" in qp_provenance and "f_lumo" in qp_provenance:
        if (getattr(args, "qp_edge_split", "model") == "anchor"
                and qp_provenance.get("edge_split_source") is None
                and entry is not None and len(entry) >= 14):
            # Delta-W models split the correction almost 50/50 (classical charging);
            # take the HOMO/LUMO split from the per-edge two-anchor curves instead.
            from qdex.hardness import anchor_edge_curves, get_cluster_size_metrics
            r_split = qp_provenance.get("cluster_radius_ang")
            if r_split is None:
                r_split = get_cluster_size_metrics(np.array(coords_ang), syms, args.material)["R_eff_hull"]
            from qdex.hardness import anchor_residual_scale
            dec, _ = anchor_residual_scale(args.material, float(r_split), dft_gap,
                                           getattr(args, "qp_residual_scaling", "econf"), args.qp_residual_power)
            edges = anchor_edge_curves(args.material, float(r_split), args.eps_out,
                                       residual_power=args.qp_residual_power, decay=dec)
            if edges is not None:
                qp_provenance.setdefault("f_homo_micro", qp_provenance["f_homo"])
                qp_provenance.setdefault("f_lumo_micro", qp_provenance["f_lumo"])
                qp_provenance["f_homo"] = float(edges["f_homo"])
                qp_provenance["f_lumo"] = float(edges["f_lumo"])
                qp_provenance["edge_split_source"] = "anchor_edge_curves"
                logger.info(f"\n  [QP Edge Split] Anchor-calibrated: HOMO {edges['f_homo']*100:.1f}% / LUMO "
                      f"{edges['f_lumo']*100:.1f}% (model's own split: {qp_provenance['f_homo_micro']*100:.1f}% / "
                      f"{qp_provenance['f_lumo_micro']*100:.1f}%)")
        f_homo = float(qp_provenance["f_homo"])
        f_lumo = float(qp_provenance["f_lumo"])
        model_name = qp_provenance.get("qp_model", "microscopic").upper()

        if entry is not None and len(entry) >= 14:
            pbe_h_mono, pbe_l_mono = entry[10], entry[11]
            gap_pbe_mono = pbe_l_mono - pbe_h_mono
            shrinkage_pbe = gap_pbe_mono - dft_gap
            true_pbe_homo = pbe_h_mono + (shrinkage_pbe * f_homo)
            true_pbe_lumo = pbe_l_mono - (shrinkage_pbe * f_lumo)
            qp_homo = true_pbe_homo - (scissor * f_homo)
            qp_lumo = true_pbe_lumo + (scissor * f_lumo)
            logger.info(f"\n  [Absolute Band Edges (IP & EA - Microscopic Wavefunction Asymmetry)]")
            logger.info(f"    Raw CP2K HOMO    : {dft_homo_raw:8.4f} eV (Floating Vacuum)")
            logger.info(f"    Modeled PBE HOMO : {true_pbe_homo:8.4f} eV (vacuum-anchored)")
            logger.info(f"    -> Shift Split   : HOMO takes {f_homo*100:.1f}%, LUMO takes {f_lumo*100:.1f}% ({model_name})")
        else:
            qp_homo = dft_homo_raw - (scissor * f_homo)
            qp_lumo = dft_lumo_raw + (scissor * f_lumo)
            logger.info(f"\n  [Absolute Band Edges (IP & EA - Microscopic Wavefunction Asymmetry)]")
            logger.info(f"    Raw CP2K HOMO    : {dft_homo_raw:8.4f} eV")
            logger.info(f"    -> Shift Split   : HOMO takes {f_homo*100:.1f}%, LUMO takes {f_lumo*100:.1f}% ({model_name})")

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
            logger.warning("  [QP Warning] Anchor frontier shifts do not bracket the PBE gap; using a symmetric edge split.")
        
        # 2. Project Absolute PBE Levels (Bypassing CP2K floating vacuum)
        # We use the computed intermediate PBE gap (dft_gap) as the physical truth
        shrinkage_pbe = gap_pbe_mono - dft_gap 
        
        true_pbe_homo = pbe_h_mono + (shrinkage_pbe * f_homo)
        true_pbe_lumo = pbe_l_mono - (shrinkage_pbe * f_lumo)
        
        # 3. Apply the Dielectric Scissor to get QP levels
        qp_homo = true_pbe_homo - (scissor * f_homo)
        qp_lumo = true_pbe_lumo + (scissor * f_lumo)
        
        logger.info(f"\n  [Absolute Band Edges (IP & EA)]")
        logger.info(f"    Raw CP2K HOMO    : {dft_homo_raw:8.4f} eV (Floating Vacuum)")
        logger.info(f"    Modeled PBE HOMO : {true_pbe_homo:8.4f} eV (anchor-reconstructed)")
        logger.info(f"    -> Shift Split   : HOMO takes {f_homo*100:.1f}%, LUMO takes {f_lumo*100:.1f}%")

    else:
        # Fallback if no 14-item monomer data is available
        f_homo, f_lumo = 0.5, 0.5
        qp_homo = dft_homo_raw - (scissor * f_homo)
        qp_lumo = dft_lumo_raw + (scissor * f_lumo)
        
        logger.info(f"\n  [Absolute Band Edges (IP & EA)]")
        logger.info(f"    Raw CP2K HOMO    : {dft_homo_raw:8.4f} eV")
        logger.info(f"    -> Shift Split   : HOMO takes 50.0%, LUMO takes 50.0% (Default)")

    if getattr(args, "periodic_enabled", False):
        logger.info(f"    QP HOMO-like     : {qp_homo:8.4f} eV (relative eigenvalue)")
        logger.info(f"    QP LUMO-like     : {qp_lumo:8.4f} eV (relative eigenvalue)")
    else:
        logger.info(f"    QP HOMO (IP)     : {qp_homo:8.4f} eV   -> IP = {-qp_homo:8.4f} eV")
        logger.info(f"    QP LUMO (EA)     : {qp_lumo:8.4f} eV   -> EA = {-qp_lumo:8.4f} eV")
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

    eps_dft_shifted = eps_shifted.copy()
    if eps_qp_active is not None:
        eps_shifted = eps_qp_active - e_fermi_raw
        if C_beta is None:
            eps_beta_shifted = eps_shifted

    logger.debug(f"  -> Energy axis shifted and target gap resolved in {time.time() - t0_gap:.4f} s")
 
    # -----------------------------------------------------------------
    # Unified S@C Computation
    # -----------------------------------------------------------------
    return _export(locals(), (
        "eps_beta_shifted", "eps_dft_shifted", "eps_shifted", "qp_homo", "qp_lumo"
    ))


def _orbital_populations(args, *, C, C_beta, S, _get_SC, tracker):
    """S @ C, orthonormality check and spin-free populations."""
    tracker.start_stage("MO Orthonormality & Populations")
    logger.info("\n--- Computing Unified S@C Population Analysis ---")
    t0_pop = time.time()
    C_dense = C.toarray() if hasattr(C, 'toarray') else C
    SC_dense = _get_SC()  # reused if the QP step already computed it
    C_dense_beta_pop, SC_dense_beta_pop, pops_beta = None, None, None
    
    if getattr(args, "skip_orthonormality_check", False):
        logger.debug("\n  [Diag] MO orthonormality check skipped (skip_orthonormality_check).")
        # SOC keeps its cheap re-orthonormalization of the active window (the safe path).
        soc_assume_orthonormal = False
        if C_beta is not None:
            C_dense_beta_pop = C_beta.toarray() if hasattr(C_beta, 'toarray') else np.asarray(C_beta)
            SC_dense_beta_pop = S @ C_dense_beta_pop
    else:
        # === DIAGNOSTIC: STRICT C^T S C ORTHONORMALITY CHECK ===
        overlap_label = "PBC" if getattr(args, "periodic_enabled", False) else "finite"
        logger.debug(f"\n  [Diag] Testing MO Orthonormality with {overlap_label} overlap (C^T S C = I) ...")
        norm_matrix = C_dense.conj().T @ SC_dense
        orth_delta = norm_matrix - np.eye(C_dense.shape[1])
        orth_err = np.linalg.norm(orth_delta)
        orth_max = np.max(np.abs(orth_delta))
        logger.info(f"[CHECK] alpha ||C†SC - I||_F = {orth_err:.3e}")
        logger.info(f"[CHECK] alpha max|C†SC - I|  = {orth_max:.3e}")
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
            logger.info(f"[CHECK] beta  ||C†SC - I||_F = {orth_err_beta:.3e}")
            logger.info(f"[CHECK] beta  max|C†SC - I|  = {orth_max_beta:.3e}")
            if orth_max_beta > args.orthonormality_tol:
                raise ValueError(
                    f"Beta MO orthonormality failure: max|C†SC-I|={orth_max_beta:.3e} exceeds "
                    f"{args.orthonormality_tol:.3e}."
                )
            soc_assume_orthonormal = soc_assume_orthonormal and orth_max_beta < 1.0e-6
            cross_err = np.linalg.norm(C_dense.conj().T @ SC_dense_beta_pop)
            logger.info(f"[CHECK] alpha/beta ||Caᵀ S Cb|| = {cross_err:.3e} (diagnostic)")
    # ===========================================================

    pops_sf = np.real(C_dense.conj() * SC_dense)
    if C_dense_beta_pop is not None:
        pops_beta = np.real(C_dense_beta_pop.conj() * SC_dense_beta_pop)
    logger.debug(f"  -> S@C projection and populations computed in {time.time() - t0_pop:.2f} s")
    return _export(locals(), (
        "C_dense", "SC_dense", "SC_dense_beta_pop", "pops_beta", "pops_sf", "soc_assume_orthonormal"
    ))


def _active_space_and_soc(args, *,
        C_act, C_beta, C_dense, S, SC_dense, SC_dense_beta_pop, coords_ang, dft_gap, e_fermi_raw, eps,
        eps_beta, eps_beta_shifted, eps_dft_shifted, eps_qp_active, eps_shifted, f_h, fb, homo_index,
        homo_index_beta, occ, occ_beta, pops_beta, pops_sf, qp_provenance, scissor, shells,
        soc_assume_orthonormal, syms, target_qp_gap, tracker):
    """BSE active space, SOC spinor subspace, population printouts and the scissor operator."""
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
        logger.info(f"\n--- Computing SOC Spinor Subspace for BSE (Small Window) ---")
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
                C_AO=C_dense, eps_Ha=(eps_qp_active if eps_qp_active is not None else eps) / HA_TO_EV, S_AO=S, 
                active_indices=bse_active_indices, gth_file=args.gth_file,
                nthreads=args.nthreads, assume_orthonormal=soc_assume_orthonormal,
                SC_AO=SC_dense, device=args.device,
            )
            bse_spinor_homo_idx = (bse_n_occ * 2) - 1
        bse_soc_E = (bse_soc_E * HA_TO_EV) - e_fermi_raw
        bse_soc_E -= (bse_soc_E[bse_spinor_homo_idx] + bse_soc_E[bse_spinor_homo_idx + 1]) / 2.0

    qp_breakdown_alpha = None
    if qp_provenance is not None and getattr(args, "qp_gap", "").lower() not in ("pbe", "none", "dft"):
        fb = qp_provenance.get("bulk_homo_fraction", anchor_bulk_homo_fraction(args.material))
        d_bulk = float(qp_provenance.get("bulk_gw_shift_ev", qp_provenance.get("bulk_shift_ev", 0.0)))
        bulk_arr = np.where(np.arange(len(eps)) <= homo_index, -fb * d_bulk, (1.0 - fb) * d_bulk)
        
        if eps_qp_active is not None:
            zn_arr = qp_provenance.get("z_factors", np.ones(len(eps)))
            d_sig_arr = qp_provenance.get("delta_sigmas", np.zeros(len(eps)))
            r_scale = float(qp_provenance.get("anchor_residual_scale", 1.0))
            if str(getattr(args, "qp_anchor_residual", None) or "off").lower() == "off":
                r_scale = 0.0
            r_h = float(qp_provenance.get("anchor_residual_homo_ev", 0.0)) / max(r_scale, 1e-12) if r_scale > 0 else 0.0
            r_l = float(qp_provenance.get("anchor_residual_lumo_ev", 0.0)) / max(r_scale, 1e-12) if r_scale > 0 else 0.0
            r_hl_arr = np.where(np.arange(len(eps)) <= homo_index, r_h, r_l)
            eps_qp_print = eps_shifted
        else:
            zn_arr = np.ones(len(eps))
            f_h = float(qp_provenance.get("f_homo", 0.5))
            f_l = float(qp_provenance.get("f_lumo", 0.5))
            pol_vac = float(qp_provenance.get("finite_size_shift_vacuum_ev", 0.0))
            pol_solv = float(qp_provenance.get("finite_size_shift_solvent_ev", pol_vac))
            pol = pol_solv if args.eps_out > 1.0 else pol_vac
            d_sig_arr = np.where(np.arange(len(eps)) <= homo_index, -f_h * pol, f_l * pol)
            r_scale = float(qp_provenance.get("residual_scale", 1.0))
            if str(getattr(args, "qp_anchor_residual", None) or "on").lower() == "off":
                r_scale = 0.0
            r_h = float(qp_provenance.get("edge_residual_homo_ev", 0.0))
            r_l = float(qp_provenance.get("edge_residual_lumo_ev", 0.0))
            r_hl_arr = np.where(np.arange(len(eps)) <= homo_index, r_h, r_l)
            tot_shift_arr = bulk_arr + zn_arr * d_sig_arr + r_hl_arr * r_scale
            eps_qp_print = eps_dft_shifted + tot_shift_arr

        qp_breakdown_alpha = {
            "eps_dft": eps_dft_shifted,
            "eps_qp": eps_qp_print,
            "bulk_shift": bulk_arr,
            "zn": zn_arr,
            "delta_sigma": d_sig_arr,
            "r_hl": r_hl_arr,
            "s_r": r_scale,
        }

    if C_beta is not None:
        pop_range = max(0, int(getattr(args, "population_print_range", 15)))
        pop_tags = getattr(args, "population_bars", None)
        logger.info("\n--- Spin-Free Alpha MO Population Analysis ---")
        print_orbital_summary(eps_shifted, occ, homo_index, pops_sf, syms, shells, is_soc=False, print_range=pop_range, population_bars=pop_tags, qp_breakdown=qp_breakdown_alpha)
        logger.info("\n--- Spin-Free Beta MO Population Analysis ---")
        print_orbital_summary(eps_beta_shifted, occ_beta, homo_index_beta, pops_beta, syms, shells, is_soc=False, print_range=pop_range, population_bars=pop_tags)
    else:
        logger.info("\n--- Spin-Free MO Population Analysis ---")
        pop_range = max(0, int(getattr(args, "population_print_range", 15)))
        pop_tags = getattr(args, "population_bars", None)
        print_orbital_summary(eps_shifted, occ, homo_index, pops_sf, syms, shells, is_soc=False, print_range=pop_range, population_bars=pop_tags, qp_breakdown=qp_breakdown_alpha)

    if args.soc_flag:
        if args.gth_file is None: sys.exit("ERROR: --gth_file is required when --soc_flag is enabled.")
        logger.info("\n--- SOC Spinor Population Analysis (BSE Active Space) ---")
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
                       
        logger.debug(f"  -> Projected populations in {time.time() - t_pop:.2f}s")
        
        soc_occ = np.zeros_like(bse_soc_E)
        soc_occ[:bse_spinor_homo_idx + 1] = 1.0
        soc_offset = (homo_index - bse_n_occ + 1) * 2 
        print_orbital_summary(
            bse_soc_E, soc_occ, bse_spinor_homo_idx, pops_soc_act, syms, shells,
            is_soc=True, offset=soc_offset, print_range=pop_range, population_bars=pop_tags
        )
        calculated_soc_gap = bse_soc_E[bse_spinor_homo_idx + 1] - bse_soc_E[bse_spinor_homo_idx]
        if eps_qp_active is not None:
            # Spinors were built from QP energies (qsgw-*), so their gap already
            # contains the QP opening.  Express it on the DFT reference so that
            # every downstream "SOC gap + scissor" is the SOC QP gap.
            calculated_soc_gap -= scissor

        # --- RIGID SCISSOR APPLICATION ---
        logger.info("\n--- Scissor Operator Application ---")
        logger.info(f"  Rigid Scissor (Computed above)  : {scissor:+8.4f} eV")
        logger.info(f"  Spin-Free DFT Gap               : {dft_gap:8.4f} eV")
        logger.info(f"  SOC-Shrunken DFT Gap            : {calculated_soc_gap:8.4f} eV")
        
        final_sf_qp_gap = dft_gap + scissor
        final_soc_qp_gap = calculated_soc_gap + scissor
        
        logger.info(f"  -> Final Spin-Free QP Gap       : {final_sf_qp_gap:8.4f} eV")
        logger.info(f"  -> Final SOC QP Gap             : {final_soc_qp_gap:8.4f} eV (D_SOC = {dft_gap - calculated_soc_gap:.4f} eV)")

    else:
        # --- RIGID SCISSOR FOR SPIN-FREE ONLY ---
        logger.info("\n--- Scissor Operator Application (Spin-Free) ---")
        logger.info(f"  Rigid Scissor (Computed above)  : {scissor:+8.4f} eV")
        logger.info(f"  Spin-Free DFT Gap               : {dft_gap:8.4f} eV")
        logger.info(f"  -> Final Spin-Free QP Gap       : {dft_gap + scissor:8.4f} eV")

    # Update Confinement Energy (Always relative to bulk)
    db_gap = MATERIAL_DB.get(args.material.upper(), MATERIAL_DB["DEFAULT"])[3]
    confinement_energy = target_qp_gap - db_gap
    return _export(locals(), (
        "C_dense_beta", "bse_n_occ", "bse_n_occ_beta", "bse_n_virt", "bse_n_virt_beta", "bse_soc_E",
        "bse_soc_U", "bse_spinor_homo_idx", "calculated_soc_gap", "compute_spinor_subspace",
        "compute_spinor_subspace_uks", "confinement_energy", "db_gap", "is_uks_sp", "soc_overlap_cache"
    ))


def _cubes_and_fuzzy(args, *,
        C_beta, C_dense, C_dense_beta, S, SC_dense, SC_dense_beta_pop, bse_n_occ, bse_soc_U,
        bse_spinor_homo_idx, compute_spinor_subspace, compute_spinor_subspace_uks, coords_ang, e_fermi_raw,
        e_homo, e_lumo, eps, eps_beta, eps_beta_shifted, eps_dft_shifted, eps_qp_active, eps_shifted,
        homo_index, homo_index_beta, is_uks_sp, occ, pops_sf, qp_homo, qp_lumo, scissor, shells,
        soc_assume_orthonormal, soc_overlap_cache, syms, tracker):
    """Optional cube files, fuzzy bands and PDOS/COOP."""
    # -----------------------------------------------------------------
    # EXCITON CUBE GENERATION: EXECUTED BEFORE FUZZY PLOTTING 
    # -----------------------------------------------------------------
    if getattr(args, 'cube', False):
        tracker.start_stage("Exciton Cube Generation")
        from qdex.exciton_cube import generate_cubes
        logger.info("\n--- Generating Cubes for MOs / Spinors ---")
        
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

        if eps_qp_active is not None:
            qp_energies_abs = eps_qp_active
            qp_energies_beta_abs = eps_qp_active
            qp_energies_rel = eps_qp_active - e_fermi_raw
            qp_energies_beta_rel = eps_qp_active - e_fermi_raw
        else:
            qp_occ_shift_abs = qp_homo - eps[homo_index]
            qp_virt_shift_abs = qp_lumo - eps_beta[homo_index_beta + 1]
            qp_energies_abs = build_qp_energies_vacuum(eps, homo_index, occ_shift=qp_occ_shift_abs, virt_shift=qp_virt_shift_abs)
            qp_energies_beta_abs = build_qp_energies_vacuum(eps_beta, homo_index_beta, occ_shift=qp_occ_shift_abs, virt_shift=qp_virt_shift_abs)
            qp_energies_rel = build_qp_energies(eps_dft_shifted, homo_index, scissor_ev=scissor)
            qp_energies_beta_rel = build_qp_energies(eps_beta_shifted, homo_index_beta, scissor_ev=scissor)
        
        fuzzy_soc_E, fuzzy_soc_E_abs, fuzzy_soc_U, fuzzy_spinor_homo_idx = None, None, None, None
        if args.soc_flag:
            logger.info(f"\n--- Computing SOC Spinor Subspace for Fuzzy Bands (|E-Ef| <= {args.soc_window:.3f} eV) ---")
            logger.info(f"  -> Alpha fuzzy SOC active MOs: {len(fuzzy_active_indices)} / {len(eps)}")
            if is_uks_sp:
                logger.info(f"  -> Beta fuzzy SOC active MOs : {len(fuzzy_active_indices_beta)} / {len(eps_beta)}")
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
            args, C_dense, S, eps_dft_shifted, occ, homo_index, e_homo, e_lumo, e_fermi_raw,
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
    return _export(locals(), (
        "C_dense_beta",
    ))


def _transition_dipoles(args, *,
        C, C_beta, C_dense_beta, bse_n_occ, bse_n_virt, compute_device, dev_obj, eps_beta, homo_index,
        homo_index_beta, is_uks_sp, shells, tracker):
    """Transition dipoles of the active space."""
    # -----------------------------------------------------------------
    # BSE CONTINUATION (Only if run_bse is True)
    # -----------------------------------------------------------------
    tracker.start_stage("Transition Dipoles")
    logger.info("\n--- Computing Transition Dipoles ---")
    if getattr(args, "periodic_enabled", False):
        logger.warning("  [Warning] Periodic mode is enabled, but transition dipoles use the finite-cell AO position operator.")
        logger.warning("  [Warning] Oscillator strengths should be treated as Gamma-only finite-cell approximations.")
    t0_dip = time.time()
    mu_ao_x, mu_ao_y, mu_ao_z = compute_dipole_ao(shells, nthreads=args.nthreads)
    logger.debug(f"  ->  Dipoles computed in {time.time() - t0_dip:.2f} s")

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
        logger.info("  [UKS Spin-Preserving] Computing alpha-alpha and beta-beta transition dipoles separately...")
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
    return _export(locals(), (
        "bse_n_occ_beta", "bse_n_virt_beta", "compute_device", "mu_ia_x", "mu_ia_y", "mu_ia_z", "spin_mode"
    ))


def _solve_excitons(args, *,
        C, C_beta, S, atom_ao_ranges, bse_n_occ, bse_n_occ_beta, bse_n_virt, bse_n_virt_beta, bse_soc_E,
        bse_soc_U, calculated_soc_gap, compute_device, confinement_energy, coords_ang, db_gap, dft_gap,
        eps_beta_shifted, eps_dft_shifted, eps_qp_active, eps_shifted, homo_index, homo_index_beta, mu_ia_x,
        mu_ia_y, mu_ia_z, occ, qp_w, scissor, shared_gamma_ao, shells, spin_mode, stda_gamma_j, stda_gamma_k,
        syms, target_qp_gap, tracker):
    """Spin-free (and SOC) exciton solvers and their analysis."""
    # Store this for the analysis printouts later
    args.qp_gap_num = target_qp_gap

    tracker.start_stage("BSE Exciton Solver (Spin-Free)")
    scissor_solver = 0.0 if eps_qp_active is not None else scissor
    eps_solver = eps_shifted
    eps_dft_solver = eps_dft_shifted

    solver_sf = ExcitonSolver(
        C=C, eps=eps_solver, occ=occ, overlap=S, atom_symbols=syms, atom_coords=np.array(coords_ang),
        atom_ao_ranges=atom_ao_ranges, homo_index=homo_index, n_occ=bse_n_occ, n_virt=bse_n_virt, 
        scissor_ev=scissor_solver, kernel=args.kernel, alpha=args.alpha, beta=args.beta,
        include_exchange=args.include_exchange,
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
        excitation_mode=args.excitation_mode,
        kernel_type=args.kernel_type,
        shared_W=(stda_gamma_j if stda_gamma_j is not None else (qp_w[0] if qp_w is not None else None)),
        shared_gamma_bare=(stda_gamma_k if stda_gamma_k is not None else shared_gamma_ao),
        selection=args.selection, selection_energy=args.selection_energy, selection_pt=args.selection_pt,
        selection_shift=(str(args.selection_shift).lower() != "off"),
        shells=shells,
        eps_dft=eps_dft_solver
    )

    if args.estimate_qp:
        # Pull the dynamically calculated QP gap correction from the solver
        calculated_gap_shift = solver_sf.ham.sigma_virt[0] - solver_sf.ham.sigma_occ[-1]
        
        if getattr(args, 'use_cohsex_gap', False):
            logger.info(f"\n  [QP] OVERRIDE: Using Pure COHSEX Gap Correction ({calculated_gap_shift:.4f} eV) instead of Tabulated GW.")
            scissor = calculated_gap_shift
        else:
            logger.info(f"\n  [QP] Note: Tabulated GW Scissor ({scissor:.4f} eV) was used as the anchor. COHSEX provided orbital dispersion.")
            
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
            C=C, eps=eps_solver, occ=occ, overlap=S, atom_symbols=syms, atom_coords=np.array(coords_ang),
            atom_ao_ranges=atom_ao_ranges, homo_index=homo_index, n_occ=bse_n_occ, n_virt=bse_n_virt, 
            scissor_ev=scissor_solver, kernel=args.kernel, alpha=args.alpha, beta=args.beta,
            include_exchange=args.include_exchange,
            include_direct_eh=args.include_direct_eh,
            estimate_qp=args.estimate_qp, material=args.material, e_thresh=args.e_thresh, f_thresh=args.f_thresh, 
            mu_ia_x=mu_ia_x, mu_ia_y=mu_ia_y, mu_ia_z=mu_ia_z, eps_out=args.eps_out,
            soc_U=bse_soc_U, soc_E=bse_soc_E, device=compute_device, 
            precomputed_sigma=precalc_sigma,
            vxc_ao_path=args.vxc_ao,
            nthreads=args.nthreads,
            spin=spin_mode,
            C_beta=C_beta, eps_beta=eps_beta_shifted, homo_index_beta=homo_index_beta,
            charge_type=args.charge_type,
            n_occ_beta=bse_n_occ_beta, n_virt_beta=bse_n_virt_beta,
            excitation_mode=args.excitation_mode,
            kernel_type=args.kernel_type,
            shared_W=(stda_gamma_j if stda_gamma_j is not None else (qp_w[0] if qp_w is not None else None)),
            shared_gamma_bare=(stda_gamma_k if stda_gamma_k is not None else shared_gamma_ao),
            shells=shells,
            eps_dft=eps_dft_solver
        )
 
        run_solver_and_analysis(solver_soc, np.array(coords_ang), syms, shells, mu_ia_x, mu_ia_y, mu_ia_z, 
                                dft_gap, scissor, confinement_energy, args, suffix="_soc", soc_gap=calculated_soc_gap, soc_U=bse_soc_U, soc_E=bse_soc_E)

    tracker.end_stage()
    tracker.print_summary(device=args.device, nthreads=args.nthreads)
    logger.info("\nAll calculations finished successfully.")
    return {}


def _export(namespace, names):
    """The named variables of a stage that exist (some are set only on certain paths)."""
    return {k: namespace[k] for k in names if k in namespace}


def _run_stage(stage, args, state):
    """Call a stage with the variables it declares (keyword-only parameters) and store what it returns."""
    params = [p.name for p in inspect.signature(stage).parameters.values() if p.kind is p.KEYWORD_ONLY]
    state.update(stage(args, **{p: state.get(p) for p in params}))


def main():
    parser = _build_parser()
    loaded = _load_arguments(parser)
    if loaded is None:
        return
    args = loaded["args"]
    state = {"parser": parser, "config_path": loaded["config_path"]}
    for stage in (_prepare_run, _read_geometry_and_basis, _read_molecular_orbitals, _quasiparticle_correction,
                  _qp_levels_and_kernel, _ip_ea_and_energy_axis, _orbital_populations, _active_space_and_soc,
                  _cubes_and_fuzzy):
        _run_stage(stage, args, state)

    tracker = state["tracker"]
    # -----------------------------------------------------------------
    # EARLY EXIT LOGIC (If run_bse is False)
    # -----------------------------------------------------------------
    run_bse = getattr(args, 'run_bse', True)
    if not run_bse:
        tracker.end_stage()
        tracker.print_summary(device=args.device, nthreads=args.nthreads)
        logger.info("\n--- BSE Calculation Skipped (run_bse: false) ---")
        logger.info("\nAll requested tasks finished successfully.")
        return

    for stage in (_transition_dipoles, _solve_excitons):
        _run_stage(stage, args, state)


if __name__ == "__main__":
    main()
