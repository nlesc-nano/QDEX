"""Layout of the QDEX YAML configuration.

The input file is organised by what each block controls::

    system:          files, material, threads, device
    environment:     dielectric environment (eps_out)
    quasiparticles:  the QP correction of the orbital energies
    integrals:       representation of the two-electron integrals
    excitations:     excited-state framework, kernel and active space
    soc:             spin-orbit coupling
    analysis:        populations, PDOS/COOP, fuzzy bands
    output:          spectra, cubes, CSV, NTOs
    periodic, auger, namd:  unchanged

Inside the new sections keys have short names (``quasiparticles.model``
instead of ``qp_gap``); the full command-line names are accepted as well.
Files written for the old layout (``physics:``, ``solver:``, ``fuzzy:``)
still work: their keys are the command-line names.
"""

# section -> {yaml key: argparse dest}.  Full dest names are always accepted too.
SECTIONS = {
    "system": {
        "mo_file": "mo_file",
        "mo_file_beta": "mo_file_beta",
        "xyz": "xyz",
        "basis_txt": "basis_txt",
        "basis_name": "basis_name",
        "material": "material",
        "gth_file": "gth_file",
        "nthreads": "nthreads",
        "device": "device",
        "cache_mos": "cache_mos",
        "skip_orthonormality_check": "skip_orthonormality_check",
        "orthonormality_tol": "orthonormality_tol",
        "log_file": "log_file",
        "inorganic_elements": "inorganic_elements",
    },
    "environment": {
        "eps_out": "eps_out",
    },
    "quasiparticles": {
        "model": "qp_gap",
        "reference": "qp_reference",
        "bulk_vertex": "bulk_vertex",
        "bulk_vertex_factor": "bulk_vertex_factor",
        "bulk_geometry": "bulk_geometry",
        "bulk_residual": "bulk_residual",
        "z": "qp_z",
        "selfenergy": "qp_selfenergy",
        "levels": "qp_levels",
        "window": "qp_window",
        "window_size": "qp_window_size",
        "solvent_term": "qp_solvent_term",
        "anchor_residual": "qp_anchor_residual",
        "anchor_calibrate": "qp_anchor_calibrate",
        "anchor_table": "qp_anchor_table",
        "residual_scaling": "qp_residual_scaling",
        "residual_power": "qp_residual_power",
        "edge_split": "qp_edge_split",
        "energy_reference": "qp_energy_reference",
        "radius": "qp_radius",
        "polarization": "qp_polarization",
        "regularization_length": "qp_regularization_length",
        "strict": "qp_strict",
        "update_orbitals": "update_orbitals",
        "estimate_qp": "estimate_qp",
        "use_cohsex_gap": "use_cohsex_gap",
        "vxc_ao": "vxc_ao",
        "dynamic_z": "dynamic_z",
    },
    "integrals": {
        "representation": "kernel_type",
        "charges": "charge_type",
        "mnok_exponent": "mnok_exponent",
        "mnok_exponent_exchange": "mnok_exponent_exchange",
        "mnok_onsite": "mnok_onsite",
        "beta": "beta",
    },
    "excitations": {
        "mode": "excitation_mode",
        "kernel": "kernel",
        "kernel_scaling": "alpha",
        "allow_inconsistent_kernel": "allow_inconsistent_kernel",
        "nhomos": "nhomos",
        "nlumos": "nlumos",
        "n_occ": "n_occ",
        "n_virt": "n_virt",
        "e_thresh": "e_thresh",
        "f_thresh": "f_thresh",
        "nroots": "nroots",
        "full_diag": "full_diag",
        "tol": "tol",
        "triplet": "triplet",
        "include_exchange": "include_exchange",
        "include_direct_eh": "include_direct_eh",
        "energy_shift": "soc",
        "selection": "selection",
        "selection_energy": "selection_energy",
        "selection_pt": "selection_pt",
        "selection_shift": "selection_shift",
        "functional": "stda_functional",
        "ax": "stda_ax",
        "stda_alpha": "stda_alpha",
        "stda_beta": "stda_beta",
        "trap_filter": "trap_filter",
        "filter_traps": "trap_filter",
    },
    "soc": {
        "enabled": "soc_flag",
        "window": "soc_window",
        "bse_window": "soc_bse_window",
        "gth_file": "gth_file",
        "gth_functional": "gth_functional",
    },
    "analysis": {k: k for k in (
        "run_fuzzy", "cif", "pdos_atoms", "coop_pairs", "population_print_range", "fuzzy_sigma",
        "pdos_sigma", "ewin", "fold_to_bz", "g_shell", "dashboard_energy_mode",
        "trap_filter", "xi_core_threshold", "xi_core", "xi_trap_threshold", "xi_trap",
        "centroid_core_elements", "f_core_min", "f_core_trap")},
    "output": {k: k for k in (
        "broadening", "sigma", "plot", "show", "cube", "cube_spacing", "disable_cpp_cube", "cube_nhomos",
        "cube_nlumos", "nbse", "bse_states", "write_csv", "csv_roots", "time", "save_xia", "nto",
        "nto_states", "nto_top", "nto_csv", "verbosity", "h5", "html", "excitations_emax",
        "excitations_min_states", "mo_cubes", "mo_cube_spacing", "spectrum_grid", "spectrum_sigmas")},
}

# Accepted spellings of a few dest names.
KEY_ALIASES = {
    "two_electron_integrals": "kernel_type",
    "2e_integrals": "kernel_type",
    "2e-integrals": "kernel_type",
}

# Sections that are not flat key/value blocks.
_NON_SCALAR = ("namd", "auger", "periodic")

# Sections of the old layout, still read with command-line key names.
LEGACY_SECTIONS = ("physics", "solver", "fuzzy")


def section_dest(section, key):
    """argparse dest of ``section.key`` or None if the section does not define it."""
    table = SECTIONS.get(section)
    norm = key.replace("-", "_")
    norm = KEY_ALIASES.get(key, KEY_ALIASES.get(norm, norm))
    if table is None:
        return norm if section in LEGACY_SECTIONS or section not in _NON_SCALAR else None
    if key in table:
        return table[key]
    if norm in table:
        return table[norm]
    if norm in table.values():
        return norm
    return None


def flatten_config(config):
    """All scalar settings of a config as {dest: value}, whatever the layout.

    Used by modules that read the YAML directly (NAMD, Auger).  Keys of the
    new sections are translated to their command-line names; keys of the old
    ``physics`` block are passed through.
    """
    flat = {}
    for section, params in (config or {}).items():
        if not isinstance(params, dict) or section in ("namd", "auger", "periodic"):
            continue
        for key, value in params.items():
            dest = section_dest(section, key)
            if dest is not None:
                flat[dest] = value
    return flat


def _owner_table():
    owner = {}
    for section, table in SECTIONS.items():
        for key, dest in table.items():
            owner.setdefault(dest, (section, key))
    return owner


def to_sections(config):
    """Rewrite a config (any layout) in the current layout.

    Keys of the old ``physics``/``solver``/``fuzzy`` blocks are moved to the
    section that owns them; ``namd``, ``auger`` and ``periodic`` are copied.
    """
    owner = _owner_table()
    out = {}
    for section, params in (config or {}).items():
        if section in _NON_SCALAR or not isinstance(params, dict):
            out[section] = params
            continue
        for key, value in params.items():
            dest = section_dest(section, key)
            if dest == "exchange":          # deprecated alias
                dest = "include_direct_eh"
            if dest not in owner:
                raise ValueError(f"'{section}.{key}' has no place in the current layout")
            sec, new_key = owner[dest]
            out.setdefault(sec, {})[new_key] = value
    order = ["system", "environment", "quasiparticles", "integrals", "excitations", "soc", "analysis",
             "output", "periodic", "auger", "namd"]
    return {k: out[k] for k in sorted(out, key=lambda k: order.index(k) if k in order else len(order))}
