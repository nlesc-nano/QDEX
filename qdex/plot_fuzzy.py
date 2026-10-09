import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
import re
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter, median_filter
import logging

logger = logging.getLogger(__name__)

def load_fuzzy(npz_path):
    d = np.load(npz_path, allow_pickle=True)
    Z = np.asarray(d["intensity"], dtype=float)
    centres = np.asarray(d["centres"], dtype=float)
    spinpol = np.asarray(d["spinpol"], dtype=float) if "spinpol" in d else None
    Zn = np.asarray(d["intensity_norm"], dtype=float) if "intensity_norm" in d else None
    vmap = np.asarray(d["soc_energy_map"], dtype=float) if "soc_energy_map" in d else None
    step = int(np.ceil((Z.size / 250_000) ** 0.5))
    if step > 1:
        Z = Z[::step, ::step]
        centres = centres[::step]
        if spinpol is not None:
            spinpol = spinpol[::step, ::step]
        if Zn is not None:
            Zn = Zn[::step, ::step]
        if vmap is not None:
            vmap = vmap[::step, ::step]
    peaks = None
    if "peak_k" in d:   # k index on the original (not downsampled) path axis, like the extent
        peaks = dict(k=np.asarray(d["peak_k"], float), energy=np.asarray(d["peak_energy"], float),
                     weight=np.asarray(d["peak_weight"], float), state=np.asarray(d["peak_state"], int),
                     soc_energy=np.asarray(d["peak_soc_energy"], float) if "peak_soc_energy" in d else None)
    trap = None
    if "state_flag" in d:
        trap = dict(energy=np.asarray(d["state_energy"], float), kpart=np.asarray(d["state_kpart"], float),
                    onband=np.asarray(d["state_onband"], float), flag=np.asarray(d["state_flag"], int))
    return dict(
        trap=trap,
        centres=centres, Z=Z, Z_norm=Zn, soc_energy_map=vmap, peaks=peaks,
        spinpol=spinpol,
        tick_positions=np.asarray(d.get("tick_positions", []), dtype=float), 
        tick_labels=[str(x) for x in d.get("tick_labels", [])],
        extent=np.asarray(d.get("extent", [0.0, float(Z.shape[1]-1), float(centres.min()), float(centres.max())]), dtype=float),
        kpath_frac=np.asarray(d["kpath_frac"], dtype=float) if "kpath_frac" in d else None,
        ewin=np.asarray(d.get("ewin", [float(centres.min()), float(centres.max())]), dtype=float)
    )

def load_pdos_csv(csv_path):
    if not os.path.exists(csv_path): return np.array([]), [], np.empty((0,0))
    df = pd.read_csv(csv_path)
    return df.iloc[:, 0].to_numpy(dtype=float), list(df.columns[1:]), df.iloc[:, 1:].to_numpy(dtype=float)

def load_dict_csv(csv_path):
    if not os.path.exists(csv_path): return np.array([]), [], {}
    df = pd.read_csv(csv_path)
    return df.iloc[:, 0].to_numpy(dtype=float), list(df.columns[1:]), {p: df[p].to_numpy(dtype=float) for p in df.columns[1:]}

def load_ipr_csv(csv_path):
    if not os.path.exists(csv_path): return np.array([]), np.array([])
    df = pd.read_csv(csv_path)
    return df.iloc[:, 0].to_numpy(dtype=float), df.iloc[:, 1].to_numpy(dtype=float)

def prepare_fuzzy_display(Z, fuzzy_display_mode="state_norm", Z_norm=None):
    """Colour data of the fuzzy map: (values, zmin, zmax, vmin_base, vmax).

    state_norm (default): each state's k-profile normalised to the same total weight
    (fuzzy_bands.state_weights), energy-smeared, on a square-root colour scale. The localised
    semicore/surface states no longer dominate the colours and the low-weight tails, which the
    logarithmic scale of 'raw' spreads over four decades, fade into the background.
    Falls back to 'raw' (log10 of the raw weights) for files without the normalised map."""
    if fuzzy_display_mode == "state_norm" and Z_norm is None:
        fuzzy_display_mode = "raw"
    if fuzzy_display_mode == "state_norm":
        Zpos = Z[Z > 1e-9]
        vmin_base = float(np.percentile(Zpos, 5)) if Zpos.size else 1e-6
        v = float(np.percentile(Z_norm, 99.7))
        v = v if v > 0 else 1.0
        return np.sqrt(np.clip(Z_norm / v, 0.0, 1.0)).astype(np.float32), 0.0, 1.0, vmin_base, v

    if fuzzy_display_mode == "raw":
        Zpos = Z[Z > 1e-9]
        vmax = float(np.percentile(Z, 99.9))
        vmin_base = float(np.percentile(Zpos, 5)) if Zpos.size else 1e-6
        Zm = Z.astype(np.float32); Zm[Zm <= 0] = np.nan
        return np.log10(Zm), float(np.log10(vmin_base)), float(np.log10(vmax)), vmin_base, vmax

    if fuzzy_display_mode == "soft_log":
        Zplot = Z.astype(float, copy=True)
        max_z = np.nanmax(Zplot)
        if not np.isfinite(max_z) or max_z <= 0:
            return np.full_like(Zplot, np.nan, dtype=float), 0.0, 1.0, 1.0, 10.0

        floor = 1e-4 * max_z
        Zplot[Zplot < floor] = np.nan
        Zplot = gaussian_filter(Zplot, sigma=(0.4, 0.0), mode="nearest")

        Zscale = np.nanpercentile(Zplot, 99.8)
        if not np.isfinite(Zscale) or Zscale <= 0:
            Zscale = max_z

        Zlog = np.log10(1.0 + 100.0 * Zplot / Zscale)
        zmin = float(np.nanpercentile(Zlog, 5))
        zmax = float(np.nanpercentile(Zlog, 99.8))
        return Zlog, zmin, zmax, floor, Zscale

    if fuzzy_display_mode == "background_subtracted":
        Zplot = Z.astype(float, copy=True)
        background = median_filter(Zplot, size=(15, 5))
        Zplot = Zplot - 0.7 * background
        Zplot[Zplot < 0] = 0.0
        Zplot = gaussian_filter(Zplot, sigma=(0.5, 0.1))
        Zplot[Zplot <= 0] = np.nan

        finite = Zplot[np.isfinite(Zplot) & (Zplot > 0)]
        if finite.size:
            zmin = float(np.log10(np.percentile(finite, 10)))
            zmax = float(np.log10(np.percentile(finite, 99.8)))
        else:
            zmin, zmax = -6.0, 0.0

        return np.log10(Zplot), zmin, zmax, 10.0**zmin, 10.0**zmax

    if fuzzy_display_mode != "column_norm":
        raise ValueError("fuzzy_display_mode must be 'state_norm', 'raw', 'column_norm', 'background_subtracted', or 'soft_log'")

    Zplot = Z.astype(np.float32, copy=True)
    col_ref = np.percentile(Zplot, 99, axis=0)
    col_ref[col_ref <= 0] = 1.0
    Zplot /= col_ref[None, :]
    Zplot[Zplot < 1e-3] = 0.0
    Zplot = gaussian_filter(Zplot, sigma=(0.7, 0.15), mode="nearest")
    Zplot[Zplot <= 0] = np.nan

    finite = Zplot[np.isfinite(Zplot) & (Zplot > 0)]
    if finite.size:
        zmin = float(np.log10(np.percentile(finite, 20)))
        zmax = float(np.log10(np.percentile(finite, 99.7)))
    else:
        zmin, zmax = -6.0, 0.0

    return np.log10(Zplot), zmin, zmax, 10.0**zmin, 10.0**zmax

def build_wireframe_agnostic(atoms):
    if len(atoms) < 2: return [], [], []
    coords = np.array([[a[1], a[2], a[3]] for a in atoms])
    tree = cKDTree(coords)
    dists, _ = tree.query(coords, k=2)
    local_radii = dists[:, 1] * 0.6
    xs, ys, zs = [], [], []
    max_r = np.max(local_radii)
    pairs = tree.query_pairs(max_r * 2.5)
    for i, j in pairs:
        d_ij = np.linalg.norm(coords[i] - coords[j])
        if d_ij <= 1.25 * (local_radii[i] + local_radii[j]):
            xs.extend([coords[i][0], coords[j][0], None])
            ys.extend([coords[i][1], coords[j][1], None])
            zs.extend([coords[i][2], coords[j][2], None])
    return xs, ys, zs

def parse_cube(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    n_atoms = int(lines[2].split()[0])
    origin = np.array([float(x) for x in lines[2].split()[1:4]])
    nx, dx, _, _ = [float(x) for x in lines[3].split()]
    ny, _, dy, _ = [float(x) for x in lines[4].split()]
    nz, _, _, dz = [float(x) for x in lines[5].split()]
    nx, ny, nz = int(nx), int(ny), int(nz)
    atoms = []
    for i in range(abs(n_atoms)):
        parts = lines[6+i].split()
        atoms.append((int(parts[0]), float(parts[1]), float(parts[2]), float(parts[3]), float(parts[4])))
    data_str = " ".join(lines[6+abs(n_atoms):])
    V = np.fromstring(data_str, sep=' ').reshape((nx, ny, nz))
    BOHR_TO_ANG = 0.529177249
    origin *= BOHR_TO_ANG
    dx *= BOHR_TO_ANG; dy *= BOHR_TO_ANG; dz *= BOHR_TO_ANG
    step_x = max(1, nx // 20); step_y = max(1, ny // 20); step_z = max(1, nz // 20)
    V = V[::step_x, ::step_y, ::step_z]
    nx_new, ny_new, nz_new = V.shape
    x = origin[0] + np.arange(nx_new) * (dx * step_x)
    y = origin[1] + np.arange(ny_new) * (dy * step_y)
    z = origin[2] + np.arange(nz_new) * (dz * step_z)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    atoms_ang = [(z, ax*BOHR_TO_ANG, ay*BOHR_TO_ANG, az*BOHR_TO_ANG) for z, charge, ax, ay, az in atoms]
    return X.flatten(), Y.flatten(), Z.flatten(), V.flatten(), atoms_ang

def generate_interactive_plot(prefix="sf", material="DEFAULT", ef=0.0, e_homo=None, e_lumo=None, normalize_coop=False, energy_label="Energy (eV)", output_html=None, fuzzy_display_mode="state_norm", bulk_bs_path=None, bulk_alignment="core_level", bulk_semicore_rel=None, bulk_overlay=True, bulk_cif=None, bulk_semicore_label=None, bulk_unfolded="auto"):
    if prefix.startswith("soc"):
        lbl = "SOC"
    elif prefix.startswith("uks"):
        lbl = "UKS"
    else:
        lbl = "Spin-Free"
    logger.info(f"  [Plotter] Generating elegant {lbl} HTML dashboard ({energy_label})...")
    
    # =========================================================================
    # PART 1: 2D DASHBOARD
    # =========================================================================
    fig = make_subplots(
        rows=1, cols=5, 
        shared_yaxes=True, 
        column_widths=[0.32, 0.19, 0.17, 0.10, 0.22],
        horizontal_spacing=0.015, 
        subplot_titles=(f"{lbl} Fuzzy Bands", "PDOS", "Trap detector", "Surf/Core", "COOP")
    )
    
    # Header above the panels, laid out in pixels from the top of the figure so that nothing
    # overlaps at any window width: one row of controls, the colour bars stacked over the fuzzy
    # panel, and one vertical legend over each of the panels it belongs to.
    fig_height, top, bottom = 960, 270, 100
    plot_h = fig_height - top - bottom

    def y_px(d):
        """Paper y of a point `d` pixels below the top of the figure."""
        return 1.0 + (top - d) / plot_h

    def x_col(col):
        return float(fig.layout[f"xaxis{'' if col == 1 else col}"].domain[0])

    def colorbar(title, row, **kw):
        """Horizontal colour bar over the fuzzy panel, title on top; row 0 or 1 of the stack."""
        return dict(title=dict(text=title, side="top", font=dict(size=16)), orientation="h", len=0.30, thickness=14,
                    x=0.0, xanchor="left", y=y_px(52 + 70 * row), yanchor="top", tickfont=dict(size=14), **kw)

    def legend(title, col):
        return dict(title=dict(text=f"<b>{title}</b>", font=dict(size=17), side="top"), orientation="v",
                    x=x_col(col), xanchor="left", y=y_px(52), yanchor="top", font=dict(size=15),
                    bgcolor="rgba(255,255,255,0)", tracegroupgap=2)

    fig.update_layout(
        template="plotly_white", paper_bgcolor="white", plot_bgcolor="white",
        height=fig_height,
        font=dict(family="Helvetica, Arial, sans-serif", size=24, color="#222"),
        margin=dict(l=100, r=40, t=top, b=bottom),
        legend=legend("Atoms", 2), legend2=legend("Localization", 4), legend3=legend("Bonds", 5),
        legend4=legend("Trap detector", 3),
    )

    for annotation in fig['layout']['annotations']: annotation['font'] = dict(size=24, family="Helvetica", color="#111")
    palette = ["#636EFA","#EF553B","#00CC96","#AB63FA","#FFA15A","#19D3F3","#FF6692"]
    
    fuzzy = load_fuzzy(f"fuzzy_data_{prefix}.npz")
    pdos_E, pdos_L, pdos_Y = load_pdos_csv(f"pdos_data_{prefix}.csv")
    coop_E, coop_P, coop_V = load_dict_csv(f"coop_data_{prefix}.csv")
    sc_E, sc_P, sc_V = load_dict_csv(f"surf_core_data_{prefix}.csv")
    ipr_E, ipr_V = load_ipr_csv(f"ipr_data_{prefix}.csv")

    Z, ewin = fuzzy["Z"], fuzzy["ewin"]
    kx = np.linspace(fuzzy["extent"][0], fuzzy["extent"][1], Z.shape[1])
    if fuzzy_display_mode == "state_norm" and fuzzy.get("Z_norm") is None:
        fuzzy_display_mode = "raw"
    Z_display, zmin_display, zmax_display, vmin_base, vmax = prepare_fuzzy_display(
        Z, fuzzy_display_mode=fuzzy_display_mode, Z_norm=fuzzy.get("Z_norm"))
    norm_mode = fuzzy_display_mode == "state_norm"
    
    fig.add_shape(type="rect", xref="x", yref="y", x0=kx[0], x1=kx[-1], y0=ewin[0], y1=ewin[1], fillcolor="black", line=dict(width=0), layer="below") 

    cbar_title = "<b>√ weight (per state)</b>" if norm_mode else "<b>log₁₀(I)</b>"
    heat = go.Heatmap(
        z=Z_display, x=kx, y=fuzzy["centres"], colorscale="Inferno",
        zmin=zmin_display, zmax=zmax_display, showscale=True,
        zsmooth="best",
        colorbar=colorbar(cbar_title, 0, **(dict(tickvals=[0, 0.5, 1]) if norm_mode else {})),
        hovertemplate="k-point %{x:.0f}<br>E=%{y:.3f} eV<br>" + ("√w" if norm_mode else "log10(I)") + "=%{z:.2f}<extra></extra>"
    )

    fig.add_trace(heat, row=1, col=1)
    view_traces = {"map": [len(fig.data) - 1], "states": [], "jmap": []}

    # j = 3/2 / 1/2 character of the spinors: <V_SOC> averaged over the states in each pixel
    peaks = fuzzy.get("peaks")
    soc_char = fuzzy.get("soc_energy_map") is not None
    V_TITLE = "<b>⟨V<sub>SOC</sub>⟩ (eV)</b>  − j=½ · + j=3/2"
    v_lim = 0.2
    if soc_char:
        v_all = peaks["soc_energy"] if peaks is not None and peaks["soc_energy"] is not None else fuzzy["soc_energy_map"]
        v_all = v_all[np.isfinite(v_all)]
        if v_all.size:
            v_lim = max(float(np.percentile(np.abs(v_all), 98)), 0.02)
        # j-character view: hue = j character (red j = 3/2, blue j = 1/2, grey mixed or no SOC),
        # brightness = the weight of the map (a translucent overlay on the map washes the colours out)
        bright = Z_display if norm_mode else np.clip((np.nan_to_num(Z_display, nan=zmin_display) - zmin_display)
                                                     / max(zmax_display - zmin_display, 1e-12), 0.0, 1.0)
        t = np.clip(np.nan_to_num(fuzzy["soc_energy_map"] / v_lim, nan=0.0), -1.0, 1.0)
        # a heatmap with a composite colour scale (an image trace would force equal axis scales):
        # n_hue bins from blue to red, each a ramp from black to its colour; value = bin + brightness
        n_hue = 21
        grey, red, blue = np.array([204, 204, 204]), np.array([255, 64, 38]), np.array([64, 140, 255])
        stops = []
        for h in range(n_hue):
            th = 2.0 * h / (n_hue - 1) - 1.0
            c = grey + max(th, 0.0) * (red - grey) + max(-th, 0.0) * (blue - grey)
            stops += [[h / n_hue, "rgb(0,0,0)"], [(h + 0.999) / n_hue, "rgb(%d,%d,%d)" % tuple(c.round().astype(int))]]
        stops.append([1.0, stops[-1][1]])
        hue = np.rint((t + 1.0) / 2.0 * (n_hue - 1))
        code = hue + 0.999 * np.clip(np.nan_to_num(bright, nan=0.0), 0.0, 1.0)
        fig.add_trace(go.Heatmap(z=code.astype(np.float32), x=kx, y=fuzzy["centres"], colorscale=stops, zmin=0.0, zmax=float(n_hue),
                                 zsmooth=False, showscale=False, visible=False, hoverinfo="skip"), row=1, col=1)
        view_traces["jmap"].append(len(fig.data) - 1)
        # colour bar of the j-character view
        fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers", showlegend=False, visible=False, hoverinfo="skip",
                                 marker=dict(color=[0.0], colorscale=[[0, "rgb(64,140,255)"], [0.5, "rgb(204,204,204)"], [1, "rgb(255,64,38)"]],
                                             cmin=-v_lim, cmax=v_lim, showscale=True,
                                             colorbar=colorbar(V_TITLE, 1, tickvals=[-v_lim, 0.0, v_lim], ticktext=[f"{-v_lim:.2f}", "0", f"+{v_lim:.2f}"]))), row=1, col=1)
        view_traces["jmap"].append(len(fig.data) - 1)

    # Dominant k of every state: one marker per k-peak at the state's own energy (no broadening)
    if peaks is not None and peaks["k"].size:
        # Kramers pairs of spinors (and degenerate MOs) give the same marker twice: keep one
        _, first = np.unique(np.stack([peaks["k"], np.round(peaks["energy"], 4)], axis=1), axis=0, return_index=True)
        keep = np.zeros(peaks["k"].size, bool); keep[first] = True
        in_win = keep & (peaks["energy"] >= ewin[0]) & (peaks["energy"] <= ewin[1])
        w = peaks["weight"]
        w_ref = max(float(np.percentile(w[in_win], 99)) if in_win.any() else float(w.max()), 2.0)
        size = 2.5 + 6.5 * np.sqrt(np.clip((w - 1.5) / (w_ref - 1.5), 0.0, 1.0))
        v_pk = peaks["soc_energy"]
        groups = [(in_win & np.isfinite(v_pk), True), (in_win & ~np.isfinite(v_pk), False)] if v_pk is not None else [(in_win, False)]
        for sel, coloured in groups:
            if not sel.any():
                continue
            if coloured:
                marker = dict(size=size[sel], color=v_pk[sel], colorscale="RdBu_r", cmin=-v_lim, cmax=v_lim,
                              line=dict(width=0.3, color="rgba(0,0,0,0.5)"), opacity=0.9, showscale=True,
                              colorbar=colorbar(V_TITLE, 1, tickvals=[-v_lim, 0.0, v_lim], ticktext=[f"{-v_lim:.2f}", "0", f"+{v_lim:.2f}"]))
                custom = np.stack([peaks["state"][sel], w[sel], v_pk[sel]], axis=1)
                hover = "state %{customdata[0]:.0f}<br>E=%{y:.3f} eV, k-point %{x:.0f}<br>weight %{customdata[1]:.1f}× mean<br>⟨V_SOC⟩=%{customdata[2]:.3f} eV<extra></extra>"
            else:
                marker = dict(size=size[sel], color="rgba(235,235,235,0.9)", line=dict(width=0.3, color="rgba(0,0,0,0.5)"))
                custom = np.stack([peaks["state"][sel], w[sel]], axis=1)
                hover = "state %{customdata[0]:.0f}<br>E=%{y:.3f} eV, k-point %{x:.0f}<br>weight %{customdata[1]:.1f}× mean<extra></extra>"
            fig.add_trace(go.Scattergl(x=peaks["k"][sel], y=peaks["energy"][sel], mode="markers", marker=marker,
                                       customdata=custom, hovertemplate=hover, showlegend=False, visible=False),
                          row=1, col=1)
            view_traces["states"].append(len(fig.data) - 1)

    if fuzzy.get("spinpol") is not None:
        spinpol = fuzzy["spinpol"].astype(float, copy=True)
        spinpol[Z <= max(vmin_base, 1e-12)] = np.nan
        fig.add_trace(
            go.Heatmap(
                z=spinpol, x=kx, y=fuzzy["centres"],
                colorscale="RdBu", zmin=-1.0, zmax=1.0,
                opacity=0.38, showscale=True, zsmooth="best",
                colorbar=colorbar("<b>spin polarisation α−β</b>", 1, tickvals=[-1, 0, 1]),
                hovertemplate="k-point %{x:.0f}<br>E=%{y:.3f} eV<br>spin pol=%{z:.2f}<extra></extra>"
            ),
            row=1, col=1
        )

    # Bulk bands overlay (CP2K reference band structure, spin-free DFT)
    has_bulk = False
    if bulk_overlay:
        try:
            from qdex.bulk_bands import get_aligned_bulk_bands, functional_label
            fl = functional_label()
            bulk_data = get_aligned_bulk_bands(
                material=material,
                bs_path=bulk_bs_path,
                alignment_mode=bulk_alignment,
                ewin=(float(ewin[0]), float(ewin[1])),
                qd_homo_rel=e_homo,
                qd_lumo_rel=e_lumo,
                qd_semicore_rel=bulk_semicore_rel,
                path_frac=fuzzy["kpath_frac"],
                soc=(lbl == "SOC"),
                cif=bulk_cif,
                qd_semicore_label=bulk_semicore_label,
                unfolded=bulk_unfolded,
            )
            if bulk_data is not None:
                if bulk_data["segments"] is not None:
                    segments = bulk_data["segments"]
                else:   # unknown path: stretch the bulk k-points over the fuzzy axis
                    segments = [(np.linspace(fuzzy["extent"][0], fuzzy["extent"][1], bulk_data["n_k"]),
                                 bulk_data["bands_aligned"])]
                if bulk_data.get("unfolded"):
                    # unfolded supercell bands: one marker per (k, band), opacity = unfolding weight
                    bulk_name = f"Bulk {fl}+SOC, unfolded" if bulk_data.get("soc") else f"Bulk {fl}, unfolded"
                    xs, ys, cs = [], [], []
                    for x_seg, b_seg, w_seg in segments:
                        mask = w_seg >= bulk_data.get("min_weight", 0.03)
                        xx = np.broadcast_to(np.asarray(x_seg, float)[:, None], b_seg.shape)
                        xs.append(xx[mask]); ys.append(b_seg[mask]); cs.append(w_seg[mask])
                    xs, ys, cs = np.concatenate(xs), np.concatenate(ys), np.clip(np.concatenate(cs), 0, 1)
                    order = np.argsort(cs)              # strongest on top
                    fig.add_trace(go.Scattergl(
                        x=xs[order], y=ys[order], mode="markers", name=bulk_name, legendgroup="bulk_bands",
                        marker=dict(size=3.5 + 3.0 * cs[order], color=[f"rgba(0,240,255,{0.15 + 0.85 * c:.2f})" for c in cs[order]],
                                    line=dict(width=0)),
                        customdata=cs[order].tolist(),   # plain list: read back by the bulk-band controls
                        hovertemplate=f"{bulk_name}: E = %{{y:.3f}} eV, weight %{{customdata:.2f}}<extra></extra>"),
                        row=1, col=1)
                    logger.info(f"  [Plotter] Overlaid unfolded bulk bands from {os.path.basename(bulk_data['source_file'])} "
                                f"({bulk_data['alignment_mode']} alignment).")
                elif bulk_data.get("soc"):
                    bulk_name = f"Bulk {fl}+SOC"
                else:
                    bulk_name = f"Bulk {fl}, no SOC" if lbl == "SOC" else f"Bulk {fl}"
                first = True
                for x_seg, b_seg in ([] if bulk_data.get("unfolded") else segments):
                    for b_i in range(b_seg.shape[1]):
                        fig.add_trace(
                            go.Scatter(
                                x=x_seg,
                                y=b_seg[:, b_i],
                                mode="lines",
                                line=dict(color="rgba(0, 240, 255, 0.85)", width=2.0),
                                name=bulk_name,
                                legendgroup="bulk_bands",
                                showlegend=first,
                                hovertemplate=f"{bulk_name}: E = %{{y:.3f}} eV<extra></extra>",
                            ),
                            row=1, col=1
                        )
                        first = False
                has_bulk = True
                if not bulk_data.get("unfolded"):
                    logger.info(f"  [Plotter] Overlaid {len(bulk_data['band_indices'])} bulk {fl} bands from "
                                f"{os.path.basename(bulk_data['source_file'])} ({bulk_data['alignment_mode']} alignment).")
        except Exception as exc:
            logger.debug(f"  [Plotter] Bulk band overlay skipped: {exc}")

    if fuzzy["tick_positions"].size:
        # FIX: Scale ticks based on the original physical width, not the downsampled pixel count
        original_width = fuzzy["extent"][1] - fuzzy["extent"][0]
        scale = (kx[-1] - kx[0]) / original_width if original_width > 0 else 1.0
        tpos_plot = kx[0] + fuzzy["tick_positions"] * scale
        
        for x in tpos_plot: fig.add_vline(x=x, line_color="rgba(255,255,255,0.4)", line_width=2, row=1, col=1)
        fig.update_xaxes(tickmode="array", tickvals=tpos_plot, ticktext=fuzzy["tick_labels"], row=1, col=1)
 
    if len(pdos_L) > 0:
        for j, lab in enumerate(pdos_L):
            fig.add_trace(go.Scatter(x=pdos_Y[:, j], y=pdos_E, mode="lines", fill="tonextx" if j > 0 else "tozerox", line=dict(width=1.0, color="rgba(0,0,0,0)"), fillcolor=palette[j % len(palette)], name=lab, showlegend=True, legend="legend", hovertemplate=f"{lab}: %{{x:.3f}}<br>E=%{{y:.3f}} eV<extra></extra>"), row=1, col=2)
        total = pdos_Y[:, -1]
        fig.add_trace(go.Scatter(x=total, y=pdos_E, mode="lines", line=dict(color="black", width=3), name="Total DOS", showlegend=False), row=1, col=2)
        fig.update_xaxes(range=[0, float(max(total.max(), 1e-12)) * 1.05], row=1, col=2)

    # Trap detector: per-state localization measures on one 0..1 axis (see the glossary)
    trap = fuzzy.get("trap")
    if len(ipr_E) > 0:
        fig.add_trace(go.Scatter(x=ipr_V, y=ipr_E, mode="markers", name="IPR", legend="legend4",
                                 marker=dict(size=7, color="#440154", symbol="circle", line=dict(width=1, color="black")),
                                 hovertemplate="IPR: %{x:.4f}<br>E: %{y:.3f} eV<extra></extra>"), row=1, col=3)
    deloc = {}
    if trap is not None:
        is_trap = trap["flag"] == 1
        for key, name, color, symbol in (("kpart", "k-particip.", "#E69F00", "diamond"),
                                         ("onband", "on-band", "#009E73", "triangle-left")):
            fig.add_trace(go.Scatter(
                x=trap[key], y=trap["energy"], mode="markers", name=name, legend="legend4",
                marker=dict(size=np.where(is_trap, 10, 7), color=color, symbol=symbol,
                            line=dict(width=np.where(is_trap, 2.5, 1), color=np.where(is_trap, "#D7263D", "black"))),
                customdata=np.where(is_trap, "trap", ""),
                hovertemplate=f"{name}: %{{x:.2f}}<br>E: %{{y:.3f}} eV<br>%{{customdata}}<extra></extra>"), row=1, col=3)
        for flag, key in ((2, "homo"), (3, "lumo")):
            hit = np.where(trap["flag"] == flag)[0]
            if hit.size:
                deloc[key] = float(trap["energy"][hit[0]])
    if len(ipr_E) > 0 or trap is not None:
        fig.update_xaxes(range=[0, 1.05], row=1, col=3)

    if len(sc_E) > 0:
        surf, core = sc_V["Surface"], sc_V["Core"]
        xs_c, ys_c, xs_s, ys_s = [], [], [], []
        for yi, cv, sv in zip(sc_E, core, surf):
            xs_c.extend([0.0, float(cv), None]); ys_c.extend([float(yi), float(yi), None])
            xs_s.extend([float(cv), float(cv+sv), None]); ys_s.extend([float(yi), float(yi), None])
        fig.add_trace(go.Scattergl(x=xs_c, y=ys_c, mode="lines", line=dict(color="#1f77b4", width=3), name="Core", legend="legend2"), row=1, col=4)
        fig.add_trace(go.Scattergl(x=xs_s, y=ys_s, mode="lines", line=dict(color="#ff7f0e", width=3), name="Surface", legend="legend2"), row=1, col=4)
        fig.update_xaxes(range=[0, 1.05], row=1, col=4)

    if len(coop_P) > 0:
        mask = (coop_E >= ewin[0]) & (coop_E <= ewin[1])
        E_sticks = coop_E[mask]
        scale, vmax_abs = 1.0, 1.0
        if normalize_coop:
            gmax = max((np.nanmax(np.abs(coop_V[p][mask])) for p in coop_P if coop_V[p][mask].size), default=0.0)
            scale = (1.0 / gmax) if gmax > 0 else 1.0
            fig.update_xaxes(range=[-1.05, 1.05], row=1, col=5)
        else:
            vmax_abs = max((np.nanmax(np.abs(coop_V[p][mask])) for p in coop_P if coop_V[p][mask].size), default=1.0)
            fig.update_xaxes(range=[-(vmax_abs*1.1), (vmax_abs*1.1)], row=1, col=5)
        for i, p in enumerate(coop_P):
            v = (coop_V[p][mask] * scale).astype(np.float32)
            if v.size == 0: continue
            xs, ys = [], []
            for yi, xv in zip(E_sticks, v):
                xs.extend([0.0, float(xv), None]); ys.extend([float(yi), float(yi), None])
            fig.add_trace(go.Scattergl(x=xs, y=ys, mode="lines", line=dict(color=palette[i % len(palette)], width=4), name=p, legend="legend3"), row=1, col=5)

    gap_mid = 0.5 * (e_homo + e_lumo) if e_homo is not None and e_lumo is not None else None
    reference_y = ef if ef is not None else gap_mid

    for col in range(1, 6):
        line_col_ref = "white" if col == 1 else "rgba(0,0,0,0.4)"
        if reference_y is not None:
            fig.add_hline(y=reference_y, line_dash="dash", line_color=line_col_ref, line_width=2.5, row=1, col=col, layer="above")
        if e_homo is not None: fig.add_hline(y=e_homo, line_dash="dot", line_color="royalblue", line_width=3, row=1, col=col, layer="above")
        if e_lumo is not None: fig.add_hline(y=e_lumo, line_dash="dot", line_color="crimson", line_width=3, row=1, col=col, layer="above")
        if "homo" in deloc: fig.add_hline(y=deloc["homo"], line_dash="solid", line_color="royalblue", line_width=2.5, row=1, col=col, layer="above")
        if "lumo" in deloc: fig.add_hline(y=deloc["lumo"], line_dash="solid", line_color="crimson", line_width=2.5, row=1, col=col, layer="above")

    if e_homo is not None and e_lumo is not None:
        gap = e_lumo - e_homo
        # in the PDOS panel, which is empty around the gap (the WebGL markers of the fuzzy panel would
        # cover them): the gap on the mid-gap line, HOMO below its line, LUMO above its line
        fig.add_annotation(x=0.97, y=reference_y, xref="x2 domain", yref="y", text=f"<b>E<sub>g</sub> = {gap:.3f} eV</b>", showarrow=False, font=dict(color="#222", size=17), xanchor="right", yanchor="middle", bgcolor="white")
        fig.add_annotation(x=0.97, y=e_homo, xref="x2 domain", yref="y", text="<b>HOMO</b>", showarrow=False, font=dict(color="royalblue", size=16), xanchor="right", yanchor="top", yshift=-4)
        fig.add_annotation(x=0.97, y=e_lumo, xref="x2 domain", yref="y", text="<b>LUMO</b>", showarrow=False, font=dict(color="crimson", size=16), xanchor="right", yanchor="bottom", yshift=4)
    if "homo" in deloc and "lumo" in deloc:
        # labels only where the delocalized edge differs from the nominal one (they would overlap)
        if e_homo is None or abs(deloc["homo"] - e_homo) > 0.02:
            fig.add_annotation(x=0.97, y=deloc["homo"], xref="x2 domain", yref="y", text="<b>deloc. HOMO</b>", showarrow=False, font=dict(color="royalblue", size=14), xanchor="right", yanchor="top", yshift=-4)
        if e_lumo is None or abs(deloc["lumo"] - e_lumo) > 0.02:
            fig.add_annotation(x=0.97, y=deloc["lumo"], xref="x2 domain", yref="y", text="<b>deloc. LUMO</b>", showarrow=False, font=dict(color="crimson", size=14), xanchor="right", yanchor="bottom", yshift=4)
        if reference_y is not None:
            fig.add_annotation(x=0.97, y=reference_y, xref="x2 domain", yref="y",
                               text=f"deloc. E<sub>g</sub> = {deloc['lumo'] - deloc['homo']:.3f} eV", showarrow=False,
                               font=dict(color="#222", size=14), xanchor="right", yanchor="top", yshift=-14, bgcolor="white")

    for col in range(1, 6):
        fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True, ticks="outside", gridcolor='rgba(0,0,0,0.1)', zeroline=False, row=1, col=col)
        fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True, ticks="outside", gridcolor='rgba(0,0,0,0.1)', zeroline=False, row=1, col=col)
        
    fig.update_yaxes(range=[ewin[0], ewin[1]], title_text=f"<b>{energy_label}</b>", title_font=dict(size=28), row=1, col=1)
    fig.update_xaxes(title_text="<b>k-Path</b>", title_font=dict(size=28), tickangle=0, tickfont=dict(size=17), row=1, col=1)
    fig.update_xaxes(title_text="<b>DOS</b>", title_font=dict(size=28), row=1, col=2)
    fig.update_xaxes(title_text="<b>0 – 1</b>", title_font=dict(size=28), tickvals=[0, 0.5, 1.0], row=1, col=3)
    fig.update_xaxes(title_text="<b>Char</b>", title_font=dict(size=28), tickvals=[0, 0.5, 1.0], row=1, col=4)
    fig.update_xaxes(title_text="<b>COOP</b>", title_font=dict(size=28), row=1, col=5)

    if norm_mode:   # saturate the colour scale at 1/s of the reference weight
        buttons = [dict(label=f"Contrast: {s}x", method="restyle", args=[{"zmax": [float(1.0 / np.sqrt(s))]}, [0]]) for s in (1, 2, 4, 8)]
    else:
        buttons = [dict(label=f"Contrast: {s}x", method="restyle", args=[{"zmin": [float(np.log10(max(vmin_base, vmax/s)))]}, [0]]) for s in (1, 10, 100, 1000)]
    # controls in the top row: contrast at the left, the views next to it (offset in pixels)
    menus = [dict(type="dropdown", buttons=buttons, x=0.0, y=y_px(6), xanchor="left", yanchor="top",
                  pad=dict(l=0, t=0), font=dict(size=16), bgcolor="#f8f9fa", bordercolor="black")]

    # Views of the fuzzy panel: smeared map, map + k-peak markers of the states, markers only,
    # and (SOC) the map coloured by the j = 3/2 / 1/2 character
    views = []
    if view_traces["states"]:
        views += [("Map", ["map"]), ("Map + states", ["map", "states"]), ("States", ["states"])]
    if view_traces["jmap"]:
        views += [("j-character", ["jmap"])] if views else [("Map", ["map"]), ("j-character", ["jmap"])]
    if views:
        managed = view_traces["map"] + view_traces["states"] + view_traces["jmap"]
        def _visible(parts):
            on = {i for part in parts for i in view_traces[part]}
            return [i in on for i in managed]
        default = 1 if view_traces["states"] else 0
        for i, vis in zip(managed, _visible(views[default][1])):
            fig.data[i].visible = vis
        menus.append(dict(type="buttons", direction="right", active=default, x=0.0, y=y_px(6), xanchor="left", yanchor="top",
                          pad=dict(l=150, t=0), font=dict(size=16), bgcolor="#f8f9fa", bordercolor="black", showactive=True,
                          buttons=[dict(label=lab, method="restyle", args=[{"visible": _visible(parts)}, managed])
                                   for lab, parts in views]))
    fig.update_layout(updatemenus=menus)

    plot_2d_html = fig.to_html(full_html=False, include_plotlyjs=False, div_id="fuzzy_2d_plot",
                               config={'responsive': True, 'displaylogo': False})
    # Bulk-band controls (HTML, above the plot): show/hide, colour, line width and opacity
    bulk_controls_html = '''
            <div class="iso-control-bar bulk-control-bar">
                <strong>Bulk bands:</strong>
                <label><input type="checkbox" id="bulkShow" checked onchange="updateBulk()"> show</label>
                <label>colour <input type="color" id="bulkColor" value="#00f0ff" oninput="updateBulk()"></label>
                <label>width <input type="range" id="bulkWidth" min="0.5" max="6" step="0.5" value="2" oninput="updateBulk()">
                    <span id="bulkWidthVal">2.0</span></label>
                <label>opacity <input type="range" id="bulkAlpha" min="0.1" max="1" step="0.05" value="0.85" oninput="updateBulk()">
                    <span id="bulkAlphaVal">0.85</span></label>
            </div>''' if has_bulk else ""



    # =========================================================================
    # PART 2: 3D MOLECULAR ORBITAL (CUBE) VISUALIZATION
    # =========================================================================
    if prefix == "soc":
        cube_files = [f for f in os.listdir('.') if f.lower().startswith('spinor_') and f.lower().endswith('.cube')]
    else:
        cube_files = [f for f in os.listdir('.') if f.lower().startswith('spatial_') and f.lower().endswith('.cube')]

    def extract_idx(name):
        name = name.upper()
        if "HOMO" in name:
            match = re.search(r'HOMO-(\d+)', name)
            return -int(match.group(1)) if match else 0
        elif "LUMO" in name:
            match = re.search(r'LUMO\+(\d+)', name)
            return 1 + (int(match.group(1)) if match else 0)
            
        # Fallback for other numbered files
        match = re.search(r'\d+', name)
        return int(match.group()) if match else 0
    
    cube_files.sort(key=extract_idx)
    
    plot_3d_html = ""
    first_iso_val = 0.0001
    
    if len(cube_files) > 0:
        # nominal frontier cubes (the four around the gap) plus the delocalized band-edge cubes
        # (spatial_MO_HOMO-20_dHOMO.cube, output.cube_nhomos_deloc / cube_nlumos_deloc) as extra panels
        deloc = [f for f in cube_files if re.search(r"_d(HOMO|LUMO)", f)]
        nominal = [f for f in cube_files if f not in deloc]
        if len(nominal) >= 4:
            mid = len(nominal) // 2
            nominal = nominal[mid-2 : mid+2]
        cube_files = nominal + deloc

        def cube_title(f):
            t = os.path.basename(f).replace('.cube', '')
            t = re.sub(r'^(spatial_|spinor_)(MO_|sp_)?', '', t).replace('_density', '')
            m = re.match(r'(.+?)_(d(?:HOMO|LUMO)[-+0-9]*)$', t)
            return f"{m.group(2)} ({m.group(1)})" if m else t
        cube_titles = [cube_title(f) for f in cube_files]
        # one row of nominal cubes, a second row with the delocalized band edges
        rows_of = [nominal, deloc] if (nominal and deloc) else [cube_files]
        n_cols = max(len(r) for r in rows_of)
        cube_pos = [(ri + 1, ci + 1) for ri, r in enumerate(rows_of) for ci in range(len(r))]
            
        fig_3d = make_subplots(
            rows=len(rows_of), cols=n_cols,
            specs=[[{'type': 'scene'} if ci < len(r) else None for ci in range(n_cols)] for r in rows_of],
            subplot_titles=cube_titles,
            horizontal_spacing=0.02, vertical_spacing=0.06
        )
        
        for idx, cfile in enumerate(cube_files):
            row_i, col_i = cube_pos[idx]
            logger.info(f"    -> Parsing and mapping {cfile} to 3D grid...")
            X_grid, Y_grid, Z_grid, V_data, atoms = parse_cube(cfile)
            
            v_max = float(np.max(V_data))
            v_min = float(np.min(V_data))
            is_density = (v_min >= -1e-6)
            
            if prefix == "soc" or is_density:
                iso_val = 0.0001
                iso_pos = min(iso_val, v_max * 0.95)
                if iso_pos < 1e-6: iso_pos = v_max * 0.50
                fig_3d.add_trace(go.Isosurface(
                    x=X_grid, y=Y_grid, z=Z_grid, value=V_data,
                    isomin=iso_pos, isomax=iso_pos, surface_count=1,
                    colorscale=[[0, 'blue'], [1, 'blue']], showscale=False,
                    caps=dict(x_show=False, y_show=False, z_show=False), opacity=0.6, name="Density"
                ), row=row_i, col=col_i)
            else:
                iso_val = max(abs(v_max), abs(v_min)) * 0.15
                iso_pos = min(iso_val, v_max * 0.95)
                iso_neg = max(-iso_val, v_min * 0.95)
                
                if iso_pos > 1e-6:
                    fig_3d.add_trace(go.Isosurface(
                        x=X_grid, y=Y_grid, z=Z_grid, value=V_data,
                        isomin=iso_pos, isomax=iso_pos, surface_count=1,
                        colorscale=[[0, 'blue'], [1, 'blue']], showscale=False,
                        caps=dict(x_show=False, y_show=False, z_show=False), opacity=0.35, name="Pos Lobe",
                        flatshading=False,
                    ), row=row_i, col=col_i)
                
                if iso_neg < -1e-6:
                    fig_3d.add_trace(go.Isosurface(
                        x=X_grid, y=Y_grid, z=Z_grid, value=V_data,
                        isomin=iso_neg, isomax=iso_neg, surface_count=1,
                        colorscale=[[0, 'red'], [1, 'red']], showscale=False,
                        caps=dict(x_show=False, y_show=False, z_show=False), opacity=0.35, name="Neg Lobe",
                        flatshading=False,
                    ), row=row_i, col=col_i)
                    
            if idx == 0: first_iso_val = iso_val
            
            b_xs, b_ys, b_zs = build_wireframe_agnostic(atoms)
            fig_3d.add_trace(go.Scatter3d(x=b_xs, y=b_ys, z=b_zs, mode='lines', line=dict(color='#555555', width=3), showlegend=False, hoverinfo='skip'), row=row_i, col=col_i)
            
            atom_colors = {1: '#FFFFFF', 6: '#777777', 7: '#0000FF', 8: '#FF0000', 16: '#CCCC00', 14: '#FFC0CB', 15: '#FFA500'}
            ax, ay, az, ac = [], [], [], []
            for a in atoms:
                ax.append(a[1]); ay.append(a[2]); az.append(a[3]); ac.append(atom_colors.get(a[0], '#A0A0A0'))
            fig_3d.add_trace(go.Scatter3d(x=ax, y=ay, z=az, mode='markers', marker=dict(size=4, color=ac, line=dict(width=1, color='black')), showlegend=False, hoverinfo='skip'), row=row_i, col=col_i)
            
        scene_config = dict(xaxis=dict(showbackground=False, showgrid=False, zeroline=False, showticklabels=False, title=''), yaxis=dict(showbackground=False, showgrid=False, zeroline=False, showticklabels=False, title=''), zaxis=dict(showbackground=False, showgrid=False, zeroline=False, showticklabels=False, title=''), aspectmode='data')
        layout_update = {f"scene{i+1}" if i > 0 else "scene": scene_config for i in range(len(cube_files))}
        
        fig_3d.update_layout(
            **layout_update, paper_bgcolor="#fdfdfd", plot_bgcolor="#fdfdfd",
            margin=dict(l=10, r=10, t=50, b=10), height=500 * len(rows_of), font=dict(family="Helvetica, Arial, sans-serif", size=24, color="#222")
        )
        
        plot_3d_html = fig_3d.to_html(full_html=False, include_plotlyjs=False, div_id="mo_3d_plot")

    # =========================================================================
    # PART 3: WRAPPING HTML TEMPLATE AND GLOSSARY WITH NATIVE JS INPUT
    # =========================================================================
    html_template = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>QDEX {lbl} Electronic Structure Analysis - {energy_label}</title>
        <script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
        <style>
            body {{ font-family: 'Helvetica', 'Arial', sans-serif; background-color: #f8f9fa; margin: 0; padding: 20px; color: #333; }}
            .dashboard-container {{ max-width: 2400px; margin: 0 auto; background: white; padding: 30px; border-radius: 12px; box-shadow: 0 4px 15px rgba(0,0,0,0.05); }}
            .plot-container {{ width: 100%; margin-bottom: 20px; }}
            .mo-container {{ width: 100%; background: #fdfdfd; border: 1px solid #eaeaea; border-radius: 8px; padding: 10px 0; margin-bottom: 40px; position: relative; }}
            .iso-control-bar {{ background: white; padding: 12px 25px; border-radius: 6px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); display: inline-flex; align-items: center; gap: 15px; border: 1px solid #eaeaea; margin: 15px 0 0 25px; z-index: 10; }}
            .iso-control-bar input {{ padding: 6px 12px; font-size: 16px; border: 1px solid #ccc; border-radius: 4px; width: 120px; }}
            .iso-control-bar button {{ padding: 8px 16px; font-size: 16px; background: #1f77b4; color: white; border: none; border-radius: 4px; cursor: pointer; transition: 0.2s; }}
            .iso-control-bar button:hover {{ background: #155d8f; }}
            .bulk-control-bar {{ margin: 0 0 10px 0; font-size: 16px; flex-wrap: wrap; }}
            .bulk-control-bar label {{ display: inline-flex; align-items: center; gap: 6px; }}
            .bulk-control-bar input[type=range] {{ width: 110px; padding: 0; }}
            .bulk-control-bar input[type=color] {{ width: 42px; height: 28px; padding: 0; border: 1px solid #ccc; }}
            .explanation-box {{ background: #fdfdfd; border: 1px solid #eaeaea; border-radius: 8px; padding: 30px; margin-top: 10px; }}
            .explanation-box h2 {{ margin-top: 0; color: #222; font-size: 26px; border-bottom: 2px solid #eee; padding-bottom: 10px; margin-bottom: 20px; }}
            .glossary-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(350px, 1fr)); gap: 25px; }}
            .glossary-item {{ background: #fff; padding: 20px; border-radius: 6px; box-shadow: 0 2px 8px rgba(0,0,0,0.04); border-left: 4px solid #1f77b4; }}
            .glossary-item > strong {{ display: block; font-size: 20px; margin-bottom: 10px; color: #111; }}
            .glossary-item p {{ margin: 0; font-size: 16px; line-height: 1.5; color: #555; }}
            .trap-indicator {{ border-left-color: #ff7f0e; }}
            .cube-indicator {{ border-left-color: #d62728; }}
        </style>
    </head>
    <body>
        <div class="dashboard-container">
            <div class="plot-container">
                {bulk_controls_html}
                {plot_2d_html}
            </div>
            
            {f'''
            <div class="mo-container">
                <div class="iso-control-bar">
                    <label for="isoInput"><strong>Isosurface Threshold (±):</strong></label>
                    <input type="number" id="isoInput" step="0.0001" value="{first_iso_val:.5f}">
                    <button onclick="updateIso()">Apply</button>
                </div>
                {plot_3d_html}
            </div>
            ''' if plot_3d_html else ''}
            
            <div class="explanation-box">
                <h2>Glossary of Electronic Structure Plots</h2>
                <div class="glossary-grid">
                    <div class="glossary-item">
                        <strong>1. Fuzzy Bands (Reciprocal Space)</strong>
                        <p>Projects the finite-size real-space Molecular Orbitals onto a bulk-like momentum ($k$) grid. A sharp, continuous band structure indicates bulk-like delocalized states, while flat, smeared lines across the Brillouin zone indicate highly localized molecules or defects.</p>
                    </div>
                    <div class="glossary-item">
                        <strong>2. PDOS (Projected Density of States)</strong>
                        <p>Shows how much specific elements (or orbitals) contribute to the overall electronic density at a given energy. The peaks correspond to available molecular orbitals.</p>
                    </div>
                    <div class="glossary-item trap-indicator">
                        <strong>3. Trap detector (IPR, k-participation, on-band weight)</strong>
                        <p>Three numbers per state, all between 0 and 1, plotted at the state's energy.</p>
                        <p><span style="color:#440154">&#9679;</span> <strong>IPR</strong> (real space): near 0 for a state spread over many atoms, near 1 for a state on one atom.</p>
                        <p><span style="color:#E69F00">&#9670;</span> <strong>k-participation</strong>: how evenly the state's fuzzy-band weight is spread over a full grid of the Brillouin zone, 1/(N&nbsp;&Sigma;<sub>k</sub>P<sub>k</sub><sup>2</sup>). A band state concentrates its weight on a few k-points (small value); a state localized in real space, such as a ligand or dangling-bond orbital, spreads it evenly over k (close to 1).</p>
                        <p><span style="color:#009E73">&#9664;</span> <strong>on-band weight</strong>: the share of the state's weight along the path that sits where a bulk band of its own side (valence for occupied, conduction for empty states, aligned on the semicore level) lies within 0.2 eV of its energy. States that follow the bulk bands score high; states in the bulk gap score low. Strongly confined conduction states of small dots also score low, because confinement lifts them above the bulk band.</p>
                        <p><strong>Band edges:</strong> a state is a trap when its k-participation is clearly above that of the most band-like state of the dot (by 0.3, and at least 0.5; very small clusters are broad in k even when delocalized), or when its on-band weight is below 0.25 while it lies more than 0.1 eV inside the aligned bulk gap. Counting inward from the gap, the first state that is neither is the <em>delocalized HOMO</em> (solid blue line) or <em>delocalized LUMO</em> (solid red line); the dotted lines are the nominal DFT HOMO and LUMO, and the states between the two (red-outlined markers) are the traps.</p>
                    </div>
                    <div class="glossary-item trap-indicator">
                        <strong>4. Surface vs. Core Character</strong>
                        <p>Evaluates whether the electron density of a specific state resides in the inner core of the nanostructure or on the outer 25% of the radius. <em>Note: A state with both a high IPR and >90% Surface Character is definitively an unpassivated surface trap.</em></p>
                    </div>
                    <div class="glossary-item">
                        <strong>5. COOP (Crystal Orbital Overlap Population)</strong>
                        <p>Quantifies the chemical bonding interactions between specific pairs of atoms. <strong>Positive</strong> values (plotted to the right) indicate stabilizing bonding interactions, while <strong>negative</strong> values (plotted to the left) indicate destabilizing anti-bonding interactions.</p>
                    </div>
                    <div class="glossary-item cube-indicator">
                        <strong>6. 3D Molecular Orbitals (MOs)</strong>
                        <p>Interactive 3D representations of the electron density limits for the key frontier orbitals (HOMO-1, HOMO, LUMO, LUMO+1). Red and blue lobes represent positive and negative phase domains of the electron wavefunctions.</p>
                    </div>
                </div>
            </div>
        </div>
        
        <script>
        function updateBulk() {{
            var plotDiv = document.getElementById('fuzzy_2d_plot');
            if (!plotDiv || !plotDiv.data) return;
            var show = document.getElementById('bulkShow').checked;
            var hex = document.getElementById('bulkColor').value;
            var width = parseFloat(document.getElementById('bulkWidth').value);
            var alpha = parseFloat(document.getElementById('bulkAlpha').value);
            document.getElementById('bulkWidthVal').textContent = width.toFixed(1);
            document.getElementById('bulkAlphaVal').textContent = alpha.toFixed(2);
            var r = parseInt(hex.substr(1, 2), 16), g = parseInt(hex.substr(3, 2), 16), b = parseInt(hex.substr(5, 2), 16);
            var rgba = function(a) {{ return 'rgba(' + r + ',' + g + ',' + b + ',' + a.toFixed(3) + ')'; }};
            var lines = [], markers = [];
            for (var i = 0; i < plotDiv.data.length; i++) {{
                if (plotDiv.data[i].legendgroup !== 'bulk_bands') continue;
                (plotDiv.data[i].mode === 'markers' ? markers : lines).push(i);
            }}
            if (lines.length)
                Plotly.restyle(plotDiv, {{'visible': show, 'line.color': rgba(alpha), 'line.width': width}}, lines);
            // unfolded bands: opacity and size follow the unfolding weight (customdata)
            markers.forEach(function(i) {{
                var cd = plotDiv.data[i].customdata;
                if (!cd || cd.length === undefined)   // base64-encoded array: use Plotly's decoded copy
                    cd = (plotDiv._fullData.find(function(t) {{ return t.index === i; }}) || {{}}).customdata;
                var w = Array.from(cd || []);
                Plotly.restyle(plotDiv, {{
                    'visible': show,
                    'marker.color': [w.map(function(c) {{ return rgba(Math.min(1, alpha * (0.15 + 0.85 * c) / 0.85)); }})],
                    'marker.size': [w.map(function(c) {{ return (width / 2) * (3.5 + 3.0 * c); }})]
                }}, [i]);
            }});
        }}

        function updateIso() {{
            var val = parseFloat(document.getElementById('isoInput').value);
            if (isNaN(val) || val <= 0) return;
            var plotDiv = document.getElementById('mo_3d_plot');
            if (!plotDiv) return;
            
            for (var i = 0; i < plotDiv.data.length; i++) {{
                var trace = plotDiv.data[i];
                if (trace.type === 'isosurface') {{
                    if (trace.isomin > 0) {{
                        Plotly.restyle(plotDiv, {{'isomin': val, 'isomax': val}}, [i]);
                    }} else if (trace.isomin < 0) {{
                        Plotly.restyle(plotDiv, {{'isomin': -val, 'isomax': -val}}, [i]);
                    }}
                }}
            }}
        }}
        </script>
    </body>
    </html>
    """
    
    out_html = output_html or f"fuzzy_dashboard_{prefix}.html"
    with open(out_html, 'w', encoding='utf-8') as f:
        f.write(html_template)
        
    logger.info(f"  [Plotter] Successfully saved elegant HTML dashboard to {out_html}")
