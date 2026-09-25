"""Generate the schematic SVG figures used in the QDEX documentation.

Usage (from the repository root)::

    python docs/figures/make_figures.py

Figures are written to ``docs/_static/figures``.  Schematics are drawn with
matplotlib so that they can be regenerated and edited in one place.  The
data-driven panels (anchor QP model, Resta screening profile) evaluate the
same constants and formulas as ``qdex.hardness``; the CdSe energy ladder uses
the numbers of the ``tests/CdSe`` benchmark run (see
``docs/validation/cdse_benchmark.rst``).
"""
import os
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Ellipse

OUT = Path(__file__).resolve().parents[1] / "_static" / "figures"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "mathtext.fontset": "dejavusans",
    "svg.fonttype": "path",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.08,
})

# Palette: occupied/hole = blue, virtual/electron = red, QP/model = violet,
# kernels = green, dynamics = amber, neutral = grey.
C_OCC = "#2563eb"
C_VIR = "#dc2626"
C_QP = "#7c3aed"
C_KER = "#059669"
C_DYN = "#d97706"
C_GREY = "#4b5563"
C_LIGHT = "#f3f4f6"
C_EDGE = "#9ca3af"
FILL = {
    "in": "#eef2ff", "gs": "#ecfeff", "soc": "#fdf4ff", "qp": "#f5f3ff",
    "ker": "#ecfdf5", "bse": "#fff7ed", "out": "#f9fafb", "dyn": "#fffbeb",
}


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.svg", facecolor="white")
    preview = os.environ.get("QDEX_FIG_PREVIEW")
    if preview:  # optional raster previews for reviewing layouts
        fig.savefig(Path(preview) / f"{name}.png", facecolor="white", dpi=90)
    plt.close(fig)
    print("wrote", OUT / f"{name}.svg")


def box(ax, x, y, w, h, text, fc="white", ec=C_GREY, fs=9, bold_first=True, lw=1.2, ha="center"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                                fc=fc, ec=ec, lw=lw))
    lines = text.split("\n")
    if bold_first and len(lines) > 1:
        ax.text(x + w / 2 if ha == "center" else x + 0.1, y + h - 0.2, lines[0], ha=ha, va="top",
                fontsize=fs + 0.5, fontweight="bold", color="#111827")
        ax.text(x + w / 2 if ha == "center" else x + 0.1, y + h - 0.72, "\n".join(lines[1:]), ha=ha,
                va="top", fontsize=fs - 0.5, color="#374151", linespacing=1.35)
    else:
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
                color="#111827", fontweight="bold" if bold_first else "normal")


def arrow(ax, p0, p1, color=C_GREY, lw=1.4, style="-|>", ls="-", rad=0.0, ms=12):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=ms, color=color, lw=lw,
                                 linestyle=ls, connectionstyle=f"arc3,rad={rad}"))


def canvas(w, h, xlim, ylim):
    fig, ax = plt.subplots(figsize=(w, h))
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.axis("off")
    return fig, ax


# ---------------------------------------------------------------------------
# 1. End-to-end pipeline
# ---------------------------------------------------------------------------
def fig_pipeline():
    fig, ax = canvas(11.5, 5.8, (0, 23), (0.9, 12.4))
    # inputs
    box(ax, 0.2, 9.6, 4.2, 2.5, "DFT inputs (CP2K)\nMO coefficients C, ε_p\ngeometry (.xyz)\nGaussian basis set\nGTH-SOC potentials (opt.)", FILL["in"], C_OCC)
    box(ax, 5.2, 9.6, 5.0, 2.5, "1 · Ground-state analysis\nAO overlap S, Löwdin S^½\nMulliken/Löwdin populations\nPDOS, COOP, IPR, fuzzy bands\norbital cubes", FILL["gs"], "#0891b2")
    box(ax, 11.0, 9.6, 5.3, 2.5, "2 · Spin–orbit coupling (opt.)\nGTH projectors → H_SOC\nactive window of MOs\n→ two-component spinors\n(soc_flag, soc_window)", FILL["soc"], "#c026d3")
    box(ax, 17.1, 9.6, 5.7, 2.5, "3 · Quasiparticle model\nscalar scissor (pbe, brus, gw)\nΔW charging (sgw-dim/resta/sbse)\neigenvalue loop (evgw-*)\norbital relaxation (qsgw-*)", FILL["qp"], C_QP)
    arrow(ax, (4.45, 10.85), (5.15, 10.85))
    arrow(ax, (10.25, 10.85), (10.95, 10.85))
    arrow(ax, (16.35, 10.85), (17.05, 10.85))

    box(ax, 17.1, 5.6, 5.7, 3.0, "4 · Interaction kernels\nbare v: MNOK γ_AB or xs (μμ|νν)\nscreened W: Resta, DIM/Thole,\nsBSE (RPA-like)\natomic charges: Mulliken / Löwdin", FILL["ker"], C_KER)
    arrow(ax, (19.95, 9.55), (19.95, 8.65))
    box(ax, 9.3, 5.6, 7.0, 3.0, "5 · Excitons: BSE / TDA\nA = ΔE_QP + 2Kˣ(v) − Kᵈ(W)\nframeworks: independent_dft / _qp,\ndiagonal_bse, full bse\nsolvers: dense eigh or Davidson", FILL["bse"], C_VIR)
    arrow(ax, (17.05, 7.1), (16.35, 7.1))

    box(ax, 0.2, 5.6, 8.3, 3.0, "6 · Exciton analysis & spectra\noscillator strengths, absorption spectrum\ntransition density, NTOs\nPlasser–Dreuw descriptors (d_eh, d_CT, σ_h, σ_e)\ncubes and HTML dashboards", FILL["out"], C_GREY)
    arrow(ax, (9.25, 7.1), (8.55, 7.1))

    box(ax, 0.2, 1.2, 7.2, 3.1, "7 · Non-adiabatic dynamics\nMD frames → overlaps, phase &\ncrossing tracking → NACs\nFSSH-EDC · DISH · PME master eq.\ncarrier cooling, populations", FILL["dyn"], C_DYN)
    box(ax, 8.1, 1.2, 7.0, 3.1, "8 · Recombination\nradiative (Einstein A), energy-gap law,\nFCWD, traps, PLQY\nAuger: eeh / hhe / XX, screened M_if\nstatic rates and ECSH hops", FILL["dyn"], C_DYN)
    box(ax, 15.8, 1.2, 7.0, 3.1, "9 · Time-resolved spectroscopy\nstate filling from P_ia(t)\nGSB / SE / ESA, ΔA(E,t)\n1S-bleach rise → τ_C\n2D vibronic maps", FILL["dyn"], C_DYN)
    arrow(ax, (12.8, 5.55), (3.8, 4.35), rad=0.08)
    arrow(ax, (12.8, 5.55), (11.6, 4.35))
    arrow(ax, (12.8, 5.55), (19.3, 4.35), rad=-0.08)
    ax.text(11.5, 12.35, "QDEX workflow: from DFT orbitals to excitons and carrier dynamics",
            ha="center", va="bottom", fontsize=12, fontweight="bold")
    save(fig, "pipeline")


# ---------------------------------------------------------------------------
# 2. QP model hierarchy
# ---------------------------------------------------------------------------
def fig_qp_hierarchy():
    rows = [
        ("pbe", "none", "—", "—", "—", "DFT reference only"),
        ("brus", "bulk exp. gap + EMA confinement", "scalar", "—", "—", "kinetic confinement only"),
        ("gw / sgw-anchor", "bulk ΔGW + κ/(R+ℓ) + anchor term", "scalar", "—", "—", "bulk ↔ monomer interpolation"),
        ("sgw-dim / sgw-resta(-pure)", "bulk ΔGW + ½ qᵀ ΔW q (edges)", "scalar", "HOMO / LUMO", "—", "static charging model"),
        ("sgw  (+ --kernel sbse)", "bulk ΔGW + ½ Σ q_A ΔW_AA", "scalar", "HOMO / LUMO", "—", "site-diagonal sBSE screening"),
        ("evgw-dim / -resta", "iterate gap ↔ screening", "scalar", "HOMO / LUMO", "—", "fixed point in the gap only"),
        ("qsgw-dim / -resta", "H_DFT + H_bulk + Z·ΔCOHSEX", "all ε_p", "all levels", "C_QP", "static, Löwdin AO basis"),
    ]
    fig, ax = canvas(12.5, 4.6, (0, 23), (0, 9.4))
    heads = ["qp_gap option", "correction", "gap", "edge split", "orbitals", "status"]
    xs = [0.2, 6.6, 13.4, 15.4, 18.5, 20.6]
    ws = [5.0, 6.7, 1.9, 2.4, 1.7, 4.4]
    y = 8.6
    for x, w, h in zip(xs, ws, heads):
        ax.text(x + 0.1, y, h, fontsize=9.5, fontweight="bold", va="center")
    ax.plot([0.2, 26.4], [8.2, 8.2], color=C_GREY, lw=1)
    cols = ["#f9fafb", "#f3f4f6"]
    for k, r in enumerate(rows):
        yy = 7.35 - k * 1.05
        ax.add_patch(FancyBboxPatch((0.2, yy - 0.45), 26.2, 0.9, boxstyle="round,pad=0.0,rounding_size=0.05",
                                    fc=cols[k % 2], ec="none"))
        depth = [0, 1, 2, 3, 3, 4, 5][k]
        ax.add_patch(FancyBboxPatch((0.25, yy - 0.4), 0.12, 0.8, boxstyle="square,pad=0", fc=C_QP,
                                    ec="none", alpha=0.25 + 0.13 * depth))
        for j, (x, t) in enumerate(zip(xs, r)):
            ax.text(x + 0.3 if j == 0 else x + 0.1, yy, t, fontsize=8.6, va="center",
                    family="DejaVu Sans Mono" if j == 0 else "DejaVu Sans",
                    color=C_QP if (j == 4 and t != "—") else "#111827")
    arrow(ax, (27.0, 7.6), (27.0, 0.9), color=C_QP, lw=1.6)
    ax.text(27.0, 0.55, "more\nstate\nresolution", fontsize=7.5, ha="center", va="top", color=C_QP)
    ax.set_xlim(0, 27.6)
    ax.set_ylim(-0.6, 9.4)
    save(fig, "qp_hierarchy")


# ---------------------------------------------------------------------------
# 3. Anchor-scaled QP model (data driven, CdSe database entry)
# ---------------------------------------------------------------------------
def fig_anchor_model():
    try:
        from qdex.hardness import MATERIAL_DB
        from qdex.constants import IMAGE_CHARGE_CONST_EV_ANG as KAPPA0
        e = MATERIAL_DB["CDSE"]
    except Exception:  # documentation build without the compiled package
        KAPPA0 = 11.52
        e = (6.2, 56.0, 6.05, 1.74, 9.5, 0.13, 0.026, 0.64, 1.91, 5.2578, -6.4196, -3.7852, -7.8177, -1.7884)
    eps_inf, gpbe, ggw, R0 = e[0], e[7], e[8], e[9]
    anchor = (e[13] - e[12]) - (e[11] - e[10])
    dbulk = ggw - gpbe
    ell, p = 1.0, 2.0
    kvac = KAPPA0 * (1 - 1 / eps_inf)
    A = anchor - dbulk - kvac / (R0 + ell)
    R = np.linspace(R0, 40, 400)

    def scissor(eps_out, R):
        k = KAPPA0 * (1 / eps_out - 1 / eps_inf)
        return dbulk + k / (R + ell) + A * (R0 / R) ** p

    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    ax.axhline(dbulk, color=C_GREY, ls="--", lw=1)
    ax.text(39.5, dbulk - 0.04, f"bulk limit ΔGW = {dbulk:.2f} eV", ha="right", va="top", fontsize=8.5, color=C_GREY)
    for eps_out, c, lab in [(1.0, C_QP, "vacuum, ε_out = 1"), (2.4, C_KER, "toluene-like, ε_out = 2.4"),
                            (eps_inf, C_DYN, f"matched, ε_out = ε_∞ = {eps_inf}")]:
        ax.plot(R, scissor(eps_out, R), color=c, lw=2, label=lab)
    Rr = np.linspace(2.5, R0, 50)
    ax.axvspan(2.5, R0, color=C_LIGHT)
    ax.text((2.5 + R0) / 2, 3.9, "R < R₀:\nclamped to R₀\n(no extrapolation)", ha="center", fontsize=7.5, color=C_GREY, va="top")
    ax.plot([R0], [anchor], "o", color="black", ms=6, zorder=5)
    ax.annotate(f"finite GW anchor\n(R₀ = {R0:.2f} Å, Δ₀ = {anchor:.2f} eV)", (R0, anchor), (R0 + 5, anchor - 0.15),
                fontsize=8, arrowprops=dict(arrowstyle="-", color="black", lw=0.8))
    Rq = 9.221
    ax.plot([Rq], [scissor(1.0, Rq)], "s", color=C_QP, ms=6, zorder=5)
    ax.annotate(f"CdSe 2 nm test (R_eff = {Rq:.2f} Å)\nscissor = {scissor(1.0, Rq):.2f} eV", (Rq, scissor(1.0, Rq)),
                (Rq + 6, scissor(1.0, Rq) + 0.45), fontsize=8, arrowprops=dict(arrowstyle="-", color=C_QP, lw=0.8))
    ax.set_xlim(2.5, 40)
    ax.set_ylim(1.0, 4.0)
    ax.set_xlabel("effective radius R (Å)")
    ax.set_ylabel("PBE → QP gap correction Δ(R) (eV)")
    ax.set_title(r"Anchor-scaled model:  $\Delta(R)=\Delta_b+\kappa_{out}/(R+\ell)+A\,(R_0/R)^p$", fontsize=10)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    save(fig, "anchor_model")


# ---------------------------------------------------------------------------
# 4. Dielectric sphere: self-polarization vs electron-hole polarization
# ---------------------------------------------------------------------------
def fig_dielectric_sphere():
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3))
    for ax, title in zip(axes, ["Charged excitation (QP gap)", "Neutral excitation (optical gap)"]):
        ax.set_xlim(-3.3, 3.3)
        ax.set_ylim(-3.1, 3.0)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.add_patch(Circle((0, 0), 2.0, fc="#fef3c7", ec=C_DYN, lw=1.8))
        ax.text(0, 2.55, title, ha="center", fontsize=10.5, fontweight="bold")
        ax.text(-1.3, -1.25, "ε_in", fontsize=11, color=C_DYN)
        ax.text(2.2, 1.9, "ε_out", fontsize=11, color=C_GREY)
    ax = axes[0]
    ax.add_patch(Circle((0.6, 0.3), 0.16, fc=C_VIR, ec="none"))
    ax.text(0.6, 0.3, "−", color="white", ha="center", va="center", fontsize=10, fontweight="bold")
    for t in np.linspace(0, 2 * np.pi, 20, endpoint=False):
        ax.text(2.0 * np.cos(t), 2.0 * np.sin(t), "−", ha="center", va="center", fontsize=12, fontweight="bold",
                color=C_VIR, alpha=0.3 + 0.7 * max(0, np.cos(t - np.arctan2(0.3, 0.6))) ** 2)
    ax.text(-3.2, 2.2, "induced surface charge\nhas the carrier's sign\n(ε_in > ε_out): repulsive", fontsize=7.5, color=C_VIR, va="top")
    ax.text(0, -2.55, "added carrier polarizes the interface → self-energy Σ_pol ∝ (1/ε_out − 1/ε_in)/R\n"
            "raises IP and lowers EA: QP gap opens by ≈ 2Σ_pol when ε_out < ε_in",
            ha="center", fontsize=8, color="#111827")
    ax = axes[1]
    ax.add_patch(Circle((0.7, 0.35), 0.16, fc=C_VIR, ec="none"))
    ax.text(0.7, 0.35, "−", color="white", ha="center", va="center", fontsize=10, fontweight="bold")
    ax.add_patch(Circle((-0.6, -0.2), 0.16, fc=C_OCC, ec="none"))
    ax.text(-0.6, -0.2, "+", color="white", ha="center", va="center", fontsize=10, fontweight="bold")
    ax.plot([-0.45, 0.55], [-0.14, 0.29], color=C_KER, lw=1.5, ls=":")
    ax.text(0.05, 0.35, "W_eh", color=C_KER, fontsize=9, ha="center")
    ax.text(0, -2.55, "e and h self-energies (+) and the induced e–h attraction (−)\n"
            "largely cancel: optical gap ≈ insensitive to ε_out (Brus; Delerue–Lannoo–Allan)",
            ha="center", fontsize=8, color="#111827")
    save(fig, "dielectric_sphere")


# ---------------------------------------------------------------------------
# 5. Screening profiles used by the kernels (data driven)
# ---------------------------------------------------------------------------
def fig_screening_profile():
    HA, BOHR = 27.211386245988, 0.52917721
    eps = 6.2
    d_nn = 2.596  # Å, median Cd–Se nearest-neighbour distance of the CdSe test cluster
    ks = np.sqrt(eps - 1) / (d_nn / BOHR)  # 1/bohr
    eta = {"Cd": 3.4994, "Se": 5.4795}
    a = 0.5 * (HA / eta["Cd"] + HA / eta["Se"])  # bohr
    r = np.linspace(0.0, 15, 600)  # Å
    rb = r / BOHR
    v = HA / np.maximum(rb, 1e-9)
    gam = HA / np.sqrt(rb ** 2 + a ** 2)
    S = 1 / eps + (1 - 1 / eps) * np.exp(-ks * rb)
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 3.9))
    ax = axes[0]
    ax.plot(r, v, color=C_GREY, ls="--", lw=1.3, label="bare 1/r")
    ax.plot(r, gam, color=C_OCC, lw=2, label="MNOK γ_CdSe(r)  (bare, damped)")
    ax.plot(r, S * gam, color=C_KER, lw=2, label="Resta-MNOK W(r) = S(r) γ(r)")
    ax.plot(r, gam / eps, color=C_DYN, lw=1.3, ls=":", label="γ(r)/ε_∞")
    ax.set_ylim(0, 12)
    ax.set_xlim(0, 15)
    ax.set_xlabel("r (Å)")
    ax.set_ylabel("interaction (eV)")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("Bare and screened pair interactions (CdSe)", fontsize=10)
    ax = axes[1]
    ax.plot(r, S, color=C_KER, lw=2)
    ax.axhline(1 / eps, color=C_DYN, ls=":", lw=1.2)
    ax.axhline(1.0, color=C_GREY, ls="--", lw=1)
    lam = 1 / ks * BOHR
    ax.axvline(lam, color=C_GREY, lw=0.8)
    ax.text(lam + 0.2, 0.62, f"λ_s = d_NN/√(ε−1) = {lam:.2f} Å", fontsize=8)
    ax.text(14.5, 1 / eps + 0.03, "1/ε_∞ (long range)", ha="right", fontsize=8, color=C_DYN)
    ax.text(14.5, 0.96, "unscreened (short range)", ha="right", va="top", fontsize=8, color=C_GREY)
    ax.set_ylim(0, 1.08)
    ax.set_xlim(0, 15)
    ax.set_xlabel("r (Å)")
    ax.set_ylabel("S(r) = W/γ")
    ax.set_title(r"$S(r)=\epsilon_\infty^{-1}+(1-\epsilon_\infty^{-1})\,e^{-k_s r}$", fontsize=10)
    fig.tight_layout()
    save(fig, "screening_profile")


# ---------------------------------------------------------------------------
# 6. BSE kernel diagrams
# ---------------------------------------------------------------------------
def fig_bse_kernels():
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    for ax in axes:
        ax.set_xlim(0, 10)
        ax.set_ylim(-0.6, 7.4)
        ax.axis("off")
    ax = axes[0]
    ax.text(5, 7.1, r"Exchange (e–h annihilation):  $+2K^x_{ia,jb}$ (singlet)", ha="center", fontsize=10, fontweight="bold")
    for x0, top, bot in [(1.8, "a", "i"), (8.2, "b", "j")]:
        ax.add_patch(Ellipse((x0, 4.4), 1.6, 2.6, fc="none", ec=C_GREY, lw=1.3))
        arrow(ax, (x0 - 0.8, 4.6), (x0 - 0.8, 4.3), color=C_VIR, ms=10)
        arrow(ax, (x0 + 0.8, 4.2), (x0 + 0.8, 4.5), color=C_OCC, ms=10)
        ax.text(x0 - 1.05, 4.4, top, color=C_VIR, fontsize=11, ha="right", va="center")
        ax.text(x0 + 1.05, 4.4, bot, color=C_OCC, fontsize=11, ha="left", va="center")
    ax.text(1.8, 2.7, r"$\rho_{ia}$", ha="center", fontsize=10)
    ax.text(8.2, 2.7, r"$\rho_{jb}$", ha="center", fontsize=10)
    ax.plot([3.0, 7.0], [4.4, 4.4], color=C_OCC, lw=2.2, ls=(0, (4, 2)))
    ax.text(5, 4.7, "bare v", ha="center", fontsize=9, color=C_OCC)
    ax.text(5, 3.9, r"($\gamma_{AB}$ or $(\mu\mu|\nu\nu)$)", ha="center", fontsize=8.5, color=C_OCC)
    ax.text(5, 1.4, r"$K^x_{ia,jb}=\sum_{AB} q^{ia}_A\,v_{AB}\,q^{jb}_B$", ha="center", fontsize=10)
    ax.text(5, 0.4, "local-field (e–h exchange) term; vanishes for triplets", ha="center", fontsize=8.5, color=C_GREY)
    ax = axes[1]
    ax.text(5, 7.1, r"Direct (e–h attraction):  $-K^d_{ia,jb}$", ha="center", fontsize=10, fontweight="bold")
    arrow(ax, (1.2, 5.6), (8.8, 5.6), color=C_VIR, lw=2, ms=13)
    arrow(ax, (8.8, 3.0), (1.2, 3.0), color=C_OCC, lw=2, ms=13)
    ax.text(1.0, 5.6, "a", color=C_VIR, fontsize=11, ha="right", va="center")
    ax.text(9.0, 5.6, "b", color=C_VIR, fontsize=11, ha="left", va="center")
    ax.text(1.0, 3.0, "i", color=C_OCC, fontsize=11, ha="right", va="center")
    ax.text(9.0, 3.0, "j", color=C_OCC, fontsize=11, ha="left", va="center")
    yy = np.linspace(3.0, 5.6, 80)
    ax.plot(5 + 0.18 * np.sin(np.linspace(0, 12 * np.pi, 80)), yy, color=C_KER, lw=2)
    ax.text(5.4, 4.3, "screened W\n(Resta / DIM / sBSE)", fontsize=8.5, color=C_KER, va="center")
    ax.text(5, 1.4, r"$K^d_{ia,jb}=\sum_{AB} q^{ij}_A\,W_{AB}\,q^{ab}_B$", ha="center", fontsize=10)
    ax.text(5, 0.4, "binds the electron–hole pair; identical for singlets and triplets", ha="center", fontsize=8.5, color=C_GREY)
    fig.text(0.5, -0.02, r"$A_{ia,jb}=(\varepsilon^{QP}_a-\varepsilon^{QP}_i)\,\delta_{ij}\delta_{ab}+2K^x_{ia,jb}-K^d_{ia,jb}$"
             "    (TDA, closed shell; in the spinor basis the exchange term has no factor 2)", ha="center", fontsize=9.5)
    save(fig, "bse_kernels")


# ---------------------------------------------------------------------------
# 7. Energy ladder for the CdSe test
# ---------------------------------------------------------------------------
def fig_energy_ladder():
    cases = [
        ("ε_out = 1 (vacuum)", 1.457, 1.270, 1.134, 3.861, 3.634),
        ("ε_out = 2.4", 1.457, 1.270, 0.476, 3.203, 2.977),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), sharey=True)
    for ax, (lab, dft, bulk, fin, qp, s1) in zip(axes, cases):
        steps = [("DFT (PBE)\nHOMO–LUMO", dft, "#9ca3af"),
                 ("+ bulk\nΔGW", dft + bulk, C_QP),
                 ("+ finite-size\npolarization", qp, C_QP),
                 ("− exciton\nbinding", s1, C_KER)]
        prev = 0
        for k, (name, val, col) in enumerate(steps):
            ax.bar(k, val - prev if k else val, bottom=prev if k else 0, color=col, width=0.62,
                   alpha=0.9 if k else 0.6)
            ax.plot([k - 0.31, k + 0.31], [val, val], color="black", lw=1)
            ax.text(k, val + 0.06, f"{val:.2f} eV", ha="center", fontsize=8.5)
            if k < 3:
                ax.plot([k + 0.31, k + 0.69], [val, val], color=C_GREY, lw=0.8, ls=":")
            prev = val
        ax.set_xticks(range(4))
        ax.set_xticklabels([s[0] for s in steps], fontsize=8)
        ax.set_title(lab, fontsize=10)
        ax.set_ylim(0, 4.4)
        ax.axhspan(2.70, 2.95, color="#fde68a", alpha=0.6, lw=0)
        ax.text(-0.35, 4.35, "shaded: exp. 1S absorption for D ≈ 1.6–2.0 nm\n(Yu et al. 2003 sizing curve, in solution)", fontsize=7.2, va="top", ha="left")
    axes[0].set_ylabel("energy (eV)")
    fig.suptitle("CdSe (Cd₆₈Se₅₅Cl₂₆, R_eff = 9.2 Å): gw model + Resta BSE — binding stays 0.23 eV for any ε_out",
                 fontsize=10)
    fig.tight_layout()
    save(fig, "cdse_energy_ladder")


# ---------------------------------------------------------------------------
# 8. The four excitation frameworks
# ---------------------------------------------------------------------------
def fig_frameworks():
    rng = np.random.default_rng(3)
    n = 12
    d_dft = np.sort(rng.uniform(1.4, 2.6, n))
    d_qp = d_dft + 2.4
    K = rng.normal(0, 0.05, (n, n))
    K = 0.5 * (K + K.T)
    np.fill_diagonal(K, -0.28)
    mats = [("independent_dft", np.diag(d_dft), "ΔE_DFT on the diagonal"),
            ("independent_qp", np.diag(d_qp), "ΔE_QP on the diagonal"),
            ("diagonal_bse", np.diag(d_qp + np.diag(K)), r"$\Delta E^{QP}+2K^x_{ia,ia}-K^d_{ia,ia}$" + "\n(no configuration mixing)"),
            ("bse", np.diag(d_qp) + K, "full A: coupled transitions\n(dense eigh or Davidson)")]
    fig, axes = plt.subplots(1, 4, figsize=(11, 3.3))
    for ax, (name, M, cap) in zip(axes, mats):
        off = M - np.diag(np.diag(M))
        show = np.where(np.eye(n, dtype=bool), 1.0, np.abs(off) / 0.12)
        ax.imshow(show, cmap="Oranges", vmin=0, vmax=1.2)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(True)
            s.set_color(C_EDGE)
        ax.set_title(name, fontsize=9.5, family="DejaVu Sans Mono")
        ax.set_xlabel(cap, fontsize=8)
    fig.suptitle("excitation_mode: what enters the transition-space matrix (rows/columns = i→a transitions)", fontsize=10)
    save(fig, "excitation_frameworks")


# ---------------------------------------------------------------------------
# 9. qsGW orbital-relaxation loop
# ---------------------------------------------------------------------------
def fig_qsgw_loop():
    fig, ax = canvas(11.5, 5.2, (0, 23), (0, 10.4))
    box(ax, 0.3, 7.4, 6.0, 2.6, "Start (Löwdin basis)\nC^L = S^½ C\nH_DFT = C^L diag(ε) C^Lᵀ\nrequires the full square MO set", FILL["gs"], "#0891b2")
    box(ax, 7.4, 7.4, 7.0, 2.6, "Density & screening\nP = 2 C_occ C_occᵀ\nW_QD (DIM or Resta)\nΔW = max(0, W_QD − W_bulk) + ΔW_solv", FILL["ker"], C_KER)
    box(ax, 15.5, 7.4, 7.2, 2.6, "Static ΔCOHSEX operator\nΣ_SEX = −½ P ∘ ΔW_AO\nΣ_COH = ½ diag(ΔW_AO)\nZ from damping model", FILL["qp"], C_QP)
    box(ax, 15.5, 2.8, 7.2, 3.2, "Effective Hamiltonian\nH = H_DFT + H_bulk\n  + Z (Σ_SEX + Σ_COH)\nH_bulk: f_H/f_L split of ΔGW_bulk\nlinear mixing (damping 0.5)", FILL["qp"], C_QP)
    box(ax, 7.4, 2.8, 7.0, 3.2, "Diagonalize\nH C_new = C_new E\nphase alignment of columns\ngap change < tol ?", FILL["bse"], C_VIR)
    box(ax, 0.3, 2.8, 6.0, 3.2, "Output\nε_QP (all levels), C_QP = S^-½ C^L\npassed to BSE with scissor = 0\nprovenance: Z, f_H, f_L,\nconverged flag, fidelities", FILL["out"], C_GREY)
    arrow(ax, (6.35, 8.7), (7.35, 8.7))
    arrow(ax, (14.45, 8.7), (15.45, 8.7))
    arrow(ax, (19.1, 7.35), (19.1, 6.05))
    arrow(ax, (15.45, 4.4), (14.45, 4.4))
    arrow(ax, (7.35, 4.4), (6.35, 4.4))
    ax.text(6.85, 4.65, "yes", fontsize=8.5, ha="center")
    arrow(ax, (10.9, 6.05), (10.9, 7.35), color=C_VIR, ls="--")
    ax.text(11.1, 6.7, "no: update C, gap", fontsize=8.5, color=C_VIR)
    ax.text(11.5, 0.9, "qsgw-* is a static COHSEX-like orbital-relaxation model on atom-block ΔW,\n"
            "not the quasiparticle self-consistent GW of Kotani–van Schilfgaarde–Faleev.",
            ha="center", fontsize=8.5, color=C_GREY, style="italic")
    save(fig, "qsgw_loop")


# ---------------------------------------------------------------------------
# 10. SOC spinor construction
# ---------------------------------------------------------------------------
def fig_soc():
    fig, ax = canvas(10.5, 4.0, (0, 21), (0, 8))
    # levels left: spatial MOs with window
    for k, e in enumerate(np.linspace(0.8, 3.0, 5)):
        ax.plot([0.8, 2.6], [e, e], color=C_OCC, lw=2)
        ax.text(1.7, e + 0.07, "↑↓", ha="center", fontsize=8, color=C_OCC)
    for e in np.linspace(4.2, 6.4, 5):
        ax.plot([0.8, 2.6], [e, e], color=C_VIR, lw=2)
    ax.add_patch(FancyBboxPatch((0.5, 1.6), 2.4, 3.9, boxstyle="round,pad=0.05", fc="none", ec="#c026d3", ls="--", lw=1.3))
    ax.text(1.7, 7.2, "spatial MOs φ_p\n(DFT, spin-free)", ha="center", fontsize=9, fontweight="bold")
    ax.text(3.05, 3.55, "soc_window", color="#c026d3", fontsize=8, rotation=90, va="center")
    arrow(ax, (3.4, 3.8), (5.0, 3.8))
    box(ax, 5.1, 2.0, 5.2, 3.7, "Spin-orbital basis\n{φ_p α, φ_p β} in window\nH = diag(ε_p) ⊗ 1₂ + H_SOC\nH_SOC from GTH projectors\n⟨φ_p σ| V_SO |φ_q σ'⟩", FILL["soc"], "#c026d3")
    arrow(ax, (10.4, 3.8), (12.0, 3.8))
    box(ax, 12.1, 2.0, 4.0, 3.7, "Diagonalize\nH U = U E\nKramers-paired\nspinors ψ_P\n= Σ U φ_p σ", FILL["soc"], "#c026d3")
    arrow(ax, (16.2, 3.8), (17.3, 3.8))
    for k, e in enumerate(np.linspace(0.8, 3.0, 5)):
        for dx in (0, 0.9):
            ax.plot([17.5 + dx, 18.2 + dx], [e + 0.05 * dx, e + 0.05 * dx], color=C_OCC, lw=2)
    for e in np.linspace(4.2, 6.4, 5):
        for dx in (0, 0.9):
            ax.plot([17.5 + dx, 18.2 + dx], [e - 0.04 * dx, e - 0.04 * dx], color=C_VIR, lw=2)
    ax.text(18.3, 7.2, "spinors ψ_P → spinor BSE\n(Kˣ without the factor 2)", ha="center", fontsize=9, fontweight="bold")
    save(fig, "soc_spinors")


# ---------------------------------------------------------------------------
# 11. NAMD workflow
# ---------------------------------------------------------------------------
def fig_namd():
    fig, ax = canvas(11, 4.6, (0, 22), (0, 9.2))
    b = [
        (0.2, "MD trajectory\nframes R(t_k)\n(ground-state or\nexcited-state MD)", FILL["in"], C_OCC),
        (4.6, "Electronic structure\nper frame (CP2K MOs)\nQP shift + BSE\n(optional excitons)", FILL["bse"], C_VIR),
        (9.0, "Overlaps, tracking\n⟨ψ_p(t)|ψ_q(t+Δt)⟩\nphase alignment,\ncrossing detection", FILL["gs"], "#0891b2"),
        (13.4, "NACs, propagation\nd_pq ≈ (O_pq − O_qp)/2Δt\nelectronic Δt ≪ nuclear Δt\nCPA (classical path)", FILL["ker"], C_KER),
        (17.8, "Hopping / kinetics\nFSSH-EDC\nDISH (decoherence-\ninduced hops) · PME", FILL["qp"], C_QP),
    ]
    for x, t, fc, ec in b:
        box(ax, x, 5.0, 4.0, 3.6, t, fc, ec, fs=8.0)
    for k in range(4):
        arrow(ax, (b[k][0] + 4.05, 6.8), (b[k + 1][0] - 0.05, 6.8))
    outs = ["populations P_p(t)\nenergy relaxation", "carrier cooling\nτ_C, k_C", "transient absorption\nΔA(E,t), 2D maps", "recombination\nradiative / Auger (ECSH)"]
    for k, t in enumerate(outs):
        box(ax, 1.0 + k * 5.3, 0.6, 4.4, 2.2, t, FILL["dyn"], C_DYN, fs=8.8, bold_first=False)
        arrow(ax, (19.8, 4.95), (3.2 + k * 5.3, 2.85), rad=0.05 * (k - 1.5), color=C_DYN)
    save(fig, "namd_workflow")


# ---------------------------------------------------------------------------
# 12. Auger processes
# ---------------------------------------------------------------------------
def _levels(ax, x0, title):
    ax.plot([x0, x0 + 3.0], [3.0, 3.0], color=C_GREY, lw=1)
    ax.plot([x0, x0 + 3.0], [1.0, 1.0], color=C_GREY, lw=1)
    ax.text(x0 + 1.5, 7.0, title, ha="center", fontsize=9.5, fontweight="bold")


def _dot(ax, x, y, kind):
    fc = C_VIR if kind == "e" else "white"
    ax.add_patch(Circle((x, y), 0.17, fc=fc, ec=C_VIR if kind == "e" else C_OCC, lw=1.6, zorder=4))


def fig_auger():
    fig, ax = canvas(11, 4.3, (0, 22), (0, 7.6))
    for x0, title in [(0.6, "eeh (X⁻ / XX)"), (8.1, "hhe (X⁺ / XX)"), (15.6, "biexciton XX")]:
        _levels(ax, x0, title)
        ax.text(x0 - 0.15, 3.0, "CB", ha="right", va="center", fontsize=8, color=C_GREY)
        ax.text(x0 - 0.15, 1.0, "VB", ha="right", va="center", fontsize=8, color=C_GREY)
    # eeh
    x = 0.6
    _dot(ax, x + 0.8, 3.0, "e"); _dot(ax, x + 1.8, 3.0, "e"); _dot(ax, x + 1.8, 1.0, "h")
    arrow(ax, (x + 1.8, 2.8), (x + 1.8, 1.2), color=C_VIR, lw=1.6)
    arrow(ax, (x + 0.8, 3.2), (x + 0.8, 6.0), color=C_VIR, lw=1.6, ls="--")
    _dot(ax, x + 0.8, 6.1, "e")
    ax.text(x + 2.0, 2.0, "e₂ + h\nrecombine", fontsize=7.5, va="center")
    ax.text(x + 1.0, 5.2, "e₁ → e′\n(hot electron)", fontsize=7.5)
    # hhe
    x = 8.1
    _dot(ax, x + 1.8, 3.0, "e"); _dot(ax, x + 0.8, 1.0, "h"); _dot(ax, x + 1.8, 1.0, "h")
    arrow(ax, (x + 1.8, 2.8), (x + 1.8, 1.2), color=C_VIR, lw=1.6)
    arrow(ax, (x + 0.8, 0.8), (x + 0.8, -0.2), color=C_OCC, lw=1.6, ls="--")
    _dot(ax, x + 0.8, -0.3, "h")
    ax.text(x + 2.0, 2.0, "e + h₂\nrecombine", fontsize=7.5, va="center")
    ax.text(x + 1.05, -0.3, "h₁ → h′ (deep hole)", fontsize=7.5, va="center")
    # XX
    x = 15.6
    for dx in (0.7, 1.4):
        _dot(ax, x + dx, 3.0, "e"); _dot(ax, x + dx, 1.0, "h")
    arrow(ax, (x + 1.4, 2.8), (x + 1.4, 1.2), color=C_VIR, lw=1.6)
    arrow(ax, (x + 0.7, 3.2), (x + 0.7, 6.0), color=C_VIR, lw=1.4, ls="--")
    arrow(ax, (x + 2.4, 0.8), (x + 2.4, -0.2), color=C_OCC, lw=1.4, ls=":")
    ax.text(x + 2.55, -0.1, "or h → h′", fontsize=7.5, va="center")
    ax.text(x + 0.85, 6.1, "e → e′", fontsize=7.5, va="center")
    ax.text(x + 1.3, 4.6, "4 e–h pairs × (eeh + hhe):\nk_XX = 2k_X⁻ + 2k_X⁺", fontsize=7.5)
    ax.text(11.0, -1.6, r"$M_{if}=V^{dir}-V^{exch}$,  $V_{e'e_1,he_2}=\sum_{AB} q^{e'e_1}_A W_{AB}\,q^{he_2}_B$;"
            r"   $\Gamma=\frac{2\pi}{\hbar}\sum_f|M_{if}|^2\,\rho(E_f)$", ha="center", fontsize=9)
    ax.set_ylim(-2.1, 7.6)
    save(fig, "auger_processes")


# ---------------------------------------------------------------------------
# 13. Transient absorption contributions
# ---------------------------------------------------------------------------
def fig_ta():
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.8), gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axes[0]
    ax.set_xlim(0, 12)
    ax.set_ylim(-0.6, 7.2)
    ax.axis("off")
    for x0, title, sign in [(0.4, "GSB / state filling", "ΔA < 0"), (4.4, "stimulated emission", "ΔA < 0"),
                            (8.4, "excited-state absorption", "ΔA > 0")]:
        for y, c in [(1.0, C_OCC), (3.2, C_VIR), (5.6, C_VIR)]:
            ax.plot([x0, x0 + 2.8], [y, y], color=c, lw=1.6)
        ax.text(x0 + 1.4, 6.7, title, ha="center", fontsize=9, fontweight="bold")
        ax.text(x0 + 1.4, -0.4, sign, ha="center", fontsize=9, color=C_KER if ">" in sign else C_VIR)
    _dot(ax, 1.3, 3.2, "e"); _dot(ax, 1.3, 1.0, "h")
    ax.plot([2.3, 2.3], [1.2, 3.0], color=C_GREY, lw=1.4, ls=":")
    ax.text(2.45, 2.1, "blocked\n(Pauli)", fontsize=7.5, va="center")
    _dot(ax, 5.8, 3.2, "e"); _dot(ax, 5.8, 1.0, "h")
    arrow(ax, (5.8, 3.0), (5.8, 1.2), color=C_VIR)
    ax.text(6.0, 2.1, "probe-\nstimulated", fontsize=7.5, va="center")
    _dot(ax, 9.3, 3.2, "e")
    arrow(ax, (9.3, 3.4), (9.3, 5.4), color=C_KER)
    ax.text(9.5, 4.3, "probe\nabsorbed", fontsize=7.5, va="center")
    ax = axes[1]
    E = np.linspace(2.4, 3.6, 400)
    t = np.array([0.05, 0.3, 1.0, 3.0])
    for k, tt in enumerate(t):
        w = 1 - np.exp(-tt / 0.6)
        g = -w * np.exp(-((E - 2.62) / 0.05) ** 2) - (1 - w) * 0.6 * np.exp(-((E - 3.2) / 0.12) ** 2) \
            + 0.15 * w * np.exp(-((E - 2.5) / 0.04) ** 2)
        ax.plot(E, g, color=plt.cm.viridis(k / 3.5), lw=1.8, label=f"t = {tt:g} ps")
    ax.axhline(0, color=C_GREY, lw=0.8)
    ax.set_xlabel("probe energy (eV)")
    ax.set_ylabel("ΔA (arb.)")
    ax.set_title("1S bleach rises as hot carriers cool", fontsize=9.5)
    ax.legend(frameon=False, fontsize=7.5)
    ax.set_yticks([])
    fig.tight_layout()
    save(fig, "ta_contributions")


# ---------------------------------------------------------------------------
# 14. Environment cancellation: QP vs optical energies vs eps_out
# ---------------------------------------------------------------------------
def fig_environment_cancellation():
    # CdSe 2 nm runs (tests/CdSe, spin-free, 25x25, Resta interior kernel)
    eps_inf, R, ell = 6.2, 9.221, 1.0
    eo = np.linspace(1.0, 6.2, 100)
    dqp = 11.52 * (1 / eo - 1 / eps_inf) / (R + ell)
    off = 11.52 * (1 / 2.4 - 1 / eps_inf) / (R + ell)
    env_eo = [1.0, 1.9, 2.24, 6.2]
    env_qp = [3.852, 3.223, 3.115, 2.727]
    env_s1 = [2.501, 2.503, 2.503, 2.500]
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.9), sharey=True)
    ax = axes[0]
    ax.plot(eo, 3.203 - off + dqp, color=C_QP, lw=2, label="QP gap")
    ax.plot(eo, 2.977 - off + dqp, color=C_VIR, lw=2, label="S₁")
    ax.plot([1.0, 2.4], [3.861, 3.203], "o", color=C_QP)
    ax.plot([1.0, 2.4], [3.634, 2.977], "o", color=C_VIR)
    ax.set_title("scissor model (qp_gap: gw): S₁ follows the QP gap", fontsize=9.5)
    ax = axes[1]
    ax.plot(env_eo, env_qp, "o-", color=C_QP, lw=2, label="QP gap")
    ax.plot(env_eo, env_s1, "o-", color=C_VIR, lw=2, label="S₁")
    ax.set_title("environment model (qp_gap: env): S₁ constant", fontsize=9.5)
    for ax in axes:
        ax.axhspan(2.70, 2.95, color="#fde68a", alpha=0.6, lw=0)
        ax.set_xlabel("ε_out")
        ax.legend(frameon=False, fontsize=8, loc="upper right")
    axes[0].set_ylabel("energy (eV)")
    axes[0].text(6.1, 2.72, "exp. 1S window", ha="right", fontsize=7.5)
    fig.tight_layout()
    save(fig, "environment_cancellation")


if __name__ == "__main__":
    fig_pipeline()
    fig_qp_hierarchy()
    fig_anchor_model()
    fig_dielectric_sphere()
    fig_screening_profile()
    fig_bse_kernels()
    fig_energy_ladder()
    fig_frameworks()
    fig_qsgw_loop()
    fig_soc()
    fig_namd()
    fig_auger()
    fig_ta()
    fig_environment_cancellation()
