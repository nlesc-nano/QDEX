import qdex.hardness as h


def test_mnok_block_reports_exponents_and_onsite():
    try:
        h.set_mnok_options(1.0, 2.0, "eta", verbose=False)
        txt = h.format_integrals_block("mnok", "mulliken", "resta", ["Cd", "Se"])
        assert "beta (direct)  : 1 (Mataga-Nishimoto)" in txt
        assert "beta (exchange): 2 (Ohno-Klopman)" in txt
        assert "gamma_AA = eta_A" in txt and "Cd 3.50 eV" in txt
        assert "bulk Resta" in txt
    finally:
        h.set_mnok_options(verbose=False)
    txt = h.format_integrals_block("mnok", "lowdin", "qp", ["Cd"], include_exchange=False)
    assert "Cd 7.00 eV" in txt and "K^x  : off" in txt


def test_stda_block_reports_grimme_parameters():
    _, _, info = h.build_stda_gammas(["Cd", "Se"], [[0, 0, 0], [0, 0, 2.6]], 0.25)
    txt = h.format_integrals_block("mnok", "lowdin", "stda", ["Cd", "Se"], stda_info=info, stda_ax_source="pbe0")
    assert "a_x            : 0.250 (pbe0)" in txt
    assert f"{0.20 + 1.83 * 0.25:.3f}" in txt and f"{1.42 + 0.48 * 0.25:.3f}" in txt
