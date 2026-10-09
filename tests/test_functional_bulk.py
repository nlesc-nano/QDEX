"""Bulk references of the functional of the orbitals (system.functional)."""
import pytest

from qdex import bulk_bands, hardness


@pytest.fixture
def bulk_dir(tmp_path, monkeypatch):
    for name in ("CdSe_zb.bs", "CdSe_zb_soc.bs", "CdSe_zb_hle17.bs", "CdSe_zb_hle17_soc.bs", "CdSe_wz.bs"):
        (tmp_path / name).write_text("")
    monkeypatch.setattr(bulk_bands, "DEFAULT_BULK_DIR", str(tmp_path))
    yield tmp_path
    bulk_bands.set_bulk_functional("pbe")


def test_band_files_follow_functional(bulk_dir):
    name = lambda p: p and p.split("/")[-1]
    assert name(bulk_bands.find_bulk_bs("CdSe", cif="CdSe_zb.cif")) == "CdSe_zb.bs"
    assert name(bulk_bands.find_bulk_bs("CdSe")) == "CdSe_zb.bs"            # PBE never picks a tagged file
    bulk_bands.set_bulk_functional("hle17")
    assert name(bulk_bands.find_bulk_bs("CdSe", cif="CdSe_zb.cif")) == "CdSe_zb_hle17.bs"
    assert name(bulk_bands.find_bulk_bs("CdSe", soc=True)) == "CdSe_zb_hle17_soc.bs"
    # no wurtzite HLE17 bands: the material's HLE17 bands (zinc blende), as for PBE; never a PBE file
    assert name(bulk_bands.find_bulk_bs("CdSe", cif="CdSe_wz.cif")) == "CdSe_zb_hle17.bs"
    assert bulk_bands.find_bulk_bs("CdTe") is None
    assert bulk_bands.functional_label() == "HLE17"


def test_qp_bulk_gaps_follow_functional(monkeypatch):
    monkeypatch.setitem(hardness.MATERIAL_DB_BULK_DFT, "hle17",
                        {"CDSE": dict(a_exp=6.05, gap_exp_lattice=1.0, gap_soc_exp_lattice=0.88, a_pbe=6.10,
                                      gap_pbe_lattice=0.9, gap_soc_pbe_lattice=0.78, qsgw_soc_exp_lattice=2.164,
                                      delta_sigma=1.284)})
    try:
        assert hardness.bulk_gw_gap_sf("CDSE") == pytest.approx(float(hardness.MATERIAL_DB["CDSE"][8]))
        hardness.set_dft_functional("hle17", verbose=False)
        assert hardness.bulk_dft_gap_exp("CDSE") == pytest.approx(1.0)
        assert hardness.bulk_gw_gap_sf("CDSE") == pytest.approx(2.284)
        with pytest.raises(ValueError):
            hardness.bulk_dft_gap_exp("CDTE")
        with pytest.raises(ValueError):
            hardness.set_dft_functional("b3lyp")
    finally:
        hardness.set_dft_functional("pbe", verbose=False)
