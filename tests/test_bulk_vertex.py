"""Vertex correction of the bulk QSGW shift (quasiparticles.bulk_vertex)."""
import pytest

import qdex.hardness as h


@pytest.fixture(autouse=True)
def _reset():
    yield
    h.set_bulk_vertex("none", 0.8, verbose=False)


def test_none_is_pure_qsgw():
    h.set_bulk_vertex("none", verbose=False)
    shift, info = h.bulk_qp_shift("CDSE", 1.457)
    assert shift == pytest.approx(2.19 - 0.62)
    assert info["bulk_vertex_fraction"] == 0.0


def test_full_applies_the_bulk_factor():
    h.set_bulk_vertex("full", 0.8, verbose=False)
    shift, info = h.bulk_qp_shift("CDSE", 2.639)
    assert shift == pytest.approx(0.8 * (2.19 - 0.62))
    assert info["bulk_vertex_correction_ev"] == pytest.approx(-0.2 * 1.57)


def test_scaled_follows_the_penn_screening_fraction():
    h.set_bulk_vertex("scaled", 0.8, verbose=False)
    d = 2.19 - 0.62
    bulk, _ = h.bulk_qp_shift("CDSE", 0.62)             # no confinement: full bulk correction
    assert bulk == pytest.approx(0.8 * d)
    shifts = [h.bulk_qp_shift("CDSE", g)[0] for g in (0.62, 1.0, 1.457, 2.639, 6.0)]
    assert all(a <= b for a, b in zip(shifts, shifts[1:]))   # smaller dots keep more of the QSGW shift
    s12, info = h.bulk_qp_shift("CDSE", 2.639)
    f = (info["bulk_vertex_eps_eff"] - 1.0) / (6.2 - 1.0)
    assert info["bulk_vertex_fraction"] == pytest.approx(f)
    assert s12 == pytest.approx(d - 0.2 * d * f)


def test_invalid_settings():
    with pytest.raises(ValueError):
        h.set_bulk_vertex("partial", verbose=False)
    with pytest.raises(ValueError):
        h.set_bulk_vertex("full", 1.5, verbose=False)
