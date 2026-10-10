from types import SimpleNamespace

from qdex.cli import _bse_soc_window_ev, _soc_window_ev


def args(**kw):
    base = dict(ewin=None, soc_window=None, soc_bse_window=None)
    base.update(kw)
    return SimpleNamespace(**base)


def test_defaults_fuzzy_one_ev_margin_bse_two():
    a = args()
    assert _soc_window_ev(a) == 6.0
    assert _bse_soc_window_ev(a) == 7.0


def test_defaults_follow_ewin():
    a = args(ewin=[-3.0, 4.0])
    assert _soc_window_ev(a) == 5.0
    assert _bse_soc_window_ev(a) == 6.0


def test_bse_window_follows_explicit_soc_window():
    a = args(soc_window=8.0)
    assert _soc_window_ev(a) == 8.0
    assert _bse_soc_window_ev(a) == 8.0


def test_explicit_bse_window_wins():
    a = args(soc_window=6.0, soc_bse_window=0.0)
    assert _soc_window_ev(a) == 6.0
    assert _bse_soc_window_ev(a) == 0.0
