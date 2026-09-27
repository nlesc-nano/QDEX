"""Console verbosity and the run log of qdex.logging_setup."""
import contextlib
import io
import logging

import qdex  # noqa: F401  (installs the default handler)
from qdex.logging_setup import attach_run_log, detach_run_log, set_verbosity

log = logging.getLogger("qdex.test")


def _emit():
    log.debug("detail")
    log.info("result")
    log.warning("careful")


def test_default_output_is_plain_stdout():
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        _emit()
    assert buf.getvalue() == "detail\nresult\ncareful\n"


def test_verbosity_levels():
    try:
        for verbosity, expected in (("normal", "result\ncareful\n"), ("quiet", "careful\n"), ("full", "detail\nresult\ncareful\n")):
            set_verbosity(verbosity)
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                _emit()
            assert buf.getvalue() == expected, verbosity
    finally:
        set_verbosity("full")


def test_run_log_gets_everything():
    console, logfile = io.StringIO(), io.StringIO()
    try:
        attach_run_log(logfile, console, "quiet")
        _emit()
    finally:
        detach_run_log()
        set_verbosity("full")
    assert console.getvalue() == "careful\n"
    assert logfile.getvalue() == "detail\nresult\ncareful\n"
