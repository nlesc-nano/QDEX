"""Output of QDEX through the standard ``logging`` module.

Every module logs to ``logging.getLogger(__name__)``, below the ``qdex`` logger. Messages are written
as they are (no prefix) to standard output, so the output reads as before. The verbosity
(``output.verbosity``, ``--verbosity``) selects what reaches the console:

* ``full`` (default): everything;
* ``normal``: no detail lines (iteration traces, timings, diagnostics; level DEBUG);
* ``quiet``: warnings and errors only.

The run log file (``system.log_file``) always receives everything.
"""
import logging
import sys

VERBOSITY = {"full": logging.DEBUG, "normal": logging.INFO, "quiet": logging.WARNING}
_ROOT = "qdex"


class _StdoutHandler(logging.StreamHandler):
    """Writes each record to ``stream``, or to the current ``sys.stdout`` when none is given.

    Resolving ``sys.stdout`` at emit time keeps ``contextlib.redirect_stdout`` and pytest capture working.
    """

    def __init__(self, stream=None):
        self._fixed = stream
        super().__init__(stream if stream is not None else sys.stdout)

    @property
    def stream(self):
        return self._fixed if self._fixed is not None else sys.stdout

    @stream.setter
    def stream(self, value):
        pass


def _handler(stream=None, level=logging.DEBUG):
    h = _StdoutHandler(stream)
    h.setFormatter(logging.Formatter("%(message)s"))
    h.setLevel(level)
    return h


def _logger():
    return logging.getLogger(_ROOT)


def install_default_handler():
    """Console output to the current sys.stdout at full verbosity (done when qdex is imported)."""
    lg = _logger()
    if getattr(lg, "_qdex_console", None) is not None:
        return
    lg.setLevel(logging.DEBUG)
    lg.propagate = False
    lg._qdex_console = _handler()
    lg._qdex_file = None
    lg.addHandler(lg._qdex_console)


def _level(verbosity):
    key = str(verbosity or "full").lower()
    if key not in VERBOSITY:
        raise ValueError(f"verbosity must be one of {', '.join(VERBOSITY)}, not '{verbosity}'")
    return VERBOSITY[key]


def set_verbosity(verbosity="full"):
    """Console verbosity: 'full', 'normal' or 'quiet'."""
    install_default_handler()
    _logger()._qdex_console.setLevel(_level(verbosity))


def attach_run_log(log_stream, console_stream, verbosity="full"):
    """Log everything to ``log_stream`` and the chosen verbosity to ``console_stream``."""
    install_default_handler()
    lg = _logger()
    detach_run_log()
    lg.removeHandler(lg._qdex_console)
    lg._qdex_console = _handler(console_stream, _level(verbosity))
    lg._qdex_file = _handler(log_stream, logging.DEBUG)
    lg.addHandler(lg._qdex_console)
    lg.addHandler(lg._qdex_file)


def detach_run_log():
    """Back to console output on the current sys.stdout (keeps the verbosity)."""
    lg = _logger()
    if getattr(lg, "_qdex_file", None) is None:
        return
    level = lg._qdex_console.level
    lg.removeHandler(lg._qdex_file)
    lg.removeHandler(lg._qdex_console)
    lg._qdex_file = None
    lg._qdex_console = _handler(level=level)
    lg.addHandler(lg._qdex_console)
