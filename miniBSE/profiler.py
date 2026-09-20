"""
Computational Resource Usage Profiler for miniBSE.
Tracks wall-clock time, memory (RAM), and GPU VRAM across calculation modules.
"""

import os
import sys
import time
import platform
import resource
from contextlib import contextmanager

# -----------------------------------------------------------------------------
# Host RAM Inspection (Zero External Dependencies)
# -----------------------------------------------------------------------------

def get_memory_info_mb():
    """
    Returns (current_rss_mb, peak_rss_mb) for the current process.
    Supports macOS (via Mach task_info) and Linux (via /proc/self/status & getrusage).
    """
    current_rss = 0.0
    peak_rss = 0.0

    if sys.platform == "darwin":
        try:
            import ctypes
            class mach_task_basic_info(ctypes.Structure):
                _fields_ = [
                    ('virtual_size', ctypes.c_uint64),
                    ('resident_size', ctypes.c_uint64),
                    ('resident_size_max', ctypes.c_uint64),
                    ('user_time_sec', ctypes.c_uint32),
                    ('user_time_usec', ctypes.c_uint32),
                    ('system_time_sec', ctypes.c_uint32),
                    ('system_time_usec', ctypes.c_uint32),
                    ('policy', ctypes.c_int32),
                    ('suspend_count', ctypes.c_int32)
                ]
            libc = ctypes.CDLL(None)
            mach_task_self = libc.mach_task_self
            mach_task_self.restype = ctypes.c_uint32
            info = mach_task_basic_info()
            count = ctypes.c_uint32(ctypes.sizeof(info) // 4)
            # MACH_TASK_BASIC_INFO = 20
            ret = libc.task_info(mach_task_self(), 20, ctypes.byref(info), ctypes.byref(count))
            if ret == 0:
                current_rss = info.resident_size / (1024.0 * 1024.0)
                peak_rss = info.resident_size_max / (1024.0 * 1024.0)
        except Exception:
            pass

    elif sys.platform.startswith("linux"):
        try:
            with open("/proc/self/status", "r") as f:
                for line in f:
                    if line.startswith("VmRSS:"):
                        current_rss = float(line.split()[1]) / 1024.0
                    elif line.startswith("VmHWM:"):
                        peak_rss = float(line.split()[1]) / 1024.0
        except Exception:
            pass

    # Fallback for peak RSS using getrusage if peak was not resolved
    if peak_rss <= 0.0:
        ru = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if sys.platform == "darwin":
            peak_rss = ru / (1024.0 * 1024.0) # bytes to MB on macOS
        else:
            peak_rss = ru / 1024.0 # KB to MB on Linux/others

    if current_rss <= 0.0:
        current_rss = peak_rss

    return current_rss, peak_rss


# -----------------------------------------------------------------------------
# GPU VRAM Inspection
# -----------------------------------------------------------------------------

def get_gpu_memory_info_mb():
    """
    Returns (current_vram_mb, peak_vram_mb) if PyTorch CUDA is initialized, else None.
    """
    try:
        import torch
        if torch.cuda.is_available() and torch.cuda.is_initialized():
            curr = torch.cuda.memory_allocated() / (1024.0 * 1024.0)
            peak = torch.cuda.max_memory_allocated() / (1024.0 * 1024.0)
            return curr, peak
    except Exception:
        pass
    return None


# -----------------------------------------------------------------------------
# Resource Tracker
# -----------------------------------------------------------------------------

class ResourceTracker:
    def __init__(self):
        self.stages = []
        self._active_stage = None
        self._stage_start_time = None
        self._stage_start_ram = None
        self._stage_start_gpu = None
        self.t0_total = time.time()

    def start_stage(self, name):
        """Starts recording metrics for a stage."""
        if self._active_stage is not None:
            self.end_stage()

        self._active_stage = name
        self._stage_start_time = time.time()
        ram_curr, _ = get_memory_info_mb()
        self._stage_start_ram = ram_curr
        gpu_info = get_gpu_memory_info_mb()
        self._stage_start_gpu = gpu_info[0] if gpu_info is not None else None

    def end_stage(self, name=None):
        """Ends recording metrics for the active stage."""
        if self._active_stage is None:
            return

        end_time = time.time()
        elapsed = end_time - self._stage_start_time
        ram_curr, ram_peak = get_memory_info_mb()
        delta_ram = ram_curr - self._stage_start_ram

        gpu_info = get_gpu_memory_info_mb()
        if gpu_info is not None and self._stage_start_gpu is not None:
            gpu_curr, gpu_peak = gpu_info
            delta_gpu = gpu_curr - self._stage_start_gpu
        else:
            gpu_curr = None
            gpu_peak = None
            delta_gpu = None

        stage_name = name if name is not None else self._active_stage
        self.stages.append({
            "name": stage_name,
            "elapsed": elapsed,
            "ram_curr": ram_curr,
            "ram_peak": ram_peak,
            "delta_ram": delta_ram,
            "gpu_curr": gpu_curr,
            "gpu_peak": gpu_peak,
            "delta_gpu": delta_gpu,
        })

        self._active_stage = None
        self._stage_start_time = None
        self._stage_start_ram = None
        self._stage_start_gpu = None

    @contextmanager
    def stage(self, name):
        """Context manager for tracking a calculation stage."""
        self.start_stage(name)
        try:
            yield
        finally:
            self.end_stage(name)

    def print_summary(self, device=None, nthreads=None):
        """
        Prints a comprehensive, formatted ASCII summary table of resource usage.
        """
        if self._active_stage is not None:
            self.end_stage()

        if not self.stages:
            return

        total_elapsed = sum(s["elapsed"] for s in self.stages)
        has_gpu = any(s["gpu_peak"] is not None for s in self.stages)
        _, overall_peak_ram = get_memory_info_mb()

        # System metadata
        os_name = f"{platform.system()} ({platform.machine()})"
        cpu_cores = nthreads if nthreads else os.cpu_count()
        dev_str = str(device) if device else "numpy (CPU)"

        print("\n" + "=" * 90)
        print(" miniBSE COMPUTATIONAL RESOURCE USAGE SUMMARY")
        print("=" * 90)
        print(f" Platform: {os_name} | Threads/Cores: {cpu_cores} | Device: {dev_str}")
        print("-" * 90)

        if has_gpu:
            header = f" {'Module / Calculation Stage':<34} {'Wall Time':>11} {'% Total':>8} {'Peak RAM':>11} {'ΔRAM':>10} {'Peak VRAM':>11}"
            print(header)
            print("-" * 90)
            for s in self.stages:
                pct = (s["elapsed"] / total_elapsed * 100.0) if total_elapsed > 0 else 0.0
                d_ram_sign = "+" if s["delta_ram"] >= 0 else ""
                d_ram_str = f"{d_ram_sign}{s['delta_ram']:.1f} MB"
                vram_str = f"{s['gpu_peak']:.1f} MB" if s["gpu_peak"] is not None else "-"
                print(f" {s['name']:<34} {s['elapsed']:>9.2f} s {pct:>7.1f}% {s['ram_peak']:>9.1f} MB {d_ram_str:>10} {vram_str:>11}")
            print("-" * 90)
            print(f" {'Total Execution Time':<34} {total_elapsed:>9.2f} s {'100.0%':>8}")
            print(f" {'Overall Peak RAM Usage':<34} {overall_peak_ram:>9.1f} MB ({overall_peak_ram / 1024.0:.2f} GB)")
            overall_gpu = max((s["gpu_peak"] for s in self.stages if s["gpu_peak"] is not None), default=None)
            if overall_gpu is not None:
                print(f" {'Overall Peak GPU VRAM':<34} {overall_gpu:>9.1f} MB ({overall_gpu / 1024.0:.2f} GB)")
        else:
            header = f" {'Module / Calculation Stage':<36} {'Wall Time':>11} {'% Total':>8} {'Peak RAM':>11} {'ΔRAM':>10}"
            print(header)
            print("-" * 90)
            for s in self.stages:
                pct = (s["elapsed"] / total_elapsed * 100.0) if total_elapsed > 0 else 0.0
                d_ram_sign = "+" if s["delta_ram"] >= 0 else ""
                d_ram_str = f"{d_ram_sign}{s['delta_ram']:.1f} MB"
                print(f" {s['name']:<36} {s['elapsed']:>9.2f} s {pct:>7.1f}% {s['ram_peak']:>9.1f} MB {d_ram_str:>10}")
            print("-" * 90)
            print(f" {'Total Execution Time':<36} {total_elapsed:>9.2f} s {'100.0%':>8}")
            print(f" {'Overall Peak RAM Usage':<36} {overall_peak_ram:>9.1f} MB ({overall_peak_ram / 1024.0:.2f} GB)")

        print("=" * 90 + "\n")
