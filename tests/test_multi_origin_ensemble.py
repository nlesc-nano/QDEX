"""
Unit tests for the Automated Multi-Origin Ensemble Sampling module in QDEX NAMD.
Tests:
  1. Energy autocorrelation function and phonon dephasing / correlation time extraction.
  2. Auto-calibration of simulation window and origin grid.
  3. Graceful fallback on short MD trajectories.
  4. Energy-referenced initial state sampling preserving excess energy Delta E_excess.
  5. Multi-origin statistical aggregation and thermal confidence intervals.
  6. End-to-end multi-origin PME and DISH execution.
"""

import os
import shutil
import tempfile
import unittest
import numpy as np

from qdex.namd.ensemble import (
    compute_energy_autocorrelation_time,
    estimate_pilot_cooling_time,
    auto_calibrate_ensemble_origins,
    load_origin_frame_data,
    sample_origin_initial_states,
    aggregate_multi_origin_results,
)
from qdex.namd.surface_hopping import run_namd_dynamics


class TestMultiOriginEnsemble(unittest.TestCase):
    def setUp(self):
        np.random.seed(42)
        self.temp_dir = tempfile.mkdtemp(prefix="qdex_test_ensemble_")

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_compute_energy_autocorrelation_time(self):
        """Test ACF extraction on a known damped oscillation."""
        dt_fs = 2.0
        times = np.arange(500) * dt_fs
        # Damped harmonic oscillation: tau = 100 fs, omega = 150 cm^-1 (period ~ 220 fs)
        tau_true = 100.0
        signal = np.exp(-times / tau_true) * np.cos(2.0 * np.pi * times / 220.0)

        tau_corr = compute_energy_autocorrelation_time(times, signal)
        self.assertTrue(np.isfinite(tau_corr))
        self.assertGreater(tau_corr, 10.0)
        self.assertLess(tau_corr, 200.0)

        # Test constant signal fallback
        flat_signal = np.ones_like(times) * 2.5
        tau_flat = compute_energy_autocorrelation_time(times, flat_signal)
        self.assertEqual(tau_flat, 200.0)

    def _create_synthetic_precomputed_data(self, n_frames=50, n_occ=4, n_virt=4, dt_fs=2.0):
        """Creates minimal synthetic precomputed NAMD files for testing."""
        precompute_dir = os.path.join(self.temp_dir, "precompute")
        os.makedirs(precompute_dir, exist_ok=True)

        i_pairs = []
        a_pairs = []
        for i in range(n_occ):
            for a in range(n_virt):
                i_pairs.append(i)
                a_pairs.append(a)
        i_pairs = np.array(i_pairs, dtype=int)
        a_pairs = np.array(a_pairs, dtype=int)
        n_pairs = len(i_pairs)

        eps_occ0 = np.linspace(-3.0, -1.0, n_occ)
        eps_virt0 = np.linspace(1.0, 3.0, n_virt)
        E0 = eps_virt0[a_pairs] - eps_occ0[i_pairs]
        f0 = np.random.uniform(0.1, 1.0, n_pairs)

        np.savez_compressed(
            os.path.join(precompute_dir, "namd_metadata.npz"),
            n_occ=n_occ,
            n_virt=n_virt,
            n_frames=n_frames,
            dt_fs=dt_fs,
            dt_nuc_fs=dt_fs,
            soc=False,
            qp_gap=float(eps_virt0[0] - eps_occ0[-1]),
            dft_gap=float(eps_virt0[0] - eps_occ0[-1])
        )

        np.savez_compressed(
            os.path.join(precompute_dir, "frame_00000.npz"),
            time_fs=0.0,
            eps_occ=eps_occ0,
            eps_virt=eps_virt0,
            i_pairs=i_pairs,
            a_pairs=a_pairs,
            E_pairs=E0,
            f_pairs=f0,
            scissor=0.0
        )

        for k in range(n_frames - 1):
            t_prev = k * dt_fs
            t_curr = (k + 1) * dt_fs
            # Introduce small thermal fluctuating shift
            fluct = 0.05 * np.sin(2.0 * np.pi * t_curr / 100.0)
            eps_occ_curr = eps_occ0 + fluct
            eps_virt_curr = eps_virt0 - fluct
            E_curr = eps_virt_curr[a_pairs] - eps_occ_curr[i_pairs]

            # Identity overlap with small coupling
            S_occ = np.eye(n_occ, dtype=np.complex128) + 1e-4 * (np.random.randn(n_occ, n_occ) + 1j * np.random.randn(n_occ, n_occ))
            S_virt = np.eye(n_virt, dtype=np.complex128) + 1e-4 * (np.random.randn(n_virt, n_virt) + 1j * np.random.randn(n_virt, n_virt))

            np.savez_compressed(
                os.path.join(precompute_dir, f"step_{k:05d}_to_{k+1:05d}.npz"),
                time_prev_fs=t_prev,
                time_curr_fs=t_curr,
                eps_occ_prev=eps_occ0 if k == 0 else eps_occ_curr - fluct,
                eps_occ_curr=eps_occ_curr,
                eps_virt_prev=eps_virt0 if k == 0 else eps_virt_curr + fluct,
                eps_virt_curr=eps_virt_curr,
                E_prev=E0 if k == 0 else E_curr,
                E_curr=E_curr,
                f_curr=f0,
                i_pairs=i_pairs,
                a_pairs=a_pairs,
                S_occ=S_occ,
                S_virt=S_virt,
                max_norm_loss=0.0
            )

        return precompute_dir, float(eps_virt0[0] - eps_occ0[-1])

    def test_auto_calibrate_short_trajectory_fallback(self):
        """Verify graceful fallback to single origin when trajectory is short."""
        precompute_dir, qp_gap = self._create_synthetic_precomputed_data(n_frames=10, dt_fs=2.0)
        # Total MD duration = 18 fs, cooling window = 50 fs -> T_span < 0
        dyn_cfg = {"window_fs": 50.0}
        calib = auto_calibrate_ensemble_origins(
            precompute_dir=precompute_dir,
            qp_gap_ev=qp_gap,
            pump_energy_ev=qp_gap + 1.0,
            dyn_cfg=dyn_cfg
        )
        self.assertTrue(calib["is_fallback"])
        self.assertEqual(calib["n_origins"], 1)
        self.assertEqual(list(calib["origin_frames"]), [0])

    def test_auto_calibrate_long_trajectory(self):
        """Verify generation of multiple independent origins on a sufficiently long trajectory."""
        precompute_dir, qp_gap = self._create_synthetic_precomputed_data(n_frames=150, dt_fs=2.0)
        # Total MD = 298 fs, set window = 80 fs, dt0 = 50 fs -> span = 218 fs -> origins ~ 5
        dyn_cfg = {"window_fs": 80.0, "dt0_fs": 50.0}
        calib = auto_calibrate_ensemble_origins(
            precompute_dir=precompute_dir,
            qp_gap_ev=qp_gap,
            pump_energy_ev=qp_gap + 1.0,
            dyn_cfg=dyn_cfg
        )
        self.assertFalse(calib["is_fallback"])
        self.assertGreater(calib["n_origins"], 2)
        # Verify origins are non-decreasing and within bounds
        frames = calib["origin_frames"]
        self.assertEqual(frames[0], 0)
        self.assertTrue(np.all(np.diff(frames) > 0))
        self.assertLessEqual(frames[-1] + calib["n_win_steps"], 150)

    def test_energy_referenced_sampling(self):
        """Verify that sample_origin_initial_states references instantaneous band gap."""
        precompute_dir, qp_gap_0 = self._create_synthetic_precomputed_data(n_frames=30, dt_fs=2.0)
        dyn_mask = np.ones(16, dtype=bool)
        pump_excess = 1.0

        # Sample at frame 0
        sampled_0, _, orig_0 = sample_origin_initial_states(
            precompute_dir=precompute_dir,
            k0=0,
            pump_excess_ev=pump_excess,
            dyn_mask=dyn_mask,
            dyn_cfg={},
            n_trajectories=10
        )
        # Sample at frame 10 (which has dynamic gap fluctuation)
        sampled_10, _, orig_10 = sample_origin_initial_states(
            precompute_dir=precompute_dir,
            k0=10,
            pump_excess_ev=pump_excess,
            dyn_mask=dyn_mask,
            dyn_cfg={},
            n_trajectories=10
        )

        self.assertAlmostEqual(orig_0["pump_energy_ev"] - orig_0["qp_gap"], pump_excess, places=6)
        self.assertAlmostEqual(orig_10["pump_energy_ev"] - orig_10["qp_gap"], pump_excess, places=6)

    def test_aggregate_multi_origin_results(self):
        """Verify grand ensemble mean and thermal standard deviation aggregation."""
        times = np.array([0.0, 2.0, 4.0])
        n_frames = len(times)
        n_states = 5

        # 3 mock origins
        res1 = {
            "mean_energy": np.array([3.0, 2.5, 2.0]),
            "mean_excess_e": np.array([0.5, 0.3, 0.1]),
            "mean_excess_h": np.array([0.5, 0.2, 0.1]),
            "populations": np.zeros((n_frames, n_states)),
            "trajectory_energies": np.ones((n_frames, 10)) * 2.0
        }
        res2 = {
            "mean_energy": np.array([3.2, 2.7, 2.2]),
            "mean_excess_e": np.array([0.6, 0.4, 0.2]),
            "mean_excess_h": np.array([0.6, 0.3, 0.2]),
            "populations": np.zeros((n_frames, n_states)),
            "trajectory_energies": np.ones((n_frames, 10)) * 2.2
        }

        agg = aggregate_multi_origin_results([res1, res2], times_fs=times)

        # Expected means
        np.testing.assert_allclose(agg["mean_energy"], [3.1, 2.6, 2.1])
        # Expected sample standard deviation with ddof=1: std([3.0, 3.2]) = sqrt(0.02) ≈ 0.14142
        self.assertAlmostEqual(agg["std_energy"][0], np.sqrt(0.02), places=5)
        # Expected concatenated trajectory shape: (3, 20)
        self.assertEqual(agg["trajectory_energies"].shape, (3, 20))

    def test_end_to_end_pme_multi_origin(self):
        """End-to-end test of run_namd_dynamics with PME and initial_conditions: multiple."""
        precompute_dir, qp_gap = self._create_synthetic_precomputed_data(n_frames=60, dt_fs=2.0)
        out_dir = os.path.join(self.temp_dir, "namd_out")
        os.makedirs(out_dir, exist_ok=True)

        config = {
            "namd": {
                "storage": {"precompute_dir": precompute_dir},
                "dynamics": {
                    "method": "pme",
                    "initial_conditions": "multiple",
                    "window_fs": 40.0,
                    "dt0_fs": 20.0,
                    "pump_energy_ev": qp_gap + 0.8,
                },
                "output": {
                    "cooling_curve_csv": os.path.join(out_dir, "carrier_cooling.csv"),
                    "populations_npz": os.path.join(out_dir, "populations.npz"),
                    "plot_cooling": False,
                }
            }
        }

        run_namd_dynamics(config)

        csv_path = os.path.join(out_dir, "carrier_cooling.csv")
        npz_path = os.path.join(out_dir, "populations.npz")
        self.assertTrue(os.path.exists(csv_path))
        self.assertTrue(os.path.exists(npz_path))

        data = np.load(npz_path)
        self.assertIn("std_energies_ev", data)
        self.assertIn("n_origins", data)
        self.assertGreater(int(data["n_origins"]), 1)

    def test_end_to_end_dish_multi_origin(self):
        """End-to-end test of run_namd_dynamics with DISH and initial_conditions: multiple."""
        precompute_dir, qp_gap = self._create_synthetic_precomputed_data(n_frames=60, dt_fs=2.0)
        out_dir = os.path.join(self.temp_dir, "namd_dish_out")
        os.makedirs(out_dir, exist_ok=True)

        config = {
            "namd": {
                "storage": {"precompute_dir": precompute_dir},
                "dynamics": {
                    "method": "dish",
                    "initial_conditions": "multiple",
                    "window_fs": 40.0,
                    "dt0_fs": 20.0,
                    "n_trajectories": 12,
                    "pump_energy_ev": qp_gap + 0.8,
                },
                "output": {
                    "cooling_curve_csv": os.path.join(out_dir, "carrier_cooling.csv"),
                    "populations_npz": os.path.join(out_dir, "populations.npz"),
                    "plot_cooling": False,
                }
            }
        }

        run_namd_dynamics(config)

        csv_path = os.path.join(out_dir, "carrier_cooling.csv")
        self.assertTrue(os.path.exists(csv_path))


if __name__ == "__main__":
    unittest.main()
