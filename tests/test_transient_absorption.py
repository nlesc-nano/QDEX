"""
Unit tests for QDEX Ultrafast Pump-Probe Transient Absorption (TA) module.
"""

import os
import unittest
import numpy as np
import tempfile
import pandas as pd

from qdex.namd.transient_absorption import (
    compute_transient_absorption,
    fit_bleach_rise_kinetics,
    plot_transient_absorption,
    export_transient_absorption_data,
)


class TestTransientAbsorption(unittest.TestCase):
    def setUp(self):
        # 3-level exciton cascade model:
        # State 0: 1S Band Edge (E = 2.0 eV), i=0, a=0, f = 1.0
        # State 1: Intermediate (E = 3.0 eV), i=1, a=1, f = 0.8
        # State 2: Hot Pumped State (E = 4.0 eV), i=2, a=2, f = 0.5
        self.n_times = 51
        self.times_fs = np.linspace(0.0, 1000.0, self.n_times)  # 0 to 1 ps
        self.E_pairs = np.array([2.0, 3.0, 4.0])
        self.f_pairs = np.array([1.0, 0.8, 0.5])
        self.i_pairs = np.array([0, 1, 2])
        self.a_pairs = np.array([0, 1, 2])

        # Cascade dynamics:
        # P_2(t) decays: exp(-t / 150)
        # P_0(t) rises: 1 - exp(-t / 250)
        tau_rise = 250.0
        p0 = 1.0 - np.exp(- self.times_fs / tau_rise)
        p2 = np.exp(- self.times_fs / 150.0)
        p1 = np.maximum(0.0, 1.0 - p0 - p2)

        self.populations = np.column_stack([p0, p1, p2])
        self.expected_tau_rise_fs = tau_rise

    def test_fit_bleach_rise_kinetics(self):
        """Verify exponential rise fit accurately extracts known cooling timescale."""
        tau_true = 300.0
        times = np.linspace(0.0, 1500.0, 60)
        # S(t) = - Delta A (positive rising signal)
        delta_A_1s = - (1.0 - np.exp(- times / tau_true))

        fit = fit_bleach_rise_kinetics(times, delta_A_1s)
        self.assertTrue(fit["success"])
        self.assertAlmostEqual(fit["tau_rise_fs"], tau_true, delta=5.0)
        self.assertAlmostEqual(fit["k_cool_ps"], 1.0 / (tau_true * 1e-3), places=2)

    def test_compute_transient_absorption_structure(self):
        """Verify 2D TA map, dimensions, and sign conventions."""
        res = compute_transient_absorption(
            times_fs=self.times_fs,
            populations=self.populations,
            E_pairs=self.E_pairs,
            f_pairs=self.f_pairs,
            i_pairs=self.i_pairs,
            a_pairs=self.a_pairs,
            sigma_ev=0.05,
            e_range=(1.5, 4.5),
            n_e_points=100,
            include_se=True,
        )

        self.assertEqual(res["delta_A"].shape, (self.n_times, 100))
        self.assertEqual(len(res["delta_A_1s"]), self.n_times)
        self.assertEqual(res["e_1s_ev"], 2.0)

        # Differential absorption is negative (bleach / stimulated emission)
        self.assertLessEqual(np.min(res["delta_A"]), 0.0)
        # Final 1S bleach should be deeply negative
        self.assertLess(res["delta_A_1s"][-1], res["delta_A_1s"][0])

        # Fitted rise time should be close to 250 fs
        fit = res["fit_results"]
        self.assertTrue(fit["success"])
        self.assertAlmostEqual(fit["tau_rise_fs"], self.expected_tau_rise_fs, delta=20.0)

    def test_plot_and_export_transient_absorption(self):
        """Verify plot generation and data export."""
        res = compute_transient_absorption(
            times_fs=self.times_fs,
            populations=self.populations,
            E_pairs=self.E_pairs,
            f_pairs=self.f_pairs,
            i_pairs=self.i_pairs,
            a_pairs=self.a_pairs,
            sigma_ev=0.05,
            e_range=(1.5, 4.5),
            n_e_points=50,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            plot_file = os.path.join(tmpdir, "test_ta_map.png")
            csv_file = os.path.join(tmpdir, "test_kinetics.csv")
            npz_file = os.path.join(tmpdir, "test_map.npz")

            plot_transient_absorption(res, plot_file=plot_file, material_name="CSPBBR3")
            self.assertTrue(os.path.exists(plot_file))
            self.assertGreater(os.path.getsize(plot_file), 1000)

            export_transient_absorption_data(res, kinetics_csv=csv_file, map_npz=npz_file)
            self.assertTrue(os.path.exists(csv_file))
            self.assertTrue(os.path.exists(npz_file))

            df = pd.read_csv(csv_file)
            self.assertEqual(len(df), self.n_times)
            self.assertIn("delta_A_1s", df.columns)
            self.assertIn("bleach_amplitude", df.columns)


if __name__ == "__main__":
    unittest.main()
