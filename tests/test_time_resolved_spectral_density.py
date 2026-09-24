"""
Unit tests for Time-Resolved Vibrational Action Spectrum / Dynamical Phonon Spectrogram J(omega, t).
"""

import os
import unittest
import tempfile
import numpy as np

from qdex.namd.analysis import (
    compute_pair_spectral_density,
    compute_time_resolved_spectral_density,
    plot_time_resolved_spectral_density,
    export_time_resolved_spectral_density,
    compute_2d_vibronic_action_map,
    plot_2d_vibronic_action_map,
    export_2d_vibronic_action_map,
)


class TestTimeResolvedSpectralDensity(unittest.TestCase):
    def setUp(self):
        np.random.seed(42)
        self.n_frames = 200
        self.dt_fs = 2.0
        self.times_fs = np.arange(self.n_frames) * self.dt_fs
        self.n_occ = 6
        self.n_virt = 6

        # Harmonic oscillation at 150 cm^-1 for orbital 1 vs 0
        c_cm_per_fs = 2.99792458e-5
        target_wn = 150.0
        period_fs = 1.0 / (target_wn * c_cm_per_fs)

        self.eps_occ_traj = np.random.randn(self.n_frames, self.n_occ) * 0.01
        self.eps_virt_traj = np.random.randn(self.n_frames, self.n_virt) * 0.01

        # Embed clean mode in virtual orbitals 0 and 1
        t_grid = self.times_fs
        self.eps_virt_traj[:, 1] = self.eps_virt_traj[:, 0] + 0.1 * np.sin(2.0 * np.pi * t_grid / period_fs)

    def test_compute_pair_spectral_density_peak_recovery(self):
        """Verify compute_pair_spectral_density recovers embedded 150 cm^-1 phonon frequency."""
        e1 = self.eps_virt_traj[:, 1]
        e0 = self.eps_virt_traj[:, 0]
        wn, psd, var_d = compute_pair_spectral_density(e1, e0, self.dt_fs, w_max_cm=300.0)

        self.assertTrue(len(wn) > 10)
        self.assertEqual(len(wn), len(psd))
        self.assertTrue(np.all(psd >= 0.0), "PSD contains negative values")

        peak_freq = wn[np.argmax(psd)]
        self.assertAlmostEqual(peak_freq, 150.0, delta=5.0,
                               msg=f"Expected peak near 150 cm^-1, got {peak_freq:.1f} cm^-1")

    def test_compute_pair_spectral_density_constant_zero(self):
        """Verify identical orbital energies produce exactly zero spectral density."""
        e_const = np.ones(self.n_frames)
        wn, psd, var_d = compute_pair_spectral_density(e_const, e_const, self.dt_fs, w_max_cm=300.0)
        self.assertAlmostEqual(var_d, 0.0, places=12)
        np.testing.assert_allclose(psd, 0.0, atol=1e-12)

    def test_compute_time_resolved_spectral_density_hops(self):
        """Verify spectrogram construction from discrete stochastic surface hops (FSSH/DISH)."""
        hops = [
            {"time_fs": 40.0, "channel": "electron", "from": 1, "to": 0, "traj": 0},
            {"time_fs": 42.0, "channel": "electron", "from": 1, "to": 0, "traj": 1},
            {"time_fs": 45.0, "channel": "electron", "from": 1, "to": 0, "traj": 2},
            {"time_fs": 120.0, "channel": "hole", "from": 3, "to": 2, "traj": 0},
        ]

        res = compute_time_resolved_spectral_density(
            times_fs=self.times_fs,
            eps_occ_traj=self.eps_occ_traj,
            eps_virt_traj=self.eps_virt_traj,
            hop_records=hops,
            n_trajectories=10,
            sigma_t_fs=10.0,
            w_max_cm=300.0,
            method_name="DISH",
        )

        J_tot = res["J_total"]
        wn = res["wavenumbers_cm"]
        self.assertEqual(J_tot.shape, (len(wn), len(self.times_fs)))
        self.assertTrue(np.all(J_tot >= 0.0))
        self.assertAlmostEqual(np.max(J_tot), 1.0, places=5)

        # Check that action density is localized near t ~ 40 fs
        t_axis = self.times_fs
        mask_early = (t_axis >= 30.0) & (t_axis <= 55.0)
        mask_late = (t_axis >= 80.0) & (t_axis <= 100.0)
        self.assertGreater(np.mean(J_tot[:, mask_early]), np.mean(J_tot[:, mask_late]) * 5.0)

    def test_compute_time_resolved_spectral_density_pme_flux(self):
        """Verify spectrogram construction from Pauli Master Equation continuous fluxes."""
        pme_flux = []
        for k, t in enumerate(self.times_fs):
            flux_v = np.zeros((self.n_virt, self.n_virt))
            flux_o = np.zeros((self.n_occ, self.n_occ))
            if 30 <= k <= 40:
                flux_v[1, 0] = 0.05
            pme_flux.append({
                "time_fs": t,
                "flux_virt": flux_v,
                "flux_occ": flux_o,
            })

        res = compute_time_resolved_spectral_density(
            times_fs=self.times_fs,
            eps_occ_traj=self.eps_occ_traj,
            eps_virt_traj=self.eps_virt_traj,
            pme_flux_records=pme_flux,
            sigma_t_fs=8.0,
            w_max_cm=300.0,
            method_name="PME",
        )

        J_tot = res["J_total"]
        self.assertTrue(np.all(J_tot >= 0.0))
        self.assertAlmostEqual(np.max(J_tot), 1.0, places=5)

        # Peak of 1D integrated spectrum should be near 150 cm^-1
        wn = res["wavenumbers_cm"]
        J_int = res["J_int_total"]
        peak_wn = wn[np.argmax(J_int)]
        self.assertAlmostEqual(peak_wn, 150.0, delta=6.0)

    def test_plot_and_export_pipeline(self):
        """Verify plotting (PNG + HTML) and data export (NPZ + CSV) execute without errors."""
        hops = [{"time_fs": 50.0, "channel": "electron", "from": 1, "to": 0, "traj": 0}]
        res = compute_time_resolved_spectral_density(
            times_fs=self.times_fs,
            eps_occ_traj=self.eps_occ_traj,
            eps_virt_traj=self.eps_virt_traj,
            hop_records=hops,
            n_trajectories=1,
            sigma_t_fs=15.0,
            w_max_cm=300.0,
            method_name="CPA-FSSH-GDC",
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            png_file = os.path.join(tmpdir, "tr_sd.png")
            html_file = os.path.join(tmpdir, "tr_sd.html")
            npz_file = os.path.join(tmpdir, "tr_sd.npz")
            csv_file = os.path.join(tmpdir, "tr_sd.csv")

            plot_time_resolved_spectral_density(
                res, plot_file=png_file, html_file=html_file, material_name="CsPbBr3", method_name="GDC"
            )
            export_time_resolved_spectral_density(
                res, output_npz=npz_file, output_csv=csv_file
            )

            self.assertTrue(os.path.exists(png_file))
            self.assertTrue(os.path.getsize(png_file) > 1000)
            self.assertTrue(os.path.exists(html_file))
            self.assertTrue(os.path.getsize(html_file) > 1000)
            self.assertTrue(os.path.exists(npz_file))
            self.assertTrue(os.path.exists(csv_file))

            # Verify exported npz content
            d = np.load(npz_file)
            self.assertIn("J_total", d.files)
            self.assertIn("wavenumbers_cm", d.files)
            self.assertIn("times_fs", d.files)

    def test_2d_vibronic_action_map_pipeline(self):
        """Verify 2D vibronic action map computation, plotting, and export."""
        # Create synthetic time-resolved data with periodic bursts at period ~200 fs (167 cm^-1)
        # and accepting mode at 25 cm^-1
        hops = []
        for t in np.arange(20.0, 380.0, 40.0):
            hops.append({"time_fs": float(t), "channel": "electron", "from": 1, "to": 0, "traj": 0})

        tr_data = compute_time_resolved_spectral_density(
            times_fs=self.times_fs,
            eps_occ_traj=self.eps_occ_traj,
            eps_virt_traj=self.eps_virt_traj,
            hop_records=hops,
            n_trajectories=1,
            sigma_t_fs=10.0,
            w_max_cm=300.0,
            method_name="DISH",
        )

        vib2d = compute_2d_vibronic_action_map(
            tr_data, n_fft=1024, detrend_mode="linear", wmax_acc=200.0, wmax_prom=300.0
        )

        self.assertIn("S_2d_total", vib2d)
        self.assertIn("wavenumbers_acc_cm", vib2d)
        self.assertIn("wavenumbers_prom_cm", vib2d)
        self.assertIn("proj_acc", vib2d)
        self.assertIn("proj_prom", vib2d)
        self.assertEqual(vib2d["S_2d_total"].shape, (len(vib2d["wavenumbers_acc_cm"]), len(vib2d["wavenumbers_prom_cm"])))

        with tempfile.TemporaryDirectory() as tmpdir:
            png_file = os.path.join(tmpdir, "map2d.png")
            html_file = os.path.join(tmpdir, "map2d.html")
            npz_file = os.path.join(tmpdir, "map2d.npz")
            csv_file = os.path.join(tmpdir, "map2d.csv")

            plot_2d_vibronic_action_map(
                vib2d, output_png=png_file, output_html=html_file, material_name="CsPbBr3", method_name="DISH"
            )
            export_2d_vibronic_action_map(
                vib2d, output_npz=npz_file, output_csv=csv_file
            )

            self.assertTrue(os.path.exists(png_file))
            self.assertTrue(os.path.getsize(png_file) > 1000)
            self.assertTrue(os.path.exists(html_file))
            self.assertTrue(os.path.getsize(html_file) > 1000)
            self.assertTrue(os.path.exists(npz_file))
            self.assertTrue(os.path.exists(csv_file))


if __name__ == "__main__":
    unittest.main()
