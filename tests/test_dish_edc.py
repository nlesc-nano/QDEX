"""
Unit tests for Decoherence-Induced Surface Hopping (DISH) and Energy-based
Decoherence Corrections (EDC) in QDEX NAMD.
"""

import os
import unittest
import numpy as np
import tempfile

from qdex.namd.integrator import (
    apply_edc_decoherence_batch,
    apply_gdc_decoherence_batch,
    step_dish_batch,
    propagate_channel_batch_strang,
    KB_EV
)
from qdex.namd.precompute import compute_trajectory_decoherence_times


class TestDishAndEdc(unittest.TestCase):
    def setUp(self):
        np.random.seed(42)
        self.n_states = 10
        self.n_traj = 50
        self.dt_fs = 2.0
        self.temp_k = 300.0
        self.beta = 1.0 / (KB_EV * self.temp_k)

        # Random normalized initial wavepackets
        self.C = (np.random.randn(self.n_states, self.n_traj) +
                  1j * np.random.randn(self.n_states, self.n_traj))
        self.C /= np.linalg.norm(self.C, axis=0, keepdims=True)

        self.active_surfaces = np.random.randint(0, self.n_states, size=self.n_traj)
        self.E_batch = np.sort(np.random.uniform(1.0, 3.0, size=(self.n_states, self.n_traj)), axis=0)

        # Dephasing matrix: off-diagonals ~25 fs, diagonals 500 fs
        self.tau_mat = np.full((self.n_states, self.n_states), 25.0)
        np.fill_diagonal(self.tau_mat, 500.0)

    def test_edc_batch_norm_conservation(self):
        """Verify apply_edc_decoherence_batch preserves exact unit norm across all trajectories."""
        C_damped = apply_edc_decoherence_batch(
            self.C,
            self.active_surfaces,
            self.E_batch,
            self.dt_fs,
            tau_mat=self.tau_mat,
            temp_k=self.temp_k
        )
        norms = np.linalg.norm(C_damped, axis=0)
        np.testing.assert_allclose(norms, 1.0, rtol=1e-10, atol=1e-10,
                                   err_msg="EDC batch damping violated norm conservation")

    def test_edc_batch_damping_rate(self):
        """Verify that non-active amplitudes are damped by exp(-dt / tau)."""
        dt = 5.0
        tau_val = 20.0
        tau_m = np.full((self.n_states, self.n_states), tau_val)
        np.fill_diagonal(tau_m, 500.0)

        # Start from equal superposition
        C_in = np.ones((self.n_states, self.n_traj), dtype=np.complex128) / np.sqrt(self.n_states)
        act = np.zeros(self.n_traj, dtype=int)  # state 0 active

        C_out = apply_edc_decoherence_batch(
            C_in,
            act,
            self.E_batch,
            dt,
            tau_mat=tau_m,
            temp_k=self.temp_k
        )
        expected_damp = np.exp(-dt / tau_val)
        actual_ratio = np.abs(C_out[1, 0]) / np.abs(C_in[1, 0])
        self.assertAlmostEqual(actual_ratio, expected_damp, places=6,
                               msg="EDC did not damp inactive amplitude by expected exponential factor")

    def test_gdc_batch_norm_conservation(self):
        """Verify apply_gdc_decoherence_batch preserves exact unit norm across all trajectories."""
        C_damped = apply_gdc_decoherence_batch(
            self.C,
            self.active_surfaces,
            self.E_batch,
            self.dt_fs,
            tau_mat=self.tau_mat,
            temp_k=self.temp_k
        )
        norms = np.linalg.norm(C_damped, axis=0)
        np.testing.assert_allclose(norms, 1.0, rtol=1e-10, atol=1e-10,
                                   err_msg="GDC batch damping violated norm conservation")

    def test_gdc_batch_damping_rate(self):
        """Verify that non-active amplitudes in GDC are damped by exp(-0.5 * (dt / tau)^2)."""
        dt = 5.0
        tau_val = 20.0
        tau_m = np.full((self.n_states, self.n_states), tau_val)
        np.fill_diagonal(tau_m, 500.0)

        # Start from equal superposition
        C_in = np.ones((self.n_states, self.n_traj), dtype=np.complex128) / np.sqrt(self.n_states)
        act = np.zeros(self.n_traj, dtype=int)  # state 0 active

        C_out = apply_gdc_decoherence_batch(
            C_in,
            act,
            self.E_batch,
            dt,
            tau_mat=tau_m,
            temp_k=self.temp_k
        )
        expected_damp = np.exp(-0.5 * (dt / tau_val) ** 2)
        actual_ratio = np.abs(C_out[1, 0]) / np.abs(C_in[1, 0])
        self.assertAlmostEqual(actual_ratio, expected_damp, places=6,
                               msg="GDC did not damp inactive amplitude by expected Gaussian factor")

    def test_dish_batch_norm_conservation(self):
        """Verify step_dish_batch preserves exact unit norm across all trajectories."""
        C_new, active_new, hopped = step_dish_batch(
            self.C,
            self.active_surfaces,
            self.E_batch,
            self.dt_fs,
            tau_mat=self.tau_mat,
            beta=self.beta,
            detailed_balance=True
        )
        norms = np.linalg.norm(C_new, axis=0)
        np.testing.assert_allclose(norms, 1.0, rtol=1e-10, atol=1e-10,
                                   err_msg="DISH step violated norm conservation")

    def test_dish_collapse_on_hop(self):
        """Verify that when a DISH hop occurs, the wavepacket collapses onto the new active state."""
        # Force a hop by setting large dt and high population on candidate state
        n_traj = 20
        C_test = np.zeros((3, n_traj), dtype=np.complex128)
        act_test = np.zeros(n_traj, dtype=int)
        # Put 99% population on state 1, 1% on state 0
        C_test[0, :] = 0.1
        C_test[1, :] = np.sqrt(0.99)
        E_down = np.zeros((3, n_traj))
        E_down[0, :] = 2.0
        E_down[1, :] = 1.0  # downward transition (no Boltzmann suppression)
        E_down[2, :] = 3.0

        tau_fast = np.full((3, 3), 1.0)
        np.fill_diagonal(tau_fast, 500.0)

        C_new, active_new, hopped = step_dish_batch(
            C_test,
            act_test,
            E_down,
            dt_fs=10.0,
            tau_mat=tau_fast,
            beta=self.beta,
            detailed_balance=True
        )

        for tr in range(n_traj):
            if hopped[tr]:
                new_k = active_new[tr]
                self.assertEqual(new_k, 1, "Hop should have chosen dominant candidate state")
                self.assertAlmostEqual(np.abs(C_new[new_k, tr]), 1.0, places=10,
                                       msg="DISH hop must collapse wavepacket to pure state")
                for other in range(3):
                    if other != new_k:
                        self.assertEqual(C_new[other, tr], 0.0, "DISH hop must zero non-active states")

    def test_dish_upward_detailed_balance(self):
        """Verify that DISH suppresses upward hops at near-zero temperature."""
        n_traj = 100
        C_test = np.zeros((2, n_traj), dtype=np.complex128)
        act_test = np.zeros(n_traj, dtype=int)  # all at ground state 0
        C_test[0, :] = 0.5
        C_test[1, :] = np.sqrt(0.75)
        C_test /= np.linalg.norm(C_test, axis=0, keepdims=True)

        E_up = np.zeros((2, n_traj))
        E_up[0, :] = 1.0
        E_up[1, :] = 2.0  # +1 eV upward transition

        # Near zero temperature: beta -> very large
        beta_cold = 1.0 / (KB_EV * 1e-4)
        tau_fast = np.full((2, 2), 1.0)
        np.fill_diagonal(tau_fast, 500.0)

        C_new, active_new, hopped = step_dish_batch(
            C_test,
            act_test,
            E_up,
            dt_fs=10.0,
            tau_mat=tau_fast,
            beta=beta_cold,
            detailed_balance=True
        )
        self.assertEqual(np.sum(hopped), 0, "No upward hops should occur at near-zero Kelvin in DISH")

    def test_compute_trajectory_decoherence_times(self):
        """Verify compute_trajectory_decoherence_times returns valid, capped matrices."""
        with tempfile.TemporaryDirectory() as tmpdir:
            n_frames = 20
            n_occ = 15
            n_virt = 8
            # Save frame_00000.npz
            np.savez_compressed(
                os.path.join(tmpdir, "frame_00000.npz"),
                eps_occ=np.sort(np.random.randn(n_occ) * 0.1 - 2.0),
                eps_virt=np.sort(np.random.randn(n_virt) * 0.1 + 1.0),
            )
            # Generate fake step files with fluctuating energies
            for k in range(n_frames - 1):
                fpath = os.path.join(tmpdir, f"step_{k:05d}_to_{k+1:05d}.npz")
                eps_occ = np.sort(np.random.randn(n_occ) * 0.1 - 2.0)
                eps_virt = np.sort(np.random.randn(n_virt) * 0.1 + 1.0)
                np.savez_compressed(
                    fpath,
                    eps_occ_curr=eps_occ,
                    eps_virt_curr=eps_virt,
                )

            tau_occ, tau_virt = compute_trajectory_decoherence_times(
                tmpdir, max_tau_fs=500.0, update_metadata=False, verbose=False
            )

            self.assertEqual(tau_occ.shape, (n_occ, n_occ))
            self.assertEqual(tau_virt.shape, (n_virt, n_virt))
            # Diagonals must equal cap
            np.testing.assert_allclose(np.diag(tau_occ), 500.0)
            np.testing.assert_allclose(np.diag(tau_virt), 500.0)
            # Off-diagonals must be positive
            self.assertTrue(np.all(tau_occ > 0.0))
            self.assertTrue(np.all(tau_virt > 0.0))


if __name__ == "__main__":
    unittest.main()
