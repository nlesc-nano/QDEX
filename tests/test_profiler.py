import unittest
import io
import sys
import numpy as np
from qdex.profiler import get_memory_info_mb, ResourceTracker

class TestProfiler(unittest.TestCase):
    def test_get_memory_info(self):
        curr, peak = get_memory_info_mb()
        self.assertGreater(curr, 0.0)
        self.assertGreater(peak, 0.0)
        self.assertGreaterEqual(peak, curr * 0.5)

    def test_resource_tracker_stages(self):
        tracker = ResourceTracker()

        with tracker.stage("Stage 1"):
            # Allocate some memory
            arr = np.ones((1000, 1000), dtype=np.float64)
            _ = arr.sum()

        tracker.start_stage("Stage 2")
        arr2 = np.ones((500, 500), dtype=np.float64)
        _ = arr2.sum()
        tracker.end_stage()

        self.assertEqual(len(tracker.stages), 2)
        self.assertEqual(tracker.stages[0]["name"], "Stage 1")
        self.assertEqual(tracker.stages[1]["name"], "Stage 2")
        self.assertGreater(tracker.stages[0]["elapsed"], 0.0)
        self.assertGreater(tracker.stages[1]["elapsed"], 0.0)

    def test_print_summary(self):
        tracker = ResourceTracker()
        with tracker.stage("Test Stage"):
            pass

        captured_output = io.StringIO()
        old_stdout = sys.stdout
        try:
            sys.stdout = captured_output
            tracker.print_summary(device="numpy", nthreads=4)
        finally:
            sys.stdout = old_stdout

        output = captured_output.getvalue()
        self.assertIn("COMPUTATIONAL RESOURCE USAGE SUMMARY", output)
        self.assertIn("Test Stage", output)
        self.assertIn("Total Execution Time", output)

if __name__ == "__main__":
    unittest.main()
