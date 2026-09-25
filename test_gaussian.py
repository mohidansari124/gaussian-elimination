"""CLI tests for the compiled solver; no third-party dependencies."""
import itertools
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent
EXE = ROOT / ("gaussian.exe" if (ROOT / "gaussian.exe").exists() else "gaussian")


class SolverTests(unittest.TestCase):
    def run_system(self, data, flags=(), success=True):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "test system.lin"
            source.write_text(data)
            result = subprocess.run([str(EXE), *flags, str(source)],
                                    capture_output=True, text=True)
            solution = source.with_suffix(".sol")
            if success:
                self.assertEqual(result.returncode, 0, result.stderr)
                lines = solution.read_text().splitlines()
                self.assertEqual(len(lines), 1)
                return list(map(float, lines[0].split()))
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(solution.exists())

    def test_known_solution_all_modes(self):
        for spp, double in itertools.product([False, True], repeat=2):
            flags = (["--spp"] if spp else []) + (["--double"] if double else [])
            with self.subTest(flags=flags):
                x = self.run_system((ROOT / "sys1.lin").read_text(), flags)
                for actual, expected in zip(x, [2, 3, -1]):
                    self.assertAlmostEqual(actual, expected, places=5)
                self.assertEqual(len(x), 3)

    def test_row_exchange(self):
        data = "2\n0 2\n1 3\n4 7\n"
        self.run_system(data, success=False)
        self.assertEqual(self.run_system(data, ["--spp"]), [1, 2])

    def test_scaled_pivot_and_precision(self):
        # Scaled ratios select row 2 despite its smaller absolute coefficient.
        data = "2\n100 1e10\n1 1\n10000000100 2\n"
        naive = self.run_system(data)
        scaled = self.run_system(data, ["--spp"])
        double = self.run_system(data, ["--double"])
        self.assertEqual(naive, [0, 1])
        self.assertEqual(scaled, [1, 1])
        for actual in double:
            self.assertAlmostEqual(actual, 1, places=6)

    def test_single_equation(self):
        self.assertEqual(self.run_system("1\n4\n10\n"), [2.5])

    def test_singular_and_invalid_input(self):
        for data in ["2\n1 2\n2 4\n3 6", "0", "1.5\n2\n3",
                     "2\n1 2", "1\n2\n3 extra", "1\nNaN\n1"]:
            for flags in [[], ["--spp"]]:
                with self.subTest(data=data, flags=flags):
                    self.run_system(data, flags, success=False)

    def test_bad_arguments(self):
        for args in [[], ["--unknown"], ["wrong.txt"], ["missing.lin"]]:
            result = subprocess.run([str(EXE), *args], capture_output=True)
            self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()
