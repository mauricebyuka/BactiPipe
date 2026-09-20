import contextlib
import io
import sys
import unittest
from unittest.mock import patch
from bactipipe import cli


class CliSurfaceTests(unittest.TestCase):
    def test_help_is_short_and_successful(self):
        output = io.StringIO()
        with patch.object(sys, "argv", ["bactipipe", "--help"]), contextlib.redirect_stdout(output), self.assertRaises(SystemExit) as exc:
            cli.main()
        self.assertEqual(exc.exception.code, 0)
        for hidden in ("plan-identify", "build-identify", "prepare-speciesfinder-db", "render-placement"):
            self.assertNotIn(hidden, output.getvalue())
        self.assertIn("update-db speciesfinder", output.getvalue())

    def test_new_and_legacy_routes(self):
        cases = [(["identify", "--place", "--help"], "place_assembly.py"),
                 (["identify", "--render", "--help"], "render_placement.py"),
                 (["update-db", "speciesfinder", "--help"], "prepare_speciesfinder_database.py"),
                 (["update-db", "speciesfinder", "--verify", "--help"], "validate_species_database.py"),
                 (["prepare-speciesfinder-db", "--help"], "prepare_speciesfinder_database.py"),
                 (["place-assembly", "--help"], "place_assembly.py"),
                 (["detect", "--organism-list"], "run_traits.py")]
        for args, script in cases:
            with self.subTest(args=args), patch.object(sys, "argv", ["bactipipe"] + args), patch.object(cli.subprocess, "run") as run:
                cli.main()
                command = run.call_args.args[0]
                self.assertEqual(command[0], sys.executable)
                self.assertTrue(command[1].endswith(script))
                self.assertEqual(command[2:], [args[-1]])

    def test_modes_are_exclusive(self):
        with patch.object(sys, "argv", ["bactipipe", "identify", "--place", "--render"]), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as exc:
            cli.main()
        self.assertEqual(exc.exception.code, 2)


if __name__ == "__main__": unittest.main()
