"""Boundary tests for the sustained 3 dB capture-window semantics."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from main_THZ_220GHz_deterministic_mechanical_tracking import (
    _dwell_window_sample_count,
    _first_true_window_start,
)


class CaptureWindowTests(unittest.TestCase):
    def test_no_lock(self) -> None:
        self.assertIsNone(_first_true_window_start(np.array([False, False]), 2))

    def test_single_point_does_not_satisfy_two_sample_dwell(self) -> None:
        self.assertIsNone(_first_true_window_start(np.array([False, True, False]), 2))

    def test_exact_window_returns_its_start(self) -> None:
        self.assertEqual(
            _first_true_window_start(np.array([False, True, True, False]), 2), 1
        )

    def test_tail_shorter_than_window_is_rejected(self) -> None:
        self.assertIsNone(_first_true_window_start(np.array([False, True, True]), 3))

    def test_drop_and_relock_returns_second_window(self) -> None:
        self.assertEqual(
            _first_true_window_start(
                np.array([True, True, False, True, True, True]), 3
            ),
            3,
        )

    def test_invalid_window_length_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            _first_true_window_start(np.array([True]), 0)

    def test_dwell_samples_cover_elapsed_span(self) -> None:
        self.assertEqual(_dwell_window_sample_count(1.5, 0.5), 4)
        self.assertEqual(_dwell_window_sample_count(1.51, 0.5), 5)
        self.assertEqual(_dwell_window_sample_count(0.0, 0.5), 1)

    def test_invalid_dwell_parameters_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            _dwell_window_sample_count(-0.1, 0.5)
        with self.assertRaises(ValueError):
            _dwell_window_sample_count(1.0, 0.0)


if __name__ == '__main__':
    unittest.main()
