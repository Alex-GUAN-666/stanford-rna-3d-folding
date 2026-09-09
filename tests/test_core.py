"""Behavioral tests for input contracts, rigid assembly, and safe ranking."""

import tempfile
from pathlib import Path
import unittest

import numpy as np

from rna_folding.geometry import kabsch_align, stitch_segments, validate_coordinates
from rna_folding.selection import Candidate, select_candidates
from rna_folding.sequence import RNASequence, load_sequences, route_model, segment_ranges


def curve(length: int) -> np.ndarray:
    t = np.arange(length, dtype=float)
    return np.column_stack((np.cos(t * 0.7), np.sin(t * 0.7), t * 0.25))


class SequenceTests(unittest.TestCase):
    def test_routing_boundaries(self):
        for length, expected in [(1, "trrosettarna"), (100, "trrosettarna"),
                                 (101, "drfold2"), (480, "drfold2"), (481, "protenix")]:
            with self.subTest(length=length):
                self.assertEqual(route_model(length), expected)
        for bad in [0, -1, 100.5, True]:
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                route_model(bad)

    def test_sequence_requires_canonical_uppercase(self):
        self.assertEqual(RNASequence("target-1.a", "ACGU").sequence, "ACGU")
        for bad in ["", "ACGT", "acgu", "ACGN", "AC GU", "ACGU\n"]:
            with self.subTest(sequence=bad), self.assertRaises(ValueError):
                RNASequence("target", bad)
        for bad in ["", "../target", "a/b", "a b", ".."]:
            with self.subTest(target=bad), self.assertRaises(ValueError):
                RNASequence(bad, "ACGU")

    def test_csv_extras_and_duplicate_rejection(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sequences.csv"
            path.write_text("target_id,sequence,extra\na,ACGU,okay\nb,UGCA,more\n")
            self.assertEqual([record.target_id for record in load_sequences(path)], ["a", "b"])
            path.write_text("target_id,sequence\na,ACGU\na,UGCA\n")
            with self.assertRaisesRegex(ValueError, "duplicate target_id"):
                load_sequences(path)
            for text in ["id,sequence\na,ACGU\n", "target_id,sequence\n",
                         "target_id,sequence\na\n"]:
                path.write_text(text)
                with self.assertRaises(ValueError):
                    load_sequences(path)

    def test_window_coverage_and_full_tail(self):
        self.assertEqual(segment_ranges(100), [(0, 100)])
        self.assertEqual(segment_ranges(551), [(0, 300), (250, 550), (251, 551)])
        for length in [301, 501, 550, 551, 1001]:
            ranges = segment_ranges(length)
            covered = np.zeros(length, dtype=bool)
            self.assertEqual(len(ranges), len(set(ranges)))
            for start, end in ranges:
                self.assertEqual(end - start, 300)
                covered[start:end] = True
            self.assertTrue(covered.all())
        for args in [(0, 300, 50), (10, 0, 0), (10, 3, 3), (10, 3, -1)]:
            with self.subTest(args=args), self.assertRaises(ValueError):
                segment_ranges(*args)


class GeometryTests(unittest.TestCase):
    def test_coordinate_validation(self):
        for bad in [[], [1, 2, 3], [[0, 1]], [[0, 0, np.nan]],
                    [[np.inf, 0, 0]], [[1e6, 0, 0]], [[-1e18, 0, 0]],
                    np.array([[1 + 8j, 2, 3]])]:
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                validate_coordinates(bad)
        with self.assertRaisesRegex(ValueError, "length"):
            validate_coordinates(curve(4), length=5)

    def test_kabsch_recovers_rotation_and_translation(self):
        reference = curve(12)
        rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        moving = reference @ rotation + [13., -8., 2.]
        np.testing.assert_allclose(kabsch_align(moving, reference), reference, atol=1e-12)

    def test_kabsch_never_fits_reflection_as_rotation(self):
        reference = np.array([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.], [0., 0., 3.]])
        reflected = reference * [-1., 1., 1.]
        aligned = kabsch_align(reflected, reference)
        self.assertGreater(float(np.linalg.norm(aligned - reference)), 0.1)

    def test_rotation_invariant_stitching(self):
        reference = curve(20)
        rotation = np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]])
        segments = [(0, reference[:12]), (8, reference[8:] @ rotation + [30., -10., 4.])]
        np.testing.assert_allclose(stitch_segments(segments[::-1], 20), reference, atol=1e-11)

    def test_stitching_rejects_missing_and_degenerate_overlap(self):
        reference = curve(20)
        with self.assertRaisesRegex(ValueError, "gap"):
            stitch_segments([(0, reference[:8]), (9, reference[9:])], 20)
        with self.assertRaisesRegex(ValueError, "at least 3"):
            stitch_segments([(0, reference[:10]), (8, reference[8:])], 20)
        line = np.column_stack((np.arange(12), np.zeros(12), np.zeros(12)))
        with self.assertRaisesRegex(ValueError, "non-collinear"):
            stitch_segments([(0, line[:8]), (4, line[4:])], 12)
        with self.assertRaisesRegex(ValueError, "complete target"):
            stitch_segments([(0, reference[:12])], 20)


class SelectionTests(unittest.TestCase):
    def setUp(self):
        self.sequence = RNASequence("target", "ACGUAC")
        self.coordinates = curve(6)

    def candidate(self, candidate_id, model="trrosettarna", energy=None, group=None):
        return Candidate(candidate_id, "target", model, self.coordinates, energy, group)

    def test_round_robin_starts_with_routed_model(self):
        candidates = [self.candidate("p1", "protenix"), self.candidate("t1"),
                      self.candidate("t2"), self.candidate("d1", "drfold2"),
                      self.candidate("p2", "protenix"), self.candidate("t3")]
        selected = select_candidates(candidates, self.sequence)
        self.assertEqual([entry.candidate_id for entry in selected], ["t1", "d1", "p1", "t2", "p2"])

    def test_energy_group_isolation_and_missing_energy_slots(self):
        candidates = [self.candidate("a_high", energy=8., group="a"),
                      self.candidate("b_low", energy=-10000., group="b"),
                      self.candidate("unscored", group="a"),
                      self.candidate("a_low", energy=1., group="a"),
                      self.candidate("ungrouped", energy=-1e10),
                      self.candidate("b_high", energy=200., group="b")]
        result = select_candidates(candidates, self.sequence, count=6)
        self.assertEqual([entry.candidate_id for entry in result],
                         ["a_low", "b_low", "unscored", "a_high", "ungrouped", "b_high"])

    def test_same_group_name_does_not_compare_models(self):
        candidates = [self.candidate("t_high", energy=1000., group="same"),
                      self.candidate("p_low", "protenix", -1000., "same"),
                      self.candidate("t_low", energy=500., group="same")]
        result = select_candidates(candidates, self.sequence, count=3)
        self.assertEqual([entry.candidate_id for entry in result], ["t_low", "p_low", "t_high"])

    def test_strict_candidate_validation(self):
        with self.assertRaisesRegex(ValueError, "at least 5"):
            select_candidates([self.candidate("one")], self.sequence)
        with self.assertRaisesRegex(ValueError, "duplicate"):
            select_candidates([self.candidate("same"), self.candidate("same")], self.sequence, 2)
        bad_target = Candidate("other", "wrong", "protenix", self.coordinates)
        with self.assertRaisesRegex(ValueError, "different target"):
            select_candidates([bad_target], self.sequence, 1)
        bad_shape = Candidate("bad", "target", "protenix", curve(5))
        with self.assertRaisesRegex(ValueError, "length"):
            select_candidates([bad_shape], self.sequence, 1)
        with self.assertRaisesRegex(ValueError, "finite"):
            select_candidates([self.candidate("bad", energy=np.nan)], self.sequence, 1)


if __name__ == "__main__":
    unittest.main()
