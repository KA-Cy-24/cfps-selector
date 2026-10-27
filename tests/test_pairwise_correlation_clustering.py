import json
import unittest

import numpy as np

from analysis.cfps_correlation.pairwise_correlation_clustering import (
    CFPSSubmissionCorrelationAdapter,
    AbstractPairwiseCorrelationClusterer,
    adjusted_rand_index,
)


def _record(record_id, sample_id, target, reagents, batch="batch-a"):
    return {
        "id": record_id,
        "created": "2026-10-27 10:00:00.000Z",
        "design_json": json.dumps({
            "sample_id": sample_id,
            "sample_type": "experimental",
            "reagent_list": {"base_buffer": 2000, "lysate": 5000, **reagents},
            "metadata": {"protein_target": target, "batch_id": batch},
        }),
    }


class TestCFPSAdapter(unittest.TestCase):
    def setUp(self):
        self.records = [
            _record("id-a", "sample-a", "sfGFP", {"r1": 100, "r3": 50, "r5": 25}),
            _record("id-b", "sample-b", "sfGFP", {"r1": 200, "r3": 100, "r5": 50}),
            _record("id-c", "sample-c", "sfGFP", {"r2": 100, "r4": 50}),
            _record("id-d", "sample-d", "sfGFP", {"r2": 50, "r4": 100, "r5": 25}),
            _record("id-e", "sample-e", "PETase", {"r1": 10, "r2": 20, "r3": 5, "r4": 10}),
            _record("id-f", "sample-f", "PETase", {"r1": 20, "r2": 40, "r3": 10, "r4": 20}),
            _record("id-g", "sample-g", "PETase", {"r2": 10, "r3": 30, "r5": 50}),
        ]

    def test_abstract_contract(self):
        self.assertTrue(issubclass(CFPSSubmissionCorrelationAdapter, AbstractPairwiseCorrelationClusterer))

    def test_extracts_sparse_json_to_zero_filled_matrix(self):
        adapter = CFPSSubmissionCorrelationAdapter(self.records, group_by="protein_target")
        X = adapter.build_feature_matrix()
        self.assertEqual(X.loc["id-a", "r2"], 0.0)
        self.assertEqual(X.loc["id-a", "r1"], 100.0)
        self.assertIn("protein_target", adapter.build_metadata().columns)

    def test_constant_features_dropped_and_cluster_fit(self):
        adapter = CFPSSubmissionCorrelationAdapter(
            self.records, correlation_method="spearman", linkage_method="average", n_clusters=2,
            group_by="protein_target", min_features=3,
        ).fit()
        self.assertIn("sfGFP", adapter.results_)
        self.assertIn("PETase", adapter.results_)
        sf = adapter.results_["sfGFP"]
        self.assertIn("base_buffer", sf.dropped_features)
        self.assertIn("lysate", sf.dropped_features)
        self.assertGreater(sf.correlation_matrix.loc["id-a", "id-b"], 0.99)
        self.assertTrue(np.allclose(sf.distance_matrix.to_numpy(), sf.distance_matrix.to_numpy().T))
        self.assertEqual(len(adapter.labels_frame()), 7)

    def test_presence_mode(self):
        adapter = CFPSSubmissionCorrelationAdapter(self.records, feature_mode="presence", n_clusters=2).fit()
        self.assertEqual(adapter.results_["sfGFP"].feature_mode, "presence")
        self.assertTrue(set(np.unique(adapter.results_["sfGFP"].feature_matrix.to_numpy())).issubset({0.0, 1.0}))

    def test_bad_volume_fails_loudly(self):
        records = [_record("bad", "bad", "sfGFP", {"r1": "not-a-volume", "r2": 10, "r3": 25})]
        with self.assertRaisesRegex(ValueError, "invalid/nonnegative-volume"):
            CFPSSubmissionCorrelationAdapter(records)

    def test_snapshot_mapping_input(self):
        snapshot = {"collections": {"cfps_designs": self.records}}
        adapter = CFPSSubmissionCorrelationAdapter(snapshot)
        self.assertEqual(len(adapter.build_feature_matrix()), len(self.records))

    def test_agent_single_and_batch_payload_shapes(self):
        sample = {
            "sample_id": "agent-sample",
            "sample_type": "experimental",
            "reagent_list": {"r1": 100, "r2": 50, "r3": 25},
            "metadata": {"protein_target": "sfGFP"},
        }
        single = CFPSSubmissionCorrelationAdapter({"sample": sample})
        self.assertEqual(single.build_feature_matrix().loc["agent-sample", "r1"], 100.0)
        batch = CFPSSubmissionCorrelationAdapter({"batch_id": "batch-agent", "samples": [sample, {**sample, "sample_id": "agent-sample-2", "reagent_list": {"r1": 0, "r2": 100, "r3": 50}}]})
        self.assertEqual(len(batch.build_feature_matrix()), 2)
        self.assertTrue((batch.build_metadata()["batch_id"] == "batch-agent").all())

    def test_adjusted_rand_index_is_label_permutation_invariant(self):
        self.assertAlmostEqual(adjusted_rand_index([1, 1, 2, 2], ["x", "x", "y", "y"]), 1.0)
        self.assertLess(adjusted_rand_index([1, 1, 2, 2], ["x", "y", "x", "y"]), 0.0)

    def test_compare_with_reports_overlap_and_ari(self):
        a = CFPSSubmissionCorrelationAdapter(self.records, correlation_method="spearman", n_clusters=2).fit()
        b = CFPSSubmissionCorrelationAdapter(self.records, correlation_method="pearson", n_clusters=2).fit()
        comparison = a.compare_with(b)
        self.assertEqual(set(comparison["group"]), {"sfGFP", "PETase"})
        self.assertTrue((comparison["n_overlap"] > 0).all())


if __name__ == "__main__":
    unittest.main()
