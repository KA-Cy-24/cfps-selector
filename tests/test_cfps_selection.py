"""Synthetic-data tests: no private/source dataset or network required."""
import copy
import unittest

from analysis.cfps_selection.selection_framework import (
    audit_gp_model_snapshot, canonical_recipe, changed_reagents,
    coverage_metrics, suggest_portfolio, audit_gp_model_link_matrices,
)


def sample_data():
    def d(rid, recipe, role, group=None):
        total = sum(recipe.values())
        return {
            "id": rid, "sampleId": rid, "author": "synthetic", "protein": "sfGFP",
            "recipeNl": {**recipe, "nuclease_free_water": 20_000-total},
            "round1Role": role, "sameRecipeAs": None,
            "group": {"letter":group} if group else None,
        }
    return {"built":"synthetic-test", "proteins": {"sfGFP":{
        "designs":[
            d("a",{"magnesium_glutamate":100,"dsbc_ecoli":100},"first run","A"),
            d("b",{"magnesium_glutamate":200,"dsbc_ecoli":100},"check run","A"),
            d("c",{"magnesium_glutamate":300,"dsbc_ecoli":100},"follows its group's runs","A"),
            d("d",{"magnesium_glutamate":300,"dsbc_ecoli":0},"independent run"),
            d("e",{"magnesium_glutamate":300,"dsbc_ecoli":100},"follows its group's runs","A"),
        ],
        "groups":[{"letter":"A","firstRun":"a","linkFromSharedUntested":.1,
                   "checkRuns":[{"id":"b","ingredient":"Magnesium Glutamate", "newInformation":44}]}],
        "sweeps":[{"reagent":"magnesium_glutamate","recommended":True,
                   "sweeps":[{"toRun":[{"id":"a"},{"id":"b"},{"id":"c"}]}]}],
        "totals":{}}}}


class SelectionTests(unittest.TestCase):
    def test_canonical_ignores_balancing_water_and_zero(self):
        self.assertEqual(canonical_recipe({"a":100,"nuclease_free_water":19900}),
                         canonical_recipe({"a":100,"b":0,"nuclease_free_water":1}))

    def test_one_factor_diff_and_sweep(self):
        data = sample_data()
        result = audit_gp_model_snapshot(data,working_quota=4)
        s = result.summary["proteins"]["sfGFP"]
        self.assertEqual(s["candidateDesigns"],5)
        self.assertEqual(s["distinctFormulations"],4)
        self.assertEqual(s["selectedUniqueFormulations"],3)
        self.assertEqual(s["exactOneFactorChecksVsFirst"],1)
        self.assertEqual(s["newExactCounterpartsForChecks"],0)
        self.assertEqual(s["cleanSweepFamilies"],1)
        self.assertEqual(s["sweepFamiliesWithThreeOrMoreLevels"],0)
        self.assertTrue(result.additional_contrasts.empty)

    def test_modes(self):
        data=sample_data()
        rows=data["proteins"]["sfGFP"]["designs"]
        baseline=suggest_portfolio(rows,quota=4,mode="gp_model_baseline")
        self.assertEqual(set(baseline["ids"]), {"a","b","d"})
        hybrid=suggest_portfolio(rows,quota=4,mode="hybrid",objective="facility")
        self.assertEqual(hybrid["count"],4)
        self.assertEqual(len(set(hybrid["ids"])),4)
        model_free=suggest_portfolio(rows,quota=3,mode="model_free",objective="maximin")
        self.assertEqual(model_free["count"],3)
        metrics=coverage_metrics(rows,set(model_free["ids"]))
        self.assertGreaterEqual(metrics["mean_candidate_to_selected_distance"],0)

    def test_fail_closed_volume(self):
        with self.assertRaises(ValueError):
            canonical_recipe({"a":-1})
        with self.assertRaises(ValueError):
            canonical_recipe({"a":"not-a-number"})

    def test_order_independent(self):
        a=sample_data(); b=copy.deepcopy(a)
        b["proteins"]["sfGFP"]["designs"].reverse()
        aa=audit_gp_model_snapshot(a,working_quota=4)
        bb=audit_gp_model_snapshot(b,working_quota=4)
        self.assertEqual(aa.summary,bb.summary)
        self.assertEqual(aa.additional_contrasts.to_dict("records"),bb.additional_contrasts.to_dict("records"))

    def test_exact_unselected_counterpart_for_confounded_check(self):
        d=sample_data()
        p=d["proteins"]["sfGFP"]
        p["designs"][0]["recipeNl"] = {"magnesium_glutamate":100,"nuclease_free_water":19900}
        p["designs"][1]["recipeNl"] = {"magnesium_glutamate":200,"dsbc_ecoli":100,"nuclease_free_water":19700}
        p["sweeps"] = []
        audit=audit_gp_model_snapshot(d,working_quota=4)
        self.assertEqual(audit.summary["proteins"]["sfGFP"]["multifactorChecksVsFirst"],1)
        self.assertEqual(audit.summary["proteins"]["sfGFP"]["newExactCounterpartsForChecks"],1)
        self.assertEqual(audit.additional_contrasts.iloc[0]["candidateId"],"c")

    def test_model_free_selection_does_not_look_at_forecasts_or_roles(self):
        d=sample_data()["proteins"]["sfGFP"]["designs"]
        before=suggest_portfolio(d,quota=3,mode="model_free")["ids"]
        changed=copy.deepcopy(d)
        for item in changed:
            item["round1Role"]="first run"
            item["forecast"]={"pBeat":float(len(item["id"])),"lnSpecCost":-1000}
        after=suggest_portfolio(changed,quota=3,mode="model_free")["ids"]
        self.assertEqual(before,after)

    def test_display_matrix_psd_warning(self):
        data={"tips":{"sfgfp":{"C":[[1,1.01],[1.01,1]],"names":["a","b"]}}}
        report=audit_gp_model_link_matrices(data)
        self.assertFalse(report["targets"]["sfGFP"]["PSD_as_stored_tolerance_1e_6"])

    def test_changed_water_only_is_empty(self):
        self.assertEqual(changed_reagents({"a":10,"nuclease_free_water":20},
                                          {"a":10,"nuclease_free_water":200}),[])


if __name__ == "__main__":
    unittest.main()
