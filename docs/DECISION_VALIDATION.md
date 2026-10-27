# Decision framework validation protocol

The following tests are automated and reproducible in the release checkout:

```bash
python -m pip install -r requirements-analysis.txt
PYTHONPATH=. python -m pytest -q tests
python -m py_compile analysis/cfps_selection/decision_framework.py tools/cfps_decision_cli.py
node --check repo_additions/static/cfps/decision/decision.js
```

Additional release checks:

1. Load the full authorized frozen GP-annotated snapshot in a **private runtime**. Confirm exactly 86/79/54 source groups and no malformed group links, and explicitly report that independent-cluster equality is not established.
2. Generate 80 unique selections per target at 4 replicates. Count 320 assigned experimental positions per target, no assignment in columns 1–4, and 40+40 unique candidates at 50%.
3. Compare the JS selector on the exact same generated private input with the Python engine at the default settings; compare immutable **selected-ID sets**, not display order.
4. Run the notebook cells headlessly against the same frozen input and verify that the interactive private HTML preview and ZIP are produced.
5. Independently inspect the current `htgaa-automation` authentication/session paths before production deployment; run the SvelteKit build in **that** repo. The standalone viewer and JS computational path do not alone verify production route integration.
6. Review every exported batch/group warning before laboratory submission. The displayed deterministic well order is provisional and should not be used as the final executable plate file.
7. Scan the source tree for hard-coded credentials. No confidential-client secret should appear in the source. The PB password and token are requested when needed rather than bundled.
