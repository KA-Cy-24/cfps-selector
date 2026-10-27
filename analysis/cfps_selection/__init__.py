"""Pre-run CFPS experiment-selection audit, coverage, and candidate selectors."""
from .selection_framework import (
    AuditResult, audit_gp_model_snapshot, audit_gp_model_link_matrices, changed_reagents, canonical_recipe,
    read_gp_model_file, coverage_metrics, suggest_portfolio, write_selection_overlay,
)
