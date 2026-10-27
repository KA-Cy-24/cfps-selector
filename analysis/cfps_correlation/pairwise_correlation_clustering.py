"""Reusable correlation-based clustering for entity-by-feature data.

The abstract class implements shared validation, pairwise correlation, conversion to
correlation dissimilarity, hierarchical clustering, output tables, and plots. A
subclass/adapter is responsible for converting a project-specific source (e.g. CFPS
submission JSON) into a numeric feature matrix and aligned metadata.

Interpretation warning
----------------------
For CFPS designs, the rows are submissions and the columns are reagent features.
The correlation between two rows measures similarity of their across-reagent *patterns*.
It does NOT measure correlation between experimental outcomes, causal relationships,
or expected information gain.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Any, Mapping, Optional

import numpy as np
import pandas as pd


def adjusted_rand_index(labels_a: Any, labels_b: Any) -> float:
    """Compute adjusted Rand index for two partitions over the same entities.

    Labels are aligned by position; callers must align entity IDs first. The score
    is invariant to arbitrary cluster label names and is intended for sensitivity
    analysis, not for choosing a preferred clustering by itself.
    """
    a = np.asarray(labels_a, dtype=object).ravel()
    b = np.asarray(labels_b, dtype=object).ravel()
    if len(a) != len(b):
        raise ValueError("Both label arrays must have the same length")
    n = len(a)
    if n < 2:
        return float("nan")
    contingency = pd.crosstab(pd.Series(a, name="a"), pd.Series(b, name="b")).to_numpy(dtype=float)
    choose2 = lambda x: x * (x - 1.0) / 2.0
    sum_cells = float(choose2(contingency).sum())
    sum_rows = float(choose2(contingency.sum(axis=1)).sum())
    sum_cols = float(choose2(contingency.sum(axis=0)).sum())
    total_pairs = float(choose2(n))
    expected = (sum_rows * sum_cols) / total_pairs
    maximum = 0.5 * (sum_rows + sum_cols)
    denominator = maximum - expected
    if abs(denominator) < np.finfo(float).eps:
        # Degenerate partitions: score identical pairwise memberships as 1, else 0.
        same = all((a[i] == a[j]) == (b[i] == b[j]) for i in range(n) for j in range(i + 1, n))
        return 1.0 if same else 0.0
    return float((sum_cells - expected) / denominator)


@dataclass
class CorrelationClusterResult:
    """Results for a single group of entities (e.g. one protein target)."""

    group_name: str
    feature_matrix: pd.DataFrame
    correlation_matrix: pd.DataFrame
    distance_matrix: pd.DataFrame
    labels: pd.Series
    linkage_matrix: np.ndarray
    leaf_order: list[str]
    dropped_entities: list[str]
    dropped_features: list[str]
    correlation_method: str
    linkage_method: str
    feature_mode: str


class AbstractPairwiseCorrelationClusterer(ABC):
    """Reusable analysis workflow for correlation-based hierarchical clustering.

    Subclasses must implement ``build_feature_matrix`` and ``build_metadata``.
    Feature matrix rows are entities and columns are numeric features. Entity indexes
    must be unique strings and must align with the metadata index. Missing feature
    values are rejected in the abstract layer: each adapter must state how its source
    encodes absence/missingness rather than having zeros guessed silently.

    ``1 - correlation`` is used as a dissimilarity. Since it is not guaranteed to be
    Euclidean, Ward/centroid/median linkage are intentionally disallowed; average or
    complete linkage are reasonable defaults for this distance.
    """

    ALLOWED_CORRELATIONS = {"pearson", "spearman"}
    ALLOWED_LINKAGES = {"average", "complete", "single", "weighted"}

    def __init__(
        self,
        *,
        correlation_method: str = "spearman",
        linkage_method: str = "average",
        n_clusters: int = 6,
        min_features: int = 3,
        group_by: Optional[str] = "protein_target",
    ) -> None:
        self.correlation_method = correlation_method.lower()
        self.linkage_method = linkage_method.lower()
        self.n_clusters = int(n_clusters)
        self.min_features = int(min_features)
        self.group_by = group_by

        if self.correlation_method not in self.ALLOWED_CORRELATIONS:
            raise ValueError(f"correlation_method must be one of {sorted(self.ALLOWED_CORRELATIONS)}")
        if self.linkage_method not in self.ALLOWED_LINKAGES:
            raise ValueError(
                f"linkage_method must be one of {sorted(self.ALLOWED_LINKAGES)}; "
                "Ward/centroid/median are not supported for 1-r dissimilarity."
            )
        if self.n_clusters < 1:
            raise ValueError("n_clusters must be >= 1")
        if self.min_features < 2:
            raise ValueError("min_features must be >= 2")

        self.results_: dict[str, CorrelationClusterResult] = {}
        self.skipped_groups_: dict[str, str] = {}
        self.feature_matrix_: Optional[pd.DataFrame] = None
        self.metadata_: Optional[pd.DataFrame] = None

    @abstractmethod
    def build_feature_matrix(self) -> pd.DataFrame:
        """Return a numeric DataFrame: entities as rows, features as columns."""
        raise NotImplementedError

    @abstractmethod
    def build_metadata(self) -> pd.DataFrame:
        """Return entity metadata indexed by the same entity IDs as the feature matrix."""
        raise NotImplementedError

    @staticmethod
    def _safe_group_name(value: Any) -> str:
        if value is None or (isinstance(value, float) and np.isnan(value)):
            return "Unknown"
        text = str(value).strip()
        return text if text else "Unknown"

    @staticmethod
    def _slug(value: str) -> str:
        slug = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_.-")
        return slug or "group"

    def validate_inputs(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Build and validate matrices; does not mutate the adapter's source records."""
        X = self.build_feature_matrix().copy()
        metadata = self.build_metadata().copy()
        if not isinstance(X, pd.DataFrame) or X.empty:
            raise ValueError("build_feature_matrix() must return a non-empty pandas DataFrame")
        if not isinstance(metadata, pd.DataFrame):
            raise TypeError("build_metadata() must return a pandas DataFrame")

        X.index = X.index.map(str)
        metadata.index = metadata.index.map(str)
        if X.index.has_duplicates:
            raise ValueError("Feature matrix index must be unique; create stable unique entity IDs in the adapter")
        if metadata.index.has_duplicates:
            raise ValueError("Metadata index must be unique")
        if not X.index.equals(metadata.index):
            missing = sorted(set(X.index) - set(metadata.index))[:5]
            extra = sorted(set(metadata.index) - set(X.index))[:5]
            raise ValueError(f"Feature/metadata indexes must match and be ordered alike. Missing metadata={missing}; extra metadata={extra}")

        X = X.apply(pd.to_numeric, errors="coerce")
        if X.columns.duplicated().any():
            raise ValueError("Feature names must be unique")
        if X.isna().any().any():
            n = int(X.isna().sum().sum())
            raise ValueError(
                f"Feature matrix contains {n} missing/non-numeric values. Resolve missingness explicitly in the adapter; "
                "do not silently turn unknown values into zero."
            )
        values = X.to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError("Feature matrix contains infinite values")
        self.feature_matrix_ = X
        self.metadata_ = metadata
        return X, metadata

    def _fit_one_group(self, group_name: str, X: pd.DataFrame) -> Optional[CorrelationClusterResult]:
        if len(X) < 2:
            self.skipped_groups_[group_name] = f"Only {len(X)} entity/entities; need at least 2."
            return None

        # Invariant features carry no discriminatory information within this group.
        variable_mask = X.nunique(axis=0, dropna=False) > 1
        dropped_features = X.columns[~variable_mask].astype(str).tolist()
        Xv = X.loc[:, variable_mask].copy()
        if Xv.shape[1] < self.min_features:
            self.skipped_groups_[group_name] = (
                f"Only {Xv.shape[1]} variable feature(s) after removing constants; "
                f"need at least {self.min_features}."
            )
            return None

        # Correlation is undefined for an entity vector with zero variance across features.
        row_sd = Xv.std(axis=1, ddof=0)
        valid_rows = row_sd > np.finfo(float).eps
        dropped_entities = Xv.index[~valid_rows].astype(str).tolist()
        Xv = Xv.loc[valid_rows]
        if len(Xv) < 2:
            self.skipped_groups_[group_name] = "Fewer than two non-constant entity vectors remain."
            return None

        correlation = Xv.T.corr(method=self.correlation_method, min_periods=self.min_features)
        correlation = correlation.clip(lower=-1.0, upper=1.0)
        correlation = (correlation + correlation.T) / 2.0
        np.fill_diagonal(correlation.values, 1.0)
        if not np.isfinite(correlation.to_numpy(dtype=float)).all():
            raise ValueError(
                f"Group {group_name!r} produced undefined pairwise correlations. "
                "Check sparse/tied features or raise feature quality before clustering."
            )

        distance = (1.0 - correlation).clip(lower=0.0, upper=2.0)
        distance = (distance + distance.T) / 2.0
        np.fill_diagonal(distance.values, 0.0)

        from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
        from scipy.spatial.distance import squareform

        condensed = squareform(distance.to_numpy(dtype=float), checks=True)
        Z = linkage(condensed, method=self.linkage_method, optimal_ordering=True)
        k = min(self.n_clusters, len(Xv))
        cluster_numbers = fcluster(Z, t=k, criterion="maxclust")
        labels = pd.Series(
            [f"C{int(v):02d}" for v in cluster_numbers], index=Xv.index, name="cluster_id"
        )
        order = leaves_list(Z)
        leaf_order = Xv.index[order].astype(str).tolist()

        return CorrelationClusterResult(
            group_name=group_name,
            feature_matrix=Xv,
            correlation_matrix=correlation,
            distance_matrix=distance,
            labels=labels,
            linkage_matrix=Z,
            leaf_order=leaf_order,
            dropped_entities=dropped_entities,
            dropped_features=dropped_features,
            correlation_method=self.correlation_method,
            linkage_method=self.linkage_method,
            feature_mode=getattr(self, "feature_mode", "adapter-defined"),
        )

    def fit(self) -> "AbstractPairwiseCorrelationClusterer":
        """Run one clustering per group (protein target by default)."""
        X, metadata = self.validate_inputs()
        self.results_ = {}
        self.skipped_groups_ = {}

        if self.group_by and self.group_by in metadata.columns:
            groups = metadata[self.group_by].map(self._safe_group_name)
            grouped_indices = [(name, metadata.index[groups == name]) for name in pd.unique(groups)]
        elif self.group_by is None:
            grouped_indices = [("All", metadata.index)]
        else:
            raise ValueError(f"group_by={self.group_by!r} is not a metadata column; available: {list(metadata.columns)}")

        for name, ids in grouped_indices:
            result = self._fit_one_group(str(name), X.loc[ids])
            if result is not None:
                self.results_[str(name)] = result

        if not self.results_:
            raise ValueError(f"No groups could be clustered. Skip reasons: {self.skipped_groups_}")
        return self

    def labels_frame(self) -> pd.DataFrame:
        """Return row-level labels plus metadata. Treat as private unless approved for sharing."""
        if not self.results_:
            raise RuntimeError("Call fit() before requesting labels")
        assert self.metadata_ is not None
        frames = []
        for group_name, result in self.results_.items():
            part = self.metadata_.loc[result.labels.index].copy()
            part.insert(0, "entity_id", part.index.astype(str))
            part["cluster_group"] = group_name
            part["cluster_id"] = result.labels.astype(str)
            part["correlation_method"] = result.correlation_method
            part["linkage_method"] = result.linkage_method
            part["feature_mode"] = result.feature_mode
            frames.append(part.reset_index(drop=True))
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    def compare_with(self, other: "AbstractPairwiseCorrelationClusterer") -> pd.DataFrame:
        """Compare this fit with another fit using adjusted Rand index by group.

        Only entities clustered by both runs are compared; the overlap count/fraction
        is reported so that an apparently good score based on a small subset is visible.
        """
        if not self.results_ or not other.results_:
            raise RuntimeError("Both clusterers must be fitted before compare_with()")
        rows = []
        for group in sorted(set(self.results_) & set(other.results_)):
            left = self.results_[group].labels.rename("left")
            right = other.results_[group].labels.rename("right")
            joined = pd.concat([left, right], axis=1, join="inner").dropna()
            n_left = len(left)
            n_right = len(right)
            rows.append({
                "group": group,
                "n_overlap": int(len(joined)),
                "n_left": int(n_left),
                "n_right": int(n_right),
                "overlap_pct_of_left": (100.0 * len(joined) / n_left) if n_left else np.nan,
                "overlap_pct_of_right": (100.0 * len(joined) / n_right) if n_right else np.nan,
                "adjusted_rand_index": adjusted_rand_index(joined["left"], joined["right"]) if len(joined) >= 2 else np.nan,
            })
        return pd.DataFrame(rows, columns=["group", "n_overlap", "n_left", "n_right", "overlap_pct_of_left", "overlap_pct_of_right", "adjusted_rand_index"])

    def cluster_sizes_frame(self) -> pd.DataFrame:
        """Aggregate cluster sizes; review before publication because groups can be small."""
        rows = []
        for group_name, result in self.results_.items():
            counts = result.labels.value_counts().sort_index()
            for cluster_id, count in counts.items():
                rows.append({"group": group_name, "cluster_id": cluster_id, "n_entities": int(count)})
        return pd.DataFrame(rows, columns=["group", "cluster_id", "n_entities"])

    def save_private_outputs(self, output_dir: str | Path) -> list[Path]:
        """Write row-level matrices/labels to a caller-specified private directory."""
        if not self.results_:
            raise RuntimeError("Call fit() before saving outputs")
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        written: list[Path] = []
        label_path = out / "cluster_labels_private.csv"
        self.labels_frame().to_csv(label_path, index=False)
        written.append(label_path)
        for name, result in self.results_.items():
            slug = self._slug(name)
            for filename, table in (
                (f"correlation_{slug}_private.csv", result.correlation_matrix),
                (f"distance_{slug}_private.csv", result.distance_matrix),
                (f"features_{slug}_private.csv", result.feature_matrix),
            ):
                path = out / filename
                table.to_csv(path, index=True)
                written.append(path)
            path = out / f"cluster_sizes_{slug}_aggregate.csv"
            result.labels.value_counts().rename_axis("cluster_id").rename("n_entities").reset_index().to_csv(path, index=False)
            written.append(path)
        return written

    def plot_correlation_heatmap(
        self,
        group_name: str,
        *,
        annotate_max_entities: int = 35,
        figsize: Optional[tuple[float, float]] = None,
        save_path: str | Path | None = None,
    ):
        """Plot pairwise correlations, ordered by hierarchical leaves and cluster."""
        import matplotlib.pyplot as plt

        result = self.results_[group_name]
        order = result.leaf_order
        corr = result.correlation_matrix.loc[order, order]
        n = len(order)
        if figsize is None:
            side = min(16.0, max(6.0, n * 0.035))
            figsize = (side, side)
        fig, ax = plt.subplots(figsize=figsize)
        image = ax.imshow(corr.to_numpy(dtype=float), vmin=-1, vmax=1, cmap="coolwarm", interpolation="nearest", aspect="auto")
        cb = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label(f"{result.correlation_method.title()} correlation (recipe vectors)")
        ax.set_title(f"Pairwise recipe correlation — {group_name}\n{n} submissions · {result.linkage_method} linkage · {result.feature_mode} features")
        if n <= annotate_max_entities:
            ax.set_xticks(range(n), labels=order, rotation=90, fontsize=7)
            ax.set_yticks(range(n), labels=order, fontsize=7)
        else:
            ax.set_xticks([])
            ax.set_yticks([])
        # Draw separators between adjacent clusters in leaf order.
        ordered_labels = result.labels.loc[order].astype(str).tolist()
        for i in range(1, n):
            if ordered_labels[i] != ordered_labels[i - 1]:
                ax.axhline(i - 0.5, color="black", linewidth=0.7, alpha=0.7)
                ax.axvline(i - 0.5, color="black", linewidth=0.7, alpha=0.7)
        fig.tight_layout()
        if save_path is not None:
            Path(save_path).parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path, dpi=220, bbox_inches="tight")
        return fig, ax

    def plot_dendrogram(
        self,
        group_name: str,
        *,
        show_labels_max_entities: int = 60,
        figsize: Optional[tuple[float, float]] = None,
        save_path: str | Path | None = None,
    ):
        """Plot a hierarchical dendrogram; suppress labels on large groups."""
        import matplotlib.pyplot as plt
        from scipy.cluster.hierarchy import dendrogram

        result = self.results_[group_name]
        n = len(result.labels)
        if figsize is None:
            figsize = (max(10.0, min(22.0, n * 0.06)), 6.5)
        fig, ax = plt.subplots(figsize=figsize)
        dendrogram(
            result.linkage_matrix,
            labels=result.feature_matrix.index.astype(str).tolist() if n <= show_labels_max_entities else None,
            leaf_rotation=90,
            leaf_font_size=7,
            ax=ax,
            no_labels=(n > show_labels_max_entities),
            color_threshold=None,
        )
        ax.set_title(f"Hierarchical clustering — {group_name}")
        ax.set_xlabel("Submission" if n <= show_labels_max_entities else "Submissions (labels suppressed)")
        ax.set_ylabel("1 − {0} correlation dissimilarity".format(result.correlation_method))
        fig.tight_layout()
        if save_path is not None:
            Path(save_path).parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path, dpi=220, bbox_inches="tight")
        return fig, ax

    def plot_cluster_sizes(
        self,
        group_name: str,
        *,
        figsize: tuple[float, float] = (7.5, 4.5),
        save_path: str | Path | None = None,
    ):
        """Plot the number of entities in each cluster."""
        import matplotlib.pyplot as plt

        result = self.results_[group_name]
        counts = result.labels.value_counts().sort_index()
        fig, ax = plt.subplots(figsize=figsize)
        ax.bar(counts.index.astype(str), counts.values)
        ax.set_title(f"Cluster sizes — {group_name}")
        ax.set_xlabel("Cluster")
        ax.set_ylabel("Number of submissions")
        fig.tight_layout()
        if save_path is not None:
            Path(save_path).parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path, dpi=220, bbox_inches="tight")
        return fig, ax


class CFPSSubmissionCorrelationAdapter(AbstractPairwiseCorrelationClusterer):
    """Adapter from the CFPS coopetition snapshot's records to recipe feature vectors.

    Accepted input shapes:
      * raw ``cfps_designs`` record list;
      * snapshot mapping with ``collections.cfps_designs``;
      * iterable PocketBase records with ``design_json`` as a dict or JSON string;
      * a single submission payload (``{"sample": ...}``) or batch payload
        (``{"samples": [...]}``) matching the coopetition agent instructions.

    In the agent/API contract, a sample has ``sample_id``, ``sample_type``,
    ``reagent_list`` (reagent ID -> volume in nL), and ``metadata``; the live PocketBase
    collection wraps this under ``design_json``. Omitted reagent keys are treated as
    zero volume because the submission contract uses a sparse reagent mapping. A key
    present with a non-numeric value is not silently treated as zero; that row's feature
    value is marked invalid and input validation will fail loudly.

    ``feature_mode`` options:
      * ``volume``: recorded nL volumes (default);
      * ``log1p_volume``: log(1 + nL), to reduce domination by large volume values;
      * ``presence``: 1 if volume > 0, otherwise 0.
    """

    FEATURE_MODES = {"volume", "log1p_volume", "presence"}

    def __init__(
        self,
        source: Any,
        *,
        feature_mode: str = "volume",
        targets: Optional[list[str]] = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(**kwargs)
        if feature_mode not in self.FEATURE_MODES:
            raise ValueError(f"feature_mode must be one of {sorted(self.FEATURE_MODES)}")
        self.feature_mode = feature_mode
        self.targets = set(targets) if targets is not None else None
        self.records = self._extract_records(source)
        self._parsed_rows = self._parse_records(self.records)
        if self.targets is not None:
            allowed = {self._normalize_target(t) for t in self.targets}
            self._parsed_rows = [r for r in self._parsed_rows if r["protein_target"] in allowed]
        if not self._parsed_rows:
            raise ValueError("No CFPS submission records remain after parsing/target filtering")

    @staticmethod
    def _parse_json(value: Any) -> Any:
        if isinstance(value, str):
            text = value.strip()
            if not text:
                return None
            try:
                import json
                return json.loads(text)
            except (TypeError, ValueError):
                return value
        return value

    @classmethod
    def _as_dict(cls, value: Any) -> dict:
        parsed = cls._parse_json(value)
        return parsed if isinstance(parsed, dict) else {}

    @classmethod
    def _extract_records(cls, source: Any) -> list[dict]:
        inherited_batch_id = None
        if isinstance(source, Mapping):
            collections = source.get("collections", source)
            if isinstance(collections, Mapping) and "cfps_designs" in collections:
                records = collections["cfps_designs"]
            elif isinstance(source.get("samples"), (list, tuple)):
                # Agent API batch payload: {experiment_hypothesis, batch_id, samples: [...]}
                records = source["samples"]
                inherited_batch_id = source.get("batch_id") or source.get("batchId")
            elif isinstance(source.get("sample"), Mapping):
                # Agent API single-sample payload: {sample: {...}}
                records = [source]
            elif isinstance(source.get("reagent_list"), Mapping):
                # Bare sample object: {sample_id, sample_type, reagent_list, metadata}
                records = [source]
            else:
                raise ValueError(
                    "Mapping input must be a snapshot containing collections['cfps_designs'], "
                    "an agent payload with 'sample'/'samples', or a bare sample object."
                )
        else:
            records = source
        if not isinstance(records, (list, tuple)):
            raise TypeError("source must be a list of records, a snapshot mapping, or an agent sample payload")
        output = []
        for rec in records:
            if not isinstance(rec, Mapping):
                continue
            item = dict(rec)
            if inherited_batch_id is not None:
                item.setdefault("batch_id", inherited_batch_id)
            output.append(item)
        return output

    @staticmethod
    def _normalize_target(value: Any) -> str:
        if value is None or not str(value).strip():
            return "Unknown"
        raw = str(value).strip()
        known = {"sfgfp": "sfGFP", "sf gfp": "sfGFP", "petase": "PETase", "reteplase": "Reteplase"}
        return known.get(raw.casefold(), raw)

    def _parse_records(self, records: list[dict]) -> list[dict]:
        parsed_rows = []
        seen_entity_ids: set[str] = set()
        for idx, record in enumerate(records):
            design = self._as_dict(record.get("design_json"))
            # Agent POST payloads wrap a single composition under `sample`; other exports
            # may keep its fields directly on the record.
            if not design and isinstance(record.get("sample"), Mapping):
                design = self._as_dict(record.get("sample"))
            if not design and isinstance(record.get("reagent_list"), Mapping):
                design = record
            metadata = self._as_dict(design.get("metadata"))
            raw_reagents = self._as_dict(design.get("reagent_list"))
            target_value = metadata.get("protein_target", design.get("protein_target", record.get("protein_target", record.get("proteinTarget"))))
            target = self._normalize_target(target_value)
            if not raw_reagents:
                continue
            base_id = record.get("id") or design.get("sample_id") or record.get("title") or f"row_{idx:06d}"
            entity_id = str(base_id)
            if entity_id in seen_entity_ids:
                entity_id = f"{entity_id}__row_{idx:06d}"
            seen_entity_ids.add(entity_id)

            values: dict[str, float] = {}
            bad_fields: list[str] = []
            for reagent_id, raw_value in raw_reagents.items():
                reagent = str(reagent_id)
                try:
                    value = float(raw_value)
                except (TypeError, ValueError):
                    bad_fields.append(reagent)
                    continue
                if not np.isfinite(value) or value < 0:
                    bad_fields.append(reagent)
                    continue
                values[reagent] = value
            if bad_fields:
                raise ValueError(
                    f"Record {entity_id!r} has invalid/nonnegative-volume violations in reagent field(s): {bad_fields[:8]}"
                )
            if not values:
                continue
            if self.feature_mode == "presence":
                values = {k: float(v > 0) for k, v in values.items()}
            elif self.feature_mode == "log1p_volume":
                values = {k: float(np.log1p(v)) for k, v in values.items()}

            parsed_rows.append({
                "entity_id": entity_id,
                "record_id": str(record.get("id") or entity_id),
                "sample_id": design.get("sample_id") or record.get("title"),
                "protein_target": target,
                "batch_id": metadata.get("batch_id") or record.get("batch_id") or record.get("batchId"),
                "sample_type": design.get("sample_type"),
                "created": record.get("created"),
                "reagents": values,
            })
        return parsed_rows

    def build_feature_matrix(self) -> pd.DataFrame:
        # Missing mapping keys represent unadded reagents and are therefore zero in this API contract.
        feature_names = sorted({name for row in self._parsed_rows for name in row["reagents"]})
        rows = []
        ids = []
        for row in self._parsed_rows:
            rows.append({name: row["reagents"].get(name, 0.0) for name in feature_names})
            ids.append(row["entity_id"])
        return pd.DataFrame(rows, index=pd.Index(ids, name="entity_id"), columns=feature_names, dtype=float)

    def build_metadata(self) -> pd.DataFrame:
        columns = ["record_id", "sample_id", "protein_target", "batch_id", "sample_type", "created"]
        rows = []
        ids = []
        for row in self._parsed_rows:
            ids.append(row["entity_id"])
            rows.append({key: row.get(key) for key in columns})
        return pd.DataFrame(rows, index=pd.Index(ids, name="entity_id"), columns=columns)

    def feature_prevalence(self) -> pd.DataFrame:
        """Return a private QC table showing how frequently each reagent is present."""
        X = self.build_feature_matrix()
        return pd.DataFrame({
            "reagent_id": X.columns.astype(str),
            "n_entities_present": (X > 0).sum(axis=0).to_numpy(dtype=int),
            "n_entities": len(X),
            "prevalence_pct": 100.0 * (X > 0).mean(axis=0).to_numpy(dtype=float),
        }).sort_values(["prevalence_pct", "reagent_id"], ascending=[False, True]).reset_index(drop=True)
