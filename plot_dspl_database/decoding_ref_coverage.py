"""Script decoding_ref-coverage.

Phase 2 scope:
- validate script settings
- load SQL data model for all configured inputs (barcodes)
- build normalized intermediate dataframes used by future plotting phase
"""

from __future__ import annotations

from dataclasses import dataclass
import sqlite3
from typing import ClassVar

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from pydantic import BaseModel, ConfigDict, field_validator, model_validator

from optional_script_utils import (
    MetadataSection,
    PlotTextSettings,
    apply_figure_title,
    build_metadata_pages,
    make_axes_grid,
    resolve_database_path,
    save_figure_page_and_png,
)
from protocols import (
    InputSpec,
    ScriptExecutionContext,
)


SUPPORTED_DECODING_REF_COVERAGE_PLOTS: set[str] = {
    "ref-coverage-at-image-decoding_distribution",
    "ref-coverage-at-ref-decoding_distribution",
    "ref-coverage-at-ref-decoding_vs_delta-g_scatter",
    "ref-coverage-at-ref-decoding_vs_delta-g_mean-in-bin",
    "underdecoded-ref-proba_vs_delta-g_mean-in-bin",
}

DEFAULT_PLOT_TITLES: dict[str, str] = {
    "ref-coverage-at-image-decoding_distribution": "Reference Coverage distribution at image decoding",
    "ref-coverage-at-ref-decoding_distribution": "Reference Coverage distribution at reference decoding (drop out at infinity)",
    "ref-coverage-at-ref-decoding_vs_delta-g_scatter": "Reference Coverage at first perfect decoding against delta G",
    "ref-coverage-at-ref-decoding_vs_delta-g_mean-in-bin": "Binned mean positive coverage over delta G intervals",
    "underdecoded-ref-proba_vs_delta-g_mean-in-bin": "Underdecoded reference probability by delta G bin",
}

DEFAULT_PLOT_EXPLANATIONS: dict[str, str] = {
    "ref-coverage-at-image-decoding_distribution": (
        "Each bar shows the normalized number of references per image-decoding coverage bin. "
        "For each item, statistics mu_cov, sigma_cov and m_cov summarize the aggregated distribution."
    ),
    "ref-coverage-at-ref-decoding_distribution": (
        "Coverage at first perfect decoding for references decoded at least once is shown on finite bins; "
        "references never decoded are grouped at infinity."
    ),
    "ref-coverage-at-ref-decoding_vs_delta-g_scatter": (
        "Each point is a reference with X=delta G and Y=mean positive coverage at first perfect decoding; "
        "drop-out references are shown on an infinity horizontal line."
    ),
    "ref-coverage-at-ref-decoding_vs_delta-g_mean-in-bin": (
        "Delta G is binned and each bar reports mean positive coverage for bins passing the minimum point threshold."
    ),
    "underdecoded-ref-proba_vs_delta-g_mean-in-bin": (
        "For each delta G bin, the probability of underdecoded references is computed as "
        "N_underdecoded / (N_underdecoded + N_accepted)."
    ),
}


class DecodingRefCoverageSettings(BaseModel):
    model_config = ConfigDict(extra="forbid")

    use_single_run_id_instead: str | None = None


class ImageDecodingDistributionSettings(PlotTextSettings):
    coverage_display_min: int = 0
    coverage_display_max: int = 50
    max_xticks: int = 16

    @field_validator("max_xticks")
    @classmethod
    def _validate_max_xticks(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("max_xticks must be > 0.")
        return value

    @model_validator(mode="after")
    def _validate_coverage_window(self) -> "ImageDecodingDistributionSettings":
        if self.coverage_display_min > self.coverage_display_max:
            raise ValueError("coverage_display_min must be <= coverage_display_max.")
        return self


class RefDecodingDistributionSettings(PlotTextSettings):
    bin_width: float = 0.10
    max_xticks: int = 18

    @field_validator("bin_width")
    @classmethod
    def _validate_bin_width(cls, value: float) -> float:
        if value <= 0:
            raise ValueError("bin_width must be > 0.")
        return value

    @field_validator("max_xticks")
    @classmethod
    def _validate_max_xticks(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("max_xticks must be > 0.")
        return value


class RefDecodingDeltaGScatterSettings(PlotTextSettings):
    visual_infinity_factor: float = 1.08

    @field_validator("visual_infinity_factor")
    @classmethod
    def _validate_visual_infinity_factor(cls, value: float) -> float:
        if value <= 0:
            raise ValueError("visual_infinity_factor must be > 0.")
        return value


class RefDecodingDeltaGMeanInBinSettings(PlotTextSettings):
    delta_g_precision: float = 0.33333
    min_points_per_bin: int = 10

    @field_validator("delta_g_precision")
    @classmethod
    def _validate_delta_g_precision(cls, value: float) -> float:
        if value <= 0:
            raise ValueError("delta_g_precision must be > 0.")
        return value

    @field_validator("min_points_per_bin")
    @classmethod
    def _validate_min_points_per_bin(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("min_points_per_bin must be > 0.")
        return value


class UnderdecodedPlotSettings(PlotTextSettings):
    max_zero_run_ratio_for_inclusion: float = 0.5
    display_linear_regression: bool = False
    delta_g_precision: float = 0.33333
    min_points_per_bin: int = 2

    @field_validator("max_zero_run_ratio_for_inclusion")
    @classmethod
    def _validate_ratio(cls, value: float) -> float:
        if value < 0.0 or value > 1.0:
            raise ValueError("max_zero_run_ratio_for_inclusion must be in [0, 1].")
        return value

    @field_validator("delta_g_precision")
    @classmethod
    def _validate_delta_g_precision(cls, value: float) -> float:
        if value <= 0:
            raise ValueError("delta_g_precision must be > 0.")
        return value

    @field_validator("min_points_per_bin")
    @classmethod
    def _validate_min_points(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("min_points_per_bin must be > 0.")
        return value


# Graphical style config is intentionally declared at module top to keep the
# same script ergonomics as legacy plotting files.
ITEM_ID_COLOR_MAP = {
    0: "#1f77b4",  # blue
    1: "#d62728",  # red
    2: "#2ca02c",  # green
    3: "#ff7f0e",  # orange
    4: "#9467bd",  # purple
    5: "#8c564b",  # brown
}

PLOT_STYLE_MAP = {
    "cluster_size_distribution": {
        "alpha": 0.85,
    },
    "delta_g_scatter": {
        "marker": "o",
        "s": 18,
        "alpha": 0.45,
        "edgecolors": "none",
    },
}


SQL_QUERY_CLUSTER_BASE = """
WITH eligible_run_item AS (
    SELECT
        m.exp_id,
        m.dec_run_id,
        m.item_id
    FROM metrics_view_run_cluster_at_decoding m
    JOIN decoding_run_label_record l
      ON l.exp_id = m.exp_id
     AND l.dec_run_id = m.dec_run_id
    WHERE l.label = ?
      AND m.item_id IN ({item_placeholders})
      AND (? IS NULL OR m.dec_run_id = ?)
    GROUP BY m.exp_id, m.dec_run_id, m.item_id
    HAVING SUM(
        CASE
            WHEN NULLIF(TRIM(CAST(m.dec_pass_id AS TEXT)), '') IS NOT NULL THEN 1
            ELSE 0
        END
    ) > 0
)
SELECT
    m.exp_id,
    m.dec_run_id,
    m.item_id,
    m.region_id,
    m.position_id,
    m.dec_pass_id,
    m.delta_g,
    m.count_used_for_consensus,
    m.count_at_first_decoding,
    m.hamming_dist_to_ref
FROM metrics_view_run_cluster_at_decoding m
JOIN eligible_run_item e
  ON e.exp_id = m.exp_id
 AND e.dec_run_id = m.dec_run_id
 AND e.item_id = m.item_id
ORDER BY m.dec_run_id ASC, m.item_id ASC, m.region_id ASC, m.position_id ASC
"""

SQL_QUERY_LABEL_RUNS_TOTAL = """
SELECT COUNT(*) AS n_runs_total
FROM (
    SELECT DISTINCT exp_id, dec_run_id
    FROM decoding_run_label_record
    WHERE label = ?
)
"""


@dataclass
class InputReadModel:
    input_spec: InputSpec
    rows_df: pd.DataFrame
    reference_level_df: pd.DataFrame
    accepted_reference_level_df: pd.DataFrame
    underdecoded_reference_level_df: pd.DataFrame
    never_decoded_reference_level_df: pd.DataFrame
    decoded_at_least_once_reference_level_df: pd.DataFrame
    n_runs_total_labeled: int
    n_runs_displayed: int
    has_delta_g: bool


class DecodingRefCoverageScript:
    script_name: ClassVar[str] = "decoding_ref-coverage"

    def __init__(self, context: ScriptExecutionContext):
        self.context = context
        self._validate_requested_plots()
        self.settings = DecodingRefCoverageSettings.model_validate(context.script.script_settings)
        self.plot_specs_by_name = {plot.name: plot for plot in context.script.plot_settings}
        self.input_models: list[InputReadModel] = []

    def _settings_for_plot(self, plot_name: str) -> PlotTextSettings:
        raw = self.plot_specs_by_name.get(plot_name)
        raw_settings = raw.settings if raw is not None else {}
        model_by_plot: dict[str, type[PlotTextSettings]] = {
            "ref-coverage-at-image-decoding_distribution": ImageDecodingDistributionSettings,
            "ref-coverage-at-ref-decoding_distribution": RefDecodingDistributionSettings,
            "ref-coverage-at-ref-decoding_vs_delta-g_scatter": RefDecodingDeltaGScatterSettings,
            "ref-coverage-at-ref-decoding_vs_delta-g_mean-in-bin": RefDecodingDeltaGMeanInBinSettings,
            "underdecoded-ref-proba_vs_delta-g_mean-in-bin": UnderdecodedPlotSettings,
        }
        settings_model = model_by_plot[plot_name]
        return settings_model.model_validate(raw_settings)

    def _underdecoded_plot_settings(self, plot_name: str) -> UnderdecodedPlotSettings:
        settings = self._settings_for_plot(plot_name)
        if not isinstance(settings, UnderdecodedPlotSettings):
            raise TypeError(f"Expected underdecoded settings for '{plot_name}'.")
        return settings

    def _input_metadata_lines(self) -> list[str]:
        per_input = []
        for model in self.input_models:
            item_names = ", ".join(item.name for item in model.input_spec.items)
            per_input.append(
                f"- {model.input_spec.name} (read_pool_id={model.input_spec.read_pool_id}, "
                f"label={model.input_spec.decoding_run_label}, runs_studied={model.n_runs_total_labeled}, "
                f"runs_displayed={model.n_runs_displayed}, items=[{item_names}], has_delta_g={model.has_delta_g})"
            )
        return per_input

    def _apply_figure_title(self, fig: Figure, plot_name: str, plot_settings: PlotTextSettings) -> None:
        title = plot_settings.title or DEFAULT_PLOT_TITLES[plot_name]
        apply_figure_title(fig, title)

    def _build_metadata_pages(self) -> list[Figure]:
        script_metadata_lines = [
            f"Output directory: {self.context.output_path / self.script_name}",
            f"Configured figures: {len(self.context.script.plot_settings)}",
            f"Script settings: {self.settings.model_dump()}",
            *self._input_metadata_lines(),
        ]
        sections: list[MetadataSection] = []
        for index, plot_spec in enumerate(self.context.script.plot_settings, start=1):
            plot_settings = self._settings_for_plot(plot_spec.name)
            title = plot_settings.title or DEFAULT_PLOT_TITLES[plot_spec.name]
            explanation = plot_settings.description or DEFAULT_PLOT_EXPLANATIONS[plot_spec.name]
            figure_metadata = dict(plot_spec.settings)
            figure_metadata.pop("title", None)
            figure_metadata.pop("description", None)
            delta_g_plot = "delta-g" in plot_spec.name
            available_inputs = [
                model.input_spec.name
                for model in self.input_models
                if not delta_g_plot or model.has_delta_g
            ]
            sections.append(
                MetadataSection(
                    title=f"{index}. {title}",
                    description=explanation,
                    lines=(
                        f"Figure settings: {figure_metadata or '{}'}",
                        f"Inputs plotted: {', '.join(available_inputs) if available_inputs else 'none'}",
                    ),
                )
            )
        return build_metadata_pages(self.script_name, script_metadata_lines, sections)

    def _validate_requested_plots(self) -> None:
        unknown_plots = [
            plot.name
            for plot in self.context.script.plot_settings
            if plot.name not in SUPPORTED_DECODING_REF_COVERAGE_PLOTS
        ]
        if unknown_plots:
            allowed = ", ".join(sorted(SUPPORTED_DECODING_REF_COVERAGE_PLOTS))
            unknown = ", ".join(sorted(unknown_plots))
            raise ValueError(
                f"Unknown plot(s) for script '{self.script_name}': {unknown}. Allowed values: {allowed}"
            )

    def run(self) -> None:
        self.input_models = [self._build_read_model_for_input(input_spec) for input_spec in self.context.inputs]

        print(f"[{self.script_name}] phase=2 completed")
        print(
            f"[{self.script_name}] configured plots="
            + ", ".join(plot.name for plot in self.context.script.plot_settings)
        )
        for input_model in self.input_models:
            print(
                f"[{self.script_name}] input={input_model.input_spec.name}, "
                f"rows={len(input_model.rows_df)}, "
                f"refs={len(input_model.reference_level_df)}, "
                f"runs(studied={input_model.n_runs_total_labeled}, displayed={input_model.n_runs_displayed}), "
                f"has_delta_g={input_model.has_delta_g}"
            )

        plot_dispatcher = {
            "ref-coverage-at-image-decoding_distribution": self._plot_ref_coverage_at_image_decoding_distribution,
            "ref-coverage-at-ref-decoding_distribution": self._plot_ref_coverage_at_ref_decoding_distribution,
            "ref-coverage-at-ref-decoding_vs_delta-g_scatter": self._plot_ref_coverage_at_ref_decoding_vs_delta_g_scatter,
            "ref-coverage-at-ref-decoding_vs_delta-g_mean-in-bin": self._plot_ref_coverage_at_ref_decoding_vs_delta_g_mean_in_bin,
            "underdecoded-ref-proba_vs_delta-g_mean-in-bin": self._plot_underdecoded_ref_proba_vs_delta_g_mean_in_bin,
        }

        output_dir = (self.context.output_path / self.script_name).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        output_pdf = output_dir / f"{self.script_name}.pdf"

        for plot_name in SUPPORTED_DECODING_REF_COVERAGE_PLOTS:
            legacy_plot_pdf = output_dir / f"{plot_name}.pdf"
            if legacy_plot_pdf.exists():
                legacy_plot_pdf.unlink()

        with PdfPages(output_pdf) as pdf:
            for metadata_fig in self._build_metadata_pages():
                pdf.savefig(metadata_fig)
                plt.close(metadata_fig)

            for plot_spec in self.context.script.plot_settings:
                plot_fn = plot_dispatcher.get(plot_spec.name)
                if plot_fn is None:
                    raise ValueError(f"No plotting implementation for '{plot_spec.name}'.")
                print(f"[{self.script_name}] rendering plot={plot_spec.name}")
                fig = plot_fn()
                save_figure_page_and_png(
                    fig=fig,
                    plot_name=plot_spec.name,
                    output_dir=output_dir,
                    pdf=pdf,
                )

    def _build_read_model_for_input(self, input_spec: InputSpec) -> InputReadModel:
        item_ids_to_plot = [item.item_id for item in input_spec.items]
        if not item_ids_to_plot:
            raise ValueError(f"Input '{input_spec.name}' has no configured items.")

        db_path = resolve_database_path(input_spec.database)
        if not db_path.exists():
            raise FileNotFoundError(f"Database not found for input '{input_spec.name}': {db_path}")

        item_placeholders = ",".join(["?"] * len(item_ids_to_plot))
        query = SQL_QUERY_CLUSTER_BASE.format(item_placeholders=item_placeholders)
        query_params = [
            input_spec.decoding_run_label,
            *item_ids_to_plot,
            self.settings.use_single_run_id_instead,
            self.settings.use_single_run_id_instead,
        ]

        with sqlite3.connect(db_path) as conn:
            rows_df = pd.read_sql_query(query, conn, params=query_params)
            n_runs_total_labeled = int(
                conn.execute(SQL_QUERY_LABEL_RUNS_TOTAL, (input_spec.decoding_run_label,)).fetchone()[0]
            )

        if rows_df.empty:
            raise ValueError(
                f"No rows returned for input '{input_spec.name}'. "
                f"Check decoding_run_label and configured items."
            )

        rows_df = rows_df.copy()
        rows_df["delta_g"] = pd.to_numeric(rows_df["delta_g"], errors="coerce")
        rows_df["count_used_for_consensus"] = pd.to_numeric(rows_df["count_used_for_consensus"], errors="coerce")
        rows_df["count_at_first_decoding"] = pd.to_numeric(rows_df["count_at_first_decoding"], errors="coerce").fillna(0)
        rows_df["hamming_dist_to_ref"] = pd.to_numeric(rows_df["hamming_dist_to_ref"], errors="coerce")

        reference_level_df = self._build_reference_level_df(rows_df)

        accepted_reference_level_df = reference_level_df[
            reference_level_df["zero_ratio"] < 0.5
        ].copy()
        underdecoded_reference_level_df = reference_level_df[
            reference_level_df["zero_ratio"] >= 0.5
        ].copy()
        never_decoded_reference_level_df = reference_level_df[
            reference_level_df["n_runs_zero"] == reference_level_df["n_runs_total"]
        ].copy()
        decoded_at_least_once_reference_level_df = reference_level_df[
            reference_level_df["n_runs_zero"] < reference_level_df["n_runs_total"]
        ].copy()

        n_runs_displayed = int(rows_df["dec_run_id"].nunique())
        has_delta_g = bool(reference_level_df["delta_g"].notna().any())

        return InputReadModel(
            input_spec=input_spec,
            rows_df=rows_df,
            reference_level_df=reference_level_df,
            accepted_reference_level_df=accepted_reference_level_df,
            underdecoded_reference_level_df=underdecoded_reference_level_df,
            never_decoded_reference_level_df=never_decoded_reference_level_df,
            decoded_at_least_once_reference_level_df=decoded_at_least_once_reference_level_df,
            n_runs_total_labeled=n_runs_total_labeled,
            n_runs_displayed=n_runs_displayed,
            has_delta_g=has_delta_g,
        )

    @staticmethod
    def _build_reference_level_df(rows_df: pd.DataFrame) -> pd.DataFrame:
        grouped = (
            rows_df.groupby(["item_id", "region_id", "position_id"], as_index=False)
            .agg(
                delta_g=("delta_g", "first"),
                n_runs_total=("count_at_first_decoding", "size"),
                n_runs_zero=("count_at_first_decoding", lambda s: int((pd.to_numeric(s, errors="coerce") <= 0).sum())),
                mean_positive_coverage=(
                    "count_at_first_decoding",
                    lambda s: pd.to_numeric(s, errors="coerce").loc[pd.to_numeric(s, errors="coerce") > 0].mean(),
                ),
                mean_count_used_for_consensus=("count_used_for_consensus", "mean"),
                mean_hamming_dist_to_ref=("hamming_dist_to_ref", "mean"),
            )
        )

        grouped["zero_ratio"] = grouped["n_runs_zero"] / grouped["n_runs_total"]
        grouped["mean_positive_coverage"] = grouped["mean_positive_coverage"].replace([np.inf, -np.inf], np.nan)
        return grouped

    @staticmethod
    def _item_name_map(input_spec: InputSpec) -> dict[int, str]:
        return {item.item_id: item.name for item in input_spec.items}

    @staticmethod
    def _make_axes_grid(n_inputs: int, height: float = 4.8) -> tuple[Figure, list[Axes]]:
        return make_axes_grid(n_inputs, height)

    def _compute_image_decoding_distribution(
        self,
        input_model: InputReadModel,
    ) -> tuple[list[int], dict[int, np.ndarray], dict[int, tuple[float, float, float]]]:
        rows_df = input_model.rows_df.copy()
        rows_df["count_used_for_consensus"] = (
            pd.to_numeric(rows_df["count_used_for_consensus"], errors="coerce").fillna(0).astype(int)
        )

        run_item_cluster_counts = (
            rows_df.groupby(["dec_run_id", "item_id", "count_used_for_consensus"])
            .size()
            .rename("n_refs")
            .reset_index()
        )

        run_item_histograms = run_item_cluster_counts.pivot_table(
            index=["dec_run_id", "item_id"],
            columns="count_used_for_consensus",
            values="n_refs",
            fill_value=0,
            aggfunc="sum",
        )
        cluster_sizes = sorted(int(size) for size in run_item_histograms.columns.tolist())
        run_item_histograms = run_item_histograms.reindex(cluster_sizes, axis=1, fill_value=0)

        normalized_by_item: dict[int, np.ndarray] = {}
        moments_by_item: dict[int, tuple[float, float, float]] = {}
        for item_id in sorted(self._item_name_map(input_model.input_spec).keys()):
            try:
                item_histograms = run_item_histograms.xs(item_id, level="item_id")
            except KeyError:
                continue
            if isinstance(item_histograms, pd.Series):
                item_histograms = item_histograms.to_frame().T

            total_refs_item = float(item_histograms.to_numpy(dtype=float).sum())
            if total_refs_item <= 0:
                normalized_distribution = np.zeros(len(cluster_sizes), dtype=float)
            else:
                aggregated_counts = item_histograms.sum(axis=0).to_numpy(dtype=float)
                normalized_distribution = aggregated_counts / total_refs_item
            normalized_by_item[item_id] = normalized_distribution

            support = np.array(cluster_sizes, dtype=float)
            total_prob = float(normalized_distribution.sum())
            if total_prob <= 0:
                moments_by_item[item_id] = (0.0, 0.0, 0.0)
            else:
                probs = normalized_distribution / total_prob
                mu_cov = float(np.sum(support * probs))
                var_cov = float(np.sum(((support - mu_cov) ** 2) * probs))
                sigma_cov = float(np.sqrt(max(var_cov, 0.0)))
                cdf = np.cumsum(probs)
                median_idx = int(np.searchsorted(cdf, 0.5, side="left"))
                median_idx = min(max(median_idx, 0), len(support) - 1)
                moments_by_item[item_id] = (mu_cov, sigma_cov, float(support[median_idx]))
        return cluster_sizes, normalized_by_item, moments_by_item

    def _plot_ref_coverage_at_image_decoding_distribution(self) -> Figure:
        plot_name = "ref-coverage-at-image-decoding_distribution"
        plot_settings = ImageDecodingDistributionSettings.model_validate(
            self.plot_specs_by_name[plot_name].settings
        )
        fig, axes = self._make_axes_grid(len(self.input_models), height=5.0)
        for ax, input_model in zip(axes, self.input_models):
            cluster_sizes, normalized_by_item, moments_by_item = self._compute_image_decoding_distribution(input_model)
            if not cluster_sizes:
                ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
                continue

            x = np.array(cluster_sizes, dtype=float)
            mask = (
                (x >= float(plot_settings.coverage_display_min))
                & (x <= float(plot_settings.coverage_display_max))
            )
            x_plot = x[mask]
            if x_plot.size == 0:
                ax.text(0.5, 0.5, "No bins in display range", ha="center", va="center", transform=ax.transAxes)
                continue

            item_ids = sorted(normalized_by_item.keys())
            bar_width = 1.0 / max(len(item_ids), 1)
            names = self._item_name_map(input_model.input_spec)
            for idx, item_id in enumerate(item_ids):
                mu_cov, sigma_cov, m_cov = moments_by_item.get(item_id, (0.0, 0.0, 0.0))
                offsets = x_plot + (idx - (len(item_ids) - 1) / 2.0) * bar_width
                ax.bar(
                    offsets,
                    normalized_by_item[item_id][mask],
                    width=bar_width,
                    color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                    alpha=PLOT_STYLE_MAP["cluster_size_distribution"]["alpha"],
                    edgecolor="none",
                    linewidth=0.0,
                    label=f"{names.get(item_id, str(item_id))} (mu={mu_cov:.1f}, sigma={sigma_cov:.1f}, m={m_cov:.1f})",
                )

            ax.set_title(input_model.input_spec.name)
            ax.set_xlabel("Coverage")
            ax.set_ylabel("Normalized number of references")
            finite_ticks = x_plot.astype(int).tolist()
            if finite_ticks:
                tick_step = max(1, int(np.ceil(len(finite_ticks) / float(plot_settings.max_xticks))))
                reduced_ticks = finite_ticks[::tick_step]
                if reduced_ticks[-1] != finite_ticks[-1]:
                    reduced_ticks.append(finite_ticks[-1])
                ax.set_xticks(reduced_ticks)
            ax.grid(axis="y", alpha=0.25, linestyle="--")
            ax.legend(fontsize=8)

        self._apply_figure_title(fig, plot_name, plot_settings)
        return fig

    def _plot_ref_coverage_at_ref_decoding_distribution(self) -> Figure:
        plot_name = "ref-coverage-at-ref-decoding_distribution"
        plot_settings = RefDecodingDistributionSettings.model_validate(
            self.plot_specs_by_name[plot_name].settings
        )
        fig, axes = self._make_axes_grid(len(self.input_models), height=5.0)
        for ax, input_model in zip(axes, self.input_models):
            included_df = input_model.decoded_at_least_once_reference_level_df
            excluded_df = input_model.never_decoded_reference_level_df
            item_names = self._item_name_map(input_model.input_spec)
            total_refs_by_item = input_model.reference_level_df.groupby("item_id").size().to_dict()

            finite_values = included_df["mean_positive_coverage"].dropna()
            max_finite = int(np.ceil(float(finite_values.max()))) if not finite_values.empty else 1
            max_bin_edge = np.ceil(max_finite / plot_settings.bin_width) * plot_settings.bin_width
            infinity_x = float(max_bin_edge + 2.0 * plot_settings.bin_width)

            rows: list[pd.DataFrame] = []
            for item_id in sorted(item_names.keys()):
                item_included = included_df[included_df["item_id"] == item_id].copy()
                item_excluded = excluded_df[excluded_df["item_id"] == item_id].copy()
                total_refs = float(total_refs_by_item.get(item_id, 0))
                if total_refs <= 0:
                    continue
                if not item_included.empty:
                    item_included = item_included[np.isfinite(item_included["mean_positive_coverage"])].copy()
                    item_included["bin_center"] = (
                        np.floor(item_included["mean_positive_coverage"] / plot_settings.bin_width)
                        * plot_settings.bin_width
                        + (plot_settings.bin_width / 2.0)
                    )
                    binned = item_included.groupby("bin_center", as_index=False).size().rename(columns={"size": "n_refs"})
                    binned["n_refs_normalized"] = binned["n_refs"] / total_refs
                    binned["item_id"] = item_id
                    binned["kind"] = "included"
                    rows.append(binned)
                if not item_excluded.empty:
                    rows.append(
                        pd.DataFrame(
                            [{
                                "bin_center": infinity_x,
                                "n_refs": len(item_excluded),
                                "n_refs_normalized": len(item_excluded) / total_refs,
                                "item_id": item_id,
                                "kind": "excluded",
                            }]
                        )
                    )

            if not rows:
                ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
                continue

            merged = pd.concat(rows, ignore_index=True)
            item_ids = sorted(merged["item_id"].unique().tolist())
            bar_width = plot_settings.bin_width / max(len(item_ids), 1)
            for idx, item_id in enumerate(item_ids):
                data = merged[merged["item_id"] == item_id].sort_values("bin_center")
                x_pos = data["bin_center"].to_numpy(dtype=float) + (idx - (len(item_ids) - 1) / 2.0) * bar_width
                mask_in = data["kind"].eq("included").to_numpy()
                mask_out = data["kind"].eq("excluded").to_numpy()
                if mask_in.any():
                    ax.bar(
                        x_pos[mask_in],
                        data.loc[data["kind"].eq("included"), "n_refs_normalized"],
                        width=bar_width,
                        color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                        alpha=0.75,
                        edgecolor="none",
                        linewidth=0.0,
                        label=f"{item_names.get(item_id, str(item_id))} decoded",
                    )
                if mask_out.any():
                    ax.bar(
                        x_pos[mask_out],
                        data.loc[data["kind"].eq("excluded"), "n_refs_normalized"],
                        width=bar_width,
                        color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                        alpha=0.95,
                        hatch="//",
                        edgecolor="none",
                        linewidth=0.0,
                        label=f"{item_names.get(item_id, str(item_id))} drop out",
                    )
            ax.axvline(infinity_x, color="black", linewidth=1.0, linestyle=":", label="infinity=drop out")
            finite_tick_candidates = sorted(
                {
                    float(x_tick)
                    for x_tick in merged.loc[
                        merged["bin_center"] < infinity_x,
                        "bin_center",
                    ].tolist()
                }
            )
            if finite_tick_candidates:
                tick_step = max(1, int(np.ceil(len(finite_tick_candidates) / float(plot_settings.max_xticks))))
                finite_ticks = finite_tick_candidates[::tick_step]
                if finite_ticks[-1] != finite_tick_candidates[-1]:
                    finite_ticks.append(finite_tick_candidates[-1])
            else:
                finite_ticks = []
            xticks = [*finite_ticks, infinity_x]
            ax.set_xticks(xticks)
            ax.set_xticklabels(
                [
                    "∞" if np.isclose(x_tick, infinity_x)
                    else (f"{x_tick:.1f}" if not np.isclose(x_tick, round(x_tick)) else str(int(round(x_tick))))
                    for x_tick in xticks
                ]
            )
            left_xlim = min(finite_ticks) - (plot_settings.bin_width / 2.0) if finite_ticks else infinity_x - plot_settings.bin_width
            ax.set_xlim(left_xlim, infinity_x + plot_settings.bin_width)
            ax.set_title(input_model.input_spec.name)
            ax.set_xlabel("Coverage at first perfect decoding")
            ax.set_ylabel("Normalized number of references in bin")
            ax.grid(True, alpha=0.25, linestyle="--")
            ax.legend(fontsize=8)

        self._apply_figure_title(fig, plot_name, plot_settings)
        return fig

    def _plot_ref_coverage_at_ref_decoding_vs_delta_g_scatter(self) -> Figure:
        plot_name = "ref-coverage-at-ref-decoding_vs_delta-g_scatter"
        plot_settings = RefDecodingDeltaGScatterSettings.model_validate(
            self.plot_specs_by_name[plot_name].settings
        )
        fig, axes = self._make_axes_grid(len(self.input_models), height=5.0)
        for ax, input_model in zip(axes, self.input_models):
            if not input_model.has_delta_g:
                ax.text(0.5, 0.5, "No delta_g data", ha="center", va="center", transform=ax.transAxes)
                ax.set_title(input_model.input_spec.name)
                continue

            item_names = self._item_name_map(input_model.input_spec)
            accepted_df = input_model.decoded_at_least_once_reference_level_df.dropna(subset=["delta_g"])
            dropped_df = input_model.never_decoded_reference_level_df.dropna(subset=["delta_g"])

            finite_means = accepted_df["mean_positive_coverage"].dropna()
            infinity_y = float(finite_means.max()) * plot_settings.visual_infinity_factor if not finite_means.empty else 1.0
            for item_id in sorted(item_names.keys()):
                item_acc = accepted_df[accepted_df["item_id"] == item_id]
                item_drop = dropped_df[dropped_df["item_id"] == item_id]
                if not item_acc.empty:
                    ax.scatter(
                        item_acc["delta_g"],
                        item_acc["mean_positive_coverage"],
                        color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                        label=f"{item_names.get(item_id, str(item_id))} included",
                        marker=PLOT_STYLE_MAP["delta_g_scatter"]["marker"],
                        s=PLOT_STYLE_MAP["delta_g_scatter"]["s"],
                        alpha=PLOT_STYLE_MAP["delta_g_scatter"]["alpha"],
                        edgecolors=PLOT_STYLE_MAP["delta_g_scatter"]["edgecolors"],
                    )
                if not item_drop.empty:
                    ax.scatter(
                        item_drop["delta_g"],
                        np.full(len(item_drop), infinity_y),
                        color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                        label=f"{item_names.get(item_id, str(item_id))} drop out",
                        marker="x",
                        s=22,
                        alpha=0.8,
                    )

            ax.axhline(infinity_y, color="black", linewidth=1.0, linestyle=":")
            ax.set_title(input_model.input_spec.name)
            ax.set_xlabel("Delta G")
            ax.set_ylabel("Reference coverage at first decoding")
            ax.grid(True, alpha=0.25, linestyle="--")
            ax.legend(fontsize=8)

        self._apply_figure_title(fig, plot_name, plot_settings)
        return fig

    def _plot_ref_coverage_at_ref_decoding_vs_delta_g_mean_in_bin(self) -> Figure:
        plot_name = "ref-coverage-at-ref-decoding_vs_delta-g_mean-in-bin"
        plot_settings = RefDecodingDeltaGMeanInBinSettings.model_validate(
            self.plot_specs_by_name[plot_name].settings
        )
        fig, axes = self._make_axes_grid(len(self.input_models), height=5.0)
        for ax, input_model in zip(axes, self.input_models):
            if not input_model.has_delta_g:
                ax.text(0.5, 0.5, "No delta_g data", ha="center", va="center", transform=ax.transAxes)
                ax.set_title(input_model.input_spec.name)
                continue

            item_names = self._item_name_map(input_model.input_spec)
            rows: list[pd.DataFrame] = []
            for item_id in sorted(item_names.keys()):
                item_df = input_model.decoded_at_least_once_reference_level_df[
                    input_model.decoded_at_least_once_reference_level_df["item_id"] == item_id
                ].dropna(subset=["delta_g", "mean_positive_coverage"]).copy()
                if item_df.empty:
                    continue
                item_df["bin_upper"] = np.ceil(item_df["delta_g"] / plot_settings.delta_g_precision) * plot_settings.delta_g_precision
                binned = (
                    item_df.groupby("bin_upper", as_index=False)
                    .agg(
                        n_points=("delta_g", "size"),
                        mean_coverage=("mean_positive_coverage", "mean"),
                    )
                )
                binned = binned[binned["n_points"] >= plot_settings.min_points_per_bin].copy()
                if binned.empty:
                    continue
                binned["bin_center"] = binned["bin_upper"] - (plot_settings.delta_g_precision / 2.0)
                binned["item_id"] = item_id
                rows.append(binned)

            if not rows:
                ax.text(0.5, 0.5, "No bins pass threshold", ha="center", va="center", transform=ax.transAxes)
                ax.set_title(input_model.input_spec.name)
                continue

            merged = pd.concat(rows, ignore_index=True)
            item_ids = sorted(merged["item_id"].unique().tolist())
            bar_width = plot_settings.delta_g_precision / max(len(item_ids), 1)
            for idx, item_id in enumerate(item_ids):
                item_binned = merged[merged["item_id"] == item_id].sort_values("bin_center")
                x_pos = item_binned["bin_center"].to_numpy(dtype=float) + (idx - (len(item_ids) - 1) / 2.0) * bar_width
                ax.bar(
                    x_pos,
                    item_binned["mean_coverage"],
                    width=bar_width,
                    color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                    alpha=0.85,
                    edgecolor="none",
                    linewidth=0.0,
                    label=f"{item_names.get(item_id, str(item_id))}",
                )

            ax.set_title(input_model.input_spec.name)
            ax.set_xlabel("Delta G")
            ax.set_ylabel("Mean coverage in bin")
            ax.grid(True, alpha=0.25, linestyle="--")
            ax.legend(fontsize=8)

        self._apply_figure_title(fig, plot_name, plot_settings)
        return fig

    def _plot_underdecoded_ref_proba_vs_delta_g_mean_in_bin(self) -> Figure:
        plot_name = "underdecoded-ref-proba_vs_delta-g_mean-in-bin"
        plot_settings = self._underdecoded_plot_settings(plot_name)
        fig, axes = self._make_axes_grid(len(self.input_models), height=5.0)
        for ax, input_model in zip(axes, self.input_models):
            if not input_model.has_delta_g:
                ax.text(0.5, 0.5, "No delta_g data", ha="center", va="center", transform=ax.transAxes)
                ax.set_title(input_model.input_spec.name)
                continue

            item_names = self._item_name_map(input_model.input_spec)
            underdecoded_threshold = plot_settings.max_zero_run_ratio_for_inclusion
            min_points_per_bin = plot_settings.min_points_per_bin

            accepted_reference_level_df = input_model.reference_level_df[
                input_model.reference_level_df["zero_ratio"] < underdecoded_threshold
            ].copy()
            underdecoded_reference_level_df = input_model.reference_level_df[
                input_model.reference_level_df["zero_ratio"] >= underdecoded_threshold
            ].copy()

            rows: list[pd.DataFrame] = []
            for item_id in sorted(item_names.keys()):
                accepted = accepted_reference_level_df[
                    accepted_reference_level_df["item_id"] == item_id
                ].dropna(subset=["delta_g"]).copy()
                underdecoded = underdecoded_reference_level_df[
                    underdecoded_reference_level_df["item_id"] == item_id
                ].dropna(subset=["delta_g"]).copy()

                if accepted.empty and underdecoded.empty:
                    continue

                if not accepted.empty:
                    accepted["bin_upper"] = np.ceil(accepted["delta_g"] / plot_settings.delta_g_precision) * plot_settings.delta_g_precision
                    accepted_counts = accepted.groupby("bin_upper", as_index=False).size().rename(columns={"size": "n_accepted"})
                else:
                    accepted_counts = pd.DataFrame(columns=["bin_upper", "n_accepted"])

                if not underdecoded.empty:
                    underdecoded["bin_upper"] = np.ceil(underdecoded["delta_g"] / plot_settings.delta_g_precision) * plot_settings.delta_g_precision
                    underdecoded_counts = underdecoded.groupby("bin_upper", as_index=False).size().rename(columns={"size": "n_underdecoded"})
                else:
                    underdecoded_counts = pd.DataFrame(columns=["bin_upper", "n_underdecoded"])

                merged = accepted_counts.merge(underdecoded_counts, on="bin_upper", how="outer").fillna(0)
                if merged.empty:
                    continue
                merged["n_accepted"] = merged["n_accepted"].astype(int)
                merged["n_underdecoded"] = merged["n_underdecoded"].astype(int)
                merged["n_total"] = merged["n_accepted"] + merged["n_underdecoded"]
                merged = merged[merged["n_total"] >= min_points_per_bin].copy()
                if merged.empty:
                    continue
                merged["underdecoded_probability"] = merged["n_underdecoded"] / merged["n_total"]
                merged["bin_center"] = merged["bin_upper"] - (plot_settings.delta_g_precision / 2.0)
                merged["item_id"] = item_id
                rows.append(merged)

            if not rows:
                ax.text(0.5, 0.5, "No bins pass threshold", ha="center", va="center", transform=ax.transAxes)
                ax.set_title(input_model.input_spec.name)
                continue

            merged_rows = pd.concat(rows, ignore_index=True)
            item_ids = sorted(merged_rows["item_id"].unique().tolist())
            bar_width = plot_settings.delta_g_precision / max(len(item_ids), 1)
            for idx, item_id in enumerate(item_ids):
                item_df = merged_rows[merged_rows["item_id"] == item_id].sort_values("bin_center")
                x_pos = item_df["bin_center"].to_numpy(dtype=float) + (idx - (len(item_ids) - 1) / 2.0) * bar_width
                ax.bar(
                    x_pos,
                    item_df["underdecoded_probability"],
                    width=bar_width,
                    color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                    alpha=0.9,
                    edgecolor="none",
                    linewidth=0.0,
                    label=f"{item_names.get(item_id, str(item_id))}",
                )

                if plot_settings.display_linear_regression and len(item_df) >= 2:
                    x_data = item_df["bin_center"].to_numpy(dtype=float)
                    y_data = item_df["underdecoded_probability"].to_numpy(dtype=float)
                    if np.std(x_data) > 0:
                        slope, intercept = np.polyfit(x_data, y_data, deg=1)
                        x_line = np.linspace(float(np.min(x_data)), float(np.max(x_data)), 200)
                        y_line = (slope * x_line) + intercept
                        ax.plot(
                            x_line,
                            y_line,
                            color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                            linestyle="-",
                            linewidth=2.0,
                            alpha=0.95,
                        )

            ax.set_title(input_model.input_spec.name)
            ax.set_xlabel("Delta G")
            ax.set_ylabel("Underdecoded reference probability")
            ax.grid(True, alpha=0.25, linestyle="--")
            ax.legend(fontsize=8)

        self._apply_figure_title(fig, plot_name, plot_settings)
        return fig
