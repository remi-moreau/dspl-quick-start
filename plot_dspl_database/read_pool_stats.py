"""Read-pool reference statistics plots backed by persisted database snapshots."""

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
from pydantic import BaseModel, ConfigDict, field_validator
from scipy.stats import spearmanr

from optional_script_utils import (
    MetadataSection,
    PlotTextSettings,
    apply_figure_title,
    build_metadata_pages,
    harmonize_axes_scales,
    make_axes_grid,
    resolve_database_path,
    save_figure_page_and_png,
)
from protocols import InputSpec, ScriptExecutionContext


SUPPORTED_READ_POOL_STATS_PLOTS: set[str] = {
    "ref-coverage-in-read-pool_distribution",
    "ref-coverage-in-read-pool_vs_delta-g_scatter",
    "ref-coverage-in-read-pool_vs_delta-g_mean-in-bin",
}

DEFAULT_PLOT_TITLES: dict[str, str] = {
    "ref-coverage-in-read-pool_distribution": "Reference count distribution over normalized read-pool coverage",
    "ref-coverage-in-read-pool_vs_delta-g_scatter": "Normalized read-pool coverage against delta G",
    "ref-coverage-in-read-pool_vs_delta-g_mean-in-bin": "Binned mean normalized read-pool coverage over delta G intervals",
}

DEFAULT_PLOT_EXPLANATIONS: dict[str, str] = {
    "ref-coverage-in-read-pool_distribution": (
        "References are assigned to common normalized-coverage bins. Bars distinguish configured items, "
        "while every reference is normalized by the mean coverage across all selected items."
    ),
    "ref-coverage-in-read-pool_vs_delta-g_scatter": (
        "Each point is one reference with X=delta G and Y=read count divided by the global mean read count. "
        "Items are distinguished visually; regression and Spearman statistics are computed globally."
    ),
    "ref-coverage-in-read-pool_vs_delta-g_mean-in-bin": (
        "Delta G is binned without item distinction. Each bar reports the global mean normalized coverage "
        "for bins passing the minimum point threshold."
    ),
}

ITEM_ID_COLOR_MAP = {
    0: "#1f77b4",
    1: "#d62728",
    2: "#2ca02c",
    3: "#ff7f0e",
    4: "#9467bd",
    5: "#8c564b",
}

SUPPORTED_SNAPSHOT_VERSION = 1

SQL_QUERY_POOL_SNAPSHOT = """
SELECT
    exp_id,
    read_pool_id,
    enc_run_id,
    n_reads_total,
    n_reads_mapped,
    n_reads_unmapped,
    n_fastq_total,
    source_max_position_exclusive,
    refreshed_at,
    snapshot_version
FROM read_pool_stats_snapshot
WHERE exp_id = ?
  AND read_pool_id = ?
"""

SQL_QUERY_REFERENCE_STATS = """
SELECT
    exp_id,
    read_pool_id,
    enc_run_id,
    item_id,
    region_id,
    position_id,
    delta_g,
    n_reads_on_reference
FROM metrics_view_pool_reference_stats
WHERE exp_id = ?
  AND read_pool_id = ?
  AND item_id IN ({item_placeholders})
ORDER BY item_id, region_id, position_id
"""


class ReadPoolStatsSettings(BaseModel):
    model_config = ConfigDict(extra="forbid")


class DistributionPlotSettings(PlotTextSettings):
    read_count_values_per_bin: int | None = None

    @field_validator("read_count_values_per_bin")
    @classmethod
    def _validate_read_count_values_per_bin(cls, value: int | None) -> int | None:
        if value is not None and value <= 0:
            raise ValueError("read_count_values_per_bin must be > 0 or null.")
        return value


class ScatterPlotSettings(PlotTextSettings):
    display_linear_regression: bool = True


class MeanInBinPlotSettings(PlotTextSettings):
    delta_g_precision: float = 1.0
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


@dataclass(frozen=True)
class PoolSnapshot:
    exp_id: str
    read_pool_id: str
    enc_run_id: str
    n_reads_total: int
    n_reads_mapped: int
    n_reads_unmapped: int
    n_fastq_total: int
    source_max_position_exclusive: int
    refreshed_at: str
    snapshot_version: int


@dataclass(frozen=True)
class CorrelationStatistics:
    n_points: int
    rho: float | None
    p_value: float | None


@dataclass(frozen=True)
class RegressionStatistics:
    n_points: int
    slope: float | None
    intercept: float | None
    r_squared: float | None


@dataclass
class InputReadModel:
    input_spec: InputSpec
    snapshot: PoolSnapshot
    reference_df: pd.DataFrame
    global_mean_coverage: float
    has_delta_g: bool


def normalize_reference_coverage(reference_df: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    if reference_df.empty:
        raise ValueError("Cannot normalize an empty reference dataset.")

    normalized_df = reference_df.copy()
    coverage = pd.to_numeric(normalized_df["n_reads_on_reference"], errors="coerce")
    if coverage.isna().any():
        raise ValueError("n_reads_on_reference contains non-numeric or NULL values.")
    if (coverage < 0).any():
        raise ValueError("n_reads_on_reference contains negative values.")

    global_mean = float(coverage.mean())
    if not np.isfinite(global_mean) or global_mean <= 0:
        raise ValueError("Global mean reference coverage is zero; normalized coverage is undefined.")

    normalized_df["n_reads_on_reference"] = coverage
    normalized_df["normalized_read_pool_coverage"] = coverage / global_mean
    return normalized_df, global_mean


def compute_distribution_bin_edges(
    reference_df: pd.DataFrame,
    global_mean_coverage: float,
    read_count_values_per_bin: int | None,
) -> np.ndarray:
    normalized_values = reference_df["normalized_read_pool_coverage"].to_numpy(dtype=float)
    if read_count_values_per_bin is None:
        bins: str | int = "fd" if len(normalized_values) >= 2 else 10
        return np.histogram_bin_edges(normalized_values, bins=bins)

    read_counts = pd.to_numeric(reference_df["n_reads_on_reference"], errors="coerce")
    if read_counts.isna().any() or (read_counts < 0).any():
        raise ValueError("Cannot build distribution bins from invalid read counts.")
    if not np.allclose(read_counts, np.round(read_counts)):
        raise ValueError("read_count_values_per_bin requires integer read counts.")
    if global_mean_coverage <= 0:
        raise ValueError("global_mean_coverage must be > 0.")

    bin_width = int(read_count_values_per_bin)
    minimum_count = int(read_counts.min())
    maximum_count = int(read_counts.max())
    first_edge = (minimum_count // bin_width) * bin_width
    last_edge_exclusive = ((maximum_count // bin_width) + 1) * bin_width
    raw_edges = np.arange(first_edge, last_edge_exclusive + bin_width, bin_width, dtype=float)
    return raw_edges / global_mean_coverage


def compute_spearman_statistics(reference_df: pd.DataFrame) -> CorrelationStatistics:
    pairs = _finite_delta_g_coverage_pairs(reference_df)
    n_points = len(pairs)
    if (
        n_points < 2
        or pairs["delta_g"].nunique() < 2
        or pairs["normalized_read_pool_coverage"].nunique() < 2
    ):
        return CorrelationStatistics(n_points, None, None)

    result = spearmanr(pairs["delta_g"], pairs["normalized_read_pool_coverage"])
    rho = float(result.statistic)
    p_value = float(result.pvalue)
    if not np.isfinite(rho) or not np.isfinite(p_value):
        return CorrelationStatistics(n_points, None, None)
    return CorrelationStatistics(n_points, rho, p_value)


def compute_regression_statistics(reference_df: pd.DataFrame) -> RegressionStatistics:
    pairs = _finite_delta_g_coverage_pairs(reference_df)
    n_points = len(pairs)
    if n_points < 2 or pairs["delta_g"].nunique() < 2:
        return RegressionStatistics(n_points, None, None, None)

    x_values = pairs["delta_g"].to_numpy(dtype=float)
    y_values = pairs["normalized_read_pool_coverage"].to_numpy(dtype=float)
    slope, intercept = np.polyfit(x_values, y_values, deg=1)
    predicted = slope * x_values + intercept
    residual_sum = float(np.sum((y_values - predicted) ** 2))
    total_sum = float(np.sum((y_values - float(np.mean(y_values))) ** 2))
    r_squared = 1.0 - residual_sum / total_sum if total_sum > 0 else 1.0
    return RegressionStatistics(
        n_points=n_points,
        slope=float(slope),
        intercept=float(intercept),
        r_squared=float(r_squared),
    )


def compute_global_delta_g_bins(
    reference_df: pd.DataFrame,
    delta_g_precision: float,
    min_points_per_bin: int,
) -> pd.DataFrame:
    pairs = _finite_delta_g_coverage_pairs(reference_df).copy()
    if pairs.empty:
        return pd.DataFrame(columns=["bin_upper", "bin_center", "n_points", "mean_normalized_coverage"])

    pairs["bin_upper"] = np.ceil(pairs["delta_g"] / delta_g_precision) * delta_g_precision
    binned = (
        pairs.groupby("bin_upper", as_index=False)
        .agg(
            n_points=("normalized_read_pool_coverage", "size"),
            mean_normalized_coverage=("normalized_read_pool_coverage", "mean"),
        )
    )
    binned = binned[binned["n_points"] >= min_points_per_bin].copy()
    binned["bin_center"] = binned["bin_upper"] - delta_g_precision / 2.0
    return binned.sort_values("bin_upper").reset_index(drop=True)


def _finite_delta_g_coverage_pairs(reference_df: pd.DataFrame) -> pd.DataFrame:
    pairs = reference_df[["delta_g", "normalized_read_pool_coverage"]].copy()
    pairs["delta_g"] = pd.to_numeric(pairs["delta_g"], errors="coerce")
    pairs["normalized_read_pool_coverage"] = pd.to_numeric(
        pairs["normalized_read_pool_coverage"], errors="coerce"
    )
    return pairs.replace([np.inf, -np.inf], np.nan).dropna()


class ReadPoolStatsScript:
    script_name: ClassVar[str] = "read_pool_stats"

    def __init__(self, context: ScriptExecutionContext):
        self.context = context
        self._validate_requested_plots()
        self.settings = ReadPoolStatsSettings.model_validate(context.script.script_settings)
        self.plot_specs_by_name = {plot.name: plot for plot in context.script.plot_settings}
        self.input_models: list[InputReadModel] = []

    def _validate_requested_plots(self) -> None:
        unknown_plots = sorted(
            plot.name
            for plot in self.context.script.plot_settings
            if plot.name not in SUPPORTED_READ_POOL_STATS_PLOTS
        )
        if unknown_plots:
            raise ValueError(
                f"Unknown plot(s) for script '{self.script_name}': {', '.join(unknown_plots)}. "
                f"Allowed values: {', '.join(sorted(SUPPORTED_READ_POOL_STATS_PLOTS))}"
            )

    def _settings_for_plot(self, plot_name: str) -> PlotTextSettings:
        settings_types: dict[str, type[PlotTextSettings]] = {
            "ref-coverage-in-read-pool_distribution": DistributionPlotSettings,
            "ref-coverage-in-read-pool_vs_delta-g_scatter": ScatterPlotSettings,
            "ref-coverage-in-read-pool_vs_delta-g_mean-in-bin": MeanInBinPlotSettings,
        }
        return settings_types[plot_name].model_validate(self.plot_specs_by_name[plot_name].settings)

    def run(self) -> None:
        self.input_models = [self._build_read_model_for_input(input_spec) for input_spec in self.context.inputs]
        output_dir = (self.context.output_path / self.script_name).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        output_pdf = output_dir / f"{self.script_name}.pdf"
        plot_dispatcher = {
            "ref-coverage-in-read-pool_distribution": self._plot_distribution,
            "ref-coverage-in-read-pool_vs_delta-g_scatter": self._plot_scatter,
            "ref-coverage-in-read-pool_vs_delta-g_mean-in-bin": self._plot_mean_in_bin,
        }

        with PdfPages(output_pdf) as pdf:
            for metadata_fig in self._build_metadata_pages():
                pdf.savefig(metadata_fig)
                plt.close(metadata_fig)
            for plot_spec in self.context.script.plot_settings:
                print(f"[{self.script_name}] rendering plot={plot_spec.name}")
                figure = plot_dispatcher[plot_spec.name]()
                save_figure_page_and_png(
                    figure,
                    plot_spec.name,
                    output_dir,
                    pdf,
                    show_figure=self.context.show_figures,
                )

    def _build_read_model_for_input(self, input_spec: InputSpec) -> InputReadModel:
        item_ids = [item.item_id for item in input_spec.items]
        if not item_ids:
            raise ValueError(f"Input '{input_spec.name}' has no configured items.")

        database_path = resolve_database_path(input_spec.database)
        if not database_path.exists():
            raise FileNotFoundError(f"Database not found for input '{input_spec.name}': {database_path}")

        placeholders = ",".join(["?"] * len(item_ids))
        connection_uri = f"file:{database_path}?mode=ro"
        try:
            with sqlite3.connect(connection_uri, uri=True) as connection:
                connection.row_factory = sqlite3.Row
                snapshot_rows = connection.execute(
                    SQL_QUERY_POOL_SNAPSHOT,
                    (input_spec.exp_id, input_spec.read_pool_id),
                ).fetchall()
                reference_df = pd.read_sql_query(
                    SQL_QUERY_REFERENCE_STATS.format(item_placeholders=placeholders),
                    connection,
                    params=[input_spec.exp_id, input_spec.read_pool_id, *item_ids],
                )
        except (sqlite3.OperationalError, pd.errors.DatabaseError) as error:
            raise ValueError(
                f"Read-pool statistics schema is unavailable for input '{input_spec.name}'. "
                "Update the database schema and refresh its read-pool statistics snapshot first. "
                f"SQLite reported: {error}"
            ) from error

        if len(snapshot_rows) != 1:
            raise ValueError(
                f"No read-pool statistics snapshot found for input '{input_spec.name}' "
                f"(exp_id={input_spec.exp_id}, read_pool_id={input_spec.read_pool_id}). "
                "Update the database and refresh this read-pool statistics snapshot first."
            )

        snapshot = self._snapshot_from_row(snapshot_rows[0])
        if snapshot.snapshot_version != SUPPORTED_SNAPSHOT_VERSION:
            raise ValueError(
                f"Unsupported read-pool snapshot version for input '{input_spec.name}': "
                f"expected {SUPPORTED_SNAPSHOT_VERSION}, found {snapshot.snapshot_version}. "
                "Update the database and refresh the snapshot."
            )
        if not snapshot.enc_run_id:
            raise ValueError(f"Read pool for input '{input_spec.name}' has no enc_run_id in its snapshot.")
        if reference_df.empty:
            raise ValueError(
                f"No reference statistics found for input '{input_spec.name}' and its configured items. "
                "Check exp_id, read_pool_id and items; refresh the database snapshot if needed."
            )

        enc_run_ids = set(reference_df["enc_run_id"].dropna().astype(str).unique())
        if enc_run_ids != {snapshot.enc_run_id}:
            raise ValueError(
                f"Inconsistent enc_run_id values for input '{input_spec.name}': "
                f"snapshot={snapshot.enc_run_id}, view={sorted(enc_run_ids)}."
            )

        reference_df["delta_g"] = pd.to_numeric(reference_df["delta_g"], errors="coerce")
        reference_df, global_mean = normalize_reference_coverage(reference_df)
        return InputReadModel(
            input_spec=input_spec,
            snapshot=snapshot,
            reference_df=reference_df,
            global_mean_coverage=global_mean,
            has_delta_g=bool(reference_df["delta_g"].notna().any()),
        )

    @staticmethod
    def _snapshot_from_row(row: sqlite3.Row) -> PoolSnapshot:
        return PoolSnapshot(
            exp_id=str(row["exp_id"]),
            read_pool_id=str(row["read_pool_id"]),
            enc_run_id=str(row["enc_run_id"] or ""),
            n_reads_total=int(row["n_reads_total"]),
            n_reads_mapped=int(row["n_reads_mapped"]),
            n_reads_unmapped=int(row["n_reads_unmapped"]),
            n_fastq_total=int(row["n_fastq_total"]),
            source_max_position_exclusive=int(row["source_max_position_exclusive"]),
            refreshed_at=str(row["refreshed_at"]),
            snapshot_version=int(row["snapshot_version"]),
        )

    def _plot_distribution(self) -> Figure:
        plot_name = "ref-coverage-in-read-pool_distribution"
        plot_settings = DistributionPlotSettings.model_validate(self.plot_specs_by_name[plot_name].settings)
        figure, axes = self._make_axes_grid()

        for axis, input_model in zip(axes, self.input_models):
            bin_edges = compute_distribution_bin_edges(
                input_model.reference_df,
                input_model.global_mean_coverage,
                plot_settings.read_count_values_per_bin,
            )
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
            bin_widths = np.diff(bin_edges)
            item_names = self._item_name_map(input_model.input_spec)
            item_ids = [item_id for item_id in item_names if (input_model.reference_df["item_id"] == item_id).any()]
            slot_widths = bin_widths * 0.88 / max(len(item_ids), 1)

            for index, item_id in enumerate(item_ids):
                item_values = input_model.reference_df.loc[
                    input_model.reference_df["item_id"] == item_id,
                    "normalized_read_pool_coverage",
                ].to_numpy(dtype=float)
                counts, _ = np.histogram(item_values, bins=bin_edges)
                offsets = (index - (len(item_ids) - 1) / 2.0) * slot_widths
                axis.bar(
                    bin_centers + offsets,
                    counts,
                    width=slot_widths,
                    color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                    alpha=0.82,
                    edgecolor="black",
                    linewidth=0.35,
                    label=f"{item_names[item_id]} (n={len(item_values)})",
                )

            axis.set_title(input_model.input_spec.name)
            axis.set_xlabel("Normalized read-pool coverage")
            axis.set_ylabel("Number of references")
            axis.grid(axis="y", alpha=0.25, linestyle="--")
            axis.legend(fontsize=8)

        self._finish_figure(figure, axes, plot_name, plot_settings)
        return figure

    def _plot_scatter(self) -> Figure:
        plot_name = "ref-coverage-in-read-pool_vs_delta-g_scatter"
        plot_settings = ScatterPlotSettings.model_validate(self.plot_specs_by_name[plot_name].settings)
        figure, axes = self._make_axes_grid()

        for axis, input_model in zip(axes, self.input_models):
            rows = input_model.reference_df.dropna(subset=["delta_g"])
            item_names = self._item_name_map(input_model.input_spec)
            for item_id, item_name in item_names.items():
                item_rows = rows[rows["item_id"] == item_id]
                if item_rows.empty:
                    continue
                axis.scatter(
                    item_rows["delta_g"],
                    item_rows["normalized_read_pool_coverage"],
                    color=ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f"),
                    marker="o",
                    s=18,
                    alpha=0.5,
                    edgecolors="none",
                    label=f"{item_name} (n={len(item_rows)})",
                )

            regression = compute_regression_statistics(rows)
            if plot_settings.display_linear_regression and regression.slope is not None and regression.intercept is not None:
                x_values = rows["delta_g"].to_numpy(dtype=float)
                x_line = np.linspace(float(np.min(x_values)), float(np.max(x_values)), 200)
                axis.plot(
                    x_line,
                    regression.slope * x_line + regression.intercept,
                    color="black",
                    linewidth=2.4,
                    label="Global linear regression",
                )

            axis.set_title(input_model.input_spec.name)
            axis.set_xlabel("Delta G")
            axis.set_ylabel("Normalized read-pool coverage")
            axis.grid(True, alpha=0.25, linestyle="--")
            if axis.has_data():
                axis.legend(fontsize=8)
            else:
                axis.text(0.5, 0.5, "No delta_g data", ha="center", va="center", transform=axis.transAxes)

        self._finish_figure(figure, axes, plot_name, plot_settings)
        return figure

    def _plot_mean_in_bin(self) -> Figure:
        plot_name = "ref-coverage-in-read-pool_vs_delta-g_mean-in-bin"
        plot_settings = MeanInBinPlotSettings.model_validate(self.plot_specs_by_name[plot_name].settings)
        figure, axes = self._make_axes_grid()

        for axis, input_model in zip(axes, self.input_models):
            binned = compute_global_delta_g_bins(
                input_model.reference_df,
                plot_settings.delta_g_precision,
                plot_settings.min_points_per_bin,
            )
            if binned.empty:
                message = "No delta_g data" if not input_model.has_delta_g else "No bins pass threshold"
                axis.text(0.5, 0.5, message, ha="center", va="center", transform=axis.transAxes)
            else:
                axis.bar(
                    binned["bin_center"],
                    binned["mean_normalized_coverage"],
                    width=plot_settings.delta_g_precision * 0.92,
                    color="#4c4c4c",
                    alpha=0.9,
                    edgecolor="black",
                    linewidth=0.4,
                )
            axis.set_title(input_model.input_spec.name)
            axis.set_xlabel("Delta G")
            axis.set_ylabel("Mean normalized read-pool coverage in bin")
            axis.grid(axis="y", alpha=0.25, linestyle="--")

        self._finish_figure(figure, axes, plot_name, plot_settings)
        return figure

    def _finish_figure(
        self,
        figure: Figure,
        axes: list,
        plot_name: str,
        plot_settings: PlotTextSettings,
    ) -> None:
        harmonize_axes_scales(
            axes,
            same_x_scale=plot_settings.same_x_scale_across_inputs,
            same_y_scale=plot_settings.same_y_scale_across_inputs,
        )
        apply_figure_title(figure, plot_settings.title or DEFAULT_PLOT_TITLES[plot_name])

    def _make_axes_grid(self) -> tuple[Figure, list[Axes]]:
        return make_axes_grid(
            len(self.input_models),
            height=self.context.figure_height,
            width_per_input=self.context.figure_width_per_input,
        )

    def _build_metadata_pages(self) -> list[Figure]:
        script_lines = [
            f"Output directory: {self.context.output_path / self.script_name}",
            f"Configured figures: {len(self.context.script.plot_settings)}",
            f"Script settings: {self.settings.model_dump()}",
            "Normalization: global mean across every reference in the configured items, including zero coverage",
        ]
        for model in self.input_models:
            snapshot = model.snapshot
            script_lines.append(
                f"- {model.input_spec.name}: exp_id={snapshot.exp_id}, read_pool_id={snapshot.read_pool_id}, "
                f"enc_run_id={snapshot.enc_run_id}, refs={len(model.reference_df)}, "
                f"zero_refs={(model.reference_df['n_reads_on_reference'] == 0).sum()}, "
                f"global_mean={model.global_mean_coverage:.6g}, refs_with_delta_g={model.reference_df['delta_g'].notna().sum()}"
            )
            script_lines.append(
                f"  snapshot: reads_total={snapshot.n_reads_total}, mapped={snapshot.n_reads_mapped}, "
                f"unmapped={snapshot.n_reads_unmapped}, fastq={snapshot.n_fastq_total}, "
                f"source_max_position_exclusive={snapshot.source_max_position_exclusive}, "
                f"version={snapshot.snapshot_version}, refreshed_at={snapshot.refreshed_at}"
            )

        sections: list[MetadataSection] = []
        for index, plot_spec in enumerate(self.context.script.plot_settings, start=1):
            settings = self._settings_for_plot(plot_spec.name)
            lines = [f"Figure settings: {settings.model_dump(exclude={'title', 'description'})}"]
            if plot_spec.name == "ref-coverage-in-read-pool_distribution":
                lines.extend(
                    f"{model.input_spec.name}: mean read-pool coverage used as normalization denominator="
                    f"{model.global_mean_coverage:.6g} reads/reference"
                    for model in self.input_models
                )
            elif plot_spec.name == "ref-coverage-in-read-pool_vs_delta-g_scatter":
                lines.extend(self._scatter_metadata_lines())
            sections.append(
                MetadataSection(
                    title=f"{index}. {settings.title or DEFAULT_PLOT_TITLES[plot_spec.name]}",
                    description=settings.description or DEFAULT_PLOT_EXPLANATIONS[plot_spec.name],
                    lines=tuple(lines),
                )
            )
        return build_metadata_pages(self.script_name, script_lines, sections)

    def _scatter_metadata_lines(self) -> list[str]:
        lines: list[str] = []
        for model in self.input_models:
            correlation = compute_spearman_statistics(model.reference_df)
            regression = compute_regression_statistics(model.reference_df)
            correlation_text = (
                "rho=N/A, p_value=N/A"
                if correlation.rho is None or correlation.p_value is None
                else f"rho={correlation.rho:.6g}, p_value={correlation.p_value:.6g}"
            )
            regression_text = (
                "slope=N/A, intercept=N/A, R2=N/A"
                if regression.slope is None or regression.intercept is None or regression.r_squared is None
                else (
                    f"slope={regression.slope:.6g}, intercept={regression.intercept:.6g}, "
                    f"R2={regression.r_squared:.6g}"
                )
            )
            lines.append(
                f"{model.input_spec.name}: n={correlation.n_points}, Spearman {correlation_text}; "
                f"linear regression {regression_text}"
            )
        return lines

    @staticmethod
    def _item_name_map(input_spec: InputSpec) -> dict[int, str]:
        return {item.item_id: item.name for item in input_spec.items}