"""Figures de metriques au fil des passes et au premier decodage."""

from __future__ import annotations

from dataclasses import dataclass
import sqlite3
from typing import ClassVar, Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from pydantic import BaseModel, ConfigDict, field_validator

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


PROGRESSION_PLOTS: dict[str, str] = {
    "psnr_vs_run-progression": "psnr",
    "hamming-dist_vs_run-progression": "hamming_distance_normalized",
    "pdpr_vs_run-progression": "perfectly_decoded_payload_ratio",
    "fpdpc_vs_run-progression": "fpdpc",
}
AT_DECODING_PLOTS: dict[str, str] = {
    "psnr-at-decoding_vs_run-number": "psnr_at_decoding",
    "hamming-dist-at-decoding_vs_run-number": "hamming_distance_normalized_at_decoding",
    "pdpr-at-decoding_vs_run-number": "perfectly_decoded_payload_ratio_at_decoding",
    "fpdpc-at-decoding_vs_run-number": "fpdpc_at_decoding",
    "coverage-at-decoding_vs_run-number": "coverage_at_decoding",
    "estimated-seq-duration-at-decoding_vs_run-number": "estimated_sequencing_duration_at_decoding",
}
SUPPORTED_DECODING_RUN_METRICS_PLOTS = set(PROGRESSION_PLOTS) | set(AT_DECODING_PLOTS)

METRIC_LABELS: dict[str, str] = {
    "psnr": "PSNR (dB)",
    "hamming_distance_normalized": "Normalized Hamming distance",
    "perfectly_decoded_payload_ratio": "Perfectly decoded payload ratio (%)",
    "fpdpc": "First-time perfectly decoded payload count",
    "psnr_at_decoding": "PSNR at decoding (dB)",
    "hamming_distance_normalized_at_decoding": "Normalized Hamming distance at decoding",
    "perfectly_decoded_payload_ratio_at_decoding": "Perfectly decoded payload ratio at decoding (%)",
    "fpdpc_at_decoding": "First-time perfectly decoded payload count at decoding",
    "coverage_at_decoding": "Coverage at decoding",
    "estimated_sequencing_duration_at_decoding": "Estimated sequencing duration at decoding (minutes)",
}

AT_DECODING_UNIT_FACTORS: dict[str, float] = {
    "perfectly_decoded_payload_ratio_at_decoding": 100.0,
    "estimated_sequencing_duration_at_decoding": 1.0 / 60.0,
}


def convert_at_decoding_to_display_units(at_decoding_df: pd.DataFrame) -> pd.DataFrame:
    converted = at_decoding_df.copy()
    for column, factor in AT_DECODING_UNIT_FACTORS.items():
        if column in converted:
            converted[column] *= factor
    return converted

X_AXIS_COLUMNS: dict[str, str] = {
    "pass_index": "dec_pass_id",
    "coverage": "coverage",
    "n_tot_reads": "N_trimmed_reads_cum_all_items",
    "estimated_sequencing_duration": "estimated_sequencing_duration",
    "run_duration": "run_duration",
}
X_AXIS_LABELS: dict[str, str] = {
    "pass_index": "Pass index",
    "coverage": "Coverage",
    "n_tot_reads": "Total trimmed reads (cumulative across items)",
    "estimated_sequencing_duration": "Estimated sequencing duration",
    "run_duration": "Run duration",
}
ITEM_ID_COLOR_MAP = {
    0: "#1f77b4",
    1: "#d62728",
    2: "#2ca02c",
    3: "#ff7f0e",
    4: "#9467bd",
    5: "#8c564b",
}

DEFAULT_PLOT_TITLES = {
    plot_name: (
        f"{METRIC_LABELS[metric_column]} over run progression"
        if plot_name in PROGRESSION_PLOTS
        else f"{METRIC_LABELS[metric_column]} versus run number"
    )
    for plot_name, metric_column in {**PROGRESSION_PLOTS, **AT_DECODING_PLOTS}.items()
}
DEFAULT_PLOT_EXPLANATIONS = {
    **{
        plot_name: (
            "Each item curve is the pointwise average of decoding runs interpolated on a shared X grid. "
            "The shaded interval is the standard error when at least two runs contribute."
        )
        for plot_name in PROGRESSION_PLOTS
    },
    **{
        plot_name: (
            "Each point is the metric value at the first successful item-decoding transition for one run. "
            "Failed decoding runs are excluded and gaps in database run numbers are preserved."
        )
        for plot_name in AT_DECODING_PLOTS
    },
}


class DecodingRunMetricsSettings(BaseModel):
    model_config = ConfigDict(extra="forbid")

    x_axis_key: Literal[
        "pass_index",
        "coverage",
        "n_tot_reads",
        "estimated_sequencing_duration",
        "run_duration",
    ] = "coverage"
    common_x_origin: float = 0.0
    grid_n_points: int = 220
    min_points_per_curve: int = 2

    @field_validator("grid_n_points", "min_points_per_curve")
    @classmethod
    def _validate_minimum_two(cls, value: int) -> int:
        if value < 2:
            raise ValueError("grid_n_points and min_points_per_curve must be >= 2.")
        return value


class RunProgressionPlotSettings(PlotTextSettings):
    plot_individual_curves: bool = False
    plot_average_curve: bool = True
    plot_se_band: bool = True
    plot_std_band: bool = False


class AtDecodingPlotSettings(PlotTextSettings):
    display_metric_mean: bool = False


SQL_QUERY_PROGRESSION = """
SELECT
    m.exp_id,
    m.dec_run_id,
    d.dec_run_name,
    d.dec_run_number,
    m.item_id,
    m.dec_pass_id,
    m.coverage,
    m.hamming_distance_normalized,
    m.perfectly_decoded_payload_ratio,
    m.fpdpc,
    m.psnr,
    p.N_trimmed_reads_cum_all_items,
    p.N_decoding_reads_cum_all_items,
    p.run_duration,
    p.estimated_sequencing_duration
FROM metrics_view_run_item_pass m
JOIN metrics_view_run_pass p
  ON p.exp_id = m.exp_id
 AND p.dec_run_id = m.dec_run_id
 AND p.dec_pass_id = m.dec_pass_id
JOIN decoding_run_label_record l
  ON l.exp_id = m.exp_id
 AND l.dec_run_id = m.dec_run_id
JOIN decoding_run d
  ON d.exp_id = m.exp_id
 AND d.dec_run_id = m.dec_run_id
WHERE l.label = ?
  AND m.item_id IN ({item_placeholders})
ORDER BY d.dec_run_name, d.dec_run_number, m.item_id, m.dec_pass_id
"""

SQL_QUERY_AT_DECODING = """
SELECT
    m.exp_id,
    m.dec_run_id,
    d.dec_run_name,
    d.dec_run_number,
    m.item_id,
    m.pass_at_decoding,
    m.psnr_at_decoding,
    m.hamming_distance_normalized_at_decoding,
    m.perfectly_decoded_payload_ratio_at_decoding,
    m.fpdpc_at_decoding,
    m.coverage_at_decoding,
    m.estimated_sequencing_duration_at_decoding
FROM metrics_view_run_item_at_decoding m
JOIN decoding_run_label_record l
  ON l.exp_id = m.exp_id
 AND l.dec_run_id = m.dec_run_id
JOIN decoding_run d
  ON d.exp_id = m.exp_id
 AND d.dec_run_id = m.dec_run_id
WHERE l.label = ?
  AND m.item_id IN ({item_placeholders})
ORDER BY d.dec_run_name, d.dec_run_number, m.item_id
"""

SQL_QUERY_LABEL_RUNS_TOTAL = """
SELECT COUNT(*)
FROM (
    SELECT DISTINCT exp_id, dec_run_id
    FROM decoding_run_label_record
    WHERE label = ?
)
"""


@dataclass
class InputReadModel:
    input_spec: InputSpec
    progression_df: pd.DataFrame
    at_decoding_df: pd.DataFrame
    n_runs_total_labeled: int
    run_name: str


class DecodingRunMetricsScript:
    script_name: ClassVar[str] = "decoding_run-metrics"

    def __init__(self, context: ScriptExecutionContext):
        self.context = context
        self._validate_requested_plots()
        self.settings = DecodingRunMetricsSettings.model_validate(context.script.script_settings)
        self.plot_specs_by_name = {plot.name: plot for plot in context.script.plot_settings}
        self.input_models: list[InputReadModel] = []

    def _validate_requested_plots(self) -> None:
        unknown_plots = sorted(
            plot.name
            for plot in self.context.script.plot_settings
            if plot.name not in SUPPORTED_DECODING_RUN_METRICS_PLOTS
        )
        if unknown_plots:
            raise ValueError(
                f"Unknown plot(s) for script '{self.script_name}': {', '.join(unknown_plots)}. "
                f"Allowed values: {', '.join(sorted(SUPPORTED_DECODING_RUN_METRICS_PLOTS))}"
            )

    def _settings_for_plot(self, plot_name: str) -> PlotTextSettings:
        model_type: type[PlotTextSettings]
        if plot_name in PROGRESSION_PLOTS:
            model_type = RunProgressionPlotSettings
        else:
            model_type = AtDecodingPlotSettings
        return model_type.model_validate(self.plot_specs_by_name[plot_name].settings)

    def run(self) -> None:
        needs_progression = any(plot.name in PROGRESSION_PLOTS for plot in self.context.script.plot_settings)
        needs_at_decoding = any(plot.name in AT_DECODING_PLOTS for plot in self.context.script.plot_settings)
        self.input_models = [
            self._build_read_model_for_input(input_spec, needs_progression, needs_at_decoding)
            for input_spec in self.context.inputs
        ]

        output_dir = (self.context.output_path / self.script_name).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        output_pdf = output_dir / f"{self.script_name}.pdf"
        plot_dispatcher = {
            **{name: self._plot_run_progression for name in PROGRESSION_PLOTS},
            **{name: self._plot_at_decoding_vs_run_number for name in AT_DECODING_PLOTS},
        }

        with PdfPages(output_pdf) as pdf:
            for metadata_fig in self._build_metadata_pages():
                pdf.savefig(metadata_fig)
                plt.close(metadata_fig)
            for plot_spec in self.context.script.plot_settings:
                print(f"[{self.script_name}] rendering plot={plot_spec.name}")
                fig = plot_dispatcher[plot_spec.name](plot_spec.name)
                save_figure_page_and_png(
                    fig,
                    plot_spec.name,
                    output_dir,
                    pdf,
                    show_figure=self.context.show_figures,
                )

    def _build_read_model_for_input(
        self,
        input_spec: InputSpec,
        needs_progression: bool,
        needs_at_decoding: bool,
    ) -> InputReadModel:
        item_ids = [item.item_id for item in input_spec.items]
        if not item_ids:
            raise ValueError(f"Input '{input_spec.name}' has no configured items.")
        database_path = resolve_database_path(input_spec.database)
        if not database_path.exists():
            raise FileNotFoundError(f"Database not found for input '{input_spec.name}': {database_path}")

        placeholders = ",".join(["?"] * len(item_ids))
        query_params = [input_spec.decoding_run_label, *item_ids]
        with sqlite3.connect(database_path) as connection:
            progression_df = (
                pd.read_sql_query(
                    SQL_QUERY_PROGRESSION.format(item_placeholders=placeholders),
                    connection,
                    params=query_params,
                )
                if needs_progression
                else pd.DataFrame()
            )
            at_decoding_df = (
                pd.read_sql_query(
                    SQL_QUERY_AT_DECODING.format(item_placeholders=placeholders),
                    connection,
                    params=query_params,
                )
                if needs_at_decoding
                else pd.DataFrame()
            )
            n_runs_total_labeled = int(
                connection.execute(SQL_QUERY_LABEL_RUNS_TOTAL, (input_spec.decoding_run_label,)).fetchone()[0]
            )

        frames = [frame for frame in (progression_df, at_decoding_df) if not frame.empty]
        if not frames:
            raise ValueError(
                f"No rows returned for input '{input_spec.name}'. Check decoding_run_label and configured items."
            )
        run_names = sorted({str(value) for frame in frames for value in frame["dec_run_name"].dropna().unique()})
        if len(run_names) != 1:
            raise ValueError(
                f"Input '{input_spec.name}' label must select exactly one dec_run_name to plot dec_run_number; "
                f"found {run_names or 'none'}."
            )

        progression_numeric = [
            "dec_run_number",
            "dec_pass_id",
            "coverage",
            "hamming_distance_normalized",
            "perfectly_decoded_payload_ratio",
            "fpdpc",
            "psnr",
            "N_trimmed_reads_cum_all_items",
            "N_decoding_reads_cum_all_items",
            "run_duration",
            "estimated_sequencing_duration",
        ]
        for column in progression_numeric:
            if column in progression_df:
                progression_df[column] = pd.to_numeric(progression_df[column], errors="coerce")
        at_decoding_numeric = [
            "dec_run_number",
            "pass_at_decoding",
            *AT_DECODING_PLOTS.values(),
        ]
        for column in at_decoding_numeric:
            if column in at_decoding_df:
                at_decoding_df[column] = pd.to_numeric(at_decoding_df[column], errors="coerce")

        if "perfectly_decoded_payload_ratio" in progression_df:
            progression_df["perfectly_decoded_payload_ratio"] *= 100.0
        at_decoding_df = convert_at_decoding_to_display_units(at_decoding_df)

        return InputReadModel(
            input_spec=input_spec,
            progression_df=progression_df,
            at_decoding_df=at_decoding_df,
            n_runs_total_labeled=n_runs_total_labeled,
            run_name=run_names[0],
        )

    def _plot_run_progression(self, plot_name: str) -> Figure:
        metric_column = PROGRESSION_PLOTS[plot_name]
        plot_settings = RunProgressionPlotSettings.model_validate(
            self.plot_specs_by_name[plot_name].settings
        )
        x_column = X_AXIS_COLUMNS[self.settings.x_axis_key]
        fig, axes = make_axes_grid(
            len(self.input_models),
            height=self.context.figure_height,
            width_per_input=self.context.figure_width_per_input,
        )

        for axis, input_model in zip(axes, self.input_models):
            item_names = self._item_name_map(input_model.input_spec)
            for item_id, item_name in item_names.items():
                item_df = input_model.progression_df[input_model.progression_df["item_id"] == item_id]
                curves: list[tuple[np.ndarray, np.ndarray]] = []
                for _, run_df in item_df.groupby(["exp_id", "dec_run_id"], sort=False):
                    curve = self._prepare_xy_curve(run_df, x_column, metric_column)
                    if curve is not None:
                        curves.append(curve)
                if not curves:
                    continue

                grid_start = max(self.settings.common_x_origin, min(float(x.min()) for x, _ in curves))
                grid_end = float(np.mean([x.max() for x, _ in curves]))
                if grid_end <= grid_start:
                    continue
                x_grid = np.linspace(grid_start, grid_end, self.settings.grid_n_points)
                y_matrix = np.vstack([self._interpolate_on_grid(x, y, x_grid) for x, y in curves])
                n_contributors = np.sum(~np.isnan(y_matrix), axis=0)
                sums = np.nansum(y_matrix, axis=0)
                y_mean = np.divide(
                    sums,
                    n_contributors,
                    out=np.full_like(sums, np.nan, dtype=float),
                    where=n_contributors > 0,
                )
                y_std = np.full_like(y_mean, np.nan)
                valid_spread = n_contributors > 1
                if np.any(valid_spread):
                    y_std[valid_spread] = np.nanstd(y_matrix[:, valid_spread], axis=0, ddof=1)
                y_se = np.divide(
                    y_std,
                    np.sqrt(n_contributors),
                    out=np.full_like(y_std, np.nan),
                    where=valid_spread,
                )
                color = ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f")

                if plot_settings.plot_individual_curves:
                    for y_values in y_matrix:
                        axis.plot(x_grid, y_values, color=color, linewidth=1.0, alpha=0.20)
                if plot_settings.plot_std_band:
                    axis.fill_between(x_grid, y_mean - y_std, y_mean + y_std, color=color, alpha=0.10)
                if plot_settings.plot_se_band:
                    axis.fill_between(x_grid, y_mean - y_se, y_mean + y_se, color=color, alpha=0.18)
                if plot_settings.plot_average_curve:
                    axis.plot(
                        x_grid,
                        y_mean,
                        color=color,
                        linewidth=2.4,
                        label=f"{item_name} AVG (n_runs={len(curves)})",
                    )

            axis.set_title(input_model.input_spec.name)
            axis.set_xlabel(X_AXIS_LABELS[self.settings.x_axis_key])
            axis.set_ylabel(METRIC_LABELS[metric_column])
            axis.grid(True, alpha=0.25)
            if axis.has_data():
                axis.legend(fontsize=8)
            else:
                axis.text(0.5, 0.5, "No valid curves", ha="center", va="center", transform=axis.transAxes)

        harmonize_axes_scales(
            axes,
            same_x_scale=plot_settings.same_x_scale_across_inputs,
            same_y_scale=plot_settings.same_y_scale_across_inputs,
        )
        apply_figure_title(fig, plot_settings.title or DEFAULT_PLOT_TITLES[plot_name])
        return fig

    def _plot_at_decoding_vs_run_number(self, plot_name: str) -> Figure:
        metric_column = AT_DECODING_PLOTS[plot_name]
        plot_settings = AtDecodingPlotSettings.model_validate(self.plot_specs_by_name[plot_name].settings)
        fig, axes = make_axes_grid(
            len(self.input_models),
            height=self.context.figure_height,
            width_per_input=self.context.figure_width_per_input,
        )

        for axis, input_model in zip(axes, self.input_models):
            item_names = self._item_name_map(input_model.input_spec)
            decoded_df = input_model.at_decoding_df.dropna(subset=["pass_at_decoding"])
            for item_id, item_name in item_names.items():
                item_df = decoded_df[decoded_df["item_id"] == item_id].dropna(
                    subset=["dec_run_number", metric_column]
                ).sort_values("dec_run_number")
                if item_df.empty:
                    continue
                color = ITEM_ID_COLOR_MAP.get(item_id, "#7f7f7f")
                axis.plot(
                    item_df["dec_run_number"],
                    item_df[metric_column],
                    color=color,
                    marker="o",
                    markersize=3.5,
                    linewidth=1.5,
                    label=item_name,
                )
                if plot_settings.display_metric_mean:
                    mean_value = float(item_df[metric_column].mean())
                    axis.axhline(
                        mean_value,
                        color=color,
                        linestyle="--",
                        linewidth=1.2,
                        alpha=0.8,
                        label=f"{item_name} mean={mean_value:.3g}",
                    )

            axis.set_title(input_model.input_spec.name)
            axis.set_xlabel("Decoding run number")
            axis.set_ylabel(METRIC_LABELS[metric_column])
            axis.grid(True, alpha=0.25)
            if axis.has_data():
                axis.legend(fontsize=8)
            else:
                axis.text(0.5, 0.5, "No decoded values", ha="center", va="center", transform=axis.transAxes)

        harmonize_axes_scales(
            axes,
            same_x_scale=plot_settings.same_x_scale_across_inputs,
            same_y_scale=plot_settings.same_y_scale_across_inputs,
        )
        apply_figure_title(fig, plot_settings.title or DEFAULT_PLOT_TITLES[plot_name])
        return fig

    def _build_metadata_pages(self) -> list[Figure]:
        script_lines = [
            f"Output directory: {self.context.output_path / self.script_name}",
            f"Configured figures: {len(self.context.script.plot_settings)}",
            f"Script settings: {self.settings.model_dump()}",
            *(
                f"- {model.input_spec.name}: label={model.input_spec.decoding_run_label}, "
                f"dec_run_name={model.run_name}, runs_labelled={model.n_runs_total_labeled}"
                for model in self.input_models
            ),
        ]
        sections: list[MetadataSection] = []
        for index, plot_spec in enumerate(self.context.script.plot_settings, start=1):
            plot_name = plot_spec.name
            settings = self._settings_for_plot(plot_name)
            lines = [f"Figure settings: {settings.model_dump(exclude={'title', 'description'})}"]
            if plot_name in PROGRESSION_PLOTS:
                for model in self.input_models:
                    metric_column = PROGRESSION_PLOTS[plot_name]
                    item_lines = []
                    for item_id, item_name in self._item_name_map(model.input_spec).items():
                        item_df = model.progression_df[model.progression_df["item_id"] == item_id]
                        valid_runs = sum(
                            self._prepare_xy_curve(run_df, X_AXIS_COLUMNS[self.settings.x_axis_key], metric_column)
                            is not None
                            for _, run_df in item_df.groupby(["exp_id", "dec_run_id"], sort=False)
                        )
                        item_lines.append(f"{item_name}: valid_curves={valid_runs}")
                    lines.append(f"{model.input_spec.name}: " + "; ".join(item_lines))
            else:
                metric_column = AT_DECODING_PLOTS[plot_name]
                for model in self.input_models:
                    decoded_df = model.at_decoding_df.dropna(subset=["pass_at_decoding"])
                    lines.append(f"{model.input_spec.name}:")
                    for item_id, item_name in self._item_name_map(model.input_spec).items():
                        values = decoded_df.loc[decoded_df["item_id"] == item_id, metric_column].dropna()
                        decoded_count = int(
                            decoded_df.loc[decoded_df["item_id"] == item_id, "dec_run_id"].nunique()
                        )
                        lines.append(
                            "  "
                            + self._format_statistics(
                                item_name,
                                values,
                                decoded_count,
                                model.n_runs_total_labeled,
                            )
                        )

            sections.append(
                MetadataSection(
                    title=f"{index}. {settings.title or DEFAULT_PLOT_TITLES[plot_name]}",
                    description=settings.description or DEFAULT_PLOT_EXPLANATIONS[plot_name],
                    lines=tuple(lines),
                )
            )
        return build_metadata_pages(self.script_name, script_lines, sections)

    def _prepare_xy_curve(
        self,
        run_item_df: pd.DataFrame,
        x_column: str,
        metric_column: str,
    ) -> tuple[np.ndarray, np.ndarray] | None:
        curve = run_item_df[[x_column, metric_column]].dropna().copy()
        if curve.empty:
            return None
        curve = curve.groupby(x_column, as_index=False).agg({metric_column: "mean"})
        curve = curve.sort_values(by=x_column)
        if len(curve) < self.settings.min_points_per_curve:
            return None
        return (
            curve[x_column].to_numpy(dtype=float),
            curve[metric_column].to_numpy(dtype=float),
        )

    @staticmethod
    def _interpolate_on_grid(
        x_values: np.ndarray,
        y_values: np.ndarray,
        x_grid: np.ndarray,
    ) -> np.ndarray:
        result = np.full_like(x_grid, np.nan, dtype=float)
        in_support = (x_grid >= x_values.min()) & (x_grid <= x_values.max())
        result[in_support] = np.interp(x_grid[in_support], x_values, y_values)
        return result

    @staticmethod
    def _format_statistics(
        item_name: str,
        values: pd.Series,
        decoded_count: int,
        labelled_count: int,
    ) -> str:
        numeric = pd.to_numeric(values, errors="coerce").dropna()
        if numeric.empty:
            return f"{item_name}: n=0, decoded/labelled={decoded_count}/{labelled_count}, stats=N/A"
        std_value = float(numeric.std(ddof=1)) if len(numeric) >= 2 else float("nan")
        std_text = f"{std_value:.4g}" if np.isfinite(std_value) else "N/A"
        return (
            f"{item_name}: n={len(numeric)}, decoded/labelled={decoded_count}/{labelled_count}, "
            f"mean={numeric.mean():.4g}, std={std_text}, min={numeric.min():.4g}, max={numeric.max():.4g}"
        )

    @staticmethod
    def _item_name_map(input_spec: InputSpec) -> dict[int, str]:
        return {item.item_id: item.name for item in input_spec.items}
