from __future__ import annotations

from pathlib import Path
import sqlite3
import sys

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import pandas as pd
import pytest
from pydantic import ValidationError

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from optional_script_utils import save_figure_page_and_png  # noqa: E402
from protocols import InputSpec, ItemSpec, PlotSpec, ScriptExecutionContext, ScriptSpec  # noqa: E402
from read_pool_stats import (  # noqa: E402
    DistributionPlotSettings,
    InputReadModel,
    MeanInBinPlotSettings,
    PoolSnapshot,
    ReadPoolStatsScript,
    ReadPoolStatsSettings,
    compute_distribution_bin_edges,
    compute_global_delta_g_bins,
    compute_pooled_mean_coverage,
    compute_regression_statistics,
    compute_shared_normalized_bin_edges,
    compute_spearman_statistics,
    normalize_reference_coverage,
)


def _reference_rows() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "item_id": [0, 0, 0, 1, 1, 1],
            "delta_g": [1.0, 2.0, 3.0, 1.0, 2.0, None],
            "n_reads_on_reference": [0, 2, 4, 0, 6, 3],
        }
    )


def _input_spec(database: Path) -> InputSpec:
    return InputSpec(
        input_id=1,
        name="Synthetic input",
        database=database,
        exp_id="test-exp",
        read_pool_id="test-pool",
        decoding_run_label="unused",
        items=[
            ItemSpec(item_id=0, name="Item zero"),
            ItemSpec(item_id=1, name="Item one"),
        ],
    )


def _script_context(
    database: Path,
    output_path: Path,
    *,
    figure_width_per_input: float = 6.0,
    figure_height: float = 5.0,
) -> ScriptExecutionContext:
    script_spec = ScriptSpec(
        name="read_pool_stats",
        plot_settings=[
            PlotSpec(name="ref-coverage-in-read-pool_distribution"),
            PlotSpec(
                name="ref-coverage-in-read-pool_vs_delta-g_scatter",
                settings={"display_linear_regression": True},
            ),
            PlotSpec(
                name="ref-coverage-in-read-pool_vs_delta-g_mean-in-bin",
                settings={"delta_g_precision": 1.0, "min_points_per_bin": 1},
            ),
        ],
        script_settings={},
    )
    return ScriptExecutionContext(
        output_path=output_path,
        figure_width_per_input=figure_width_per_input,
        figure_height=figure_height,
        inputs=[_input_spec(database)],
        script=script_spec,
    )


def _create_snapshot_database(database: Path, *, include_snapshot: bool = True) -> None:
    with sqlite3.connect(database) as connection:
        connection.executescript(
            """
            CREATE TABLE read_pool_stats_snapshot (
                exp_id TEXT NOT NULL,
                read_pool_id TEXT NOT NULL,
                enc_run_id TEXT,
                n_reads_total INTEGER NOT NULL,
                n_reads_mapped INTEGER NOT NULL,
                n_reads_unmapped INTEGER NOT NULL,
                n_fastq_total INTEGER NOT NULL,
                source_max_position_exclusive INTEGER NOT NULL,
                refreshed_at TEXT NOT NULL,
                snapshot_version INTEGER NOT NULL,
                PRIMARY KEY (exp_id, read_pool_id)
            );

            CREATE TABLE pool_reference_stats_source (
                exp_id TEXT NOT NULL,
                read_pool_id TEXT NOT NULL,
                enc_run_id TEXT NOT NULL,
                item_id INTEGER NOT NULL,
                region_id INTEGER NOT NULL,
                position_id INTEGER NOT NULL,
                delta_g REAL,
                n_reads_on_reference INTEGER NOT NULL
            );

            CREATE VIEW metrics_view_pool_reference_stats AS
            SELECT * FROM pool_reference_stats_source;
            """
        )
        if include_snapshot:
            connection.execute(
                """
                INSERT INTO read_pool_stats_snapshot VALUES
                    ('test-exp', 'test-pool', 'enc-run', 18, 15, 3, 2, 18,
                     '2026-08-11 10:00:00', 1)
                """
            )
            rows = _reference_rows()
            connection.executemany(
                """
                INSERT INTO pool_reference_stats_source VALUES
                    ('test-exp', 'test-pool', 'enc-run', ?, 0, ?, ?, ?)
                """,
                [
                    (int(row.item_id), index, row.delta_g, int(row.n_reads_on_reference))
                    for index, row in enumerate(rows.itertuples(index=False))
                ],
            )


def test_settings_are_strict_and_validate_binning() -> None:
    with pytest.raises(ValidationError):
        ReadPoolStatsSettings.model_validate({"analysis_mode": "global"})
    with pytest.raises(ValidationError):
        MeanInBinPlotSettings.model_validate({"delta_g_precision": 0})
    with pytest.raises(ValidationError):
        MeanInBinPlotSettings.model_validate({"min_points_per_bin": 0})
    with pytest.raises(ValidationError):
        DistributionPlotSettings.model_validate({"read_count_values_per_bin": 0})
    with pytest.raises(ValidationError, match="requires read_count_values_per_bin"):
        DistributionPlotSettings.model_validate(
            {"bin_width_mode": "shared_normalized", "read_count_values_per_bin": None}
        )
    with pytest.raises(ValidationError):
        _script_context(Path("unused.db"), Path("plots"), figure_width_per_input=0)


def test_normalization_is_global_across_items_and_includes_zeros() -> None:
    normalized, global_mean = normalize_reference_coverage(_reference_rows())

    assert global_mean == pytest.approx(2.5)
    assert normalized["normalized_read_pool_coverage"].tolist() == pytest.approx(
        [0.0, 0.8, 1.6, 0.0, 2.4, 1.2]
    )
    assert normalized.groupby("item_id")["normalized_read_pool_coverage"].mean().tolist() == pytest.approx(
        [0.8, 1.2]
    )


def test_normalization_fails_when_global_mean_is_zero() -> None:
    rows = pd.DataFrame({"item_id": [0, 1], "n_reads_on_reference": [0, 0]})

    with pytest.raises(ValueError, match="Global mean reference coverage is zero"):
        normalize_reference_coverage(rows)


def test_distribution_can_group_multiple_integer_read_counts_per_bin() -> None:
    normalized, global_mean = normalize_reference_coverage(_reference_rows())

    bin_edges = compute_distribution_bin_edges(
        normalized,
        global_mean_coverage=global_mean,
        read_count_values_per_bin=3,
    )

    assert bin_edges.tolist() == pytest.approx([0.0, 1.2, 2.4, 3.6])


def test_shared_bins_use_reference_weighted_pooled_mean_and_strict_upper_edge() -> None:
    first, first_mean = normalize_reference_coverage(
        pd.DataFrame({"n_reads_on_reference": [0, 10]})
    )
    second, second_mean = normalize_reference_coverage(
        pd.DataFrame({"n_reads_on_reference": [100]})
    )

    pooled_mean = compute_pooled_mean_coverage([first, second])
    bin_edges = compute_shared_normalized_bin_edges(
        [first, second],
        pooled_mean_coverage=pooled_mean,
        read_count_values_per_bin=55,
    )

    assert first_mean == pytest.approx(5.0)
    assert second_mean == pytest.approx(100.0)
    assert pooled_mean == pytest.approx((2 * first_mean + second_mean) / 3)
    assert bin_edges.tolist() == pytest.approx([0.0, 1.5, 3.0])
    assert bin_edges[-1] > max(first["normalized_read_pool_coverage"].max(), 1.0)
    assert first_mean * (bin_edges[1] - bin_edges[0]) == pytest.approx(7.5)
    assert second_mean * (bin_edges[1] - bin_edges[0]) == pytest.approx(150.0)


def test_shared_bins_add_an_upper_bin_when_maximum_is_exactly_on_an_edge() -> None:
    normalized = pd.DataFrame(
        {
            "n_reads_on_reference": [0, 10, 20],
            "normalized_read_pool_coverage": [0.0, 1.0, 2.0],
        }
    )

    bin_edges = compute_shared_normalized_bin_edges(
        [normalized],
        pooled_mean_coverage=10.0,
        read_count_values_per_bin=10,
    )

    assert bin_edges.tolist() == pytest.approx([0.0, 1.0, 2.0, 3.0])


def test_scatter_statistics_are_global_and_ignore_missing_delta_g() -> None:
    normalized, _ = normalize_reference_coverage(_reference_rows())

    correlation = compute_spearman_statistics(normalized)
    regression = compute_regression_statistics(normalized)

    assert correlation.n_points == 5
    assert correlation.rho is not None
    assert correlation.p_value is not None
    assert regression.n_points == 5
    assert regression.slope is not None
    assert regression.intercept is not None
    assert regression.r_squared is not None


def test_delta_g_bins_are_item_blind() -> None:
    normalized, _ = normalize_reference_coverage(_reference_rows())

    binned = compute_global_delta_g_bins(normalized, delta_g_precision=1.0, min_points_per_bin=1)

    assert binned["n_points"].tolist() == [2, 2, 1]
    assert binned["mean_normalized_coverage"].tolist() == pytest.approx([0.0, 1.6, 1.6])


def test_missing_snapshot_fails_with_update_guidance(tmp_path: Path) -> None:
    database = tmp_path / "missing-snapshot.db"
    _create_snapshot_database(database, include_snapshot=False)
    script = ReadPoolStatsScript(_script_context(database, tmp_path / "plots"))

    with pytest.raises(ValueError, match="refresh this read-pool statistics snapshot"):
        script._build_read_model_for_input(_input_spec(database))


def test_missing_view_fails_with_schema_update_guidance(tmp_path: Path) -> None:
    database = tmp_path / "missing-view.db"
    _create_snapshot_database(database)
    with sqlite3.connect(database) as connection:
        connection.execute("DROP VIEW metrics_view_pool_reference_stats")
    script = ReadPoolStatsScript(_script_context(database, tmp_path / "plots"))

    with pytest.raises(ValueError, match="Update the database schema"):
        script._build_read_model_for_input(_input_spec(database))


def test_read_model_excludes_items_not_selected_in_configuration(tmp_path: Path) -> None:
    database = tmp_path / "configured-items.db"
    _create_snapshot_database(database)
    with sqlite3.connect(database) as connection:
        connection.execute(
            """
            INSERT INTO pool_reference_stats_source VALUES
                ('test-exp', 'test-pool', 'enc-run', 99, 0, 99, NULL, 50)
            """
        )
    script = ReadPoolStatsScript(_script_context(database, tmp_path / "plots"))

    input_model = script._build_read_model_for_input(_input_spec(database))

    assert set(input_model.reference_df["item_id"]) == {0, 1}


def test_read_pool_figure_size_is_configurable_for_one_input(tmp_path: Path) -> None:
    database = tmp_path / "sized-figure.db"
    _create_snapshot_database(database)
    context = _script_context(
        database,
        tmp_path / "plots",
        figure_width_per_input=9.0,
        figure_height=5.5,
    )
    script = ReadPoolStatsScript(context)
    script.input_models = [script._build_read_model_for_input(_input_spec(database))]

    figure = script._plot_distribution()

    try:
        assert figure.get_size_inches().tolist() == pytest.approx([9.0, 5.5])
    finally:
        plt.close(figure)


def test_distribution_normalizes_each_item_population_on_y_axis(tmp_path: Path) -> None:
    database = tmp_path / "normalized-y.db"
    _create_snapshot_database(database)
    script = ReadPoolStatsScript(_script_context(database, tmp_path / "plots"))
    script.input_models = [script._build_read_model_for_input(_input_spec(database))]

    figure = script._plot_distribution()

    try:
        axis = figure.axes[0]
        assert sum(bar.get_height() for bar in axis.patches) == pytest.approx(2.0)
        assert axis.get_ylabel() == "Fraction of item references"
    finally:
        plt.close(figure)


def test_shared_normalized_mode_uses_common_bar_widths_and_reports_effective_widths() -> None:
    item_ids = [0, 0, 1, 1]
    input_specs = [_input_spec(Path(f"unused-{index}.db")) for index in range(2)]
    input_specs[0].name = "Low coverage"
    input_specs[1].name = "High coverage"
    reference_dfs = []
    means = []
    for raw_counts in ([0, 10, 20, 30], [0, 100, 200, 300]):
        normalized_df, mean = normalize_reference_coverage(
            pd.DataFrame(
                {
                    "item_id": item_ids,
                    "delta_g": [None] * 4,
                    "n_reads_on_reference": raw_counts,
                }
            )
        )
        reference_dfs.append(normalized_df)
        means.append(mean)
    snapshot = PoolSnapshot(
        exp_id="test-exp",
        read_pool_id="test-pool",
        enc_run_id="enc-run",
        n_reads_total=1,
        n_reads_mapped=1,
        n_reads_unmapped=0,
        n_fastq_total=1,
        source_max_position_exclusive=1,
        refreshed_at="2026-08-12",
        snapshot_version=1,
    )
    context = ScriptExecutionContext(
        output_path=Path("plots"),
        inputs=input_specs,
        script=ScriptSpec(
            name="read_pool_stats",
            plot_settings=[
                PlotSpec(
                    name="ref-coverage-in-read-pool_distribution",
                    settings={
                        "read_count_values_per_bin": 50,
                        "bin_width_mode": "shared_normalized",
                        "bin_width_reference": "pooled",
                    },
                )
            ],
        ),
    )
    script = ReadPoolStatsScript(context)
    script.input_models = [
        InputReadModel(
            input_spec=input_spec,
            snapshot=snapshot,
            reference_df=reference_df,
            global_mean_coverage=mean,
            has_delta_g=False,
        )
        for input_spec, reference_df, mean in zip(input_specs, reference_dfs, means)
    ]

    figure = script._plot_distribution()
    metadata_pages = script._build_metadata_pages()
    rendered_text = "\n".join(
        text.get_text() for page in metadata_pages for text in page.texts
    )

    try:
        first_widths = [bar.get_width() for bar in figure.axes[0].patches]
        second_widths = [bar.get_width() for bar in figure.axes[1].patches]
        assert first_widths == pytest.approx(second_widths)
        assert figure.axes[0].get_xlim() == pytest.approx(figure.axes[1].get_xlim())
        for axis in figure.axes:
            bars_per_item = len(axis.patches) // 2
            assert sum(bar.get_height() for bar in axis.patches[:bars_per_item]) == pytest.approx(1.0)
            assert sum(bar.get_height() for bar in axis.patches[bars_per_item:]) == pytest.approx(1.0)
        pooled_mean = compute_pooled_mean_coverage(reference_dfs)
        normalized_width = 50 / pooled_mean
        assert f"pooled_mean={pooled_mean:.6g}" in rendered_text
        assert f"effective raw-read width={means[0] * normalized_width:.6g}" in rendered_text
        assert f"effective raw-read width={means[1] * normalized_width:.6g}" in rendered_text
        assert "Low coverage / Item zero: N_item=2" in rendered_text
        assert "High coverage / Item one: N_item=2" in rendered_text
    finally:
        plt.close(figure)
        for page in metadata_pages:
            plt.close(page)


def test_distribution_metadata_reports_normalization_denominator(tmp_path: Path) -> None:
    database = tmp_path / "metadata.db"
    _create_snapshot_database(database)
    script = ReadPoolStatsScript(_script_context(database, tmp_path / "plots"))
    script.input_models = [script._build_read_model_for_input(_input_spec(database))]

    pages = script._build_metadata_pages()
    rendered_text = "\n".join(text.get_text() for page in pages for text in page.texts)

    try:
        assert (
            "Synthetic input: mean read-pool coverage used as normalization denominator="
            "2.5 reads/reference"
        ) in rendered_text
    finally:
        for page in pages:
            plt.close(page)


def test_interactive_resize_happens_before_png_and_pdf_export(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    figure = plt.figure(figsize=(6.0, 5.0))
    events: list[tuple[str, tuple[float, float]]] = []

    def resize_on_show(*, block: bool) -> None:
        assert block is True
        figure.set_size_inches(8.0, 6.0)
        events.append(("show", tuple(figure.get_size_inches())))

    def record_png(*args: object, **kwargs: object) -> None:
        events.append(("png", tuple(figure.get_size_inches())))

    class RecordingPdf:
        def savefig(self, saved_figure: Figure) -> None:
            events.append(("pdf", tuple(saved_figure.get_size_inches())))

    monkeypatch.setattr(plt, "show", resize_on_show)
    monkeypatch.setattr(figure, "savefig", record_png)

    save_figure_page_and_png(
        figure,
        "interactive",
        tmp_path,
        RecordingPdf(),  # type: ignore[arg-type]
        show_figure=True,
    )

    assert events == [
        ("show", (8.0, 6.0)),
        ("png", (8.0, 6.0)),
        ("pdf", (8.0, 6.0)),
    ]


def test_script_writes_one_pdf_and_three_pngs_from_temporary_sqlite(tmp_path: Path) -> None:
    database = tmp_path / "snapshot.db"
    output_path = tmp_path / "plots"
    _create_snapshot_database(database)
    script = ReadPoolStatsScript(_script_context(database, output_path))

    script.run()

    script_output = output_path / "read_pool_stats"
    expected_outputs = {
        "read_pool_stats.pdf",
        "ref-coverage-in-read-pool_distribution.png",
        "ref-coverage-in-read-pool_vs_delta-g_scatter.png",
        "ref-coverage-in-read-pool_vs_delta-g_mean-in-bin.png",
    }
    assert {path.name for path in script_output.iterdir()} == expected_outputs
    assert all((script_output / filename).stat().st_size > 0 for filename in expected_outputs)