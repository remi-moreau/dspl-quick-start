from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from decoding_ref_coverage import (  # noqa: E402
    DecodingRefCoverageScript,
    InputReadModel,
    UnderdecodedPlotSettings,
    compute_infinity_bin_center,
    compute_spearman_statistics,
    compute_underdecoded_probability_bins,
    compute_visual_infinity_level,
    format_distribution_statistics,
)
from protocols import (  # noqa: E402
    InputSpec,
    ItemSpec,
    PlotSpec,
    ScriptExecutionContext,
    ScriptSpec,
)


def _coverage_input_model(
    name: str,
    finite_coverages: list[float],
    *,
    drop_out_count: int = 0,
) -> InputReadModel:
    input_spec = InputSpec(
        input_id=1,
        name=name,
        database=Path("unused.db"),
        exp_id="test-exp",
        read_pool_id="test-pool",
        decoding_run_label="test-label",
        items=[ItemSpec(item_id=0, name="Item zero")],
    )
    decoded_df = pd.DataFrame(
        {
            "item_id": [0] * len(finite_coverages),
            "mean_positive_coverage": finite_coverages,
        }
    )
    never_decoded_df = pd.DataFrame({"item_id": [0] * drop_out_count})
    reference_df = pd.concat(
        [decoded_df, never_decoded_df.assign(mean_positive_coverage=np.nan)],
        ignore_index=True,
    )
    rows_df = pd.DataFrame(
        {
            "dec_run_id": ["run"] * len(finite_coverages),
            "item_id": [0] * len(finite_coverages),
            "count_used_for_consensus": finite_coverages,
        }
    )
    empty_df = pd.DataFrame()
    return InputReadModel(
        input_spec=input_spec,
        rows_df=rows_df,
        reference_level_df=reference_df,
        accepted_reference_level_df=empty_df,
        underdecoded_reference_level_df=empty_df,
        never_decoded_reference_level_df=never_decoded_df,
        decoded_at_least_once_reference_level_df=decoded_df,
        n_runs_total_labeled=1,
        n_runs_displayed=1,
        has_delta_g=False,
    )


def _coverage_script(
    plot_name: str,
    input_models: list[InputReadModel],
    settings: dict[str, object] | None = None,
) -> DecodingRefCoverageScript:
    context = ScriptExecutionContext(
        output_path=Path("plots"),
        inputs=[model.input_spec for model in input_models],
        script=ScriptSpec(
            name="decoding_ref-coverage",
            plot_settings=[PlotSpec(name=plot_name, settings=settings or {})],
        ),
    )
    script = DecodingRefCoverageScript(context)
    script.input_models = input_models
    return script


def test_infinity_is_one_bin_after_last_finite_center() -> None:
    infinity_x = compute_infinity_bin_center(pd.Series([1.01, 1.26, 1.77]), bin_width=0.1)

    assert infinity_x == pytest.approx(1.85)


def test_shared_infinity_uses_the_global_finite_maximum() -> None:
    all_inputs = pd.concat(
        [pd.Series([1.01, 1.77]), pd.Series([1.20, 2.34])],
        ignore_index=True,
    )

    assert compute_infinity_bin_center(all_inputs, bin_width=0.1) == pytest.approx(2.45)
    assert compute_visual_infinity_level(all_inputs, factor=1.08) == pytest.approx(2.5272)


def test_spearman_statistics_use_only_finite_pairs() -> None:
    rows = pd.DataFrame(
        {
            "delta_g": [1.0, 2.0, 3.0, None],
            "mean_positive_coverage": [10.0, 20.0, 30.0, 40.0],
        }
    )

    n_points, rho, p_value = compute_spearman_statistics(rows)

    assert n_points == 3
    assert rho == pytest.approx(1.0)
    assert p_value == pytest.approx(0.0)


def test_spearman_statistics_are_unavailable_for_constant_values() -> None:
    rows = pd.DataFrame(
        {
            "delta_g": [1.0, 2.0, 3.0],
            "mean_positive_coverage": [5.0, 5.0, 5.0],
        }
    )

    assert compute_spearman_statistics(rows) == (3, None, None)


def test_distribution_statistics_distinguish_population_sigma_and_sample_std() -> None:
    statistics = format_distribution_statistics(pd.Series([1.0, 2.0, 3.0]))

    assert statistics == "n=3, mu=2, sigma=0.8165, median=2, std=1"


def test_distribution_statistics_are_metadata_and_not_image_legend() -> None:
    input_model = _coverage_input_model("Input A", [1.0, 2.0, 3.0], drop_out_count=2)
    image_script = _coverage_script(
        "ref-coverage-at-image-decoding_distribution",
        [input_model],
    )
    ref_script = _coverage_script(
        "ref-coverage-at-ref-decoding_distribution",
        [input_model],
    )

    image_lines = image_script._distribution_metadata_lines(
        "ref-coverage-at-image-decoding_distribution"
    )
    ref_lines = ref_script._distribution_metadata_lines(
        "ref-coverage-at-ref-decoding_distribution"
    )
    figure = image_script._plot_ref_coverage_at_image_decoding_distribution()

    try:
        assert image_lines[-1] == "  Item zero: n=3, mu=2, sigma=0.8165, median=2, std=1"
        assert ref_lines[-1].endswith("median=2, std=1, drop_outs=2")
        assert figure.axes[0].get_legend_handles_labels()[1] == ["Item zero"]
    finally:
        plt.close(figure)


def test_ref_distribution_shared_x_scale_uses_identical_ticks() -> None:
    input_models = [
        _coverage_input_model("Short range", [1.2, 2.2], drop_out_count=1),
        _coverage_input_model("Long range", [1.2, 10.2, 20.2], drop_out_count=1),
    ]
    script = _coverage_script(
        "ref-coverage-at-ref-decoding_distribution",
        input_models,
        {
            "same_x_scale_across_inputs": True,
            "bin_width": 1.0,
            "max_xticks": 4,
        },
    )

    figure = script._plot_ref_coverage_at_ref_decoding_distribution()

    try:
        first_ticks = figure.axes[0].get_xticks()
        second_ticks = figure.axes[1].get_xticks()
        assert first_ticks.tolist() == pytest.approx(second_ticks.tolist())
        assert max(first_ticks[:-1]) > 2.2
        assert figure.axes[0].get_xlim() == pytest.approx(figure.axes[1].get_xlim())
    finally:
        plt.close(figure)


def test_delta_g_metadata_reports_global_spearman_across_items() -> None:
    item_zero_reference_df = pd.DataFrame(
        {
            "item_id": [0, 0, 0, 0, 0, 0],
            "delta_g": [0.2, 0.3, 1.2, 1.3, 2.2, 2.3],
            "zero_ratio": [0.0, 0.1, 0.2, 0.8, 0.9, 1.0],
        }
    )
    item_one_reference_df = item_zero_reference_df.copy()
    item_one_reference_df["item_id"] = 1
    reference_df = pd.concat(
        [item_zero_reference_df, item_one_reference_df],
        ignore_index=True,
    )
    settings = UnderdecodedPlotSettings(
        max_zero_run_ratio_for_inclusion=0.5,
        delta_g_precision=1.0,
        min_points_per_bin=2,
    )
    binned = compute_underdecoded_probability_bins(
        reference_df,
        item_ids=[0, 1],
        max_zero_run_ratio_for_inclusion=settings.max_zero_run_ratio_for_inclusion,
        delta_g_precision=settings.delta_g_precision,
        min_points_per_bin=settings.min_points_per_bin,
    )
    input_spec = InputSpec(
        input_id=1,
        name="Synthetic input",
        database=Path("unused.db"),
        exp_id="test-exp",
        read_pool_id="test-pool",
        decoding_run_label="test-label",
        items=[
            ItemSpec(item_id=0, name="Item zero"),
            ItemSpec(item_id=1, name="Item one"),
        ],
    )
    empty_df = pd.DataFrame()
    decoded_df = pd.DataFrame(
        {
            "item_id": [0, 0, 1, 1],
            "delta_g": [1.0, 2.0, 3.0, 4.0],
            "mean_positive_coverage": [10.0, 20.0, 30.0, 40.0],
        }
    )
    input_model = InputReadModel(
        input_spec=input_spec,
        rows_df=empty_df,
        reference_level_df=reference_df,
        accepted_reference_level_df=empty_df,
        underdecoded_reference_level_df=empty_df,
        never_decoded_reference_level_df=empty_df,
        decoded_at_least_once_reference_level_df=decoded_df,
        n_runs_total_labeled=0,
        n_runs_displayed=0,
        has_delta_g=True,
    )
    script = object.__new__(DecodingRefCoverageScript)
    script.input_models = [input_model]

    scatter_metadata_lines = script._spearman_metadata_lines()
    underdecoded_metadata_lines = script._underdecoded_probability_spearman_metadata_lines(settings)

    assert binned["underdecoded_probability"].tolist() == pytest.approx(
        [0.0, 0.5, 1.0, 0.0, 0.5, 1.0]
    )
    assert scatter_metadata_lines == [
        "Global Spearman - Synthetic input: n=4, rho=1, p_value=0"
    ]
    assert underdecoded_metadata_lines == [
        "Global Spearman (delta G bin center vs underdecoded probability) - "
        "Synthetic input: n_bins=6, rho=1, p_value=0"
    ]