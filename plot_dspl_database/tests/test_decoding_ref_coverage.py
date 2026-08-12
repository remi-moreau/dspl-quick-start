from __future__ import annotations

from pathlib import Path
import sys

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
)
from protocols import InputSpec, ItemSpec  # noqa: E402


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