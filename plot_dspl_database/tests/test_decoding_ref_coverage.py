from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from decoding_ref_coverage import (  # noqa: E402
    compute_infinity_bin_center,
    compute_spearman_statistics,
    compute_visual_infinity_level,
)


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