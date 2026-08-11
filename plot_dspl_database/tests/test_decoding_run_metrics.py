from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from pydantic import ValidationError

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from decoding_run_metrics import (  # noqa: E402
    DecodingRunMetricsScript,
    DecodingRunMetricsSettings,
    convert_at_decoding_to_display_units,
)
from optional_script_utils import (  # noqa: E402
    MetadataSection,
    build_metadata_pages,
    harmonize_axes_scales,
)
from protocols import PlotSpec, ScriptExecutionContext, ScriptSpec  # noqa: E402


def _script_with_settings(**settings: object) -> DecodingRunMetricsScript:
    context = ScriptExecutionContext(
        output_path=Path("plots"),
        inputs=[],
        script=ScriptSpec(
            name="decoding_run-metrics",
            plot_settings=[PlotSpec(name="psnr_vs_run-progression")],
            script_settings=settings,
        ),
    )
    return DecodingRunMetricsScript(context)


def test_script_settings_are_strict_and_validate_curve_sizes() -> None:
    with pytest.raises(ValidationError):
        DecodingRunMetricsSettings.model_validate({"unknown": True})
    with pytest.raises(ValidationError):
        DecodingRunMetricsSettings.model_validate({"grid_n_points": 1})
    with pytest.raises(ValidationError):
        DecodingRunMetricsSettings.model_validate({"min_points_per_curve": 1})


def test_prepare_curve_filters_only_the_requested_metric() -> None:
    script = _script_with_settings()
    rows = pd.DataFrame(
        {
            "coverage": [0.0, 1.0, 2.0],
            "psnr": [np.nan, 20.0, 21.0],
            "hamming_distance_normalized": [0.4, 0.2, 0.1],
        }
    )

    hamming_curve = script._prepare_xy_curve(
        rows,
        "coverage",
        "hamming_distance_normalized",
    )
    psnr_curve = script._prepare_xy_curve(rows, "coverage", "psnr")

    assert hamming_curve is not None
    assert hamming_curve[0].tolist() == [0.0, 1.0, 2.0]
    assert psnr_curve is not None
    assert psnr_curve[0].tolist() == [1.0, 2.0]


def test_interpolation_does_not_extrapolate_outside_run_support() -> None:
    result = DecodingRunMetricsScript._interpolate_on_grid(
        np.array([1.0, 2.0]),
        np.array([10.0, 20.0]),
        np.array([0.0, 1.0, 1.5, 2.0, 3.0]),
    )

    assert np.isnan(result[0])
    assert result[1:4].tolist() == [10.0, 15.0, 20.0]
    assert np.isnan(result[4])


def test_at_decoding_statistics_use_sample_std_and_report_denominator() -> None:
    line = DecodingRunMetricsScript._format_statistics(
        "JPEGDNA-control",
        pd.Series([10.0, 14.0]),
        decoded_count=2,
        labelled_count=3,
    )

    assert "n=2" in line
    assert "decoded/labelled=2/3" in line
    assert "mean=12" in line
    assert "std=2.828" in line
    assert "min=10" in line
    assert "max=14" in line


def test_at_decoding_display_units_convert_duration_and_preserve_coverage() -> None:
    source = pd.DataFrame(
        {
            "estimated_sequencing_duration_at_decoding": [60.0, 150.0],
            "perfectly_decoded_payload_ratio_at_decoding": [0.25, 1.0],
            "coverage_at_decoding": [3.5, 7.0],
        }
    )

    converted = convert_at_decoding_to_display_units(source)

    assert converted["estimated_sequencing_duration_at_decoding"].tolist() == [1.0, 2.5]
    assert converted["perfectly_decoded_payload_ratio_at_decoding"].tolist() == [25.0, 100.0]
    assert converted["coverage_at_decoding"].tolist() == [3.5, 7.0]
    assert source["estimated_sequencing_duration_at_decoding"].tolist() == [60.0, 150.0]


def test_metadata_builder_creates_multiple_pages_without_dropping_sections() -> None:
    sections = [
        MetadataSection(
            title=f"Figure {index}",
            description="Description",
            lines=tuple(f"line {line}" for line in range(8)),
        )
        for index in range(8)
    ]

    pages = build_metadata_pages("decoding_run-metrics", ["script line"], sections)
    rendered_text = "\n".join(text.get_text() for page in pages for text in page.texts)

    try:
        assert len(pages) > 1
        for index in range(8):
            assert f"Figure {index}" in rendered_text
    finally:
        for page in pages:
            plt.close(page)


def test_axes_scales_can_be_harmonized_independently() -> None:
    figure, axes = plt.subplots(1, 2)
    axes[0].plot([0.0, 1.0], [10.0, 20.0])
    axes[1].plot([2.0, 4.0], [-5.0, 5.0])
    original_y_limits = [axis.get_ylim() for axis in axes]

    harmonize_axes_scales(axes, same_x_scale=True, same_y_scale=False)

    try:
        assert axes[0].get_xlim() == axes[1].get_xlim()
        assert axes[0].get_ylim() == original_y_limits[0]
        assert axes[1].get_ylim() == original_y_limits[1]
    finally:
        plt.close(figure)
