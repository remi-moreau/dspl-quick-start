from __future__ import annotations

from pathlib import Path
import sys

import pytest
import yaml
from pypdf import PdfReader, PdfWriter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import main as main_module  # noqa: E402
from main import _parse_args  # noqa: E402
from pdf_artifacts import MERGED_PDF_FILENAME  # noqa: E402
from protocols import PlotConfig  # noqa: E402


def _minimal_raw_config() -> dict:
    return {
        "output_path": "plots",
        "inputs": [
            {
                0: {
                    "name": "Input",
                    "database": "database.db",
                    "exp_id": "experiment",
                    "read_pool_id": "pool",
                    "decoding_run_label": "run",
                    "items": [{0: {"name": "Item"}}],
                }
            }
        ],
        "scripts": [
            {
                "read_pool_stats": {
                    "plot_settings": ["ref-coverage-in-read-pool_distribution"]
                }
            }
        ],
    }


def test_config_path_is_required() -> None:
    with pytest.raises(SystemExit) as error:
        _parse_args([])

    assert error.value.code == 2


def test_config_path_is_parsed() -> None:
    args = _parse_args(["custom.yml"])

    assert args.config == Path("custom.yml")


def test_pdf_merge_is_disabled_by_default() -> None:
    config = PlotConfig.model_validate(_minimal_raw_config())

    assert config.merge_script_pdfs is False


def test_pdf_merge_can_be_enabled() -> None:
    raw_config = _minimal_raw_config()
    raw_config["merge_script_pdfs"] = True

    config = PlotConfig.model_validate(raw_config)

    assert config.merge_script_pdfs is True


def test_all_plot_config_variants_declare_boolean_pdf_merge() -> None:
    repository_root = Path(__file__).resolve().parents[2]
    config_paths = sorted(repository_root.glob("plot_config_*.yml"))

    assert config_paths
    for config_path in config_paths:
        raw_config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
        config = PlotConfig.model_validate(raw_config)
        assert isinstance(raw_config.get("merge_script_pdfs"), bool), config_path.name
        assert config.merge_script_pdfs is raw_config["merge_script_pdfs"]


@pytest.mark.parametrize("merge_script_pdfs", [False, True])
def test_main_optionally_merges_script_pdfs_in_yaml_order(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    merge_script_pdfs: bool,
) -> None:
    class FakePlotScript:
        def __init__(self, context):
            self.context = context

        def run(self) -> None:
            page_widths = {
                "decoding_ref-coverage": 101,
                "read_pool_stats": 202,
            }
            output_dir = self.context.output_path / self.context.script.name
            output_dir.mkdir(parents=True, exist_ok=True)
            writer = PdfWriter()
            writer.add_blank_page(
                width=page_widths[self.context.script.name],
                height=300,
            )
            with (output_dir / f"{self.context.script.name}.pdf").open("wb") as file:
                writer.write(file)
            writer.close()

    raw_config = _minimal_raw_config()
    raw_config["merge_script_pdfs"] = merge_script_pdfs
    raw_config["scripts"] = [
        {
            "decoding_ref-coverage": {
                "plot_settings": ["ref-coverage-at-ref-decoding_distribution"]
            }
        },
        {
            "read_pool_stats": {
                "plot_settings": ["ref-coverage-in-read-pool_distribution"]
            }
        },
    ]
    config_path = tmp_path / "plot_config.yml"
    config_path.write_text(yaml.safe_dump(raw_config), encoding="utf-8")
    monkeypatch.setattr(
        main_module,
        "SCRIPT_REGISTRY",
        {
            "decoding_ref-coverage": FakePlotScript,
            "read_pool_stats": FakePlotScript,
        },
    )
    monkeypatch.setattr(sys, "argv", ["main.py", str(config_path)])

    main_module.main()

    output_path = tmp_path / "plots"
    merged_pdf_path = output_path / MERGED_PDF_FILENAME
    assert merged_pdf_path.exists() is merge_script_pdfs
    assert (output_path / "decoding_ref-coverage" / "decoding_ref-coverage.pdf").is_file()
    assert (output_path / "read_pool_stats" / "read_pool_stats.pdf").is_file()
    if merge_script_pdfs:
        reader = PdfReader(merged_pdf_path)
        page_widths = [float(page.mediabox.width) for page in reader.pages]
        assert page_widths == [101.0, 202.0]