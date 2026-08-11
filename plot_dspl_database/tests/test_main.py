from __future__ import annotations

from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from main import _parse_args  # noqa: E402


def test_config_path_is_required() -> None:
    with pytest.raises(SystemExit) as error:
        _parse_args([])

    assert error.value.code == 2


def test_config_path_is_parsed() -> None:
    args = _parse_args(["custom.yml"])

    assert args.config == Path("custom.yml")