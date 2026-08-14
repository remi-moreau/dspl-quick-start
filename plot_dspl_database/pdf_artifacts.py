"""Composition des artefacts PDF produits par les scripts de plotting."""

from __future__ import annotations

from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Sequence


def merge_pdf_files(input_paths: Sequence[Path], output_path: Path) -> None:
	if not input_paths:
		raise ValueError("At least one input PDF is required for merging.")

	missing_paths = [path for path in input_paths if not path.is_file()]
	if missing_paths:
		missing = ", ".join(str(path) for path in missing_paths)
		raise FileNotFoundError(f"Cannot merge missing script PDF(s): {missing}")

	from pypdf import PdfWriter

	output_path.parent.mkdir(parents=True, exist_ok=True)
	writer = PdfWriter()
	for input_path in input_paths:
		writer.append(input_path)

	temporary_path: Path | None = None
	try:
		with NamedTemporaryFile(
			mode="wb",
			dir=output_path.parent,
			prefix=f".{output_path.name}.",
			suffix=".tmp",
			delete=False,
		) as temporary_file:
			temporary_path = Path(temporary_file.name)
			writer.write(temporary_file)
		temporary_path.replace(output_path)
	finally:
		writer.close()
		if temporary_path is not None and temporary_path.exists():
			temporary_path.unlink()