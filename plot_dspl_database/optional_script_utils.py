"""Utilitaires optionnels partagés par les scripts de plotting DSPL."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from pydantic import BaseModel, ConfigDict


class PlotTextSettings(BaseModel):
	model_config = ConfigDict(extra="forbid")

	title: str | None = None
	description: str | None = None


@dataclass(frozen=True)
class MetadataSection:
	title: str
	description: str
	lines: tuple[str, ...] = ()


def resolve_database_path(database_path: Path) -> Path:
	if database_path.is_absolute():
		return database_path
	repository_root = Path(__file__).resolve().parent.parent
	return (repository_root / database_path).resolve()


def make_axes_grid(n_inputs: int, height: float = 4.8) -> tuple[Figure, list[Axes]]:
	if n_inputs <= 0:
		raise ValueError("n_inputs must be > 0.")
	fig, axes = plt.subplots(1, n_inputs, figsize=(6.0 * n_inputs, height), sharey=False)
	if n_inputs == 1:
		axes = [axes]
	return fig, list(axes)


def apply_figure_title(fig: Figure, title: str) -> None:
	fig.suptitle(title, fontsize=14, y=0.98)
	fig.tight_layout(rect=(0.02, 0.02, 0.98, 0.94))


def save_figure_page_and_png(
	fig: Figure,
	plot_name: str,
	output_dir: Path,
	pdf: PdfPages,
) -> None:
	fig.savefig(output_dir / f"{plot_name}.png", dpi=250)
	pdf.savefig(fig)
	plt.close(fig)


def build_metadata_pages(
	document_title: str,
	script_lines: Iterable[str],
	sections: Iterable[MetadataSection],
) -> list[Figure]:
	blocks: list[tuple[str, str, tuple[str, ...]]] = [
		("Script metadata", "", tuple(script_lines)),
		*((section.title, section.description, section.lines) for section in sections),
	]
	pages: list[Figure] = []
	page: Figure | None = None
	y_position = 0.0

	def start_page() -> tuple[Figure, float]:
		new_page = plt.figure(figsize=(11.69, 16.53))
		new_page.patch.set_facecolor("white")
		new_page.text(0.06, 0.965, document_title, fontsize=22, fontweight="bold", va="top")
		return new_page, 0.915

	for heading, description, lines in blocks:
		required_height = 0.055 + (0.030 if description else 0.0) + 0.024 * len(lines)
		if page is None or y_position - required_height < 0.055:
			if page is not None:
				pages.append(page)
			page, y_position = start_page()

		page.text(0.06, y_position, heading, fontsize=14, fontweight="bold", va="top")
		y_position -= 0.028
		if description:
			page.text(0.075, y_position, description, fontsize=9.5, va="top", wrap=True)
			y_position -= 0.034
		for line in lines:
			page.text(0.075, y_position, line, fontsize=9, va="top", wrap=True)
			y_position -= 0.024
		y_position -= 0.014

	if page is not None:
		pages.append(page)
	return pages
