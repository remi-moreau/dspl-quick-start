"""Utilitaires optionnels partagés par les scripts de plotting DSPL."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from textwrap import wrap
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
	same_x_scale_across_inputs: bool = False
	same_y_scale_across_inputs: bool = False


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


def make_axes_grid(
	n_inputs: int,
	height: float = 4.8,
	width_per_input: float = 6.0,
) -> tuple[Figure, list[Axes]]:
	if n_inputs <= 0:
		raise ValueError("n_inputs must be > 0.")
	if height <= 0:
		raise ValueError("height must be > 0.")
	if width_per_input <= 0:
		raise ValueError("width_per_input must be > 0.")
	fig, axes = plt.subplots(
		1,
		n_inputs,
		figsize=(width_per_input * n_inputs, height),
		sharey=False,
	)
	if n_inputs == 1:
		axes = [axes]
	return fig, list(axes)


def apply_figure_title(fig: Figure, title: str) -> None:
	fig.suptitle(title, fontsize=14, y=0.98)
	fig.tight_layout(rect=(0.02, 0.02, 0.98, 0.94))


def harmonize_axes_scales(
	axes: Iterable[Axes],
	*,
	same_x_scale: bool,
	same_y_scale: bool,
) -> None:
	populated_axes = [axis for axis in axes if axis.has_data()]
	if len(populated_axes) < 2:
		return

	if same_x_scale:
		x_limits = [axis.get_xlim() for axis in populated_axes]
		shared_x_limits = (
			min(lower for lower, _ in x_limits),
			max(upper for _, upper in x_limits),
		)
		for axis in populated_axes:
			axis.set_xlim(shared_x_limits)

	if same_y_scale:
		y_limits = [axis.get_ylim() for axis in populated_axes]
		shared_y_limits = (
			min(lower for lower, _ in y_limits),
			max(upper for _, upper in y_limits),
		)
		for axis in populated_axes:
			axis.set_ylim(shared_y_limits)


def save_figure_page_and_png(
	fig: Figure,
	plot_name: str,
	output_dir: Path,
	pdf: PdfPages,
	*,
	show_figure: bool = False,
) -> None:
	if show_figure:
		plt.figure(fig.number)
		plt.show(block=True)
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

	def wrapped_lines(text: str, width: int) -> tuple[str, ...]:
		return tuple(
			wrap(
				text,
				width=width,
				break_long_words=False,
				break_on_hyphens=False,
			)
		) or ("",)

	def start_page() -> tuple[Figure, float]:
		new_page = plt.figure(figsize=(11.69, 16.53))
		new_page.patch.set_facecolor("white")
		new_page.text(0.06, 0.965, document_title, fontsize=22, fontweight="bold", va="top")
		return new_page, 0.915

	for heading, description, lines in blocks:
		description_lines = wrapped_lines(description, 145) if description else ()
		content_lines = tuple(
			wrapped_line
			for line in lines
			for wrapped_line in wrapped_lines(line, 155)
		)
		if page is None or y_position - 0.077 < 0.055:
			if page is not None:
				pages.append(page)
			page, y_position = start_page()

		page.text(0.06, y_position, heading, fontsize=14, fontweight="bold", va="top")
		y_position -= 0.028
		for line in description_lines:
			if y_position - 0.022 < 0.055:
				pages.append(page)
				page, y_position = start_page()
				page.text(0.06, y_position, f"{heading} (continued)", fontsize=14, fontweight="bold", va="top")
				y_position -= 0.028
			page.text(0.075, y_position, line, fontsize=9.5, va="top")
			y_position -= 0.022
		if description_lines:
			y_position -= 0.012
		for line in content_lines:
			if y_position - 0.024 < 0.055:
				pages.append(page)
				page, y_position = start_page()
				page.text(0.06, y_position, f"{heading} (continued)", fontsize=14, fontweight="bold", va="top")
				y_position -= 0.028
			page.text(0.075, y_position, line, fontsize=9, va="top")
			y_position -= 0.024
		y_position -= 0.014

	if page is not None:
		pages.append(page)
	return pages
