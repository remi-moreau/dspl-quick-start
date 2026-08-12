"""Contracts metier pour l'orchestrateur de plots DSPL.

Ce module contient:
- les modeles pydantic pour valider/normaliser plot_config.yml
- les protocoles d'interface que chaque script de plotting doit implementer
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, ClassVar, Protocol

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


SUPPORTED_SCRIPTS: set[str] = {
	"decoding_run-metrics",
	"decoding_ref-coverage",
	"read_pool_stats",
}


class ItemSpec(BaseModel):
	model_config = ConfigDict(extra="forbid")

	item_id: int
	name: str


class InputSpec(BaseModel):
	model_config = ConfigDict(extra="forbid")

	input_id: int
	name: str
	database: Path
	exp_id: str
	read_pool_id: str
	decoding_run_label: str
	items: list[ItemSpec]


class PlotSpec(BaseModel):
	model_config = ConfigDict(extra="forbid")

	name: str
	settings: dict[str, Any] = Field(default_factory=dict)


class ScriptSpec(BaseModel):
	model_config = ConfigDict(extra="forbid")

	name: str
	plot_settings: list[PlotSpec]
	script_settings: dict[str, Any] = Field(default_factory=dict)

	@field_validator("name")
	@classmethod
	def _validate_known_script_name(cls, script_name: str) -> str:
		if script_name not in SUPPORTED_SCRIPTS:
			allowed = ", ".join(sorted(SUPPORTED_SCRIPTS))
			raise ValueError(f"Unknown script '{script_name}'. Allowed values: {allowed}")
		return script_name

	@field_validator("plot_settings")
	@classmethod
	def _validate_plot_settings(cls, plot_settings: list[PlotSpec]) -> list[PlotSpec]:
		if not plot_settings:
			raise ValueError("plot_settings must contain at least one figure.")

		plot_names = [plot.name for plot in plot_settings]
		duplicate_names = sorted({name for name in plot_names if plot_names.count(name) > 1})
		if duplicate_names:
			raise ValueError(
				"plot_settings contains duplicate figure name(s): "
				+ ", ".join(duplicate_names)
			)
		return plot_settings


class PlotConfig(BaseModel):
	model_config = ConfigDict(extra="forbid")

	output_path: Path
	merge_script_pdfs: bool = False
	show_figures: bool = False
	figure_width_per_input: float = 6.0
	figure_height: float = 5.0
	inputs: list[InputSpec]
	scripts: list[ScriptSpec]

	@field_validator("figure_width_per_input", "figure_height")
	@classmethod
	def _validate_figure_dimension(cls, value: float) -> float:
		if value <= 0:
			raise ValueError("Global figure dimensions must be > 0.")
		return value

	@model_validator(mode="before")
	@classmethod
	def _normalize_yaml_shape(cls, data: Any) -> Any:
		if not isinstance(data, dict):
			raise ValueError("Top-level config must be a mapping.")

		normalized = dict(data)
		normalized["inputs"] = cls._normalize_inputs(data.get("inputs"))
		normalized["scripts"] = cls._normalize_scripts(data.get("scripts"))
		return normalized

	@staticmethod
	def _normalize_inputs(raw_inputs: Any) -> list[dict[str, Any]]:
		if not isinstance(raw_inputs, list):
			raise ValueError("'inputs' must be a list.")

		normalized_inputs: list[dict[str, Any]] = []
		for raw_input in raw_inputs:
			input_id_raw, payload = PlotConfig._extract_anchor_and_payload(
				raw_entry=raw_input,
				entry_kind="input",
				payload_field_names={
					"name",
					"database",
					"exp_id",
					"read_pool_id",
					"items",
					"decoding_run_label",
				},
			)

			normalized_payload = dict(payload)
			normalized_payload["input_id"] = int(input_id_raw)
			normalized_payload["items"] = PlotConfig._normalize_items(payload.get("items"))
			normalized_inputs.append(normalized_payload)
		return normalized_inputs

	@staticmethod
	def _normalize_items(raw_items: Any) -> list[dict[str, Any]]:
		if not isinstance(raw_items, list):
			raise ValueError("Each input 'items' field must be a list.")

		normalized_items: list[dict[str, Any]] = []
		for raw_item in raw_items:
			item_id_raw, payload = PlotConfig._extract_anchor_and_payload(
				raw_entry=raw_item,
				entry_kind="item",
				payload_field_names={"name"},
			)
			normalized_items.append({"item_id": int(item_id_raw), **payload})
		return normalized_items

	@staticmethod
	def _normalize_scripts(raw_scripts: Any) -> list[dict[str, Any]]:
		if not isinstance(raw_scripts, list):
			raise ValueError("'scripts' must be a list.")

		normalized_scripts: list[dict[str, Any]] = []
		for raw_script in raw_scripts:
			script_name, payload = PlotConfig._extract_anchor_and_payload(
				raw_entry=raw_script,
				entry_kind="script",
				payload_field_names={"plot_settings", "script_settings"},
			)

			normalized_scripts.append(
				{
					"name": str(script_name),
					"plot_settings": PlotConfig._normalize_plot_settings(payload.get("plot_settings")),
					"script_settings": payload.get("script_settings") or {},
				}
			)
		return normalized_scripts

	@staticmethod
	def _normalize_plot_settings(raw_plot_settings: Any) -> list[dict[str, Any]]:
		if not isinstance(raw_plot_settings, list):
			raise ValueError("Each script 'plot_settings' field must be a list.")

		normalized_plot_settings: list[dict[str, Any]] = []
		for raw_plot in raw_plot_settings:
			if isinstance(raw_plot, str):
				normalized_plot_settings.append({"name": raw_plot, "settings": {}})
				continue

			plot_name, plot_payload = PlotConfig._extract_anchor_and_payload(
				raw_entry=raw_plot,
				entry_kind="plot",
				payload_field_names={"settings"},
			)

			normalized_plot_settings.append({"name": str(plot_name), "settings": plot_payload})
		return normalized_plot_settings

	@staticmethod
	def _extract_anchor_and_payload(
		raw_entry: Any,
		entry_kind: str,
		payload_field_names: set[str] | None = None,
	) -> tuple[Any, dict[str, Any]]:
		if not isinstance(raw_entry, dict) or not raw_entry:
			raise ValueError(
				f"Each {entry_kind} entry must be a mapping with an anchor key: '- <id_or_name>:'."
			)

		if len(raw_entry) == 1:
			anchor_key, payload = next(iter(raw_entry.items()))
			if payload is None:
				return anchor_key, {}
			if not isinstance(payload, dict):
				raise ValueError(f"Each {entry_kind} payload must be a mapping.")
			return anchor_key, dict(payload)

		# Supports YAML style:
		# - 1:
		#   name: ...
		#   ...
		payload_field_names = payload_field_names or set()
		anchor_candidates = [
			key
			for key, value in raw_entry.items()
			if key not in payload_field_names and value is None
		]
		if len(anchor_candidates) != 1:
			raise ValueError(
				f"Each {entry_kind} entry must contain exactly one anchor key with null value."
			)

		anchor_key = anchor_candidates[0]
		payload = {key: value for key, value in raw_entry.items() if key != anchor_key}
		return anchor_key, payload


class ScriptExecutionContext(BaseModel):
	model_config = ConfigDict(extra="forbid")

	output_path: Path
	show_figures: bool = False
	figure_width_per_input: float = 6.0
	figure_height: float = 5.0
	inputs: list[InputSpec]
	script: ScriptSpec

	@field_validator("figure_width_per_input", "figure_height")
	@classmethod
	def _validate_figure_dimension(cls, value: float) -> float:
		if value <= 0:
			raise ValueError("Global figure dimensions must be > 0.")
		return value


class PlotScriptProtocol(Protocol):
	script_name: ClassVar[str]

	def __init__(self, context: ScriptExecutionContext):
		"""Initialize script with validated orchestration context."""

	def run(self) -> None:
		"""Run all figures required by the script."""