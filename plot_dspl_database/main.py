"""Orchestrateur principal des scripts de plotting DSPL."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence, Type

import yaml

from decoding_ref_coverage import DecodingRefCoverageScript
from decoding_run_metrics import DecodingRunMetricsScript
from protocols import PlotConfig, PlotScriptProtocol, ScriptExecutionContext


SCRIPT_REGISTRY: dict[str, Type[PlotScriptProtocol]] = {
	DecodingRefCoverageScript.script_name: DecodingRefCoverageScript,
	DecodingRunMetricsScript.script_name: DecodingRunMetricsScript,
}


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
	parser = argparse.ArgumentParser(description="Run DSPL plotting scripts from a YAML config.")
	parser.add_argument(
		"config",
		type=Path,
		help="Path to the required YAML config.",
	)
	return parser.parse_args(argv)


def _load_config_file(config_path: Path) -> dict:
	if not config_path.exists():
		raise FileNotFoundError(f"Config file not found: {config_path}")

	with config_path.open("r", encoding="utf-8") as file:
		content = yaml.safe_load(file)

	if not isinstance(content, dict):
		raise ValueError("YAML root must be a mapping.")
	return content


def main() -> None:
	args = _parse_args()
	raw_config = _load_config_file(args.config)
	validated_config = PlotConfig.model_validate(raw_config)

	config_dir = args.config.resolve().parent
	output_path = (config_dir / validated_config.output_path).resolve()

	for script_spec in validated_config.scripts:
		script_type = SCRIPT_REGISTRY.get(script_spec.name)
		if script_type is None:
			raise ValueError(f"No registered script implementation for '{script_spec.name}'.")

		context = ScriptExecutionContext(
			output_path=output_path,
			inputs=validated_config.inputs,
			script=script_spec,
		)

		script_instance = script_type(context)
		print(f"[orchestrator] running script={script_spec.name}")
		script_instance.run()


if __name__ == "__main__":
	main()